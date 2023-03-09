### Packages
library(MASS)
library(tidyverse)
library(survey)
#devtools::install_github("csterbenz1/KBAL", ref = "cat_kernel")
library(kbal)
library(parallel)
library(knitr)
library(glmnet)
library(tictoc)

###### SET PARAMS  ###############
set.seed(9345876)

if(detectCores() > 10) {
  path_data = "/home/csterbenz/Data/"
  cores_saved = 10
} else if(detectCores() != 4) {
  path_data= "/Users/Ciara_1/Dropbox/kpop/Updated/application/data/" 
  cores_saved = 6
} else {
    path_data= "/Users/Ciara/Dropbox/kpop/Updated/application/data/" 
    cores_saved = 2
}

POPW = FALSE
eval_kpop = TRUE
TEST = FALSE # to run with a linear kernel so it's way faster; UPDATE: errors catch this as mistake and prevent
tolerance = 1e-4
maxit = 500
#both for runtime
increment = 5
#for more complex dgp
# min_num_dims = 150
# max_num_dims = 500 
min_num_dims = NULL
max_num_dims = NULL 
SAVE = TRUE
##### Central Params to adjust
n_sample = 2500
simple_selection_model = TRUE
nsims = (detectCores()-cores_saved)*13
nsims
###################### Formulas ################
formula_rake_demos_noeduc <- ~recode_age_bucket + recode_female + recode_race + 
  recode_region + recode_pid_3way

#updated to include 6 way edu
formula_rake_demos_weduc <- ~recode_age_bucket + recode_female + 
  recode_race + recode_region + recode_educ + recode_pid_3way

formula_rake_all_vars <- ~recode_age_bucket + recode_female + 
  recode_race + recode_region + recode_pid_3way + recode_educ +
  recode_income_5way + recode_relig_6way + recode_born + recode_attndch_4way

formula_ps <- ~recode_age_3way + recode_female + recode_race +
  recode_region + recode_educ_wh_3way + recode_pid_3way

formula_ps_reduc <- ~recode_age_3way + recode_female + 
  recode_race + recode_region + recode_pid_3way  +
  recode_income_3way + recode_born + recode_educ_wh_3way

# #let's coarsen income and religion and attend church for all:
formula_ps_all <- ~recode_age_3way + recode_female +
  recode_race + recode_region + recode_pid_3way + recode_educ_3way +
  recode_income_3way + recode_relig_6way + recode_born + recode_attndch_bin


#### (4) Selection model (also rake_demo_truth)
# selection_model = as.formula(~recode_female:recode_pid_3way +
#                                  recode_age:recode_pid_3way + 
#                                  recode_race:recode_region:recode_educ_wh_3way:recode_pid_3way +
#                                  recode_age + 
#                                  I(recode_age^2))

create_targets <- function (target_design, target_formula) {
  target_mf <- model.frame(target_formula, model.frame(target_design))
  target_mm <- model.matrix(target_formula, target_mf)
  wts <- weights(target_design)
  
  return(colSums(target_mm * wts) / sum(wts))
}

### Post-stratification function
## For now assumes that strata variable is already created and in
## the data set and called "strata"
postStrat <- function(survey, pop_counts, pop_w_col, strata_pass) {
  survey_counts <- survey %>%
    group_by(!!as.symbol(strata_pass)) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    mutate(w_survey = n / sum(n))
  
  pop_counts <- pop_counts %>%
    rename(w_pop = matches(pop_w_col))
  
  post_strat <- pop_counts %>%
    left_join(survey_counts, by = strata_pass) %>%
    filter(!is.na(w_survey)) %>%
    ## Normalizes back to 1 after dropping
    ## empty cells
    mutate(w_pop = w_pop * 1/sum(w_pop),
           w = w_pop / w_survey) %>%
    dplyr::select(!!as.symbol(strata_pass), w)
  
  survey <- survey %>%
    left_join(post_strat)
  
  return(survey)
}



############# Load Data #####################
#these data have been cleaned already see app_modeled for how it was done
## Load Pew 
pew <- readRDS(paste0(path_data, "pew_lasso_061021.rds"))

#adjusting size of pew here:
# n_sample = nrow(pew)
# samp_units = sample(1:nrow(pew), n_sample, replace = F)
# pew = pew[samp_units, ]


### Load Target Data
cces <- readRDS(paste0(path_data, "cces_lasso_061021.rds"))

#################### Define Outcome #########
#cces$outcome = cces$diff_cces_on_cces
#cces = cces %>% mutate(outcome = rbinom(nrow(cces), 1, mod_cces_on_cces_pR))


#this is wrong lol
# cces$outcome =plogis((var(cces$diff_cces_on_cces) + abs(cces$diff_cces_on_cces))*sign(cces$diff_cces_on_cces))
# cces$outcome =cces$diff_cces_on_cces + rnorm(nrow(cces), 0 , sqrt(var(cces$diff_cces_on_cces)))
# #so this is actually a bit annoying bc the difference can be on a -1 to 1 scale but plogis pushes us to 0-1 reducing the var
# #so let's remake the difference with increased var:
# cces$pD_upvar = plogis(cces$mod_cces_on_cces_pD + rnorm(nrow(cces), 0 , sqrt(var(cces$mod_cces_on_cces_pD))))
# cces$pR_upvar = plogis(cces$mod_cces_on_cces_pR + rnorm(nrow(cces), 0 , sqrt(var(cces$mod_cces_on_cces_pR))))
# cces$outcome = cces$pD_upvar - cces$pR_upvar
# mean(cces$outcome)
# # cces$diff_cces_on_cces[1:10]
# # cces$mod_cces_on_cces_pD[1:10] - cces$mod_cces_on_cces_pR[1:10]
# # cces$outcome[1:10]
# var(cces$outcome)
# var(cces$diff_cces_on_cces)
# range(cces$diff_cces_on_cces)
# range(cces$outcome)
cces$outcome = cces$diff_cces_on_cces
#cces$outcome = cces$diff_cces_on_cces

######### Make STRATA variable in CCES and Pew ############
cces <- bind_cols(cces, cces %>% 
                    unite("strata", all.vars(formula_ps), remove = FALSE) %>%
                    unite("strata_reduc", all.vars(formula_ps_reduc), 
                          remove = FALSE) %>%
                    unite("strata_all", all.vars(formula_ps_all)) %>%
                    dplyr::select(strata, strata_reduc, strata_all))

pew <- bind_cols(pew, pew %>% 
                   unite("strata", all.vars(formula_ps), remove = FALSE) %>%
                   unite("strata_reduc", all.vars(formula_ps_reduc), 
                         remove = FALSE) %>%
                   unite("strata_all", all.vars(formula_ps_all)) %>%
                   dplyr::select(strata, strata_reduc, strata_all))


 ##################### LASSO: Selection #############################


#income, religion, pid x race
if(simple_selection_model) {
    #first attempt to make this worse
    selection_model = as.formula(~recode_pid_3way:poly(recode_age, 2) + 
                                     recode_female:recode_pid_3way)
    #first attempt to make this worse
    #center age then square
    pew = pew %>% mutate(centered_age = scale(recode_age, scale = F))
    cces = cces %>% mutate(centered_age = scale(recode_age, scale = F))
    
    selection_model = as.formula(~recode_pid_3way*poly(centered_age, 2) + 
                                     recode_female*recode_pid_3way)
    

} else {
    selection_model = as.formula(~recode_female:recode_pid_3way + 
                                     recode_age:recode_pid_3way +
                                     #adds a bit mroe bias to edu+d
                                     recode_pid_race + 
                                     recode_race_reg_wh_educ +
                                     recode_educ_wh_3way +
                                     poly(recode_age, 3))
}

# Stack data with S = 1 indicating Pew
stack_data <- data.frame(bind_rows(pew, cces), 
                         S = c(rep(1, nrow(pew)), rep(0, nrow(cces))))


mod <- model.matrix(selection_model, data = stack_data)
nrow(mod)
## Remove columns where Pew missing strata
ncol(mod)
mod <- mod[, apply(mod[stack_data$S == 1, ], 2, sum) != 0]
## Remove columns where CCES missing Strata
mod <- mod[, apply(mod[stack_data$S == 0, ], 2, sum) != 0]
ncol(mod)

#lambda from 10 fold default CV
#had to remove the intercept in mod and have glmnet add it cause other wise even with intercept false it was adding two intercepts
lasso_lambda <- cv.glmnet(x= mod[,-1], 
                          y = as.matrix(stack_data$S),
                          alpha = 1,
                          family = "binomial",
                          intercept = TRUE)
lasso_lambda$lambda.min
lasso_include <- glmnet(x= mod[,-1], 
                        y = as.matrix(stack_data$S),
                        alpha = 1,
                        lambda = lasso_lambda$lambda.min,
                        weights = if(POPW){ c(rep(1, nrow(pew)),
                                              cces$commonweight_vv_post)} else {NULL},
                        family = "binomial",
                        intercept = TRUE)

lasso_include_coefs <- coef(lasso_include) 
res <- as.matrix(lasso_include_coefs)
#View(res)
ncol(mod)
sum(lasso_include_coefs == 0)

lasso_include_coefs
#probability of being in pew for cces 
lasso_pinclude = predict(lasso_include,
                         s= lasso_lambda$lambda.min,
                         type = "response",
                         newx = mod[stack_data$S == 0,-1])

coefs <- as.matrix(lasso_include_coefs)
coefs

p_include <- lasso_pinclude
sum(p_include)


########### Sampling Inclusion Results
summary(p_include)

#p_include = p_include*(n_sample/sum(p_include))
#sum(p_include)

cor(cces$diff_cces_on_cces, p_include)
cor(cces$outcome, p_include)

#cor(cces$outcome, cces$recode_age)
#cor(cces$outcome, cces$recode_age*cces$recode_age)
cor(cces$outcome, cces$centered_age)
cor(cces$outcome, cces$centered_age*cces$centered_age)

cor(cces$diff_cces_on_cces, cces$centered_age*cces$centered_age)
cor(cces$diff_cces_on_cces, cces$centered_age)

#################### Targets ###################
if(POPW) {
  cces_svy <- svydesign(ids = ~1, weights = ~commonweight_vv_post, data = cces)
} else {
  cces_svy <- svydesign(ids = ~1, data = cces)
}
#margin_sim = svymean(~diff_cces_on_cces, cces_svy)[1]* 100
margin_sim = svymean(~outcome, cces_svy)[1]* 100

targets_rake_demos_noeduc <- create_targets(cces_svy,
                                            formula_rake_demos_noeduc)
targets_rake_demos_weduc <- create_targets(cces_svy, formula_rake_demos_weduc)
targets_rake_all_vars <- create_targets(cces_svy, 
                                        formula_rake_all_vars)
targets_demo_truth <- create_targets(cces_svy, selection_model)


## Make table of Population Counts for post-stratification for manual ps function
cces_counts <- cces %>%
  group_by(strata) %>%
  summarize(n = if(!POPW) {n()} else {sum(commonweight_vv_post, na.rm = TRUE)}) %>%
  ungroup() %>%
  mutate(w = n / sum(n, na.rm = TRUE))

cces_reduc_counts <- cces %>%
  group_by(strata_reduc) %>%
  summarize(n = if(!POPW) {n()} else {sum(commonweight_vv_post, na.rm = TRUE)}) %>%
  ungroup() %>%
  mutate(w_reduc = n / sum(n, na.rm = TRUE)) 


cces_all_counts <- cces %>% 
  group_by(strata_all) %>% 
  summarize(n = if(!POPW) {n()} else {sum(commonweight_vv_post, na.rm = TRUE)}) %>%
  ungroup() %>%
  mutate(w_all = n / sum(n, na.rm = TRUE)) 


########################### RUN SIMS ##################################
#save mem
# rm(pew, lasso_pinclude, lasso_include,lasso_lambda,
#    stack_data, mod)

cces <- cces %>% mutate(recode_agesq = recode_age^2/ mean(recode_age^2),
                        recode_agecubed = recode_age^3/ mean(recode_age^3))

margins_formula <- ~recode_vote_2016 + 
  mod_cces_on_cces_pR + mod_cces_on_cces_pD + mod_cces_on_cces_pO+
  diff_cces_on_cces + margin_cces_on_cces + 
  recode_pid_3way + 
  recode_female + recode_race +recode_region + recode_educ + recode_relig_6way + 
  recode_born + recode_attndch_4way + recode_income_5way +
  recode_age_bucket + recode_age + recode_agesq + recode_agecubed + recode_age_factor + 
  recode_race_educ_reg + recode_educ_wh_3way + 
  recode_educ_pid_race +
  recode_pid_race + 
  recode_educ_pid +
  recode_race_reg_wh_educ + 
  recode_midwest_edu_race +
  recode_midwest_wh_edu

est_mean <- function(outcome, design) {
  svymean(as.formula(paste0("~", outcome)), design, na.rm = TRUE)[1]
}


#########################################
############## variance calc ###########
## Variance functions
var_fixed <- function(Y, weights, pop_size) {
    ## note: needs weights that sum to population total
    if(round(sum(weights)) != pop_size) { weights = weights*pop_size/sum(weights)}
    return(Hmisc::wtd.var(Y, weights))
}

## kott (14) (under poisson)
var_quasi <- function(weights, residuals, pop_size) {
    #moving from kott 14 w sum w =N to weights that sum to 1 + var of total to var of mean:
    #sum^n (w_i^2 - 1/N_pop w_i)e_i^2 
    return(sum((weights^2 - (weights / pop_size))*residuals^2))
}

## kott (15) linearization
var_linear <- function(weights, residuals, sample_size) {
    #moving from kott 14 w sum w =N to weights that sum to 1 + var of total to var of mean:
    # n/(n-1) sum^n (w_i*e_i)^2 - (1/n-1) [sum^n] *using this for now
    # approx = sum^n (w_i*e_i)^2 - (1/n) [sum^n]
    n = sample_size
    return((n/(n-1))*sum((weights * residuals)^2) - 1/n * sum(weights * residuals)^2)
}

## chad
var_chad <- function(weights, residuals) {
    return(sum(weights^2 * residuals^2))
}

## calculate all variances
calc_SEs <- function(Y, residuals, pop_size, weights, sample_size) {
    if(round(sum(weights)) != 1 ) {
        weights = weights/sum(weights)
    }
    return(data.frame(SE_fixed = sqrt(var_fixed(Y, weights, pop_size) / length(Y)),
                      SE_quasi = sqrt(var_quasi(weights, residuals, pop_size)),
                      SE_linear = sqrt(var_linear(weights, residuals, sample_size)),
                      SE_chad = sqrt(var_chad(weights, residuals))))
}

#rescale p_include to be diff n:
p_include = p_include*(n_sample/sum(p_include))
sum(p_include)

if(detectCores() < 20) {
    nsims = 500
    eval_kpop = FALSE
    SAVE = FALSE
}
coverage_eval = T
nsims = 1000
n_sample = 500
eval_kpop =F
#RAW ORIGINAL OUTCOME
cces$outcome = cces$diff_cces_on_cces
#adding noise but not rescaling
#cces$outcome = cces$diff_cces_on_cces + rnorm(nrow(cces),0, sd = sqrt(var(cces$diff_cces_on_cces)))
var(cces$outcome)
jank_rescale <- function(a, b, x) {
    min = min(x)
    max = max(x)
    ((x - min)/(max-min))*(b-a) + a
}
test = jank_rescale(-1, 1,x =  cces$outcome)
range(cces$outcome)
range(test)
var(test)
#cces$outcome = test
#var(cces$diff_cces_on_cces)
# cces$pD_upvar = plogis(cces$mod_cces_on_cces_pD + rnorm(nrow(cces), 0 , sqrt(var(cces$mod_cces_on_cces_pD))))
# cces$pR_upvar = plogis(cces$mod_cces_on_cces_pR + rnorm(nrow(cces), 0 , sqrt(var(cces$mod_cces_on_cces_pR))))
# cces$outcome = cces$pD_upvar - cces$pR_upvar
# cces$pD_upvar = (cces$mod_cces_on_cces_pD + rnorm(nrow(cces), 0 , sqrt(var(cces$mod_cces_on_cces_pD))))
# cces$pR_upvar = (cces$mod_cces_on_cces_pR + rnorm(nrow(cces), 0 , sqrt(var(cces$mod_cces_on_cces_pR))))
# cces$outcome = jank_rescale(-1,1,cces$pD_upvar - cces$pR_upvar)
var(cces$outcome)


system.time({
  sims <- mclapply(1:nsims, function(nsim) {
    #
    cat(paste("=====================  SIM:",nsim, 
              "===================== \n"))
      
    if(coverage_eval) {
        sample = sample.int(nrow(cces), size = n_sample, replace = F)
        survey_sim <- cces[sample, ]
    } else {
        sample <- rbinom(nrow(cces), 1, p_include)
        survey_sim <- cces[sample == 1, ]
    }
    
    survey_design <- suppressWarnings(svydesign(ids = ~1, data = survey_sim))
    
    ############################################
    ## Unweighted estimate
    ############################################
    
    unweighted <- est_mean("outcome", survey_design)
    ############################################
    ## Sample size
    ############################################
    n <- if(!coverage_eval) {sum(sample) } else { length(sample)}
    
    ############################################
    ## Raking on demographics (no education)
    ############################################
    rake_demos_noeduc_svyd <- calibrate(design = survey_design,
                                        formula = formula_rake_demos_noeduc,
                                        population = targets_rake_demos_noeduc,
                                        calfun = "raking")
    rake_demos_noeduc <- est_mean("outcome", rake_demos_noeduc_svyd)
    
    #SEs
    residuals = residuals(lm(update(formula_rake_demos_noeduc, outcome ~ .), 
                             data = rake_demos_noeduc_svyd$variables))
    rake_demos_noeduc_se <- calc_SEs(Y = rake_demos_noeduc_svyd$variables$outcome, 
                                     residuals = residuals, 
                                     pop_size = nrow(cces), 
                                     sample_size = sum(sample),
                                     weights = weights(rake_demos_noeduc_svyd))
    names(rake_demos_noeduc_se) = paste0("rake_demos_noeduc_", names(rake_demos_noeduc_se))
    
    if(coverage_eval) {
        rake_demos_noeduc_se_SVY = data.frame(rake_demos_noeduc_se_SVY = data.frame(svymean(~outcome, rake_demos_noeduc_svyd,
                                                                                            na.rm = TRUE))[1,2])
    }
    
    ############################################
    #### Raking on demographics (with education)
    ############################################
    
    rake_demos_weduc_svyd <- try(calibrate(design = survey_design,
                                       formula = formula_rake_demos_weduc,
                                       population = targets_rake_demos_weduc,
                                       calfun = "raking"), silent = T)

    rake_demos_weduc <- tryCatch(est_mean("outcome", rake_demos_weduc_svyd), 
                                error = function(e) NA)
    
    
    #SEs
    if(class(rake_demos_weduc_svyd)[1] == "try-error") {
        rake_demos_weduc_se <- data.frame(SE_fixed = NA,
                                          SE_quasi = NA, 
                                          SE_linear = NA, 
                                          SE_chad = NA) 
        names(rake_demos_weduc_se) = paste0("rake_demos_weduc_", names(rake_demos_weduc_se))
        
    } else {
        residuals = residuals(lm(update(formula_rake_demos_weduc, outcome ~ .), 
                                 data = rake_demos_weduc_svyd$variables))
        rake_demos_weduc_se <- calc_SEs(Y = rake_demos_weduc_svyd$variables$outcome, 
                                        residuals = residuals, 
                                        pop_size = nrow(cces),
                                        sample_size = sum(sample),
                                        weights = weights(rake_demos_weduc_svyd))
        names(rake_demos_weduc_se) = paste0("rake_demos_weduc_", names(rake_demos_weduc_se))
        
    }
    
    
    
    ############################################
    #### Raking on everything
    ############################################
    rake_all_svyd <- try(calibrate(design = survey_design,
                               formula = formula_rake_all_vars,
                               population = targets_rake_all_vars,
                               calfun = "raking"))
    
    rake_all <- tryCatch(est_mean("outcome", rake_all_svyd), 
                         error = function(e) NA)
    #SEs
    if(class(rake_all_svyd)[1] == "try-error") {
        rake_all_se <- data.frame(SE_fixed = NA,
                                          SE_quasi = NA, 
                                          SE_linear = NA, 
                                          SE_chad = NA) 
        names(rake_all_se) = paste0("rake_all_", names(rake_all_se))
        
    } else {
        residuals = residuals(lm(update(formula_rake_all_vars, outcome ~ .), 
                                 data = rake_all_svyd$variables))
        rake_all_se <- calc_SEs(Y = rake_all_svyd$variables$outcome, 
                                residuals = residuals, 
                                pop_size = nrow(cces), 
                                sample_size = sum(sample),
                                weights = weights(rake_all_svyd))
        names(rake_all_se) = paste0("rake_all_", names(rake_all_se))
        
    }
    
    
    ############################################
    ## Post-stratification: Old Formula
    ############################################
    post_stratification_svyd = svydesign(~1, data = postStrat(survey_sim, 
                                                              cces_counts, "w", 
                                                              strata_pass = "strata"),
                                         weights = ~w)
    
    post_stratification <- est_mean("outcome", post_stratification_svyd)
    
    #SEs
    post_stratification_se <- data.frame(post_strat_SE_svy = data.frame(svymean(~outcome, post_stratification_svyd, na.rm = TRUE))[1,2])
    
    ############################################
    ## Post-stratification: Reduced
    ############################################
    post_strat_reduc_svyd = svydesign(~1, data = postStrat(survey_sim, 
                                                           cces_reduc_counts, "w_reduc", 
                                                           strata_pass = "strata_reduc"),
                                      weights = ~w)
    
    post_strat_reduc <- est_mean("outcome", post_strat_reduc_svyd)
    
    #SEs
    post_strat_reduc_se <- data.frame(post_strat_reduc_SE_svy = data.frame(svymean(~outcome, post_strat_reduc_svyd, na.rm = TRUE))[1,2])
    
    
    
    ############################################
    ## Post-stratification: All
    ############################################
    post_strat_all_svyd = svydesign(~1, data = postStrat(survey_sim, 
                                                         cces_all_counts, "w_all", 
                                                         strata_pass = "strata_all"),
                                    weights = ~w)
    
    post_strat_all <- est_mean("outcome", post_strat_all_svyd)
    
    #SEs
    post_strat_all_se <- data.frame(post_strat_all_SE_svy = data.frame(svymean(~outcome, post_strat_all_svyd, na.rm = TRUE))[1,2])
    
    ############################################
    #### Raking on true model
    ############################################
    #very messy error catching for the moment just a stop gap to see how things look
    rake_truth_svyd <- try(calibrate(design = survey_design,
                                     formula = selection_model,
                                     population = targets_demo_truth,
                                     calfun = "raking",
                                     epsilon = .009), silent = T)

    
    rake_truth <- tryCatch(est_mean("outcome", rake_truth_svyd), 
                           error = function(e) NA)
    truth_margins <- tryCatch(svymean(margins_formula, rake_truth_svyd), 
                           error = function(e) NA)
    
    #SEs
    lambdas <- 10^seq(3, -2, by = -.1)
    x <- model.matrix(update(selection_model, outcome ~ .),
                      data = rake_truth_svyd$variables)
    fit <- glmnet(x, 
                  rake_truth_svyd$variables$outcome, alpha = 0, lambda = lambdas)
    cv_fit <- cv.glmnet(x, rake_truth_svyd$variables$outcome, alpha = 0, lambda = lambdas)
    opt_lambda <- cv_fit$lambda.min
    fit <- cv_fit$glmnet.fit
    
    residuals = rake_truth_svyd$variables$outcome - predict(fit, s = opt_lambda, newx = x)
    
    rake_truth_se <- tryCatch(calc_SEs(Y = rake_truth_svyd$variables$outcome,
                                       residuals = residuals,
                                       pop_size = nrow(cces),
                                       sample_size = sum(sample),
                                       weights = weights(rake_truth_svyd)), error = function(e) NA)
    
    if(length(rake_truth_se) == 1) {
        rake_truth_se <- data.frame(SE_fixed = NA, 
                                    SE_quasi = NA, 
                                    SE_linear = NA, 
                                    SE_chad = NA)
    }
    names(rake_truth_se) = tryCatch(paste0("rake_truth_", names(rake_truth_se)), error = function(e) NA)
    #why not still use lm?? OH bc i guess the selection model still has unused coefs?
    #but then like ok fine, lets only fit with the nonzero coefs
    #ohhh ok yeah it's annoying bc they're not one hot encoded so you can't easily auto select them
    #so erin's way makes sense but for now let me just try and replicate.. no it's annoying ok whatever
    #wait but then why is it ridge? are the coefs the same? like no right?
    #waaaait shit but wait a second shouldnt we then only be raking on the non-zero coefs? but like you cant do that ok ok nvmd
    #im still confused...
    
    #HT and Hayek
    p_sample <- if(!coverage_eval) {
        as.matrix((p_include[sample==1]))
    } else {
        as.matrix(rep(length(sample)/nrow(cces), length(sample)))
    }
    
    #HT
    ht_truth = if(!coverage_eval) {
        sum((cces[sample ==1, "outcome"]/p_sample))/nrow(cces) 
        } else {  sum((cces[sample, "outcome"]/p_sample))/nrow(cces) }
    
    #Hayek
    hayek_truth = if(!coverage_eval) {
        sum((cces[sample ==1, "outcome"]/p_sample))/sum(1/p_sample)
    } else { sum((cces[sample, "outcome"]/p_sample))/sum(1/p_sample)}
  
    
    ############################################
    ## Kpop: Categorical Data + b = argmax V(K)
    ############################################
    if(eval_kpop) {

      # Select the covariates for use in Kbal: updated cat data no cont age
      #one-hot coded for cat kernel
      kbal_data <- bind_rows(survey_sim %>% dplyr::select(recode_age_bucket,
                                                          recode_female,
                                                          recode_race,
                                                          recode_region,
                                                          recode_pid_3way,
                                                          recode_educ,

                                                          recode_income_5way,
                                                          recode_relig_6way,
                                                          recode_born,
                                                          recode_attndch_4way),
                             cces %>% dplyr::select(recode_age_bucket,
                                                    recode_female,
                                                    recode_race,
                                                    recode_region,
                                                    recode_pid_3way,
                                                    recode_educ,

                                                    recode_income_5way,
                                                    recode_relig_6way,
                                                    recode_born,
                                                    recode_attndch_4way))
      kbal_data_sampled <- c(rep(1, nrow(survey_sim)), rep(0, nrow(cces)))
      ##### Demos Constraint
      rake_demos_constraint <- bind_rows(survey_sim %>% dplyr::select(recode_age_bucket,
                                                                      recode_female,
                                                                      recode_race,
                                                                      recode_region,
                                                                      recode_pid_3way),
                                         cces %>% dplyr::select(recode_age_bucket,
                                                                recode_female,
                                                                recode_race,
                                                                recode_region,
                                                                recode_pid_3way))%>%
        model.matrix(as.formula("~."), .)

      rake_demos_constraint <- rake_demos_constraint[,-1]
      rake_demos_constraint <- scale(rake_demos_constraint)


      rake_demos_wedu_constraint <- bind_rows(survey_sim %>% dplyr::select(recode_age_bucket,
                                                                           recode_female,
                                                                           recode_race,
                                                                           recode_region,
                                                                           recode_pid_3way,
                                                                           recode_educ),
                                              cces %>% dplyr::select(recode_age_bucket,
                                                                     recode_female,
                                                                     recode_race,
                                                                     recode_region,
                                                                     recode_pid_3way,
                                                                     recode_educ))%>%
        model.matrix(as.formula("~."), .)

      rake_demos_wedu_constraint <- rake_demos_wedu_constraint[,-1]
      rake_demos_wedu_constraint <- scale(rake_demos_wedu_constraint)


      rake_all_constraint <- bind_rows(survey_sim %>% dplyr::select(recode_age_bucket,
                                                                    recode_female,
                                                                    recode_race,
                                                                    recode_region,
                                                                    recode_pid_3way,
                                                                    recode_educ,

                                                                    recode_income_5way,
                                                                    recode_relig_6way,
                                                                    recode_born,
                                                                    recode_attndch_4way),
                                       cces %>% dplyr::select(recode_age_bucket,
                                                              recode_female,
                                                              recode_race,
                                                              recode_region,
                                                              recode_pid_3way,
                                                              recode_educ,
                                                              recode_income_5way,
                                                              recode_relig_6way,
                                                              recode_born,
                                                              recode_attndch_4way)) %>%
        model.matrix(as.formula("~."), .)

      rake_all_constraint <- rake_all_constraint[,-1]
      rake_all_constraint <- scale(rake_all_constraint)
      b_manual = NA
      
      kpop <- lapply(1:length(b_manual), function(i) {
        
        #cat("===================== Running Kbal ==================\n")
        
        
        #### DEFAULT ######
        cat(paste("b:", b_manual[i], "nsim:", nsim, "DEFAULT", "\n"))
        kbal_est <- kbal(allx=kbal_data,
                         sampled = kbal_data_sampled,
                         #b = b_manual[i],
                         cat_data = TRUE,
                         incrementby = increment,
                         meanfirst = FALSE,
                         ebal.tol = tolerance,
                         ebal.maxit = maxit,
                         minnumdims = min_num_dims,
                         maxnumdims = max_num_dims,
                         linkernel = if(TEST){TRUE} else{ FALSE},
                         sampledinpop = FALSE,
                         fullSVD = TRUE)
        
        kpop_svyd <- svydesign(~1, data = survey_sim,
                               weights = kbal_est$w[kbal_data_sampled ==1])
        
        kpop <- est_mean("outcome", kpop_svyd)
        b_kpop = kbal_est$b
        #save memory by saving only the svd to re use
        svdK = kbal_est$svdK 
        numdims = kbal_est$numdims
        biasbound_r = kbal_est$biasbound_ratio
        biasbound = kbal_est$biasbound_opt
        
        ##### Kpop SEs
        kpop <- tryCatch(est_mean("outcome", kpop_svyd), error = function(e) NA)
        
        lambdas <- 10^seq(3, -2, by = -.1)
        x <- model.matrix( ~ ., data = as.data.frame(kbal_est$svdK$v[, 1:kbal_est$numdims]))
        fit <- glmnet(x, 
                      kpop_svyd$variables$outcome, alpha = 0, lambda = lambdas)
        cv_fit <- cv.glmnet(x, kpop_svyd$variables$outcome, alpha = 0, lambda = lambdas)
        opt_lambda <- cv_fit$lambda.min
        fit <- cv_fit$glmnet.fit
        
        residuals = kpop_svyd$variables$outcome - predict(fit, s = opt_lambda, newx = x)
        
        # lm_full = lm(kpop_svyd$variables$outcome ~ kbal_est$svdK$v[, 1:kbal_est$numdims])
        # full_res = residuals(lm_full)
        # cbind(residuals, full_res)
        kpop_se <- tryCatch(calc_SEs(Y = kpop_svyd$variables$outcome,
                                     residuals = residuals,
                                     pop_size = nrow(cces),
                                     sample_size = sum(sample),
                                     weights = weights(kpop_svyd)), error = function(e) NA)
       
        if(length(kpop_se) == 1) {
            kpop_se <- data.frame(SE_fixed = NA, 
                                  SE_quasi = NA, 
                                  SE_linear = NA, 
                                  SE_chad = NA)
        }
        names(kpop_se) = tryCatch(paste0("kpop_", names(kpop_se)), error = function(e) NA)
        
        #CONVERGED
        dist_record = data.frame(t(kbal_est$dist_record))
        min_converged = dist_record[which.min(dist_record[dist_record$Ebal.Convergence ==1,"BiasBound"]), "Dims"]
        
        rm(kbal_est)
        
        #### CONVG ####
        cat(paste("b:", b_manual[i], "nsim:", nsim, "CONV", "\n"))
        if(is.null(min_converged) | length(min_converged) ==0) {
          kpop_svyd_conv <- "dn converge"
          kpop_conv <- "dn converge"
          
          numdims_conv = "dn converge"
          biasbound_r_conv = "dn converge"
          biasbound_conv = "dn converge"
          kpop_conv_se = data.frame(SE_fixed = NA, 
                                    SE_quasi = NA, 
                                    SE_linear = NA, 
                                    SE_chad = NA)
          
        } else {
          kbal_est_conv <- kbal(allx=kbal_data,
                                K.svd = svdK,
                                sampled = kbal_data_sampled,
                                numdims = min_converged,
                                ebal.tol = tolerance,
                                ebal.maxit = maxit,
                                minnumdims = min_num_dims,
                                maxnumdims = max_num_dims,
                                scale_data = FALSE,
                                drop_MC = FALSE,
                                incrementby = increment,
                                meanfirst = FALSE,
                                sampledinpop = FALSE,
                                ebal.convergence = TRUE)
          kpop_svyd_conv <- svydesign(~1, data = survey_sim,
                                      weights = kbal_est_conv$w[kbal_data_sampled ==1])
          kpop_conv <- est_mean("outcome", kpop_svyd_conv)
          
          numdims_conv = kbal_est_conv$numdims
          biasbound_r_conv = kbal_est_conv$biasbound_ratio
          biasbound_conv = kbal_est_conv$biasbound_opt
          
          #SEs
          x <- model.matrix( ~ .,
                             data = as.data.frame(kbal_est_conv$svdK$v[, 1:kbal_est_conv$numdims]))
          fit <- glmnet(x, 
                        kpop_svyd_conv$variables$outcome, alpha = 0, lambda = lambdas)
          cv_fit <- cv.glmnet(x, kpop_svyd_conv$variables$outcome, alpha = 0, 
                              lambda = lambdas)
          opt_lambda <- cv_fit$lambda.min
          fit <- cv_fit$glmnet.fit
          residuals = kpop_svyd_conv$variables$outcome - predict(fit, s = opt_lambda, 
                                                                           newx = x)
          kpop_conv_se <- tryCatch(calc_SEs(Y = kpop_svyd_conv$variables$outcome,
                                       residuals = residuals,
                                       pop_size = nrow(cces),
                                       sample_size = sum(sample),
                                       weights = weights(kpop_svyd_conv)), error = function(e) NA)
          if(length(kpop_conv_se) == 1) {
              kpop_conv_se <- data.frame(SE_fixed = NA, 
                                    SE_quasi = NA, 
                                    SE_linear = NA, 
                                    SE_chad = NA)
          }
          names(kpop_conv_se) = tryCatch(paste0("kpop_conv_", names(kpop_se)), error = function(e) NA)
         
          rm(kbal_est_conv) 
        }
        
        
        ####### MF #######
        cat(paste("b:", b_manual[i], "nsim:", nsim, "MEANFIRST", "\n"))
        #### BROKEN HERE
        kbal_mf_est <- kbal(K.svd = svdK,
                            cat_data = T,
                            allx=kbal_data,
                            sampled = kbal_data_sampled,
                            ebal.tol = tolerance,
                            ebal.maxit = maxit,
                            minnumdims = min_num_dims,
                            maxnumdims = max_num_dims,
                            incrementby = increment,
                            meanfirst = TRUE,
                            sampledinpop = FALSE)
        
        kpop_mf_svyd <- svydesign(~1, data = survey_sim, 
                                  weights = kbal_mf_est$w[kbal_data_sampled ==1])
        
        kpop_mf <- est_mean("outcome", kpop_mf_svyd)
        
        mfnumdims = kbal_mf_est$numdims
        mf_appended_dims = kbal_mf_est$meanfirst_dims
        if(is.null(mf_appended_dims)) {mf_appended_dims = c("dn_converge")}
        biasbound_r_mf = kbal_mf_est$biasbound_ratio
        biasbound_mf = kbal_mf_est$biasbound_opt
        
        if(is.null(mfnumdims)) {
            mfnumdims = c("dn_converge") 
            kpop_mf_se = data.frame(SE_fixed = NA, 
                                    SE_quasi = NA, 
                                    SE_linear = NA, 
                                    SE_chad = NA)
        } else {
            
            #SEs X = V[, 1:numdims]
            V <- model.matrix( ~ .,
                               data = as.data.frame(kbal_mf_est$svdK$v[, 1:kbal_mf_est$numdims]))
            #binding mf cols for sample units to V
            X <- cbind(kbal_mf_est$appended_constraint_cols[kbal_data_sampled, ], V)
            
            fit <- glmnet(X, 
                          kpop_mf_svyd$variables$outcome, alpha = 0, lambda = lambdas)
            
            cv_fit <- cv.glmnet(X, kpop_mf_svyd$variables$outcome, alpha = 0, 
                                lambda = lambdas)
            opt_lambda <- cv_fit$lambda.min
            fit <- cv_fit$glmnet.fit
            residuals = kpop_mf_svyd$variables$outcome - predict(fit, s = opt_lambda, 
                                                                             newx = X)
            
            
            kpop_mf_se <- tryCatch(calc_SEs(Y = kpop_mf_svyd$variables$outcome,
                                              residuals = residuals,
                                              pop_size = nrow(cces),
                                            sample_size = sum(sample),
                                              weights = weights(kpop_mf_svyd)), 
                                   error = function(e) NA)
            if(length(kpop_mf_se) == 1) {
                kpop_mf_se <- data.frame(SE_fixed = NA, 
                                      SE_quasi = NA, 
                                      SE_linear = NA, 
                                      SE_chad = NA)
            }
            names(kpop_mf_se) = tryCatch(paste0("kpop_mf_", names(kpop_mf_se)),
                                         error = function(e) NA)
            
        }

        rm(kbal_mf_est)
        
        #########demos constraint method:
        cat(paste("b:", b_manual[i], "nsim:", nsim, "CONSTR", "\n"))
        kbal_demos_est <- kbal(K.svd = svdK,
                               allx=kbal_data,
                               #cat_data = TRUE,
                               sampled = kbal_data_sampled,
                               ebal.tol = tolerance,
                               ebal.maxit = maxit,
                               minnumdims = min_num_dims,
                               maxnumdims = max_num_dims,
                               scale_data = FALSE,
                               drop_MC = FALSE,
                               incrementby = increment,
                               #scaling these
                               constraint = rake_demos_constraint,
                               meanfirst = FALSE,
                               sampledinpop = FALSE)
        kpop_demos_svyd <- svydesign(~1, data = survey_sim, 
                                     weights = kbal_demos_est$w[kbal_data_sampled ==1])
        
        kpop_demos <- est_mean("outcome", kpop_demos_svyd)
        
        numdims_demos = kbal_demos_est$numdims
        if(is.null(numdims_demos)) {
            numdims_demos = c("dn_converge") 
            kpop_demos_se <- data.frame(SE_fixed = NA, 
                                     SE_quasi = NA, 
                                     SE_linear = NA, 
                                     SE_chad = NA)
        } else {
            #SEs X = V[, 1:numdims]
            V <- model.matrix( ~ .,
                               data = as.data.frame(kbal_demos_est$svdK$v[, 1:kbal_demos_est$numdims]))
            #binding constraint cols to
            X <- cbind(kbal_demos_est$appended_constraint_cols[kbal_data_sampled, ], V)
            
            fit <- glmnet(X, 
                          kpop_demos_svyd$variables$outcome, alpha = 0, lambda = lambdas)
            
            cv_fit <- cv.glmnet(X, kpop_demos_svyd$variables$outcome, alpha = 0, 
                                lambda = lambdas)
            opt_lambda <- cv_fit$lambda.min
            fit <- cv_fit$glmnet.fit
            residuals = kpop_demos_svyd$variables$outcome - predict(fit, s = opt_lambda, 
                                                                              newx = X)
            
            
            kpop_demos_se <- tryCatch(calc_SEs(Y = kpop_demos_svyd$variables$outcome,
                                               residuals = residuals,
                                               pop_size = nrow(cces),
                                               sample_size = sum(sample),
                                               weights = weights(kpop_demos_svyd)), 
                                      error = function(e) NA)
            if(length(kpop_demos_se) == 1) {
                kpop_demos_se <- data.frame(SE_fixed = NA, 
                                            SE_quasi = NA, 
                                            SE_linear = NA, 
                                            SE_chad = NA)
            }
            names(kpop_demos_se) = tryCatch(paste0("kpop_demos_", names(kpop_mf_se)),
                                            error = function(e) NA)
        }
        biasbound_r_demos = kbal_demos_est$biasbound_ratio
        biasbound_demos = kbal_demos_est$biasbound_opt
        
        rm(kbal_demos_est)
        
        
        #########demos + educ constraint method:
        cat(paste("b:", b_manual[i], "nsim:", nsim, "CONSTR", "\n"))
        kbal_demos_wedu_est <- kbal(K.svd = svdK,
                                    allx=kbal_data,
                                    cat_data = TRUE,
                                    sampled = kbal_data_sampled,
                                    ebal.tol = tolerance,
                                    ebal.maxit = maxit,
                                    minnumdims = min_num_dims,
                                    maxnumdims = max_num_dims,
                                    scale_data = FALSE,
                                    drop_MC = FALSE,
                                    incrementby = increment,
                                    #scaling these
                                    constraint = rake_demos_wedu_constraint,
                                    meanfirst = FALSE,
                                    sampledinpop = FALSE)
        kpop_demos_wedu_svyd <- svydesign(~1, data = survey_sim, 
                                          weights = kbal_demos_wedu_est$w[kbal_data_sampled ==1])
        
        kpop_demos_wedu <- est_mean("outcome", kpop_demos_wedu_svyd)
        
        numdims_demos_wedu = kbal_demos_wedu_est$numdims
        if(is.null(numdims_demos_wedu)) {
            numdims_demos_wedu = c("dn_converge")
            kpop_demos_wedu_se <- data.frame(SE_fixed = NA, 
                                        SE_quasi = NA, 
                                        SE_linear = NA, 
                                        SE_chad = NA)
        } else {
            #SEs X = V[, 1:numdims]
            V <- model.matrix( ~ .,
                               data = as.data.frame(kbal_demos_wedu_est$svdK$v[, 1:kbal_demos_wedu_est$numdims]))
            #binding constraint cols to
            X <- cbind(kbal_demos_wedu_est$appended_constraint_cols[kbal_data_sampled, ], V)
            
            fit <- glmnet(X, 
                          kpop_demos_wedu_svyd$variables$outcome, alpha = 0, lambda = lambdas)
            
            cv_fit <- cv.glmnet(X, kpop_demos_wedu_svyd$variables$outcome, alpha = 0, 
                                lambda = lambdas)
            opt_lambda <- cv_fit$lambda.min
            fit <- cv_fit$glmnet.fit
            residuals = kpop_demos_wedu_svyd$variables$outcome - predict(fit, s = opt_lambda, 
                                                                              newx = X)
            
            kpop_demos_wedu_se <- tryCatch(calc_SEs(Y = kpop_demos_wedu_svyd$variables$outcome,
                                               residuals = residuals,
                                               pop_size = nrow(cces),
                                               sample_size = sum(sample),
                                               weights = weights(kpop_demos_wedu_svyd)), 
                                      error = function(e) NA)
            if(length(kpop_demos_wedu_se) == 1) {
                kpop_demos_wedu_se <- data.frame(SE_fixed = NA, 
                                            SE_quasi = NA, 
                                            SE_linear = NA, 
                                            SE_chad = NA)
            }
            names(kpop_demos_wedu_se) = tryCatch(paste0("kpop_demos_wedu_", names(kpop_mf_se)),
                                            error = function(e) NA)
            
        }
        biasbound_r_demos_wedu = kbal_demos_wedu_est$biasbound_ratio
        biasbound_demos_wedu = kbal_demos_wedu_est$biasbound_opt
        
        rm(kbal_demos_wedu_est)
        
        
        #########all constraint method:
        cat(paste("b:", b_manual[i], "nsim:", nsim, "CONSTR", "\n"))
        kbal_all_est <- kbal(K.svd = svdK,
                             allx=kbal_data,
                             #cat_data = TRUE,
                             sampled = kbal_data_sampled,
                             ebal.tol = tolerance,
                             ebal.maxit = maxit,
                             minnumdims = min_num_dims,
                             maxnumdims = max_num_dims,
                             scale_data = FALSE,
                             drop_MC = FALSE,
                             incrementby = increment,
                             #scaling these
                             constraint = rake_all_constraint,
                             meanfirst = FALSE,
                             sampledinpop = FALSE)
        kpop_all_svyd <- svydesign(~1, data = survey_sim, 
                                   weights = kbal_all_est$w[kbal_data_sampled ==1])
        
        kpop_all <- est_mean("outcome", kpop_all_svyd)
        
        numdims_all = kbal_all_est$numdims
        if(is.null(numdims_all)) {
            numdims_all = c("dn_converge")
            numdims_all_se <- data.frame(SE_fixed = NA, 
                                             SE_quasi = NA, 
                                             SE_linear = NA, 
                                             SE_chad = NA)
        } else {
            #SEs X = V[, 1:numdims]
            V <- model.matrix( ~ .,
                               data = as.data.frame(kbal_all_est$svdK$v[, 1:kbal_all_est$numdims]))
            #binding constraint cols to
            X <- cbind(kbal_all_est$appended_constraint_cols[kbal_data_sampled, ], V)
            
            fit <- glmnet(X, 
                          kpop_all_svyd$variables$outcome, alpha = 0, lambda = lambdas)
            
            cv_fit <- cv.glmnet(X, kpop_all_svyd$variables$outcome, alpha = 0, 
                                lambda = lambdas)
            opt_lambda <- cv_fit$lambda.min
            fit <- cv_fit$glmnet.fit
            residuals = kpop_all_svyd$variables$outcome - predict(fit, s = opt_lambda, 
                                                                                   newx = X)
            
            kpop_all_se <- tryCatch(calc_SEs(Y = kpop_all_svyd$variables$outcome,
                                                    residuals = residuals,
                                                    pop_size = nrow(cces),
                                                    sample_size = sum(sample),
                                                    weights = weights(kpop_all_svyd)), 
                                           error = function(e) NA)
            if(length(kpop_demos_wedu_se) == 1) {
                kpop_demos_wedu_se <- data.frame(SE_fixed = NA, 
                                                 SE_quasi = NA, 
                                                 SE_linear = NA, 
                                                 SE_chad = NA)
            }
            names(kpop_all_se) = tryCatch(paste0("kpop_all_", names(kpop_mf_se)),
                                                 error = function(e) NA)
            
        }
        
        biasbound_r_all = kbal_all_est$biasbound_ratio
        biasbound_all = kbal_all_est$biasbound_opt
        
        rm(kbal_all_est)
        
        rm(svdK)
        
        ##### return
        out = list()
        b_out = b_kpop
        b = b_out
        out$sims = data.frame(b_out,
                              kpop,
                              kpop_mf,
                              kpop_conv,
                              kpop_demos,
                              kpop_demos_wedu,
                              kpop_all,
                              bb = biasbound,
                              bbr = biasbound_r,
                              bb_conv = biasbound_conv,
                              bbr_conv = biasbound_r_conv,
                              bb_mf = biasbound_mf,
                              bbr_mf = biasbound_r_mf,
                              bb_demos = biasbound_demos,
                              bbr_demos = biasbound_r_demos,
                              bb_demos_wedu = biasbound_demos_wedu,
                              bbr_demos_wedu = biasbound_r_demos_wedu,
                              bb_all = biasbound_all,
                              bbr_all = biasbound_r_all,
                              numdims,
                              numdims_conv,
                              mfnumdims, 
                              mf_appended_dims, 
                              numdims_demos,
                              numdims_demos_wedu,
                              numdims_all)
        
        #Starndard Errors:
        out$SEs = data.frame(rake_demos_noeduc_se,
                             rake_demos_weduc_se,
                             rake_all_se,
                             post_stratification_se,
                             post_strat_reduc_se,
                             post_strat_all_se,
                             rake_truth_se,
                             kpop_se,
                             kpop_conv_se,
                             kpop_mf_se,
                             kpop_demos_se,
                             kpop_demos_wedu_se,
                             kpop_all_se)
        
        
        
        #weights
        out$weights = list(b = b_out,
                           kpop_w = weights(kpop_svyd),
                           kpop_w_conv = weights(kpop_svyd_conv),
                           kpop_mf_w = weights(kpop_mf_svyd), 
                           kpop_demos_w = weights(kpop_demos_svyd),
                           kpop_demos_wedu_w = weights(kpop_demos_wedu_svyd),
                           kpop_all_w = weights(kpop_all_svyd))
        
        ######## Kpop Margins ########
        
        out$km <- round(cbind(b = b_out/100,
                              kpop = svymean(margins_formula, kpop_svyd),
                              kpop_conv = svymean(margins_formula, kpop_svyd_conv),
                              kpop_mf = svymean(margins_formula, kpop_mf_svyd),
                              kpop_demos = svymean(margins_formula, kpop_demos_svyd),
                              kpop_demos_wedu = svymean(margins_formula, kpop_demos_wedu_svyd),
                              kpop_all = svymean(margins_formula, kpop_all_svyd)) * 100,
                        4)
        
        rm(kpop_svyd, kpop_mf_svyd, kpop_svyd_conv, kpop_demos_svyd,
           kpop_demos_wedu_svyd, kpop_all_svyd)
        
        return(out)
      })
      
    }
    
    ############################################ OUTPUT
    out = list()
    if(eval_kpop) {
      out$sims = cbind(nsim,
                       n,
                       unweighted,
                       rake_demos_noeduc,
                       rake_demos_weduc,
                       rake_all,
                       post_stratification,
                       post_strat_reduc,
                       post_strat_all,
                       rake_truth,
                       ht_truth, 
                       hayek_truth,
                       lapply(kpop, `[[`,1) %>% bind_rows())
      
      out$SEs = lapply(kpop, `[[`,2) %>% bind_rows()
      out$weights = lapply(kpop, `[[`,3) %>% bind_cols()
      
      margin <- round(cbind(sample = svymean(margins_formula, survey_design),
                            cces =  svymean(margins_formula, cces_svy),
                            rake_demos_noeduc = svymean(margins_formula,
                                                        rake_demos_noeduc_svyd),
                            rake_demos_weduc = svymean(margins_formula,
                                                       rake_demos_weduc_svyd),
                            rake_all = svymean(margins_formula,
                                               rake_all_svyd),
                            post_stratification = svymean(margins_formula,
                                                          post_stratification_svyd),
                            post_strat_reduc = svymean(margins_formula,
                                                       post_strat_reduc_svyd),
                            post_strat_all = svymean(margins_formula,
                                                     post_strat_all_svyd),
                            rake_truth = truth_margins) * 100, 5)
      
      #these are just means so let's not multiply by 100
      margin["recode_age",] <- margin["recode_age",]/100 
      margin["recode_agesq",] <- margin["recode_agesq",]/100
      margin["recode_agecubed",] <- margin["recode_agecubed",]/100
      
      
      margin = cbind(margin,
                     as.data.frame(sapply(kpop, `[`, 4)))
      
      margin <- margin[,grepl("km.kpop", colnames(margin))|
                         !grepl("km.", colnames(margin)) ]
      
      # colnames(margin)[grepl("km.kpop", colnames(margin))] <- unlist(lapply(b_out, function(x) c(paste0("kpop_b", round(x,3)),
      #                                                                                        paste0("kpop_cvg_b", round(x,3)), 
      #                                                                                        paste0("kpop_mf_b", round(x,3)), 
      #                                                                                        paste0("kpop_demos_b", round(x,3)),
      #                                                                                        paste0("kpop_demos_wedu_b", round(x,3)) ,
      #                                                                                        paste0("kpop_all_b", round(x,3)) ) ))
    } else {
     ##### NO KPOP EVAL
    if(!coverage_eval) {
        margin <- round(cbind(sample = svymean(margins_formula, survey_design),
                              cces =  svymean(margins_formula, cces_svy),
                              rake_demos_noeduc = svymean(margins_formula,
                                                          rake_demos_noeduc_svyd),
                              rake_demos_weduc = svymean(margins_formula,
                                                         rake_demos_weduc_svyd),
                              rake_all = svymean(margins_formula,
                                                 rake_all_svyd),
                              post_stratification = svymean(margins_formula,
                                                            post_stratification_svyd),
                              post_strat_reduc = svymean(margins_formula,
                                                         post_strat_reduc_svyd),
                              post_strat_all = svymean(margins_formula,
                                                       post_strat_all_svyd),

                              rake_truth = truth_margins) * 100, 5)

        #these are just means so let's not multiply by 100
        margin["recode_age",] <- margin["recode_age",]/100
        margin["recode_agesq",] <- margin["recode_agesq",]/100
        margin["recode_agecubed",] <- margin["recode_agecubed",]/100
    } else {
        margin = NULL
    }
      
      
      out$sims = data.frame(nsim, n,
                            unweighted,
                            rake_demos_noeduc,
                            rake_demos_weduc,
                            rake_all,
                            post_stratification,
                            post_strat_reduc,
                            post_strat_all,
                            rake_truth,
                            ht_truth, 
                            hayek_truth)
      
      out$SEs = data.frame(rake_demos_noeduc_se,
                           rake_demos_weduc_se,
                           rake_all_se,
                           post_stratification_se,
                           post_strat_reduc_se,
                           post_strat_all_se,
                           rake_truth_se)
      if(coverage_eval) {
          out$SEs = data.frame(rake_demos_noeduc_se,
                     rake_demos_weduc_se,
                     rake_demos_noeduc_se_SVY,
                     rake_all_se,
                     post_stratification_se,
                     post_strat_reduc_se,
                     post_strat_all_se,
                     rake_truth_se)
      }

  
    } 
    
    out$margins = margin
    
    return(out)
    
  }, mc.cores = detectCores() - cores_saved) 
})

good = which(lapply(sims, function (x) return(class(x))) == "list")
length(good)

if(SAVE) {
    save(sims, tolerance, maxit, increment, min_num_dims,
     file = paste0("./cat_sims_lump_modeled_outcome_nodiag_",POPW, "_m",maxit, "_t",tolerance, "_inc",
                   increment, "mindims", min_num_dims,
                   Sys.Date(),
                   # str_sub(gsub("[[:space:]|[:punct:]]", "_",
                   #              gsub("[:alpha:]", "",
                   #                   Sys.time())),
                   #         start = 1, end = -3),
                   "_nsims", length(good),
                   ".RData"))
}
sims = sims[good]


#combines all weights across rows but can group by b to get them per iteration
if(eval_kpop) { 
    weights <- lapply(sims, `[[`, 3) %>% bind_rows() 
    margins <- lapply(sims, `[[`, 4) 
} else { 
    weights = NULL
    margins = lapply(sims, `[[`, 3) 
}

if(SAVE) {
    save(est, SEs, weights, margins, tolerance, maxit, increment, min_num_dims,
         file = paste0("./cat_sims_modeled_outcome_nodiag_",POPW, "_m",maxit, "_t",tolerance, "_inc",
                       increment, "mindims", min_num_dims,
                       Sys.Date(),
                       "_nsims", length(good),
                       ".RData")) 
}


##################### eval coverage ####################
colnames(SEs)
est
#bootstrap SE
apply(est[, -c(1,2)], 2, sd)
apply(est[, -c(1,2)], 2, function(x) quantile(x, probs = c(0.025, 0.975)))

SE_fixed = SEs[grepl("SE_fixed$", colnames(SEs))]
SE_linear = SEs[grepl("SE_linear$", colnames(SEs))]
SE_quasi = SEs[grepl("SE_quasi$", colnames(SEs))]
SE_chad= SEs[grepl("SE_chad$", colnames(SEs))]

SE_svy= SEs[grepl("SVY", colnames(SEs))]

coverage <- function(SE, x_bar, truth =NULL, crit_val= qnorm(0.975)) {
    if(is.null(truth)) {
        truth = svymean(~outcome, cces_svy)[1]
    }
    x_upper = x_bar +  (SE*crit_val)
    x_lower = x_bar - (SE*crit_val)
    contains_truth = matrix(NA, ncol = ncol(SE), nrow = 1)
    for(i in 1:ncol(x_upper)) {
        contains_truth[,i] = sum((truth <= x_upper[,i] & truth >= x_lower[,i]))/nrow(SE)
    }
    colnames(contains_truth) = colnames(x_bar)
    return(contains_truth)
}
est3 = est[grepl("rake", colnames(est))]
nrow(est3)
#est3 = na.omit(est3)
#nrow(est3)
#SE_chad = na.omit(SE_chad)
#unbiasedness?
mean(cces$outcome) -colMeans(est[, -c(1,2)])


# eval coverage of diff SEs
all_SE_coverage <- function(sims, drop_NA = F, truth = NULL) {
    est <- lapply(sims, `[[`, 1) %>% bind_rows()
    SEs <- lapply(sims, `[[`, 2) %>% bind_rows()
    if(drop_NA) {
        est_temp = na.omit(est)
        n_drop = nrow(est) - nrow(est_temp)
        est = est_temp
        SEs = na.omit(SEs)
    }
    
    est_c = est[grepl("rake", colnames(est))]
    SE_fixed = SEs[grepl("SE_fixed$", colnames(SEs))]
    SE_linear = SEs[grepl("SE_linear$", colnames(SEs))]
    SE_quasi = SEs[grepl("SE_quasi$", colnames(SEs))]
    SE_chad= SEs[grepl("SE_chad$", colnames(SEs))]
    
    SE_svy= SEs[grepl("SVY", colnames(SEs))]
    if(ncol(SE_svy) != 0){
        #just making sure we're getting the estimates for the same SEs that we output from svy obj which currently is demos_noedu
        search = gsub("_se_SVY","", colnames(SE_svy))
        grepl(search, colnames(est_c))
        s = coverage(SE_svy, est_c[,grepl(search, colnames(est_c))])
        s1 = rep(NA, ncol(SE_fixed))
        s1[grepl(search, colnames(est_c))] = s
        #colnames(s1) = colnames(SE_fixed)
    }
    coverage_out = rbind(coverage(SE_fixed, est_c, truth = truth),
                coverage(SE_linear, est_c, truth = truth),
                coverage(SE_quasi, est_c,truth = truth),
                coverage(SE_chad,est_c,truth = truth), 
                s1)
    rownames(coverage_out) = c("SE_fixed", "SE_linear", "SE_quasi", "SE_chad", "SE_svy")
    if(drop_NA) {
        out = list()
        out$n_drop = n_drop
        out$coverage = coverage_out
    } else {
        out = coverage_out
    }
    return(out)
}

#sims_n500_diff = sims
#sims_n500_diff_upvar = sims
all_SE_coverage(sims_n500_diff, truth = mean(cces$diff_cces_on_cces))
#upvar = cces$outcome
all_SE_coverage(sims_n500_diff_upvar, drop_NA = T, truth = mean(upvar))
#upvar_scaled = cces$outcome
sims_n500_diff_upvar_scaled = sims
all_SE_coverage(sims_n500_diff_upvar_scaled,drop_NA = T, truth = mean(upvar_scaled))
#plogis PD_upvar - plog pR_upvar
#plog  = cces$outcome
#sims_n500_plog = sims
all_SE_coverage(sims_n500_plog,drop_NA = T, truth = mean(plog))
man_pD_pR = cces$outcome
sims_n500_man_pD_pR = sims
all_SE_coverage(sims_n500_man_pD_pR,drop_NA = T, truth = mean(man_pD_pR))


#lets check the bias
est <- lapply(sims, `[[`, 1) %>% bind_rows()
colMeans(est[, -c(1,2)] - mean(cces$outcome))

colMeans(SE_fixed)
colMeans(SE_linear)
colMeans(SE_svy)
#on par with bootstrap, a little bigger in bootstrap 0.008 vs 0.007
apply(est[, -c(1,2, 7,8,9)], 2, sd)


###############


########## RES ################
plot = est
nrow(est)
#estimates
colMeans(est)

####### Complex DGPP + no POPW + n= 500




################## results plots
#### ported over from tex
if(eval_kpop) {
    plot_lasso_margin <- plot %>% 
        dplyr::select(unweighted, 
                      rake_demos_noeduc,
                      rake_demos_weduc,
                      rake_all,
                      #post_stratification,
                      post_strat_reduc,
                      #post_strat_all,
                      rake_truth,
                      kpop, 
                      #kpop_conv,
                      #kpop_mf, 
                      kpop_demos,
                      kpop_demos_wedu,
                      kpop_all, 
                      ht_truth,
                      hayek_truth) %>% 
        pivot_longer(everything(),
                     names_to = "estimator", 
                     values_to = "margin") %>%
        mutate(margin = margin * 100,
               estimator_name = factor(case_when(estimator == "kpop" ~ "kpop",
                                                 estimator == "kpop_mf" ~ "kpop aMF (All)",
                                                 # estimator == "kpop_conv" ~ "kpop Converged",
                                                 estimator == "kpop_demos" ~ "kpop+MF:\n (Demos)",
                                                 estimator == "kpop_demos_wedu" ~ "kpop+MF:\n (Demos+Edu)",
                                                 estimator == "kpop_all" ~ "kpop+MF:\n (All)",
                                                 estimator == "rake_demos_noeduc" ~ "Mean Calibration:\n (Demos)",
                                                 estimator == "rake_demos_weduc" ~  "Mean Calibration:\n (Demos+Edu)",
                                                 estimator == "rake_all" ~ "Mean Calibration:\n (All)",
                                                 estimator == "rake_truth" ~ "Mean Calibration:\n True Selection\nModel",
                                                 estimator == "post_stratification" ~ "Post-Strat Prev",
                                                 estimator == "post_strat_reduc" ~ "Post-Stratification:\n (Reduc)",
                                                 estimator == "post_strat_all" ~ "Post-Strat All",
                                                 estimator == "unweighted" ~ "Unweighted",
                                                 estimator == "ht_truth" ~ "Horvitz-Thompson",
                                                 estimator == "hayek_truth" ~ "Hayek"),
                                       levels = c("Unweighted", 
                                                  "Mean Calibration:\n (Demos)",
                                                  "Mean Calibration:\n (Demos+Edu)",
                                                  "Mean Calibration:\n (All)",
                                                  #"Post-Strat Prev", 
                                                  "Post-Stratification:\n (Reduc)", 
                                                  "Post-Strat All",
                                                  "kpop",
                                                  # "kpop Converged",
                                                  #"kpop aMF (All)",
                                                  "kpop+MF:\n (Demos)",
                                                  "kpop+MF:\n (Demos+Edu)",
                                                  "kpop+MF:\n (All)",
                                                  "Mean Calibration:\n True Selection\nModel",
                                                  "Horvitz-Thompson",
                                                  "Hayek"
                                       ) ) )
    
} else {
    plot_lasso_margin <- plot %>% 
        dplyr::select(unweighted, 
                      rake_demos_noeduc,
                      rake_demos_weduc,
                      rake_all,
                      post_stratification,
                      post_strat_reduc,
                      post_strat_all,
                      rake_truth, 
                      ht_truth, 
                      hayek_truth) %>% 
        pivot_longer(everything(),
                     names_to = "estimator", 
                     values_to = "margin") %>%
        mutate(margin = margin * 100,
               estimator_name = factor(case_when(estimator == "kpop" ~ "kpop",
                                                 estimator == "kpop_mf" ~ "kpop aMF (All)",
                                                 # estimator == "kpop_conv" ~ "kpop Converged",
                                                 estimator == "kpop_demos" ~ "kpop+MF:\n (Demos)",
                                                 estimator == "kpop_demos_wedu" ~ "kpop+MF:\n (Demos+Edu)",
                                                 estimator == "kpop_all" ~ "kpop+MF:\n (All)",
                                                 estimator == "rake_demos_noeduc" ~ "Mean Calibration:\n (Demos)",
                                                 estimator == "rake_demos_weduc" ~  "Mean Calibration:\n (Demos+Edu)",
                                                 estimator == "rake_all" ~ "Mean Calibration:\n (All)",
                                                 estimator == "rake_truth" ~ "Mean Calibration:\n True Selection\nModel",
                                                 estimator == "post_stratification" ~ "Post-Strat Prev",
                                                 estimator == "post_strat_reduc" ~ "Post-Stratification:\n (Reduc)",
                                                 estimator == "post_strat_all" ~ "Post-Strat All",
                                                 estimator == "unweighted" ~ "Unweighted",
                                                 
                                                 estimator == "ht_truth" ~ "Horvitz-Thompson",
                                                 estimator == "hayek_truth" ~ "Hayek"),
                                       levels = c("Unweighted", 
                                                  "Mean Calibration:\n (Demos)",
                                                  "Mean Calibration:\n (Demos+Edu)",
                                                  "Mean Calibration:\n (All)",
                                                   "Post-Strat Prev", 
                                                  "Post-Stratification:\n (Reduc)", 
                                                  "Post-Strat All",
                                                  "kpop",
                                                  #"kpop Converged",
                                                  #"kpop aMF (All)",
                                                  "kpop+MF:\n (Demos)",
                                                  "kpop+MF:\n (Demos+Edu)",
                                                  "kpop+MF:\n (All)",
                                                  "Mean Calibration:\n True Selection\nModel",
                                                  "Horvitz-Thompson",
                                                  "Hayek"
                                       )))
    
}

#target:
#margin_sim = svymean(~outcome, cces_svy)[1]* 100
#### Box Plot
ggplot(data = plot_lasso_margin,
       aes(x = estimator_name, y = margin)) +
    geom_boxplot(alpha = 0.2) +
    geom_hline(yintercept = margin_sim) +
    theme_bw() +
    # ggtitle(paste0("Simulation Results: lasso ",
    #                " ", "pop weights: ", POPW, "\n",length(good), " sims, b=maxvar")) +
    xlab("") +
    ylab("Modeled Vote Margin") +
    annotate(geom = "text", x = 0.85, y = margin_sim+0.25, size = 2.7, angle = 90,
             label = "True Target\nPopulation\nMargin", hjust = 0) +
    ggtitle(paste0(nrow(est)," sims w/avg n_samp =", round(mean(est$n)), " simple DGP ", simple_selection_model)) +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))


### table
table = plot_lasso_margin %>% 
    mutate(estimator_name = gsub("\n", " ", estimator_name)) %>%
    group_by(estimator_name) %>%
    summarize(
        Bias = mean(margin - margin_sim),
        SE = sd(margin),
        MSE = mean((margin - margin_sim)^2)
    ) %>%
    mutate(
        Bias_Reduc = 1- Bias / Bias[estimator_name == "Unweighted"]
    ) %>%
    arrange(MSE)
table
### kable
knitr::kable(
    plot_lasso_margin %>% 
        mutate(estimator_name = gsub("\n", " ", estimator_name)) %>%
        group_by(estimator_name) %>%
        summarize(
            Bias = mean(margin - margin_sim),
            SE = sd(margin),
            MSE = mean((margin - margin_sim)^2)
        ) %>%
        mutate(
            Bias_Reduc = 1- Bias / Bias[estimator_name == "Unweighted"]
        ) %>%
        arrange(MSE), caption = "Simulation Results (numeric)",
    booktabs = T,
    format = "latex", digits = 2)
