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
  cores_saved = 14
} else {
  path_data= "/Users/Ciara_1/Dropbox/kpop/Updated/application/data/" 
  cores_saved = 6
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

##### Central Params to adjust
n_sample = 500
simple_selection_model = TRUE
nsims = (detectCores()-cores_saved)*1
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
samp_units = sample(1:nrow(pew), n_sample, replace = F)
pew = pew[samp_units, ]


### Load Target Data
cces <- readRDS(paste0(path_data, "cces_lasso_061021.rds"))


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
# Stack data with S = 1 indicating Pew
stack_data <- data.frame(bind_rows(pew, cces), 
                         S = c(rep(1, nrow(pew)), rep(0, nrow(cces))))


#income, religion, pid x race
if(simple_selection_model) {
    selection_model = as.formula(~recode_pid_3way:poly(recode_age, 2) + 
                                     recode_female:recode_pid_3way)
} else {
    selection_model = as.formula(~recode_female:recode_pid_3way + 
                                     recode_age:recode_pid_3way +
                                     #adds a bit mroe bias to edu+d
                                     recode_pid_race + 
                                     recode_race_reg_wh_educ +
                                     recode_educ_wh_3way +
                                     poly(recode_age, 3))
}



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


#probability of being in pew for cces 
lasso_pinclude = predict(lasso_include,
                         s= lasso_lambda$lambda.min,
                         type = "response",
                         newx = mod[stack_data$S == 0,-1])

coefs <- as.matrix(lasso_include_coefs)
#coefs

p_include <- lasso_pinclude
sum(p_include)


########### Sampling Inclusion Results
summary(p_include)

cor(cces$mod_cces_on_cces_pD, p_include)
cor(cces$mod_cces_on_cces_pR, p_include)
cor(cces$diff_cces_on_cces, p_include)


#################### Targets ###################
if(POPW) {
  cces_svy <- svydesign(ids = ~1, weights = ~commonweight_vv_post, data = cces)
} else {
  cces_svy <- svydesign(ids = ~1, data = cces)
}
margin_sim = svymean(~diff_cces_on_cces, cces_svy)[1]* 100

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
    #if(sum(weights) != pop_size) { weights = weights*pop_size/sum(weights)}
    return(Hmisc::wtd.var(Y, weights * pop_size))
}

## kott (14) (under poisson)
var_quasi <- function(weights, residuals, pop_size) {
    return(sum((weights^2 - weights / pop_size) * residuals^2))
}

## kott (15) linearization
var_linear <- function(weights, residuals, pop_size) {
    return(sum((weights * residuals)^2) - 1/pop_size * sum(weights * residuals)^2)
}

## chad
var_chad <- function(weights, residuals) {
    return(sum(weights^2 * residuals^2))
}

## calculate all variances
calc_SEs <- function(Y, residuals, pop_size, weights) {
    return(data.frame(SE_fixed = sqrt(var_fixed(Y, weights, pop_size) / length(Y)),
                      SE_quasi = sqrt(var_quasi(weights, residuals, pop_size)),
                      SE_linear = sqrt(var_linear(weights, residuals, pop_size)),
                      SE_chad = sqrt(var_chad(weights, residuals))))
}


system.time({
  sims <- mclapply(1:nsims, function(nsim) {
    #
    cat(paste("=====================  SIM:",nsim, 
              "===================== \n"))
    sample <- rbinom(nrow(cces), 1, p_include)
    
    survey_sim <- cces[sample == 1, ]
    survey_design <- suppressWarnings(svydesign(ids = ~1, data = survey_sim))
    
    ############################################
    ## Unweighted estimate
    ############################################
    
    unweighted <- est_mean("diff_cces_on_cces", survey_design)
    ############################################
    ## Sample size
    ############################################
    n <- sum(sample)
    
    ############################################
    ## Raking on demographics (no education)
    ############################################
    rake_demos_noeduc_svyd <- calibrate(design = survey_design,
                                        formula = formula_rake_demos_noeduc,
                                        population = targets_rake_demos_noeduc,
                                        calfun = "raking")
    rake_demos_noeduc <- est_mean("diff_cces_on_cces", rake_demos_noeduc_svyd)
    
    #SEs
    residuals = residuals(lm(update(formula_rake_demos_noeduc, diff_cces_on_cces ~ .), 
                             data = rake_demos_noeduc_svyd$variables))
    rake_demos_noeduc_se <- calc_SEs(Y = rake_demos_noeduc_svyd$variables$diff_cces_on_cces, 
                                     residuals = residuals, 
                                     pop_size = nrow(cces), 
                                     weights = weights(rake_demos_noeduc_svyd))
    names(rake_demos_noeduc_se) = paste0("rake_demos_noeduc_", names(rake_demos_noeduc_se))
    
    
    ############################################
    #### Raking on demographics (with education)
    ############################################
    rake_demos_weduc_svyd <- calibrate(design = survey_design,
                                       formula = formula_rake_demos_weduc,
                                       population = targets_rake_demos_weduc,
                                       calfun = "raking")
    
    rake_demos_weduc <- est_mean("diff_cces_on_cces", rake_demos_weduc_svyd)
    
    
    #SEs
    residuals = residuals(lm(update(formula_rake_demos_weduc, diff_cces_on_cces ~ .), 
                             data = rake_demos_weduc_svyd$variables))
    rake_demos_weduc_se <- calc_SEs(Y = rake_demos_weduc_svyd$variables$diff_cces_on_cces, 
                                    residuals = residuals, 
                                    pop_size = nrow(cces), 
                                    weights = weights(rake_demos_weduc_svyd))
    names(rake_demos_weduc_se) = paste0("rake_demos_weduc_", names(rake_demos_weduc_se))
    
    
    ############################################
    #### Raking on everything
    ############################################
    rake_all_svyd <- calibrate(design = survey_design,
                               formula = formula_rake_all_vars,
                               population = targets_rake_all_vars,
                               calfun = "raking")
    
    rake_all <- est_mean("diff_cces_on_cces", rake_all_svyd)
    
    #SEs
    residuals = residuals(lm(update(formula_rake_all_vars, diff_cces_on_cces ~ .), 
                             data = rake_all_svyd$variables))
    rake_all_se <- calc_SEs(Y = rake_all_svyd$variables$diff_cces_on_cces, 
                                    residuals = residuals, 
                                    pop_size = nrow(cces), 
                                    weights = weights(rake_all_svyd))
    names(rake_all_se) = paste0("rake_all_", names(rake_all_se))
    
    
    
    ############################################
    ## Post-stratification: Old Formula
    ############################################
    post_stratification_svyd = svydesign(~1, data = postStrat(survey_sim, 
                                                              cces_counts, "w", 
                                                              strata_pass = "strata"),
                                         weights = ~w)
    
    post_stratification <- est_mean("diff_cces_on_cces", post_stratification_svyd)
    
    #SEs
    residuals = residuals(lm(diff_cces_on_cces ~ strata,
                             data = post_stratification_svyd$variables))
    post_stratification_se <- calc_SEs(Y = post_stratification_svyd$variables$diff_cces_on_cces, 
                                       residuals = residuals, 
                                       pop_size = nrow(cces), 
                                       weights = weights(post_stratification_svyd))
    names(post_stratification_se) = paste0("post_strat_", names(post_stratification_se))
    
    ############################################
    ## Post-stratification: Reduced
    ############################################
    post_strat_reduc_svyd = svydesign(~1, data = postStrat(survey_sim, 
                                                           cces_reduc_counts, "w_reduc", 
                                                           strata_pass = "strata_reduc"),
                                      weights = ~w)
    
    post_strat_reduc <- est_mean("diff_cces_on_cces", post_strat_reduc_svyd)
    
    #SEs
    residuals = residuals(lm(diff_cces_on_cces ~ strata_reduc,
                             data = post_strat_reduc_svyd$variables))
    post_strat_reduc_se <- calc_SEs(Y = post_strat_reduc_svyd$variables$diff_cces_on_cces, 
                                       residuals = residuals, 
                                       pop_size = nrow(cces), 
                                       weights = weights(post_strat_reduc_svyd))
    names(post_strat_reduc_se) = paste0("post_strat_reduc_", names(post_strat_reduc_se))
    
    ############################################
    ## Post-stratification: All
    ############################################
    post_strat_all_svyd = svydesign(~1, data = postStrat(survey_sim, 
                                                         cces_all_counts, "w_all", 
                                                         strata_pass = "strata_all"),
                                    weights = ~w)
    
    post_strat_all <- est_mean("diff_cces_on_cces", post_strat_all_svyd)
    
    #SEs
    residuals = residuals(lm(diff_cces_on_cces ~ strata_all,
                             data = post_strat_all_svyd$variables))
    post_strat_all_se <- calc_SEs(Y = post_strat_all_svyd$variables$diff_cces_on_cces, 
                                       residuals = residuals, 
                                       pop_size = nrow(cces), 
                                       weights = weights(post_strat_all_svyd))
    names(post_strat_all_se) = paste0("post_strat_all_", names(post_strat_all_se))
    
    ############################################
    #### Raking on true model
    ############################################
    #very messy error catching for the moment just a stop gap to see how things look
    rake_truth_svyd <- try(calibrate(design = survey_design,
                                     formula = selection_model,
                                     population = targets_demo_truth,
                                     calfun = "raking",
                                     epsilon = .009), silent = T)

    
    rake_truth <- tryCatch(est_mean("diff_cces_on_cces", rake_truth_svyd), 
                           error = function(e) NA)
    truth_margins <- tryCatch(svymean(margins_formula, rake_truth_svyd), 
                           error = function(e) NA)
    
    #SEs
    lambdas <- 10^seq(3, -2, by = -.1)
    x <- model.matrix(update(selection_model, diff_cces_on_cces ~ .),
                      data = rake_truth_svyd$variables)
    fit <- glmnet(x, 
                  rake_truth_svyd$variables$diff_cces_on_cces, alpha = 0, lambda = lambdas)
    cv_fit <- cv.glmnet(x, rake_truth_svyd$variables$diff_cces_on_cces, alpha = 0, lambda = lambdas)
    opt_lambda <- cv_fit$lambda.min
    fit <- cv_fit$glmnet.fit
    
    residuals = rake_truth_svyd$variables$diff_cces_on_cces - predict(fit, s = opt_lambda, newx = x)
    
    rake_truth_se <- tryCatch(calc_SEs(Y = rake_truth_svyd$variables$diff_cces_on_cces,
                                       residuals = residuals,
                                       pop_size = nrow(cces),
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
    p_sample <- as.matrix((p_include[sample==1]))
    
    #HT
    ht_truth = sum((cces[sample ==1, "diff_cces_on_cces"]/p_sample))/nrow(cces)
    
    #Hayek
    hayek_truth = sum((cces[sample ==1, "diff_cces_on_cces"]/p_sample))/sum(1/p_sample)
  
    
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
        
        kpop <- est_mean("diff_cces_on_cces", kpop_svyd)
        b_kpop = kbal_est$b
        #save memory by saving only the svd to re use
        svdK = kbal_est$svdK 
        numdims = kbal_est$numdims
        biasbound_r = kbal_est$biasbound_ratio
        biasbound = kbal_est$biasbound_opt
        
        ##### Kpop SEs
        kpop <- tryCatch(est_mean("diff_cces_on_cces", kpop_svyd), error = function(e) NA)
        
        lambdas <- 10^seq(3, -2, by = -.1)
        x <- model.matrix( ~ ., data = as.data.frame(kbal_est$svdK$v[, 1:kbal_est$numdims]))
        fit <- glmnet(x, 
                      kpop_svyd$variables$diff_cces_on_cces, alpha = 0, lambda = lambdas)
        cv_fit <- cv.glmnet(x, kpop_svyd$variables$diff_cces_on_cces, alpha = 0, lambda = lambdas)
        opt_lambda <- cv_fit$lambda.min
        fit <- cv_fit$glmnet.fit
        
        residuals = kpop_svyd$variables$diff_cces_on_cces - predict(fit, s = opt_lambda, newx = x)
        
        # lm_full = lm(kpop_svyd$variables$diff_cces_on_cces ~ kbal_est$svdK$v[, 1:kbal_est$numdims])
        # full_res = residuals(lm_full)
        # cbind(residuals, full_res)
        kpop_se <- tryCatch(calc_SEs(Y = kpop_svyd$variables$diff_cces_on_cces,
                                     residuals = residuals,
                                     pop_size = nrow(cces),
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
          kpop_conv <- est_mean("diff_cces_on_cces", kpop_svyd_conv)
          
          numdims_conv = kbal_est_conv$numdims
          biasbound_r_conv = kbal_est_conv$biasbound_ratio
          biasbound_conv = kbal_est_conv$biasbound_opt
          
          #SEs
          x <- model.matrix( ~ .,
                             data = as.data.frame(kbal_est_conv$svdK$v[, 1:kbal_est_conv$numdims]))
          fit <- glmnet(x, 
                        kpop_svyd_conv$variables$diff_cces_on_cces, alpha = 0, lambda = lambdas)
          cv_fit <- cv.glmnet(x, kpop_svyd_conv$variables$diff_cces_on_cces, alpha = 0, 
                              lambda = lambdas)
          opt_lambda <- cv_fit$lambda.min
          fit <- cv_fit$glmnet.fit
          residuals = kpop_svyd_conv$variables$diff_cces_on_cces - predict(fit, s = opt_lambda, 
                                                                           newx = x)
          kpop_conv_se <- tryCatch(calc_SEs(Y = kpop_svyd_conv$variables$diff_cces_on_cces,
                                       residuals = residuals,
                                       pop_size = nrow(cces),
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
        
        kpop_mf <- est_mean("diff_cces_on_cces", kpop_mf_svyd)
        
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
                          kpop_mf_svyd$variables$diff_cces_on_cces, alpha = 0, lambda = lambdas)
            
            cv_fit <- cv.glmnet(X, kpop_mf_svyd$variables$diff_cces_on_cces, alpha = 0, 
                                lambda = lambdas)
            opt_lambda <- cv_fit$lambda.min
            fit <- cv_fit$glmnet.fit
            residuals = kpop_mf_svyd$variables$diff_cces_on_cces - predict(fit, s = opt_lambda, 
                                                                             newx = X)
            
            
            kpop_mf_se <- tryCatch(calc_SEs(Y = kpop_mf_svyd$variables$diff_cces_on_cces,
                                              residuals = residuals,
                                              pop_size = nrow(cces),
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
        
        kpop_demos <- est_mean("diff_cces_on_cces", kpop_demos_svyd)
        
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
                          kpop_demos_svyd$variables$diff_cces_on_cces, alpha = 0, lambda = lambdas)
            
            cv_fit <- cv.glmnet(X, kpop_demos_svyd$variables$diff_cces_on_cces, alpha = 0, 
                                lambda = lambdas)
            opt_lambda <- cv_fit$lambda.min
            fit <- cv_fit$glmnet.fit
            residuals = kpop_demos_svyd$variables$diff_cces_on_cces - predict(fit, s = opt_lambda, 
                                                                              newx = X)
            
            
            kpop_demos_se <- tryCatch(calc_SEs(Y = kpop_demos_svyd$variables$diff_cces_on_cces,
                                               residuals = residuals,
                                               pop_size = nrow(cces),
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
        
        kpop_demos_wedu <- est_mean("diff_cces_on_cces", kpop_demos_wedu_svyd)
        
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
                          kpop_demos_wedu_svyd$variables$diff_cces_on_cces, alpha = 0, lambda = lambdas)
            
            cv_fit <- cv.glmnet(X, kpop_demos_wedu_svyd$variables$diff_cces_on_cces, alpha = 0, 
                                lambda = lambdas)
            opt_lambda <- cv_fit$lambda.min
            fit <- cv_fit$glmnet.fit
            residuals = kpop_demos_wedu_svyd$variables$diff_cces_on_cces - predict(fit, s = opt_lambda, 
                                                                              newx = X)
            
            kpop_demos_wedu_se <- tryCatch(calc_SEs(Y = kpop_demos_wedu_svyd$variables$diff_cces_on_cces,
                                               residuals = residuals,
                                               pop_size = nrow(cces),
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
        
        kpop_all <- est_mean("diff_cces_on_cces", kpop_all_svyd)
        
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
                          kpop_all_svyd$variables$diff_cces_on_cces, alpha = 0, lambda = lambdas)
            
            cv_fit <- cv.glmnet(X, kpop_all_svyd$variables$diff_cces_on_cces, alpha = 0, 
                                lambda = lambdas)
            opt_lambda <- cv_fit$lambda.min
            fit <- cv_fit$glmnet.fit
            residuals = kpop_all_svyd$variables$diff_cces_on_cces - predict(fit, s = opt_lambda, 
                                                                                   newx = X)
            
            kpop_all_se <- tryCatch(calc_SEs(Y = kpop_all_svyd$variables$diff_cces_on_cces,
                                                    residuals = residuals,
                                                    pop_size = nrow(cces),
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
      
    } else { # XXXX THIS HAS NOT BEEN UPDATED AND WILL BREAK IF eval_kpop = F
      kpop <- NA
      kpop_mf<- NA
      kpop_conv<- NA
      kpop_demos <- NA
      kpop_demos_wedu <- NA
      kpop_all <- NA
      
      numdims = NA
      mfnumdims = NA
      mf_appended_dims = NA
      numdims_conv = NA
      numdims_demos = NA
      numdims_demos_wedu = NA
      numdims_all = NA
      biasbound_r = NA
      biasbound = NA
      bb_conv = NA
      bbr_conv = NA
      bb_mf = NA
      bbr_mf = NA
      bb_demos = NA
      bbr_demos = NA
      bb_demos_wedu = NA
      bbr_demos_wedu = NA
      bb_all = NA
      bbr_all = NA
      
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
      
      colnames(margin)[grepl("km.kpop", colnames(margin))] <- unlist(lapply(b, function(x) c(paste0("kpop_b", round(x,3)),
                                                                                             paste0("kpop_cvg_b", round(x,3)), 
                                                                                             paste0("kpop_mf_b", round(x,3)), 
                                                                                             paste0("kpop_demos_b", round(x,3)),
                                                                                             paste0("kpop_demos_wedu_b", round(x,3)) ,
                                                                                             paste0("kpop_all_b", round(x,3)) ) ))
    } else {
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
    }
    
    out$margins = margin
    
    return(out)
    
  }, mc.cores = detectCores() - cores_saved) 
})

good = which(lapply(sims, function (x) return(class(x))) == "list")
length(good)

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

est <- lapply(sims, `[[`, 1) %>% bind_rows()
SEs <- lapply(sims, `[[`, 2) %>% bind_rows()

#combines all weights across rows but can group by b to get them per iteration
weights <- lapply(sims, `[[`, 3) %>% bind_rows()
margins <- lapply(sims, `[[`, 4) 

save(est, SEs, weights, margins, tolerance, maxit, increment, min_num_dims,
     file = paste0("./cat_sims_modeled_outcome_nodiag_",POPW, "_m",maxit, "_t",tolerance, "_inc",
                   increment, "mindims", min_num_dims,
                   Sys.Date(),
                   "_nsims", length(good),
                   ".RData"))
