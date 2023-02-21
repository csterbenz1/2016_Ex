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
coverage_eval = TRUE
noise = 1
TEST = FALSE # to run with a linear kernel so it's way faster; UPDATE: errors catch this as mistake and prevent
tolerance = 1e-4
maxit = 500
#both for runtime
increment = 5
min_num_dims = NULL
max_num_dims = NULL
SAVE = TRUE
##### Central Params to adjust
n_sample = 500
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

if(coverage_eval) {
    formula_ps <- ~recode_age_3way + recode_female + recode_pid_3way
}

formula_ps_reduc <- ~recode_age_3way + recode_female +
  recode_race + recode_region + recode_pid_3way  +
  recode_income_3way + recode_born + recode_educ_wh_3way

# #let's coarsen income and religion and attend church for all:
formula_ps_all <- ~recode_age_3way + recode_female +
  recode_race + recode_region + recode_pid_3way + recode_educ_3way +
  recode_income_3way + recode_relig_6way + recode_born + recode_attndch_bin

create_targets <- function (target_design, target_formula) {
  target_mf <- model.frame(target_formula, model.frame(target_design))
  target_mm <- model.matrix(target_formula, target_mf)
  wts <- weights(target_design)

  return(colSums(target_mm * wts) / sum(wts))
}

manual_rescale <- function(dat, a =0, b= 1) {
    rescaled = (b-a)*( dat - min(dat))/(max(dat) - min(dat)) + a
    return(rescaled)
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


### Load Target Data
cces <- readRDS(paste0(path_data, "cces_lasso_061021.rds"))


######### Make STRATA variable in CCES and Pew ############
cces <- bind_cols(cces, cces %>%
                    unite("strata", all.vars(formula_ps), remove = FALSE) %>%
                    unite("strata_reduc", all.vars(formula_ps_reduc),
                          remove = FALSE) %>%
                    unite("strata_all", all.vars(formula_ps_all)) %>%
                    dplyr::select(strata, strata_reduc, strata_all))
#weirdly this is producing still bias w ps on the correct formula apparently maybe because we are dropping the one category  with unite so im going to try this manual way:
cces = cces %>% mutate(strata = paste(recode_pid_3way,recode_female, recode_age_bucket, sep = "_"))

pew <- bind_cols(pew, pew %>%
                   unite("strata", all.vars(formula_ps), remove = FALSE) %>%
                   unite("strata_reduc", all.vars(formula_ps_reduc),
                         remove = FALSE) %>%
                   unite("strata_all", all.vars(formula_ps_all)) %>%
                   dplyr::select(strata, strata_reduc, strata_all))

pew = pew %>% mutate(strata = paste(recode_pid_3way, recode_female, recode_age_bucket, sep = "_"))


 ##################### LASSO: Selection #############################


#income, religion, pid x race
if(!coverage_eval) {

    if(simple_selection_model) {
        #first attempt to make this worse
        selection_model = as.formula(~recode_pid_3way:poly(recode_age, 2) +
                                         recode_female:recode_pid_3way)
        #first attempt to make this worse
        #center age then square
        pew = pew %>% mutate(centered_age = scale(recode_age, scale = F))
        cces = cces %>% mutate(centered_age = scale(recode_age, scale = F))

        selection_model = as.formula(~recode_pid_3way*poly(centered_age, 2) +
                                         #recode_age:recode_pid_3way +
                                         recode_female*recode_pid_3way)

        #selection_model = as.formula(~recode_pid_3way + recode_age_bucket + recode_female)

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

    # #imbalance:
    # mod_p <- mod[stack_data$S ==1,]
    # mod_c <- mod[stack_data$S ==0,]
    # nrow(mod_p)
    # nrow(mod_c)
    # imb = cbind(pew = colMeans(mod_p), cces=colMeans(mod_c))
    # imb = cbind(imb, diff = imb[,2]- imb[,1])
    # (round(imb, 3))
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
    # df = as.data.frame(as.matrix(lasso_include_coefs))
    #
    # kable(round(df,3), format = "latex", booktabs = T)
    # #probability of being in pew for cces
    lasso_pinclude = predict(lasso_include,
                             s= lasso_lambda$lambda.min,
                             type = "response",
                             newx = mod[stack_data$S == 0,-1])

    p_include <- lasso_pinclude
    sum(p_include)
    #check manual: yep
    # colnames(mod[stack_data$S == 0,])
    # xb = mod[stack_data$S == 0,] %*% lasso_include_coefs
    # xb[1:10]
    # dim(xb)
    # p_inc_man = plogis(xb[,1])
    # cbind(p_include, p_inc_man)[1:10,]

    cor(p_include, cces$outcome)
    cor(p_include, cces$mod_cces_on_cces_pD)
    cor(p_include, cces$mod_cces_on_cces_pR)

    ########### Sampling Inclusion Results
    summary(p_include)
    sum(p_include)
    mean(p_include)
    p_include_raw = p_include


    intercept_shift = 0
    p_include = p_include*(n_sample/sum(p_include)) + intercept_shift
    mean(p_include)
    summary(p_include)
    sum(p_include)

    #################### Define Outcome #########
    cces$outcome = cces$diff_cces_on_cces
} else {


    ########## DESIGN SELECTION MODEL: Specifying Coefs Directly
    #coefs: pid, age, gender
    selection_model = as.formula(~recode_pid_3way + recode_age_bucket + recode_female)
    cces_expanded = model.matrix(selection_model, data = cces)
    #needs to be n x p X p x 1 -> coef matrix is p x 1
    coefs = matrix(NA,nrow = ncol(cces_expanded), ncol =1 )
    rownames(coefs) = colnames(cces_expanded)
    coefs[,1] = c(-5, #intercept
                  .1, #selection of indep pos
                  .4, #selection of R pos
                  -.2, #36-50,
                  .3, #51-64,
                  .8, #65+,
                  .4 #male pos
                  )
    xbeta = cces_expanded %*% coefs
    p_include = plogis(xbeta)
    sum(p_include)
    summary(p_include)
    #two options to deal with sample size from 3k or so to 500: 1) rescale manaully 2) adjust coefs
    #since we can directly adjust the ceofs i will just do that (making the intercept quite neg)
    
    #out = data.frame(raw = c(summary(p_include), sum(p_include)) )
    
    #rownames(out) = c("Min", "25%", "Median", "Mean", "75%", "Max", "Sum")
    
    # kable(out, format = "latex", booktabs = T, 
    #       caption = "Sample Inclusion Probabilities")

    #################### DESIGN OUTCOME MODEL ##################
    #let's first try modeling pR and pD and subtracting and then doing bernoulli; same vars just diff coefs

    coefs_outcome = coefs
    #pR
    # coefs_outcome[,1] = c(-1, #intercept
    #                       -.1, #  indep pos
    #                       3, #  R pos
    #                       -.1, #36-50,
    #                       -.3, #51-64,
    #                       1, #65+,xbeta_outcome
    #                       .6#male pos
    # )
    # xbeta_outcome = cces_expanded %*% coefs_outcome
    # summary(xbeta_outcome)
    # p_R = xbeta_outcome
    # #p_R = manual_rescale(xbeta_outcome)
    # mean(p_R)
    # #sum(p_R)/nrow(cces)
    # #pD
    coefs_outcome[,1] = c(1, #intercept
                          .5, #  indep pos
                           -3, #  R pos
                          .4,#36-50,
                          .5, #51-64,
                          -.6, #65+,
                          -.3 #male pos
    )
    xbeta_outcome = cces_expanded %*% coefs_outcome
    # p_D = xbeta_outcome
    # mean(p_D)
    # cbind(p_R, p_D)[1:15,]
    # mean(p_D)
    # vote_diff = p_D - p_R
    # mean(vote_diff)
    # summary(vote_diff)
    # vote_diff = manual_rescale(p_D) - manual_rescale(p_R)
    # summary(vote_diff)
    # #vote_diff = manual_rescale(vote_diff)
   
    
    # .14:
    #     coefs_outcome[,1] = c(0, #intercept
    #                           2.5, #  indep pos
    #                           -3, #  R pos
    #                           .2,#36-50,
    #                           .5, #51-64,
    #                           -.8, #65+,
    #                           -.6 #male pos
    #     )
    
    #13.7 and wild:
    # c(0, #intercept
    #   2.5, #  indep pos
    #   -3, #  R pos
    #   .2,#36-50,
    #   .5, #51-64,
    #   -.8, #65+,
    #   -1.1 #male pos
    # )
    #.11
    # coefs_outcome[,1] = c(0, #intercept
    #                       2.5, #  indep pos
    #                       -3, #  R pos
    #                       .2,#36-50,
    #                       .5, #51-64,
    #                       -.3, #65+,
    #                       -1.1 #male pos
    # )
    #mean 0.06!
    # coefs_outcome[,1] = c(0, #intercept
    #                       2.5, #  indep pos
    #                       -3, #  R pos
    #                       -.5,#36-50,
    #                       -.5, #51-64,
    #                       .4, #65+,
    #                       -1.1 #male pos
    # )
    #coefs_outcome = coefs_outcome/2
    xbeta_outcome = cces_expanded %*% coefs_outcome
    vote_diff = xbeta_outcome
    #rescaleL
    #vote_diff = manual_rescale(vote_diff, -1,1)
    mean(vote_diff)
    #plot(density(vote_diff))
    summary(vote_diff)
    
    #add a little noise to lower the R^2 from literally 1 and SE estimators kind of break or get coverage of 1
    #let's not add any noise yet bc that makes the whole dist normal for the r^2 to be reasonable
    cat(paste("Adding sd(outcome)*",noise, "\n"))
    vote_diff = vote_diff + rnorm(nrow(cces), mean = 0, sd = sd(vote_diff)*noise)
    vote_diff = manual_rescale(vote_diff, -1,1)
    plot(density(vote_diff))
    summary(vote_diff)
    s = summary(lm(vote_diff ~ recode_pid_3way + recode_age_bucket + recode_female,data = cces))
    s$adj.r.squared
    #mean(vote_diff)
    #the mean is still a little high... how is the mean set
    
    # sum(p_D + p_R >1 )/nrow(cces)
    # vote_diff = p_D - p_R
    # mean(vote_diff)
    # s = summary(lm(vote_diff ~ recode_pid_3way + recode_age_bucket + recode_female,data = cces))
    # s$adj.r.squared
    # #oof let's try and add some noise
    # vote_diff = p_D - p_R
    # mean(vote_diff)
    # vote_diff = vote_diff + rnorm(nrow(cces), mean = 0, sd = sd(vote_diff))
    # mean(vote_diff)
    #svymean(~recode_vote_2016, cces_svy)[1] - svymean(~recode_vote_2016, cces_svy)[3]
    # test = vote_diff + rnorm(nrow(cces), 1, sd = sd(vote_diff))
    # var(test)
    # var(vote_diff)
    s = summary(lm(vote_diff ~ recode_pid_3way + recode_age_bucket + recode_female,data = cces))
    cat(paste("R^2 outcome is", round(s$adj.r.squared,3), "\n"))
    cat(paste("Mean outcome (target) is", round(mean(vote_diff)*100,3)))
    cat(paste("\nCorr of sampling prob and outcome ", round(cor(vote_diff, p_include),3)))
    cces$outcome = vote_diff

}


#################### Targets ###################
if(POPW) {
  cces_svy <- svydesign(ids = ~1, weights = ~commonweight_vv_post, data = cces)
} else {
  cces_svy <- svydesign(ids = ~1, data = cces)
}
margin_sim = svymean(~outcome, cces_svy)[1]* 100
margin_sim
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
#rm(pew, lasso_pinclude, lasso_include,lasso_lambda,
#   stack_data, mod)

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


if(detectCores() < 20) {
    nsims = 1000
    eval_kpop = F
    SAVE = FALSE
}
nsims
sum(p_include)
SAVE
coverage_eval
eval_kpop
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
    
    unweighted <- est_mean("outcome", survey_design)
    ############################################
    ## Sample size
    ############################################\
    n <- sum(sample) 
    
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
                                     maxit = 100,
                                     calfun = "raking"), silent = T)
    
    if(class(rake_truth_svyd)[1] == "try-error") {
        rake_truth_svyd <- try(calibrate(design = survey_design,
                                         formula = selection_model,
                                         population = targets_demo_truth,
                                         calfun = "raking",
                                         epsilon = .009), silent = T)
    }
    
    rake_truth <- tryCatch(est_mean("outcome", rake_truth_svyd), 
                           error = function(e) NA)
    truth_margins <- tryCatch(svymean(margins_formula, rake_truth_svyd), 
                           error = function(e) NA)
    
    #SEs
    lambdas <- 10^seq(3, -2, by = -.1)
    x <- model.matrix(update(selection_model, outcome ~ .),
                      data = rake_truth_svyd$variables)[, -1]
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
    
    #HT and Hayek
    p_sample <- as.matrix((p_include[sample==1]))
    
    #HT
    ht_truth = sum((cces[sample ==1, "outcome"]/p_sample))/nrow(cces) 
    #Hayek
    hayek_truth = sum((cces[sample ==1, "outcome"]/p_sample))/sum(1/p_sample)
    
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
        x <- as.matrix(data.frame(kbal_dims = kbal_est$svdK$v[, 1:kbal_est$numdims]))
        cv_fit <- cv.glmnet(x, kpop_svyd$variables$outcome, alpha = 0, lambda = lambdas)
        residuals = kpop_svyd$variables$outcome - predict(cv_fit$glmnet.fit, s = cv_fit$lambda.min, newx = x)
        
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
          x <- as.matrix(data.frame(kbal_dims = kbal_est_conv$svdK$v[, 1:kbal_est_conv$numdims]))
          cv_fit <- cv.glmnet(x, kpop_svyd_conv$variables$outcome, alpha = 0, 
                              lambda = lambdas)
          fit <- cv_fit$glmnet.fit
          residuals = kpop_svyd_conv$variables$outcome - predict(cv_fit$glmnet.fit, s = cv_fit$lambda.min, 
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
          names(kpop_conv_se) = tryCatch(paste0("kpop_conv_", names(kpop_conv_se)), error = function(e) NA)
         
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
            #after much hair puling i have determined the following
            #1. it is best to allow glmnet to add the intercept (it's not receptive to turning it off if you have it already bc you use model matrix; also bc the penalty.factor will always regularize the intercept when it's intercept when its internally aded)
            #2. glmnet is entirely redundant when you run cv.glmnet, but it is very important to predict with cv_glmnet$glmnet.fit + s = lambda.min otherwise you get different answers, BUT no errors (eg predict on the raw cv.glmnet obj)
            #3. penalty.factor will rescale the lambdas internally except when the user specifies the lambda sequence directly
            #4. you do not need to anticipate the addition of the intercept in penalty.factor when it's added internally
            #5. FOR GODS SAKE dont be an idiot and subset with only a vector name (kbal_data_sampled), make it a logical ==1!!!!! 
            #SEs X = V[, 1:numdims]
            V <-  data.frame(kbal_dims = kbal_mf_est$svdK$v[, c(1:kbal_mf_est$numdims)])
            #binding mf cols for sample units to V
            X <- as.matrix(cbind(kbal_mf_est$appended_constraint_cols[kbal_data_sampled==1, ], V))

            cv_fit <- cv.glmnet(X, kpop_mf_svyd$variables$outcome, alpha = 0, lambda = lambdas,
                                penalty.factor = c(rep(0, kbal_mf_est$meanfirst_dims), rep(1, kbal_mf_est$numdims)))

            residuals = kpop_mf_svyd$variables$outcome - predict(cv_fit$glmnet.fit, s = cv_fit$lambda.min, 
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
            V <-  data.frame(kbal_dims = kbal_demos_est$svdK$v[, c(1:kbal_demos_est$numdims)])
            X <- as.matrix(cbind(kbal_demos_est$appended_constraint_cols[kbal_data_sampled==1, ], V))
            
            cv_fit <- cv.glmnet(X, kpop_demos_svyd$variables$outcome, alpha = 0, lambda = lambdas,
                                penalty.factor = c(rep(0, ncol(kbal_demos_est$appended_constraint_cols)), rep(1, kbal_demos_est$numdims)))
            
            residuals =  kpop_demos_svyd$variables$outcome - predict(cv_fit$glmnet.fit, s = cv_fit$lambda.min, 
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
            names(kpop_demos_se) = tryCatch(paste0("kpop_demos_", names(kpop_demos_se)),
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
            V <-  data.frame(kbal_dims = kbal_demos_wedu_est$svdK$v[, c(1:kbal_demos_wedu_est$numdims)])
            X <- as.matrix(cbind(kbal_demos_wedu_est$appended_constraint_cols[kbal_data_sampled==1, ], V))
            
            cv_fit <- cv.glmnet(X, kpop_demos_wedu_svyd$variables$outcome, alpha = 0, lambda = lambdas,
                                penalty.factor = c(rep(0, ncol(kbal_demos_wedu_est$appended_constraint_cols)),
                                                   rep(1, kbal_demos_wedu_est$numdims)))
            
            residuals =  kpop_demos_wedu_svyd$variables$outcome - predict(cv_fit$glmnet.fit, s = cv_fit$lambda.min, 
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
            names(kpop_demos_wedu_se) = tryCatch(paste0("kpop_demos_wedu_", names(kpop_demos_wedu_se)),
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
            V <-  data.frame(kbal_dims = kbal_all_est$svdK$v[, c(1:kbal_all_est$numdims)])
            X <- as.matrix(cbind(kbal_all_est$appended_constraint_cols[kbal_data_sampled==1, ], V))
            
            cv_fit <- cv.glmnet(X, kpop_all_svyd$variables$outcome, alpha = 0, lambda = lambdas,
                                penalty.factor = c(rep(0, ncol(kbal_all_est$appended_constraint_cols)),
                                                   rep(1, kbal_all_est$numdims)))
            
            residuals =  kpop_all_svyd$variables$outcome - predict(cv_fit$glmnet.fit, s = cv_fit$lambda.min, 
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
            names(kpop_all_se) = tryCatch(paste0("kpop_all_", names(kpop_all_se)),
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
        if(coverage_eval) {
            out$SEs = data.frame(rake_demos_noeduc_se,
                                 rake_demos_noeduc_se_SVY,
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
        } else {
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
        }
        
        
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

outcome = cces$outcome
good = which(lapply(sims, function (x) return(class(x))) == "list")
length(good)

if(SAVE) {
    save(sims, outcome,  tolerance, maxit, increment, min_num_dims,
         file = paste0("./coverage_kpop", eval_kpop,
                       Sys.Date(),
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
    if(!coverage_eval) {margins = lapply(sims, `[[`, 3) }
}


##################### eval coverage ####################


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

# eval coverage of diff SEs
all_SE_coverage <- function(sims, drop_NA = F, truth = NULL, methods = c("rake|kpop")) {
    est <- lapply(sims, `[[`, 1) %>% bind_rows()
    SEs <- lapply(sims, `[[`, 2) %>% bind_rows()
    if(drop_NA) {
        est_temp = na.omit(est)
        n_drop = nrow(est) - nrow(est_temp)
        est = est_temp
        SEs = na.omit(SEs)
    }
    
    est_c = est[grepl(methods, colnames(est))]
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
        coverage_out = rbind(coverage(SE_fixed, est_c, truth = truth),
                             coverage(SE_linear, est_c, truth = truth),
                             coverage(SE_quasi, est_c,truth = truth),
                             coverage(SE_chad,est_c,truth = truth), 
                             s1)
        rownames(coverage_out) = c("SE_fixed", "SE_linear", "SE_quasi", "SE_chad", "SE_svy")
    } else {
        coverage_out = rbind(coverage(SE_fixed, est_c, truth = truth),
                             coverage(SE_linear, est_c, truth = truth),
                             coverage(SE_quasi, est_c,truth = truth),
                             coverage(SE_chad,est_c,truth = truth))
        rownames(coverage_out) = c("SE_fixed", "SE_linear", "SE_quasi", "SE_chad")
    }
    
    
    if(drop_NA) {
        out = list()
        out$n_drop = n_drop
        out$coverage = coverage_out
    } else {
        out = coverage_out
    }
    return(out)
}


#res: 
good = which(lapply(sims, function (x) return(class(x))) == "list")
length(good)

est <- lapply(sims, `[[`, 1) %>% bind_rows()
cols = if(eval_kpop) { c(3:7,10:12,14:19)} else {c(3:7,10:12)}
bias = colMeans(est[, cols])
bias = bias - mean(cces$outcome)
bias = data.frame(bias = t(t(bias))*100)
bias = bias %>% arrange(desc(abs(bias)))
round(bias,3)
# kable(round(bias, 3), format = "latex", booktabs = T, 
#       caption = paste0("Bias \\textbf{in Percent} across ", length(good), " sims: All Methods (Target = ", round(mean(outcome),3)*100, ")"))

SE_coverage = all_SE_coverage(sims, truth = mean(outcome))
SE_coverage

SEs = lapply(sims, `[[`, 2) %>% bind_rows()
colMeans(SEs)*100

# kable(round(SE_coverage[, c(1:4)], 3), format = "latex", booktabs = T, 
#       caption = paste("SE Coverage Results: Raking", length(good), "sims"))
# kable(round(SE_coverage[, c(5:ncol(SE_coverage))], 3), format = "latex", booktabs = T, 
#       caption = paste("SE Coverage Results: Kpop", length(good), "sims"))


if(SAVE) {
    save(sims, outcome, SE_coverage, bias,  tolerance, maxit, increment, min_num_dims,
         file = paste0("./coverage_res_kpop", eval_kpop,
                       Sys.Date(),
                       "_nsims", length(good),
                       ".RData"))
}

########## RES ################
plot = est


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
gg_out = ggplot(data = plot_lasso_margin,
                aes(x = estimator_name, y = margin)) +
    geom_boxplot(alpha = 0.2) +
    geom_hline(yintercept = margin_sim) +
    theme_bw() +
    xlab("") +
    ylab("Modeled Vote Margin") +
    annotate(geom = "text", x = 0.85, y = margin_sim+0.25, size = 2.7, angle = 90,
             label = "True Target\nPopulation\nMargin", hjust = 0) +
    ggtitle(paste0(nrow(est)," sims w/avg n_samp =", round(mean(est$n)), " simple DGP ", simple_selection_model)) +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))

gg_out
### table
table = plot_lasso_margin %>% 
    mutate(estimator_name = gsub("\n", " ", estimator_name)) %>%
    group_by(estimator_name) %>%
    summarize(
        Bias = mean(margin - margin_sim),
        SE_boot= sd(margin),
        MSE = mean((margin - margin_sim)^2)
    ) %>%
    mutate(
        Bias_Reduc = 1- Bias / Bias[estimator_name == "Unweighted"]
    ) %>%
    arrange(MSE)
table
colMeans(SEs)*100


if(SAVE) {
    save(sims, outcome, SE_coverage, bias,table, gg_out,
         tolerance, maxit, increment, min_num_dims,
         file = paste0("./coverage_res_all_kpop", eval_kpop,
                       Sys.Date(),
                       "_nsims", length(good),
                       ".RData"))
}
