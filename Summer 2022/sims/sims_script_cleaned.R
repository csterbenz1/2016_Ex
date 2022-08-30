### Packages
library(MASS)
library(tidyverse)
library(survey)
library(kbal)
library(parallel)
library(knitr)
library(glmnet)
library(tictoc)

###### SET PARAMS  ###############
set.seed(9345876)

if(detectCores() > 10) {
  path_data = "/home/csterbenz/Data/"
} else {
  path_data= "/Users/Ciara/Dropbox/kpop/Updated/application/data/" 
}

POPW = FALSE
cores_saved = 14
nsims = (detectCores()-cores_saved)*15
eval_kpop = TRUE
TEST = FALSE # to run with a linear kernel so it's way faster
tolerance = 1e-4
maxit = 500
#both for runtime
increment = 5
min_num_dims = 150
max_num_dims = 500 

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
selection_model = as.formula(~recode_female:recode_pid_3way + recode_age:recode_pid_3way + recode_race:recode_region:recode_educ_wh_3way:recode_pid_3way + recode_age + I(recode_age^2))

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
selection_model = as.formula(~recode_female:recode_pid_3way + 
                               recode_age:recode_pid_3way +
                               #adds a bit mroe bias to edu+d
                               recode_pid_race + 
                               recode_race_reg_wh_educ +
                               recode_educ_wh_3way +
                               poly(recode_age, 3))


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
rm(pew, lasso_pinclude, lasso_include,lasso_lambda,
   stack_data, mod)

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

#function to get argmax V(k) b value
var_K= function(b, n_d, diag_length){
  d <- n_d[,1] %>% pull()
  n_d <- as.vector(n_d[,2] %>% pull())
  #REMOVING DIAGONAL 0 COUNTS FROM MAXIMIZATION CONSIDERATION
  n_d[1] <- n_d[1] - diag_length
  p_d <- n_d/ sum(n_d) 
  
  mean_k = sum(exp(-1*d/b)*p_d)
  var_k = sum((exp(-1*d/b)-mean_k)^2 * p_d)
  return(var_k)
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
    ############################################
    #### Raking on demographics (with education)
    ############################################
    rake_demos_weduc_svyd <- calibrate(design = survey_design,
                                       formula = formula_rake_demos_weduc,
                                       population = targets_rake_demos_weduc,
                                       calfun = "raking")
    
    rake_demos_weduc <- est_mean("diff_cces_on_cces", rake_demos_weduc_svyd)
    
    
    ############################################
    #### Raking on everything
    ############################################
    rake_all_svyd <- calibrate(design = survey_design,
                               formula = formula_rake_all_vars,
                               population = targets_rake_all_vars,
                               calfun = "raking")
    
    rake_all <- est_mean("diff_cces_on_cces", rake_all_svyd)
    
    
    ############################################
    ## Post-stratification: Old Formula
    ############################################
    post_stratification_svyd = svydesign(~1, data = postStrat(survey_sim, 
                                                              cces_counts, "w", 
                                                              strata_pass = "strata"),
                                         weights = ~w)
    
    post_stratification <- est_mean("diff_cces_on_cces", post_stratification_svyd)
    
    ############################################
    ## Post-stratification: Reduced
    ############################################
    post_strat_reduc_svyd = svydesign(~1, data = postStrat(survey_sim, 
                                                           cces_reduc_counts, "w_reduc", 
                                                           strata_pass = "strata_reduc"),
                                      weights = ~w)
    
    post_strat_reduc <- est_mean("diff_cces_on_cces", post_strat_reduc_svyd)
    
    ############################################
    ## Post-stratification: All
    ############################################
    post_strat_all_svyd = svydesign(~1, data = postStrat(survey_sim, 
                                                         cces_all_counts, "w_all", 
                                                         strata_pass = "strata_all"),
                                    weights = ~w)
    
    post_strat_all <- est_mean("diff_cces_on_cces", post_strat_all_svyd)
    
    ############################################
    #### Raking on true model
    ############################################
    #very messy error catching for the moment just a stop gap to see how things look
    rake_truth_svyd <- try(calibrate(design = survey_design,
                                     formula = selection_model,
                                     population = targets_demo_truth,
                                     calfun = "raking",
                                     epsilon = .009), silent = T)
    
    if ("try-error"%in%class(rake_truth_svyd)[1]){
      rake_truth <- NA
      truth_margins = NA
    } else {
      rake_truth <- est_mean("diff_cces_on_cces", rake_truth_svyd)
      truth_margins = svymean(margins_formula, rake_truth_svyd)
    }
    
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
      kbal_data <- data.frame(lapply(kbal_data, as.factor))
      kbal_data <- model.matrix(~ ., kbal_data,
                                contrasts.arg = lapply(kbal_data,
                                                       contrasts, contrasts=FALSE))
      kbal_data <- kbal_data[, -1]
      kbal_data_sampled <- c(rep(1, nrow(survey_sim)), rep(0, nrow(cces)))

      K <- makeK(kbal_data, b=2, useasbases = kbal_data_sampled,
                 linkernel = FALSE, scale = FALSE)
      raw_counts <- -log(K)

      #### now get alt b's
      #run optim: dplyr much faster than table
      n_d <- data.frame(diff = c(raw_counts)) %>% group_by(diff) %>% summarise(n())
      res = optimize(var_K, n_d= n_d, diag_length = length(diag(K)),
                     interval=c(0,2000), maximum=TRUE)
      #note with one hot encoding this needs to be * by 2
      b_maxvar <- res$maximum
      
      #bmean:
      #b_mean <- mean(as.vector(raw_counts))

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

      rm(K, res, raw_counts)
      #option to check more b values
      #b <- c(b_mean, b_maxvar*2, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048)
      #b <- c(b_mean, b_maxvar*2)
      b <- c(b_maxvar*2)
      
      kpop <- lapply(1:length(b), function(i) {
        
        #cat("===================== Running Kbal ==================\n")
        
        
        #### DEFAULT ######
        cat(paste("b:", b[i], "nsim:", nsim, "DEFAULT", "\n"))
        kbal_est <- kbal(allx=kbal_data,
                         sampled = kbal_data_sampled,
                         b = b[i],
                         scale_data = FALSE,
                         drop_multicollin = FALSE,
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
        
        #save memory by saving only the svd to re use
        svdK = kbal_est$svdK 
        numdims = kbal_est$numdims
        biasbound_r = kbal_est$biasbound.orig/kbal_est$biasbound.opt
        biasbound = kbal_est$biasbound.opt
        
        #CONVERGED
        dist_record = data.frame(t(kbal_est$dist.record))
        min_converged = dist_record[which.min(dist_record[dist_record$Ebal.Convergence ==1,"BiasBound"]), "Dims"]
        
        rm(kbal_est)
        
        #### CONVG ####
        cat(paste("b:", b[i], "nsim:", nsim, "CONV", "\n"))
        if(is.null(min_converged) | length(min_converged) ==0) {
          kpop_svyd_conv <- "dn converge"
          kpop_conv <- "dn converge"
          
          numdims_conv = "dn converge"
          biasbound_r_conv = "dn converge"
          biasbound_conv = "dn converge"
          
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
                                drop_multicollin = FALSE,
                                incrementby = increment,
                                meanfirst = FALSE,
                                sampledinpop = FALSE,
                                ebal.convergence = TRUE)
          kpop_svyd_conv <- svydesign(~1, data = survey_sim,
                                      weights = kbal_est_conv$w[kbal_data_sampled ==1])
          kpop_conv <- est_mean("diff_cces_on_cces", kpop_svyd_conv)
          
          numdims_conv = kbal_est_conv$numdims
          biasbound_r_conv = kbal_est_conv$biasbound.orig/kbal_est_conv$biasbound.opt
          biasbound_conv = kbal_est_conv$biasbound.opt
          rm(kbal_est_conv)
          
        }
        
        
        
        
        ####### MF #######
        cat(paste("b:", b[i], "nsim:", nsim, "MEANFIRST", "\n"))
        kbal_mf_est <- kbal(K.svd = svdK,
                            allx=kbal_data,
                            sampled = kbal_data_sampled,
                            ebal.tol = tolerance,
                            ebal.maxit = maxit,
                            minnumdims = min_num_dims,
                            maxnumdims = max_num_dims,
                            scale_data = FALSE,
                            drop_multicollin = FALSE,
                            incrementby = increment,
                            meanfirst = TRUE,
                            sampledinpop = FALSE)
        
        kpop_mf_svyd <- svydesign(~1, data = survey_sim, 
                                  weights = kbal_mf_est$w[kbal_data_sampled ==1])
        
        kpop_mf <- est_mean("diff_cces_on_cces", kpop_mf_svyd)
        
        mfnumdims = kbal_mf_est$numdims
        if(is.null(mfnumdims)) {mfnumdims = c("dn_converge") }
        mf_appended_dims = kbal_mf_est$meanfirst.dims
        if(is.null(mf_appended_dims)) {mf_appended_dims = c("dn_converge")}
        biasbound_r_mf = kbal_mf_est$biasbound.orig/kbal_mf_est$biasbound.opt
        biasbound_mf = kbal_mf_est$biasbound.opt
        
        rm(kbal_mf_est)
        
        
        #########demos constraint method:
        cat(paste("b:", b[i], "nsim:", nsim, "CONSTR", "\n"))
        kbal_demos_est <- kbal(K.svd = svdK,
                               allx=kbal_data,
                               sampled = kbal_data_sampled,
                               ebal.tol = tolerance,
                               ebal.maxit = maxit,
                               minnumdims = min_num_dims,
                               maxnumdims = max_num_dims,
                               scale_data = FALSE,
                               drop_multicollin = FALSE,
                               incrementby = increment,
                               #scaling these
                               constraint = rake_demos_constraint,
                               meanfirst = FALSE,
                               sampledinpop = FALSE)
        kpop_demos_svyd <- svydesign(~1, data = survey_sim, 
                                     weights = kbal_demos_est$w[kbal_data_sampled ==1])
        
        kpop_demos <- est_mean("diff_cces_on_cces", kpop_demos_svyd)
        
        numdims_demos = kbal_demos_est$numdims
        if(is.null(numdims_demos)) {numdims_demos = c("dn_converge") }
        biasbound_r_demos = kbal_demos_est$biasbound.orig/kbal_demos_est$biasbound.opt
        biasbound_demos = kbal_demos_est$biasbound.opt
        
        rm(kbal_demos_est)
        
        
        #########demos + educ constraint method:
        cat(paste("b:", b[i], "nsim:", nsim, "CONSTR", "\n"))
        kbal_demos_wedu_est <- kbal(K.svd = svdK,
                                    allx=kbal_data,
                                    sampled = kbal_data_sampled,
                                    ebal.tol = tolerance,
                                    ebal.maxit = maxit,
                                    minnumdims = min_num_dims,
                                    maxnumdims = max_num_dims,
                                    scale_data = FALSE,
                                    drop_multicollin = FALSE,
                                    incrementby = increment,
                                    #scaling these
                                    constraint = rake_demos_wedu_constraint,
                                    meanfirst = FALSE,
                                    sampledinpop = FALSE)
        kpop_demos_wedu_svyd <- svydesign(~1, data = survey_sim, 
                                          weights = kbal_demos_wedu_est$w[kbal_data_sampled ==1])
        
        kpop_demos_wedu <- est_mean("diff_cces_on_cces", kpop_demos_wedu_svyd)
        
        numdims_demos_wedu = kbal_demos_wedu_est$numdims
        if(is.null(numdims_demos_wedu)) {numdims_demos_wedu = c("dn_converge") }
        biasbound_r_demos_wedu = kbal_demos_wedu_est$biasbound.orig/kbal_demos_wedu_est$biasbound.opt
        biasbound_demos_wedu = kbal_demos_wedu_est$biasbound.opt
        
        rm(kbal_demos_wedu_est)
        
        
        #########all constraint method:
        cat(paste("b:", b[i], "nsim:", nsim, "CONSTR", "\n"))
        kbal_all_est <- kbal(K.svd = svdK,
                             allx=kbal_data,
                             sampled = kbal_data_sampled,
                             ebal.tol = tolerance,
                             ebal.maxit = maxit,
                             minnumdims = min_num_dims,
                             maxnumdims = max_num_dims,
                             scale_data = FALSE,
                             drop_multicollin = FALSE,
                             incrementby = increment,
                             #scaling these
                             constraint = rake_all_constraint,
                             meanfirst = FALSE,
                             sampledinpop = FALSE)
        kpop_all_svyd <- svydesign(~1, data = survey_sim, 
                                   weights = kbal_all_est$w[kbal_data_sampled ==1])
        
        kpop_all <- est_mean("diff_cces_on_cces", kpop_all_svyd)
        
        numdims_all = kbal_all_est$numdims
        if(is.null(numdims_all)) {numdims_all = c("dn_converge") }
        biasbound_r_all = kbal_all_est$biasbound.orig/kbal_all_est$biasbound.opt
        biasbound_all = kbal_all_est$biasbound.opt
        
        rm(kbal_all_est)
        
        rm(svdK)
        
        ##### return
        out = list()
        out$sims = data.frame(b[i],
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
        #weights
        out$weights = list(b[i],
                           kpop_w = weights(kpop_svyd),
                           kpop_w_conv = weights(kpop_svyd_conv),
                           kpop_mf_w = weights(kpop_mf_svyd), 
                           kpop_demos_w = weights(kpop_demos_svyd),
                           kpop_demos_wedu_w = weights(kpop_demos_wedu_svyd),
                           kpop_all_w = weights(kpop_all_svyd))
        
        ######## Kpop Margins ########
        
        out$km <- round(cbind(b = b[i]/100,
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
      
    } else {
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
                     as.data.frame(sapply(kpop, `[`, 3)))
      
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

save(sims,tolerance, maxit, increment, min_num_dims,
     file = paste0("./cat_sims_modeled_outcome_nodiag_",POPW, "_m",maxit, "_t",tolerance, "_inc",
                   increment, "mindims", min_num_dims,
                   Sys.Date(),
                   # str_sub(gsub("[[:space:]|[:punct:]]", "_",
                   #              gsub("[:alpha:]", "",
                   #                   Sys.time())),
                   #         start = 1, end = -3),
                   "_nsims", length(good),
                   ".RData"))


################# OUTCOME ESTIMATES ##################
#### load runs:
#load("./cat_sims_modeled_outcomeFALSE_m500_t1e-04_inc5mindims1502021-06-04_nsims176.RData")

good = which(lapply(sims, function (x) return(class(x))) == "list")
length(good)

plot <- lapply(sims, `[[`, 1)
#plot[-good]
#some of them may fail (rake truth dn converge usually)
plot <- plot[good] %>% bind_rows()

mean(plot$n)

plot_lasso_margin <- plot %>% 
  dplyr::select(unweighted, 
                kpop,
                kpop_conv,
                kpop_mf,
                kpop_demos,
                kpop_demos_wedu,
                kpop_all,
                rake_demos_noeduc,
                rake_demos_weduc,
                rake_all,
                #post_stratification,
                post_strat_reduc,
                #post_strat_all,
                rake_truth, 
                ht_truth, 
                hayek_truth
               ) %>% 
  pivot_longer(everything(),
               names_to = "estimator", 
               values_to = "margin") %>%
  mutate(#margin = (0.5 - margin) * 2 * 100,
    margin =  margin* 100,
    estimator_name = factor(case_when(
      estimator == "kpop" ~ "KPop",
      estimator == "kpop_mf" ~ "Kpop aMF (All)",
      estimator == "kpop_conv" ~ "Kpop Converged",
      estimator == "kpop_demos" ~ "Kpop MF (Demos)",
      estimator == "kpop_demos_wedu" ~ "Kpop MF (Demos+Edu)",
      estimator == "kpop_all" ~ "Kpop MF (All)",
      estimator == "rake_demos_noeduc" ~ "Raking\nDemos (No Educ)",
      estimator == "rake_demos_weduc" ~ "Raking\nDemos (w/ Educ)",
      estimator == "rake_demos_noeduc" ~ "Raking\nDemos (No Educ)",
      estimator == "rake_demos_weduc" ~ "Raking\nDemos (w/ Educ)",
      estimator == "rake_all" ~ "Raking\n All",
      estimator == "rake_truth" ~ "Raking\nTrue Selection\nModel",
      estimator == "post_stratification" ~ "Post-Strat Prev",
      estimator == "post_strat_reduc" ~ "Post-Strat Reduc",
      estimator == "post_strat_all" ~ "Post-Strat All",
      estimator == "unweighted" ~ "Unweighted", 
      estimator == "ebal_truth" ~ "Ebal Lasso Truth", 
      estimator == "ht_truth" ~ "Horvitz-Thompson Truth",
      estimator == "hayek_truth" ~ "Hayek Truth" ),
      levels = c("Unweighted", 
                 "Raking\nDemos (No Educ)",
                 "Raking\nDemos (w/ Educ)",
                 "Raking\n All",
                 "Post-Strat Prev", 
                 "Post-Strat Reduc", 
                 "Post-Strat All",
                 "Raking\nTrue Selection\nModel",
                 "Ebal Lasso Truth",
                 "Horvitz-Thompson Truth", 
                 "Hayek Truth",
                 "KPop",
                 "Kpop Converged",
                 "Kpop aMF (All)",
                 "Kpop MF (Demos)",
                 "Kpop MF (Demos+Edu)",
                 "Kpop MF (All)"
      )))

##BOX PLOT
ggplot(data = plot_lasso_margin,
       aes(x = estimator_name, y = margin)) +
  geom_boxplot(alpha = 0.2) +
  geom_hline(yintercept = margin_sim) +
  theme_bw() +
  ggtitle(paste0("Sim Results: (lasso) Testing",
                 " ", "pop weights: ", POPW, "\n",length(good), " sims, b=maxvar")) +
  xlab("Estimator") +
  ylab("Lasso Simulated Vote Margin") +
  annotate(geom = "text", x = 0.5, y = margin_sim,
           label = "True\nTarget\nPopulation\nMargin", hjust = 0) +
  theme(panel.grid.major.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "./sims_modout_fullrun.pdf", height = 4, width = 6)

table <- plot_lasso_margin %>% 
  mutate(estimator_name = gsub("\n|-", " ", estimator_name)) %>%
  group_by(estimator_name) %>%
  summarize(#truth = margin_sim,
    #M_est = mean(margin, na.rm = T),
    MSE = mean((margin - margin_sim)^2, na.rm = T),
    Bias = mean(margin - margin_sim, na.rm = T),
    SE = sd(margin, na.rm = T),
    MAE = mean(abs(margin - margin_sim), na.rm = T)
  ) %>%
  mutate(#MSE_ratio = mse / mse[estimator_name == "Unweighted"],
    Bias_reduc = 1- Bias / Bias[estimator_name == "Unweighted"]#,
    #se_ratio = se / se[estimator_name == "Unweighted"],
    #mae_ratio = mae / mae[estimator_name == "Unweighted"]
  ) %>%
  arrange(MSE)
table
#kable(table, format = "latex", digits = 3, caption = "lasso mod", booktabs = T)


sum(is.na(plot_lasso_margin %>% filter(estimator == "ebal_truth") %>% dplyr::select(margin)))








################# MARGINS #############################

margins <- lapply(sims, `[`, 2)
margins = margins[good]


test <- as.data.frame(margins) 
pop_weighted_avg = TRUE

#other methods margins
#cols are nsim results, rows are each individual level of each variable
#to summarize perfromance across simulations we can take an average of these
#might also want to take the sd to get some sense of how much variance we get
margins_base <- data.frame(
  sample = test[,grepl("margins.sample", colnames(test))] %>% rowMeans(),
  cces = test[,grepl("margins.cces", colnames(test))] %>% rowMeans(),
  rake_demos_noeduc = test[,grepl("margins.rake_demos_noeduc",
                                  colnames(test))] %>% rowMeans(),
  rake_demos_weduc = test[,grepl("margins.rake_demos_weduc",
                                 colnames(test))] %>% rowMeans(),
  rake_all = test[,grepl("margins.rake_all",
                         colnames(test))] %>% rowMeans(),
  post_stratification = test[,grepl("margins.post_stratification",
                                    colnames(test))] %>% rowMeans(),
  post_strat_reduc = test[,grepl("margins.post_strat_reduc",
                                 colnames(test))] %>% rowMeans(),
  post_strat_all = test[,grepl("margins.post_strat_all",
                               colnames(test))] %>% rowMeans(),
  rake_truth = test[,grepl("margins.rake_truth", colnames(test))] %>% rowMeans())

rownames(margins_base) <- sapply(rownames(test), function(x) { 
  if(grepl("recode", x)) {substr(x, 8,nchar(x))} else{x} } )


margins_avg <- data.frame(
  kpop = test[,grepl("margins.kpop_b", colnames(test))] %>% rowMeans(), 
  kpop_conv = test[,grepl("margins.kpop_cvg", colnames(test))] %>% rowMeans(), 
  kpop_mf = test[,grepl("margins.kpop_mf", colnames(test))] %>% rowMeans(),
  kpop_demos = test[,grepl("margins.kpop_demos", colnames(test))] %>% rowMeans(),
  kpop_demos_wedu = test[,grepl("margins.kpop_demos_wedu", colnames(test))] %>%
    rowMeans(),
  kpop_all = test[,grepl("margins.kpop_all", colnames(test))] %>% rowMeans()
)

rownames(margins_avg) <- sapply(rownames(test), function(x) { 
  if(grepl("recode", x)) {substr(x, 8,nchar(x))} else{x} } )

margins_total <- cbind(margins_base, margins_avg)

margins_diff <- margins_total %>% mutate(
  sample = cces - sample,
  kpop = cces - kpop,
  kpop_conv = cces - kpop_conv,
  kpop_mf = cces - kpop_mf,
  kpop_demos = cces - kpop_demos,
  kpop_demos_wedu = cces - kpop_demos_wedu,
  kpop_all = cces - kpop_all,
  rake_demos_noeduc = cces - rake_demos_noeduc,
  rake_demos_weduc = cces - rake_demos_weduc,
  rake_all = cces - rake_all,
  post_stratification = cces - post_stratification,
  post_strat_reduc = cces - post_strat_reduc,
  post_strat_all = cces - post_strat_all,
  rake_truth = cces - rake_truth,
  ebal_truth = cces - ebal_truth)
rownames(margins_diff) <- rownames(margins_base)


if(pop_weighted_avg) {
  w_margins_diff <- margins_diff*(margins_base[, "cces"]/100)
  #dont's square cces:
  w_margins_diff[, "cces"] <- margins_diff[, "cces"]
  #manugally fix the one level margins: then we'll just report the abs diff on these
  w_margins_diff["age",] <- margins_diff["age",]
  w_margins_diff["agesq",] <- margins_diff["agesq",]
  w_margins_diff["agecubed",] <- margins_diff["agecubed",]
  # w_margins_diff["vote_2016_sim",] <- margins_diff["vote_2016_sim",]
  w_margins_diff["diff_cces_on_cces",] <- margins_diff["diff_cces_on_cces",]
  w_margins_diff["mod_cces_on_cces_pR",] <- margins_diff["mod_cces_on_cces_pR",]
  w_margins_diff["mod_cces_on_cces_pO",] <- margins_diff["mod_cces_on_cces_pO",]
  w_margins_diff["mod_cces_on_cces_pD",] <- margins_diff["mod_cces_on_cces_pD",]
  
  
  wmabserr <- w_margins_diff %>% 
    mutate(var = sapply(rownames(margins_diff), 
                        function(y) {
                          substr(y, 1, 
                                 sapply(rownames(margins_diff), 
                                        function(x) {
                                          if(grepl("income",x)) {
                                            11 #income_5way
                                            #this is so ridiculous WHY (sunk cost ugh)
                                          } else if(grepl("mod", x) |grepl("sim", x) |
                                                    (grepl("age",x) & 
                                                     !grepl("bucket",x) &
                                                     !grepl("factor",x)) ){
                                            nchar(x) 
                                            #THIS IS SO STUPID
                                          } else if(grepl("factor",x)) {
                                            10 #age_factor
                                          } else {
                                            (gregexpr("[A-Z]|[1-9][1-9] |[1-9][1-9]\\+|[1-9]\\-", x)[[1]][1]-1)
                                          }
                                        })[y])
                        })) %>%
    group_by(var) %>% summarize(sample = mean(abs(sample)), 
                                kpop = mean(abs(kpop)),
                                #kpop_conv = mean(abs(kpop_conv)),
                                #kpop_mf = mean(abs(kpop_mf)),
                                kpop_demos = mean(abs(kpop_demos)),
                                rake_demos_noeduc = mean(abs(rake_demos_noeduc)),
                                kpop_demos_wedu = mean(abs(kpop_demos_wedu)),
                                rake_demos_weduc = mean(abs(rake_demos_weduc)), 
                                kpop_all = mean(abs(kpop_all)),
                                rake_all = mean(abs(rake_all)), 
                                post_stratification =mean(abs(post_stratification)),
                                post_strat_reduc =mean(abs(post_strat_reduc)),
                                post_strat_all =mean(abs(post_strat_all)),
                                rake_truth = mean(abs(rake_truth)) )
  wmabserr$var[wmabserr$var == "educ"]<- "educ_6way"
  
  
} else {
  #combining the factored age bins into a weighted L1 distance
  #weighting the difference by the proportion in that bin in cces
  #then we throw in n (number of bins) so that when i take the mean later we get the correct normalization
  margins_diff[(grep("age_factor", rownames(margins_diff))), ] <-
    (margins_diff[(grep("age_factor", 
                        rownames(margins_diff))),
    ])* (margins_total[grep("age_factor", rownames(margins_total)), "cces"]/100)*length(levels(cces$recode_age_factor))
  
  
  wmabserr <- margins_diff %>% 
    mutate(var = sapply(rownames(margins_diff), 
                        function(y) {
                          substr(y, 1, 
                                 sapply(rownames(margins_diff), 
                                        function(x) {
                                          if(grepl("income",x)) {
                                            11 #income_5way
                                            #this is so ridiculous WHY (sunk cost ugh)
                                            #keeping all of the outcomes and age margins
                                            #that are not buckets or factored
                                          } else if(grepl("mod", x) | grepl("sim", x) |
                                                    (grepl("age",x) & 
                                                     !grepl("bucket",x) &
                                                     !grepl("factor",x)) ){
                                            nchar(x) 
                                            #THIS IS SO STUPID
                                          } else if(grepl("factor",x)) {
                                            10 #age_factor
                                          } else {
                                            (gregexpr("[A-Z]|[1-9][1-9] |[1-9][1-9]\\+|[1-9]\\-", x)[[1]][1]-1)
                                          }
                                        })[y])
                        })) %>% 
    group_by(var)  %>% 
    summarize(sample = mean(abs(sample)), 
              #kpop = mean(abs(kpop)),
              #kpop_conv = mean(abs(kpop_conv)),
              #kpop_mf = mean(abs(kpop_mf)),
              #kpop_demos = mean(abs(kpop_demos)),
              rake_demos_noeduc = mean(abs(rake_demos_noeduc)),
              #kpop_demos_wedu = mean(abs(kpop_demos_wedu)),
              rake_demos_weduc = mean(abs(rake_demos_weduc)), 
              #kpop_all = mean(abs(kpop_all)),
              rake_all = mean(abs(rake_all)), 
              post_stratification =mean(abs(post_stratification)),
              post_strat_reduc =mean(abs(post_strat_reduc)),
              post_strat_all =mean(abs(post_strat_all)),
              rake_truth = mean(abs(rake_truth)),
              ebal_truth= mean(ebal_truth))
  
  wmabserr$var[wmabserr$var == "educ"]<- "educ_6way"
  
}

cbind(wmabserr[,1], round(wmabserr[,-1], 2)))

