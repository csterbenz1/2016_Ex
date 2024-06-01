### Packages
library(tidyverse)

library(survey)
library(parallel)
#devtools::install_github("csterbenz1/KBAL", ref = "cat_kernel") 
library(kbal)
library(glmnet)

if(detectCores() > 10) {
    path_data = "/nas/home/ciaras/kpop/data/"
} else {
    path_data= "/Users/Ciara_1/Dropbox/kpop/Updated/application/data/" 
}
PRE = TRUE
TEST = FALSE
DEBUG = TRUE
POPW = TRUE
SAVE = T
tolerance = 1e-6
maxit = 500
increment = 1
min_num_dims = 1
max_num_dims = 500
#use the manually specified range of lambdas in the ridge residuQ                                                                                           1Q2121111       21     alization or allow glmnet to choose internally?
manual_lambda = FALSE 
#T=lambda as that which minimizes cverror in residualization; F= 1 sd from min choice
lambda_min = FALSE 

formula_rake_demos_noeduc <- ~recode_age_bucket + recode_female + recode_race +
    recode_region + recode_pid_3way

#updated to include 6 way edu
formula_rake_demos_weduc <- ~recode_age_bucket + recode_female +
    recode_race + recode_region + recode_educ + recode_pid_3way

formula_rake_all_vars <- ~recode_age_bucket + recode_female +
    recode_race + recode_region + recode_pid_3way + recode_educ +
    recode_income_5way + recode_relig_6way + recode_born + recode_attndch_4way

############################ Standard Error Functions ######################
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
    return((n/(n-1))*sum((weights * residuals)^2) - 1/(n-1) * sum(weights * residuals)^2)
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

############################## Load Data ##########################
## SURVEY DATA (PEW)
### Load
if(PRE) {
    pew <- readRDS(paste0(path_data, "pew_lasso_061021.rds"))
} else {
    pew <- readRDS(paste0(path_data, "pew_post_CS_q4.rds"))
}

## AUXILIARY INFORMATION (CCES)
### Load
cces <- readRDS(paste0(path_data, "cces_lasso_061021.rds"))

if(TEST) {
    rs_cces = sample.int(nrow(cces),  (nrow(cces)/10)*1)
    cces =cces[rs_cces,]
    rs_pew = sample.int(nrow(pew),  (nrow(pew)/10)*1)
    pew = pew[rs_pew,]
} 

cces <- cces%>% 
    mutate(commonweight_vv_post = commonweight_vv_post/ mean(commonweight_vv_post))


####### DEFINE OUTCOME: important to do this to match what plot you're interested in bc to generate the SEs we need to run reg on the appropriate outcome
#outcome is now projected cces modeled vote margin
pew = pew %>% mutate(outcome = diff_cces_on_pew)
cces = cces %>% mutate(outcome = diff_cces_on_cces)


kbal_data <- bind_rows(pew %>% dplyr::select(recode_age_bucket,
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

kbal_data_sampled <- c(rep(1, nrow(pew)), rep(0, nrow(cces)))


##### expanding b:
vote_diff <- quote((recode_vote_2016Democrat - recode_vote_2016Republican) /(recode_vote_2016Democrat + recode_vote_2016Republican + recode_vote_2016Other)) 


########################### RUN KPOP ####################################################
#### kpop Default
kbal_est <- kbal(allx=kbal_data,
             sampled = kbal_data_sampled,
             cat_data = TRUE,
             incrementby = increment,
             meanfirst = FALSE,
             ebal.tol = tolerance,
             ebal.maxit = maxit,
             population.w = cces$commonweight_vv_post,
             minnumdims = min_num_dims,
             maxnumdims = max_num_dims,
             sampledinpop = FALSE,
             fullSVD = TRUE)

kpop_svyd <- svydesign(~1, data = pew,
                       weights = kbal_est$w[kbal_data_sampled ==1])
kpop <- svymean(~outcome, kpop_svyd,na.rm = TRUE)*100
b_kpop = kbal_est$b

if(DEBUG) {
    #path = "/Users/Ciara_1/Dropbox/kpop/2023/application/weights/DEBUG_"
    path = "./DEBUG_"
    save(kbal_est,
         file = paste0(path, "full_kbal_obj_", 
                       Sys.Date(), ".Rdata"))
}

l1_orig = ifelse(!is.null(kbal_est$L1_orig),kbal_est$L1_orig, NA)
l1 = ifelse(!is.null(kbal_est$L1_opt),kbal_est$L1_opt, NA)

#save memory by saving only the svd to re use
svdK = kbal_est$svdK 
numdims = kbal_est$numdims
biasbound_r = kbal_est$biasbound_ratio
biasbound = kbal_est$biasbound_opt

##### Kpop SEs
lambdas <- if(manual_lambda) { 10^seq(3, -2, by = -.1) } else {NULL}

x <- as.matrix(data.frame(kbal_dims = kbal_est$svdK$v[, 1:kbal_est$numdims]))
cv_fit <- cv.glmnet(x, kpop_svyd$variables$outcome, alpha = 0, lambda = lambdas)
lambda_pass = if(lambda_min) { cv_fit$lambda.min} else {cv_fit$lambda.1se}
residuals = kpop_svyd$variables$outcome - predict(cv_fit$glmnet.fit,
                                                  s = lambda_pass, newx = x)
res_kpop = data.frame(min = min(residuals), 
                      perc_25 = quantile(residuals, .25), 
                      mean = mean(residuals),
                      perc_75 = quantile(residuals, .75),
                      var = var(residuals))
kpop_se <- tryCatch(calc_SEs(Y = kpop_svyd$variables$outcome,
                             residuals = residuals,
                             pop_size = nrow(cces),
                             sample_size = nrow(pew),
                             weights = weights(kpop_svyd)), error = function(e) NA)

if(length(kpop_se) == 1) {
    kpop_se <- data.frame(SE_fixed = NA, 
                          SE_quasi = NA, 
                          SE_linear = NA, 
                          SE_chad = NA)
}
#names(kpop_se) = tryCatch(paste0("kpop_", names(kpop_se)), error = function(e) NA)
rownames(kpop_se) = "kpop"

######################################################
#KPOP CONVERGED
dist_record = data.frame(t(kbal_est$dist_record))
min_converged = dist_record[which.min(dist_record[dist_record$Ebal.Convergence ==1,"BiasBound"]), "Dims"]

rm(kbal_est)
if(is.null(min_converged) | length(min_converged) ==0) {
    kpop_conv_svyd <- "dn converge"
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
                          population.w = cces$commonweight_vv_post,
                          scale_data = FALSE,
                          drop_MC = FALSE,
                          incrementby = increment,
                          meanfirst = FALSE,
                          sampledinpop = FALSE,
                          ebal.convergence = TRUE)
    
    kpop_conv_svyd <- svydesign(~1, data = pew,
                           weights = kbal_est_conv$w[kbal_data_sampled ==1])
    kpop_conv <- svymean(~outcome, kpop_conv_svyd,na.rm = TRUE)*100
    
    numdims_conv = kbal_est_conv$numdims
    biasbound_r_conv = kbal_est_conv$biasbound_ratio
    biasbound_conv = kbal_est_conv$biasbound_opt
    l1_conv = ifelse(!is.null(kbal_est_conv$L1_opt), kbal_est_conv$L1_opt, NA)
    #SEs
    x <- as.matrix(data.frame(kbal_dims = kbal_est_conv$svdK$v[, 1:kbal_est_conv$numdims]))
    cv_fit <- cv.glmnet(x, kpop_conv_svyd$variables$outcome, alpha = 0, 
                        lambda = lambdas)
    fit <- cv_fit$glmnet.fit
    lambda_pass = if(lambda_min) { cv_fit$lambda.min} else {cv_fit$lambda.1se}
    residuals = kpop_conv_svyd$variables$outcome - predict(cv_fit$glmnet.fit, 
                                                           s = lambda_pass, 
                                                           newx = x)
    res_kpop_conv = data.frame(min = min(residuals), 
                               perc_25 = quantile(residuals, .25), 
                               mean = mean(residuals),
                               perc_75 = quantile(residuals, .75),
                               var = var(residuals))
    kpop_conv_se <- tryCatch(calc_SEs(Y = kpop_conv_svyd$variables$outcome,
                                      residuals = residuals,
                                      pop_size = nrow(cces),
                                      sample_size = nrow(pew),
                                      weights = weights(kpop_conv_svyd)), error = function(e) NA)
    if(length(kpop_conv_se) == 1) {
        kpop_conv_se <- data.frame(SE_fixed = NA, 
                                   SE_quasi = NA, 
                                   SE_linear = NA, 
                                   SE_chad = NA)
    }
    #names(kpop_conv_se) = tryCatch(paste0("kpop_conv_", names(kpop_conv_se)), error = function(e) NA)
    #KRLS SEs are exactly the same for coverged
    rownames(kpop_conv_se) = "kpop_conv"
    
    rm(kbal_est_conv) 
}

########################################################################

####### MF #######
kbal_mf_est <- kbal(K.svd = svdK,
                    cat_data = T,
                    allx=kbal_data,
                    sampled = kbal_data_sampled,
                    ebal.tol = tolerance,
                    ebal.maxit = maxit,
                    minnumdims = min_num_dims,
                    maxnumdims = max_num_dims,
                    population.w = cces$commonweight_vv_post,
                    incrementby = increment,
                    meanfirst = TRUE,
                    sampledinpop = FALSE)

kpop_mf_svyd <- svydesign(~1, data = pew, 
                          weights = kbal_mf_est$w[kbal_data_sampled ==1])

kpop_mf <- svymean(~outcome, kpop_mf_svyd, na.rm = TRUE)*100

mfnumdims = kbal_mf_est$numdims
mf_appended_dims = kbal_mf_est$meanfirst_dims
if(is.null(mf_appended_dims)) {mf_appended_dims = c(NA)}
biasbound_r_mf = kbal_mf_est$biasbound_ratio
biasbound_mf = kbal_mf_est$biasbound_opt
l1_mf = ifelse(!is.null(kbal_mf_est$L1_opt), kbal_mf_est$L1_opt, NA)

mfnumdims = kbal_mf_est$numdims
if(is.null(mfnumdims)) {
    mfnumdims = c(NA) 
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
    
    cv_fit <- cv.glmnet(X, kpop_mf_svyd$variables$outcome, alpha = 0,
                        lambda = lambdas, 
                        penalty.factor = c(rep(0, kbal_mf_est$meanfirst_dims), 
                                           rep(1, kbal_mf_est$numdims)))
    lambda_pass = if(lambda_min) { cv_fit$lambda.min} else {cv_fit$lambda.1se}
    residuals = kpop_mf_svyd$variables$outcome - predict(cv_fit$glmnet.fit, 
                                                         s = lambda_pass, 
                                                         newx = X)
    res_kpop_mf = data.frame(min = min(residuals), 
                             perc_25 = quantile(residuals, .25), 
                             mean = mean(residuals),
                             perc_75 = quantile(residuals, .75),
                             var = var(residuals))
    
    kpop_mf_se <- tryCatch(calc_SEs(Y = kpop_mf_svyd$variables$outcome,
                                    residuals = residuals,
                                    pop_size = nrow(cces),
                                    sample_size = nrow(pew),
                                    weights = weights(kpop_mf_svyd)), 
                           error = function(e) NA)
    if(length(kpop_mf_se) == 1) {
        kpop_mf_se <- data.frame(SE_fixed = NA, 
                                 SE_quasi = NA, 
                                 SE_linear = NA, 
                                 SE_chad = NA)
    }
    # names(kpop_mf_se) = tryCatch(paste0("kpop_mf_", names(kpop_mf_se)),
    #                              error = function(e) NA)
    rownames(kpop_mf_se) = "kpop_mf"
}

rm(kbal_mf_est)

#########KPOP + demos constraint method:
kbal_demos_est <- kbal(K.svd = svdK,
                       allx=kbal_data,
                       cat_data = TRUE,
                       sampled = kbal_data_sampled,
                       ebal.tol = tolerance,
                       ebal.maxit = maxit,
                       minnumdims = min_num_dims,
                       maxnumdims = max_num_dims,
                       population.w = cces$commonweight_vv_post,
                       scale_data = FALSE,
                       drop_MC = FALSE,
                       incrementby = increment,
                       #scaling these
                       #constraint = rake_demos_constraint,
                       #scale_constraint = TRUE,
                       meanfirst = TRUE,
                       mf_columns = all.vars(formula_rake_demos_noeduc),
                       sampledinpop = FALSE)

kpop_demos_svyd <- svydesign(~1, data = pew, 
                             weights = kbal_demos_est$w[kbal_data_sampled ==1])

kpop_demos <- svymean(~outcome, kpop_demos_svyd, na.rm = TRUE)*100
numdims_demos = kbal_demos_est$numdims
l1_demos = ifelse(!is.null(kbal_demos_est$L1_opt), kbal_demos_est$L1_opt, NA)

if(is.null(numdims_demos)) {
    numdims_demos = c(NA) 
    kpop_demos_se <- data.frame(SE_fixed = NA, 
                                SE_quasi = NA, 
                                SE_linear = NA, 
                                SE_chad = NA)
} else {
    V <-  data.frame(kbal_dims = kbal_demos_est$svdK$v[, c(1:kbal_demos_est$numdims)])
    X <- as.matrix(cbind(kbal_demos_est$appended_constraint_cols[kbal_data_sampled==1, ], V))
    
    cv_fit <- cv.glmnet(X, kpop_demos_svyd$variables$outcome, alpha = 0, lambda = lambdas,
                        penalty.factor = c(rep(0, ncol(kbal_demos_est$appended_constraint_cols)), rep(1, kbal_demos_est$numdims)))
    
    lambda_pass = if(lambda_min) { cv_fit$lambda.min} else {cv_fit$lambda.1se}
    residuals =  kpop_demos_svyd$variables$outcome - predict(cv_fit$glmnet.fit,
                                                             s = lambda_pass, 
                                                             newx = X)
    res_kpop_demos = data.frame(min = min(residuals), 
                                perc_25 = quantile(residuals, .25), 
                                mean = mean(residuals),
                                perc_75 = quantile(residuals, .75),
                                var = var(residuals))
    
    kpop_demos_se <- tryCatch(calc_SEs(Y = kpop_demos_svyd$variables$outcome,
                                       residuals = residuals,
                                       pop_size = nrow(cces),
                                       sample_size = nrow(pew),
                                       weights = weights(kpop_demos_svyd)), 
                              error = function(e) NA)
    if(length(kpop_demos_se) == 1) {
        kpop_demos_se <- data.frame(SE_fixed = NA, 
                                    SE_quasi = NA, 
                                    SE_linear = NA, 
                                    SE_chad = NA)
    }
    # names(kpop_demos_se) = tryCatch(paste0("kpop_demos_", names(kpop_demos_se)),
    #                                 error = function(e) NA)
    rownames(kpop_demos_se) = "kpop_demos"
}
biasbound_r_demos = kbal_demos_est$biasbound_ratio
biasbound_demos = kbal_demos_est$biasbound_opt

rm(kbal_demos_est)


#########KPOP demos + educ constraint method:
kbal_demos_wedu_est <- kbal(K.svd = svdK,
                            allx=kbal_data,
                            cat_data = TRUE,
                            sampled = kbal_data_sampled,
                            ebal.tol = tolerance,
                            ebal.maxit = maxit,
                            minnumdims = min_num_dims,
                            maxnumdims = max_num_dims,
                            population.w = cces$commonweight_vv_post,
                            scale_data = FALSE,
                            drop_MC = FALSE,
                            incrementby = increment,
                            #scaling these
                            #constraint = rake_demos_wedu_constraint,
                            #scale_constraint = TRUE,
                            meanfirst = TRUE,
                            mf_columns = all.vars(formula_rake_demos_weduc),
                            sampledinpop = FALSE)

kpop_demos_wedu_svyd <- svydesign(~1, data = pew, 
                             weights = kbal_demos_wedu_est$w[kbal_data_sampled ==1])
kpop_demos_wedu <- svymean(~outcome, kpop_demos_wedu_svyd, na.rm = TRUE)*100
numdims_demos_wedu = kbal_demos_wedu_est$numdims
l1_demos_weduc = ifelse(!is.null(kbal_demos_wedu_est$L1_opt), kbal_demos_wedu_est$L1_opt, NA)

if(is.null(numdims_demos_wedu)) {
    numdims_demos_wedu = c(NA)
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
    
    lambda_pass = if(lambda_min) { cv_fit$lambda.min} else {cv_fit$lambda.1se}
    residuals =  kpop_demos_wedu_svyd$variables$outcome - predict(cv_fit$glmnet.fit,
                                                                  s = lambda_pass, 
                                                                  newx = X)
    res_kpop_demos_wedu = data.frame(min = min(residuals), 
                                     perc_25 = quantile(residuals, .25), 
                                     mean = mean(residuals),
                                     perc_75 = quantile(residuals, .75),
                                     var = var(residuals))
    kpop_demos_wedu_se <- tryCatch(calc_SEs(Y = kpop_demos_wedu_svyd$variables$outcome,
                                            residuals = residuals,
                                            pop_size = nrow(cces),
                                            sample_size = nrow(pew),
                                            weights = weights(kpop_demos_wedu_svyd)), 
                                   error = function(e) NA)
    if(length(kpop_demos_wedu_se) == 1) {
        kpop_demos_wedu_se <- data.frame(SE_fixed = NA, 
                                         SE_quasi = NA, 
                                         SE_linear = NA, 
                                         SE_chad = NA)
    }
    # names(kpop_demos_wedu_se) = tryCatch(paste0("kpop_demos_wedu_", names(kpop_demos_wedu_se)),
    #                                      error = function(e) NA)
    rownames(kpop_demos_wedu_se) = "kpop_demos_wedu"
    
}
biasbound_r_demos_wedu = kbal_demos_wedu_est$biasbound_ratio
biasbound_demos_wedu = kbal_demos_wedu_est$biasbound_opt

rm(kbal_demos_wedu_est)


#########KPOP + all constraint method:
kbal_all_est <- kbal(K.svd = svdK,
                     allx=kbal_data,
                     cat_data = TRUE,
                     sampled = kbal_data_sampled,
                     ebal.tol = tolerance,
                     ebal.maxit = maxit,
                     minnumdims = min_num_dims,
                     maxnumdims = max_num_dims,
                     population.w = cces$commonweight_vv_post,
                     scale_data = FALSE,
                     drop_MC = FALSE,
                     incrementby = increment,
                     #scaling these
                     #constraint = rake_all_constraint,
                     #scale_constraint = TRUE,
                     meanfirst = TRUE,
                     mf_columns = all.vars(formula_rake_all_vars),
                     sampledinpop = FALSE)

kpop_all_svyd <- svydesign(~1, data = pew, 
                           weights = kbal_all_est$w[kbal_data_sampled ==1])
kpop_all <- svymean(~outcome, kpop_all_svyd, na.rm = TRUE)*100
numdims_all = kbal_all_est$numdims
l1_all = ifelse(!is.null(kbal_all_est$L1_opt),  kbal_all_est$L1_opt, NA)

if(is.null(numdims_all)) {
    numdims_all = c(NA)
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
    
    lambda_pass = if(lambda_min) { cv_fit$lambda.min} else {cv_fit$lambda.1se}
    residuals =  kpop_all_svyd$variables$outcome - predict(cv_fit$glmnet.fit,
                                                           s = lambda_pass, 
                                                           newx = X)
    res_kpop_all = data.frame(min = min(residuals), 
                              perc_25 = quantile(residuals, .25), 
                              mean = mean(residuals),
                              perc_75 = quantile(residuals, .75),
                              var = var(residuals))
    kpop_all_se <- tryCatch(calc_SEs(Y = kpop_all_svyd$variables$outcome,
                                     residuals = residuals,
                                     pop_size = nrow(cces),
                                     sample_size = nrow(pew),
                                     weights = weights(kpop_all_svyd)), 
                            error = function(e) NA)
    if(length(kpop_all_se) == 1) {
        kpop_all_se <- data.frame(SE_fixed = NA, 
                                         SE_quasi = NA, 
                                         SE_linear = NA, 
                                         SE_chad = NA)
    }
    # names(kpop_all_se) = tryCatch(paste0("kpop_all_", names(kpop_all_se)),
    #                               error = function(e) NA)
    rownames(kpop_all_se) = "kpop_all"
    
}

biasbound_r_all = kbal_all_est$biasbound_ratio
biasbound_all = kbal_all_est$biasbound_opt
rm(kbal_all_est)

rm(svdK)



##### return
out = list()
b_out = b_kpop

out$est = data.frame(b_out,
                     tolerance = tolerance, 
                     maxit = maxit,
                     POPW = POPW,
                      kpop = kpop,
                      kpop_mf = kpop_mf,
                      kpop_conv = kpop_conv,
                      kpop_demos = kpop_demos,
                      kpop_demos_wedu = kpop_demos_wedu,
                      kpop_all = kpop_all,
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
                      numdims_all,
                      l1_orig ,
                      l1,
                      l1_conv,
                      l1_mf,
                      l1_demos,
                      l1_demos_weduc,
                      l1_all)

#Standard Errors:
out$SEs = rbind(kpop_se,
                     kpop_conv_se,
                     kpop_mf_se,
                     kpop_demos_se,
                     kpop_demos_wedu_se,
                     kpop_all_se)

#weights
out$weights = list(b = b_out,
                   kpop_w = weights(kpop_svyd),
                   kpop_w_conv = weights(kpop_conv_svyd),
                   kpop_mf_w = weights(kpop_mf_svyd), 
                   kpop_demos_w = weights(kpop_demos_svyd),
                   kpop_demos_wedu_w = weights(kpop_demos_wedu_svyd),
                   kpop_all_w = weights(kpop_all_svyd))

#residuals
out$residuals = rbind(b = b_out,
                      kpop_res = res_kpop,
                      kpop_w_conv = res_kpop_conv,
                      kpop_mf_w = res_kpop_mf,
                      kpop_demos_w = res_kpop_demos,
                      kpop_demos_wedu_w = res_kpop_demos_wedu,
                      kpop_all_w = res_kpop_all
)


if(SAVE) {
    if(TEST) {
        #path = "/Users/Ciara_1/Dropbox/kpop/2023/application/weights/TEST_"
        path = "./TEST_"
        out$test_sample = list(rs_cces = rs_cces, 
                          rs_pew = rs_pew)
    } else {
        #path = "/Users/Ciara_1/Dropbox/kpop/2023/application/weights/"
        path = "./"
    }
    save(out, tolerance, maxit, POPW,min_num_dims, max_num_dims, increment,TEST,
         file = paste0(path, "app_update_june", 
                       Sys.Date(), ".Rdata"))
}
 