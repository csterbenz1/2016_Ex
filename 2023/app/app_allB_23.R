### Packages
library(tidyverse)

library(survey)
library(parallel)
#devtools::install_github("csterbenz1/KBAL", ref = "cat_kernel") 
library(kbal)
library(glmnet)
#range of b:
b <- c(2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048)
if(detectCores() > 10) {
    path_data = "/nas/home/ciaras/kpop/data/"
    cores_used = length(b)
} else {
    path_data= "/Users/Ciara_1/Dropbox/kpop/Updated/application/data/" 
    cores_used = detectCores()
}
PRE = TRUE
TEST = FALSE
POPW = TRUE
SAVE = TRUE
tolerance = 1e-6
maxit = 500
increment = 1
min_num_dims = 1
max_num_dims = 500
#use the manually specified range of lambdas in the ridge residualization or allow glmnet to choose internally?
manual_lambda = FALSE 
#T=lambda as that which minimizes cverror in residualization; F= 1 sd from min choice
lambda_min = FALSE 



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

cces <- cces%>% 
  mutate(commonweight_vv_post = commonweight_vv_post/ mean(commonweight_vv_post))

if(TEST) {
    rs_cces = sample.int(nrow(cces),  (nrow(cces)/10)*1)
    cces =cces[rs_cces,]
    rs_pew = sample.int(nrow(pew),  (nrow(pew)/10)*1)
    pew = pew[rs_pew,]
    increment = 10
    b = c(2, 4)
} 

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


cces <- cces%>% 
    mutate(commonweight_vv_post = commonweight_vv_post/ mean(commonweight_vv_post))


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

b
cat_app <- mclapply(1:length(b), function(i) {
  
  cat("==========================================================================\n")
  cat(b[i], "\n")
  cat("################# KPOP DEFAULT:", b[i],"####################\n")
  kbal_est <- kbal(allx=kbal_data,
               sampled = kbal_data_sampled,
               b = b[i],
               cat_data = TRUE,
               incrementby = increment,
               meanfirst = FALSE,
               ebal.tol = tolerance,
               ebal.maxit = maxit,
               population.w = if(POPW) {cces$commonweight_vv_post} else {NULL},
               minnumdims = min_num_dims,
               maxnumdims = max_num_dims,
               sampledinpop = FALSE,
               fullSVD = TRUE)
  
  kpop_svyd <- svydesign(~1, data = pew,
                        weights = kbal_est$w[kbal_data_sampled ==1])
  kpop <- svymean(~outcome, kpop_svyd,na.rm = TRUE)*100
  b_kpop = kbal_est$b
  
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
  
  rownames(kpop_se) = "kpop"
  
  dist_record = data.frame(t(kbal_est$dist_record))
  min_converged = dist_record[which.min(dist_record[dist_record$Ebal.Convergence ==1,"BiasBound"]), "Dims"]
  rm(kbal_est)
  
  if(is.null(min_converged) | length(min_converged) ==0) {
    cat(paste("KPOP DEFAULT HAS NO CONVERGENT DIMS FOR", b[i], "\n"))
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
    cat("################# KPOP CONV:", b[i],"####################\n")
    kbal_est_conv <- kbal(allx=kbal_data,
                          K.svd = svdK,
                          sampled = kbal_data_sampled,
                          numdims = min_converged,
                          ebal.tol = tolerance,
                          ebal.maxit = maxit,
                          minnumdims = min_num_dims,
                          maxnumdims = max_num_dims,
                          population.w = if(POPW) {cces$commonweight_vv_post} else {NULL},
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
    
    rownames(kpop_conv_se) = "kpop_conv"
    
    rm(kbal_est_conv)
    
  }

  ####### MF #######
  cat("################# KPOP aMF:", b[i],"####################\n")
  kbal_mf_est <- kbal(K.svd = svdK,
                      cat_data = T,
                      allx=kbal_data,
                      sampled = kbal_data_sampled,
                      ebal.tol = tolerance,
                      ebal.maxit = maxit,
                      minnumdims = min_num_dims,
                      maxnumdims = max_num_dims,
                      population.w = if(POPW) {cces$commonweight_vv_post} else {NULL},
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
      res_kpop_mf = NA
  } else {
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
      rownames(kpop_mf_se) = "kpop_mf"
  }
  
  rm(kbal_mf_est)
 

  out = list()
  b_out = b_kpop
  
  out$est = data.frame(b_out,
                       tolerance = tolerance, 
                       maxit = maxit,
                       POPW = POPW,
                       kpop = kpop,
                       kpop_mf = kpop_mf,
                       bb = biasbound,
                       bbr = biasbound_r,
                       bb_conv = biasbound_conv,
                       bbr_conv = biasbound_r_conv,
                       bb_mf = biasbound_mf,
                       bbr_mf = biasbound_r_mf,
                       numdims,
                       numdims_conv,
                       mfnumdims, 
                       mf_appended_dims, 
                       l1_orig ,
                       l1,
                       l1_mf)
  
  #Standard Errors:
  out$SEs = rbind(kpop_se,
                  kpop_conv_se,
                  kpop_mf_se)
  
  #weights
  out$weights = list(b = b_out,
                     kpop_w = weights(kpop_svyd),
                     kpop_w_conv = weights(kpop_conv_svyd),
                     kpop_mf_w = weights(kpop_mf_svyd))
  
  #residuals
  out$residuals = rbind(b =b_out,
                        kpop_res = res_kpop,
                        kpop_res_conv = res_kpop_conv,
                        kpop_res_mf = res_kpop_mf)
  return(out)
  
},  mc.cores = cores_used)


popw_lab = ifelse(POPW, "POPW", NULL)


if(SAVE) {
    if(TEST) {
        path = "/Users/Ciara_1/Dropbox/kpop/2023/application/weights/TEST_"
        out$test_sample = list(rs_cces = rs_cces, 
                               rs_pew = rs_pew)
    } else {
        path = "./"
    }
    cat("SAVING TO: ", path, "\n")
    save(cat_app, b, tolerance, maxit, POPW, min_num_dims, max_num_dims,increment, TEST,
         file = paste0(path, "app_allB", popw_lab, "_",
                       Sys.Date(), ".Rdata"))
}
