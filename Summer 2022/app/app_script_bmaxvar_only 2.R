
### Packages
library(tidyverse)
library(survey)
library(parallel)
#devtools::install_github("csterbenz1/KBAL", ref = "cat_kernel") 
library(kbal)

if(detectCores() > 10) {
  path_data= "/home/csterbenz/Data/"
} else {
  path_data= "/Users/Ciara/Dropbox/kpop/Updated/application/data/" 
}
PRE = TRUE
TEST = FALSE
POPW = TRUE
tolerance = 1e-6
maxit = 500
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
  subset =sample.int(nrow(cces),  (nrow(cces)/10)*1)
} else {
  subset = c(1:nrow(cces))
}


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
                       cces[subset, ] %>% dplyr::select(recode_age_bucket,
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
                          contrasts.arg = lapply(kbal_data[,], contrasts, contrasts=FALSE))
kbal_data <- kbal_data[, -1]
#colnames(kbal_data)
#nrow(kbal_data)

kbal_data_sampled <- c(rep(1, nrow(pew)), rep(0, nrow(cces[subset,])))

#categorical kernel + b range:
#get raw counts:
K <- makeK(kbal_data, b=2, useasbases = kbal_data_sampled,
           linkernel = FALSE, scale = FALSE)

raw_counts <- -log(K)
#check this:
# K[1:3, 1:3]
# exp(-sum((kbal_data[1, ] - kbal_data[2,])^2)/2)
# exp(-sum((kbal_data[2, ] - kbal_data[3,])^2)/2)
# raw_counts[1:3 ,1:3]
# sum((kbal_data[1, ] - kbal_data[2,])^2)
# sum((kbal_data[1, ] - kbal_data[3,])^2)
# 


#the b you pass in/the b we get out when optimizing has this factor of 2 issue
#b=2 K corresponds to var_K b =1
# i think this is ok because we mulitply the outputted b by 2 for the max var:

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

n_d <- data.frame(diff = c(raw_counts)) %>% group_by(diff) %>% summarise(n())
res = optimize(var_K, n_d, length(diag(K)),
               interval=c(0,2000), maximum=TRUE)


b_maxvar <- res$maximum
#b_maxvar*2
#check:
# K_max <- makeK(kbal_data, b=b_maxvar*2, useasbases = kbal_data_sampled,
#                linkernel = FALSE, scale = FALSE)
# diag(K_max) <- 0
# var(as.vector(K_max))
# res$objective



# K_maxvar <- makeK(kbal_data, b=b_maxvar*2, useasbases = kbal_data_sampled,
#                    linkernel = FALSE, scale = FALSE)
# hist(as.vector(K_maxvar))
# summary(as.vector(K_maxvar))
# dat <- data.frame(raw_counts = c(0:max(raw_counts))) %>% mutate(  K = exp(-raw_counts/b_maxvar))
# ggplot(dat, aes(x= raw_counts, y = K)) +
#     geom_point() +
#     scale_x_continuous(breaks = c(dat$raw_counts))
# hist(as.vector(K))
# hist(as.vector(raw_counts))

#we want raw_counts/b_maxvar = double_coutns/(2*b_maxvar)

#bmean:
b_mean <- mean(as.vector(raw_counts))
# we want: raw_counts/(.5*mean) = double_counts/mean

##### Demos Constraint
rake_demos_constraint <- bind_rows(pew %>% dplyr::select(recode_age_bucket,
                                                         recode_female,
                                                         recode_race,
                                                         recode_region,
                                                         recode_pid_3way),
                                   cces[subset, ] %>% dplyr::select(recode_age_bucket,
                                                                    recode_female,
                                                                    recode_race,
                                                                    recode_region,
                                                                    recode_pid_3way))%>%
  model.matrix(as.formula("~."), .)

rake_demos_constraint <- rake_demos_constraint[,-1]
rake_demos_constraint <- scale(rake_demos_constraint)


rake_demos_wedu_constraint <- bind_rows(pew %>% dplyr::select(recode_age_bucket,
                                                         recode_female,
                                                         recode_race,
                                                         recode_region,
                                                         recode_pid_3way,
                                                         recode_educ),
                                   cces[subset, ] %>% dplyr::select(recode_age_bucket,
                                                                    recode_female,
                                                                    recode_race,
                                                                    recode_region,
                                                                    recode_pid_3way,
                                                                    recode_educ))%>%
  model.matrix(as.formula("~."), .)

rake_demos_wedu_constraint <- rake_demos_wedu_constraint[,-1]
rake_demos_wedu_constraint <- scale(rake_demos_wedu_constraint)


rake_all_constraint <- bind_rows(pew %>% dplyr::select(recode_age_bucket,
                                                       recode_female,
                                                       recode_race,
                                                       recode_region,
                                                       recode_pid_3way,
                                                       recode_educ,
                                                       
                                                       recode_income_5way,
                                                       recode_relig_6way,
                                                       recode_born,
                                                       recode_attndch_4way),
                                 cces[subset, ] %>% dplyr::select(recode_age_bucket,
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

#check for collinearity
# qr_all <- qr(rake_all_constraint)
# qr_all$rank
# ncol(rake_all_constraint)


##### expanding b:
vote_diff <- quote((recode_vote_2016Democrat - recode_vote_2016Republican) /(recode_vote_2016Democrat + recode_vote_2016Republican + recode_vote_2016Other)) 

#expanding b
b <- c( b_maxvar*2)
rm(K, raw_counts)
cat_app <- mclapply(1:length(b), function(i) {
  
  cat("###############################################################\n")
  cat(b[i], "\n")
  kbal <- kbal(allx=kbal_data,
               sampled = kbal_data_sampled,
               b = b[i],
               scale_data = FALSE,
               drop_multicollin = FALSE,
               incrementby = 1,
               meanfirst = FALSE,
               population.w = if(POPW) {cces$commonweight_vv_post} else {NULL},
               ebal.tol = tolerance,
               ebal.maxit = maxit,
               maxnumdims = 500,
               #linkernel = if(TEST){TRUE} else{ FALSE},
               sampledinpop = FALSE,
               fullSVD = TRUE)
  
  
  dist_record = data.frame(t(kbal$dist.record))
  min_converged = dist_record[which.min(dist_record[dist_record$Ebal.Convergence ==1,
                                                    "BiasBound"]), "Dims"]
  
  if(is.null(min_converged) | length(min_converged) ==0) {
    cat(paste("KPOP DEFAULT HAS NO CONVERGENT DIMS FOR", b[i], "\n"))
    kbal_conv <- kbal(allx=kbal_data,
                      K.svd = kbal$svdK,
                      K= kbal$K,
                      sampled = kbal_data_sampled,
                      scale_data = FALSE,
                      drop_multicollin = FALSE,
                      incrementby = 1,
                      maxnumdims = 100,
                      population.w = if(POPW) {cces$commonweight_vv_post} else {NULL},
                      ebal.tol = tolerance,
                      ebal.maxit = maxit,
                      meanfirst = FALSE,
                      sampledinpop = FALSE,
                      ebal.convergence = TRUE)
  } else {
    kbal_conv <- kbal(allx=kbal_data,
                      K.svd = kbal$svdK,
                      K= kbal$K,
                      sampled = kbal_data_sampled,
                      numdims = min_converged,
                      ebal.tol = tolerance,
                      ebal.maxit = maxit,
                      population.w = if(POPW) {cces$commonweight_vv_post} else {NULL},
                      maxnumdims = 500,
                      scale_data = FALSE,
                      drop_multicollin = FALSE,
                      incrementby = 1,
                      meanfirst = FALSE,
                      sampledinpop = FALSE,
                      ebal.convergence = TRUE)
  }
  
  kbal_mf <- kbal(K.svd = kbal$svdK,
                  K= kbal$K,
                  allx=kbal_data,
                  sampled = kbal_data_sampled,
                  ebal.tol = tolerance,
                  ebal.maxit = maxit,
                  maxnumdims = 500,
                  population.w = if(POPW) {cces$commonweight_vv_post} else {NULL},
                  scale_data = FALSE,
                  drop_multicollin = FALSE,
                  incrementby = 1,
                  meanfirst = TRUE,
                  sampledinpop = FALSE)
  
  kbal_demos <- kbal(K.svd = kbal$svdK,
                     K= kbal$K,
                     allx=kbal_data,
                     fullSVD = TRUE,
                     sampled = kbal_data_sampled,
                     ebal.tol = tolerance,
                     ebal.maxit = maxit,
                     population.w = if(POPW) {cces$commonweight_vv_post} else {NULL},
                     maxnumdims = 500,
                     scale_data = FALSE,
                     drop_multicollin = FALSE,
                     incrementby = 1,
                     #scaling these
                     constraint = rake_demos_constraint,
                     scale_constraint = TRUE,
                     meanfirst = FALSE,
                     sampledinpop = FALSE)
  
  kbal_demos_weduc <- kbal(K.svd = kbal$svdK,
                           K= kbal$K,
                           allx=kbal_data,
                           fullSVD = TRUE,
                           sampled = kbal_data_sampled,
                           ebal.tol = tolerance,
                           ebal.maxit = maxit,
                           population.w = if(POPW) {cces$commonweight_vv_post} else {NULL},
                           maxnumdims = 500,
                           scale_data = FALSE,
                           drop_multicollin = FALSE,
                           incrementby = 1,
                           #scaling these above
                           constraint = rake_demos_wedu_constraint,
                           scale_constraint = TRUE,
                           meanfirst = FALSE,
                           sampledinpop = FALSE)
  
  
  kbal_full <- kbal(K.svd = kbal$svdK,
                    K= kbal$K,
                    allx=kbal_data,
                    fullSVD = TRUE,
                    sampled = kbal_data_sampled,
                    ebal.tol = tolerance,
                    ebal.maxit = maxit,
                    scale_data = FALSE,
                    population.w = if(POPW) {cces$commonweight_vv_post} else {NULL},
                    maxnumdims = 500,
                    drop_multicollin = FALSE,
                    incrementby = 1,
                    constraint = rake_all_constraint,
                    scale_constraint = TRUE,
                    meanfirst = FALSE,
                    sampledinpop = FALSE)
  
  out <- list()
  # out$kernel <- list(b = b[i],
  #                    K = kbal$K)
  
  out$weights <- list(b = b[i],
                      popw = POPW,
                      tolerance = tolerance, 
                      maxit = maxit, 
                      kpop_w = kbal$w,
                      kpop_w_conv = kbal_conv$w,
                      kpop_mf_w = kbal_mf$w, 
                      kpop_demos_w = kbal_demos$w,
                      kpop_demos_weduc_w = kbal_demos_weduc$w,
                      kpop_all_w = kbal_full$w
  )
  out$out_cat <- data.frame(b = b[i],
                            tolerance = tolerance, 
                            popw = POPW,
                            maxit = maxit,
                            bb_orig = kbal$biasbound.orig,
                            bb_opt =  kbal$biasbound.opt,
                            est = svycontrast(svymean(~recode_vote_2016,
                                                      svydesign(~1, data = pew,
                                                                weights = kbal$w[kbal_data_sampled ==1]),
                                                      na.rm = TRUE), vote_diff)*100,
                            bb_orig_conv = kbal_conv$biasbound.orig,
                            bb_opt_conv = kbal_conv$biasbound.opt,
                            est_conv = svycontrast(svymean(~recode_vote_2016,
                                                           svydesign(~1, data = pew,
                                                                     weights = kbal_conv$w[kbal_data_sampled ==1]),
                                                           na.rm = TRUE), vote_diff)*100,
                            bb_orig_mf = kbal_mf$biasbound.orig,
                            bb_opt_mf = kbal_mf$biasbound.opt,
                            est_mf = svycontrast(svymean(~recode_vote_2016,
                                                         svydesign(~1, data = pew,
                                                                   weights = kbal_mf$w[kbal_data_sampled ==1]),
                                                         na.rm = TRUE), vote_diff)*100,
                            bb_orig_demos =  kbal_demos$biasbound.orig,
                            bb_opt_demos =  kbal_demos$biasbound.opt,
                            est_demos = svycontrast(svymean(~recode_vote_2016,
                                                            svydesign(~1, data = pew,
                                                                      weights = kbal_demos$w[kbal_data_sampled ==1]),
                                                            na.rm = TRUE), vote_diff)*100,
                            bb_orig_demos_weduc =  kbal_demos_weduc$biasbound.orig,
                            bb_opt_demos_weduc =  kbal_demos_weduc$biasbound.opt,
                            est_demos_weduc = svycontrast(svymean(~recode_vote_2016,
                                                                  svydesign(~1, data = pew,
                                                                            weights = kbal_demos_weduc$w[kbal_data_sampled ==1]),
                                                                  na.rm = TRUE), vote_diff)*100,
                            
                            bb_orig_all = kbal_full$biasbound.orig,
                            bb_opt_all =  kbal_full$biasbound.opt,
                            est_all = svycontrast(svymean(~recode_vote_2016,
                                                           svydesign(~1, data = pew,
                                                                     weights = kbal_full$w[kbal_data_sampled ==1]),
                                                           na.rm = TRUE), vote_diff)*100,
                            
                            numdims = ifelse(!is.null(kbal$numdims),kbal$numdims, NA),
                            numdims_conv = ifelse(!is.null(kbal_conv$numdims), 
                                                  kbal_conv$numdims, NA),
                            numdims_mf = ifelse(!is.null(kbal_mf$numdims), 
                                                kbal_mf$numdims, NA),
                            mf_append_dims = ifelse(!is.null(kbal_mf$meanfirst.dims),
                                                    kbal_mf$meanfirst.dims, NA),
                            numdims_demos = ifelse(!is.null(kbal_demos$numdims),
                                                   kbal_demos$numdims, NA),
                            numdims_demos_weduc = ifelse(!is.null(kbal_demos_weduc$numdims),
                                                         kbal_demos_weduc$numdims, NA),
                            numdims_all = ifelse(!is.null(kbal_full$numdims),
                                                  kbal_full$numdims, NA), 
                            l1_orig = ifelse(!is.null(kbal$L1.orig),kbal$L1.orig, NA),
                            l1 = ifelse(!is.null(kbal$L1.opt),kbal$L1.opt, NA),
                            l1_conv = ifelse(!is.null(kbal_conv$L1.opt), 
                                             kbal_conv$L1.opt, NA),
                            l1_mf = ifelse(!is.null(kbal_mf$L1.opt), 
                                           kbal_mf$L1.opt, NA),
                            l1_demos = ifelse(!is.null(kbal_demos$L1.opt),
                                              kbal_demos$L1.opt, NA),
                            l1_demos_weduc = ifelse(!is.null(kbal_demos_weduc$L1.opt),
                                                    kbal_demos_weduc$L1.opt, NA),
                            l1_all = ifelse(!is.null(kbal_full$L1.opt),
                                             kbal_full$L1.opt, NA)
                            
  ) %>% 
    mutate(bbr = bb_orig/bb_opt, 
           bbr_conv = bb_orig_conv/bb_opt_conv,
           bbr_mf = bb_orig_mf/bb_opt_mf,
           bbr_demos = bb_orig_demos/bb_opt_demos,
           bbr_demos_weduc = bb_orig_demos_weduc/bb_opt_demos_weduc,
           bbr_all = bb_orig_all/bb_opt_all
    )
  
  rm(kbal, kbal_demos,kbal_demos_wedu, kbal_full, kbal_mf, kbal_conv)
  return(out)
  
},  mc.cores = length(b))


save(cat_app, tolerance, maxit,
     file = paste0("bmaxvar_nodiag_cat_lasso061021_w", POPW, "_m",maxit, "_t",tolerance, "_",
                   Sys.Date(), ".Rdata"))
out_cat <-lapply(lapply(cat_app, `[`, 2), '[[', 1)
out_cat <- out_cat %>% bind_rows()
weights_cat <- lapply(lapply(cat_app, `[`, 1), '[[', 1)
save(out_cat, weights_cat, tolerance, file = paste0("bmaxvar_cat_lasso061021_w", POPW, "_m",maxit, "_t",tolerance, "_",
                                                    Sys.Date(), ".Rdata") )
print(out_cat)