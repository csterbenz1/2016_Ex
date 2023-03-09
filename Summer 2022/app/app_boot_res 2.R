#processing results of cluster boostrap

### Packages
library(tidyverse)
library(survey)
library(parallel)
#devtools::install_github("csterbenz1/KBAL", ref = "cat_kernel") 
library(kbal)
library(dplyr)

if(detectCores() > 10) {
    path_data= "/home/csterbenz/Data/"
} else {
    path_data= "/Users/Ciara_1/Dropbox/kpop/Updated/application/data/" 
}

POPW = TRUE

############################## Load Data ##########################
## SURVEY DATA (PEW)
pew <- readRDS(paste0(path_data, "pew_lasso_061021.rds"))

## AUXILIARY INFORMATION (CCES)
### Load
cces <- readRDS(paste0(path_data, "cces_lasso_061021.rds"))

cces <- cces%>% 
    mutate(commonweight_vv_post = commonweight_vv_post/ mean(commonweight_vv_post))

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

######### load weights for results
path_weights = "/Users/Ciara_1/Dropbox/kpop/Updated/application/weights/"
load(paste0(path_weights,"bmaxvar_cat_par_wTRUE_m500_t1e-06_2021-03-30.Rdata"))
index = 1
b_num = out_cat$b
la = weights_cat %>% bind_cols()
la


kpop_w <- weights_cat[[index]]$kpop_w
kpop_conv_w <- weights_cat[[index]]$kpop_w_conv
kpop_mf_w <- weights_cat[[index]]$kpop_mf_w
kpop_demos_w <- weights_cat[[index]]$kpop_contraint_w
kpop_demos_wedu_w <- weights_cat[[index]]$kpop_contraint_weduc_w
kpop_all_w <- weights_cat[[index]]$kpop_full_constraint_w


if(POPW) {
    cces_svy <- svydesign(ids = ~1, weights = ~commonweight_vv_post,
                          data = cces)
} else {
    cces_svy <- suppressWarnings(svydesign(ids = ~1,
                                           data = cces))
}

kpop <- svydesign(~1, data = pew,
                  weights = kpop_w[kbal_data_sampled ==1])

kpop_conv <- svydesign(~1, data = pew,
                       weights = kpop_conv_w[kbal_data_sampled ==1])

kpop_mf <-  svydesign(~1, data = pew,
                      weights = kpop_mf_w[kbal_data_sampled ==1])

kpop_demos <-  svydesign(~1, data = pew,
                         weights = kpop_demos_w[kbal_data_sampled ==1])
kpop_demos_wedu <-  svydesign(~1, data = pew,
                              weights = kpop_demos_wedu_w[kbal_data_sampled ==1])

kpop_all <-  svydesign(~1, data = pew,
                       weights = kpop_all_w[kbal_data_sampled ==1])

############### estimates
vote_diff <- quote((recode_vote_2016Democrat - recode_vote_2016Republican) /(recode_vote_2016Democrat + recode_vote_2016Republican + recode_vote_2016Other)) 
kpop_real = svycontrast(svymean(~recode_vote_2016, kpop,  na.rm = TRUE), vote_diff)

var_pew = "diff_cces_on_pew"
kpop = svymean(as.formula(paste0("~", var_pew)), kpop, 
                na.rm = TRUE)
kpop_conv = svymean(as.formula(paste0("~", var_pew)), kpop_conv,
                    na.rm = TRUE)
kpop_amf = svymean(as.formula(paste0("~", var_pew)), kpop_mf, 
                   na.rm = TRUE)
kpop_demos_mf =svymean(as.formula(paste0("~", var_pew)), kpop_demos,
                       na.rm = TRUE)
kpop_demos_wedu_mf =svymean(as.formula(paste0("~", var_pew)), 
                            kpop_demos_wedu, na.rm = TRUE)
kpop_all_mf = svymean(as.formula(paste0("~", var_pew)), 
                      kpop_all, na.rm = TRUE)

#matches Appendix table C.3
kpop
kpop_conv
kpop_demos_mf
kpop_demos_wedu_mf
kpop_all_mf

#### but what is this newer file with no diag???
load(paste0(path_weights,"bmaxvar_cat_nodiag_lasso061021_wTRUE_m500_t1e-06_2021-07-07.Rdata"))
#this file us made via app_script_bmaxvar_only.R
#this does not replicate to the same results, though all that should have changed is the best b
la = weights_cat %>% bind_cols()
la
#but it's the same so wtf?


################# ANYWAY: Bootstrap
load("./kbal_boot_popwTRUE_m500_t1e-06_incby5_2022-07-28_nbootreps252.Rdata")
colnames(estimates) <- c("boot_real_est", "boot_real_SE", "boot_diff_est", "boot_diff_SE")
colMeans(estimates)[c(1,3)]
kpop*100
kpop_real*100

###### variance of estimates
apply(estimates, 2, sd)[c(1,3)]

var_w = weights %>% group_by(iter) %>% summarise(var_w = var(boot_w))
var_w
mean(var_w$var_w)
summary(var_w$var_w)
sd(var_w$var_w)
