## WEIGHTING A 2016 PRE-ELECTION SURVEY

## SETUP

### Packages
library(tidyverse)
library(survey)

#to be absolutely sure this is installing the right version:
#these should be the same rn
#devtools::install_github("csterbenz1/KBAL", ref = "pop_weights") 
#devtools::install_github("chadhazlett/KBAL", ref = "master") 
library(kbal)
#setwd("/Users/Ciara/Documents/Cloud Documents/Hazlett:Hartman RA/Code/2016_reweighting_example/")



############################## NON KBAL METHODS #######################

### Function for creating targets from auxiliary information and formula
create_targets <- function (target_design, target_formula) {
    target_mf <- model.frame(target_formula, model.frame(target_design))
    target_mm <- model.matrix(target_formula, target_mf)
    wts <- weights(target_design)
    colSums(target_mm * wts) / sum(wts)
}



## SURVEY DATA (PEW)
### Load
#This needs to be loaded from dropbox/local not github bc it's too big
pew <- readRDS("../data/pew.rds")
### Make survey design 
pew_srs <- svydesign(ids = ~1, data = pew)


### Unweighted survey estimates of presidential vote
#function to get vote margin
vote_contrast <- quote((recode_vote_2016Democrat - recode_vote_2016Republican) /
                           (recode_vote_2016Democrat + recode_vote_2016Republican))
#svymean(~recode_vote_2016, design = pew_srs)
#svycontrast(svymean(~recode_vote_2016, pew_srs), vote_contrast)

### Auxiliary variables balance
# svymean(~recode_female, design = pew_srs)
# svymean(~recode_age_bucket, design = pew_srs)
# svymean(~recode_race, design = pew_srs)
# svymean(~recode_region, design = pew_srs)
# svymean(~recode_educ, design = pew_srs)

## AUXILIARY INFORMATION (CCES)
### Load
#This needs to be loaded from dropbox/local not github bc it's too big
cces <- readRDS("../data/cces.rds")
### Drop invalid cases
cces <- cces %>%
    filter((CC16_401 == "I definitely voted in the General Election.") &
               !is.na(commonweight_vv_post))

### Make survey design
cces_awt <- svydesign(ids = ~1, weights = ~commonweight_vv_post, data = cces)

### unweighted CCES
cces_nowt <- svydesign(ids = ~1, data = cces)

### Presidential vote estimates from weighted CCES
#### National
# svymean(~recode_vote_2016, design = cces_awt, na.rm = TRUE)
# svycontrast(svymean(~recode_vote_2016, cces_awt, na.rm = TRUE), vote_contrast)
# #### State
# svyby(~I(as.numeric(recode_vote_2016 == "Democrat")), ~recode_inputstate,
#       design = cces_awt, svymean, na.rm = TRUE, keep.var = FALSE)
# 
# ### Auxiliary variables balance in cces weighted
# svymean(~recode_female, design = cces_awt)
# svymean(~recode_age_bucket, design = cces_awt)
# svymean(~recode_race, design = cces_awt)
# svymean(~recode_region, design = cces_awt)
# svymean(~recode_educ, design = cces_awt)


###### POPULATION TARGETS
### Formulas for auxiliary vector

# NO PID 
#### (1) Marginal distributions of age, female, race, and region
formula_1 <- ~recode_age_bucket + recode_female + recode_race + recode_region 

#### (2) Marginal distributions of age, female, race, region, and education
formula_2 <- ~recode_age_bucket + recode_female + recode_race + recode_region + 
    recode_educ 

#### (3) Joint distribution of age (coarsened), female, race, region, and
#### education (coarsened)
formula_3 <- ~recode_age_3way * recode_female * recode_race *
    recode_region * recode_educ_3way

#### (4) Marginal distributions of age, female, and race and, among whites,
#### the joint distribution of region and education
#### (recode_race_educ_reg = race * educ * reg if race == "white" and
#### race * reg otherwise)
formula_4 <- ~recode_age_bucket + recode_female + recode_race_educ_reg


# WITH PID
formula_1_pid <- ~recode_age_bucket + recode_female + recode_race + recode_region + recode_pid_3way
formula_2_pid <- ~recode_age_bucket + recode_female + recode_race + recode_region + 
    recode_educ + recode_pid_3way
formula_3_pid <- ~recode_age_3way * recode_female * recode_race *
    recode_region * recode_educ_3way * recode_pid_3way
formula_4_pid <- ~recode_age_bucket + recode_female + recode_race_educ_reg + recode_pid_3way
#add interaction with pid
formula_5_pid <- ~recode_age_bucket:recode_pid_3way + recode_female:recode_pid_3way+
    recode_race_educ_reg:recode_pid_3way
#sum(svytable(formula_5_pid, pew_srs) == 0) 

### Population targets
targets_1 <- create_targets(cces_awt, formula_1)
targets_2 <- create_targets(cces_awt, formula_2)
targets_3 <- create_targets(cces_awt, formula_3) # will have to modify below
targets_4 <- create_targets(cces_awt, formula_4)

targets_1_pid <- create_targets(cces_awt, formula_1_pid)
targets_2_pid <- create_targets(cces_awt, formula_2_pid)
targets_3_pid <- create_targets(cces_awt, formula_3_pid) # will have to modify below
targets_4_pid <- create_targets(cces_awt, formula_4_pid)

targets_5_pid <- create_targets(cces_awt, formula_5_pid)

#### Weighted survey designs
# NO PID
#### (1)
pew_lwt_1 <- calibrate(design = pew_srs,
                       formula = formula_1,
                       population = targets_1,
                       calfun = "raking")

#### (2)
pew_lwt_2 <- calibrate(design = pew_srs,
                       formula = formula_2,
                       population = targets_2,
                       calfun = "raking")

#### (3)
##### Can't compute below because some cells are empty.
try(pew_lwt_3 <- calibrate(design = pew_srs,
                           formula = formula_3,
                           population = targets_3,
                           calfun = "linear"),
    silent = T)
##### So instead we use the `postStratify` function with `partial = TRUE`, which
##### ignores empty cells.
formula_3_ps <- as.formula(str_replace_all(formula_3, "\\*", "+"))
targets_3_ps <- svytable(formula = formula_3_ps, design = cces_awt)
sum(svytable(formula_3_ps, pew_srs) == 0) # 265 empty cells
pew_ps_3 <- postStratify(design = pew_srs,
                         strata = formula_3_ps,
                         population = targets_3_ps,
                         partial = TRUE) 

#### (4) our "post hoc" best estimate
pew_lwt_4 <- calibrate(design = pew_srs,
                       formula = formula_4,
                       population = targets_4,
                       calfun = "raking")



## W PID
#### (1)
pew_lwt_1_pid <- calibrate(design = pew_srs,
                       formula = formula_1_pid,
                       population = targets_1_pid,
                       calfun = "raking")

#### (2)
pew_lwt_2_pid <- calibrate(design = pew_srs,
                       formula = formula_2_pid,
                       population = targets_2_pid,
                       calfun = "raking")

#### (3)
##### Still can't compute below because some cells are empty.
try(pew_lwt_3_pid <- calibrate(design = pew_srs,
                           formula = formula_3_pid,
                           population = targets_3_pid,
                           calfun = "linear"),
    silent = T)
#Try using postStratify() with partial = TRUE to ignore empty cells but we end up with 56 NA weights
formula_3_ps_pid <- as.formula(str_replace_all(formula_3_pid, "\\*", "+"))
targets_3_ps_pid <- svytable(formula = formula_3_ps_pid, design = cces_awt)
sum(svytable(formula_3_ps_pid, pew_srs) == 0) # 265 empty cells -> 870 with pid
pew_ps_3_pid <- postStratify(design = pew_srs,
                         strata = formula_3_ps_pid,
                         population = targets_3_ps_pid,
                         partial = TRUE)

#this gives 56 na's 
sum(is.na((weights(pew_ps_3_pid ))))

#### (4) our "post-hoc" 
pew_lwt_4_pid <- calibrate(design = pew_srs,
                       formula = formula_4_pid,
                       population = targets_4_pid,
                       calfun = "raking")



###(5) not working with pid interacted with the huge interaction
pew_lwt_5_pid <- calibrate(design = pew_srs,
                           formula = formula_5_pid,
                           population = targets_5_pid,
                           calfun = "raking")
#how many empty cells? riiight need to make this target in cces too
formula_5_ps_pid <- as.formula(str_replace_all(formula_5_pid, "\\:", "+"))
targets_5_ps_pid <- svytable(formula = formula_5_ps_pid, design = cces_awt)
sum(svytable(formula_5_ps_pid, pew_srs) == 0)  #325 :()
############################## KBAL METHODS ##########################


##### KBAL DATA PREP
#filling in NAs in age with mean est. in cces unweighted
#sum(is.na(pew$recode_age)) #37 missing
#sum(is.na(cces$recode_age)) #none missing
avg_age <- mean(cces$recode_age, na.rm = TRUE)
pew$recode_age[is.na(pew$recode_age)] <- avg_age
cces$recode_age[is.na(cces$recode_age)] <- avg_age


# W PID - including all vars used in pew4 + 3way pid
kbal_data <- rbind(pew %>% select(recode_age,
                                  recode_female,
                                  recode_race, 
                                  recode_region, 
                                  recode_pid_3way,
                                  recode_educ, 
                                  recode_age_bucket),
                   cces %>% select(recode_age, 
                                   recode_female, 
                                   recode_race, 
                                   recode_region,
                                   recode_pid_3way, 
                                   recode_educ,
                                   recode_age_bucket)) %>%
  model.matrix(as.formula("~. - 1"), .)

# NO PID: selecting all vars in pew4 + NO PID
kbal_data_nopid <- rbind(pew %>% select(recode_age, 
                                        recode_female, 
                                        recode_race, 
                                        recode_region, 
                                        recode_educ,
                                        recode_age_bucket),
                         cces %>% select(recode_age,
                                         recode_female,
                                         recode_race, 
                                         recode_region,
                                         recode_educ,
                                         recode_age_bucket)) %>%
    model.matrix(as.formula("~. - 1"), .)

kbal_data_sampled <- c(rep(1, nrow(pew)), rep(0, nrow(cces)))

## NOTE: the weights for the cces are in cces$commonweight_vv_post

 ########################### Load KBal Runs ###############
#these are enormous so try to avoid loading if possible
#this is needed for choice of b comparisons, and to get numdims out for plots
#but these are saved in bb_dat.Rdata and numdims_XXXX.Rdata
#can use weights alone and survey designs and numdims which are much smaller
#so shouldnt need to load these huge objects
#loading old runs
#load("mf_T_wPid_varyb.Rdata")
#load("mf_F_wPid_varyb.Rdata")

#load("mf_T_NoPid_varyb.Rdata")
#load("mf_F_NoPid_varyb.Rdata")
#load("frankenstein_erin.Rdata")

################## KPOP + W PID #########################
####### KPOP + MF = T + W PID

#these runs saved here: mf_T_wPid_varyb.Rdata 
### KPOP + MF=T + PID + b = 2 + Conv
# tictoc::tic()
# kbal_mf_b2x_f <- kbal(allx=kbal_data,
#                       sampled = kbal_data_sampled,
#                       b = 2 * ncol(kbal_data),
#                       incrementby = 1,
#                       fullSVD = TRUE,
#                       meanfirst = TRUE,
#                       population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
#                       sampledinpop = FALSE,
#                       ebal.convergence = TRUE)
# tictoc::toc()
# 
# 
# # 
# # ### KPOP + MF=T + PID + b = 1
# kbal_mf_b1x_f <- kbal(allx=kbal_data,
#                       sampled = kbal_data_sampled,
#                       b = 1 * ncol(kbal_data),
#                       incrementby = 1,
#                       fullSVD = TRUE,
#                       meanfirst = TRUE,
#                       population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
#                       sampledinpop = FALSE,
#                       ebal.convergence = TRUE)
# # 
# # 
# # 
# # ### KPOP + MF=T + PID + b = .5
# kbal_mf_b.5x_f <- kbal(allx=kbal_data,
#                        sampled = kbal_data_sampled,
#                        b = 0.5 * ncol(kbal_data),
#                        fullSVD = TRUE,
#                        meanfirst = TRUE,
#                        incrementby = 1,
#                        population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
#                        sampledinpop = FALSE,
#                        ebal.convergence = TRUE)
# 
# 
# ### KPOP + MF=T + PID + b = .25 
# kbal_mf_b.25x_f <- kbal(allx=kbal_data,
#                         sampled = kbal_data_sampled,
#                         b = 0.25 * ncol(kbal_data),
#                         incrementby = 1,
#                         fullSVD = TRUE,
#                         meanfirst = TRUE,
#                         population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
#                         sampledinpop = FALSE,
#                         ebal.convergence = TRUE)
# # 
# 
# ### KPOP + MF=T + PID + b = .125
# kbal_mf_b.125x_f <- kbal(allx=kbal_data,
#                          sampled = kbal_data_sampled,
#                          b = 0.125 * ncol(kbal_data),
#                          incrementby = 1,
#                          fullSVD = TRUE,
#                          meanfirst = TRUE,
#                          population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
#                          sampledinpop = FALSE,
#                          ebal.convergence = TRUE)
# 
# save(kbal_mf_b.125x_f, kbal_mf_b.25x_f, kbal_mf_b.5x_f,kbal_mf_b1x_f,kbal_mf_b2x_f,
#      file = "mf_T_wPID_f.Rdata")
# 
# ############### KPOP + MF = F + PID (No Conv)
# #these runs saved here: mf_F_wPid_varyb.Rdata
# 
# kbal_b2x_f <- kbal(allx=kbal_data,
#                    sampled = kbal_data_sampled,
#                    b = 2 * ncol(kbal_data),
#                    incrementby = 1,
#                    fullSVD = TRUE,
#                    meanfirst = FALSE,
#                    population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
#                    sampledinpop = FALSE,
#                    ebal.convergence = TRUE)
# 
# 
# # 
# # ### KPOP + MF=T + PID + b = 1
# kbal_b1x_f <- kbal(allx=kbal_data,
#                    sampled = kbal_data_sampled,
#                    b = 1 * ncol(kbal_data),
#                    incrementby = 1,
#                    fullSVD = TRUE,
#                    meanfirst = FALSE,
#                    population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
#                    sampledinpop = FALSE,
#                    ebal.convergence = TRUE)
# # 
# # 
# # 
# # ### KPOP + MF=T + PID + b = .5
# kbal_b.5x_f <- kbal(allx=kbal_data,
#                     sampled = kbal_data_sampled,
#                     b = 0.5 * ncol(kbal_data),
#                     fullSVD = TRUE,
#                     meanfirst = FALSE,
#                     incrementby = 1,
#                     population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
#                     sampledinpop = FALSE,
#                     ebal.convergence = TRUE)
# 
# 
# ### KPOP + MF=T + PID + b = .25 
# kbal_b.25x_f <- kbal(allx=kbal_data,
#                      sampled = kbal_data_sampled,
#                      b = 0.25 * ncol(kbal_data),
#                      incrementby = 1,
#                      fullSVD = TRUE,
#                      meanfirst = FALSE,
#                      population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
#                      sampledinpop = FALSE,
#                      ebal.convergence = TRUE)
# # 
# 
# ### KPOP + MF=T + PID + b = .125
# kbal_b.125x_f <- kbal(allx=kbal_data,
#                       sampled = kbal_data_sampled,
#                       b = 0.125 * ncol(kbal_data),
#                       incrementby = 1,
#                       fullSVD = TRUE,
#                       meanfirst = FALSE,
#                       population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
#                       sampledinpop = FALSE,
#                       ebal.convergence = TRUE)
# 
# save(kbal_b.125x_f,kbal_b.25x_f,kbal_b.5x_f,kbal_b1x_f,kbal_b2x_f, 
#      file = "mf_F_wPid_f.Rdata")
# 
# #####test
# 
# ################## KPOP + NO PID #########################
# 
# ######## KPOP + MF=F + NO PID + No Conv
# 
# #these runs are saved in: mf_F_NoPid_varyb.Rdata
# 
# ##### Kbal + NO PID + MF = FALSE + b = 2 + No Conv 
# kbalout_nopid_b2x_f <- kbal(allx=kbal_data_nopid,
#                             sampled = kbal_data_sampled,
#                             b = 2 * ncol(kbal_data_nopid),
#                             incrementby = 1,
#                             fullSVD = TRUE,
#                             population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
#                             sampledinpop = FALSE)
# ##### Kbal + NO PID + MF = FALSE + b = 1 + No Conv 
# kbalout_nopid_b1x_f <- kbal(allx=kbal_data_nopid,
#                             sampled = kbal_data_sampled,
#                             b = 1 * ncol(kbal_data_nopid),
#                             incrementby = 1,
#                             fullSVD = TRUE,
#                             population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
#                             sampledinpop = FALSE)
# 
# #### (2) Kbal + NO PID + MF = FALSE + b=0.5x + NoConv
# kbalout_nopid_b.5x_f <- kbal(allx=kbal_data_nopid,
#                              sampled = kbal_data_sampled,
#                              b = 0.5 * ncol(kbal_data_nopid),
#                              incrementby = 1,
#                              fullSVD = TRUE,
#                              population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
#                              sampledinpop = FALSE)
# 
# #### (1) Kbal + NO PID + MF=FALSE + b=0.25x + NoConv
# kbalout_nopid_b.25x_f <- kbal(allx=kbal_data_nopid,
#                               sampled = kbal_data_sampled,
#                               b = 0.25 * ncol(kbal_data_nopid),
#                               incrementby = 1,
#                               fullSVD = TRUE,
#                               population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
#                               sampledinpop = FALSE)
# 
# #### (3) Kbal + NO PID + MF = FALSE + b=0.125x + NoConv
# kbalout_nopid_b.125x_f <- kbal(allx=kbal_data_nopid,
#                                sampled = kbal_data_sampled,
#                                b = 0.125 * ncol(kbal_data_nopid),
#                                incrementby = 1,
#                                fullSVD = TRUE,
#                                population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
#                                sampledinpop = FALSE)
# 
# #### (4) Kbal + NO PID + MF = FALSE + b=0.0625x + NoConv
# # kbalout_nopid_b.0625x_f <- kbal(allx=kbal_data_nopid,
# #                                 sampled = kbal_data_sampled,
# #                                 b = 0.0625 * ncol(kbal_data_nopid),
# #                                 incrementby = 1,
# #                                 fullSVD = TRUE,
# #                 population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
# #                 sampledinpop = FALSE)
# 
# save( kbalout_nopid_b.125x_f,kbalout_nopid_b.25x_f,
#       kbalout_nopid_b.5x_f, kbalout_nopid_b1x_f, kbalout_nopid_b2x_f,
#       file = "mf_F_NoPid_f.Rdata")
# 
# ############ KPOP + MF=TRUE + NO PID +  Conv
# 
# #these runs stored in: mf_T_NoPid_varyb.Rdata
# #### (1) Kbal + NO PID + MF = TRUE + b=2x  + Conv
# kbalout_mf_nopid_b2x_f <- kbal(allx=kbal_data_nopid,
#                                sampled = kbal_data_sampled,
#                                b = 2 * ncol(kbal_data_nopid),
#                                meanfirst = TRUE,
#                                fullSVD = TRUE,
#                                population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
#                                sampledinpop = FALSE,
#                                ebal.convergence = TRUE)
# 
# #### (2) Kbal + NO PID + MF = TRUE + b=1x  + Conv
# kbalout_mf_nopid_b1x_f <- kbal(allx=kbal_data_nopid,
#                                sampled = kbal_data_sampled,
#                                b = 1 * ncol(kbal_data_nopid),
#                                meanfirst = TRUE,
#                                fullSVD = TRUE,
#                                population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
#                                sampledinpop = FALSE,
#                                ebal.convergence = TRUE)
# 
# #### (3) Kbal + NO PID + MF = TRUE + b=0.5x  + Conv
# kbalout_mf_nopid_b.5x_f <- kbal(allx=kbal_data_nopid,
#                                 sampled = kbal_data_sampled,
#                                 b = 0.5 * ncol(kbal_data_nopid),
#                                 meanfirst = TRUE,
#                                 fullSVD = TRUE,
#                                 population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
#                                 sampledinpop = FALSE,
#                                 ebal.convergence = TRUE)
# 
# #### (4)
# kbalout_mf_nopid_b.25x_f <- kbal(allx=kbal_data_nopid,
#                                  sampled = kbal_data_sampled,
#                                  b = 0.25 * ncol(kbal_data_nopid),
#                                  meanfirst = TRUE,
#                                  fullSVD = TRUE,
#                                  population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
#                                  sampledinpop = FALSE,
#                                  ebal.convergence = TRUE)
# 
# #### (5) Kbal + NO PID + MF = TRUE + b=0.125x  + Conv
# kbalout_mf_nopid_b.125x_f <- kbal(allx=kbal_data_nopid,
#                                   sampled = kbal_data_sampled,
#                                   b = 0.125 * ncol(kbal_data_nopid),
#                                   meanfirst = TRUE,
#                                   fullSVD = TRUE,
#                                   population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
#                                   sampledinpop = FALSE,
#                                   ebal.convergence = TRUE)
# 
# save(kbalout_mf_nopid_b2x_f, kbalout_mf_nopid_b1x_f, kbalout_mf_nopid_b.5x_f,
#      kbalout_mf_nopid_b.25x_f, kbalout_mf_nopid_b.125x_f,
#      file = "mf_T_NoPid_f.Rdata")
# 
# ############################## Choice of B ############################
# 
# 
# 
# #### WITH PID
# 
# ##### (1) Kbal + MF = FALSE + PID
# bb_comp <- data.frame(biasbound_orig = c(kbal_b2x_f$biasbound.orig,
#                                          kbal_b1x_f$biasbound.orig,
#                                          kbal_b.5x_f$biasbound.orig,
#                                          kbal_b.25x_f$biasbound.orig,
#                                          kbal_b.125x_f$biasbound.orig),
#                       biasbound_opt = c(kbal_b2x_f$biasbound.opt,
#                                         kbal_b1x_f$biasbound.opt,
#                                         kbal_b.5x_f$biasbound.opt,
#                                         kbal_b.25x_f$biasbound.opt,
#                                         kbal_b.125x_f$biasbound.opt),
#                       L1_orig = c(kbal_b2x_f$L1.orig,
#                                   kbal_b1x_f$L1.orig,
#                                   kbal_b.5x_f$L1.orig,
#                                   kbal_b.25x_f$L1.orig,
#                                   kbal_b.125x_f$L1.orig),
#                       L1_opt = c(kbal_b2x_f$L1.opt,
#                                  kbal_b1x_f$L1.opt,
#                                  kbal_b.5x_f$L1.opt,
#                                  kbal_b.25x_f$L1.opt,
#                                  kbal_b.125x_f$L1.opt)
# )
# 
# bb_comp <- bb_comp %>% mutate(bb_ratio = biasbound_orig/biasbound_opt,
#                               L1_ratio = L1_orig/L1_opt) %>%
#     select(bb_ratio,L1_ratio, biasbound_orig, biasbound_opt,  L1_orig, L1_opt)
# 
# rownames(bb_comp) <- c("b=2x",
#                        "b=1x",
#                        "b=.5x",
#                        "b=.25x",
#                        "b=.125x")
# 
# # kable(bb_comp,
# #       format = "latex",
# #       caption = "KPOP + MF=FALSE + W PID: Comparison of Bias bound and L1 distance by choice of b",
# #       col.names = c("Bias Bound Ratio", "L1 Ratio",
# #                     "Bias Bound Orig", "Bias Bound Opt",
# #                     "L1 Orig", "L1 Opt"),
# #       booktabs = T, digits= 4) %>%
# #     kable_styling(position = "center", latex_options = "hold_position")
# 
# ##### (2) Kbal + MF = TRUE + PID
# bb_comp_mf <- data.frame(biasbound_orig = c(kbal_mf_b2x_f$biasbound.orig,
#                                             kbal_mf_b1x_f$biasbound.orig,
#                                             kbal_mf_b.5x_f$biasbound.orig,
#                                             kbal_mf_b.25x_f$biasbound.orig,
#                                             kbal_mf_b.125x_f$biasbound.orig),
#                          biasbound_opt = c(kbal_mf_b2x_f$biasbound.opt,
#                                            kbal_mf_b1x_f$biasbound.opt,
#                                            kbal_mf_b.5x_f$biasbound.opt,
#                                            kbal_mf_b.25x_f$biasbound.opt,
#                                            kbal_mf_b.125x_f$biasbound.opt),
#                          L1_orig = c(kbal_mf_b2x_f$L1.orig,
#                                      kbal_mf_b1x_f$L1.orig,
#                                      kbal_mf_b.5x_f$L1.orig,
#                                      kbal_mf_b.25x_f$L1.orig,
#                                      kbal_mf_b.125x_f$L1.orig),
#                          L1_opt = c(kbal_mf_b2x_f$L1.opt,
#                                     kbal_mf_b1x_f$L1.opt,
#                                     kbal_mf_b.5x_f$L1.opt,
#                                     kbal_mf_b.25x_f$L1.opt,
#                                     kbal_mf_b.125x_f$L1.opt)
# )
# 
# bb_comp_mf <- bb_comp_mf %>% mutate(bb_ratio = biasbound_orig/biasbound_opt,
#                                     L1_ratio = L1_orig/L1_opt) %>%
#     select(bb_ratio,L1_ratio, biasbound_orig, biasbound_opt,  L1_orig, L1_opt)
# 
# rownames(bb_comp_mf) <- c("b=2x",
#                           "b=1x",
#                           "b=.5x",
#                           "b=.25x",
#                           "b=.125x"
# )
# 
# 
# # library(kableExtra)
# # kable(bb_comp_mf,
# #       format = "latex",
# #       caption = "KPOP + MF=TRUE + W PID: Comparison of Bias bound and L1 distance by choice of b",
# #       col.names = c("Bias Bound Ratio", "L1 Ratio",
# #                     "Bias Bound Orig", "Bias Bound Opt",
# #                     "L1 Orig", "L1 Opt"),
# #       booktabs = T, digits= 4) %>%
# #     kable_styling(position = "center", latex_options = "hold_position")
# 
# ########## NO PID ###############
# 
# ##### (3) Kbal + MF = FALSE + NO PID
# bb_comp_nopid <- data.frame(biasbound_orig = c(kbalout_nopid_b2x_f$biasbound.orig,
#                                                kbalout_nopid_b1x_f$biasbound.orig,
#                                                kbalout_nopid_b.5x_f$biasbound.orig,
#                                                kbalout_nopid_b.25x_f$biasbound.orig,
#                                                kbalout_nopid_b.125x_f$biasbound.orig),
#                             #kbalout_nopid_b.0625x$biasbound.orig),
#                             biasbound_opt = c(kbalout_nopid_b2x_f$biasbound.opt,
#                                               kbalout_nopid_b1x_f$biasbound.opt,
#                                               kbalout_nopid_b.5x_f$biasbound.opt,
#                                               kbalout_nopid_b.25x_f$biasbound.opt,
#                                               kbalout_nopid_b.125x_f$biasbound.opt),
#                             # kbalout_nopid_b.0625x$biasbound.opt),
#                             L1_orig = c(kbalout_nopid_b2x_f$L1.orig,
#                                         kbalout_nopid_b1x_f$L1.orig,
#                                         kbalout_nopid_b.5x_f$L1.orig,
#                                         kbalout_nopid_b.25x_f$L1.orig,
#                                         kbalout_nopid_b.125x_f$L1.orig),
#                             # kbalout_nopid_b.0625x$L1.orig),
#                             L1_opt = c(kbalout_nopid_b2x_f$L1.opt,
#                                        kbalout_nopid_b1x_f$L1.opt,
#                                        kbalout_nopid_b.5x_f$L1.opt,
#                                        kbalout_nopid_b.25x_f$L1.opt,
#                                        kbalout_nopid_b.125x_f$L1.opt)
#                             #kbalout_nopid_b.0625x$L1.opt)
# )
# 
# bb_comp_nopid <- bb_comp_nopid %>% mutate(bb_ratio = biasbound_orig/biasbound_opt,
#                                           L1_ratio = L1_orig/L1_opt) %>%
#     select(bb_ratio,L1_ratio, biasbound_orig, biasbound_opt,  L1_orig, L1_opt)
# 
# rownames(bb_comp_nopid) <- c("b=2x",
#                              "b=1x",
#                              "b=.5x",
#                              "b=.25x",
#                              "b=.125x")
# #"b=.0625x")
# save(bb_comp_nopid, bb_comp, bb_comp_mf, file = "bb_dat_full.Rdata")
# # kable(bb_comp_nopid,
# #       format = "latex",
# #       caption = "KPOP + MF=FALSE + NO PID: Comparison of Bias bound and L1 distance by choice of b",
# #       col.names = c("Bias Bound Ratio", "L1 Ratio",
# #                     "Bias Bound Orig", "Bias Bound Opt",
# #                     "L1 Orig", "L1 Opt"),
# #       booktabs = T, digits= 4) %>%
# #     kable_styling(position = "center", latex_options = "hold_position")
# 
# 
# 
# ##### (4) Kbal + MF = TRUE + NO PID: only ran for b.25 and b.125 so far
# bb_comp_mf_nopid <- data.frame(biasbound_orig = c(kbalout_mf_nopid_b2x_f$biasbound.orig,
#                                                   kbalout_mf_nopid_b1x_f$biasbound.orig,
#                                                   kbalout_mf_nopid_b.5x_f$biasbound.orig,
#                                                   kbalout_mf_nopid_b.25x_f$biasbound.orig,
#                                                   kbalout_mf_nopid_b.125x_f$biasbound.orig),
#                                biasbound_opt = c(kbalout_mf_nopid_b2x_f$biasbound.opt,
#                                                  kbalout_mf_nopid_b1x_f$biasbound.opt,
#                                                  kbalout_mf_nopid_b.5x_f$biasbound.opt,
#                                                  kbalout_mf_nopid_b.25x_f$biasbound.opt,
#                                                  kbalout_mf_nopid_b.125x_f$biasbound.opt),
#                                L1_orig = c(kbalout_mf_nopid_b2x_f$L1.orig,
#                                            kbalout_mf_nopid_b1x_f$L1.orig,
#                                            kbalout_mf_nopid_b.5x_f$L1.orig,
#                                            kbalout_mf_nopid_b.25x_f$L1.orig,
#                                            kbalout_mf_nopid_b.125x_f$L1.orig),
#                                L1_opt = c(kbalout_mf_nopid_b2x_f$L1.opt,
#                                           kbalout_mf_nopid_b1x_f$L1.opt,
#                                           kbalout_mf_nopid_b.5x_f$L1.opt,
#                                           kbalout_mf_nopid_b.25x_f$L1.opt,
#                                           kbalout_mf_nopid_b.125x_f$L1.opt)
# )
# # 
# bb_comp_mf_nopid <- bb_comp_mf_nopid %>% mutate(bb_ratio = biasbound_orig/biasbound_opt,
#                                                 L1_ratio = L1_orig/L1_opt) %>%
#     select(bb_ratio,L1_ratio, biasbound_orig, biasbound_opt,  L1_orig, L1_opt)
# rownames(bb_comp_mf_nopid) <- c("b=2x",
#                                 "b=1x",
#                                 "b=.5x",
#                                 "b=.25x",
#                                 "b=.125x")
# 
# save(bb_comp_nopid, bb_comp_mf_nopid, bb_comp, bb_comp_mf, file = "bb_dat_full.Rdata")
# # kable(bb_comp_mf_nopid,
# #       format = "latex",
# #       caption = "KPOP + MF=TRUE + NO PID: Comparison of Bias bound and L1 distance by choice of b",
# #       col.names = c("Bias Bound Ratio", "L1 Ratio",
# #                     "Bias Bound Orig", "Bias Bound Opt",
# #                     "L1 Orig", "L1 Opt"),
# #       booktabs = T, digits= 4) %>%
# #     kable_styling(position = "center", latex_options = "hold_position")
# 
# #save(bb_comp, bb_comp_mf, bb_comp_nopid, bb_comp_mf_nopid,file = "bb_dat.Rdata")
# 
# 
# ############################# Save all but K and svdK ##############
# kbal_b2x_s <- kbal_b2x_f[!(names(kbal_b2x_f) %in% c("K", "svdK"))]
# kbal_b1x_s <- kbal_b1x_f[!(names(kbal_b1x_f) %in% c("K", "svdK"))]
# kbal_b.5x_s <- kbal_b.5x_f[!(names(kbal_b.5x_f) %in% c("K", "svdK"))]
# kbal_b.25x_s <- kbal_b.25x_f[!(names(kbal_b.25x_f) %in% c("K", "svdK"))]
# kbal_b.125x_s <- kbal_b.125x_f[!(names(kbal_b.125x_f) %in% c("K", "svdK"))]
# 
# 
# kbal_mf_b2x_s <- kbal_mf_b2x_f[!(names(kbal_mf_b2x_f) %in% c("K", "svdK"))]
# kbal_mf_b1x_s <- kbal_mf_b1x_f[!(names(kbal_mf_b1x_f) %in% c("K", "svdK"))]
# kbal_mf_b.5x_s <- kbal_mf_b.5x_f[!(names(kbal_mf_b.5x_f) %in% c("K", "svdK"))]
# kbal_mf_b.25x_s <- kbal_mf_b.25x_f[!(names(kbal_mf_b.25x_f) %in% c("K", "svdK"))]
# kbal_mf_b.125x_s <- kbal_mf_b.125x_f[!(names(kbal_mf_b.125x_f) %in% c("K", "svdK"))]
# 
# 
# kbalout_nopid_b2x_s <- kbalout_nopid_b2x_f[!(names(kbalout_nopid_b2x_f) %in% c("K", "svdK"))]
# kbalout_nopid_b1x_s <- kbalout_nopid_b1x_f[!(names(kbalout_nopid_b1x_f) %in% c("K", "svdK"))]
# kbalout_nopid_b.5x_s <- kbalout_nopid_b.5x_f[!(names(kbalout_nopid_b.5x_f) %in% c("K", "svdK"))]
# kbalout_nopid_b.25x_s <- kbalout_nopid_b.25x_f[!(names(kbalout_nopid_b.25x_f) %in% c("K", "svdK"))]
# kbalout_nopid_b.125x_s <- kbalout_nopid_b.125x_f[!(names(kbalout_nopid_b.125x_f) %in% c("K", "svdK"))]
# 
# kbalout_mf_nopid_b2x_s <- kbalout_mf_nopid_b2x_f[!(names(kbalout_mf_nopid_b2x_f) %in% c("K", "svdK"))]
# kbalout_mf_nopid_b1x_s <- kbalout_mf_nopid_b1x_f[!(names(kbalout_mf_nopid_b1x_f) %in% c("K", "svdK"))]
# kbalout_mf_nopid_b.5x_s <- kbalout_mf_nopid_b.5x_f[!(names(kbalout_mf_nopid_b.5x_f) %in% c("K", "svdK"))]
# kbalout_mf_nopid_b.25x_s <- kbalout_mf_nopid_b.25x_f[!(names(kbalout_mf_nopid_b.25x_f) %in% c("K", "svdK"))]
# kbalout_mf_nopid_b.125x_s <- kbalout_mf_nopid_b.125x_f[!(names(kbalout_mf_nopid_b.125x_f) %in% c("K", "svdK"))]
# 
# 
# save( kbal_b2x_s,kbal_b1x_s,kbal_b.5x_s,kbal_b.25x_s,kbal_b.125x_s,
#       kbal_mf_b2x_s,kbal_mf_b1x_s,kbal_mf_b.5x_s, kbal_mf_b.25x_s, kbal_mf_b.125x_s,
#       kbalout_nopid_b2x_s,kbalout_nopid_b1x_s, kbalout_nopid_b.5x_s, kbalout_nopid_b.25x_s, kbalout_nopid_b.125x_s, 
#       kbalout_mf_nopid_b2x_s, kbalout_mf_nopid_b1x_s, kbalout_mf_nopid_b.5x_s, kbalout_mf_nopid_b.25x_s,
#       kbalout_mf_nopid_b.125x_s,
#       file = "fullruns_noKsvd.Rdata")
# 
# ######################################### Saving numdims
# #load("mf_F_wPid_varyb.Rdata")
# #load("mf_T_wPid_varyb.Rdata")
# numdims_wPid <- data.frame(numdims_mf_F_wPid = c(kbal_b.125x_f$numdims,
#                                                  kbal_b.25x_f$numdims,
#                                                  kbal_b.5x_f$numdims,
#                                                  kbal_b1x_f$numdims,
#                                                  kbal_b2x_f$numdims),
#                            numdims_mf_T_wPid = c(kbal_mf_b.125x_f$numdims,
#                                                  kbal_mf_b.25x_f$numdims,
#                                                  kbal_mf_b.5x_f$numdims,
#                                                  kbal_mf_b1x_f$numdims,
#                                                  kbal_mf_b2x_f$numdims),
#                            mfdims_mf_T_wPid = c(kbal_mf_b.125x_f$meanfirst.dims,
#                                                 kbal_mf_b.25x_f$meanfirst.dims,
#                                                 kbal_mf_b.5x_f$meanfirst.dims,
#                                                 kbal_mf_b1x_f$meanfirst.dims,
#                                                 kbal_mf_b2x_f$meanfirst.dims))
# rownames(numdims_wPid) <- c("b=.125x",
#                             "b=.25x",
#                             "b=.5x",
#                             "b=1x",
#                             "b=2x")
# save(numdims_wPid,
#      file = "numdims_wPid_full.Rdata")
# #load("mf_F_NoPid_varyb.Rdata")
# #load("mf_T_NoPid_varyb.Rdata")
# 
# numdims_NoPid <- data.frame(
#     numdims_mf_F_NoPid = c(#kbalout_nopid_b.0625x$numdims,
#         kbalout_nopid_b.125x_f$numdims,
#         kbalout_nopid_b.25x_f$numdims,
#         kbalout_nopid_b.5x_f$numdims,
#         kbalout_nopid_b1x_f$numdims,
#         kbalout_nopid_b2x_f$numdims),
#     numdims_mf_T_NoPid = c(#NA,
#         kbalout_mf_nopid_b.125x_f$numdims,
#         kbalout_mf_nopid_b.25x_f$numdims,
#         kbalout_mf_nopid_b.5x_f$numdims,
#         kbalout_mf_nopid_b1x_f$numdims,
#         kbalout_mf_nopid_b2x_f$numdims),
#     mfdims_mf_T_NoPid = c(#NA,
#         kbalout_mf_nopid_b.125x_f$meanfirst.dims,
#         kbalout_mf_nopid_b.25x_f$meanfirst.dims,
#         kbalout_mf_nopid_b.5x_f$meanfirst.dims,
#         kbalout_mf_nopid_b1x_f$meanfirst.dims,
#         kbalout_mf_nopid_b2x_f$meanfirst.dims)
# )
# 
# rownames(numdims_NoPid) <- c(#"b=.0625",
#     "b=.125x",
#     "b=.25x",
#     "b=.5x",
#     "b=1x",
#     "b=2x")
# save(numdims_NoPid,
#      file = "numdims_NoPid_full.Rdata")
# 
# 
# ############################# GET WEIGHTS #############################
# 
# 
# ##################### Survey Designs
# #####(1) Survey Designs from W PID runs
# ## (1.1) KPOP + MF =FALSE + W PID
# kbal_wt_b2x <- svydesign(~1, data = pew,
#                          weights = kbal_b2x_f$w[kbal_data_sampled ==1])
# kbal_wt_b1x <- svydesign(~1, data = pew,
#                          weights = kbal_b1x_f$w[kbal_data_sampled ==1])
# kbal_wt_b.5x <- svydesign(~1, data = pew,
#                           weights = kbal_b.5x_f$w[kbal_data_sampled ==1])
# kbal_wt_b.25x <- svydesign(~1, data = pew,
#                            weights = kbal_b.25x_f$w[kbal_data_sampled ==1])
# kbal_wt_b.125x <- svydesign(~1, data = pew,
#                             weights = kbal_b.125x_f$w[kbal_data_sampled ==1])
# 
# ## (1.2) KPOP + MF = TRUE + W PID
# kbal_mf_wt_b2x <- svydesign(~1, data = pew,
#                             weights = kbal_mf_b2x_f$w[kbal_data_sampled ==1])
# kbal_mf_wt_b1x <- svydesign(~1, data = pew,
#                             weights = kbal_mf_b1x_f$w[kbal_data_sampled ==1])
# kbal_mf_wt_b.5x <- svydesign(~1, data = pew,
#                              weights = kbal_mf_b.5x_f$w[kbal_data_sampled ==1])
# kbal_mf_wt_b.25x <- svydesign(~1, data = pew,
#                               weights = kbal_mf_b.25x_f$w[kbal_data_sampled ==1])
# kbal_mf_wt_b.125x <- svydesign(~1, data = pew,
#                                weights = kbal_mf_b.125x_f$w[kbal_data_sampled ==1])
# 
# save(kbal_wt_b2x, kbal_wt_b1x, kbal_wt_b.5x, kbal_wt_b.25x, kbal_wt_b.125x,
#      kbal_mf_wt_b2x, kbal_mf_wt_b1x, kbal_mf_wt_b.5x, kbal_mf_wt_b.25x,
#      kbal_mf_wt_b.125x,
#      file = "surveys_wPID_full.Rdata")
# 
# #### (2) Survey Desings from NO PID runs
# ## (2.1) KPOP + MF = FALSE + NO PID
# kbal_wt_b2x_nopid <- svydesign(~1, data = pew,
#                                weights = kbalout_nopid_b2x_f$w[kbal_data_sampled ==1])
# kbal_wt_b1x_nopid <- svydesign(~1, data = pew,
#                                weights = kbalout_nopid_b1x_f$w[kbal_data_sampled ==1])
# kbal_wt_b.5x_nopid <- svydesign(~1, data = pew,
#                                 weights = kbalout_nopid_b.5x_f$w[kbal_data_sampled ==1])
# kbal_wt_b.25x_nopid <- svydesign(~1, data = pew,
#                                  weights = kbalout_nopid_b.25x_f$w[kbal_data_sampled ==1])
# kbal_wt_b.125x_nopid <- svydesign(~1, data = pew,
#                                   weights = kbalout_nopid_b.125x_f$w[kbal_data_sampled ==1])
# # kbal_wt_b.0625x_nopid <- svydesign(~1, data = pew,
# #                                 weights = kbalout_nopid_b.0625x$w[kbal_data_sampled ==1])
# # save(kbal_wt_b2x_nopid, kbal_wt_b1x_nopid, kbal_wt_b.5x_nopid, 
# #           kbal_wt_b.25x_nopid, kbal_wt_b.125x_nopid, 
# #           file = "surveys_NoPID_mf_F_full.Rdata")
# 
# # (2.2) KPOP + MF=TRUE + NO PID
# kbal_wt_mf_nopid_b2x <- svydesign(~1, data = pew,
#                                   weights = kbalout_mf_nopid_b2x_f$w[kbal_data_sampled ==1])
# kbal_wt_mf_nopid_b1x <- svydesign(~1, data = pew,
#                                   weights = kbalout_mf_nopid_b1x_f$w[kbal_data_sampled ==1])
# kbal_wt_mf_nopid_b.5x <- svydesign(~1, data = pew,
#                                    weights = kbalout_mf_nopid_b.5x_f$w[kbal_data_sampled ==1])
# kbal_wt_mf_nopid_b.25x <- svydesign(~1, data = pew,
#                                     weights = kbalout_mf_nopid_b.25x_f$w[kbal_data_sampled ==1])
# kbal_wt_mf_nopid_b.125x <- svydesign(~1, data = pew,
#                                      weights = kbalout_mf_nopid_b.125x_f$w[kbal_data_sampled ==1])
# 
# 
# save(kbal_wt_b2x_nopid, kbal_wt_b1x_nopid, kbal_wt_b.5x_nopid,
#      kbal_wt_b.25x_nopid, kbal_wt_b.125x_nopid, #kbal_wt_b.0625x_nopid,
#      kbal_wt_mf_nopid_b2x, kbal_wt_mf_nopid_b1x, kbal_wt_mf_nopid_b.5x,
#      kbal_wt_mf_nopid_b.25x, kbal_wt_mf_nopid_b.125x,
#      file = "surveys_NoPID_full.Rdata")

###### MISC:
#(1) mf runs where did not require convergence
# kbal_wt_mf_nc_nopid_b.125 <- svydesign(~1, data = pew, 
#          weights = kbalout_mf_nebal_nopid_b.125x$w[kbal_data_sampled ==1])

#kbal_wt_mf_nc_nopid_b.25x <- svydesign(~1, data = pew, 
#                           weights = kbalout_mf_nebal_nopid$w[kbal_data_sampled ==1])

#(2) Weights from Frankenstein Runs
#from: U = opt svd(X)u (16) + svd(K with b=.25)(140) + Conv
# kbal_wt_frank <- svydesign(~1, data = pew, 
#                            weights = kbal_frank$w[kbal_data_sampled ==1])
# #from: U = opt svd(X)u (16) + svd(K with b=.25)(140) + NoConv
# kbal_wt_frank2 <- svydesign(~1, data = pew,
#                             weights = kbal_frank2$w[kbal_data_sampled ==1])


#Save Weights: scale so sum to n_sampled
# wts_wPid <- data.frame(
#     wt1 = weights(pew_lwt_1) / mean(weights(pew_lwt_1)),
#     wt1_pid = weights(pew_lwt_1_pid ) / mean(weights(pew_lwt_1_pid )) ,
#     wt2 = weights(pew_lwt_2) / mean(weights(pew_lwt_2)),
#     wt2_pid = weights(pew_lwt_2_pid ) / mean(weights(pew_lwt_2_pid )),
#     wt3 = weights(pew_ps_3) / mean(weights(pew_ps_3)),
#     #wt3_pid = weights(pew_ps_3_pid ) / mean(weights(pew_ps_3_pid )), #NAs cause issue
#     wt4 = weights(pew_lwt_4) / mean(weights(pew_lwt_4)),
#     wt4_pid = weights(pew_lwt_4_pid ) / mean(weights(pew_lwt_4_pid )),
#     
#     wtkbal_b2x = weights(kbal_wt_b2x) / mean(weights(kbal_wt_b2x)),
#     wtkbal_b1x = weights(kbal_wt_b1x) / mean(weights(kbal_wt_b1x)),
#     wtkbal_b.5x = weights(kbal_wt_b.5x) / mean(weights(kbal_wt_b.5x)),
#     wtkbal_b.25x = weights(kbal_wt_b.25x) / mean(weights(kbal_wt_b.25x)),
#     wtkbal_b.125x = weights(kbal_wt_b.125x) / mean(weights(kbal_wt_b.125x)),
#     
#     wtkbal_mf_b2x = weights(kbal_mf_wt_b2x) / mean(weights(kbal_mf_wt_b2x)),
#     wtkbal_mf_b1x = weights(kbal_mf_wt_b1x) / mean(weights(kbal_mf_wt_b1x)),
#     wtkbal_mf_b.5x = weights(kbal_mf_wt_b.5x) / mean(weights(kbal_mf_wt_b.5x)),
#     wtkbal_mf_b.25x = weights(kbal_mf_wt_b.25x) / mean(weights(kbal_mf_wt_b.25x)),
#     wtkbal_mf_b.125x = weights(kbal_mf_wt_b.125x) / mean(weights(kbal_mf_wt_b.125x)))
# 
# save(wts_wPid, file = "weights_wPid_full.Rdata")

#misc  
#wtkbal_frank_b.25x = weights(kbal_wt_frank) / mean(weights(kbal_wt_frank)),
#wtkbal_frank2_nc_b.25x = weights(kbal_wt_frank2) / mean(weights(kbal_wt_frank2)),


# wts_NoPid <- data.frame(
#     wt1 = weights(pew_lwt_1) / mean(weights(pew_lwt_1)),
#     wt1_pid = weights(pew_lwt_1_pid ) / mean(weights(pew_lwt_1_pid )) ,
#     wt2 = weights(pew_lwt_2) / mean(weights(pew_lwt_2)),
#     wt2_pid = weights(pew_lwt_2_pid ) / mean(weights(pew_lwt_2_pid )),
#     wt3 = weights(pew_ps_3) / mean(weights(pew_ps_3)),
#     #wt3_pid = weights(pew_ps_3_pid ) / mean(weights(pew_ps_3_pid )), #NAs cause issue
#     wt4 = weights(pew_lwt_4) / mean(weights(pew_lwt_4)),
#     wt4_pid = weights(pew_lwt_4_pid ) / mean(weights(pew_lwt_4_pid )),
#     
#     wtkbal_nopid_b2x = weights(kbal_wt_b2x_nopid) /
#         mean(weights(kbal_wt_b2x_nopid)),
#     wtkbal_nopid_b1x = weights(kbal_wt_b1x_nopid) /
#         mean(weights(kbal_wt_b1x_nopid)),
#     wtkbal_nopid_b.5x = weights(kbal_wt_b.5x_nopid) /
#         mean(weights(kbal_wt_b.5x_nopid)),
#     wtkbal_nopid_b.25x = weights(kbal_wt_b.25x_nopid)/
#         mean(weights(kbal_wt_b.25x_nopid)),
#     wtkbal_nopid_b.125x = weights(kbal_wt_b.125x_nopid)/
#         mean(weights(kbal_wt_b.125x_nopid)),
#     
#     wtkbal_nopid_mf_b2x = weights(kbal_wt_mf_nopid_b2x) /
#         mean(weights(kbal_wt_mf_nopid_b2x)),
#     wtkbal_nopid_mf_b1x = weights(kbal_wt_mf_nopid_b1x)/
#         mean(weights(kbal_wt_mf_nopid_b1x)),
#     wtkbal_nopid_mf_b.5x = weights(kbal_wt_mf_nopid_b.5x) /
#         mean(weights(kbal_wt_mf_nopid_b.5x)),
#     wtkbal_nopid_mf_b.25x = weights(kbal_wt_mf_nopid_b.25x) /
#         mean(weights(kbal_wt_mf_nopid_b.25x)),
#     wtkbal_nopid_mf_b.125x = weights(kbal_wt_mf_nopid_b.125x) /
#         mean(weights(kbal_wt_mf_nopid_b.125x)) )
# 
# save(wts_NoPid, file = "weights_NoPid_full.Rdata")
#sapply(wts, summary)
#sapply(wts, sd)


########################### Load Runs for Plots ############################

#NOTE: can load the weights and the survey designs rather than the huge kabl objects:
#then do not need to run code in get weights section


#for 500 max dims:

# load("cleaned data/500 Max SVD/surveys_NoPID.Rdata")
# load("cleaned data/500 Max SVD/surveys_wPID.Rdata")
# 
# load("cleaned data/500 Max SVD/weights_NoPid.Rdata")
# load("cleaned data/500 Max SVD/weights_wPid.Rdata")
# 
# load("cleaned data/500 Max SVD/numdims_wPid_varyb.Rdata")
# load("cleaned data/500 Max SVD/numdims_NoPid_varyb.Rdata")

#for Full SVD max dims:

#for full plot all you need to load is
load("cleaned data/Full SVD/comp_df_table_full.Rdata")

#full svd runs with out K and ksvd are here:
load("cleaned data/Full SVD/fullruns_noKsvd.Rdata")

#but if you want to build it internally use:
load("cleaned data/Full SVD/surveys_NoPID_full.Rdata")
load("cleaned data/Full SVD/surveys_wPID_full.Rdata")

load("cleaned data/Full SVD/weights_NoPid_full.Rdata")
load("cleaned data/Full SVD/weights_wPid_full.Rdata")

load("cleaned data/Full Max SVD/numdims_wPid_full.Rdata")
load("cleaned data/Full Max SVD/numdims_NoPid_full.Rdata")



################################ PLOTS ##########################
## Presidential vote estimates


####full SVD


### Actual results
#this needs to be loaded from local/dropbox 
pres <- readRDS("../data/election.rds")
natl_margin <- pres %>%
    summarise(margin = (sum(demtotal) - sum(reptotal)) /
                  (sum(demtotal) + sum(reptotal))) %>%
    as.numeric()
#natl_dvote <- sum(pres$demtotal)/ sum(pres$totalvotes)
#natl_margin

### Compare estimates
#run this whole thing to build the dataframe from the survey designs saved in 
#survey_xxxx files 
#OR just load the comp_df dataframe from the full svd runs saved in 
load("cleaned data/Full SVD/comp_df_table_full.Rdata")

# comp_df <- data.frame(
#     CCES = svycontrast(svymean(~recode_vote_2016, cces_awt, na.rm = TRUE),
#                        vote_contrast),
#     Pew_0 = svycontrast(svymean(~recode_vote_2016, pew_srs, na.rm = TRUE),
#                         vote_contrast),
#     Pew_1 = svycontrast(svymean(~recode_vote_2016, pew_lwt_1, na.rm = TRUE),
#                         vote_contrast),
#     Pew_1_pid = svycontrast(svymean(~recode_vote_2016, pew_lwt_1_pid, na.rm = TRUE),
#                         vote_contrast),
#     Pew_2 = svycontrast(svymean(~recode_vote_2016, pew_lwt_2, na.rm = TRUE),
#                         vote_contrast),
#     Pew_2_pid = svycontrast(svymean(~recode_vote_2016, pew_lwt_2_pid, na.rm = TRUE),
#                         vote_contrast),
#     Pew_3 = svycontrast(svymean(~recode_vote_2016, pew_ps_3, na.rm = TRUE),
#                        vote_contrast),
#     Pew_4 = svycontrast(svymean(~recode_vote_2016, pew_lwt_4, na.rm = TRUE),
#                         vote_contrast),
#     Pew_4_pid = svycontrast(svymean(~recode_vote_2016, pew_lwt_4_pid, na.rm = TRUE),
#                         vote_contrast),
#     
#     #b = 2
#     ## no meanfirst
#     Pew_kbal_b2x_nopid = svycontrast(svymean(~recode_vote_2016,
#                                               kbal_wt_b2x_nopid, na.rm = TRUE),
#                                       vote_contrast),
#     Pew_kbal_b2x = svycontrast(svymean(~recode_vote_2016, kbal_wt_b2x,
#                                         na.rm = TRUE), vote_contrast),
#     ## meanfirst
#     Pew_kbal_mf_b2x_nopid = svycontrast(svymean(~recode_vote_2016,
#                                                 kbal_wt_mf_nopid_b2x, na.rm = TRUE),
#                                      vote_contrast),
#     Pew_kbal_mf_b2x = svycontrast(svymean(~recode_vote_2016, kbal_mf_wt_b2x,
#                                        na.rm = TRUE), vote_contrast),
#     #b = 1
#     ## no meanfirst
#     Pew_kbal_b1x_nopid = svycontrast(svymean(~recode_vote_2016,
#                                              kbal_wt_b1x_nopid, na.rm = TRUE),
#                                      vote_contrast),
#     Pew_kbal_b1x = svycontrast(svymean(~recode_vote_2016, kbal_wt_b1x,
#                                        na.rm = TRUE), vote_contrast),
#     ## meanfirst 
#     Pew_kbal_mf_b1x_nopid = svycontrast(svymean(~recode_vote_2016,
#                                                 kbal_wt_mf_nopid_b1x, na.rm = TRUE),
#                                         vote_contrast),
#     Pew_kbal_mf_b1x = svycontrast(svymean(~recode_vote_2016, kbal_mf_wt_b1x,
#                                           na.rm = TRUE), vote_contrast),
#     
#     #b = 0.5
#     ##no meanfirst
#     Pew_kbal_b.5x_nopid = svycontrast(svymean(~recode_vote_2016,
#                                               kbal_wt_b.5x_nopid, na.rm = TRUE),
#                                       vote_contrast),
#     Pew_kbal_b.5x = svycontrast(svymean(~recode_vote_2016, kbal_wt_b.5x,
#                                         na.rm = TRUE), vote_contrast),
#     ## meanfirst
#     Pew_kbal_mf_b.5x_nopid = svycontrast(svymean(~recode_vote_2016,
#                                                 kbal_wt_mf_nopid_b.5x, na.rm = TRUE),
#                                         vote_contrast),
#     Pew_kbal_mf_b.5x = svycontrast(svymean(~recode_vote_2016, kbal_mf_wt_b.5x,
#                                           na.rm = TRUE), vote_contrast),
#     
#     #b = 0.25
#     ## no meanfirst
#     Pew_kbal_b.25x_nopid = svycontrast(svymean(~recode_vote_2016,
#                                                kbal_wt_b.25x_nopid, na.rm = TRUE),
#                                        vote_contrast),
#     Pew_kbal_b.25x = svycontrast(svymean(~recode_vote_2016, kbal_wt_b.25x, 
#                                          na.rm = TRUE), vote_contrast),
#     ## meanfirst 
#     Pew_kbal_mf_b.25x_nopid = svycontrast(svymean(~recode_vote_2016,
#                                                 kbal_wt_mf_nopid_b.25x, na.rm = TRUE),
#                                         vote_contrast),
#     Pew_kbal_mf_b.25x = svycontrast(svymean(~recode_vote_2016, kbal_mf_wt_b.25x,
#                                           na.rm = TRUE), vote_contrast),
#     
#     #b = 0.125
#     ## nomeanfirst
#     Pew_kbal_b.125x_nopid = svycontrast(svymean(~recode_vote_2016,
#                                                 kbal_wt_b.125x_nopid, na.rm = TRUE),
#                                         vote_contrast),
#     Pew_kbal_b.125x = svycontrast(svymean(~recode_vote_2016, kbal_wt_b.125x, 
#                                           na.rm = TRUE), vote_contrast),
#     ## meanfirst
#     Pew_kbal_mf_b.125x_nopid = svycontrast(svymean(~recode_vote_2016,
#                                                 kbal_wt_mf_nopid_b.125x, na.rm = TRUE),
#                                         vote_contrast),
#     Pew_kbal_mf_b.125x = svycontrast(svymean(~recode_vote_2016, kbal_mf_wt_b.125x,
#                                           na.rm = TRUE), vote_contrast),
# 
#     #b =0.0625 (only did for nopid + mf = F)
#     Pew_kbal_b.0625x_nopid = svycontrast(svymean(~recode_vote_2016,
#                                                 kbal_wt_b.0625x_nopid, na.rm = TRUE),
#                                         vote_contrast)
#    
# ) %>%
#     pivot_longer(cols = everything(),
#                  names_to = c("source", ".value"), 
#                  names_pattern = "(.*)\\.(.*)") %>%
#     rename(est = nlcon) %>%
#     mutate(err = est - natl_margin, 
#            source = str_replace(source, "_", " ")) %>%
#     mutate(err_target = est - est[source == "CCES"]) %>%
#     mutate(source = factor(source,  
#                          levels = c("CCES", 
#                                     "Pew 0", 
#                                     "Pew 1",
#                                     "Pew 1_pid",
#                                     "Pew 2",
#                                     "Pew 2_pid",
#                                     "Pew 3",
#                                     "Pew 4",
#                                     "Pew 4_pid" ,
#                                     
#                                     "Pew kbal_b2x_nopid", 
#                                     "Pew kbal_b2x",
#                                     "Pew kbal_mf_b2x_nopid", 
#                                     "Pew kbal_mf_b2x",
#                                     
#                                     "Pew kbal_b1x_nopid", 
#                                     "Pew kbal_b1x",
#                                     "Pew kbal_mf_b1x_nopid", 
#                                     "Pew kbal_mf_b1x",
#                                     
#                                     "Pew kbal_b.5x_nopid", 
#                                     "Pew kbal_b.5x",
#                                     "Pew kbal_mf_b.5x_nopid", 
#                                     "Pew kbal_mf_b.5x",
#                                     
#                                     "Pew kbal_b.25x_nopid", 
#                                     "Pew kbal_b.25x",
#                                     "Pew kbal_mf_b.25x_nopid", 
#                                     "Pew kbal_mf_b.25x",
#                                     
#                                     "Pew kbal_b.125x_nopid", 
#                                     "Pew kbal_b.125x",
#                                     "Pew kbal_mf_b.125x_nopid", 
#                                     "Pew kbal_mf_b.125x",
#                                     
#                                     "Pew kbal_b.0625x_nopid"
#                          ),
#                          labels = c("CCES\n(Target)", 
#                                     "Pew\nUnweighted", 
#                                     "Pew\nRaking (demographics)", 
#                                     "Pew + PID \nRaking (demographics)", 
#                                     "Pew\nRaking (demographics +education)",
#                                     "Pew + PID \nRaking (demographics +education)",
#                                     "Pew\nPost-Stratification",
#                                     "Pew\nRaking (Post-Hoc)", 
#                                     "Pew + PID \nRaking (Post-Hoc)" , 
#                                     #b=2
#                                     paste0("Pew KPOP: b=2x
#                                            (dimKu = ",
#                                            numdims_NoPid["b=2x","numdims_mf_F_NoPid"]
#                                            ,")"),
#                                     paste0("Pew KPOP + PID: b=2x
#                                            (dimKu = ",
#                                            numdims_wPid["b=2x", "numdims_mf_F_wPid"]
#                                            ,")"),
#                                     paste0("Pew KPOP + MF: b=2x
#                                            (dimXu = ",
#                                            numdims_NoPid["b=2x","mfdims_mf_T_NoPid"],
#                                            ", dimKu = ",
#                                            numdims_NoPid["b=2x","numdims_mf_T_NoPid"],
#                                            ")"),
#                                     paste0("Pew KPOP + MF + PID: b=2x
#                                            (dimXu = ",
#                                            numdims_wPid["b=2x", "mfdims_mf_T_wPid"],
#                                            ", dimKu = ",
#                                            numdims_wPid["b=2x", "numdims_mf_T_wPid"]
#                                            ,")"),
#                                     #b=1
#                                     paste0("Pew KPOP: b=1x
#                                            (dimKu = ",
#                                            numdims_NoPid["b=1x","numdims_mf_F_NoPid"]
#                                            ,")"),
#                                     paste0("Pew KPOP + PID: b=1x
#                                            (dimKu = ",
#                                            numdims_wPid["b=1x", "numdims_mf_F_wPid"]
#                                            ,")"),
#                                     paste0("Pew KPOP + MF: b=1x
#                                            (dimXu = ",
#                                            numdims_NoPid["b=1x","mfdims_mf_T_NoPid"],
#                                            ", dimKu = ",
#                                            numdims_NoPid["b=1x","numdims_mf_T_NoPid"],
#                                            ")"),
#                                     paste0("Pew KPOP + MF + PID: b=1x
#                                            (dimXu = ",
#                                            numdims_wPid["b=1x", "mfdims_mf_T_wPid"],
#                                            ", dimKu = ",
#                                            numdims_wPid["b=1x", "numdims_mf_T_wPid"]
#                                            ,")"),
#                                     #b = 0.5
#                                     paste0("Pew KPOP: b=0.5x
#                                            (dimKu = ",
#                                            numdims_NoPid["b=.5x","numdims_mf_F_NoPid"]
#                                            ,")"),
#                                     paste0("Pew KPOP + PID: b=0.5x
#                                            (dimKu = ",
#                                            numdims_wPid["b=.5x", "numdims_mf_F_wPid"]
#                                            ,")"),
#                                     paste0("Pew KPOP + MF: b=0.5x
#                                            (dimXu = ",
#                                            numdims_NoPid["b=.5x","mfdims_mf_T_NoPid"],
#                                            ", dimKu = ",
#                                            numdims_NoPid["b=.5x","numdims_mf_T_NoPid"],
#                                            ")"),
#                                     paste0("Pew KPOP + MF + PID: b=0.5x
#                                            (dimXu = ",
#                                            numdims_wPid["b=.5x", "mfdims_mf_T_wPid"],
#                                            ", dimKu = ",
#                                            numdims_wPid["b=.5x", "numdims_mf_T_wPid"]
#                                            ,")"),
#                                     #b = 0.25
#                                     paste0("Pew KPOP: b=0.25x
#                                            (dimKu = ",
#                                            numdims_NoPid["b=.25x","numdims_mf_F_NoPid"]
#                                            ,")"),
#                                     paste0("Pew KPOP + PID: b=0.25x
#                                            (dimKu = ",
#                                            numdims_wPid["b=.25x", "numdims_mf_F_wPid"]
#                                            ,")"),
#                                     paste0("Pew KPOP + MF: b=0.25x
#                                            (dimXu = ",
#                                            numdims_NoPid["b=.25x","mfdims_mf_T_NoPid"],
#                                            ", dimKu = ",
#                                            numdims_NoPid["b=.25x","numdims_mf_T_NoPid"],
#                                            ")"),
#                                     paste0("Pew KPOP + MF + PID: b=0.25x
#                                            (dimXu = ",
#                                            numdims_wPid["b=.25x", "mfdims_mf_T_wPid"],
#                                            ", dimKu = ",
#                                            numdims_wPid["b=.25x", "numdims_mf_T_wPid"]
#                                            ,")"),
#                                     #b = .125
#                                     paste0("Pew KPOP: b=0.125x
#                                            (dimKu = ",
#                                            numdims_NoPid["b=.125x","numdims_mf_F_NoPid"]
#                                            ,")"),
#                                     paste0("Pew KPOP + PID: b=0.125x
#                                            (dimKu = ",
#                                            numdims_wPid["b=.125x", "numdims_mf_F_wPid"]
#                                            ,")"),
#                                     paste0("Pew KPOP + MF: b=0.125x
#                                            (dimXu = ",
#                                            numdims_NoPid["b=.125x","mfdims_mf_T_NoPid"],
#                                            ", dimKu = ",
#                                            numdims_NoPid["b=.125x","numdims_mf_T_NoPid"],
#                                            ")"),
#                                     paste0("Pew KPOP + MF + PID: b=0.125x
#                                            (dimXu = ",
#                                            numdims_wPid["b=.125x", "mfdims_mf_T_wPid"],
#                                            ", dimKu = ",
#                                            numdims_wPid["b=.125x", "numdims_mf_T_wPid"]
#                                            ,")"),
#                                     #b = 0.0625 for nopid mf=F only
#                                     paste0("Pew KPOP: b=0.0625x
#                                            (dimKu = ",
#                                            numdims_NoPid["b=.0625","numdims_mf_F_NoPid"]
#                                            ,")")
#                                     )
#                          ))
# comp_df$PID <- NA
# comp_df$PID[grep("PID", levels(comp_df$source))] <- "PID"
# comp_df$PID[-grep("PID", levels(comp_df$source)) ]<- "No PID"
# comp_df$PID[1:2]<- "Not Weighted" #making CCES target different color
# 
# comp_df$MF <- NA
# comp_df$MF[grep("MF", levels(comp_df$source))] <- "MF"
# comp_df$MF[-grep("MF", levels(comp_df$source)) ]<- "No MF"
# comp_df$MF[1:9]<- "NA" #making CCES target different color


# comp_df %>%  #filter(PID != "No PID") %>%
#     ggplot() +
#     aes(x = source, y = est, ymin = est - 1.96*SE, ymax = est + 1.96*SE, color = as.factor(PID),shape=as.factor(MF)) +
#     
#     geom_hline(yintercept = c(0, natl_margin, comp_df$est[comp_df$source == "CCES\n(Target)"]),
#                linetype = c("solid", "dashed", "longdash"),
#                color = c("black", "gray", "black")) +
#     geom_pointrange() +
#     scale_y_continuous(breaks = c(natl_margin, seq(-.05, .1, .05)),
#                        minor_breaks = NULL,
#                        labels = scales::percent_format()) +
#     theme_bw() +
#     theme(plot.title = element_text(hjust = 0.5)) +
#     
#     labs(x = NULL, y = "Estimated Margin (95% CI)") +
#     ggtitle("Estimates of Clinton National Popular Vote Margin") +
#     theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
#     annotate(geom = "text", x = 31.25, y = natl_margin, label = "True\nNational\nMargin", hjust = -0.1, angle = -90, color = "gray") +
#     annotate(geom = "text", x = 31.35, y = comp_df$est[comp_df$source == "CCES\n(Target)"], label = "CCES Estimated\n  Margin", hjust = -0.1, angle = 90) +
#     theme(legend.title = element_blank())

#ggsave("plots/fullSVD.pdf", width = 12, height = 6)


######################## Cleaned Plots #######################
#################################################
# Clean Plot:
##### (1) w PID: MF and no MF at optimal B + all raking/stratifying

#get optimal b's:
load("Biasbound Tables/bb_dat_wEst_full.Rdata")

### MF = F + Pid: 0.5 is best
#best in terms of bias bound ratio
rownames(bb_comp)[which(abs(bb_comp$bb_ratio) == max(abs(bb_comp$bb_ratio)))]

### MF = T + Pid: 2 is best
#best in terms of bias bound ratio
rownames(bb_comp_mf)[which(abs(bb_comp_mf$bb_ratio) == max(abs(bb_comp_mf$bb_ratio)))]

#load big df: 
load("cleaned data/Full SVD/comp_df_table_full.Rdata")

comp_df$b <- NA
comp_df[grep("b", levels(comp_df$source)), "b"] <- c(rep(2,4), 
                                                rep(1,4),
                                                rep(.5, 4),
                                                rep(.25,4),
                                                rep(.125,4))

# a verrrrrrryyy stupid way to do this
comp_clean <- comp_df %>% filter(PID != "No PID" & b %in% c(NA, 2,0.5)) %>% 
    filter(!(source %in% c("Pew KPOP + PID: b=2x\n                                           (dimKu = 13)","Pew KPOP + MF + PID: b=0.5x\n                                           (dimXu = 18, dimKu = 17)"))) %>% 
    mutate(source_clean = factor(source, 
                                 levels = c(as.character(comp_clean$source)),
                                 labels = c("CCES\n(Target)", 
                            "Pew\nUnweighted", 
                            "Pew Raking\n(demographics)",
                            "Pew Raking\n(demographics + education)",
                            "Pew Raking\n(Post-Hoc)", 
                            "Pew KPOP + MF",
                            "Pew KPOP\n")))


#### !!! issue bc still do not have working post-stratificaiton with pid

comp_clean %>%
    ggplot() +
    aes(x = source_clean, y = est, ymin = est - 1.96*SE, ymax = est + 1.96*SE) +

    geom_hline(yintercept = c(0, natl_margin, comp_clean$est[comp_clean$source == "CCES\n(Target)"]),
               linetype = c("solid", "dashed", "longdash"),
               color = c("black", "gray", "black")) +
    geom_pointrange() +
    scale_y_continuous(breaks = c(natl_margin, seq(-.05, .1, .05)),
                       minor_breaks = NULL,
                       labels = scales::percent_format()) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +

    labs(x = NULL, y = "Estimated Margin (95% CI)") +
    ggtitle("Estimates of Clinton National Popular Vote Margin") +
    theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
    annotate(geom = "text", x = 8.25, y = natl_margin, label = "True\nNational\nMargin", hjust = -0.1, angle = -90, color = "gray") +
    annotate(geom = "text", x = 8.35, y = comp_clean$est[comp_clean$source == "CCES\n(Target)"], label = "CCES Estimated\n  Margin", hjust = -0.1, angle = 90) +
    theme(legend.title = element_blank())

ggsave("plots/cleaned_v1.pdf", width = 12, height = 6)

