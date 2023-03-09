## KPOP RUNS FOR 2016 ELECTION

## SETUP

### Packages
library(tidyverse)
library(survey)
#devtools::install_github("chadhazlett/KBAL", ref = "master") 
library(kbal)

path_save = "/pewavg_nobucketfix/"
path_data= "./"
############################## Load Data ##########################
## SURVEY DATA (PEW)
### Load
pew <- readRDS(paste0(path_data, "pew.rds"))
### Make survey design 
pew_srs <- svydesign(ids = ~1, data = pew)


### Unweighted survey estimates of presidential vote
#function to get vote margin
vote_contrast <- quote((recode_vote_2016Democrat - recode_vote_2016Republican) /
                           (recode_vote_2016Democrat + recode_vote_2016Republican))

## AUXILIARY INFORMATION (CCES)
### Load
cces <- readRDS(paste0(path_data, "cces.rds"))
### Drop invalid cases
cces <- cces %>%
    filter((CC16_401 == "I definitely voted in the General Election.") &
               !is.na(commonweight_vv_post))

### Make survey design
cces_awt <- svydesign(ids = ~1, weights = ~commonweight_vv_post, data = cces)


############################## Run Kbal ##########################


##### Kbal Data Prep:
#filling in NAs in age with mean est. in cces unweighted
sum(is.na(pew$recode_age))
pew$recode_age[is.na(pew$recode_age)] <-  mean(pew$recode_age, na.rm = TRUE)
#adjust age bucket assignments:
#pew$recode_age_bucket[is.na(pew$recode_age)] <- "18 to 35"


# W PID included
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

# NO PID: 
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



########################## (1) KPOP + W PID #########################

##################### (1.1) KPOP + MF = F + PID 

#b=2
kbal_b2x_f <- kbal(allx=kbal_data,
                   sampled = kbal_data_sampled,
                   b = 2 * ncol(kbal_data),
                   incrementby = 1,
                   fullSVD = TRUE,
                   meanfirst = FALSE,
        population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
                   sampledinpop = FALSE)
#b = 1
kbal_b1x_f <- kbal(allx=kbal_data,
                   sampled = kbal_data_sampled,
                   b = 1 * ncol(kbal_data),
                   incrementby = 1,
                   fullSVD = TRUE,
                   meanfirst = FALSE,
        population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
                   sampledinpop = FALSE)
#b = .5
kbal_b.5x_f <- kbal(allx=kbal_data,
                    sampled = kbal_data_sampled,
                    b = 0.5 * ncol(kbal_data),
                    fullSVD = TRUE,
                    meanfirst = FALSE,
                    incrementby = 1,
         population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
                    sampledinpop = FALSE)
#b = .25
kbal_b.25x_f <- kbal(allx=kbal_data,
                     sampled = kbal_data_sampled,
                     b = 0.25 * ncol(kbal_data),
                     incrementby = 1,
                     fullSVD = TRUE,
                     meanfirst = FALSE,
    population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
                     sampledinpop = FALSE)
#b = .125
kbal_b.125x_f <- kbal(allx=kbal_data,
                      sampled = kbal_data_sampled,
                      b = 0.125 * ncol(kbal_data),
                      incrementby = 1,
                      fullSVD = TRUE,
                      meanfirst = FALSE,
        population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
                      sampledinpop = FALSE)

save(kbal_b.125x_f,kbal_b.25x_f,kbal_b.5x_f,kbal_b1x_f,kbal_b2x_f,
     file = paste0(path_save,"mf_F_wPid.Rdata"))

######################## (1.2) KPOP + MF = T + W PID 

#b=2
kbal_mf_b2x_f <- kbal(allx=kbal_data,
                      sampled = kbal_data_sampled,
                      b = 2 * ncol(kbal_data),
                      incrementby = 1,
                      fullSVD = TRUE,
                      meanfirst = TRUE,
                      ebal.convergence = TRUE,
        population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
                      sampledinpop = FALSE)
# b = 1
kbal_mf_b1x_f <- kbal(allx=kbal_data,
                      sampled = kbal_data_sampled,
                      b = 1 * ncol(kbal_data),
                      incrementby = 1,
                      fullSVD = TRUE,
                      meanfirst = TRUE,
                      ebal.convergence = TRUE,
        population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
                      sampledinpop = FALSE)
#b = .5
kbal_mf_b.5x_f <- kbal(allx=kbal_data,
                       sampled = kbal_data_sampled,
                       b = 0.5 * ncol(kbal_data),
                       fullSVD = TRUE,
                       meanfirst = TRUE,
                       ebal.convergence = TRUE,
                       incrementby = 1,
        population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
                       sampledinpop = FALSE)
#b = .25
kbal_mf_b.25x_f <- kbal(allx=kbal_data,
                        sampled = kbal_data_sampled,
                        b = 0.25 * ncol(kbal_data),
                        incrementby = 1,
                        fullSVD = TRUE,
                        meanfirst = TRUE,
                        ebal.convergence = TRUE,
    population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
                        sampledinpop = FALSE)
#b = .125
kbal_mf_b.125x_f <- kbal(allx=kbal_data,
                         sampled = kbal_data_sampled,
                         b = 0.125 * ncol(kbal_data),
                         incrementby = 1,
                         fullSVD = TRUE,
                         meanfirst = TRUE,
                         ebal.convergence = TRUE,
        population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
                         sampledinpop = FALSE)

save(kbal_mf_b.125x_f, kbal_mf_b.25x_f, kbal_mf_b.5x_f,kbal_mf_b1x_f,kbal_mf_b2x_f,
     file = paste0(path_save, "mf_T_wPID.Rdata"))


############################## (2) KPOP + NO PID ##################################

########################### (2.1) KPOP + MF=F + NO PID 

#b = 2
kbal_nopid_b2x_f <- kbal(allx=kbal_data_nopid,
                            sampled = kbal_data_sampled,
                            b = 2 * ncol(kbal_data_nopid),
                            incrementby = 1,
                            fullSVD = TRUE,
                            meanfirst = FALSE,
            population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
                            sampledinpop = FALSE)
#b = 1
kbal_nopid_b1x_f <- kbal(allx=kbal_data_nopid,
                            sampled = kbal_data_sampled,
                            b = 1 * ncol(kbal_data_nopid),
                            incrementby = 1,
                            fullSVD = TRUE,
                            meanfirst = FALSE,
            population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
                            sampledinpop = FALSE)
#b=0.5
kbal_nopid_b.5x_f <- kbal(allx=kbal_data_nopid,
                             sampled = kbal_data_sampled,
                             b = 0.5 * ncol(kbal_data_nopid),
                             incrementby = 1,
                             fullSVD = TRUE,
                             meanfirst = FALSE,
            population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
                             sampledinpop = FALSE)
#b=0.25
kbal_nopid_b.25x_f <- kbal(allx=kbal_data_nopid,
                              sampled = kbal_data_sampled,
                              b = 0.25 * ncol(kbal_data_nopid),
                              incrementby = 1,
                              fullSVD = TRUE,
                              meanfirst = FALSE,
            population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
                              sampledinpop = FALSE)

# b=0.125
kbal_nopid_b.125x_f <- kbal(allx=kbal_data_nopid,
                               sampled = kbal_data_sampled,
                               b = 0.125 * ncol(kbal_data_nopid),
                               incrementby = 1,
                               fullSVD = TRUE,
                               meanfirst = FALSE,
        population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
                               sampledinpop = FALSE)


save( kbal_nopid_b.125x_f,kbal_nopid_b.25x_f,
      kbal_nopid_b.5x_f, kbal_nopid_b1x_f, kbal_nopid_b2x_f,
      file = paste0(path_save, "mf_F_NoPid.Rdata"))

########################## (2.2) KPOP + MF=TRUE + NO PID 

#b=2x
kbal_mf_nopid_b2x_f <- kbal(allx=kbal_data_nopid,
                            sampled = kbal_data_sampled,
                            b = 2 * ncol(kbal_data_nopid),
                            meanfirst = TRUE,
                            ebal.convergence = TRUE,
                            fullSVD = TRUE,
            population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
                            sampledinpop = FALSE)
#b=1x 
kbal_mf_nopid_b1x_f <- kbal(allx=kbal_data_nopid,
                               sampled = kbal_data_sampled,
                               b = 1 * ncol(kbal_data_nopid),
                               meanfirst = TRUE, 
                               ebal.convergence = TRUE,
                               fullSVD = TRUE,
            population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
                               sampledinpop = FALSE)
#b=0.5x
kbal_mf_nopid_b.5x_f <- kbal(allx=kbal_data_nopid,
                                sampled = kbal_data_sampled,
                                b = 0.5 * ncol(kbal_data_nopid),
                                meanfirst = TRUE,
                                ebal.convergence = TRUE,
                                fullSVD = TRUE,
            population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
                                sampledinpop = FALSE)
#b=0.25
kbal_mf_nopid_b.25x_f <- kbal(allx=kbal_data_nopid,
                              sampled = kbal_data_sampled,
                              b = 0.25 * ncol(kbal_data_nopid),
                              meanfirst = TRUE,
                              ebal.convergence = TRUE,
                              fullSVD = TRUE,
            population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
                              sampledinpop = FALSE)
#b=0.125x
kbal_mf_nopid_b.125x_f <- kbal(allx=kbal_data_nopid,
                               sampled = kbal_data_sampled,
                               b = 0.125 * ncol(kbal_data_nopid),
                               meanfirst = TRUE,
                               ebal.convergence = TRUE,
                               fullSVD = TRUE,
            population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
                               sampledinpop = FALSE)

save(kbal_mf_nopid_b2x_f, kbal_mf_nopid_b1x_f, kbal_mf_nopid_b.5x_f,
     kbal_mf_nopid_b.25x_f, kbal_mf_nopid_b.125x_f,
     file = paste0(path_save,"mf_T_NoPid.Rdata"))





############################## Choice of B ############################

####################### (1) KPOP + PID #######################

##### (1.1) Kbal + MF = FALSE + PID
bb_comp <- data.frame(biasbound_orig = c(kbal_b2x_f$biasbound.orig,
                                         kbal_b1x_f$biasbound.orig,
                                         kbal_b.5x_f$biasbound.orig,
                                         kbal_b.25x_f$biasbound.orig,
                                         kbal_b.125x_f$biasbound.orig),
                      biasbound_opt = c(kbal_b2x_f$biasbound.opt,
                                        kbal_b1x_f$biasbound.opt,
                                        kbal_b.5x_f$biasbound.opt,
                                        kbal_b.25x_f$biasbound.opt,
                                        kbal_b.125x_f$biasbound.opt),
                      L1_orig = c(kbal_b2x_f$L1.orig,
                                  kbal_b1x_f$L1.orig,
                                  kbal_b.5x_f$L1.orig,
                                  kbal_b.25x_f$L1.orig,
                                  kbal_b.125x_f$L1.orig),
                      L1_opt = c(kbal_b2x_f$L1.opt,
                                 kbal_b1x_f$L1.opt,
                                 kbal_b.5x_f$L1.opt,
                                 kbal_b.25x_f$L1.opt,
                                 kbal_b.125x_f$L1.opt)
)

bb_comp <- bb_comp %>% mutate(bb_ratio = biasbound_orig/biasbound_opt,
                              L1_ratio = L1_orig/L1_opt) %>%
    select(bb_ratio,L1_ratio, biasbound_orig, biasbound_opt,  L1_orig, L1_opt)

rownames(bb_comp) <- c("b=2x",
                       "b=1x",
                       "b=.5x",
                       "b=.25x",
                       "b=.125x")

##### (1.2) Kbal + MF = TRUE + PID
bb_comp_mf <- data.frame(biasbound_orig = c(kbal_mf_b2x_f$biasbound.orig,
                                            kbal_mf_b1x_f$biasbound.orig,
                                            kbal_mf_b.5x_f$biasbound.orig,
                                            kbal_mf_b.25x_f$biasbound.orig,
                                            kbal_mf_b.125x_f$biasbound.orig),
                         biasbound_opt = c(kbal_mf_b2x_f$biasbound.opt,
                                           kbal_mf_b1x_f$biasbound.opt,
                                           kbal_mf_b.5x_f$biasbound.opt,
                                           kbal_mf_b.25x_f$biasbound.opt,
                                           kbal_mf_b.125x_f$biasbound.opt),
                         L1_orig = c(kbal_mf_b2x_f$L1.orig,
                                     kbal_mf_b1x_f$L1.orig,
                                     kbal_mf_b.5x_f$L1.orig,
                                     kbal_mf_b.25x_f$L1.orig,
                                     kbal_mf_b.125x_f$L1.orig),
                         L1_opt = c(kbal_mf_b2x_f$L1.opt,
                                    kbal_mf_b1x_f$L1.opt,
                                    kbal_mf_b.5x_f$L1.opt,
                                    kbal_mf_b.25x_f$L1.opt,
                                    kbal_mf_b.125x_f$L1.opt)
)

bb_comp_mf <- bb_comp_mf %>% mutate(bb_ratio = biasbound_orig/biasbound_opt,
                                    L1_ratio = L1_orig/L1_opt) %>%
    select(bb_ratio,L1_ratio, biasbound_orig, biasbound_opt,  L1_orig, L1_opt)

rownames(bb_comp_mf) <- c("b=2x",
                          "b=1x",
                          "b=.5x",
                          "b=.25x",
                          "b=.125x"
)


############################ (2) KPOP + NO PID ###################################

##### (2.1) Kbal + MF = FALSE + NO PID
bb_comp_nopid <- data.frame(biasbound_orig = c(kbal_nopid_b2x_f$biasbound.orig,
                                               kbal_nopid_b1x_f$biasbound.orig,
                                               kbal_nopid_b.5x_f$biasbound.orig,
                                               kbal_nopid_b.25x_f$biasbound.orig,
                                               kbal_nopid_b.125x_f$biasbound.orig),
                            biasbound_opt = c(kbal_nopid_b2x_f$biasbound.opt,
                                              kbal_nopid_b1x_f$biasbound.opt,
                                              kbal_nopid_b.5x_f$biasbound.opt,
                                              kbal_nopid_b.25x_f$biasbound.opt,
                                              kbal_nopid_b.125x_f$biasbound.opt),
                            L1_orig = c(kbal_nopid_b2x_f$L1.orig,
                                        kbal_nopid_b1x_f$L1.orig,
                                        kbal_nopid_b.5x_f$L1.orig,
                                        kbal_nopid_b.25x_f$L1.orig,
                                        kbal_nopid_b.125x_f$L1.orig),
                            L1_opt = c(kbal_nopid_b2x_f$L1.opt,
                                       kbal_nopid_b1x_f$L1.opt,
                                       kbal_nopid_b.5x_f$L1.opt,
                                       kbal_nopid_b.25x_f$L1.opt,
                                       kbal_nopid_b.125x_f$L1.opt)
)

bb_comp_nopid <- bb_comp_nopid %>% mutate(bb_ratio = biasbound_orig/biasbound_opt,
                                          L1_ratio = L1_orig/L1_opt) %>%
    select(bb_ratio,L1_ratio, biasbound_orig, biasbound_opt,  L1_orig, L1_opt)

rownames(bb_comp_nopid) <- c("b=2x",
                             "b=1x",
                             "b=.5x",
                             "b=.25x",
                             "b=.125x")

############ (2.2) Kbal + MF = TRUE + NO PID

bb_comp_mf_nopid <- data.frame(biasbound_orig = c(kbal_mf_nopid_b2x_f$biasbound.orig,
                                                  kbal_mf_nopid_b1x_f$biasbound.orig,
                                                  kbal_mf_nopid_b.5x_f$biasbound.orig,
                                                  kbal_mf_nopid_b.25x_f$biasbound.orig,
                                                  kbal_mf_nopid_b.125x_f$biasbound.orig),
                               biasbound_opt = c(kbal_mf_nopid_b2x_f$biasbound.opt,
                                                 kbal_mf_nopid_b1x_f$biasbound.opt,
                                                 kbal_mf_nopid_b.5x_f$biasbound.opt,
                                                 kbal_mf_nopid_b.25x_f$biasbound.opt,
                                                 kbal_mf_nopid_b.125x_f$biasbound.opt),
                               L1_orig = c(kbal_mf_nopid_b2x_f$L1.orig,
                                           kbal_mf_nopid_b1x_f$L1.orig,
                                           kbal_mf_nopid_b.5x_f$L1.orig,
                                           kbal_mf_nopid_b.25x_f$L1.orig,
                                           kbal_mf_nopid_b.125x_f$L1.orig),
                               L1_opt = c(kbal_mf_nopid_b2x_f$L1.opt,
                                          kbal_mf_nopid_b1x_f$L1.opt,
                                          kbal_mf_nopid_b.5x_f$L1.opt,
                                          kbal_mf_nopid_b.25x_f$L1.opt,
                                          kbal_mf_nopid_b.125x_f$L1.opt)
)
#
bb_comp_mf_nopid <- bb_comp_mf_nopid %>% mutate(bb_ratio = biasbound_orig/biasbound_opt,
                                                L1_ratio = L1_orig/L1_opt) %>%
    select(bb_ratio,L1_ratio, biasbound_orig, biasbound_opt,  L1_orig, L1_opt)
rownames(bb_comp_mf_nopid) <- c("b=2x",
                                "b=1x",
                                "b=.5x",
                                "b=.25x",
                                "b=.125x")

save(bb_comp_nopid, bb_comp_mf_nopid, bb_comp, bb_comp_mf, 
     file = paste0(path_save, "biasbound_comp.Rdata"))


############################# Save all but K and svdK #####################
kbal_b2x_s <- kbal_b2x_f[!(names(kbal_b2x_f) %in% c("K", "svdK"))]
kbal_b1x_s <- kbal_b1x_f[!(names(kbal_b1x_f) %in% c("K", "svdK"))]
kbal_b.5x_s <- kbal_b.5x_f[!(names(kbal_b.5x_f) %in% c("K", "svdK"))]
kbal_b.25x_s <- kbal_b.25x_f[!(names(kbal_b.25x_f) %in% c("K", "svdK"))]
kbal_b.125x_s <- kbal_b.125x_f[!(names(kbal_b.125x_f) %in% c("K", "svdK"))]

kbal_mf_b2x_s <- kbal_mf_b2x_f[!(names(kbal_mf_b2x_f) %in% c("K", "svdK"))]
kbal_mf_b1x_s <- kbal_mf_b1x_f[!(names(kbal_mf_b1x_f) %in% c("K", "svdK"))]
kbal_mf_b.5x_s <- kbal_mf_b.5x_f[!(names(kbal_mf_b.5x_f) %in% c("K", "svdK"))]
kbal_mf_b.25x_s <- kbal_mf_b.25x_f[!(names(kbal_mf_b.25x_f) %in% c("K", "svdK"))]
kbal_mf_b.125x_s <- kbal_mf_b.125x_f[!(names(kbal_mf_b.125x_f) %in% c("K", "svdK"))]


kbal_nopid_b2x_s <- kbal_nopid_b2x_f[!(names(kbal_nopid_b2x_f) %in% c("K", "svdK"))]
kbal_nopid_b1x_s <- kbal_nopid_b1x_f[!(names(kbal_nopid_b1x_f) %in% c("K", "svdK"))]
kbal_nopid_b.5x_s <- kbal_nopid_b.5x_f[!(names(kbal_nopid_b.5x_f) %in% c("K", "svdK"))]
kbal_nopid_b.25x_s <- kbal_nopid_b.25x_f[!(names(kbal_nopid_b.25x_f) %in% c("K", "svdK"))]
kbal_nopid_b.125x_s <- kbal_nopid_b.125x_f[!(names(kbal_nopid_b.125x_f) %in% c("K", "svdK"))]

kbal_mf_nopid_b2x_s <- kbal_mf_nopid_b2x_f[!(names(kbal_mf_nopid_b2x_f) %in% c("K", "svdK"))]
kbal_mf_nopid_b1x_s <- kbal_mf_nopid_b1x_f[!(names(kbal_mf_nopid_b1x_f) %in% c("K", "svdK"))]
kbal_mf_nopid_b.5x_s <- kbal_mf_nopid_b.5x_f[!(names(kbal_mf_nopid_b.5x_f) %in% c("K", "svdK"))]
kbal_mf_nopid_b.25x_s <- kbal_mf_nopid_b.25x_f[!(names(kbal_mf_nopid_b.25x_f) %in% c("K", "svdK"))]
kbal_mf_nopid_b.125x_s <- kbal_mf_nopid_b.125x_f[!(names(kbal_mf_nopid_b.125x_f) %in% c("K", "svdK"))]


save( kbal_b2x_s,kbal_b1x_s,kbal_b.5x_s,kbal_b.25x_s,kbal_b.125x_s,
      kbal_mf_b2x_s,kbal_mf_b1x_s,kbal_mf_b.5x_s, kbal_mf_b.25x_s, kbal_mf_b.125x_s,
      kbal_nopid_b2x_s,kbal_nopid_b1x_s,
      kbal_nopid_b.5x_s, kbal_nopid_b.25x_s, kbal_nopid_b.125x_s,
      kbal_mf_nopid_b2x_s, kbal_mf_nopid_b1x_s, kbal_mf_nopid_b.5x_s, 
      kbal_mf_nopid_b.25x_s,kbal_mf_nopid_b.125x_s,
      file = paste0(path_save, "kbal_out_noKsvd.Rdata"))



############################# Survey Designs #######################

###############  (1) KPOP + PID ##############
## (1.1) KPOP + MF =FALSE + W PID
kbal_wt_b2x <- svydesign(~1, data = pew,
                         weights = kbal_b2x_f$w[kbal_data_sampled ==1])
kbal_wt_b1x <- svydesign(~1, data = pew,
                         weights = kbal_b1x_f$w[kbal_data_sampled ==1])
kbal_wt_b.5x <- svydesign(~1, data = pew,
                          weights = kbal_b.5x_f$w[kbal_data_sampled ==1])
kbal_wt_b.25x <- svydesign(~1, data = pew,
                           weights = kbal_b.25x_f$w[kbal_data_sampled ==1])
kbal_wt_b.125x <- svydesign(~1, data = pew,
                            weights = kbal_b.125x_f$w[kbal_data_sampled ==1])

## (1.2) KPOP + MF = TRUE + W PID
kbal_mf_wt_b2x <- svydesign(~1, data = pew,
                            weights = kbal_mf_b2x_f$w[kbal_data_sampled ==1])
kbal_mf_wt_b1x <- svydesign(~1, data = pew,
                            weights = kbal_mf_b1x_f$w[kbal_data_sampled ==1])
kbal_mf_wt_b.5x <- svydesign(~1, data = pew,
                             weights = kbal_mf_b.5x_f$w[kbal_data_sampled ==1])
kbal_mf_wt_b.25x <- svydesign(~1, data = pew,
                              weights = kbal_mf_b.25x_f$w[kbal_data_sampled ==1])
kbal_mf_wt_b.125x <- svydesign(~1, data = pew,
                               weights = kbal_mf_b.125x_f$w[kbal_data_sampled ==1])

################### (2) KPOP + NO PID ##############

######## (2.1) KPOP + MF = FALSE + NO PID
kbal_wt_b2x_nopid <- svydesign(~1, data = pew,
                               weights = kbal_nopid_b2x_f$w[kbal_data_sampled ==1])
kbal_wt_b1x_nopid <- svydesign(~1, data = pew,
                               weights = kbal_nopid_b1x_f$w[kbal_data_sampled ==1])
kbal_wt_b.5x_nopid <- svydesign(~1, data = pew,
                                weights = kbal_nopid_b.5x_f$w[kbal_data_sampled ==1])
kbal_wt_b.25x_nopid <- svydesign(~1, data = pew,
                                 weights = kbal_nopid_b.25x_f$w[kbal_data_sampled ==1])
kbal_wt_b.125x_nopid <- svydesign(~1, data = pew,
                                  weights = kbal_nopid_b.125x_f$w[kbal_data_sampled ==1])

############### (2.2) KPOP + MF=TRUE + NO PID
kbal_wt_mf_nopid_b2x <- svydesign(~1, data = pew,
                                  weights = kbal_mf_nopid_b2x_f$w[kbal_data_sampled ==1])
kbal_wt_mf_nopid_b1x <- svydesign(~1, data = pew,
                                  weights = kbal_mf_nopid_b1x_f$w[kbal_data_sampled ==1])
kbal_wt_mf_nopid_b.5x <- svydesign(~1, data = pew,
                                 weights = kbal_mf_nopid_b.5x_f$w[kbal_data_sampled ==1])
kbal_wt_mf_nopid_b.25x <- svydesign(~1, data = pew,
                            weights = kbal_mf_nopid_b.25x_f$w[kbal_data_sampled ==1])
kbal_wt_mf_nopid_b.125x <- svydesign(~1, data = pew,
                            weights = kbal_mf_nopid_b.125x_f$w[kbal_data_sampled ==1])



################################ Summarize Results ##########################
## Presidential vote estimates
### Actual results
pres <- readRDS(paste0(path_data,"election.rds"))
natl_margin <- pres %>%
    summarise(margin = (sum(demtotal) - sum(reptotal)) /
                  (sum(demtotal) + sum(reptotal))) %>%
    as.numeric()


comp_df <- data.frame(
    CCES = svycontrast(svymean(~recode_vote_2016, 
                               cces_awt, na.rm = TRUE),
                       vote_contrast),
    Pew_unweighted = svycontrast(svymean(~recode_vote_2016, 
                                     pew_srs, na.rm = TRUE),
                             vote_contrast),
    #with PID
    #b = 2
    Pew_kbal_b2x = svycontrast(svymean(~recode_vote_2016, kbal_wt_b2x,
                                        na.rm = TRUE), vote_contrast),
    Pew_kbal_mf_b2x = svycontrast(svymean(~recode_vote_2016, kbal_mf_wt_b2x,
                                       na.rm = TRUE), vote_contrast),
    #b = 1
    Pew_kbal_b1x = svycontrast(svymean(~recode_vote_2016, kbal_wt_b1x,
                                       na.rm = TRUE), vote_contrast),
    Pew_kbal_mf_b1x = svycontrast(svymean(~recode_vote_2016, kbal_mf_wt_b1x,
                                          na.rm = TRUE), vote_contrast),

    #b = 0.5
    Pew_kbal_b.5x = svycontrast(svymean(~recode_vote_2016, kbal_wt_b.5x,
                                        na.rm = TRUE), vote_contrast),
    Pew_kbal_mf_b.5x = svycontrast(svymean(~recode_vote_2016, kbal_mf_wt_b.5x,
                                          na.rm = TRUE), vote_contrast),

    #b = 0.25
    Pew_kbal_b.25x = svycontrast(svymean(~recode_vote_2016, kbal_wt_b.25x,
                                         na.rm = TRUE), vote_contrast),
    Pew_kbal_mf_b.25x = svycontrast(svymean(~recode_vote_2016, kbal_mf_wt_b.25x,
                                          na.rm = TRUE), vote_contrast),

    #b = 0.125
    Pew_kbal_b.125x = svycontrast(svymean(~recode_vote_2016, kbal_wt_b.125x,
                                          na.rm = TRUE), vote_contrast),
    Pew_kbal_mf_b.125x = svycontrast(svymean(~recode_vote_2016, kbal_mf_wt_b.125x,
                                          na.rm = TRUE), vote_contrast),
    
    #without pid
    #b = 2
    Pew_kbal_b2x_nopid = svycontrast(svymean(~recode_vote_2016,
                                             kbal_wt_b2x_nopid, na.rm = TRUE),
                                     vote_contrast),
    Pew_kbal_mf_b2x_nopid = svycontrast(svymean(~recode_vote_2016,
                                                kbal_wt_mf_nopid_b2x, na.rm = TRUE),
                                        vote_contrast),
    #b=1
    Pew_kbal_b1x_nopid = svycontrast(svymean(~recode_vote_2016,
                                             kbal_wt_b1x_nopid, na.rm = TRUE),
                                     vote_contrast),
    Pew_kbal_mf_b1x_nopid = svycontrast(svymean(~recode_vote_2016,
                                                kbal_wt_mf_nopid_b1x, na.rm = TRUE),
                                        vote_contrast),
    #b=.5
    Pew_kbal_b.5x_nopid = svycontrast(svymean(~recode_vote_2016,
                                              kbal_wt_b.5x_nopid, na.rm = TRUE),
                                      vote_contrast),
    Pew_kbal_mf_b.5x_nopid = svycontrast(svymean(~recode_vote_2016,
                                                 kbal_wt_mf_nopid_b.5x, na.rm = TRUE),
                                         vote_contrast),
    #b=0.25
    Pew_kbal_b.25x_nopid = svycontrast(svymean(~recode_vote_2016,
                                               kbal_wt_b.25x_nopid, na.rm = TRUE),
                                       vote_contrast),
    Pew_kbal_mf_b.25x_nopid = svycontrast(svymean(~recode_vote_2016,
                                                  kbal_wt_mf_nopid_b.25x, na.rm = TRUE),
                                          vote_contrast),
    #b=.125
    Pew_kbal_b.125x_nopid = svycontrast(svymean(~recode_vote_2016,
                                                kbal_wt_b.125x_nopid, na.rm = TRUE),
                                        vote_contrast),
    Pew_kbal_mf_b.125x_nopid = svycontrast(svymean(~recode_vote_2016,
                                                   kbal_wt_mf_nopid_b.125x, na.rm = TRUE),
                                           vote_contrast) 
) %>%
    pivot_longer(cols = everything(),
                 names_to = c("source", ".value"),
                 names_pattern = "(.*)\\.(.*)") %>%
    rename(est = nlcon) %>%
    mutate(err = est - natl_margin,
           source = str_replace(source, "_", " ")) %>%
    mutate(err_target = est - est[source == "CCES"]) %>%
    mutate(source = factor(source,
                         levels = c("CCES",
                                    "Pew unweighted",

                                    "Pew kbal_b2x",
                                    "Pew kbal_mf_b2x",

                                    "Pew kbal_b1x",
                                    "Pew kbal_mf_b1x",

                                    "Pew kbal_b.5x",
                                    "Pew kbal_mf_b.5x",

                                    "Pew kbal_b.25x",
                                    "Pew kbal_mf_b.25x",

                                    "Pew kbal_b.125x",
                                    "Pew kbal_mf_b.125x",

                                    "Pew kbal_b2x_nopid",
                                    "Pew kbal_mf_b2x_nopid",

                                    "Pew kbal_b1x_nopid",
                                    "Pew kbal_mf_b1x_nopid",

                                    "Pew kbal_b.5x_nopid",
                                    "Pew kbal_mf_b.5x_nopid",

                                    "Pew kbal_b.25x_nopid",
                                    "Pew kbal_mf_b.25x_nopid",

                                    "Pew kbal_b.125x_nopid",
                                    "Pew kbal_mf_b.125x_nopid"),
                         labels = c("Target", "Unweighted",
                                    "Kpop b2x",
                                    "Kpop + MF +  b2x",
                                    
                                    "Kpop b1x",
                                    "Kpop + MF +  b1x",
                                    
                                    "Kpop b.5x",
                                    "Kpop + MF +  b.5x",
                                    
                                    "Kpop b.25x",
                                    "Kpop + MF +  b.25x",
                                    
                                    "Kpop b.125x",
                                    "Kpop + MF +  b.125x",
                                    
                                    "Kpop NoPID b2x",
                                    "kbal NoPID + MF b2x",
                                    
                                    "Kpop NoPID b1x",
                                    "Kpop NoPID + MF b1x",
                                    
                                    "Kpop NoPID b.5x",
                                    "Kpop NoPID + MF b2x",
                                    
                                    "Kpop NoPID b.25x",
                                    "Kpop NoPID + MF b.25x",
                                    
                                    "Kpop NoPID b.125x",
                                    "Kpop NoPID + MF b.125x")
                            ))
comp_df$PID <- NA
comp_df$PID[grep("PID", levels(comp_df$source))] <- "No PID"
comp_df$PID[-grep("PID", levels(comp_df$source)) ]<- "PID"
comp_df$PID[1:2]<- "Not Weighted" #making CCES target different color

comp_df$MF <- NA
comp_df$MF[grep("MF", levels(comp_df$source))] <- "MF"
comp_df$MF[-grep("MF", levels(comp_df$source)) ]<- "No MF"
comp_df$MF[1:2]<- "Not Weighted"


comp_df$b <- NA
comp_df[grep("b", levels(comp_df$source)), "b"] <- c(rep(2,2), 
                                                     rep(1,2),
                                                     rep(.5, 2),
                                                     rep(.25,2),
                                                     rep(.125,2),
                                                     rep(2,2), 
                                                     rep(1,2),
                                                     rep(.5, 2),
                                                     rep(.25,2),
                                                     rep(.125,2))
comp_df$b[1:2]<- "Not Weighted"

comp_df
save(comp_df, file = paste0(path_save, "comp_df.Rdata"))

######################## Cleaned Plots #######################
comp_df %>%  #filter(PID != "No PID") %>%
    ggplot() +
    aes(x = source, y = est, ymin = est - 1.96*SE, ymax = est + 1.96*SE, 
        shape = as.factor(PID), color=as.factor(MF)) +

    geom_hline(yintercept = c(0, natl_margin, comp_df$est[comp_df$source == "Target"]),
               linetype = c("solid", "dashed", "longdash"),
               color = c("black", "gray", "black")) +
    geom_pointrange() +
    scale_y_continuous(breaks = c(natl_margin, seq(-.05, .1, .05)),
                       minor_breaks = NULL,
                       labels = scales::percent_format()) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) +

    labs(x = NULL, y = "Estimated Margin (95% CI)") +
    ggtitle("Estimates of Clinton National Popular Vote Margin (Full SVD)") +
    theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
    annotate(geom = "text", x = 33.25, y = natl_margin, label = "", hjust = -0.1, angle = -90, color = "gray") +
    annotate(geom = "text", x = 32.25, y = natl_margin, label = "True\nNational\nMargin", hjust = -0.1, angle = -90, color = "gray") +
    annotate(geom = "text", x = 32.35, y = comp_df$est[comp_df$source == "CCES\n(Target)"], label = "CCES Estimated\n  Margin", hjust = -0.1, angle = 90) +
    theme(legend.title = element_blank())

ggsave(paste0(path_save, "allruns.pdf"), width = 12, height = 6)
