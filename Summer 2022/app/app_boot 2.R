
### Packages
library(tidyverse)
library(survey)
library(parallel)
#devtools::install_github("csterbenz1/KBAL", ref = "cat_kernel") 
library(kbal)

if(detectCores() > 10) {
    path_data= "/home/csterbenz/Data/"
} else {
    path_data= "/Users/Ciara_1/Dropbox/kpop/Updated/application/data/" 
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
# kbal_data <- data.frame(lapply(kbal_data, as.factor))
# kbal_data <- model.matrix(~ ., kbal_data,
#                           contrasts.arg = lapply(kbal_data[,], contrasts, contrasts=FALSE))
# kbal_data <- kbal_data[, -1]
#colnames(kbal_data)
#nrow(kbal_data)

kbal_data_sampled <- c(rep(1, nrow(pew)), rep(0, nrow(cces[subset,])))



##### Original Run:
#wait did we run with cont age and a mixed kernel ?/?? ever??? shit no i think not? 
# I should probably do that
kbal <- kbal(allx=kbal_data,
             sampled = kbal_data_sampled,
             cat_data = TRUE,
             population.w = if(POPW) {cces$commonweight_vv_post} else {NULL},
             ebal.tol = tolerance,
             ebal.maxit = maxit,
             maxnumdims = 500,
             sampledinpop = FALSE,
             fullSVD = TRUE)

colnames(kbal_data)
colnames(kbal$onehot_data)
kbal$K[1:5, 1:5]

#im literally being so stupid.
#ok so the actual variance of K with double counts is 2x b or b?
#is 2x b?
#check, var of K with b
#var of K with 2b

#var of K without double counts
#imagine two categories, 4 units, 2 are the same 2 are diff on on, 2 are diff on both





#this is different than previous cluster run?? 
#b is the same
#numdims diff (127 now vs 159 previously)
#biasbound diff (totally diff scale) (0.00696 old vs 0.0112)
#biasbound orig also diff (0.0393 now vs 0.0446) 
# :((( wtffffff
#check if force 159 dims
kbal_test = kbal(allx=kbal_data,
                sampled = kbal_data_sampled,
                cat_data = TRUE,
                population.w = if(POPW) {cces$commonweight_vv_post} else {NULL},
                ebal.tol = tolerance,
                ebal.maxit = maxit,
                numdims = 159,
                sampledinpop = FALSE,
                fullSVD = TRUE)

K_right <- makeK(kbal_data, b=kbal$b, useasbases = kbal_data_sampled,
                 linkernel = FALSE, scale = FALSE)

K_right[1:5, 1:5]
#i absolutely dont understand what is going wrong and i want to die!

#ok yes so this is the right one
kbal_test2 = kbal(allx=kbal_data,
                  K = K_right,
                 sampled = kbal_data_sampled,
                 cat_data = TRUE,
                 population.w = if(POPW) {cces$commonweight_vv_post} else {NULL},
                 ebal.tol = tolerance,
                 ebal.maxit = maxit,
                 sampledinpop = FALSE,
                 fullSVD = TRUE)

kbal_wrong = kbal
kbal = kbal_test2

##### Estimate
vote_diff <- quote((recode_vote_2016Democrat - recode_vote_2016Republican) /(recode_vote_2016Democrat + recode_vote_2016Republican + recode_vote_2016Other)) 

org_kbalreal_est = svycontrast(svymean(~recode_vote_2016,
                              svydesign(~1, data = pew, 
                                        weights = kbal$w[kbal_data_sampled ==1]), 
                              na.rm = TRUE), vote_diff)*100
org_real_est
org_est =  svymean(~diff_cces_on_pew, 
                   svydesign(~1, data = pew, weights = kbal$w[kbal_data_sampled ==1]),
                   na.rm = TRUE)*100
org_est

#boot func:
boot_kbal <- function(kbal_out,
                      dat = NULL, 
                      sampled,
                      n_boots = 500, 
                      search_range = NULL, increment_by = 1, num_dims = NULL, 
                      maxit = 500, tolerance = 1e-6, convergence = TRUE) {

    boot_rep = replicate(n_boots, {
        
        #this is dangerous generally, but assume the bases are sampled units
        boot_samp = sort(sample(1:ncol(kbal_out$K), ncol(kbal_out$K), replace = TRUE))
        pop_rows = nrow(kbal_out$K) - ncol(kbal_out$K)
        boot_K = kbal_out$K[c(1:pop_rows, boot_samp + pop_rows), boot_samp]
        #dangerous bc might have sample rows first so just pass in
        #sampled = c(rep(0, pop_rows), rep(1, length(boot_samp)))        

        kbal_boot =  kbal::kbal(allx= dat,
                                K = boot_K,
                                sampled = sampled,
                                cat_data = TRUE,
                                ebal.convergence = convergence,
                                fullSVD = TRUE,
                                incrementby = increment_by,
                                population.w = if(POPW) {cces$commonweight_vv_post} else {NULL},
                                ebal.tol = tolerance,
                                ebal.maxit = maxit,
                                numdims = num_dims,
                                minnumdims = search_range[1],
                                maxnumdims = search_range[2],
                                sampledinpop = FALSE)
        
        weights = data.frame(boot_w = kbal_boot$w) 
        b_boot = NA
        estimates = data.frame(
            boot_real_est = 
                                   svycontrast(
            svymean(~recode_vote_2016,
                    svydesign(~1, 
                              data = pew, 
                              weights = kbal_boot$w[sampled ==1]), na.rm = TRUE), vote_diff)*100, 
            boot_est = svymean(~diff_cces_on_pew, 
                               svydesign(~1, data = pew, weights = kbal_boot$w[sampled ==1]),
                               na.rm = TRUE)*100)
        weights = cbind(iter = 0, weights)
        return(list(estimates = estimates, 
                    weights = weights, 
                    b_boot =  b_boot, 
                    boot_samp = boot_samp))
        
    }, simplify = F)
    
    estimates = lapply(boot_rep, `[[`, 1) %>% bind_rows()
    weights = lapply(boot_rep, `[[`, 2) %>% bind_rows()
    weights$iter = rep(c(1:n_boots), each = nrow(dat))
    b_boot = unlist(lapply(boot_rep, `[[`, 3))
    boot_samp =  unlist(lapply(boot_rep, `[[`, 4))
   
    return(list(estimates = estimates, 
                weights = weights, 
                b_boot =  b_boot, 
                boot_samp = boot_samp))
}



system.time({app_boot = boot_kbal(kbal_out = kbal, 
                                  dat= kbal_data, 
                                  sampled = kbal_data_sampled,
                                  n_boots = 2,
                                  num_dims = kbal$numdims)})
app_boot







############################# parallelize

cores_saved = 5
cores_used = detectCores()-cores_saved
n_boots = 2
nsims = cores_used*n_boots


##make params global for all cores
kbal_out = kbal
dat= kbal_data
sampled = kbal_data_sampled 
search_range = NULL
increment_by = 5
num_dims = NULL
maxit = 500
tolerance = 1e-6
convergence = TRUE
numdims = 159


parallel_boot_kbal <- mclapply(1:cores_used, function(core_i) {
    
    cat("###############################################################\n")
    cat(paste("===================== core:",core_i, 
              "===================== \n"))
    
    #each core will run replicate for n_boots
    boot_rep = replicate(n_boots, {
        
        #this is dangerous generally, but assume the bases are sampled units
        boot_samp = sort(sample(1:ncol(kbal_out$K), ncol(kbal_out$K), replace = TRUE))
        pop_rows = nrow(kbal_out$K) - ncol(kbal_out$K)
        boot_K = kbal_out$K[c(1:pop_rows, boot_samp + pop_rows), boot_samp]
        
        kbal_boot =  kbal::kbal(allx= dat,
                                K = boot_K,
                                sampled = sampled,
                                cat_data = TRUE,
                                ebal.convergence = convergence,
                                fullSVD = TRUE,
                                incrementby = increment_by,
                                population.w = if(POPW) {cces$commonweight_vv_post} else {NULL},
                                ebal.tol = tolerance,
                                ebal.maxit = maxit,
                                numdims = num_dims,
                                minnumdims = search_range[1],
                                maxnumdims = search_range[2],
                                sampledinpop = FALSE)
        
        weights = data.frame(boot_w = kbal_boot$w) 
        b_boot = NA
        estimates = data.frame(boot_real_est = 
                                   svycontrast(
                                       svymean(~diff_cess_on_pew,
                                               svydesign(~1, 
                                                         data = pew, 
                                                         weights = kbal_boot$w[sampled ==1]),
                                               na.rm = TRUE), vote_diff)*100, 
                               boot_est = svymean(~diff_cces_on_pew, 
                                                  svydesign(~1, data = pew, weights = kbal_boot$w[sampled ==1]),
                                                  na.rm = TRUE)*100)
        weights = cbind(iter = 0, weights)
        return(list(estimates = estimates, 
                    weights = weights, 
                    b_boot =  b_boot, 
                    boot_samp = boot_samp))
        
    }, simplify = F)
    
    estimates = lapply(boot_rep, `[[`, 1) %>% bind_rows()
    weights = lapply(boot_rep, `[[`, 2) %>% bind_rows()
    weights$iter = rep(c((1:n_boots) + n_boots*(core_i -1)), each = nrow(dat))
    b_boot = unlist(lapply(boot_rep, `[[`, 3))
    boot_samp =  unlist(lapply(boot_rep, `[[`, 4))
    
    return(list(core = core_i,
                estimates = estimates, 
                weights = weights, 
                b_boot =  b_boot, 
                boot_samp = boot_samp))
    
    
    
}, mc.cores = cores_used)





cores_saved = 5
cores_used = detectCores()-cores_saved
n_boots = 2
nsims = cores_used*n_boots

kbal_out = kbal
dat= kbal_data
sampled = kbal_data_sampled 
search_range = NULL
increment_by = 1
num_dims = NULL
maxit = 500
tolerance = 1e-6
convergence = TRUE
numdims = 159

parallel2_boot_kbal <- mclapply(1:nsims, function(sim_i) {
    
    cat("###############################################################\n")
    cat(paste("===================== SIM:",sim_i, 
              "===================== \n"))
    
    
    #this is dangerous generally, but assume the bases are sampled units
    boot_samp = sort(sample(1:ncol(kbal_out$K), ncol(kbal_out$K), replace = TRUE))
    pop_rows = nrow(kbal_out$K) - ncol(kbal_out$K)
    boot_K = kbal_out$K[c(1:pop_rows, boot_samp + pop_rows), boot_samp]
    
    kbal_boot =  kbal::kbal(allx= dat,
                            K = boot_K,
                            sampled = sampled,
                            cat_data = TRUE,
                            ebal.convergence = convergence,
                            fullSVD = TRUE,
                            incrementby = increment_by,
                            population.w = if(POPW) {cces$commonweight_vv_post} else {NULL},
                            ebal.tol = tolerance,
                            ebal.maxit = maxit,
                            numdims = num_dims,
                            minnumdims = search_range[1],
                            maxnumdims = search_range[2],
                            sampledinpop = FALSE)
    
    weights = data.frame(boot_w = kbal_boot$w) 
    b_boot = NA
    estimates = data.frame(boot_real_est = 
                               svycontrast(
                                   svymean(~diff_cess_on_pew,
                                           svydesign(~1, 
                                                     data = pew, 
                                                     weights = kbal_boot$w[sampled ==1]),
                                           na.rm = TRUE), vote_diff)*100, 
                           boot_est = svymean(~diff_cces_on_pew, 
                                              svydesign(~1, data = pew, weights = kbal_boot$w[sampled ==1]),
                                              na.rm = TRUE)*100)
    weights = cbind(iter = sim_i, weights)
    
    return(list(estimates = estimates, 
                weights = weights, 
                b_boot =  b_boot, 
                boot_samp = boot_samp))

}, mc.cores = cores_used)


#clean outside of mclappy:
estimates = lapply(parallel2_boot_kbal, `[[`, 1) %>% bind_rows()
weights = lapply(parallel2_boot_kbal, `[[`, 2) %>% bind_rows()
b_boot = unlist(lapply(parallel2_boot_kbal, `[[`, 3))
boot_samp =  unlist(lapply(parallel2_boot_kbal, `[[`, 4))

























boot_kbal <- function(kbal_out,
                      dat = NULL, 
                      sampled = NULL, 
                      n_boots = 500, 
                      search_range = NULL, increment_by = 1, num_dims = NULL, 
                      maxit = 500, tolerance = 1e-6, convergence = TRUE) {
    
    
    boot_rep = replicate(n_boots, {
        
        #this is dangerous generally, but assume the bases are sampled units
        boot_samp = sort(sample(1:ncol(kbal_out$K), ncol(kbal_out$K), replace = TRUE))
        pop_rows = nrow(kbal_out$K) - ncol(kbal_out$K)
        boot_K = kbal_out$K[c(1:pop_rows, boot_samp + pop_rows), boot_samp]
        sampled = c(rep(0, pop_rows), rep(1, length(boot_samp)))        
        
        kbal_boot =  kbal::kbal(allx= dat,
                                K = boot_K,
                                sampled = sampled,
                                cat_data = TRUE,
                                ebal.convergence = convergence,
                                fullSVD = TRUE,
                                incrementby = incremebt_by,
                                population.w = if(POPW) {cces$commonweight_vv_post} else {NULL},
                                ebal.tol = tolerance,
                                ebal.maxit = maxit,
                                numdims = num_dims,
                                minnumdims = search_range[1],
                                maxnumdims = search_range[2],
                                sampledinpop = FALSE)
        
        weights = data.frame(boot_w = kbal_boot$w) 
        b_boot = NA
        estimates = data.frame(boot_real_est = 
                                   svycontrast(
                                       svymean(~recode_vote_2016,
                                               svydesign(~1, 
                                                         data = pew, 
                                                         weights = kbal_boot$w[sampled ==1]), na.rm = TRUE), vote_diff)*100, 
                    boot_est = svymean(~diff_cces_on_pew, 
                                       svydesign(~1, data = pew, weights = kbal_boot$w[sampled ==1]),
                                       na.rm = TRUE)*100)
        weights = cbind(iter = 0, weights)
        return(list(estimates = estimates, 
                    weights = weights, 
                    b_boot =  b_boot, 
                    boot_samp = boot_samp))
        
    }, simplify = F)
    
    estimates = lapply(boot_rep, `[[`, 1) %>% bind_rows()
    weights = lapply(boot_rep, `[[`, 2) %>% bind_rows()
    weights$iter = rep(c(1:n_boots), each = nrow(dat))
    b_boot = unlist(lapply(boot_rep, `[[`, 3))
    boot_samp =  unlist(lapply(boot_rep, `[[`, 4))
    
    return(list(estimates = estimates, 
                weights = weights, 
                b_boot =  b_boot, 
                boot_samp = boot_samp))
}















