
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
increment_by = 5

cores_saved = 12
cores_used = detectCores()-cores_saved
n_boots = 14
nsims = cores_used*n_boots
nsims
#numdims = 159
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
#kbal_data <- data.frame(lapply(kbal_data, as.factor))
# kbal_data <- model.matrix(~ ., kbal_data,
#                           contrasts.arg = lapply(kbal_data[,], contrasts, contrasts=FALSE))
# kbal_data <- kbal_data[, -1]
#colnames(kbal_data)
#nrow(kbal_data)

kbal_data_sampled <- c(rep(1, nrow(pew)), rep(0, nrow(cces[subset,])))


#still issue w internal kbal code but for now get around by passing in correct b
#saving time here bc i know the optimal dims are 159 don't need to research
kbal = kbal(allx=kbal_data,
                 sampled = kbal_data_sampled,
                 cat_data = TRUE,
                 numdims = 159,
                 population.w = if(POPW) {cces$commonweight_vv_post} else {NULL},
                 ebal.tol = tolerance,
                 incrementby = increment_by,
                 ebal.maxit = maxit,
                 sampledinpop = FALSE,
                 fullSVD = TRUE)
kbal$K[1:5, 1:5]
colnames(kbal_data)
colnames(kbal$onehot_data)


##### Estimate
vote_diff <- quote((recode_vote_2016Democrat - recode_vote_2016Republican) /(recode_vote_2016Democrat + recode_vote_2016Republican + recode_vote_2016Other)) 

org_real_est = svycontrast(svymean(~recode_vote_2016,
                              svydesign(~1, data = pew, 
                                        weights = kbal$w[kbal_data_sampled ==1]), 
                              na.rm = TRUE), vote_diff)*100
org_real_est
org_est =  svymean(~diff_cces_on_pew, 
                   svydesign(~1, data = pew, weights = kbal$w[kbal_data_sampled ==1]),
                   na.rm = TRUE)*100
org_est

############################# parallelize


kbal_out = kbal
dat = kbal_data
sampled = kbal_data_sampled
#based on previous run max was like 296, min was 126 or so
search_range = c(100, 340)
increment_by = 5
num_dims = NULL
maxit = 500
tolerance = 1e-6
convergence = TRUE


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
    dist_record = data.frame(t(kbal_boot$dist_record))
    ebal_status = dist_record[which(dist_record[,"Dims"] == kbal_boot$numdims), "Ebal.Convergence"]
    min_converged = dist_record[which.min(dist_record[dist_record$Ebal.Convergence ==1,"BiasBound"]), "Dims"]
    kbal_out = data.frame(numdims = kbal_boot$numdims, 
                          biasbound = kbal_boot$biasbound_opt,
                          biasbound_r = kbal_boot$biasbound_ratio,
                          l1_org = kbal_boot$L1_orig,
                          l1_opt = kbal_boot$L1_opt,
                          bb_orig = kbal_boot$biasbound_orig, 
                          ebal_status = ebal_status)
    b_boot = NA
    estimates = data.frame(boot_real_est = 
                               svycontrast(
                                   svymean(~recode_vote_2016,
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
                kbal_out = kbal_out,
                b_boot =  b_boot, 
                boot_samp = boot_samp))

}, mc.cores = cores_used)


save(parallel2_boot_kbal, file = "./fullobj_save.Rdata")
good = which(lapply(parallel2_boot_kbal, function (x) return(class(x))) == "list")
length(good)

#clean outside of mclappy:
estimates = lapply(parallel2_boot_kbal, `[[`, 1) %>% bind_rows()
weights = lapply(parallel2_boot_kbal, `[[`, 2) %>% bind_rows()
kbal_out = unlist(lapply(parallel2_boot_kbal, `[[`, 3))
kbal_out[1:7]
#ok so nsims X 7 should fix it?
fix = NULL
for(i in 1:length(good)) {
    fix = rbind(fix, kbal_out[(1+(i-1)*7):(i*7)])
}
kbal_out = fix
b_boot = unlist(lapply(parallel2_boot_kbal, `[[`, 4))
#also need to fix this: we have resampled 2052
boot_samp =  unlist(lapply(parallel2_boot_kbal, `[[`, 5))
#restarts to the next sample at 2052
boot_samp[2052:2054]
fix = NULL
for(i in 1:length(good)) {
    fix = cbind(fix, boot_samp[(1+(i-1)*2052):(i*2052)])
}
dim(fix)
#so now sims are cols
colnames(fix) = t(map_chr(1:length(good), function(i) paste0("bootrep_", i)))
boot_samp = fix
save(estimates,weights,kbal_out, b_boot, boot_samp, 
     file = paste0("./kbal_boot",POPW, "_m",maxit, "_t",tolerance, "_inc",
                   increment_by, #"mindims", min_num_dims,
                   Sys.Date(),
                   # str_sub(gsub("[[:space:]|[:punct:]]", "_",
                   #              gsub("[:alpha:]", "",
                   #                   Sys.time())),
                   #         start = 1, end = -3),
                   "_nbootreps", length(good),
                   ".RData"))




################################################################
###############################################################
#processing results from cluster:
setwd("/Users/Ciara_1/Documents/Cloud Documents/Hazlett:Hartman RA/2016 Election/Summer 2022/app")
load("./kbal_boot_popwTRUE_m500_t1e-06_incby5_2022-07-28_nbootreps252.Rdata")
orig_est
orig_real_est
colMeans(estimates)
apply(estimates, 2, var)

var_w = weights %>% group_by(iter) %>% summarise(var_w = var(boot_w))
var_w