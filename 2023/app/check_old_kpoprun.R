
#Check against old results
####################################################################################################
# what the fuck is going on
#0. check the plots are displaying the right info
#yeah seems ok?
load("/Users/Ciara_1/Dropbox/kpop/Updated/application/weights/bmaxvar_cat_nodiag_lasso061021_wTRUE_m500_t1e-06_2021-07-07.Rdata")
path_weights = "/Users/Ciara_1/Dropbox/kpop/2023/application/weights/"
weights_file = "app_update_raw2023-04-06.Rdata"
load(paste0(path_weights,weights_file))

#1. look at weights now and before
orig_w = weights_cat[[1]]$kpop_w
new_w = out$weights$kpop_w

length(orig_w)
sum(orig_w != 1)
comp = cbind(orig = orig_w[1:2052], new = new_w) 
comp = as.data.frame(comp) %>% mutate(diff = orig- new)
#shit so yeah they look pretty different
summary(comp$diff)
mean(abs(comp$diff))

#2. check num dims
old_est = t(out_cat)
new_est = t(out$est)
old_est
new_est

#same
out_cat$bb_opt
out$est$bb

#same
out_cat$bb_orig
out$est$bb*out$est$bbr

#yep
out_cat$b
out$est$b_out*2
#yep
out$est$l1
out_cat$l1
#yep
out_cat$numdims
out$est$numdims
#yep
out_cat$bb_orig
out$est$bb*out$est$bbr


#2. check kernel is the same in both
load("/Users/Ciara_1/Dropbox/kpop/2023/application/weights/DEBUG_full_kbal_obj_2023-04-06.Rdata")
kbal_est$w[1:10]
new_w[1:10]
orig_w[1:10]


#old: code  from app_scriped_bmaxvar_only_2 directly
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
K <- makeK(kbal_data, b=b_maxvar*2, useasbases = kbal_data_sampled,
           linkernel = FALSE, scale = FALSE)
#yep
kbal_est$K[1:5, 1:5]
K[1:5, 1:5]


#### compare ests:

t(out$est)
#run app-res_21.Rmd 
#yep they all match
test$comp_df_diff
out$est$kpop.mean
out$est$kpop_demos.mean
out$est$kpop_demos_wedu.mean
out$est$kpop_all.mean

 