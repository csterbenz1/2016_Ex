
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
kbal_est$K[1:10,1:10]
kbal_est$onehot_data[1:10, 1:10]
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
kbal_data[1:10, 1:3]
kbal_est$onehot_data[1:10, 1:3]

K <- makeK(kbal_data, b=2, useasbases = kbal_data_sampled,
           linkernel = FALSE, scale = FALSE)

raw_counts <- -log(K)
#check this:
K[1:3, 1:3]
exp(-sum((kbal_data[1, ] - kbal_data[2,])^2)/2)
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
b_maxvar
K <- makeK(kbal_data, b=b_maxvar*2, useasbases = kbal_data_sampled,
           linkernel = FALSE, scale = FALSE)
res$objective
var(as.vector(K))
#yep
kbal_est$K[1:5, 1:5]
K[1:5, 1:5]
exp(-sum((kbal_data[1, ] - kbal_data[2,])^2)/b_maxvar)
exp(-sum((kbal_data[1, ] - kbal_data[2,])^2)/(b_maxvar*2))

raw_counts[1:10]
K_b1 <- makeK(kbal_data, b=1,
           useasbases = kbal_data_sampled,
           linkernel = FALSE, scale = FALSE)
raw_counts_b1 <- -log(K_b1)
raw_counts_b1[1:10] == 2*raw_counts[1:10]

n_d_b1 <- data.frame(diff = c(raw_counts_b1)) %>% group_by(diff) %>% summarise(n())
res = optimize(var_K, n_d_b1, length(diag(K)),
               interval=c(0,2000), maximum=TRUE)

res$maximum
res$objective

kbal_est$K[1:5, 1:5]
K_b1_bestb <- makeK(kbal_data, b=res$maximum,
              useasbases = kbal_data_sampled,
              linkernel = FALSE, scale = FALSE)
K_b1_bestb[1:10, 1:10]

K_b1_bestb <- makeK(kbal_est$onehot_data, b=res$maximum,
                    useasbases = kbal_data_sampled,
                    linkernel = FALSE, scale = FALSE)
K_b1_bestb[1:10, 1:10]


### sqrt(2)
kbal_data_sqrt2 = kbal_est$onehot_data*sqrt(1/2)
K_bsqrt2 <- makeK(kbal_data_sqrt2, b=1,
                    useasbases = kbal_data_sampled,
                    linkernel = FALSE, scale = FALSE)
raw_counts_sqrt2 = -log(K_bsqrt2)
n_d_b1_sqrt <- data.frame(diff = c(raw_counts_sqrt2)) %>% group_by(diff) %>% summarise(n())
res_sqrt = optimize(var_K, n_d_b1_sqrt, length(diag(K)),
               interval=c(0,2000), maximum=TRUE)

res_sqrt$maximum
res_sqrt$objective

K_maxvar_sqrt = makeK(kbal_data_sqrt2, b=res_sqrt$maximum,
                      useasbases = kbal_data_sampled,
                      linkernel = FALSE, scale = FALSE)
K_maxvar_sqrt[1:10, 1:10]
kbal_est$K[1:10,1:10]

#### compare ests:

t(out$est)
#run app-res_21.Rmd 
#yep they all match
test$comp_df_diff
out$est$kpop.mean
out$est$kpop_demos.mean
out$est$kpop_demos_wedu.mean
out$est$kpop_all.mean

 


#makeK with b = 2 to compensate for double counts -> so that's single counts then we get maxvarK b for single counts
#but we have K with odubles so we need to multiply it by two:
#alternatively try K-1 wiht bmaxvar and it should be the same
K_1 <- makeK(kbal_data, b=1, useasbases = kbal_data_sampled,
             linkernel = FALSE, scale = FALSE)
#we have sqrt(2) factor in allx so that when we square it inside the the kernel, we divide double counts correctly


#so the following should be the same:
#1. K double counts with one hot encoding + b = 2 
#2. K with sqrt(2) doubel counts + b = 1
#3. K coming out of kbal() 


#2.
#so since we put sqrt 2 in data, the optimal b will be twice the size so to match the true run we need to 
#divide that b by two
K1_sqrt = makeK(kbal_data*sqrt(2), b=kbal_est$b/2, useasbases = kbal_data_sampled,
           linkernel = FALSE, scale = FALSE)
#3. this still works even though we we are double one hot encoidng making quadruple counts, we just increase the b by a factor of 2 as one can see
K2 = kbal(kbal_data, sampled = kbal_data_sampled, cat_data = T, linkernel = T)
K2$b
kbal_est$b
#1. xby 2 kbal_est$b bc the kbal_data here is already one hot encoded so it's redundantly one-hot encoding which means we get 4 counts and we need to multiply b by 2 to match the original output
K2_indep = makeK(kbal_data, b=2*kbal_est$b, useasbases = kbal_data_sampled,
                 linkernel = FALSE, scale = FALSE)


K1_sqrt[1:5, 1:5]

K2_indep[1:5, 1:5]
K2$K[1:5, 1:5]
kbal_est$K[1:5, 1:5]


#this is just me figuring out the double one hot encoding commented above
#why are they not one-hot encoding the same way
# kbal_est$onehot_data[1:10,1:2]
# K2$onehot_data[1:10,1:2]
# #but this corresponds to their b values being different
# K2$b
# kbal_est$b
# ncol(kbal_est$onehot_data)
# ncol(K2$onehot_data)
# #oh shit ofc bc i passed K2 alreadty one hot encoded data... so yeah it still works when you double everything that's good


#fix this and do everyhting without the double one hot:

kbal_data_raw <- bind_rows(pew %>% dplyr::select(recode_age_bucket,
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

#2.can't replicate teh sqrt(.5 ) thing here bc of the annoying one hot but see b_tinkering to see it work
one_hot_kbal_data = one_hot(kbal_data_raw)
K1_sqrt = makeK(one_hot_kbal_data, b=sqrt(0.5), useasbases = kbal_data_sampled,
                linkernel = FALSE, scale = FALSE)
K1_sqrt[1:5, 1:5]
#3. 
K2_corr = kbal(kbal_data_raw, sampled = kbal_data_sampled, cat_data = T, linkernel = T)
K2_corr$K[1:5, 1:5]
kbal_est$K[1:5, 1:5]
K2_corr$b
kbal_est$b


