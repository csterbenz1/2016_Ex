##################### Mini Example in Paper #####################
#Summer 2022
#exploring how bad bootstrapping would be
suppressWarnings(library(tidyverse))

##########################################
##################### Set up ##############

####minimal unit version 4 +8
pop <- data.frame(
    republican =  c(rep(1,2), rep(0,2)),
    white = c(1,0,1,0),
    support = c(.8, rep(.2,3)))
unique(pop) # white republicans have higher support
table(pop[,1:2])/nrow(pop)
mean(pop$support)  # Target value 0.35

#sample build : correct margins/means, but wrong joint distribution
samp <- data.frame( republican = c(rep(1, 4), rep(0,4)),
                    white = c(rep(1,3), rep(0,1), rep(1,1), rep(0,3)),
                    support = c(rep(.8,3), rep(.2,5)))
table(samp[,1:2])/nrow(samp)



##########################################
################ Reweight #################

dat <- as.matrix(rbind(pop,samp))
X <- dat[,1:2]

#1k: sampled <- c(rep(0,800), rep(1,80))
sampled <- c(rep(0,4), rep(1,8))
#### Raking:
# Run ebal (treatment = population units = 1-sampled)
ebal_out <- kbal::ebalance_custom(Treatment = 1-sampled,
                                  X=X,
                                  constraint.tolerance=1e-6,
                                  print.level=-1)

#already have balanceon margins so have even weights
length(unique(ebal_out$w))
unique(cbind(samp, e_bal_weight = ebal_out$w))
#estimate does not change .425 as we expect
weighted.mean(samp[,3], w = ebal_out$w) 

########## Kbal

#NB: the way to force your df to be factors depends on your version of R, this is for 3.6
#if you have 4.0 try: replacing the first line with: 
#onehot_data <- data.frame(lapply(X, as.factor))
onehot_data <- data.frame(apply(X,2, as.factor))
onehot_data <- model.matrix(~ ., onehot_data,
                            contrasts.arg = lapply(onehot_data, contrasts, contrasts=FALSE))
onehot_data <- onehot_data[, -1]
colnames(onehot_data)

# K <- kbal::makeK(onehot_data, b=2, useasbases = rep(1, nrow(onehot_data)),
#            linkernel = FALSE, scale = FALSE)

#OOPS we should proabbly be using the same bases that we use in the kbal() run
K <- kbal::makeK(onehot_data, b=2, useasbases = sampled,
                 linkernel = FALSE, scale = FALSE)

raw_counts <- -log(K)

#2. now we find the b that will max the variance of this Kernel
var_K= function(b, n_d){
    p_d <- n_d/ sum(n_d) 
    d = 0:(length(n_d)-1)
    mean_k = sum(exp(-1*d/b)*p_d)
    var_k = sum((exp(-1*d/b)-mean_k)^2 * p_d)
    return(var_k)
}
#run optim:
#we need to know the frequences each number of "misses" and dplyr is much faster than table
n_d <- data.frame(diff = c(raw_counts)) %>% group_by(diff) %>% summarise(n())
n_d
#visually we can see these
hist(as.vector(raw_counts))

n_d <- as.vector(n_d[,2] %>% pull())
#if we don't want to count the diagonal we could go back just adjust accordingly:
n_d[1] <- n_d[1] - ncol(K)

res = optimize(var_K, n_d = n_d,
               interval=c(0,2000), maximum=TRUE)
b_maxvar <- res$maximum
b_maxvar
#yields this much variance in K
res$objective


#3. Run kbal with this choice of b
#Now: Kernel balancing for weighting to a population (i.e. kpop) -------
#updating with catkernel
kbalout = kbal::kbal(allx= X,
                     useasbases= sampled, #rep(1,nrow(X)) for all units
                     sampled = sampled,
                     cat_data = TRUE,
                     ebal.convergence = TRUE,
                     scale_data = FALSE,
                     drop_MC = FALSE,
                     b = b_maxvar*2,
                     fullSVD = TRUE,
                     sampledinpop = FALSE)
# The weights now vary:
plot(kbalout$w[sampled ==1], pch=16)

unique(c(kbalout$K))

# And produce correct estimate:
#update produces same correct estimate, and same weights (minus some arbitrary rounding)
weighted.mean(samp$support, w = kbalout$w[sampled==1])

#down weight non white republicans as we expect
unique(cbind(samp[,-3], k_bal_weight = kbalout$w[sampled==1]))


###################### first attempt at boostrap with small sample #######
counter = 0
boot_kbal <- function(kbal_out, dat = NULL, sampled = NULL, 
                      n_boots = 1000, search_range = NULL, increment_by = 1) {
    
    if(is.null(search_range)) {
        if(ncol(kbal_out$dist_record) < 30) { 
            search_range = c(1, kbal_out$dist_record["Dims", ncol(kbal_out$dist_record)])
        }
        else {
            search_range = c(kbal_out$numdims -15, kbal_out$numdims + 15)
        }
    }
    
    boot_rep = replicate(n_boots, {
        #method 1: most direct = sample units
        
        #method 2: sample... cols of K? no rows of K?
        #we want to resample ... oh god do we resample population units or just sample units???
        #im being fucking stupid. this is what happens when you take a month off
        #we want to resample... from the sample as if it were the pop... so we definitely need to resample sample units, but do we also need to resample pop units? no? not unless we want to incorperate the uncertainty in the pop which... there should be none, it's only bc we're using cces as an approximation of the us popualtion... so in theory actually resampling pop units would probably be a good idea with the cces but with this example it's stupid
        #ok so for now just resample sample units, which are rows of K:
        #OH GOD NO but the sample determines the bases.... so we have to also adjust the columns?
        #fuck wait but each kernel distance is the distance from that sample unit to that corresponding unit... wait so no the rows are pop + sample, and the columns... are... only sample
        #so we def have to reselect the columns, but we need to keep all the pop units and then resample the sample rows
        
        #this is dangerous generally, but assume the bases are sampled units
        boot_samp = sort(sample(1:ncol(kbal_out$K), ncol(kbal_out$K), replace = TRUE))
        pop_rows = nrow(kbal_out$K) - ncol(kbal_out$K)
        boot_K = kbal_out$K[c(1:pop_rows, boot_samp + pop_rows), boot_samp]
        X = dat
        # diag(boot_K)
        # #ok so this is now right with rows, but I now have to also resample the columns
        # boot_K[1:4,]
        # K[1:4, ]
        # boot_K[5:9,]
        # K[5,]
        # boot_K[10:12,]
        # K[9:10,]
        
        #ugh noooo getting all the arguments is such a fucking pain
        #boot with same b
        kbal_boot =  kbal::kbal(allx= X,
                                K = boot_K,
                                sampled = sampled,
                                cat_data = TRUE,
                                ebal.convergence = TRUE,
                                #scale_data = FALSE,
                                #drop_MC = FALSE,
                                #b = b_maxvar*2,
                                fullSVD = TRUE,
                                sampledinpop = FALSE)
        
        #boot w b search
        raw_counts <- -log(boot_K)
        n_d <- data.frame(diff = c(raw_counts)) %>% group_by(diff) %>% summarise(n())
        n_d <- as.vector(n_d[,2] %>% pull())
        #if we don't want to count the diagonal we could go back just adjust accordingly:
        n_d[1] <- n_d[1] - ncol(boot_K)
        
        res = optimize(var_K, n_d = n_d,
                       interval=c(0,2000), maximum=TRUE)
        
        # kbal_boot_b = kbal::kbal(allx= X,
        #                          K = boot_K,
        #                          sampled = sampled,
        #                          cat_data = TRUE,
        #                          ebal.convergence = TRUE,
        #                          scale_data = FALSE,
        #                          drop_MC = FALSE,
        #                          b = b_maxvar_boot*2,
        #                          fullSVD = TRUE,
        #                          sampledinpop = FALSE)
        
        
        weights = data.frame(boot_w = kbal_boot$w) #, boot_bootb = kbal_boot_b$w)
        b_boot = res$maximum
        estimates = data.frame(boot_est = weighted.mean(samp$support[boot_samp],
                                                        w = kbal_boot$w[sampled==1])) #, 
                               # boot_b_est = weighted.mean(samp$support[boot_samp],
                               #                            w = kbal_boot_b$w[sampled==1]))
        weights = cbind(iter = counter, weights)
        counter = counter + 1
        return(list(estimates = estimates, weights = weights, b_boot =  b_boot, 
                    boot_samp = boot_samp))
        
    }, simplify = F)
    
    
    estimates = lapply(boot_rep, `[[`, 1) %>% bind_rows()
    weights = lapply(boot_rep, `[[`, 2) %>% bind_rows()
    boot_samp =  unlist(lapply(boot_rep, `[[`, 4))
    #something is messed up here
    b_boot = test #try(unlist(lapply(test, `[[`, 3)))
    return(list(estimates = estimates, weights = weights, b_boot =  b_boot, 
                boot_samp = boot_samp))
}


test = boot_kbal(kbalout,n_boot = 500)
summary(test$estimates)
hist(test$estimates[,1])
var(test$estimates)


#stupid replicate doesn't modify past it's own scope so the iter doesn't work, go fix it
test$weights$iter = rep(c(1:500), each =  length(kbalout$w))
w_var = test$weights %>% group_by(iter) %>% summarise(w_var = var(boot_w))
#variance of variance of the weights....
var(w_var$w_var)












#######################################################################


############## 880 unit version and walking through more slowly:

#1000 unit version
pop <- data.frame(
    republican =  c(rep(1,400), rep(0,400)),
    white = c(rep(1,200), rep(0,200), rep(1,200), rep(0,200)),
    support = c(rep(.8,200), rep(.2,600)))
unique(pop) # white republicans have higher support
table(pop[,1:2])/nrow(pop)
mean(pop$support)  # Target value 0.35

#sample build : correct margins/means, but wrong joint distribution
# 100 unit version
samp <- data.frame( republican = c(rep(1, 40), rep(0,40)),
                    white = c(rep(1,30), rep(0,10), rep(1,10), rep(0,30)),
                    support = c(rep(.8,30), rep(.2,50)))
unique(samp)
table(samp[,1:2])/nrow(samp)
mean(samp$support) #we get .425 because white republicans over represented




###run kbal

dat <- as.matrix(rbind(pop,samp))
X <- dat[,1:2]

kbalout_2 = kbal::kbal(allx= X,
                     useasbases= c(rep(0, nrow(pop)), rep(1, nrow(samp))), 
                     sampled = c(rep(0, nrow(pop)), rep(1, nrow(samp))),
                     cat_data = TRUE,
                     ebal.convergence = TRUE,
                     scale_data = FALSE,
                     drop_MC = FALSE,
                     fullSVD = TRUE,
                     sampledinpop = FALSE)
# The weights now vary:
plot(kbalout_2$w[sampled ==1], pch=16)
plot(kbalout$w[sampled ==1], pch=16)

unique(c(kbalout$K))
unique(c(kbalout_2$K))
weighted.mean(samp$support, w = kbalout_2$w[c(rep(0, nrow(pop)), rep(1, nrow(samp)))==1])


test = boot_kbal(kbalout_2, dat = X, 
                 sampled =  c(rep(0, nrow(pop)), rep(1, nrow(samp))), 
                 n_boots = 500)
#check re-sampling dist
hist(test$boot_samp)

summary(test$estimates)
hist(test$estimates[,1])
#well shit???? why is this so good??
table(round(test$estimates[,1], 5))
var(test$estimates)

test$weights$iter = rep(c(1:500), each =  nrow(X))
w_var = test$weights %>% group_by(iter) %>% summarise(w_var = var(boot_w))
#why is this so low??
hist(w_var$w_var)
#variance of variance of the weights....



#why is this so good??? I should check that this is the same as actually just rebuilding K every time and doing a full bootstrap




boot_kbal_full <- function(pop, samp,
                           samp_outcome,
                      n_boots = 1000,
                      search_range = NULL,
                      increment_by = 1) {
    
    boot_rep = replicate(n_boots, {
        #method 1: most direct = sample units
        boot_samp = sort(sample(1:nrow(samp), nrow(samp), replace = TRUE))
        boot_dat = rbind(pop, samp[boot_samp,])
        sampled = c(rep(0, nrow(pop)), rep(1, nrow(samp)))
        
        kbal_boot =  kbal::kbal(allx= boot_dat,
                                sampled = sampled,
                                cat_data = TRUE,
                                ebal.convergence = TRUE,
                                fullSVD = TRUE,
                                incrementby = increment_by,
                                minnumdims = search_range[1],
                                maxnumdims = search_range[2],
                                sampledinpop = FALSE)
        
        weights = data.frame(boot_w = kbal_boot$w)
        b_boot = kbal_boot$b
        estimates = data.frame(boot_est = weighted.mean(samp_outcome[boot_samp],
                                                        w = kbal_boot$w[sampled==1]))
        weights = cbind(iter = counter, weights)
        counter = counter + 1
        return(list(estimates = estimates, weights = weights, b_boot =  b_boot, 
                    boot_samp = boot_samp))
        
    }, simplify = F)
    
    
    estimates = lapply(boot_rep, `[[`, 1) %>% bind_rows()
    weights = lapply(boot_rep, `[[`, 2) %>% bind_rows()
    weights$iter = rep(c(1:n_boots), each =  nrow(pop) + nrow(samp))
    b_boot = unlist(lapply(boot_rep, `[[`, 3))
    boot_samp =  unlist(lapply(boot_rep, `[[`, 4))
    
    return(list(estimates = estimates, 
                weights = weights, 
                b_boot =  b_boot, 
                boot_samp = boot_samp))
}



system.time({test_full = boot_kbal_full(pop = pop[,1:2], 
                           samp = samp[,1:2],
                           samp_outcome = samp[,3],
                           n_boots = 500)})


#check re-sampling dist
hist(test_full$boot_samp)

#variance in the estimates
summary(test_full$estimates)
hist(test_full$estimates[,1])
table(round(test_full$estimates[,1], 5))
var(test_full$estimates)

#variance in the weights
w_var = test_full$weights %>% group_by(iter) %>% summarise(w_var = var(boot_w))
hist(w_var$w_var)
#variance of variance of the weights
var(w_var$w_var)


#variance in best b (max var)
summary(test_full$b_boot)
var(test_full$b_boot)


########################################################################
########################################################################
#now that i've confirmed it's the same to resample the rows in K I'll try boostrapping the application


