##################### functiions
library(tidyverse)
coverage <- function(SE, x_bar, truth, crit_val= qnorm(0.975)) {
    x_upper = x_bar +  (SE*crit_val)
    x_lower = x_bar - (SE*crit_val)
    contains_truth = matrix(NA, ncol = ncol(SE), nrow = 1)
    for(i in 1:ncol(x_upper)) {
        contains_truth[,i] = sum((truth <= x_upper[,i] & truth >= x_lower[,i]))/nrow(SE)
    }
    colnames(contains_truth) = colnames(x_bar)
    return(contains_truth)
}

# eval coverage of diff SEs
all_SE_coverage <- function(sims, adjust_bias = FALSE, drop_NA = F, truth = NULL, methods = c("rake|kpop")) {
    
    est <- lapply(sims, `[[`, 1) %>% bind_rows()
    if(adjust_bias) {
        temp = est[grepl(methods, colnames(est))]
        avg_bias = colMeans(temp) - truth
        # avg_bias
        # la %>% group_by(estimator) %>% summarise(unique(avg_bias)/100)
        # avg_bias
        # colnames(test)
        bias_adj = temp
        for(i in 1:ncol(temp)) {
            bias_adj[, i] = temp[,i] - avg_bias[i]
        }
        est[grepl(methods, colnames(est))] = bias_adj
    }
    
    SEs <- lapply(sims, `[[`, 2) %>% bind_rows()
    est_c = est[grepl(methods, colnames(est))]
    SEs = SEs[grepl(methods, colnames(SEs))]
    
    
    if(drop_NA) {
        n_drop = NULL
        coverage_out = NULL
        for(i in 1:ncol(est_c)) {
            # drop = which(is.na(est_c[,i]))
            # est_temp = est_c[-drop, ]
            est_temp = na.omit(est_c[,i])
            n_drop = c(n_drop, nrow(est) - length(est_temp))
            if(i ==1) {
                names(n_drop) = colnames(est_c)[i]
            } else {
                names(n_drop)[i] = colnames(est_c)[i]
            }
            SEs_temp = na.omit(SEs[grepl(paste0(colnames(est_c)[i],"_SE"), colnames(SEs))])
            
            SE_fixed = SEs_temp[grepl("SE_fixed$", colnames(SEs_temp))]
            SE_linear = SEs_temp[grepl("SE_linear$", colnames(SEs_temp))]
            SE_quasi = SEs_temp[grepl("SE_quasi$", colnames(SEs_temp))]
            SE_chad= SEs_temp[grepl("SE_chad$", colnames(SEs_temp))]
            # SE_svy= SEs_temp[grepl("SVY", colnames(SEs_temp))]
            # search = gsub("_se_SVY","", colnames(SE_svy))
            
            coverage_out = cbind(coverage_out, rbind(coverage(SE_fixed, est_temp, truth = truth),
                                                     coverage(SE_linear, est_temp, truth = truth),
                                                     coverage(SE_quasi, est_temp,truth = truth),
                                                     coverage(SE_chad,est_temp,truth = truth)))
            rownames(coverage_out) = c("SE_fixed", "SE_linear", "SE_quasi", "SE_chad")
            colnames(coverage_out)[i] = colnames(est_c)[i]
            
        }
    } else {
        est_c = est[grepl(methods, colnames(est))]
        SEs = SEs[grepl(methods, colnames(SEs))]
        SE_fixed = SEs[grepl("SE_fixed$", colnames(SEs))]
        SE_linear = SEs[grepl("SE_linear$", colnames(SEs))]
        SE_quasi = SEs[grepl("SE_quasi$", colnames(SEs))]
        SE_chad= SEs[grepl("SE_chad$", colnames(SEs))]
        
        SE_svy= SEs[grepl("SVY", colnames(SEs))]
        if(ncol(SE_svy) != 0){
            #just making sure we're getting the estimates for the same SEs that we output from svy obj which currently is demos_noedu
            search = gsub("_se_SVY","", colnames(SE_svy))
            grepl(search, colnames(est_c))
            s = coverage(SE_svy, truth = truth, est_c[,grepl(search, colnames(est_c))])
            s1 = rep(NA, ncol(SE_fixed))
            s1[grepl(search, colnames(est_c))] = s
            #colnames(s1) = colnames(SE_fixed)
            coverage_out = rbind(coverage(SE_fixed, est_c, truth = truth),
                                 coverage(SE_linear, est_c, truth = truth),
                                 coverage(SE_quasi, est_c,truth = truth),
                                 coverage(SE_chad,est_c,truth = truth), 
                                 s1)
            rownames(coverage_out) = c("SE_fixed", "SE_linear", "SE_quasi", "SE_chad", "SE_svy")
        } else {
            coverage_out = rbind(coverage(SE_fixed, est_c, truth = truth),
                                 coverage(SE_linear, est_c, truth = truth),
                                 coverage(SE_quasi, est_c,truth = truth),
                                 coverage(SE_chad,est_c,truth = truth))
            rownames(coverage_out) = c("SE_fixed", "SE_linear", "SE_quasi", "SE_chad")
        }
        
    }
    
    if(drop_NA) {
        out = list()
        out$n_drop = n_drop
        out$coverage = coverage_out
    } else {
        out = coverage_out
    }
    if(adjust_bias) {
        out$bias_adj_est = est
    }
    return(out)
}

empirical_SEs <- function(sims, eval_kpop = T, return_svy_package = F, na_rm = F) {
    SEs = lapply(sims, `[[`, 2) %>% bind_rows()
    est <- lapply(sims, `[[`, 1) %>% bind_rows()
    cols = grepl("kpop|rake|post_stratification|unweighted", colnames(est))
    #cols = if(eval_kpop) { c(3:7,10:12,14:19)} else {c(3:7,10:12)}
    est =est[, cols]
    
    #avg SEs
    avg_SE = colMeans(SEs, na.rm = na_rm)
    avg_SE_out = rbind(avg_SE[grepl("SE_fixed$", names(avg_SE))],
                       avg_SE[grepl("SE_linear$", names(avg_SE))],
                       avg_SE[grepl("SE_quasi$", names(avg_SE))],
                       avg_SE[grepl("SE_chad", names(avg_SE))],
                       avg_SE[grepl("SVY", names(avg_SE))])
    rownames(avg_SE_out) = c("SE_fixed", "SE_linear", "SE_quasi", "SE_chad", "SE_SVY") 
    colnames(avg_SE_out) = gsub("_SE_fixed", "", colnames(avg_SE_out))
    
    #bootstrapped SEs
    boot_SE = t(as.matrix(apply(est, 2, sd)))
    SE_boot = boot_SE[, colnames(boot_SE) %in% colnames(avg_SE_out)]
    #stupid fix for non names kpop mf Ses in some sims
    good_cnames = grepl("kpop|rake|post_stratification|unweighted", colnames(avg_SE_out))
    good_cnames = colnames(avg_SE_out)[good_cnames]
    emp_SEs = rbind(avg_SE_out, SE_boot[good_cnames])
    if(rownames(emp_SEs)[nrow(emp_SEs)] == "") { rownames(emp_SEs)[nrow(emp_SEs)] = "SE_boot" }
    
    if(!return_svy_package) {
        avg_SE_out = avg_SE_out[-which(rownames(avg_SE_out) == "SE_SVY"), ]
        emp_SEs = emp_SEs[-which(rownames(emp_SEs) == "SE_SVY"), ]
    }
    
    return(list(emp_SEs =emp_SEs, 
                boot_SE = boot_SE,
                avg_SE = avg_SE_out) )    
}



##############################################################################################################################
### Lambda 1se + auto choice of  lamdbas r^2

r5_file = "~/Documents/Cloud Documents/Hazlett:Hartman RA/2016 Election/2023 Data/noscale_kpopTRUE_noise1_on2023-02-15_nsims494.RData"
load(r5_file)
good = which(lapply(sims, function (x) return(class(x))) == "list")
#1. look at SEs:
est <- lapply(sims, `[[`, 1) %>% bind_rows()
SEs <- lapply(sims, `[[`, 2) %>% bind_rows()
methods = "rake|kpop|unweighted"
SEs = SEs[grepl(methods, colnames(SEs))]
colnames(SEs)
r5_emp_SEs = suppressWarnings(empirical_SEs(sims = sims[good], na_rm = T))
r5_emp_SEs$emp_SEs

SE_fixed = SEs[grepl("SE_fixed$", colnames(SEs))] %>% 
    pivot_longer(everything(),names_to = "estimator", 
                 values_to = "SE_fixed") %>% mutate(estimator = gsub("_SE_fixed", "", estimator))
SE_linear = SEs[grepl("SE_linear$", colnames(SEs))] %>% 
    pivot_longer(everything(),names_to = "estimator", 
                 values_to = "SE_linear") %>% mutate(estimator = gsub("_SE_linear", "", estimator))
SE_quasi = SEs[grepl("SE_quasi$", colnames(SEs))]%>% 
    pivot_longer(everything(),names_to = "estimator", 
                 values_to = "SE_quasi") %>% mutate(estimator = gsub("_SE_quasi", "", estimator))
SE_chad= SEs[grepl("SE_chad$", colnames(SEs))]%>% 
    pivot_longer(everything(),names_to = "estimator", 
                 values_to = "SE_chad") %>% mutate(estimator = gsub("_SE_chad", "", estimator))
tot = cbind(SE_fixed, SE_linear[,2], SE_quasi[,2], SE_chad[,2])


cols = grepl("kpop|unweighted|rake", colnames(est))
SE_boot =  t(as.matrix(apply(est, 2, sd)))
SE_boot = SE_boot[, cols]
Boot = data.frame(estimator=names(SE_boot),
                  SE_Boot = c(SE_boot))
Boot
method ="kpop_all|mf|truth"
sub = tot %>% filter(grepl(method, estimator) | estimator == "kpop")
sub2 = sub %>% pivot_longer(cols= c(2:5), names_to = "SE_type") %>% rename(SE = value)

gg_t2 = ggplot(sub2) + 
    geom_density(aes(x = SE, color = SE_type, fill = SE_type ), alpha = .2) +  
    geom_vline(aes(xintercept =SE_Boot),
               data = Boot %>% filter(grepl(method, estimator)| estimator == "kpop")) +
    annotate(geom = "label",
             x= .005,
             y=Inf, vjust = 1.4,
             color = "Black",
             parse = T,
             label =  "sd(Y_hat)") + 
    theme_bw()  + 
    xlab("SE Estimate") + 
    facet_grid(cols = vars(estimator))
gg_t2







#2. observing the large right tales, let's look at what the sample looks like and the weights look like for these large SEs
kpop_all_SE = SEs[, grepl("kpop_mf", colnames(SEs))]
apply(kpop_all_SE, 2, summary)
#lets get those in the top 75% percentile
#start by arranging this will be a bitch if we dont do it by column
kpop_mf_chad = kpop_all_SE$kpop_mf_SE_chad %>% sort()
#fuck this... does quantile impute the quantile? i guess
which(kpop_mf_chad == as.numeric(quantile(kpop_mf_chad,.75 )))
#fuck that let's just do 
length(kpop_mf_chad)*.75
#ok so now we want the 370:494 and we now want to go grab those units in teh unordered data
big_SE = which(kpop_all_SE$kpop_mf_SE_chad %in% kpop_mf_chad[370:494])




######### plot again by big vs small
####
SEs <- lapply(sims, `[[`, 2) %>% bind_rows()
SEs = SEs[big_SE,]
SE_fixed = SEs[grepl("SE_fixed$", colnames(SEs))] %>% 
    pivot_longer(everything(),names_to = "estimator", 
                 values_to = "SE_fixed") %>% mutate(estimator = gsub("_SE_fixed", "", estimator))
SE_linear = SEs[grepl("SE_linear$", colnames(SEs))] %>% 
    pivot_longer(everything(),names_to = "estimator", 
                 values_to = "SE_linear") %>% mutate(estimator = gsub("_SE_linear", "", estimator))
SE_quasi = SEs[grepl("SE_quasi$", colnames(SEs))]%>% 
    pivot_longer(everything(),names_to = "estimator", 
                 values_to = "SE_quasi") %>% mutate(estimator = gsub("_SE_quasi", "", estimator))
SE_chad= SEs[grepl("SE_chad$", colnames(SEs))]%>% 
    pivot_longer(everything(),names_to = "estimator", 
                 values_to = "SE_chad") %>% mutate(estimator = gsub("_SE_chad", "", estimator))
tot_big = cbind(SE_fixed, SE_linear[,2], SE_quasi[,2], SE_chad[,2])


SEs <- lapply(sims, `[[`, 2) %>% bind_rows()
SEs = SEs[-big_SE,]
SE_fixed = SEs[grepl("SE_fixed$", colnames(SEs))] %>% 
    pivot_longer(everything(),names_to = "estimator", 
                 values_to = "SE_fixed") %>% mutate(estimator = gsub("_SE_fixed", "", estimator))
SE_linear = SEs[grepl("SE_linear$", colnames(SEs))] %>% 
    pivot_longer(everything(),names_to = "estimator", 
                 values_to = "SE_linear") %>% mutate(estimator = gsub("_SE_linear", "", estimator))
SE_quasi = SEs[grepl("SE_quasi$", colnames(SEs))]%>% 
    pivot_longer(everything(),names_to = "estimator", 
                 values_to = "SE_quasi") %>% mutate(estimator = gsub("_SE_quasi", "", estimator))
SE_chad= SEs[grepl("SE_chad$", colnames(SEs))]%>% 
    pivot_longer(everything(),names_to = "estimator", 
                 values_to = "SE_chad") %>% mutate(estimator = gsub("_SE_chad", "", estimator))
tot_s = cbind(SE_fixed, SE_linear[,2], SE_quasi[,2], SE_chad[,2])

method ="kpop_all|mf|truth"
sub = tot_big %>% filter(grepl(method, estimator) | estimator == "kpop")
big = sub %>% pivot_longer(cols= c(2:5), names_to = "SE_type") %>% rename(SE = value)
sub = tot_s %>% filter(grepl(method, estimator) | estimator == "kpop")
small = sub %>% pivot_longer(cols= c(2:5), names_to = "SE_type") %>% rename(SE = value)

#if you want a closer look at those tails...
gg_big = ggplot(big) + 
    geom_density(aes(x = SE, color = SE_type, fill = SE_type ), alpha = .2) +  
    geom_vline(aes(xintercept =SE_Boot),
               data = Boot %>% filter(grepl(method, estimator)| estimator == "kpop")) +
    annotate(geom = "label",
             x= .005,
             y=Inf, vjust = 1.4,
             color = "Black",
             parse = T,
             label =  "sd(Y_hat)") + 
    theme_bw()  + 
    xlab("SE Estimate") + 
    facet_grid(cols = vars(estimator))
gg_small = ggplot(small) + 
    geom_density(aes(x = SE, color = SE_type, fill = SE_type ), alpha = .2) +  
    geom_vline(aes(xintercept =SE_Boot),
               data = Boot %>% filter(grepl(method, estimator)| estimator == "kpop")) +
    annotate(geom = "label",
             x= .005,
             y=Inf, vjust = 1.4,
             color = "Black",
             parse = T,
             label =  "sd(Y_hat)") + 
    theme_bw()  + 
    xlab("SE Estimate") + 
    facet_grid(cols = vars(estimator))



###### start by looking at the sample for those
b_SE = sims[big_SE]
s_SE = sims[-big_SE]
length(s_SE)
length(b_SE)

samp_check <- lapply(b_SE, `[[`, 5)
samp_check = samp_check %>% bind_rows()
sum(samp_check$bad_sample)
samp_check

s_samp_check <- lapply(s_SE, `[[`, 5)
s_samp_check= s_samp_check %>% bind_rows()
sum(s_samp_check$bad_sample)
s_samp_check

# <= 5%
mean(samp_check$check.leq_5pp)
mean(s_samp_check$check.leq_5pp)

summary(samp_check$check.leq_5pp)
summary(s_samp_check$check.leq_5pp)
#verrrrrry slightly worse samples but not remarkably


#what about by category, well shit i did not print that out


#could also check numdims
big_est = est[big_SE,]
s_est = est[-big_SE,]
#interesting on avg the bigger SEs have slkightlyu fewer numdims and lower variance in the choice
#that is the polar opposite of what id expect....bc if the weights have more variance, presumably we have more numdims?
#but it's the opposite, maybe weightsare  more constrained with more numdims> 
cbind(big_SE = colMeans(big_est[,grepl("numdims", colnames(est))]), s_SE = colMeans(s_est[,grepl("numdims", colnames(est))]))
cbind(big_SE = apply(big_est[,grepl("numdims", colnames(est))],2, sd),
      s_SE = apply(s_est[,grepl("numdims", colnames(est))], 2, sd))





#3. look at var(w) and max and min weights 
weights = lapply(sims, `[[`, 3) #%>% bind_rows()
b_weights = weights[big_SE]
s_weights = weights[-big_SE]


#mean summary stats
t = lapply(b_weights, function(x) data.frame(apply(x, 2, summary)))
t = t %>% bind_rows()
nrow(t)
min = colMeans(t[grepl("Min", rownames(t)),])
max = colMeans(t[grepl("Max", rownames(t)),])
mean = colMeans(t[grepl("Mean", rownames(t)),])
perc_25 = colMeans(t[grepl("1st", rownames(t)),])
perc_75 = colMeans(t[grepl("3rd", rownames(t)),])    
avg_w = rbind(min = min, perc_25 = perc_25, mean = mean, perc_75 = perc_75, max = max)
avg_w

t = lapply(s_weights, function(x) data.frame(apply(x, 2, summary)))
t = t %>% bind_rows()
nrow(t)
min = colMeans(t[grepl("Min", rownames(t)),])
max = colMeans(t[grepl("Max", rownames(t)),])
mean = colMeans(t[grepl("Mean", rownames(t)),])
perc_25 = colMeans(t[grepl("1st", rownames(t)),])
perc_75 = colMeans(t[grepl("3rd", rownames(t)),])    
avg_w_s = rbind(min = min, perc_25 = perc_25, mean = mean, perc_75 = perc_75, max = max)

#so yes the weights are bigger in the big Ses but not wildly so?
t(avg_w)
t(avg_w_s)
#bigger max, smaller min, smaller 25% and 75% 
#interesting that the larger weights have a larger 75%, though its in the hundreths so its fairly similar
t(avg_w) - t(avg_w_s)
#rounding helps: so min and 25% are slightly lower, but notable that max is in the tenths here
round(t(avg_w) - t(avg_w_s),3)

#lol ok but like we probably just want sd?? lol, indeed look slike sd of weights is bigger for big SEs
t = lapply(b_weights, function(x) apply(x, 2, sd)) %>% bind_rows()
s = lapply(s_weights, function(x) apply(x, 2, sd)) %>% bind_rows()
rbind(colMeans(s), colMeans(t))
#so yes the big SEs have weights that are in general bigger: by an avg of + below sds
colMeans(t)-  colMeans(s)

# as you'd expdct the weights are bigger, the variance of the weights are bigger... but the samples don't look as insane? maybe i didn't look closely enough there andshoudl go back and check taht again ?d 