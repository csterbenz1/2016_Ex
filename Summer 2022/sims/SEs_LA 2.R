




plot_lasso_margin_r5 
margin_sim = mean(outcome)*100
table_r5
la = plot_lasso_margin_r5 %>% 
    mutate(estimator_name = gsub("\n", " ", estimator_name)) %>%
    group_by(estimator_name) %>%
    mutate(avg_bias = mean(margin - margin_sim)
    ) %>% ungroup()
nrow(la)
#check that worked
la %>% group_by(estimator_name) %>% summarise(n_distinct(avg_bias))


la = la %>% mutate(bias_adj = margin + sign(avg_bias)*avg_bias )
#check that worked
la %>% group_by(estimator) %>% summarise(check = round(mean(bias_adj - margin_sim), 6))
#now pass to coverage:
all_SE_coverage <- function(sims, adjust_bias = FALSE, drop_NA = F, truth = NULL, methods = c("rake|kpop")) {

    est <- lapply(sims, `[[`, 1) %>% bind_rows()
    if(bias_adj) {
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


SE_coverage_bias_adj = all_SE_coverage(sims[good], adjust_bias = T, truth = mean(outcome), drop_NA = T)
SE_coverage = all_SE_coverage(sims[good], adjust_bias = F, truth = mean(outcome), drop_NA = T)
SE_coverage_bias_adj$coverage
SE_coverage$coverage
diff = SE_coverage_bias_adj$coverage - SE_coverage$coverage
#so they almost all increase as expected which is good, not by a huge lot tho
#so the bias is constributing but coverage is sill undercoverage on kpop all and kpop mf

diff





####
r5_emp_SEs = suppressWarnings(empirical_SEs(sims = sims[good], na_rm = T))
# load("~/Documents/Cloud Documents/Hazlett:Hartman RA/2016 Election/2016_Ex/Summer 2022/sims/noscale_kpopTRUE_noise1_on2023-02-01_nsims468.RData")
load("~/Documents/Cloud Documents/Hazlett:Hartman RA/2016 Election/2016_Ex/Summer 2022/sims/noscale_kpopTRUE_noise1_on2023-02-15_nsims494.RData")
SEs <- lapply(sims, `[[`, 2) %>% bind_rows()
methods = "rake|kpop|unweighted"
SEs = SEs[grepl(methods, colnames(SEs))]
colnames(SEs)
r5_emp_SEs = suppressWarnings(empirical_SEs(sims = sims[good], na_rm = T))
r5_emp_SEs$emp_SEs
#boot_SE = t(as.matrix(apply(est, 2, sd)))
#SE_boot = boot_SE[, colnames(boot_SE) %in% colnames(avg_SE_out)]
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

#let's jsut see if the're in the right order
# sum(SE_fixed$estimator != SE_linear$estimator)
# sum(SE_linear$estimator != SE_quasi$estimator)
# sum(SE_quasi$estimator != SE_chad$estimator)
tot = cbind(SE_fixed, SE_linear[,2], SE_quasi[,2], SE_chad[,2])
nrow(tot)
head(tot)
cols = grepl("kpop|unweighted|rake", colnames(est))
SE_boot =  t(as.matrix(apply(est, 2, sd)))
SE_boot = SE_boot[, cols]
Boot = data.frame(estimator =names(SE_boot),
                  SE_Boot = c(SE_boot))
Boot
# ggplot(tot %>% filter(grepl("kpop_all|mf|truth", estimator))) + 
#     geom_histogram(aes(x = SE_fixed, color = estimator, fill = estimator), alpha = .4 ) + 
#     #geom_density(aes(x = SE_fixed, color = estimator) ) + 
#     geom_vline(aes(xintercept =SE_Boot, color = estimator), 
#                data = Boot %>% filter(grepl("kpop_all|mf|truth", estimator))) + 
#     ggtitle("Distribution of Fixed SE Estimates by Method") +
#     theme_bw()

# method ="kpop_all|mf|truth"
# g1 = ggplot(tot %>% filter(grepl(method, estimator))) + 
#     #geom_histogram(aes(x = SE_fixed, color = estimator, fill = estimator), alpha = .4 ) + 
#     geom_density(aes(x = SE_fixed, color=estimator), alpha = .2 ) + 
#     geom_vline(aes(xintercept =SE_Boot, color = estimator), 
#                data = Boot %>% filter(grepl(method, estimator))) + 
#     ggtitle("Distribution of Fixed SE Estimates by Method") +
#     theme_bw()
# 
# g2 = ggplot(tot%>% filter(grepl(method, estimator))) + 
#     geom_density(aes(x = SE_linear, color = estimator) ) + 
#     geom_vline(aes(xintercept =SE_Boot, color = estimator), 
#                data = Boot %>% filter(grepl(method, estimator))) + 
#     ggtitle("Distribution of Lineariation SE Estimates by Method") +
#     theme_bw()
# 
# 
# g3 = ggplot(tot%>% filter(grepl(method, estimator))) + 
#     geom_density(aes(x = SE_quasi, color = estimator) ) + 
#     geom_vline(aes(xintercept =SE_Boot, color = estimator), 
#                data = Boot %>% filter(grepl(method, estimator))) + 
#     ggtitle("Distribution of Quasi-Poisson SE Estimates by Method") +
#     theme_bw()
# 
# g4 = ggplot(tot%>% filter(grepl(method, estimator))) + 
#     geom_density(aes(x = SE_chad, color = estimator) ) + 
#     geom_vline(aes(xintercept =SE_Boot, color = estimator),
#                data = Boot %>% filter(grepl(method, estimator))) + 
#     ggtitle("Distribution of Chad's SE Estimates by Method") +
#     theme_bw()


method ="kpop_all|mf|truth"
sub = tot %>% filter(grepl(method, estimator))
# gg_t = ggplot(sub) +
#     geom_density(aes(x = SE_fixed, color = "Fixed", fill = "Fixed"), alpha = .1 ) + 
#     geom_density(aes(x = SE_linear, color = "Linearization", fill = "Linearization"), alpha = .1 ) + 
#     geom_density(aes(x = SE_quasi, color = "Quasi-Poisson", fill = "Quasi-Poisson"), alpha = .1 ) + 
#     geom_density(aes(x = SE_chad, color = "Chad", fill = "Chad"), alpha = .1) + 
#     geom_vline(aes(xintercept =SE_Boot),
#                 data = Boot %>% filter(grepl(method, estimator))) +
#     # geom_vline(aes(xintercept =SE_Boot, color = "Fixed"),
#     #            data = Boot %>% filter(grepl(method, estimator))) +
#     # geom_vline(aes(xintercept =SE_Boot, color = "Quasi-Poisson"),
#     #            data = Boot %>% filter(grepl(method, estimator)) ) +
#     # geom_vline(aes(xintercept =SE_Boot, color = "Linearization"),
#     #            data = Boot %>% filter(grepl(method, estimator)) ) +
#     #geom_vline(aes(xintercept =SE_Boot), data = Boot %>% filter(grepl(method, estimator))) + 
#     scale_fill_continuous(position = element_blank()) + 
#     theme_bw()  + 
#     xlab("SE Estimate") + 
#     facet_grid(cols = vars(estimator))
# gg_t


method ="kpop_all|mf|truth|rake_demos"
sub = tot %>% filter(grepl(method, estimator) | estimator == "kpop")
sub
sub2 = sub %>% pivot_longer(cols= c(2:5), names_to = "SE_type") %>% rename(SE = value)

gg_t2 = ggplot(sub2) + 
    geom_density(aes(x = SE, color = SE_type, fill = SE_type ), alpha = .2) +  
    geom_vline(aes(xintercept =SE_Boot),
               data = Boot %>% filter(grepl(method, estimator))) +
    annotate(geom = "label",
             x= .005,
             y=Inf, vjust = 1.4,
             color = "Black",
             label =  "Bootstrap SE") + 
    theme_bw()  + 
    xlab("SE Estimate") + 
    facet_grid(cols = vars(estimator))
gg_t2
                     




################################################
##### lambda 1se + auto selection of lambdas
r5_lambda_1se_file = "~/Documents/Cloud Documents/Hazlett:Hartman RA/2016 Election/2016_Ex/Summer 2022/sims/noscale_kpopTRUE_noise1_on2023-02-15_nsims494.RData"
# "~/Documents/Cloud Documents/Hazlett:Hartman RA/2016 Election/2016_Ex/Summer 2022/sims/noscale_kpopTRUE_noise1_on2023-02-15_nsims481.RData"
good = which(lapply(sims, function (x) return(class(x))) == "list")
length(good)
emp_SE_1se = suppressWarnings(empirical_SEs(sims = sims[good], na_rm = T))
emp_SE_1se$emp_SEs


###################### lambda 1se + manual selection of lambdas
load("~/Documents/Cloud Documents/Hazlett:Hartman RA/2016 Election/2016_Ex/Summer 2022/sims/res_kpopTRUElambdaminFALSEmanTRUE_noise1_on2023-02-15_nsims12.RData")

good = which(lapply(sims, function (x) return(class(x))) == "list")
#rm(SE_coverage)=
SE_coverage = all_SE_coverage(sims[good], adjust_bias = F, truth = mean(outcome), drop_NA = T)
SE_coverage_1se_man = SE_coverage$coverage
SE_coverage_1se_man
emp_SE_1se_man = suppressWarnings(empirical_SEs(sims = sims[good], na_rm = T))
emp_SE_1se_man$emp_SEs


####################### lambda min + auto
load("~/Documents/Cloud Documents/Hazlett:Hartman RA/2016 Election/2016_Ex/Summer 2022/sims/res_kpopTRUElambdaminTRUEmanFALSE_noise1_on2023-02-15_nsims10.RData")
good = which(lapply(sims, function (x) return(class(x))) == "list")
length(good)
#rm(SE_coverage)=
SE_coverage = all_SE_coverage(sims[good], adjust_bias = F, truth = mean(outcome), drop_NA = T)
SE_coverage_min = SE_coverage$coverage
SE_coverage_min
emp_SE_min = suppressWarnings(empirical_SEs(sims = sims[good], na_rm = T))
emp_SE_min$emp_SEs



