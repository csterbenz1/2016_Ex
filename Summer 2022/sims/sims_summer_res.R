###### Sims out:
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

POPW = FALSE
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




################# Load results
setwd('/Users/Ciara_1/Documents/Cloud Documents/Hazlett:Hartman RA/2016 Election/2016_Ex/Summer 2022/sims')
####### SIMPLE DGP + no POPW + n= 500
load("./cat_sims_modeled_outcome_nodiag_FALSE_m500_t1e-04_inc5mindims2022-09-13_nsims494.RData")
plot = est
nrow(est)
#estimates
colMeans(est)

####### Complex DGPP + no POPW + n= 500




################## results plots
#### ported over from tex
if(POPW) {
    cces_svy <- svydesign(ids = ~1, weights = ~commonweight_vv_post, data = cces)
} else {
    cces_svy <- svydesign(ids = ~1, data = cces)
}

plot_lasso_margin <- plot %>% 
    dplyr::select(unweighted, 
                  rake_demos_noeduc,
                  rake_demos_weduc,
                  rake_all,
                  #post_stratification,
                  post_strat_reduc,
                  #post_strat_all,
                  rake_truth#,
                  #kpop, 
                  #kpop_conv,
                  #kpop_mf, 
                  #kpop_demos,
                  #kpop_demos_wedu,
                  #kpop_all
                  ) %>% 
    pivot_longer(everything(),
                 names_to = "estimator", 
                 values_to = "margin") %>%
    mutate(margin = margin * 100,
           estimator_name = factor(case_when(estimator == "kpop" ~ "kpop",
                                             estimator == "kpop_mf" ~ "kpop aMF (All)",
                                             # estimator == "kpop_conv" ~ "kpop Converged",
                                             estimator == "kpop_demos" ~ "kpop+MF:\n (Demos)",
                                             estimator == "kpop_demos_wedu" ~ "kpop+MF:\n (Demos+Edu)",
                                             estimator == "kpop_all" ~ "kpop+MF:\n (All)",
                                             estimator == "rake_demos_noeduc" ~ "Mean Calibration:\n (Demos)",
                                             estimator == "rake_demos_weduc" ~  "Mean Calibration:\n (Demos+Edu)",
                                             estimator == "rake_all" ~ "Mean Calibration:\n (All)",
                                             estimator == "rake_truth" ~ "Mean Calibration:\n True Selection\nModel",
                                             #estimator == "post_stratification" ~ "Post-Strat Prev",
                                             estimator == "post_strat_reduc" ~ "Post-Stratification:\n (Reduc)",
                                             #estimator == "post_strat_all" ~ "Post-Strat All",
                                             estimator == "unweighted" ~ "Unweighted"),
                                   levels = c("Unweighted", 
                                              "Mean Calibration:\n (Demos)",
                                              "Mean Calibration:\n (Demos+Edu)",
                                              "Mean Calibration:\n (All)",
                                              #  "Post-Strat Prev", 
                                              "Post-Stratification:\n (Reduc)", 
                                              # "Post-Strat All",
                                              "kpop",
                                              # "kpop Converged",
                                              #"kpop aMF (All)",
                                              "kpop+MF:\n (Demos)",
                                              "kpop+MF:\n (Demos+Edu)",
                                              "kpop+MF:\n (All)",
                                              "Mean Calibration:\n True Selection\nModel"
                                   )))

#target:
margin_sim = svymean(~diff_cces_on_cces, cces_svy)[1]* 100
#### Box Plot
ggplot(data = plot_lasso_margin,
       aes(x = estimator_name, y = margin)) +
    geom_boxplot(alpha = 0.2) +
    geom_hline(yintercept = margin_sim) +
    theme_bw() +
    # ggtitle(paste0("Simulation Results: lasso ",
    #                " ", "pop weights: ", POPW, "\n",length(good), " sims, b=maxvar")) +
    xlab("") +
    ylab("Modeled Vote Margin") +
    annotate(geom = "text", x = 0.85, y = margin_sim+0.25, size = 2.7, angle = 90,
             label = "True Target\nPopulation\nMargin", hjust = 0) +
    ggtitle(paste0("n_samp =", round(mean(est$n)), " simple DGP ", simple_selection_model)) +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))


### table
plot_lasso_margin %>% 
    mutate(estimator_name = gsub("\n", " ", estimator_name)) %>%
    group_by(estimator_name) %>%
    summarize(
        Bias = mean(margin - margin_sim),
        SE = sd(margin),
        MSE = mean((margin - margin_sim)^2)
    ) %>%
    mutate(
        Bias_Reduc = 1- Bias / Bias[estimator_name == "Unweighted"]
    ) %>%
    arrange(MSE)

### kable
knitr::kable(
    plot_lasso_margin %>% 
        mutate(estimator_name = gsub("\n", " ", estimator_name)) %>%
        group_by(estimator_name) %>%
        summarize(
            Bias = mean(margin - margin_sim),
            SE = sd(margin),
            MSE = mean((margin - margin_sim)^2)
        ) %>%
        mutate(
            Bias_Reduc = 1- Bias / Bias[estimator_name == "Unweighted"]
        ) %>%
        arrange(MSE), caption = "Simulation Results (numeric)",
    booktabs = T,
    format = "latex", digits = 2)


############################
# SEs
SEs
colMeans(SEs)
colnames(SEs)

#by SE type
SE_fixed = SEs[grepl("SE_fixed$", colnames(SEs))]
SE_linear = SEs[grepl("SE_linear$", colnames(SEs))]
SE_quasi = SEs[grepl("SE_quasi$", colnames(SEs))]
SE_chad= SEs[grepl("SE_chad$", colnames(SEs))]

colMeans(SE_fixed)

#by method
SE_rake_demos_noeduc = SEs[grepl("^rake_demos_noeduc", colnames(SEs))]
SE_rake_demos_weduc = SEs[grepl("^rake_demos_weduc", colnames(SEs))]
SE_rake_all = SEs[grepl("^rake_all", colnames(SEs))]
SE_ps = SEs[grepl("^post_strat", colnames(SEs))]
SE_ps_reduc = SEs[grepl("^post_strat_reduc", colnames(SEs))]
SE_ps_all = SEs[grepl("^post_strat_all", colnames(SEs))]
SE_rake_truth = SEs[grepl("^rake_truth", colnames(SEs))]

SE_kpop = SEs[grepl("^kpop_SE", colnames(SEs))]
SE_kpop_conv = SEs[grepl("^kpop_conv", colnames(SEs))]
SE_kpop_mf = SEs[grepl("^kpop_mf", colnames(SEs))]
SE_kpop_demos = SEs[grepl("^kpop_demos_kpop", colnames(SEs))]
SE_kpop_demos_wedu = SEs[grepl("^kpop_demos_wedu", colnames(SEs))]
SE_kpop_all = SEs[grepl("^kpop_dall", colnames(SEs))]


#compare seems off
colMeans(SE_kpop)
colMeans(SE_rake_demos_noeduc)
