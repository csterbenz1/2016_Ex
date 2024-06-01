### Packages
library(MASS)
library(tidyverse)
library(survey)
#devtools::install_github("csterbenz1/KBAL", ref = "cat_kernel")
library(kbal)
library(parallel)
library(knitr)
library(glmnet)
#library(tictoc)

###### SET PARAMS  ###############
set.seed(9345876)

if(detectCores() > 10) {
    path_data = "/nas/home/ciaras/kpop/data/"
    cores_saved = 10
    eval_kpop = T
} else if(detectCores() != 4) {
    path_data= "/Users/Ciara_1/Dropbox/kpop/Updated/application/data/"
    cores_saved = 3
} else {
    path_data= "/Users/Ciara/Dropbox/kpop/Updated/application/data/"
    cores_saved = 2
}
options(dplyr.print_max = 1e9)
#fit kpop and others to cces with populations weights?
POPW = FALSE
# to run with a linear kernel so it's way faster; UPDATE: errors catch this as mistake and prevent
TEST = FALSE 
#ebal tolerance and max iterations for kpop
tolerance = 1e-4
maxit = 500
#adjust these both for runtime
increment = 5
min_num_dims = NULL
max_num_dims = NULL
eval_kpop = TRUE
up_nonlin = 1
kpop_constraints = T
SAVE = TRUE #save .Rdata results?
#Need to adjust accordingly to machine for adequate number of sims
#for cluster
#nsims = (detectCores()-cores_saved)*14*2
#First test with 10 sims
#test runtime w/  1 iteration per core
#nsims = (detectCores()-cores_saved)*1
#full run nsims approx 1000
nsims = (detectCores()-cores_saved)*6
nsims

##### Central Params to adjust

#add bernoulli noise by drawing binary from p(S=1)?
bern = FALSE 
#sd(y)*noise; 1-> r^2 = .5; sqrt(2) -> r^2 = .33; 1/2*sqrt(2) -> r^2 = .66;
noise = 1 
#adjusts sample size by dividing p(S) by scalar pS_denom (i.e. pS = plogis(XBeta)/pS_denom)
pS_denom = 30
#use the manually specified range of lambdas in the ridge residualization or allow glmnet to choose internally?
manual_lambda = FALSE 
#T=lambda as that which minimizes cverror in residualization; F= 1 sd from min choice
lambda_min = FALSE 


###################### Formulas ################
formula_rake_demos_noeduc <- ~recode_age_bucket + recode_female + recode_race +
  recode_region + recode_pid_3way

#updated to include 6 way edu
formula_rake_demos_weduc <- ~recode_age_bucket + recode_female +
  recode_race + recode_region + recode_educ + recode_pid_3way

formula_rake_all_vars <- ~recode_age_bucket + recode_female +
  recode_race + recode_region + recode_pid_3way + recode_educ +
  recode_income_5way + recode_relig_6way + recode_born + recode_attndch_4way

create_targets <- function (target_design, target_formula) {
  target_mf <- model.frame(target_formula, model.frame(target_design))
  target_mm <- model.matrix(target_formula, target_mf)
  wts <- weights(target_design)

  return(colSums(target_mm * wts) / sum(wts))
}

### Post-stratification function
## For now assumes that strata variable is already created and in
## the data set and called "strata"
postStrat <- function(survey, pop_counts, pop_w_col, strata_pass, warn = T) {
  survey_counts <- survey %>%
    group_by(!!as.symbol(strata_pass)) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    mutate(w_survey = n / sum(n))

  pop_counts <- pop_counts %>%
    rename(w_pop = matches(pop_w_col))
  
  if(warn == T & nrow(survey_counts) !=  nrow(pop_counts)) {
      missing_strat = pop_counts[! (( pop_counts[, strata_pass]%>% pull()) %in% (survey_counts[, strata_pass]%>% pull() )), strata_pass]
      warning(paste("Strata in Pop not found in Sample. Dropping", 
                    sum(pop_counts[(pop_counts[, strata_pass] %>% pull()) %in% 
                                       (missing_strat %>% pull()),"n" ]), 
                    "empty cells\n"), immediate.  =T )
  } 
  post_strat <- pop_counts %>%
    left_join(survey_counts, by = strata_pass) %>%
    filter(!is.na(w_survey)) %>%
    ## Normalizes back to 1 after dropping
    ## empty cells
    mutate(w_pop = w_pop * 1/sum(w_pop),
           w = w_pop / w_survey) %>%
    dplyr::select(!!as.symbol(strata_pass), w)

  survey <- survey %>%
    left_join(post_strat)

  return(survey)
}

check_sample <- function(sample, selection_model) {
    check = model.matrix(selection_model, data = sample)
    check = colSums(check)
    fail = check[which(check ==0)]
    fail_bin = length(check[which(check ==0)]) > 0 
    check_prop = check/nrow(sample)
    
    return(list(samp_prop = check_prop,
                fail  = fail, 
                fail_bin = fail_bin,
                counts = check))
}

check_sample_outcome <- function(sample, pop = cces, selection_model, interaction_cols, interaction_cols_2 = NULL) {
    
    vars = all.vars(selection_model)
    var = NULL
    counts = NULL
    prop = NULL
    u_outcome = NULL

    #uninteracted variables
    run_counts <- function(data) {
        for(i in 1:length(vars)) {
            t = data %>% group_by_at(vars[i]) %>%
                summarise(n = n(), 
                          avg_outcome = mean(outcome)) %>%
                mutate(prop = round(n/nrow(sample),4))
            var = c(var, as.character(t[,1] %>% pull()))
            counts = c(counts, t$n)
            prop = c(prop,  t$prop)
            u_outcome = c(u_outcome, t$avg_outcome)
        }
        #interactions
        t = suppressMessages(data %>% group_by_at(interaction_cols) %>% 
                                 summarise(n = n(),
                                           avg_outcome = mean(outcome)) %>%
                                 mutate(prop = round(n/nrow(sample), 4)))
        interaction = apply(t, 1,  function(r) paste(r[1],r[2], collapse = "_"))
        counts = data.frame(var  = c(var, interaction),
                            n = c(counts, t$n), 
                            prop = c(prop, t$prop),
                            avg_outcome = c(u_outcome, t$avg_outcome))
        
        if(!is.null(interaction_cols_2)) {
            t2 = suppressMessages(data %>% group_by_at(interaction_cols_2) %>% 
                                      summarise(n = n(),
                                                avg_outcome = mean(outcome)) %>%
                                      mutate(prop = round(n/nrow(sample), 4)))
            interaction = apply(t2, 1,  function(r) paste(r[1],r[2], collapse = "_"))
            append = cbind(data.frame(var = interaction), t2[, - c(1,2)])
            counts = rbind(counts, append)
        }
    }
    samp_counts = run_counts(sample)
    pop_counts = run_counts(pop)
    
    if(nrow(pop_counts) > nrow(samp_counts)) {
        samp_counts = rbind(samp_counts, data.frame(var = pop_counts$var[!(pop_counts$var %in% samp_counts$var)],
                                                  n = rep(0,sum(!(pop_counts$var %in% samp_counts$var))),
                                                  prop = rep(0,sum(!(pop_counts$var %in% samp_counts$var))),
                                                  avg_outcome = rep(-99,sum(!(pop_counts$var %in% samp_counts$var)))))
    }
    
    
    fail = sum(samp_counts$n == 0)
    bad = sum(samp_counts$prop <= 0.05)
    bad_strata = data.frame(strata = as.character(samp_counts$var[samp_counts$prop <= 0.05]), prop = samp_counts$prop[samp_counts$prop <= 0.05])
    v_bad = sum(samp_counts$prop <= 0.01)
    v_bad_strata = data.frame(strata = as.character(samp_counts$var[samp_counts$prop <= 0.01]), prop = samp_counts$prop[samp_counts$prop <= 0.01])
    counts$var[samp_counts$prop <= 0.01]
    
    samp_counts = samp_counts %>% mutate(leq_5pp = as.numeric(prop <= 0.05),
                               leq_1pp = as.numeric(prop <= 0.01))
    
   
    return(list(samp_counts = samp_counts,
                fail = fail, 
                bad = bad, 
                v_bad = v_bad, 
                bad_strata = bad_strata,
                v_bad_strata = v_bad_strata))
}

check_outcome <- function(outcome) {
    beyond_support = (min(outcome) <0 | max(outcome) > 1)
    return(beyond_support)
}

bound_outcome <- function(outcome, coefs, increment = 1, increment_intercept = .01,  denom = 10, noise, cces_expanded, silent = T) {
   
    fail = check_outcome(outcome)
    while(fail) {
        denom = denom + increment
        if(!silent) { cat(denom, ": ") }
        coefs_use = coefs/denom
        outcome = cces_expanded %*% coefs_use
        
        outcome = outcome + rnorm(nrow(cces_expanded), mean = 0, sd = sd(outcome)*noise)
        summary(outcome)
        if(max(outcome) <=1 & min(outcome) <0 ) {
            coefs[1] = coefs[1] + increment_intercept
            if(!silent) { cat("\nmoving intercept up", coefs[1],  "\n") }
            denom = denom - increment
        }
        if(max(outcome) >1 & min(outcome) >=0 ) {
            coefs[1] = coefs[1] - increment_intercept
            if(!silent) { cat("\nmoving intercept down", coefs[1],  "\n") }
            denom = denom - increment
        }
        fail = check_outcome(outcome)
        if(!silent) { cat(round(min(outcome),2), round(max(outcome),2),  "\n") }
    }
    if(!silent) { cat(paste("Min denom:", denom)) }
    return(list(outcome = outcome, coefs = coefs_use, denom = denom))
    
}


############# Load Data #####################
#these data have been cleaned already see app_modeled for how it was done
### Load Target Data
cces <- readRDS(paste0(path_data, "cces_lasso_061021.rds"))
cces$recode_age_bucket = as.character(cces$recode_age_bucket)
cces$recode_age_3way= as.character(cces$recode_age_3way)

cces <- cces %>% mutate(recode_agesq = recode_age^2/ mean(recode_age^2),
                        recode_agecubed = recode_age^3/ mean(recode_age^3))

 ##################### LASSO: Selection #############################

selection_model = as.formula(~recode_pid_3way + recode_female +
                                 recode_age_bucket 
                                     + recode_educ
                                     + recode_race
                                     + recode_born 
                                     + recode_born:recode_age_bucket
                                     + recode_pid_3way:recode_age_bucket
        )
inter = c("recode_pid_3way", "recode_age_bucket")
inter_2 = c("recode_born", "recode_age_bucket")
cces_expanded = model.matrix(selection_model, data = cces)
coefs = matrix(NA,nrow = ncol(cces_expanded), ncol =1 )
rownames(coefs) = colnames(cces_expanded)
coefs
#(p(S)) for negative bias select non dem voters
#not enough raking bias
up_nonlin = 1
coefs[,1] = c(1, #intercept -5 w race
              2, #selection of indep pos
              2, #selection of R pos
              .5, #male
              .15, #36-50,
              .2, #51-64,
              .2, #65+,
              # .7, #college
              # -1 , #post-grad
              
              .8, #hs grad
            -.6, #some college
            -.2, #2 year collge
             -1, #4 year college
            -1.5, #post grad
              .5,#hispanic
              .3,#other
              .7,#white
              2, #bornagain
              1*up_nonlin,#bornagain x 36-50
              1.5*up_nonlin, #bornagain x 51-64
              2*up_nonlin, #bornagain x 65+
              .3*up_nonlin,#ind x 36-50
              .5*up_nonlin, #rep x 36-50,
              1*up_nonlin, #ind x 51-64,
              1*up_nonlin, #rep x 51-64,
              -.2*up_nonlin, #ind x 65+
              2*up_nonlin #rep x 65+
)
#overly extreme
# coefs[,1] = c(-4,
#               .5, #selection of indep pos
#               1.5, #selection of R pos
#               .5, #male
#               2, #36-50,
#               1.6, #51-64,
#               1, #65+,
#               .8, #hs grad
#               -.6, #some college
#               -.2, #2 year collge
#               # -1, #4 year college
#               # -1.5, #post grad
#               .5,#hispanic
#               .3,#other
#               .7,#white
#               2, #bornagain
#               1*up_nonlin,#bornagain x 36-50
#               1.5*up_nonlin, #bornagain x 51-64
#               2*up_nonlin, #bornagain x 65+
#               -2.4*up_nonlin,#ind x 36-50
#               -2.8*up_nonlin, #rep x 36-50,
#               .7*up_nonlin, #ind x 51-64,
#               -2*up_nonlin, #rep x 51-64,
#               -1.5*up_nonlin, #ind x 65+
#               1.1*up_nonlin #rep x 65+
# )
up_nonlin_pos = 1
up_nonlin_neg = 1
coefs[,1] = c(-0.5, #intercept -5 w race #-.5
              1.5, #selection of indep pos
              2, #selection of R pos
              .5, #male
              .35, #36-50,
              .5, #51-64,
              1, #65+,
              # .7, #college
              # -1 , #post-grad
              .8, #hs grad
              -.6, #some college
              -.2, #2 year collge
              -1, #4 year college
              -1.5, #post grad
              .5,#hispanic
              .3,#other
              .7,#white
              1, #bornagain
              -1.7*up_nonlin_neg,#bornagain x 36-50
              -1.5*up_nonlin_pos, #bornagain x 51-64
              -2.3*up_nonlin_pos, #bornagain x 65+
              -.9*up_nonlin_neg,#ind x 36-50
              -.9*up_nonlin_pos, #rep x 36-50,
              -1.5*up_nonlin_neg, #ind x 51-64,
              -1*up_nonlin_pos, #rep x 51-64,
              -2*up_nonlin_neg, #ind x 65+
              -5.1*up_nonlin_pos #rep x 65+ #pos 2.1
)
xbeta = cces_expanded %*% coefs
p_include = plogis(xbeta)
sum(p_include)
#pS_denom = 65
pS_denom = 45 #- .5 -> 4.1 
#pS_denom = 50
p_include = p_include/pS_denom
sum(p_include)
#plot(density(p_include))
ggplot(data = data.frame(pS = p_include)) +
    geom_density(aes(x = pS)) + 
    theme_bw() #+ 
    #xlim(.01, .017)
ggplot(data = data.frame(pS = p_include)) +
    geom_density(aes(x = pS)) + 
    theme_bw() + 
    xlim(.01, .017)
summary(p_include)
min(p_include)
unique(cces[which(p_include == max(p_include)), all.vars(selection_model) ])


# gg_dat = data.frame(Selection_Probability = p_include,
#                    Pid = cces$recode_pid_3way,
#                    age = cces$recode_age_bucket,
#                    born = cces$recode_born,
#                    Outcome_pD = NA)
# ggplot(gg_dat) +
#     geom_density(aes(x= Selection_Probability, color = age, linetype = Pid)) +
#     ggtitle("Distribution of Selection Probabilities by Party") +
#     theme_bw()



# LA:tinkering

# formula_ps <- selection_model
# cces <- bind_cols(cces, cces %>%
#                       unite("strata", all.vars(formula_ps), remove = FALSE) %>%
#                       unite("strata_reduc", all.vars(formula_ps_reduc),
#                             remove = FALSE) %>%
#                       unite("strata_all", all.vars(formula_ps_all)) %>%
#                       dplyr::select(strata, strata_reduc, strata_all))
# sample <- rbinom(nrow(cces), 1, p_include)
# sum(sample)
# survey_sim <- cces[sample == 1, ]
# 
# check_sample_outcome(survey_sim, selection_model, inter)
# #check_sample_outcome(survey_sim, selection_model, inter_2)
# n_distinct(cces$strata) - n_distinct(survey_sim$strata)
# dropped_strata = unique(cces$strata)[which(!(unique(cces$strata) %in%
#                                                  unique(survey_sim$strata)))]
# cces %>% filter(strata %in% dropped_strata) %>% count()
# 
# la = data.frame(pS = p_include, pid = cces$recode_pid_3way) %>%
#     group_by(pid) %>% summarise(max = round(max(pS)*100,2 ))
# # la
# #
# # check_sample_outcome(survey_sim,selection_model, inter)
# gg_dat = data.frame(Selection_Probability = p_include,
#                     Pid = cces$recode_pid_3way,
#                     age = cces$recode_age_bucket,
#                     born = cces$recode_born,
#                     Outcome_pD = NA)
# gg_p_include_pid = ggplot(gg_dat) +
#     geom_density(aes(x= Selection_Probability, color = Pid)) +
#     annotate(geom = "label", x=quantile(p_include,.25),y=Inf, vjust = 1,
#              color = "red",
#              label =  paste0("Dem Max P(S)= ", la[la$pid =="Dem", "max"], "%")) +
#     annotate(geom = "label",x=quantile(p_include,.25), y=Inf, vjust =3,
#              label= paste0("Ind Max P(S)= ", la[la$pid =="Ind", "max"], "%" ),
#              color = "green") +
#     annotate(geom = "label",x=quantile(p_include,.25), y=Inf, vjust = 5,
#              label= paste0("Rep Max P(S)= ", la[la$pid =="Rep", "max"], "%" ),
#              color = "blue") +
#     # geom_text(x=.10, y=250, label= paste0("Dem Max P(S)=", la[la$pid =="Dem", "max"] ),
#     #           color = "red") +
#     #  geom_text(x=.10, y=200, label= paste0("Ind Max P(S)=", la[la$pid =="Ind", "max"] ),
#     #           color = "green") +
#     #  geom_text(x=.10, y=150, label= paste0("Rep Max P(S)=", la[la$pid =="Rep", "max"] ),
#     #           color = "Blue") +
#     ggtitle("Distribution of Selection Probabilities by Party") +
#     theme_bw()
# gg_p_include_pid
# 
# ggplot(gg_dat) +
#     geom_density(aes(x= Selection_Probability, color = age, linetype = Pid)) +
#     ggtitle("Distribution of Selection Probabilities by Party") +
#     theme_bw()

# 
#most likely to select: rep, 65+, male, protestants
# cces[which(p_include == max(p_include)),
#      c("recode_pid_3way", "recode_age_bucket", "recode_female", "recode_relig_6way")]
# #summary(p_include)

#################### DESIGN OUTCOME MODEL ##################
#p(D)
#this -coefs method gets good bias on ps and other raking but rake all is p good :/
coefs_outcome = -coefs
coefs_outcome = coefs
coefs_outcome[1] = 9
cor(p_include, cces_expanded %*% coefs_outcome)
if(!bern) {
    cat(paste("Adding sd(outcome)*",round(noise, 3), "\n")) 
    set.seed(1383904)
    bound = bound_outcome(outcome = cces_expanded %*% coefs_outcome,
                          coefs = coefs_outcome,
                          denom = 14.5,
                          cces_expanded = cces_expanded,
                          noise = noise, silent = F)
    coefs_outcome = bound$coefs
    xbeta_outcome = bound$outcome
    beyond_support = check_outcome(xbeta_outcome)

    if(min(xbeta_outcome) <0 | max(xbeta_outcome) > 1) {
        warning("Outcome beyond prob support for some units when noise is added",
                                                                immediate. = T)
    }
} else {
    cat(paste("Adding bernoulli noise",noise, "\n"))  
    set.seed(1383904)
    xbeta_outcome = rbinom(nrow(cces), 1, xbeta_outcome) 
    cat(paste("Corr of S (one sample draw) and Y", round(cor(sample, xbeta_outcome),3)))
}
summary(xbeta_outcome)
# cor(xbeta_outcome, p_include)
# cces[which(xbeta_outcome == min(xbeta_outcome)),
#      c("recode_pid_3way", "recode_age_bucket", "recode_female")]
# cces[which(xbeta_outcome == max(xbeta_outcome)),
#      c("recode_pid_3way", "recode_age_bucket", "recode_female")]
# 
#plot(density(xbeta_outcome))
cat(paste("Mean outcome w/noise is", round(mean(xbeta_outcome)*100,3), "\n"))

cat(paste("Range of outcome w/noise is\n"))
cat(paste(summary(xbeta_outcome), "\n"))
s = summary(lm(update(selection_model, xbeta_outcome ~ .),data = cces))
R2_outcome = s$adj.r.squared
cat(paste("R^2 outcome is", round(s$adj.r.squared,3), "\n"))
cat(paste("Mean scaled outcome (target) is", round(mean(xbeta_outcome)*100,3)))
cat(paste("\nCorr of sampling prob and outcome ", round(cor(xbeta_outcome, p_include),3)))
cces$outcome = xbeta_outcome

test = rbinom(nrow(cces), 1, p_include)
(mean(cces[test==1,]$outcome) - mean(xbeta_outcome))*100

# #LA: 
# gg_dat = data.frame(Selection_Probability = p_include,
#                     Pid = cces$recode_pid_3way,
#                     age = cces$recode_age_bucket,
#                     Outcome_pD = cces$outcome)
# gg_p_include_outcome = ggplot(gg_dat) +
#     geom_point(aes(x= Selection_Probability, y= Outcome_pD, color = age, shape = Pid)) +
#     theme_bw()
# gg_p_include_outcome



######### Make STRATA variable in CCES ############
formula_ps <- selection_model
cces <- bind_cols(cces, cces %>%
                      unite("strata", all.vars(formula_ps), remove = FALSE) %>%
                      dplyr::select(strata))


#################### Targets ###################
if(POPW) {
  cces_svy <- svydesign(ids = ~1, weights = ~commonweight_vv_post, data = cces)
} else {
  cces_svy <- suppressWarnings(svydesign(ids = ~1, data = cces))
}
margin_sim = svymean(~outcome, cces_svy)[1]* 100
margin_sim
targets_rake_demos_noeduc <- create_targets(cces_svy,
                                            formula_rake_demos_noeduc)
targets_rake_demos_weduc <- create_targets(cces_svy, formula_rake_demos_weduc)
targets_rake_all_vars <- create_targets(cces_svy,
                                        formula_rake_all_vars)
targets_demo_truth <- create_targets(cces_svy, selection_model)


## Make table of Population Counts for post-stratification for manual ps function
cces_counts <- cces %>%
  group_by(strata) %>%
  summarize(n = if(!POPW) {n()} else {sum(commonweight_vv_post, na.rm = TRUE)}) %>%
  ungroup() %>%
  mutate(w = n / sum(n, na.rm = TRUE))



########################### PREP SIMS ##################################

margins_formula <- ~recode_vote_2016 + 
  mod_cces_on_cces_pR + mod_cces_on_cces_pD + mod_cces_on_cces_pO+
  diff_cces_on_cces + margin_cces_on_cces + 
  recode_pid_3way + 
  recode_female + recode_race +recode_region + recode_educ + recode_relig_6way + 
  recode_born + recode_attndch_4way + recode_income_5way +
  recode_age_bucket + recode_age + recode_agesq + recode_agecubed + recode_age_factor + 
  recode_race_educ_reg + recode_educ_wh_3way + 
  recode_educ_pid_race +
  recode_pid_race + 
  recode_educ_pid +
  recode_race_reg_wh_educ + 
  recode_midwest_edu_race +
  recode_midwest_wh_edu

est_mean <- function(outcome, design) {
  svymean(as.formula(paste0("~", outcome)), design, na.rm = TRUE)[1]
}

#########################################
############## variance calc ###########
## Variance functions
var_fixed <- function(Y, weights, pop_size) {
    ## note: needs weights that sum to population total
    if(round(sum(weights)) != pop_size) { weights = weights*pop_size/sum(weights)}
    return(Hmisc::wtd.var(Y, weights))
}

## kott (14) (under poisson)
var_quasi <- function(weights, residuals, pop_size) {
    #moving from kott 14 w sum w =N to weights that sum to 1 + var of total to var of mean:
    #sum^n (w_i^2 - 1/N_pop w_i)e_i^2 
    return(sum((weights^2 - (weights / pop_size))*residuals^2))
}

## kott (15) linearization
var_linear <- function(weights, residuals, sample_size) {
    #moving from kott 14 w sum w =N to weights that sum to 1 + var of total to var of mean:
    # n/(n-1) sum^n (w_i*e_i)^2 - (1/n-1) [sum^n] *using this for now
    # approx = sum^n (w_i*e_i)^2 - (1/n) [sum^n]
    n = sample_size
    return((n/(n-1))*sum((weights * residuals)^2) - 1/n * sum(weights * residuals)^2)
}

## chad
var_chad <- function(weights, residuals) {
    return(sum(weights^2 * residuals^2))
}

## calculate all variances
calc_SEs <- function(Y, residuals, pop_size, weights, sample_size) {
    if(round(sum(weights)) != 1 ) {
        weights = weights/sum(weights)
    }
    return(data.frame(SE_fixed = sqrt(var_fixed(Y, weights, pop_size) / length(Y)),
                      SE_quasi = sqrt(var_quasi(weights, residuals, pop_size)),
                      SE_linear = sqrt(var_linear(weights, residuals, sample_size)),
                      SE_chad = sqrt(var_chad(weights, residuals))))
}

########################### RUN SIMS ##################################
if(detectCores() < 20) {
    nsims = 100
    eval_kpop = F
    SAVE = F
    kpop_constraints = F
    run_local = T
} else {
    run_local = F
}
#rstudio is dumb:
rstudio_para_cat <- function(...){
    system(sprintf('echo "\n%s\n"', paste0(..., collapse="")))
}


#double check params are as expected:
nsims
sum(p_include)
SAVE
eval_kpop
run_local
up_nonlin
manual_lambda
pryr::mem_used()
rm(bound, cces_expanded, s, xbeta, xbeta_outcome)
pryr::mem_used()
# lambda_min
system.time({
    sims <- mclapply(1:nsims, function(nsim) {
        #
        system.time({
            cat(paste("=====================  SIM:",nsim, 
                      "===================== \n"))
            if(run_local) {
                rstudio_para_cat(paste("=====================  SIM:",nsim, 
                                       "====================="))
            }
            
            sample <- rbinom(nrow(cces), 1, p_include)
            survey_sim <- cces[sample == 1, ]
            
            survey_design <- suppressWarnings(svydesign(ids = ~1, data = survey_sim))
            
            ########################### check sample ##########################
            check_s = check_sample(survey_sim, selection_model)
            bad_sample = check_s$fail_bin
            if(bad_sample ==1 & run_local) {
                rstudio_para_cat(c("SAMPLING FAILURE: ", names(check_s$fail), " strata has ", as.character(check_s$fail), " sampled units\n"))
            }
            check_2 = check_sample_outcome(survey_sim, pop = cces, selection_model, interaction_cols = inter, interaction_cols_2 = inter_2)
            check_nums = c(leq_5pp = check_2$bad,
                           leq_1pp = check_2$v_bad, 
                           fail = check_2$fail)
            s = survey_sim %>% group_by(recode_pid_3way,recode_female, recode_age_bucket,
                                        recode_educ_3way) %>% count() %>%
                mutate(n_s = round(n/nrow(survey_sim), 3))
            c = cces %>% group_by(recode_pid_3way, recode_female, recode_age_bucket,
                                  recode_educ_3way) %>% count() %>%
                mutate(n_c = round(n/nrow(cces), 3))
            count = nrow(c) - nrow(s)
            
            ############################################
            ## Unweighted estimate
            ############################################
            
            unweighted <- est_mean("outcome", survey_design)
           
            ############################################
            ## Sample size
            ############################################\
            n <- sum(sample) 
            
            ############################################
            ## Raking on demographics (no education)
            ############################################
            rake_demos_noeduc_svyd <- try(calibrate(design = survey_design,
                                                    formula = formula_rake_demos_noeduc,
                                                    population = targets_rake_demos_noeduc,
                                                    calfun = "raking"), silent = T)
            
            rake_demos_noeduc <- tryCatch(est_mean("outcome", rake_demos_noeduc_svyd), error = function(e) NA)
            
            #SEs
            if(class(rake_demos_noeduc_svyd)[1] == "try-error") {
                rake_demos_noeduc_se <- data.frame(SE_fixed = NA,
                                                   SE_quasi = NA, 
                                                   SE_linear = NA, 
                                                   SE_chad = NA) 
                names(rake_demos_noeduc_se) = paste0("rake_demos_noeduc_", names(rake_demos_noeduc_se))
                rake_demos_noeduc_se_SVY = data.frame(rake_demos_noeduc_se_SVY  = NA)
            } else {
                residuals = residuals(lm(update(formula_rake_demos_noeduc, outcome ~ .), 
                                         data = rake_demos_noeduc_svyd$variables))
                
                res_rake_demos_noeduc = data.frame(min = min(residuals), 
                                                   perc_25 = quantile(residuals, .25), 
                                                   mean = mean(residuals),
                                                   perc_75 = quantile(residuals, .75),
                                                   var = var(residuals))
                
                rake_demos_noeduc_se <- calc_SEs(Y = rake_demos_noeduc_svyd$variables$outcome, 
                                                 residuals = residuals, 
                                                 pop_size = nrow(cces), 
                                                 sample_size = sum(sample),
                                                 weights = weights(rake_demos_noeduc_svyd))
                names(rake_demos_noeduc_se) = paste0("rake_demos_noeduc_", names(rake_demos_noeduc_se))
                
                
                rake_demos_noeduc_se_SVY = data.frame(rake_demos_noeduc_se_SVY = data.frame(svymean(~outcome, 
                                                                                                    rake_demos_noeduc_svyd, 
                                                                                                    na.rm = TRUE))[1,2])
                
            }
            ############################################
            #### Raking on demographics (with education)
            ############################################
            
            rake_demos_weduc_svyd <- try(calibrate(design = survey_design,
                                                   formula = formula_rake_demos_weduc,
                                                   population = targets_rake_demos_weduc,
                                                   calfun = "raking"), silent = T)
            
            rake_demos_weduc <- tryCatch(est_mean("outcome", rake_demos_weduc_svyd), 
                                         error = function(e) NA)
            
            
            #SEs
            if(class(rake_demos_weduc_svyd)[1] == "try-error") {
                rake_demos_weduc_se <- data.frame(SE_fixed = NA,
                                                  SE_quasi = NA, 
                                                  SE_linear = NA, 
                                                  SE_chad = NA) 
                names(rake_demos_weduc_se) = paste0("rake_demos_weduc_", names(rake_demos_weduc_se))
                
            } else {
                residuals = residuals(lm(update(formula_rake_demos_weduc, outcome ~ .), 
                                         data = rake_demos_weduc_svyd$variables))
                res_rake_demos_wedu = data.frame(min = min(residuals), 
                                                 perc_25 = quantile(residuals, .25), 
                                                 mean = mean(residuals),
                                                 perc_75 = quantile(residuals, .75),
                                                 var = var(residuals))
                rake_demos_weduc_se <- calc_SEs(Y = rake_demos_weduc_svyd$variables$outcome, 
                                                residuals = residuals, 
                                                pop_size = nrow(cces),
                                                sample_size = sum(sample),
                                                weights = weights(rake_demos_weduc_svyd))
                names(rake_demos_weduc_se) = paste0("rake_demos_weduc_", names(rake_demos_weduc_se))
                
            }
            
            ############################################
            #### Raking on everything
            ############################################
            rake_all_svyd <- try(calibrate(design = survey_design,
                                           formula = formula_rake_all_vars,
                                           population = targets_rake_all_vars,
                                           calfun = "raking"))
            
            rake_all <- tryCatch(est_mean("outcome", rake_all_svyd), 
                                 error = function(e) NA)
            #SEs
            if(class(rake_all_svyd)[1] == "try-error") {
                rake_all_se <- data.frame(SE_fixed = NA,
                                          SE_quasi = NA, 
                                          SE_linear = NA, 
                                          SE_chad = NA) 
                names(rake_all_se) = paste0("rake_all_", names(rake_all_se))
                
            } else {
                residuals = residuals(lm(update(formula_rake_all_vars, outcome ~ .), 
                                         data = rake_all_svyd$variables))
                res_rake_all = data.frame(min = min(residuals), 
                                          perc_25 = quantile(residuals, .25), 
                                          mean = mean(residuals),
                                          perc_75 = quantile(residuals, .75),
                                          var = var(residuals))
                rake_all_se <- calc_SEs(Y = rake_all_svyd$variables$outcome, 
                                        residuals = residuals, 
                                        pop_size = nrow(cces), 
                                        sample_size = sum(sample),
                                        weights = weights(rake_all_svyd))
                names(rake_all_se) = paste0("rake_all_", names(rake_all_se))
                
            }
            
            ############################################
            ## Post-stratification: Truth
            ############################################
            
            #track empty cells:
            #this subsets cces strata to only those in survey_sim
            missing_strata <- unique(cces$strata)[!(unique(cces$strata) %in%
                                                        unique(survey_sim$strata))]
            cat(round(length(missing_strata)/ length(unique(cces$strata)),3),
                "% cces original strata missing from sample, ",
                " and", cces %>% filter(strata %in% missing_strata) %>% summarise(n()) %>% pull(), "/", nrow(cces), "units\n" )
            dropped_cells = cces %>% filter(strata %in% missing_strata) %>% group_by(strata) %>% count()
            #dropped_cells = data.frame(sum = sum(dropped_cells$n), strata = paste(dropped_cells$strata, collapse = " | "))
            dropped_cells = sum(dropped_cells$n)
            
            #dropped_cells_all = data.frame(sum = sum(dropped_cells_all$n), strata = paste(dropped_cells_all$strata, collapse = " | "))
            
            #note that we no longer have the issue of pew having strata that cces doesn bc
            #we use survey_sim a sample of cces so it will never have diff strata
            
            post_stratification_svyd = svydesign(~1, data = postStrat(survey_sim, 
                                                                      cces_counts, "w", 
                                                                      strata_pass = "strata", 
                                                                      warn = F),
                                                 weights = ~w)
            
            post_stratification <- est_mean("outcome", post_stratification_svyd)
            
            #SEs
            post_stratification_se <- data.frame(post_strat_SE_svy = data.frame(svymean(~outcome, post_stratification_svyd, na.rm = TRUE))[1,2])
            
            
            ############################################
            #### Raking on true model
            ############################################
            #very messy error catching for the moment just a stop gap to see how things look
            rake_truth_svyd <- try(calibrate(design = survey_design,
                                             formula = selection_model,
                                             population = targets_demo_truth,
                                             maxit = 100,
                                             calfun = "raking"), silent = T)
            
            if(class(rake_truth_svyd)[1] == "try-error") {
                rake_truth_svyd <- try(calibrate(design = survey_design,
                                                 formula = selection_model,
                                                 population = targets_demo_truth,
                                                 calfun = "raking",
                                                 maxit = 100,
                                                 epsilon = .009), silent = T)
            }
            
            if(class(rake_truth_svyd)[1] == "try-error") {
                
                rake_truth_se <- data.frame(SE_fixed = NA,
                                            SE_quasi = NA, 
                                            SE_linear = NA, 
                                            SE_chad = NA) 
                names(rake_truth_se) = paste0("rake_truth_", names(rake_truth_se))
                res_rake_truth = data.frame(min = NA, 
                                            perc_25 = NA, 
                                            mean = NA,
                                            perc_75 = NA,
                                            var = NA)
                
            } else {
                lambdas <- 10^seq(3, -2, by = -.1)
                x <- model.matrix(update(selection_model, outcome ~ .),
                                  data = rake_truth_svyd$variables)[, -1]
                fit <- glmnet(x, 
                              rake_truth_svyd$variables$outcome, alpha = 0, lambda = lambdas)
                cv_fit <- cv.glmnet(x, rake_truth_svyd$variables$outcome, alpha = 0, lambda = lambdas)
                opt_lambda <- cv_fit$lambda.min
                fit <- cv_fit$glmnet.fit
                
                residuals = rake_truth_svyd$variables$outcome - predict(fit, s = opt_lambda, newx = x)
                res_rake_truth = data.frame(min = min(residuals), 
                                            perc_25 = quantile(residuals, .25), 
                                            mean = mean(residuals),
                                            perc_75 = quantile(residuals, .75),
                                            var = var(residuals))
                rake_truth_se <- tryCatch(calc_SEs(Y = rake_truth_svyd$variables$outcome,
                                                   residuals = residuals,
                                                   pop_size = nrow(cces),
                                                   sample_size = sum(sample),
                                                   weights = weights(rake_truth_svyd)), error = function(e) NA)
                names(rake_truth_se) = paste0("rake_truth_", names(rake_truth_se))
            }
            
            rake_truth <- tryCatch(est_mean("outcome", rake_truth_svyd), 
                                   error = function(e) NA)
            truth_margins <- tryCatch(svymean(margins_formula, rake_truth_svyd), 
                                      error = function(e) NA)
            
            #HT and Hayek
            p_sample <- as.matrix((p_include[sample==1]))
            
            #HT
            ht_truth = sum((cces[sample ==1, "outcome"]/p_sample))/nrow(cces) 
            #Hayek
            hayek_truth = sum((cces[sample ==1, "outcome"]/p_sample))/sum(1/p_sample)
            
            ############################################
            ## Kpop: Categorical Data + b = argmax V(K)
            ############################################
            if(eval_kpop) {
                
                # Select the covariates for use in Kbal: updated cat data no cont age
                #one-hot coded for cat kernel
                kbal_data <- bind_rows(survey_sim %>% dplyr::select(recode_age_bucket,
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
                
                kbal_data_reduc <- bind_rows(survey_sim %>% dplyr::select(all.vars(selection_model)),
                                             cces %>% dplyr::select(all.vars(selection_model)))
                
                kbal_data_sampled <- c(rep(1, nrow(survey_sim)), rep(0, nrow(cces)))
                ##### Demos Constraint
                rake_demos_constraint <- bind_rows(survey_sim %>% dplyr::select(recode_age_bucket,
                                                                                recode_female,
                                                                                recode_race,
                                                                                recode_region,
                                                                                recode_pid_3way),
                                                   cces %>% dplyr::select(recode_age_bucket,
                                                                          recode_female,
                                                                          recode_race,
                                                                          recode_region,
                                                                          recode_pid_3way))%>%
                    model.matrix(as.formula("~."), .)
                
                rake_demos_constraint <- rake_demos_constraint[,-1]
                rake_demos_constraint <- scale(rake_demos_constraint)
                
                
                rake_demos_wedu_constraint <- bind_rows(survey_sim %>% dplyr::select(recode_age_bucket,
                                                                                     recode_female,
                                                                                     recode_race,
                                                                                     recode_region,
                                                                                     recode_pid_3way,
                                                                                     recode_educ),
                                                        cces %>% dplyr::select(recode_age_bucket,
                                                                               recode_female,
                                                                               recode_race,
                                                                               recode_region,
                                                                               recode_pid_3way,
                                                                               recode_educ))%>%
                    model.matrix(as.formula("~."), .)
                
                rake_demos_wedu_constraint <- rake_demos_wedu_constraint[,-1]
                rake_demos_wedu_constraint <- scale(rake_demos_wedu_constraint)
                
                
                rake_all_constraint <- bind_rows(survey_sim %>% dplyr::select(recode_age_bucket,
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
                                                                        recode_attndch_4way)) %>%
                    model.matrix(as.formula("~."), .)
                
                rake_all_constraint <- rake_all_constraint[,-1]
                rake_all_constraint <- scale(rake_all_constraint)
                
                
                #cat("===================== Running Kbal ==================\n")
                
                #### DEFAULT ######
                #rstudio_para_cat(c("nsim: ", nsim, " DEFAULT"))
                cat(paste("nsim:", nsim, "DEFAULT", "\n"))
                kbal_est <- kbal(allx=kbal_data,
                                 sampled = kbal_data_sampled,
                                 #b = b_manual[i],
                                 cat_data = TRUE,
                                 incrementby = increment,
                                 meanfirst = FALSE,
                                 ebal.tol = tolerance,
                                 ebal.maxit = maxit,
                                 minnumdims = min_num_dims,
                                 maxnumdims = max_num_dims,
                                 linkernel =  FALSE,
                                 sampledinpop = FALSE,
                                 fullSVD = TRUE)
                
                kpop_svyd <- svydesign(~1, data = survey_sim,
                                       weights = kbal_est$w[kbal_data_sampled ==1])
                
                kpop <- est_mean("outcome", kpop_svyd)
                b_kpop = kbal_est$b
                #save memory by saving only the svd to re use
                svdK = kbal_est$svdK 
                numdims = kbal_est$numdims
                biasbound_r = kbal_est$biasbound_ratio
                biasbound = kbal_est$biasbound_opt
                
                ##### Kpop SEs
                kpop <- tryCatch(est_mean("outcome", kpop_svyd), error = function(e) NA)
                
                lambdas <- if(manual_lambda) { 10^seq(3, -2, by = -.1) } else {NULL}
                
                x <- as.matrix(data.frame(kbal_dims = kbal_est$svdK$v[, 1:kbal_est$numdims]))
                cv_fit <- cv.glmnet(x, kpop_svyd$variables$outcome, alpha = 0, lambda = lambdas)
                lambda_pass = if(lambda_min) { cv_fit$lambda.min} else {cv_fit$lambda.1se}
                
                residuals = kpop_svyd$variables$outcome - predict(cv_fit$glmnet.fit,
                                                                  s = lambda_pass, newx = x)
                res_kpop = data.frame(min = min(residuals), 
                                      perc_25 = quantile(residuals, .25), 
                                      mean = mean(residuals),
                                      perc_75 = quantile(residuals, .75),
                                      var = var(residuals))
                kpop_se <- tryCatch(calc_SEs(Y = kpop_svyd$variables$outcome,
                                             residuals = residuals,
                                             pop_size = nrow(cces),
                                             sample_size = sum(sample),
                                             weights = weights(kpop_svyd)), error = function(e) NA)
                
                if(length(kpop_se) == 1) {
                    kpop_se <- data.frame(SE_fixed = NA, 
                                          SE_quasi = NA, 
                                          SE_linear = NA, 
                                          SE_chad = NA)
                }
                names(kpop_se) = tryCatch(paste0("kpop_", names(kpop_se)), error = function(e) NA)
                
                #CONVERGED
                dist_record = data.frame(t(kbal_est$dist_record))
                min_converged = dist_record[which.min(dist_record[dist_record$Ebal.Convergence ==1,"BiasBound"]), "Dims"]
                
                rm(kbal_est, residuals, x, cv_fit)
                
                #CONVERGED
                #### CONVG ####
                cat(paste("nsim:", nsim, " CONV", "\n"))
                if(is.null(min_converged) | length(min_converged) ==0) {
                    kpop_svyd_conv <- "dn converge"
                    kpop_conv <- "dn converge"
                    
                    numdims_conv = "dn converge"
                    biasbound_r_conv = "dn converge"
                    biasbound_conv = "dn converge"
                    kpop_conv_se = data.frame(SE_fixed = NA, 
                                              SE_quasi = NA, 
                                              SE_linear = NA, 
                                              SE_chad = NA)
                    res_kpop_conv = data.frame(min = NA, 
                                               perc_25 = NA, 
                                               mean =NA,
                                               perc_75 = NA,
                                               var = NA)
                    
                } else {
                    kbal_est_conv <- kbal(allx=kbal_data,
                                          K.svd = svdK,
                                          sampled = kbal_data_sampled,
                                          numdims = min_converged,
                                          ebal.tol = tolerance,
                                          ebal.maxit = maxit,
                                          minnumdims = min_num_dims,
                                          maxnumdims = max_num_dims,
                                          scale_data = FALSE,
                                          drop_MC = FALSE,
                                          incrementby = increment,
                                          meanfirst = FALSE,
                                          sampledinpop = FALSE,
                                          ebal.convergence = TRUE)
                    kpop_svyd_conv <- svydesign(~1, data = survey_sim,
                                                weights = kbal_est_conv$w[kbal_data_sampled ==1])
                    kpop_conv <- est_mean("outcome", kpop_svyd_conv)
                    
                    numdims_conv = kbal_est_conv$numdims
                    biasbound_r_conv = kbal_est_conv$biasbound_ratio
                    biasbound_conv = kbal_est_conv$biasbound_opt
                    
                    #SEs
                    x <- as.matrix(data.frame(kbal_dims = kbal_est_conv$svdK$v[, 1:kbal_est_conv$numdims]))
                    cv_fit <- cv.glmnet(x, kpop_svyd_conv$variables$outcome, alpha = 0, 
                                        lambda = lambdas)
                    fit <- cv_fit$glmnet.fit
                    lambda_pass = if(lambda_min) { cv_fit$lambda.min} else {cv_fit$lambda.1se}
                    residuals = kpop_svyd_conv$variables$outcome - predict(cv_fit$glmnet.fit, 
                                                                           s = lambda_pass, 
                                                                           newx = x)
                    res_kpop_conv = data.frame(min = min(residuals), 
                                               perc_25 = quantile(residuals, .25), 
                                               mean = mean(residuals),
                                               perc_75 = quantile(residuals, .75),
                                               var = var(residuals))
                    kpop_conv_se <- tryCatch(calc_SEs(Y = kpop_svyd_conv$variables$outcome,
                                                      residuals = residuals,
                                                      pop_size = nrow(cces),
                                                      sample_size = sum(sample),
                                                      weights = weights(kpop_svyd_conv)), error = function(e) NA)
                    if(length(kpop_conv_se) == 1) {
                        kpop_conv_se <- data.frame(SE_fixed = NA, 
                                                   SE_quasi = NA, 
                                                   SE_linear = NA, 
                                                   SE_chad = NA)
                    }
                    names(kpop_conv_se) = tryCatch(paste0("kpop_conv_", names(kpop_conv_se)), error = function(e) NA)
                    #KRLS SEs are exactly the same for coverged
                    rm(kbal_est_conv, residuals, x, cv_fit) 
                }
                
                
                ####### MF #######
                #rstudio_para_cat(c("nsim: ", nsim, " + aMF"))
                cat(paste("nsim:", nsim, "aMEANFIRST", "\n"))
                #### BROKEN HERE
                kbal_mf_est <- kbal(K.svd = svdK,
                                    cat_data = T,
                                    allx=kbal_data,
                                    sampled = kbal_data_sampled,
                                    ebal.tol = tolerance,
                                    ebal.maxit = maxit,
                                    minnumdims = min_num_dims,
                                    maxnumdims = max_num_dims,
                                    incrementby = increment,
                                    meanfirst = TRUE,
                                    sampledinpop = FALSE)
                
                kpop_mf_svyd <- svydesign(~1, data = survey_sim, 
                                          weights = kbal_mf_est$w[kbal_data_sampled ==1])
                
                kpop_mf <- est_mean("outcome", kpop_mf_svyd)
                
                mfnumdims = kbal_mf_est$numdims
                mf_appended_dims = kbal_mf_est$meanfirst_dims
                if(is.null(mf_appended_dims)) {mf_appended_dims = c(NA)}
                biasbound_r_mf = kbal_mf_est$biasbound_ratio
                biasbound_mf = kbal_mf_est$biasbound_opt
                
                if(is.null(mfnumdims)) {
                    mfnumdims = c(NA) 
                    kpop_mf_se = data.frame(SE_fixed = NA, 
                                            SE_quasi = NA, 
                                            SE_linear = NA, 
                                            SE_chad = NA)
                } else {
                    
                    V <-  data.frame(kbal_dims = kbal_mf_est$svdK$v[, c(1:kbal_mf_est$numdims)])
                    #binding mf cols for sample units to V
                    X <- as.matrix(cbind(kbal_mf_est$appended_constraint_cols[kbal_data_sampled==1, ], V))
                    
                    cv_fit <- cv.glmnet(X, kpop_mf_svyd$variables$outcome, alpha = 0,
                                        lambda = lambdas, 
                                        penalty.factor = c(rep(0, kbal_mf_est$meanfirst_dims), 
                                                           rep(1, kbal_mf_est$numdims)))
                    lambda_pass = if(lambda_min) { cv_fit$lambda.min} else {cv_fit$lambda.1se}
                    residuals = kpop_mf_svyd$variables$outcome - predict(cv_fit$glmnet.fit, 
                                                                         s = lambda_pass, 
                                                                         newx = X)
                    res_kpop_mf = data.frame(min = min(residuals), 
                                             perc_25 = quantile(residuals, .25), 
                                             mean = mean(residuals),
                                             perc_75 = quantile(residuals, .75),
                                             var = var(residuals))
                    
                    kpop_mf_se <- tryCatch(calc_SEs(Y = kpop_mf_svyd$variables$outcome,
                                                    residuals = residuals,
                                                    pop_size = nrow(cces),
                                                    sample_size = sum(sample),
                                                    weights = weights(kpop_mf_svyd)), 
                                           error = function(e) NA)
                    if(length(kpop_mf_se) == 1) {
                        kpop_mf_se <- data.frame(SE_fixed = NA, 
                                                 SE_quasi = NA, 
                                                 SE_linear = NA, 
                                                 SE_chad = NA)
                    }
                    names(kpop_mf_se) = tryCatch(paste0("kpop_mf_", names(kpop_mf_se)),
                                                 error = function(e) NA)
                    
                }
                
                rm(kbal_mf_est, residuals,X, V, cv_fit)
                
                if(kpop_constraints) {
                    #########demos constraint method:
                    cat(paste("nsim:", nsim, "CONSTR", "\n"))
                    kbal_demos_est <- kbal(K.svd = svdK,
                                           allx=kbal_data,
                                           cat_data = TRUE,
                                           sampled = kbal_data_sampled,
                                           ebal.tol = tolerance,
                                           ebal.maxit = maxit,
                                           minnumdims = min_num_dims,
                                           maxnumdims = max_num_dims,
                                           scale_data = FALSE,
                                           drop_MC = FALSE,
                                           incrementby = increment,
                                           #scaling these
                                           meanfirst= TRUE,
                                           mf_columns = c(all.vars(formula_rake_demos_noeduc)),
                                           #constraint = rake_demos_constraint,
                                           #meanfirst = FALSE,
                                           sampledinpop = FALSE)
                    kpop_demos_svyd <- svydesign(~1, data = survey_sim, 
                                                 weights = kbal_demos_est$w[kbal_data_sampled ==1])
                    
                    kpop_demos <- est_mean("outcome", kpop_demos_svyd)
                    
                    numdims_demos = kbal_demos_est$numdims
                    if(is.null(numdims_demos)) {
                        numdims_demos = c(NA) 
                        kpop_demos_se <- data.frame(SE_fixed = NA, 
                                                    SE_quasi = NA, 
                                                    SE_linear = NA, 
                                                    SE_chad = NA)
                    } else {
                        V <-  data.frame(kbal_dims = kbal_demos_est$svdK$v[, c(1:kbal_demos_est$numdims)])
                        X <- as.matrix(cbind(kbal_demos_est$appended_constraint_cols[kbal_data_sampled==1, ], V))
                        
                        cv_fit <- cv.glmnet(X, kpop_demos_svyd$variables$outcome, alpha = 0, lambda = lambdas,
                                            penalty.factor = c(rep(0, ncol(kbal_demos_est$appended_constraint_cols)), rep(1, kbal_demos_est$numdims)))
                        
                        lambda_pass = if(lambda_min) { cv_fit$lambda.min} else {cv_fit$lambda.1se}
                        residuals =  kpop_demos_svyd$variables$outcome - predict(cv_fit$glmnet.fit,
                                                                                 s = lambda_pass, 
                                                                                 newx = X)
                        res_kpop_demos = data.frame(min = min(residuals), 
                                                    perc_25 = quantile(residuals, .25), 
                                                    mean = mean(residuals),
                                                    perc_75 = quantile(residuals, .75),
                                                    var = var(residuals))
                        
                        kpop_demos_se <- tryCatch(calc_SEs(Y = kpop_demos_svyd$variables$outcome,
                                                           residuals = residuals,
                                                           pop_size = nrow(cces),
                                                           sample_size = sum(sample),
                                                           weights = weights(kpop_demos_svyd)), 
                                                  error = function(e) NA)
                        if(length(kpop_demos_se) == 1) {
                            kpop_demos_se <- data.frame(SE_fixed = NA, 
                                                        SE_quasi = NA, 
                                                        SE_linear = NA, 
                                                        SE_chad = NA)
                        }
                        names(kpop_demos_se) = tryCatch(paste0("kpop_demos_", names(kpop_demos_se)),
                                                        error = function(e) NA)
                    }
                    biasbound_r_demos = kbal_demos_est$biasbound_ratio
                    biasbound_demos = kbal_demos_est$biasbound_opt
                    
                    rm(kbal_demos_est, residuals, X, V, cv_fit)
                    
                    
                    #########demos + educ constraint method:
                    cat(paste("nsim:", nsim, "CONSTR", "\n"))
                    kbal_demos_wedu_est <- kbal(K.svd = svdK,
                                                allx=kbal_data,
                                                cat_data = TRUE,
                                                sampled = kbal_data_sampled,
                                                ebal.tol = tolerance,
                                                ebal.maxit = maxit,
                                                minnumdims = min_num_dims,
                                                maxnumdims = max_num_dims,
                                                scale_data = FALSE,
                                                drop_MC = FALSE,
                                                incrementby = increment,
                                                meanfirst= TRUE,
                                                mf_columns = all.vars(formula_rake_demos_weduc),
                                                #scaling these
                                                #constraint = rake_demos_wedu_constraint,
                                                #meanfirst = FALSE,
                                                sampledinpop = FALSE)
                    kpop_demos_wedu_svyd <- svydesign(~1, data = survey_sim, 
                                                      weights = kbal_demos_wedu_est$w[kbal_data_sampled ==1])
                    
                    kpop_demos_wedu <- est_mean("outcome", kpop_demos_wedu_svyd)
                    
                    numdims_demos_wedu = kbal_demos_wedu_est$numdims
                    if(is.null(numdims_demos_wedu)) {
                        numdims_demos_wedu = c(NA)
                        kpop_demos_wedu_se <- data.frame(SE_fixed = NA, 
                                                         SE_quasi = NA, 
                                                         SE_linear = NA, 
                                                         SE_chad = NA)
                    } else {
                        V <-  data.frame(kbal_dims = kbal_demos_wedu_est$svdK$v[, c(1:kbal_demos_wedu_est$numdims)])
                        X <- as.matrix(cbind(kbal_demos_wedu_est$appended_constraint_cols[kbal_data_sampled==1, ], V))
                        
                        cv_fit <- cv.glmnet(X, kpop_demos_wedu_svyd$variables$outcome, alpha = 0, lambda = lambdas,
                                            penalty.factor = c(rep(0, ncol(kbal_demos_wedu_est$appended_constraint_cols)),
                                                               rep(1, kbal_demos_wedu_est$numdims)))
                        
                        lambda_pass = if(lambda_min) { cv_fit$lambda.min} else {cv_fit$lambda.1se}
                        residuals =  kpop_demos_wedu_svyd$variables$outcome - predict(cv_fit$glmnet.fit,
                                                                                      s = lambda_pass, 
                                                                                      newx = X)
                        res_kpop_demos_wedu = data.frame(min = min(residuals), 
                                                         perc_25 = quantile(residuals, .25), 
                                                         mean = mean(residuals),
                                                         perc_75 = quantile(residuals, .75),
                                                         var = var(residuals))
                        kpop_demos_wedu_se <- tryCatch(calc_SEs(Y = kpop_demos_wedu_svyd$variables$outcome,
                                                                residuals = residuals,
                                                                pop_size = nrow(cces),
                                                                sample_size = sum(sample),
                                                                weights = weights(kpop_demos_wedu_svyd)), 
                                                       error = function(e) NA)
                        if(length(kpop_demos_wedu_se) == 1) {
                            kpop_demos_wedu_se <- data.frame(SE_fixed = NA, 
                                                             SE_quasi = NA, 
                                                             SE_linear = NA, 
                                                             SE_chad = NA)
                        }
                        names(kpop_demos_wedu_se) = tryCatch(paste0("kpop_demos_wedu_", names(kpop_demos_wedu_se)),
                                                             error = function(e) NA)
                        
                    }
                    biasbound_r_demos_wedu = kbal_demos_wedu_est$biasbound_ratio
                    biasbound_demos_wedu = kbal_demos_wedu_est$biasbound_opt
                    
                    rm(kbal_demos_wedu_est, residuals, X, V, cv_fit)
                    
                    
                    #########all constraint method:
                    cat(paste("nsim:", nsim, "CONSTR", "\n"))
                    kbal_all_est <- kbal(K.svd = svdK,
                                         allx=kbal_data,
                                         #cat_data = TRUE,
                                         sampled = kbal_data_sampled,
                                         ebal.tol = tolerance,
                                         ebal.maxit = maxit,
                                         minnumdims = min_num_dims,
                                         maxnumdims = max_num_dims,
                                         scale_data = FALSE,
                                         drop_MC = FALSE,
                                         incrementby = increment,
                                         #scaling these
                                         meanfirst= TRUE,
                                         mf_columns = all.vars(formula_rake_all_vars),
                                         #constraint = rake_all_constraint,
                                         #meanfirst = FALSE,
                                         sampledinpop = FALSE)
                    kpop_all_svyd <- svydesign(~1, data = survey_sim, 
                                               weights = kbal_all_est$w[kbal_data_sampled ==1])
                    
                    kpop_all <- est_mean("outcome", kpop_all_svyd)
                    
                    numdims_all = kbal_all_est$numdims
                    if(is.null(numdims_all)) {
                        numdims_all = c(NA)
                        numdims_all_se <- data.frame(SE_fixed = NA, 
                                                     SE_quasi = NA, 
                                                     SE_linear = NA, 
                                                     SE_chad = NA)
                    } else {
                        V <-  data.frame(kbal_dims = kbal_all_est$svdK$v[, c(1:kbal_all_est$numdims)])
                        X <- as.matrix(cbind(kbal_all_est$appended_constraint_cols[kbal_data_sampled==1, ], V))
                        
                        cv_fit <- cv.glmnet(X, kpop_all_svyd$variables$outcome, alpha = 0, lambda = lambdas,
                                            penalty.factor = c(rep(0, ncol(kbal_all_est$appended_constraint_cols)),
                                                               rep(1, kbal_all_est$numdims)))
                        
                        lambda_pass = if(lambda_min) { cv_fit$lambda.min} else {cv_fit$lambda.1se}
                        residuals =  kpop_all_svyd$variables$outcome - predict(cv_fit$glmnet.fit,
                                                                               s = lambda_pass, 
                                                                               newx = X)
                        res_kpop_all = data.frame(min = min(residuals), 
                                                  perc_25 = quantile(residuals, .25), 
                                                  mean = mean(residuals),
                                                  perc_75 = quantile(residuals, .75),
                                                  var = var(residuals))
                        kpop_all_se <- tryCatch(calc_SEs(Y = kpop_all_svyd$variables$outcome,
                                                         residuals = residuals,
                                                         pop_size = nrow(cces),
                                                         sample_size = sum(sample),
                                                         weights = weights(kpop_all_svyd)), 
                                                error = function(e) NA)
                        if(length(kpop_demos_wedu_se) == 1) {
                            kpop_demos_wedu_se <- data.frame(SE_fixed = NA, 
                                                             SE_quasi = NA, 
                                                             SE_linear = NA, 
                                                             SE_chad = NA)
                        }
                        names(kpop_all_se) = tryCatch(paste0("kpop_all_", names(kpop_all_se)),
                                                      error = function(e) NA)
                        
                    }
                    
                    biasbound_r_all = kbal_all_est$biasbound_ratio
                    biasbound_all = kbal_all_est$biasbound_opt
                    
                    rm(kbal_all_est, residuals, X,V, cv_fit)
                } else {
                    kpop_demos <- NA
                    numdims_demos = c(NA) 
                    kpop_demos_se <- data.frame(SE_fixed = NA, 
                                                SE_quasi = NA, 
                                                SE_linear = NA, 
                                                SE_chad = NA)
                    res_kpop_demos = data.frame(min = NA, 
                                                perc_25 = NA, 
                                                mean = NA,
                                                perc_75 = NA,
                                                var = NA)
                    biasbound_r_demos = NA
                    biasbound_demos = NA
                    
                    kpop_demos_wedu <- NA
                    numdims_demos_wedu = c(NA) 
                    kpop_demos_wedu_se <- data.frame(SE_fixed = NA, 
                                                     SE_quasi = NA, 
                                                     SE_linear = NA, 
                                                     SE_chad = NA)
                    res_kpop_demos_wedu =  data.frame(min = NA, 
                                                      perc_25 = NA, 
                                                      mean = NA,
                                                      perc_75 = NA,
                                                      var = NA)
                    biasbound_r_demos_wedu = NA
                    biasbound_demos_wedu = NA
                    
                    kpop_all <- NA
                    numdims_all = c(NA) 
                    kpop_all_se <- data.frame(SE_fixed = NA, 
                                              SE_quasi = NA, 
                                              SE_linear = NA, 
                                              SE_chad = NA)
                    res_kpop_all =  data.frame(min = NA, 
                                               perc_25 = NA, 
                                               mean = NA,
                                               perc_75 = NA,
                                               var = NA)
                    biasbound_r_all = NA
                    biasbound_all = NA
                    
                    #for weights
                    kpop_demos_svyd = survey_design
                    kpop_demos_wedu_svyd = survey_design
                    kpop_all_svyd = survey_design
                }
                
                
                rm(svdK)
                
                ##### return
                kpop_res = list()
                b_out = b_kpop
                # b_reduc = b_kpop_reduc
                kpop_res$sims = data.frame(b_out,
                                           kpop = kpop,
                                           kpop_mf = kpop_mf,
                                           kpop_conv = kpop_conv,
                                           kpop_demos = kpop_demos,
                                           kpop_demos_wedu = kpop_demos_wedu,
                                           kpop_all = kpop_all,
                                           bb = biasbound,
                                           bbr = biasbound_r,
                                           bb_conv = biasbound_conv,
                                           bbr_conv = biasbound_r_conv,
                                           bb_mf = biasbound_mf,
                                           bbr_mf = biasbound_r_mf,
                                           bb_demos = biasbound_demos,
                                           bbr_demos = biasbound_r_demos,
                                           bb_demos_wedu = biasbound_demos_wedu,
                                           bbr_demos_wedu = biasbound_r_demos_wedu,
                                           bb_all = biasbound_all,
                                           bbr_all = biasbound_r_all,
                                           numdims,
                                           numdims_conv,
                                           mfnumdims, 
                                           mf_appended_dims, 
                                           numdims_demos,
                                           numdims_demos_wedu,
                                           numdims_all)
                
                #Starndard Errors:
                
                kpop_res$SEs = data.frame(rake_demos_noeduc_se,
                                          rake_demos_noeduc_se_SVY,
                                          rake_demos_weduc_se,
                                          rake_all_se,
                                          post_stratification_se,
                                          rake_truth_se,
                                          kpop_se,
                                          kpop_conv_se,
                                          kpop_mf_se,
                                          kpop_demos_se,
                                          kpop_demos_wedu_se,
                                          kpop_all_se)
                
                #weights
                kpop_res$weights = list(b = b_out,
                                        kpop_w = weights(kpop_svyd),
                                        #kpop_w_reduc = weights(kpop_svyd_reduc),
                                        kpop_w_conv = weights(kpop_svyd_conv),
                                        kpop_mf_w = weights(kpop_mf_svyd), 
                                        kpop_demos_w = weights(kpop_demos_svyd),
                                        kpop_demos_wedu_w = weights(kpop_demos_wedu_svyd),
                                        kpop_all_w = weights(kpop_all_svyd))
                
                #residuals
                kpop_res$residuals = rbind(b = b_out,
                                           kpop = res_kpop,
                                           # kpop_w_reduc = res_kpop_reduc,
                                           kpop_conv = res_kpop_conv,
                                           kpop_mf = res_kpop_mf,
                                           kpop_demos = res_kpop_demos,
                                           kpop_demos_wedu = res_kpop_demos_wedu,
                                           kpop_all = res_kpop_all,
                                           rake_truth = res_rake_truth,
                                           rake_demos = res_rake_demos_noeduc,
                                           rake_demos_wedu = res_rake_demos_wedu,
                                           rake_all = res_rake_all
                )
                
                ######## Kpop Margins ########
                
                kpop_res$km <- round(cbind(b = b_out/100,
                                           kpop = svymean(margins_formula, kpop_svyd),
                                           # kpop_reduc = svymean(margins_formula, kpop_svyd_reduc),
                                           kpop_conv = svymean(margins_formula, kpop_svyd_conv),
                                           kpop_mf = svymean(margins_formula, kpop_mf_svyd),
                                           kpop_demos = svymean(margins_formula, kpop_demos_svyd),
                                           kpop_demos_wedu = svymean(margins_formula, kpop_demos_wedu_svyd),
                                           kpop_all = svymean(margins_formula, kpop_all_svyd)) * 100,
                                     4)
                
                rm(kpop_svyd, 
                   #kpop_svyd_reduc,
                   kpop_mf_svyd, kpop_svyd_conv, kpop_demos_svyd,
                   kpop_demos_wedu_svyd, kpop_all_svyd)
                
            }
            
            ############################################ OUTPUT
            out = list()
            if(eval_kpop) {
                out$sims = cbind(nsim,
                                 n,
                                 unweighted,
                                 rake_demos_noeduc,
                                 rake_demos_weduc,
                                 rake_all,
                                 post_stratification,
                                 rake_truth,
                                 ht_truth, 
                                 hayek_truth,
                                 kpop_res$sims)
                
                out$SEs = kpop_res$SEs
                out$weights = kpop_res$weights
                out$residuals = kpop_res$residuals
                out$dropped_cells = c(dropped_cells = dropped_cells)
                
                out$sample = c(drop_ps = count, 
                               bad_sample = bad_sample, 
                               #chaningn temporarily to the fuller view from check_nums
                               check = check_nums)
                out$samp_counts = check_2$counts
                
                margin <- round(cbind(sample = svymean(margins_formula, survey_design),
                                      cces =  svymean(margins_formula, cces_svy),
                                      rake_demos_noeduc = svymean(margins_formula,
                                                                  rake_demos_noeduc_svyd),
                                      rake_demos_weduc = svymean(margins_formula,
                                                                 rake_demos_weduc_svyd),
                                      rake_all = svymean(margins_formula,
                                                         rake_all_svyd),
                                      post_stratification = svymean(margins_formula,
                                                                    post_stratification_svyd),
                                      rake_truth = truth_margins) * 100, 5)
                
                #these are just means so let's not multiply by 100
                margin["recode_age",] <- margin["recode_age",]/100 
                margin["recode_agesq",] <- margin["recode_agesq",]/100
                margin["recode_agecubed",] <- margin["recode_agecubed",]/100
                
                margin = cbind(margin,
                               kpop_res$km)
                
                margin <- margin[,grepl("km.kpop_res", colnames(margin))|
                                     !grepl("km.", colnames(margin)) ]
                
            } else {
                margin <- round(cbind(sample = svymean(margins_formula, survey_design),
                                      cces =  svymean(margins_formula, cces_svy),
                                      rake_demos_noeduc = svymean(margins_formula,
                                                                  rake_demos_noeduc_svyd),
                                      rake_demos_weduc = svymean(margins_formula,
                                                                 rake_demos_weduc_svyd),
                                      rake_all = svymean(margins_formula,
                                                         rake_all_svyd),
                                      post_stratification = svymean(margins_formula,
                                                                    post_stratification_svyd),
                                      rake_truth = truth_margins) * 100, 5)
                
                #these are just means so let's not multiply by 100
                margin["recode_age",] <- margin["recode_age",]/100
                margin["recode_agesq",] <- margin["recode_agesq",]/100
                margin["recode_agecubed",] <- margin["recode_agecubed",]/100
                
                out$sims = data.frame(nsim,
                                      n,
                                      unweighted,
                                      rake_demos_noeduc,
                                      rake_demos_weduc,
                                      rake_all,
                                      post_stratification,
                                      rake_truth,
                                      ht_truth, 
                                      hayek_truth)
                
                # out$SEs = data.frame(rake_demos_noeduc_se,
                #                      rake_demos_weduc_se,
                #                      rake_all_se,
                #                      rake_truth_se,
                #                      post_stratification_se)
                # 
                
                out$SEs = data.frame(rake_demos_noeduc_se,
                                     rake_demos_weduc_se,
                                     rake_demos_noeduc_se_SVY,
                                     rake_all_se,
                                     rake_truth_se,
                                     post_stratification_se)
                
                
                out$residuals = rbind(b = NULL,
                                      rake_truth = res_rake_truth,
                                      rake_demos = res_rake_demos_noeduc,
                                      rake_demos_wedu = res_rake_demos_wedu ,
                                      rake_all = res_rake_all
                )
                
                out$dropped_cells = c(dropped_cells = dropped_cells)
                
                out$sample = c(drop_ps = count, 
                               bad_sample = bad_sample, 
                               check = check_nums)
                out$samp_counts = check_2$counts
            } 
            
            out$margins = margin
            
        })
        return(out)
        
    }, mc.cores = detectCores() - cores_saved) 
    
})



################################## Clean Sims ##################
outcome = cces$outcome      
good = which(lapply(sims, function (x) return(class(x))) == "list")
length(good)

#sims[-good]
if(SAVE) {
    save(sims, outcome, tolerance, maxit, increment, min_num_dims, noise, R2_outcome, eval_kpop,
         coefs, coefs_outcome, selection_model, p_include, pS_denom, manual_lambda, lambda_min,
         file = paste0("./res_kpop", eval_kpop, "lambdamin", lambda_min, "man", manual_lambda,
                       "_noise", noise, "_on",
                       Sys.Date(),
                       "_nsims", length(good),
                       ".RData"))
}
sims = sims[good]

#combines all weights across rows but can group by b to get them per iteration
if(eval_kpop) { 
    weights <- lapply(sims, `[[`, 3) %>% bind_rows() 
    residuals <- lapply(sims, `[[`, 4) %>% bind_rows()
    ps_dropped = lapply(sims, `[[`, 5) %>% bind_rows()
    check_samp = lapply(sims, `[[`, 6) %>% bind_rows()
    margins <- lapply(sims, `[[`, 7) 
} else { 
    weights = NULL
    #SEs <- lapply(sims, `[[`, 2) %>% bind_rows()
    residuals <- lapply(sims, `[[`, 3) %>% bind_rows()
    ps_dropped = lapply(sims, `[[`, 4) %>% bind_rows()
    check_samp = lapply(sims, `[[`, 5) %>% bind_rows()
    margins = lapply(sims, `[[`, 6) 
}

##################### eval coverage ####################
coverage <- function(SE, x_bar, truth =NULL, crit_val= qnorm(0.975)) {
    if(is.null(truth)) {
        truth = svymean(~outcome, cces_svy)[1]
    }
    x_upper = x_bar + (SE*crit_val)
    x_lower = x_bar - (SE*crit_val)
    contains_truth = matrix(NA, ncol = ncol(SE), nrow = 1)
    for(i in 1:ncol(x_upper)) {
        contains_truth[,i] = sum((truth <= x_upper[,i] & truth >= x_lower[,i]))/nrow(SE)
    }
    colnames(contains_truth) = colnames(x_bar)
    return(contains_truth)
}

# eval coverage of diff SEs
all_SE_coverage <- function(sims, drop_NA = F, truth = NULL, methods = c("rake|kpop")) {
    est <- lapply(sims, `[[`, 1) %>% bind_rows()
    SEs <- lapply(sims, `[[`, 2) %>% bind_rows()
    est_c = est[grepl(methods, colnames(est))]
    SEs = SEs[grepl(methods, colnames(SEs))]
    # a bit of a pain to drop NAs colwise and get coverage rather than
    # dropping all NA rows and then getting coverage 
    #(unfairly drops rows in for methods that don't have NAs) 
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
            
            
            SEs_temp = na.omit(SEs[grepl(colnames(est_c)[i], colnames(SEs))])
            
            SE_fixed = SEs_temp[grepl("SE_fixed$", colnames(SEs_temp))]
            SE_linear = SEs_temp[grepl("SE_linear$", colnames(SEs_temp))]
            SE_quasi = SEs_temp[grepl("SE_quasi$", colnames(SEs_temp))]
            SE_chad= SEs_temp[grepl("SE_chad$", colnames(SEs_temp))]
            # SE_svy= SEs_temp[grepl("SVY", colnames(SEs_temp))]
            # search = gsub("_se_SVY","", colnames(SE_svy))
            
            coverage_out = cbind(coverage_out, 
                                 rbind(coverage(SE_fixed, est_temp, truth = truth),
                                       coverage(SE_linear, est_temp, truth = truth),
                                       coverage(SE_quasi, est_temp, truth = truth),
                                       coverage(SE_chad, est_temp, truth = truth)))
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
            s = coverage(SE_svy, est_c[,grepl(search, colnames(est_c))])
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
    return(out)
}



######## Bias
est <- lapply(sims, `[[`, 1) %>% bind_rows()
colMeans(est)
est = est[grepl(c("rake|kpop|post|unweighted|h"), colnames(est))]
apply(est, 2, function(x) sum(is.na(x)))
bias = colMeans(est, na.rm = T)
bias = bias - mean(outcome)
bias = data.frame(bias = t(t(bias))*100)
bias = bias %>% arrange(desc(abs(bias)))
round(bias,3)
# kable(round(bias, 3), format = "latex", booktabs = T, 
#       caption = paste0("Bias \\textbf{in Percent} across ", length(good), " sims: All Methods (Target = ", round(mean(outcome),3)*100, ")"))


########## SEs
SE_coverage = all_SE_coverage(sims, truth = mean(outcome), drop_NA = T)
SE_coverage

SEs = lapply(sims, `[[`, 2) %>% bind_rows()
est <- lapply(sims, `[[`, 1) %>% bind_rows()
cols = if(eval_kpop) { c(3:7,10:12,14:20)} else {c(3:ncol(est))}
est =est[, cols]

empirical_SEs <- function(sims, eval_kpop = T, na_rm = F) {
    SEs = lapply(sims, `[[`, 2) %>% bind_rows()
    est <- lapply(sims, `[[`, 1) %>% bind_rows()
    cols = if(eval_kpop) { c(3:7,10:12,14:20)} else {c(3:ncol(est))}
    #cols = if(eval_kpop) { c(3:7,10:12,14:20)} else {c(3:7,10:12)}
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
    avg_SE_out
    
    
    #avg_SE_out = cbind(unweighted = NA, avg_SE_out)
    
    #bootstrapped SEs
    boot_SE = t(as.matrix(apply(est, 2, sd)))
    SE_boot = boot_SE[, colnames(boot_SE) %in% colnames(avg_SE_out)]
    emp_SEs = rbind(avg_SE_out, SE_boot)
    
    return(list(emp_SEs =emp_SEs, 
                boot_SE = boot_SE,
                avg_SE = avg_SE_out) )    
}

emp_SE = empirical_SEs(sims = sims, eval_kpop = eval_kpop, na_rm = T)
emp_SE

########## Dropped Cells

colMeans(ps_dropped)


############# check sample

#samp_check
sum(samp_check$bad_sample)
sum(samp_check$check.fail)
#on average sample has 8/24 cats important in the selection model that have <5% units in it
mean(samp_check$check.leq_5pp)
#on average sample has 3/24 cats important in the selection model that have <1% units in it
mean(samp_check$check.leq_1pp)
counts %>% filter(leq_5pp ==1) %>% group_by(var) %>% summarise(n = n(),
                                                               avg_prop = mean(prop))
########## RES ################
est <- lapply(sims, `[[`, 1) %>% bind_rows()
plot = est


################## results plots
#### ported over from tex
if(eval_kpop) {
    plot_lasso_margin <- plot %>% 
        dplyr::select(unweighted, 
                      rake_demos_noeduc,
                      rake_demos_weduc,
                      rake_all,
                      post_stratification,
                      #post_strat_reduc,
                      #post_strat_all,
                      rake_truth,
                      kpop, 
                      #kpop_reduc, 
                      #kpop_conv,
                      #kpop_mf, 
                      kpop_demos,
                      kpop_demos_wedu,
                      kpop_all, 
                      ht_truth,
                      hayek_truth) %>% 
        pivot_longer(everything(),
                     names_to = "estimator", 
                     values_to = "margin") %>%
        mutate(margin = margin * 100,
               estimator_name = factor(case_when(estimator == "kpop" ~ "kpop",
                                                 #estimator == "kpop_reduc" ~ "kpop\n REDUC",
                                                 estimator == "kpop_mf" ~ "kpop aMF (All)",
                                                 # estimator == "kpop_conv" ~ "kpop Converged",
                                                 estimator == "kpop_demos" ~ "kpop+MF:\n (Demos)",
                                                 estimator == "kpop_demos_wedu" ~ "kpop+MF:\n (Demos+Edu)",
                                                 estimator == "kpop_all" ~ "kpop+MF:\n (All)",
                                                 estimator == "rake_demos_noeduc" ~ "Mean Calibration:\n (Demos)",
                                                 estimator == "rake_demos_weduc" ~  "Mean Calibration:\n (Demos+Edu)",
                                                 estimator == "rake_all" ~ "Mean Calibration:\n (All)",
                                                 estimator == "rake_truth" ~ "Mean Calibration:\n True Selection\nModel",
                                                 estimator == "post_stratification" ~ "Post-Strat Truth",
                                                 estimator == "post_strat_reduc" ~ "Post-Stratification:\n (Reduc)",
                                                 estimator == "post_strat_all" ~ "Post-Strat All",
                                                 estimator == "unweighted" ~ "Unweighted",
                                                 estimator == "ht_truth" ~ "Horvitz-Thompson",
                                                 estimator == "hayek_truth" ~ "Hayek"),
                                       levels = c("Unweighted", 
                                                  "Mean Calibration:\n (Demos)",
                                                  "Mean Calibration:\n (Demos+Edu)",
                                                  "Mean Calibration:\n (All)",
                                                  "Post-Strat Truth", 
                                                  #"Post-Stratification:\n (Reduc)", 
                                                  #"Post-Strat All",
                                                  "kpop",
                                                  "kpop\n REDUC",
                                                  # "kpop Converged",
                                                  #"kpop aMF (All)",
                                                  "kpop+MF:\n (Demos)",
                                                  "kpop+MF:\n (Demos+Edu)",
                                                  "kpop+MF:\n (All)",
                                                  "Mean Calibration:\n True Selection\nModel",
                                                  "Horvitz-Thompson",
                                                  "Hayek"
                                       ) ) )
    
} else {
    plot_lasso_margin <- plot %>% 
        dplyr::select(unweighted, 
                      rake_demos_noeduc,
                      rake_demos_weduc,
                      rake_all,
                      post_stratification,
                      post_strat_reduc,
                      post_strat_all,
                      rake_truth, 
                      ht_truth, 
                      hayek_truth) %>% 
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
                                                 estimator == "post_stratification" ~ "Post-Strat",
                                                 estimator == "post_strat_reduc" ~ "Post-Stratification:\n (Reduc)",
                                                 estimator == "post_strat_all" ~ "Post-Strat All",
                                                 estimator == "unweighted" ~ "Unweighted",
                                                 
                                                 estimator == "ht_truth" ~ "Horvitz-Thompson",
                                                 estimator == "hayek_truth" ~ "Hayek"),
                                       levels = c("Unweighted", 
                                                  "Mean Calibration:\n (Demos)",
                                                  "Mean Calibration:\n (Demos+Edu)",
                                                  "Mean Calibration:\n (All)",
                                                  "Post-Strat", 
                                                  "Post-Stratification:\n (Reduc)", 
                                                  "Post-Strat All",
                                                  "kpop",
                                                  #"kpop Converged",
                                                  #"kpop aMF (All)",
                                                  "kpop+MF:\n (Demos)",
                                                  "kpop+MF:\n (Demos+Edu)",
                                                  "kpop+MF:\n (All)",
                                                  "Mean Calibration:\n True Selection\nModel",
                                                  "Horvitz-Thompson",
                                                  "Hayek"
                                       )))
    
}

#target:
#margin_sim = svymean(~outcome, cces_svy)[1]* 100
#### Box Plot
#options(dplyr.print_max = 1e9)
gg_out = ggplot(data = plot_lasso_margin,
                aes(x = estimator_name, y = margin)) +
    geom_boxplot(alpha = 0.2) +
    geom_hline(yintercept = mean(outcome)*100) +
    theme_bw() +
    xlab("") +
    ylab("Modeled Vote Margin") +
    annotate(geom = "text", x = 0.85, y = mean(outcome)*100+0.25, size = 2.7, angle = 90,
             label = "True Target\nPopulation\nMargin", hjust = 0) +
    ggtitle(paste0(nrow(est)," sims w/avg n_samp =", round(mean(est$n)))) +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))

gg_out
### table
table = plot_lasso_margin %>% 
    mutate(estimator_name = gsub("\n", " ", estimator_name)) %>%
    group_by(estimator_name) %>%
    summarize(
        Bias = mean(margin - margin_sim),
        SE_boot= sd(margin),
        MSE = mean((margin - margin_sim)^2)
    ) %>%
    mutate(
        Bias_Reduc = 1- Bias / Bias[estimator_name == "Unweighted"]
    ) %>%
    arrange(MSE)
table

    if(SAVE) {
        save(sims, outcome, SE_coverage, bias, table, plot_lasso_margin, noise,eval_kpop,emp_SE,
             tolerance, maxit, increment, min_num_dims,
             coefs, coefs_outcome, selection_model, p_include, pS_denom, manual_lambda, lambda_min,
             file = paste0("./res_kpop", eval_kpop, "lambdamin", lambda_min, "man", manual_lambda, 
                           "_noise", noise, "_on",
                           Sys.Date(), 
                           "_nsims", length(good),
                           ".RData"))
    }