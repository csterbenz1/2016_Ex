#Poking around in the application and simulations too see what's going on here

#Possible differenes
# 1. 3way educ interaction in sims versus 6 level educ in app
# 2. sims drops NA vote_2016 (81) as well as any missing (40), application does not


############################# Application ######################
library(tidyverse)
library(survey)
#devtools::install_github("chadhazlett/KBAL", ref = "master") 
library(kbal)

path_data = "/Users/Ciara/Dropbox/kpop/application/data/"
############################## Load Data ##########################
pew <- readRDS(paste0(path_data, "pew_CS.rds"))

pew_srs <- svydesign(ids = ~1, data = pew)

cces <- readRDS(paste0(path_data, "cces_CS.rds"))
#XXXX: drop NA vote_2016 ??
cces <- cces %>%
    filter((CC16_401 == "I definitely voted in the General Election.") &
               !is.na(commonweight_vv_post))
cces_awt <- svydesign(ids = ~1, weights = ~commonweight_vv_post, data = cces)


############################## Run Kbal ##########################


#chad's thing
N=20
a=as.factor(sample(c(1,2,3), size=N, replace=TRUE))
b=as.factor(sample(c(4,5,6), size=N, replace=TRUE))
c = b
(m=model.matrix(~a*b))

(m=model.matrix(~b*c))



# W PID included
kbal_data <- rbind(pew %>% select(recode_age,
                                  recode_female,
                                  recode_race, 
                                  recode_region, 
                                  recode_pid_3way,
                                  recode_educ_3way),
                   cces %>% select(recode_age, 
                                   recode_female, 
                                   recode_race, 
                                   recode_region,
                                   recode_pid_3way, 
                                   recode_educ_3way)) %>%
    model.matrix(as.formula("~."), .)

kbal_data <- kbal_data[,-1]


lm(kbal_data[,1] ~ kbal_data[,-1])

kbal_data_sampled <- c(rep(1, nrow(pew)), rep(0, nrow(cces)))
##### collinearity LA: ONLY an issue with w educ_3way 
#but like wtf it finds that the two most highly correlated cols are Other and White
#with educ_3way no split = all other races, then white x 3 way level of educ so
#I would expect that like one of these categories would be highly correl with race
#wtf does using the interaction lead us to find that race among itself is highly correlated?
#the correlation between whtie and other is still .5 just as with the educ_3way version but we somehow have full rank what the fuck and its still the max but now things aren't multicolinear
#sooo using correlations to check for multicolinearity is the worst plan?

allx <- kbal_data
colnames(kbal_data)
qr(allx)$rank
#AHA removing 5 (white) as well as 11, 12, 13 which are: 3way educ_ college, no college, and post_grad also gets the right rank back to 12. So the correlation method would seem to point in the wrong direction here
#it must be that white is just barely passing as indep from other tho?
qr(allx[,-13])$rank
qr(allx[,-5])$rank
ncol(allx)
#ok so let's see is educ3way redundant withitself or?
#ok so it's only redundant with white! so that makes sense which is why that's what it's dropping though....
#ok yeah so shit fuck the correlation method does not work. the issue is one col of educ_3way is redundant with white which makes perfect sense since its an intrection and probaby like post grad is rare enough that its not independent or whatever
qr(allx[,c(5,11:13)])$rank
qr(allx[,c(11:13)])$rank

qr_X = qr(allx)
allx_update = allx
dropped_cols = NULL
if(qr_X$rank < ncol(allx)) {
    warning("\"allx\" contains collinear columns. Dropping these columns", 
            immediate. = TRUE)
    multicollin = TRUE
}
while(multicollin == TRUE){
    cor = cor(allx_update)
    diag(cor) = 0
    cor[lower.tri(cor)] = 0
    cor = abs(cor)
    View(cor)
    #row of raceOther and col raceWhite has .5 corr
    drop = which(cor == max(cor), arr.ind  =TRUE)[1,2]
    dropped_cols = c(dropped_cols, colnames(allx_update)[drop])
    allx_update = allx_update[,-drop]
    #yeah ok wtf this isn't working? we drop the highest corr but then we still are rank deficient even though we started out only 1 rank deff and we're still 1 rank def?
    #ok so if we drop the COLUMN and not the row then we do drop immediately the right column
    #but droppign white wtf
    if(qr(allx_update)$rank == ncol(allx_update)) {multicollin = FALSE}
    
}
ncol(allx_update)
dropped_cols

#OK SO THIS IS ACTUALLY OK AS LONG AS WE DROP THE COLUMN AND NOT THE ROW... IS THAT ALWAYS THE CASE?... still misleading bc its identifying the correlation between white and other not white and educ_3way... so thats probably just dumb luck and the correlation method really doesnt work


#alright soooo that's a problem, for now let's manually fix it by dropping one column ourselves
qr(kbal_data[,-13])$rank
#drop postgrad
kbal_data <- kbal_data[,-13]
ncol(kbal_data)
qr(kbal_data)$rank

########################## (1) KPOP + W PID #########################

# Estimate
vote_margin <- quote((recode_vote_2016Democrat - recode_vote_2016Republican) /
                         (recode_vote_2016Democrat + recode_vote_2016Republican))



##################### kbal with eudc_3way but not dropping too many 

kbal_b.5x_sims_drop_postgrad <- kbal(allx=kbal_data,
                    sampled = kbal_data_sampled,
                    b = 0.5 * ncol(kbal_data),
                    fullSVD = TRUE,
                    meanfirst = FALSE,
                    incrementby = 5,
        population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
                    sampledinpop = FALSE)

svy_sims_drop_postgrad <- svydesign(~1, data = pew,
                        weights = kbal_b.5x_sims_drop_postgrad$w[kbal_data_sampled ==1])
#without dropping too many of the redundant columns and using educ_3way we still get an estimate of about 7.8 % :(
svycontrast(svymean(~recode_vote_2016, svy_sims_drop_postgrad,
                    na.rm = TRUE), vote_margin)*100


##################### kbal with educ_3way dropping too many

kbal_data_mcissue <- rbind(pew %>% select(recode_age,
                                  recode_female,
                                  recode_race, 
                                  recode_region, 
                                  recode_pid_3way,
                                  recode_educ_3way),
                   cces %>% select(recode_age, 
                                   recode_female, 
                                   recode_race, 
                                   recode_region,
                                   recode_pid_3way, 
                                   recode_educ_3way)) %>%
    model.matrix(as.formula("~."), .)
kbal_data_mcissue <- kbal_data_mcissue[,-1]


kbal_b.5x_sims_mcissue <- kbal(allx=kbal_data_mcissue,
                                     sampled = kbal_data_sampled,
                                     b = 0.5 * ncol(kbal_data_mcissue),
                                     fullSVD = TRUE,
                                     meanfirst = FALSE,
                                     incrementby = 5,
            population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
                                     sampledinpop = FALSE)

svy_sims_mcissue <- svydesign(~1, data = pew,
                            weights = kbal_b.5x_sims_mcissue$w[kbal_data_sampled ==1])

#still high up at 7.2% so...
svycontrast(svymean(~recode_vote_2016, svy_sims_mcissue,
                    na.rm = TRUE), vote_margin)*100



############## Kbal with educ 6 level fully (no multicollin)
kbal_data_application <- rbind(pew %>% select(recode_age,
                                  recode_female,
                                  recode_race, 
                                  recode_region, 
                                  recode_pid_3way,
                                  recode_educ),
                   cces %>% select(recode_age, 
                                   recode_female, 
                                   recode_race, 
                                   recode_region,
                                   recode_pid_3way, 
                                   recode_educ)) %>%
    model.matrix(as.formula("~."), .)
kbal_data_application <- kbal_data_application[,-1]



kbal_b.5x_application <- kbal(allx=kbal_data_application,
                               sampled = kbal_data_sampled,
                               b = 0.5 * ncol(kbal_data_application),
                               fullSVD = TRUE,
                               meanfirst = FALSE,
                               incrementby = 5,
             population.w = cces$commonweight_vv_post / mean(cces$commonweight_vv_post),
                               sampledinpop = FALSE)

svy_application <- svydesign(~1, data = pew,
                                    weights = kbal_b.5x_application$w[kbal_data_sampled ==1])

#still high up at 5.6%
svycontrast(svymean(~recode_vote_2016, svy_application,
                    na.rm = TRUE), vote_margin)*100



#####################################################################################
#####################################################################################
################################# Simulations #####################################
#rm(list = c("cces", "pew", "cces_awt", "kbal_data"))
nsims = 2

# Set up 4 Methods

#### (1) Marginal distributions of age, female, race, and region
formula_demo_rake <- ~recode_age_3way + recode_female + recode_race + recode_region + recode_pid_3way
#### (2) Marginal distributions of age, female, race, region, and education
formula_demo_educ_rake <- ~recode_age_3way + recode_female + recode_race + recode_region +
    recode_educ_3way + recode_pid_3way
#### (3) Post-stratification strata
formula_ps <- ~recode_age_3way + recode_female + recode_race + recode_region + recode_educ_3way + recode_pid_3way
#### (4) Selection model (also rake_demo_truth)
selection_model = as.formula(~recode_female:recode_pid_3way + recode_age:recode_pid_3way + recode_race:recode_region:recode_educ_3way:recode_pid_3way + recode_age + I(recode_age^2))



create_targets <- function (target_design, target_formula) {
    target_mf <- model.frame(target_formula, model.frame(target_design))
    target_mm <- model.matrix(target_formula, target_mf)
    wts <- weights(target_design)
    
    return(colSums(target_mm * wts) / sum(wts))
}

### Post-stratification function
## For now assumes that strata variable is already created and in
## the data set and called "strata"
postStrat <- function(survey, pop_counts, pop_w_col) {
    survey_counts <- survey %>%
        group_by(strata) %>%
        summarize(n = n()) %>%
        ungroup() %>%
        mutate(w_survey = n / sum(n))
    
    pop_counts <- pop_counts %>%
        rename(w_pop = matches(pop_w_col))
    
    post_strat <- pop_counts %>%
        left_join(survey_counts, by = "strata") %>%
        filter(!is.na(w_survey)) %>%
        ## Normalizes back to 1 after dropping
        ## empty cells
        mutate(w_pop = w_pop * 1/sum(w_pop),
               w = w_pop / w_survey) %>%
        select(strata, w)
    
    survey <- survey %>%
        left_join(post_strat)
    
    return(survey)
}

## Function for getting the results
est_margin <- function(outcome, design) {
    svymean(as.formula(paste0("~", outcome)), design, na.rm = TRUE)[1]
}


######### Data 
pew <- readRDS(paste0(path_data, "pew_CS.rds"))
#XXXXX DIFFERENCE 1: dropping missing: 40 total
sum(pew$missing ==1)
pew <- pew %>% filter(missing == 0)

### Load Target Data
cces <- readRDS(paste0(path_data, "cces_CS.rds"))

# XXXXX DIFFERENCE 2: dropping NA 2016 votes: 360 total
# XXX ALSO normalizing weights make sure kbal is not normalizing them again
sum(is.na(cces$recode_vote_2016))
cces <- cces %>%
    filter((CC16_401 == "I definitely voted in the General Election.") &
               !is.na(commonweight_vv_post) &
               !is.na(recode_vote_2016)) %>%
    ## renormalize cces weights
    mutate(commonweight_vv_post = commonweight_vv_post/ mean(commonweight_vv_post))


## Make "strata" variable in CCES and Pew
cces <- bind_cols(cces, cces %>% 
                      unite("strata", all.vars(formula_ps), remove = FALSE) %>%
                      unite("strata_wage", c(all.vars(formula_ps), "recode_age"), 
                            remove = FALSE) %>%
                      select(strata, strata_wage))

pew <- bind_cols(pew, pew %>% 
                     unite("strata", all.vars(formula_ps), remove = FALSE) %>%
                     unite("strata_wage", c(all.vars(formula_ps), "recode_age"),
                           remove = FALSE) %>%
                     select(strata, strata_wage))

######## p_include: prob of being sampled in pew
# Stack data with S = 1 indicating Pew
stack_data <- data.frame(bind_rows(pew, cces), 
                         S = c(rep(1, nrow(pew)), rep(0, nrow(cces))))
mod <- model.matrix(selection_model, data = stack_data)
## Remove columns where Pew missing strata
mod <- mod[, apply(mod[stack_data$S == 1, ], 2, sum) != 0]
## Remove columns where CCES missing Strata
mod <- mod[, apply(mod[stack_data$S == 0, ], 2, sum) != 0]

## Run the model
S_model <- glm(S ~ . - 1, data = data.frame(S = stack_data$S, mod),
               family = binomial(link = "probit"),
               weights = c(rep(1, nrow(pew)), cces$commonweight_vv_post))

coefs <- S_model$coefficients

## Get inclusion probabilities for CCES respondents
p_include = predict(S_model, type = "response")[stack_data$S == 0]
summary(p_include)



######### p_vote: prob of voting rep
cces$recode_vote_2016_binary <- as.numeric(cces$recode_vote_2016 == "Republican")

## Drop the "Other" responses so we can model two-way vote share 
Y_keep <- stack_data$recode_vote_2016 != "Other" & stack_data$S == 0

## Model two-party vote share
### XXXX non integer number of successes? in binomial?
Y_model <- glm(recode_vote_2016_binary ~ . - 1, 
               data = data.frame(recode_vote_2016_binary = cces$recode_vote_2016_binary[cces$recode_vote_2016 != "Other"], mod[Y_keep, ]),
               family = binomial(link = "probit"),
               weights = cces$commonweight_vv_post[cces$recode_vote_2016 != "Other"])

#ONLY REMOVING NAS FOR NOW (so pnorm works below for p_vote)
remove <- bind_rows(is.na(Y_model$coefficients)) %>% colSums(na.rm = TRUE) %>% unlist()
## Don't remove intercept
remove[1] <- 0
mod <- mod[, !remove]

## Rerun model
Y_model <- glm(recode_vote_2016_binary ~ . - 1, 
               data = data.frame(recode_vote_2016_binary = cces$recode_vote_2016_binary[cces$recode_vote_2016 != "Other"], mod[Y_keep, ]), 
               family = binomial(link = "probit"), 
               weights = cces$commonweight_vv_post[cces$recode_vote_2016 != "Other"])

Y_coefs <- Y_model$coefficients

## probability vote republican
p_vote = pnorm(mod %*% Y_coefs)[stack_data$S == 0]

## Add simulated results to CCES data
cces <- cces %>%
    mutate(recode_vote_2016_sim = p_vote)

summary(p_vote)
cor(p_vote, p_include)


####### Build Surveys and Targets
## Make Population Survey Design
## using the CCES, with CCES weights
cces_awt <- svydesign(ids = ~1, weights = ~commonweight_vv_post, data = cces)

## Create population targets for raking from models defined above
targets_demo_rake <- create_targets(cces_awt, formula_demo_rake)
targets_demo_educ_rake <- create_targets(cces_awt, formula_demo_educ_rake)
targets_demo_truth <- create_targets(cces_awt, selection_model)

## Make table of Population Counts for post-stratification
cces_counts <- cces %>%
    group_by(strata) %>%
    summarize(n = n(),
              n_wtd = sum(commonweight_vv_post, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(w = n / sum(n, na.rm = TRUE),
           w_wtd = n_wtd / sum(n_wtd, na.rm = TRUE))

## True (CCES weighted) Margin of the simulated outcome
truth = (0.5 - svymean(~recode_vote_2016_sim, svydesign(~1, data = cces, weights = ~commonweight_vv_post))[1]) * 2 * 100
truth


nsims = 2
set.seed(9345876)
#my.seed <- 9345876
#set.seed( my.seed, kind = "L'Ecuyer-CMRG" )
system.time({
    sims <- mclapply(1:nsims, function(x) {
        
        vote_margin <- quote((recode_vote_2016Democrat - recode_vote_2016Republican) /
                                 (recode_vote_2016Democrat + recode_vote_2016Republican))
        
        ## Draw a random survey from the DGP
        sample <- rbinom(nrow(cces), 1, p_include)
        survey_sim <- cces[sample == 1, c("recode_vote_2016_sim", "recode_age", 
                                          "recode_vote_2016", "strata",
                                          all.vars(formula_ps))]
        survey_design <- svydesign(ids = ~1, data = survey_sim)
        
        ############################################
        ## Unweighted estimate
        ############################################
        
        unweighted <- est_margin("recode_vote_2016_sim", survey_design)
        #(0.5 - unweighted) * 2 * 100
        
        mean(survey_sim$recode_vote_2016_sim) #same as est_margin
        
        unweighted_real <- svycontrast(svymean(~recode_vote_2016, survey_design,
                            na.rm = TRUE), vote_margin)
        
        ############################################
        ## Sample size
        ############################################
        n <- sum(sample)
        
        ############################################
        ## Raking on demographics (no education)
        ############################################
        rake_demos_noeduc_svyd <- calibrate(design = survey_design,
                                            formula = formula_demo_rake,
                                            population = targets_demo_rake,
                                            calfun = "raking")
        
        rake_demos_noeduc <- est_margin("recode_vote_2016_sim", rake_demos_noeduc_svyd)
        #(0.5 - rake_demos_noeduc) * 2 * 100
        
        
        svy_demos_noeduc <- svydesign(ids = ~1, data = survey_sim, 
                                      weights = weights(rake_demos_noeduc_svyd))
        #wildly off 10.85
        rake_demos_noeduc_real <- svycontrast(svymean(~recode_vote_2016, svy_demos_noeduc,
                            na.rm = TRUE), vote_margin)
        
        ############################################
        #### Raking on demographics (with education)
        ############################################
        rake_demos_weduc_svyd <- calibrate(design = survey_design,
                                           formula = formula_demo_educ_rake,
                                           population = targets_demo_educ_rake,
                                           calfun = "raking")
        
        rake_demos_weduc <- est_margin("recode_vote_2016_sim", rake_demos_weduc_svyd)
        #(0.5 - rake_demos_weduc) * 2 * 100
        
        
        svy_demos_weduc <- svydesign(ids = ~1, data = survey_sim, 
                                      weights = weights(rake_demos_weduc_svyd))
        #does ok 3.9
        rake_demos_weduc_real <- svycontrast(svymean(~recode_vote_2016, svy_demos_weduc,
                            na.rm = TRUE), vote_margin)
        
        ############################################
        ## Post-stratification
        ############################################
        post_stratification_svyd = svydesign(~1, 
                                             data = postStrat(survey_sim,
                                                              cces_counts, "w_wtd"), 
                                             weights = ~w)
        
        post_stratification <- est_margin("recode_vote_2016_sim", 
                                          post_stratification_svyd)
        #(0.5 - post_stratification) * 2 * 100
        
        
        svy_ps <- svydesign(ids = ~1, data = survey_sim, 
                                     weights = weights(post_stratification_svyd))
        
        #4.9 ? this is much higher than application
        ps_real <- svycontrast(svymean(~recode_vote_2016, svy_ps,
                            na.rm = TRUE), vote_margin)
    
        ############################################
        ## Kpop goes here
        ############################################
        
        ## Select the covariates for use in Kbal
        ## We might want to specify these as a simulation parameter above
        
        kbal_data <- bind_rows(survey_sim %>% select(recode_age,
                                                     recode_female,
                                                     recode_race,
                                                     recode_region,
                                                     recode_pid_3way,
                                                     recode_educ_3way),
                               cces %>% select(recode_age,
                                               recode_female,
                                               recode_race,
                                               recode_region,
                                               recode_pid_3way,
                                               recode_educ_3way)) %>%
            model.matrix(as.formula("~."), .)
        
        #remove intercept this way
        kbal_data <- kbal_data[,-1]
        
        #another collinearity issue in kbal's way of finding the problematic column is ng
        #right now
        #manually remove it (one of the three educ_3way among whites is rdundant, I remove
        #postgrad)
        kbal_data <- kbal_data[,-13]
        
        kbal_data_sampled <- c(rep(1, nrow(survey_sim)), rep(0, nrow(cces)))
        
        kbal_est <- kbal(allx=kbal_data, 
                         sampled = kbal_data_sampled,
                         b = 0.5 * ncol(kbal_data),
                         incrementby = 5, 
                         meanfirst = FALSE,
                         population.w = cces$commonweight_vv_post,
                         sampledinpop = FALSE,
                         fullSVD = TRUE)
        
        kpop_svyd <- svydesign(~1, data = survey_sim,
                               weights = kbal_est$w[kbal_data_sampled ==1])
        
        kpop <- est_margin("recode_vote_2016_sim", kpop_svyd)
        #not good: -0.0635
        #(0.5 - kpop) * 2 * 100
        
        
        #coming in low again around 1.15
        kpop_real <- svycontrast(svymean(~recode_vote_2016, kpop_svyd,
                            na.rm = TRUE), vote_margin)
        
        svdK = kbal_est$svdK # Save the SVD(K) to re-use...
        
        
        ####### with mc
        # kbal_data_withmc <- bind_rows(survey_sim %>% select(recode_age,
        #                                              recode_female,
        #                                              recode_race,
        #                                              recode_region,
        #                                              recode_pid_3way,
        #                                              recode_educ_3way),
        #                        cces %>% select(recode_age,
        #                                        recode_female,
        #                                        recode_race,
        #                                        recode_region,
        #                                        recode_pid_3way,
        #                                        recode_educ_3way)) %>%
        #     model.matrix(as.formula("~."), .)
        # 
        # #remove intercept this way
        # kbal_data_withmc <- kbal_data_withmc[,-1]
        # 
        # kbal_est_wmc <- kbal(allx=kbal_data_withmc, 
        #                  sampled = kbal_data_sampled,
        #                  b = 0.5 * ncol(kbal_data_withmc),
        #                  incrementby = 5, 
        #                  meanfirst = FALSE,
        #                  population.w = cces$commonweight_vv_post,
        #                  sampledinpop = FALSE,
        #                  fullSVD = TRUE)
        # 
        # kpop_svyd_wmc <- svydesign(~1, data = survey_sim,
        #                            weights = kbal_est_wmc$w[kbal_data_sampled ==1])
        # 
        # kpop_wmc <- est_margin("recode_vote_2016_sim", kpop_svyd_wmc)
        # #3.34 better but not great
        # (0.5 - kpop_wmc) * 2 * 100
        # 
        # 
        # #much better ar 2.49 :( so looks like the mc issue helped us
        # kpop_wmc_real <- svycontrast(svymean(~recode_vote_2016, kpop_svyd_wmc,
        #                     na.rm = TRUE), vote_margin)*100
        # 
        
        ############################################
        ## Kpop with mean first 
        ############################################
        
        ## Select the covariates for use in Kbal
        ## We might want to specify these as a simulation parameter above
        
        # Use SVD already computed in prior one.
        
        kbal_mf_est <- kbal(K.svd = svdK,
                            allx=kbal_data,
                            sampled = kbal_data_sampled,
                            b = .5 * ncol(kbal_data), 
                            incrementby = 5, meanfirst = TRUE,
                            population.w = cces$commonweight_vv_post,
                            sampledinpop = FALSE)
        
        kpop_mf_svy <- svydesign(~1, data = survey_sim,
                             weights = kbal_mf_est$w[kbal_data_sampled ==1])
        
        kpop_mf <- est_margin("recode_vote_2016_sim", kpop_mf_svy)
        #3.34  so this is much better than kpop alone
        #(0.5 - kpop_mf) * 2 * 100
        
        #real is not far off from sim at 3.86
        kpop_mf_real <- svycontrast(svymean(~recode_vote_2016, kpop_mf_svy,
                            na.rm = TRUE), vote_margin)
        #
        
        ############################################
        #### Raking on true model
        ############################################
        rake_truth_svyd <- calibrate(design = survey_design,
                                     formula = selection_model,
                                     population = targets_demo_truth,
                                     calfun = "raking",
                                     ## Giving a little extra tolerance
                                     ## Because was not converging
                                     epsilon = 0.005)
        
        rake_truth <- est_margin("recode_vote_2016_sim", rake_truth_svyd)
        #right on the money as we expect
        #(0.5 - rake_truth)*200
        
        
        svy_truth <- svydesign(ids = ~1, data = survey_sim, 
                                      weights = weights(rake_truth_svyd))
        
        #woah oops 4.2???
        rake_truth_real <- svycontrast(svymean(~recode_vote_2016, svy_truth,
                            na.rm = TRUE), vote_margin)
        ############################################
        #### Check Margin
        ############################################
        
        margins <- round(cbind(cces = svymean(formula_demo_educ_rake, cces_awt),
                    rake_demos_noeduc = svymean(formula_demo_educ_rake,
                                                rake_demos_noeduc_svyd),
                    rake_demos_weduc = svymean(formula_demo_educ_rake, 
                                               rake_demos_weduc_svyd),
                    post_stratification = svymean(formula_demo_educ_rake,
                                                  post_stratification_svyd),
                    kpop = svymean(formula_demo_educ_rake, kpop_svyd),
                    kpop_mf = svymean(formula_demo_educ_rake, kpop_mf_svy),
                    rake_truth = svymean(formula_demo_educ_rake, rake_truth_svyd))
                    * 100, 2)
        # View(as.data.frame(wtf) %>% mutate(err_truth = cces - rake_truth,
        #                               err_rake_noeduc = cces- rake_demos_noeduc,
        #                               err_rake_weduc = cces - rake_demos_weduc,
        #                               err_ps = cces - post_stratification,
        #                err_kpop = cces - kpop,
        #                err_kpop_mf = cces - kpop_mf) %>% select(err_truth,
        #                                                         err_rake_noeduc,
        #                                                         err_rake_weduc,
        #                                                         err_ps,
        #                                                      err_kpop, 
        #                                                      err_kpop_mf))
        #>??? UHH WUT why i ps so off and kpop but demos w educ spot on
        ############################################
        
        out <- list(sims = data.frame(unweighted,
                                      n, 
                                      rake_demos_noeduc, 
                                      rake_demos_weduc, 
                                      post_stratification, 
                                      rake_truth,
                                      kpop,
                                      kpop_mf,
                                      numdims = kbal_est$numdims), 
                    real_outcome = data.frame(unweighted_real,
                                              n,
                                              rake_demos_noeduc_real, 
                                              rake_demos_weduc_real, 
                                              ps_real, 
                                              rake_truth_real,
                                              kpop_real,
                                              kpop_mf_real))
        
        #data.frame(unweighted, n, rake_demos_noeduc, rake_demos_weduc, post_stratification, rake_truth, kpop, kpop_mf, numdims = kbal_est$numdims)
        
        return(out)
        
    }, mc.cores = detectCores() - 1) 
})

# save(sims, file = paste0("./sims_results_updated", 
#                          str_sub(gsub("[[:space:]|[:punct:]]", "_", gsub("[:alpha:]", "", Sys.time())), start = 1, end = -3), "_nsims", nsims, ".RData"))


sims_org <- rbind(sims[[1]]$sims, sims[[2]]$sims) %>% bind_rows()
real <- rbind(sims[[1]]$real_outcome, sims[[2]]$real_outcome)



#### plot
plot_data <- sims_org %>% pivot_longer(everything(), names_to = "estimator", values_to = "margin") %>%
    filter(estimator != "n" & estimator != "numdims") %>%
    mutate(margin = (0.5 - margin) * 2 * 100,
           estimator_name = factor(case_when(estimator == "kpop" ~ "KPop",
                                             estimator == "kpop_mf" ~ "KPop mean first",
                                             estimator == "rake_demos_noeduc" ~ "Raking\nDemos (No Educ)",
                                             estimator == "rake_demos_weduc" ~ "Raking\nDemos (w/ Educ)",
                                             estimator == "rake_truth" ~ "Raking\nTrue Selection\nModel",
                                             estimator == "post_stratification" ~ "Post-Stratification",
                                             estimator == "unweighted" ~ "Unweighted"),
                                   levels = c("Unweighted", "Raking\nDemos (No Educ)", "Raking\nDemos (w/ Educ)", "Post-Stratification", "KPop", "KPop mean first","Raking\nTrue Selection\nModel")))

ggplot(data = plot_data,
       aes(x = estimator_name, y = margin)) +
    geom_boxplot(alpha = 0.2) +
    geom_hline(yintercept = truth) +
    theme_bw() +
    ggtitle("Simulation Results") +
    xlab("Estimator") +
    ylab("Simulated Vote Margin") +
    annotate(geom = "text", x = 0.5, y = truth, label = "True\nTarget\nPopulation\nMargin", hjust = 0) +
    theme(panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1))
