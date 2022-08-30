#confirm raking is the same as ebal:


library(tidyverse)
library(survey)
library(kbal)

path_data= "/Users/Ciara/Dropbox/kpop/application/data/"


# AUXILIARY INFORMATION (CCES)
cces <- readRDS(paste0(path_data, "cces_new.rds"))

### Drop invalid cases
##### SOON DROP NAS
cces <- cces %>%
    filter((CC16_401 == "I definitely voted in the General Election.") &
               !is.na(commonweight_vv_post) 
    ) %>% 
    mutate(commonweight_vv_post = commonweight_vv_post/ mean(commonweight_vv_post))

pew <- readRDS(paste0(path_data, "pew_new.rds"))
# #in pre: eliminate don't know plan to vote respondents:
pew <- pew %>%
    filter(plan1 %in% c("Plan to vote", "Already voted"))


cces <- cces %>% mutate(
    #36 levels
    recode_educ_pid_race = as.factor(paste(recode_educ_3way,
                                           recode_pid_3way, 
                                           recode_race, sep = ", ")),
    #12 levels
    recode_pid_race = as.factor(paste(recode_pid_3way, 
                                      recode_race, sep = ", ")),
    #9 levels
    recode_educ_pid = as.factor(paste(recode_educ_3way,
                                      recode_pid_3way, sep = ", ")),
    recode_agesq = recode_age*recode_age,
    recode_agecubed = recode_age*recode_age*recode_age,
    #13 levles
    recode_midwest_edu_race = as.factor(
        case_when(recode_region == "Midwest" ~ paste(recode_region,
                                                     recode_race,
                                                     recode_educ_3way, 
                                                     sep = ", "),
                  TRUE ~ "No Split")), 
    #4 levels ( low educated whites in midwest only)
    recode_midwest_wh_edu = factor(case_when(
        (recode_race != "White" | recode_region != "Midwest") ~ "No Split", 
        TRUE ~ as.character(recode_educ_3way)), 
        levels = c("No Split", "No College", "College", "Post-grad")), 
    recode_age_factor = factor(case_when(recode_age <92 ~ as.character(recode_age), 
                                         TRUE ~ "92+")) )


pew<- pew %>% mutate(
    #this has 36 levles
    recode_educ_pid_race = as.factor(paste(recode_educ_3way,
                                           recode_pid_3way, 
                                           recode_race, sep = ", ")),
    #12 levels
    recode_pid_race = as.factor(paste(recode_pid_3way, 
                                      recode_race, sep = ", ")),
    #this has 9 levles
    recode_educ_pid = as.factor(paste(recode_educ_3way,
                                      recode_pid_3way, sep = ", ")),
    recode_agesq = recode_age*recode_age,
    recode_agecubed = recode_age*recode_age*recode_age,
    
    #this has 13 levels
    recode_midwest_edu_race = as.factor(
        case_when(recode_region == "Midwest" ~ paste(recode_region,
                                                     recode_race,
                                                     recode_educ_3way, 
                                                     sep = ", "),
                  TRUE ~ "No Split")),
    #this has 4 levels
    recode_midwest_wh_edu = factor(case_when(
        (recode_race != "White" | recode_region != "Midwest") ~ "No Split", 
        TRUE ~ as.character(recode_educ_3way)), 
        levels = c("No Split", "No College", "College", "Post-grad")),
    recode_age_factor = factor(case_when(recode_age <92 ~ as.character(recode_age), 
                                         TRUE ~ "92+")) )



########## RUN EBAL
kbal_data <- kbal_data <- bind_rows(pew %>% dplyr::select(recode_age,
                                                                 recode_age_bucket,
                                                                 recode_female,
                                                                 recode_race,
                                                                 recode_region,
                                                                 recode_pid_3way,
                                                                 recode_educ,
                                                                 
                                                                 recode_income_5way,
                                                                 recode_relig_6way,
                                                                 recode_born,
                                                                 recode_attndch_4way),
                                    cces %>% dplyr::select(recode_age,
                                                               recode_age_bucket,
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
kbal_data <- kbal_data[,-1]

    
X = kbal_data[, sapply(colnames(kbal_data), function(x) {
    grepl("recode_age_bucket", x) | 
        grepl("recode_female", x) |
        grepl("recode_race", x) |
        grepl("recode_region", x) |
        grepl("recode_pid_3way", x) }) ]
X_all <-  kbal_data[, sapply(colnames(kbal_data), function(x) {
    grepl("recode_age_bucket", x) | 
        grepl("recode_female", x) |
        grepl("recode_race", x) |
        grepl("recode_region", x) |
        grepl("recode_pid_3way", x) | 
        grepl("recode_educ", x) |
        grepl("recode_income_5way", x) | 
        grepl("recode_relig_6way", x) | 
        grepl("recode_born", x) | 
        grepl("recode_attndch_4way", x)  }) ]

target_ebal = c(rep(0, nrow(pew)), rep(1, nrow(cces)) )

ebal_demos_noeduc <- ebalance_custom(Treatment = target_ebal,
                               X = X,
                               base.weight = NULL,
                               norm.constant  = NULL,
                               coefs = NULL ,
                               max.iterations = 200,
                               constraint.tolerance = 1e-6,
                               print.level=0)

ebal_all_vars <- ebalance_custom(Treatment = target_ebal,
                                     X = X_all,
                                     base.weight = NULL,
                                     norm.constant  = NULL,
                                     coefs = NULL ,
                                     max.iterations = 200,
                                     constraint.tolerance = 1e-6,
                                     print.level=0)

############### RUN RAKE #########
formula_rake_demos_noeduc <- ~recode_age_bucket + recode_female + recode_race + 
    recode_region + recode_pid_3way

formula_rake_all_vars <- ~recode_age_bucket + recode_female + 
    recode_race + recode_region + recode_pid_3way + recode_educ +
    recode_income_5way + recode_relig_6way + recode_born + recode_attndch_4way

cces_nowt <- svydesign(ids = ~1, data = cces)
pew_srs <- svydesign(ids = ~1, data = pew)

create_targets <- function(target_design, target_formula) {
    target_mf <- model.frame(target_formula, model.frame(target_design))
    target_mm <- model.matrix(target_formula, target_mf)
    wts <- weights(target_design)
    colSums(target_mm * wts) / sum(wts)
}

targets_rake_demos_noeduc <- create_targets(cces_nowt, 
                                            formula_rake_demos_noeduc)
targets_rake_all_vars <- create_targets(cces_nowt, 
                                       formula_rake_all_vars)

rake_demos_noeduc <- calibrate(design = pew_srs,
                               formula = formula_rake_demos_noeduc,
                               population = targets_rake_demos_noeduc,
                               calfun = "raking")

rake_all_vars <- calibrate(design = pew_srs,
                           formula = formula_rake_all_vars,
                           population = targets_rake_all_vars,
                           calfun = "raking")




# compare weights ---------------------------------------------------------
#samp <- round((runif(10)*10), 0)
sum(ebal_demos_noeduc$w[1:nrow(pew)])
length(ebal_demos_noeduc$w)
#our ebal weights sum to N... is that what we want?
e_w <- ebal_demos_noeduc$w
sum(e_w/mean(e_w))

r_w <- weights(rake_demos_noeduc)
sum(r_w)
r_w[1:10]
sum(r_w/mean(r_w))

#ok no so this seems right now just comapre them noramized the same
e_w <- e_w/(sum(e_w))
sum(e_w)
cbind(e_w[1:10], r_w[1:10])
sum(round(e_w - r_w, 5))
#OK SO IT IS DOING THE EXACT SAME THANK GOD
#can also check the rake all weights
e_w_all <- ebal_all_vars$w/sum(ebal_all_vars$w)
r_w_all <- weights(rake_all_vars)
sum(round(e_w_all - r_w_all, 5))
#kewwwwl 


e_w_all <- ebal_all_vars$w/sum(ebal_all_vars$w)
r_w_all <- weights(rake_all_vars)

cbind(e_w_all[1:5], r_w_all[1:5], weights_ebal[1:5])



####### compare on margins:  ######
#jsut for extra assurance that svydesign is correctly reweighting the weights so it's not an issue
rake_demos_noeduc <- svydesign(~1, data = pew, 
                               weights = weights(rake_demos_noeduc))
ebal_demos_noeduc <- svydesign(~1, data = pew, 
                               weights = ebal_demos_noeduc$w)
rake_all_vars <- svydesign(~1, data = pew, 
                               weights = weights(rake_all_vars))
ebal_all_vars <- svydesign(~1, data = pew, 
                               weights = ebal_all_vars$w)



margins_formula <- ~recode_vote_2016 +
    recode_pid_3way + 
    recode_female + recode_race +recode_region + recode_educ +recode_relig_6way+
    recode_born + recode_attndch_4way + recode_income_5way +
    recode_age_bucket+recode_age+recode_agesq+recode_agecubed+recode_age_factor+
    recode_race_educ_reg + recode_educ_wh_3way + 
    recode_educ_pid_race +
    recode_pid_race + 
    recode_educ_pid +
    recode_midwest_edu_race +
    recode_midwest_wh_edu

margins_formula_cces <- ~recode_vote_2016 +
    recode_pid_3way + 
    recode_female + recode_race +recode_region + recode_educ + recode_relig_6way + 
    recode_born + recode_attndch_4way + recode_income_5way +
    recode_age_bucket + recode_age +recode_agesq+recode_agecubed+recode_age_factor +
    recode_race_educ_reg + recode_educ_wh_3way + 
    recode_educ_pid_race +
    recode_pid_race + 
    recode_educ_pid +
    recode_midwest_edu_race +
    recode_midwest_wh_edu


margins <- round(cbind(
    pew = svymean(margins_formula, pew_srs),
    cces_margins = svymean(margins_formula_cces, cces_nowt, 
                           na.rm = TRUE),
    ebal_demos_noeduc = svymean(margins_formula, ebal_demos_noeduc),
    rake_demos_noeduc = svymean(margins_formula, rake_demos_noeduc),
    ebal_all_vars = svymean(margins_formula, ebal_all_vars),
    rake_all_vars = svymean(margins_formula, rake_all_vars)) * 100, 5)
#we multiplied the recode_age margins by 100 unnecc above bc it' just a straight mean
#so let's undo that

rownames(margins) <- sapply(rownames(margins), 
                            function(x) (
                                if(grepl("recode", x)) {
                                    substr(x, 8,nchar(x) )
                                } else {
                                    x
                                }
                            ) )
margins["age",] <- margins["age",]/100 
margins["agesq",] <- margins["agesq",]/100
margins["agecubed",] <- margins["agecubed",]/100


margins_diff <- as.data.frame(margins) %>% 
    mutate(pew = cces_margins - pew,
           ebal_demos_noeduc = cces_margins - ebal_demos_noeduc,
           rake_demos_noeduc = cces_margins- rake_demos_noeduc,
           ebal_all_vars = cces_margins- ebal_all_vars,
           rake_all_vars = cces_margins - rake_all_vars)
rownames(margins_diff) <- rownames(margins)

#get cces proportions:

#weighted average here: compensating for the fact we *100 by 100 above and
#we will take the straight average below so /100*n
#we already have p-q, so here we are just q*(p-q)
# p = margins[(grep("age_factor", rownames(margins))), ]
# q = margins[grep("age_factor", rownames(margins)), "cces_margins"]
# colMeans(abs(p-q)*q/100*length(q))
# colMeans(abs((p-q)*q)/100*length(q))

margins_diff[(grep("age_factor", rownames(margins_diff))), ] <-
    (margins_diff[(grep("age_factor", 
                        rownames(margins_diff))),
                  ])* (margins[grep("age_factor", rownames(margins)), "cces_margins"]/100)*length(levels(cces$recode_age_factor))
#### summarize
#HAHAHHAHAHAHA i dont know how to regex but it works. this is horrific I should switch to str_extract now that i belatedly discovered it
wmabserr <- margins_diff %>% 
    mutate(var = sapply(rownames(margins_diff), 
                        function(y) {
                            substr(y, 1, 
                                   sapply(rownames(margins_diff), 
                                          function(x) {
                                              if(grepl("income",x)) {
                                                  11 #income_5way
                                                  #this is so ridiculous WHY (sunk cost ugh)
                                              } else if(grepl("xgb",x) |
                                                        (grepl("age",x) & 
                                                         !grepl("bucket",x) &
                                                         !grepl("factor",x)) ){
                                                  nchar(x) #keep all xgb stuff
                                                  #THIS IS SO STUPID
                                              } else if(grepl("factor",x)) {
                                                  10 #age_factor
                                              } else {
                                                  (gregexpr("[A-Z]|[1-9][1-9] |[1-9][1-9]\\+|[1-9]\\-", x)[[1]][1]-1)
                                              }
                                          })[y])
                        })) %>%
    group_by(var) %>% summarize(pew_unweighted = mean(abs(pew)), 
                                ebal_demos_noeduc = mean(abs(ebal_demos_noeduc)), 
                                rake_demos_noeduc = mean(abs(rake_demos_noeduc)),
                                ebal_all_vars = mean(abs(ebal_all_vars)),
                                rake_all_vars = mean(abs(rake_all_vars))  )
#for kable aethetics later this helps
wmabserr$var[wmabserr$var == "educ"]<- "educ_6way"

wmabserr_nooutcomes <- wmabserr %>% 
    filter(!(var %in% c("agecubed") |grepl("xgb", var) | grepl("vote", var))) %>%
    arrange(nchar(var))
if(is.null(order)) {
    order = c("female", "pid_3way", "race", "region", "educ_6way", "income_5way", 
              "relig_6way", 
              "pid_race", "educ_pid",
              "educ_pid_race", "race_educ_reg", "educ_wh_3way", "midwest_wh_edu",
              "midwest_edu_race",
              "agesq", "age_factor")
}

wmabserr_nooutcomes <- wmabserr_nooutcomes[as.numeric(sapply(order, function(x) which(wmabserr_nooutcomes$var == x))),]

View(wmabserr)
View(margins_diff)
