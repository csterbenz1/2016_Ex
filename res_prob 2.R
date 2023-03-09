#####################
##### Defien outcome
#cces = cces %>% mutate(outcome = rbinom(nrow(cces), 1, mod_cces_on_cces_pR))
k = 1
cces$outcome = plogis((k*var(cces$diff_cces_on_cces) + abs(cces$diff_cces_on_cces))*sign(cces$diff_cces_on_cces))
#cces$outcome = cces$diff_cces_on_cces


### Look at correlations and scale of p_include
summary(p_include)
sum(p_include)
#p_include = p_include*(n_sample/sum(p_include))
#sum(p_include)

cor(cces$diff_cces_on_cces, p_include)
cor(cces$outcome, p_include)

#cor(cces$outcome, cces$recode_age)
#cor(cces$outcome, cces$recode_age*cces$recode_age)
cor(cces$outcome, cces$centered_age)
cor(cces$outcome, cces$centered_age*cces$centered_age)

#################### Targets ###################
if(POPW) {
    cces_svy <- svydesign(ids = ~1, weights = ~commonweight_vv_post, data = cces)
} else {
    cces_svy <- svydesign(ids = ~1, data = cces)
}
#margin_sim = svymean(~diff_cces_on_cces, cces_svy)[1]* 100
margin_sim = svymean(~outcome, cces_svy)[1]* 100
margin_sim
############
#Create Sample
sample <- rbinom(nrow(cces), 1, p_include)

survey_sim <- cces[sample == 1, ]
survey_design <- suppressWarnings(svydesign(ids = ~1, data = survey_sim))

unweighted <- est_mean("outcome", survey_design)
unweighted*100
#the unweighted is not that biased for the plogis(diff)
unweighted*100 - margin_sim

n <- sum(sample)
n
sum(sample)
sum(p_include)

####### Raking Svys
rake_demos_noeduc_svyd <- calibrate(design = survey_design,
                                    formula = formula_rake_demos_noeduc,
                                    population = targets_rake_demos_noeduc,
                                    calfun = "raking")
est_mean("outcome", rake_demos_noeduc_svyd)*100 - margin_sim

rake_demos_weduc_svyd <- calibrate(design = survey_design,
                                   formula = formula_rake_demos_weduc,
                                   population = targets_rake_demos_weduc,
                                   calfun = "raking")
rake_demos_weduc <- est_mean("outcome", rake_demos_weduc_svyd)
rake_demos_weduc*100 - margin_sim

rake_all_svyd <- calibrate(design = survey_design,
                           formula = formula_rake_all_vars,
                           population = targets_rake_all_vars,
                           calfun = "raking")
rake_all <- est_mean("outcome", rake_all_svyd)
rake_all*100 - margin_sim


####### R^2s
rake_demos_noeduc_lm = summary(lm(update(formula_rake_demos_noeduc, outcome ~ .), 
                                  data = rake_demos_noeduc_svyd$variables))
rake_demos_weduc_lm = summary(lm(update(formula_rake_demos_weduc, outcome ~ .), 
                                 data = rake_demos_weduc_svyd$variables))
rake_all_lm = summary(lm(update(formula_rake_all_vars, outcome ~ .), 
                         data = rake_all_svyd$variables))

R2 = data.frame((rake_demos_noeduc_lm)$adj.r.squared,
                (rake_demos_weduc_lm)$adj.r.squared,
                (rake_all_lm)$adj.r.squared)

R2



################### Look at Coefs

#within the sample
Y_s = lm(outcome ~ recode_age_bucket + recode_female + recode_race + recode_region + recode_pid_3way + 
           recode_pid_3way * poly(centered_age, 2) + recode_female * recode_pid_3way, data = rake_all_svyd$variables)
#summary(Y)

S_s = lm(p_include[sample ==1] ~ recode_age_bucket + recode_female + recode_race + recode_region + recode_pid_3way + 
           recode_pid_3way * poly(centered_age, 2) + recode_female * recode_pid_3way, data = rake_all_svyd$variables)
#summary(S)

round(cbind(summary(Y_s)$coefficients[,c(1,4)], summary(S_s)$coefficients[,c(1,4)]),3)



#### entire pop;
#rbinom(nrow(cces),1,cces$mod_cces_on_cces_pR)
Y = lm(outcome ~ recode_age_bucket + recode_female + recode_race + recode_region + recode_pid_3way + 
           recode_pid_3way * poly(centered_age, 2) + recode_female * recode_pid_3way, data = cces)
#summary(Y)

S = lm(p_include ~ recode_age_bucket + recode_female + recode_race + recode_region + recode_pid_3way + 
           recode_pid_3way * poly(centered_age, 2) + recode_female * recode_pid_3way, data = cces)
#summary(S)

round(cbind(summary(Y)$coefficients[,c(1,4)], summary(S)$coefficients[,c(1,4)]),3)
round(cbind(summary(Y)$coefficients[,c(1)], summary(S)$coefficients[,c(1)]),3)
round(cbind(summary(Y)$coefficients[,4], summary(S)$coefficients[,4]),3)
lasso_include_coefs

update(formula_rake_demos_noeduc, update( selection_model, outcome ~ .))

update(update( outcome ~ ., formula_rake_demos_noeduc), outcoem ~ . + selection_model)


update(selection_model, formula_rake_demos_noeduc + .)
reformulate(formula_rake_demos_noeduc, names(lasso_include_coefs))