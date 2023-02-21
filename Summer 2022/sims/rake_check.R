data(api)
dclus1 <- svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)
rclus1 <- as.svrepdesign(dclus1)

svymean(~api00, rclus1)
svytotal(~enroll, rclus1)

## population marginal totals for each stratum
pop.types <- data.frame(stype=c("E","H","M"), Freq=c(4421,755,1018))
pop.schwide <- data.frame(sch.wide=c("No","Yes"), Freq=c(1072,5122))

rclus1r <- rake(rclus1, list(~stype,~sch.wide), list(pop.types, pop.schwide))



pop.types <- data.frame(stype=c("E","H","M"), Freq=c(4421,755,1018))


####

pop.pid = data.frame(table(cces$recode_pid_3way))
colnames(pop.pid)[1] = "recode_pid_3way"
pop.female = data.frame(table(cces$recode_female))
colnames(pop.female)[1] = "recode_female"
pop.age = data.frame(table(cces$recode_age_bucket))
colnames(pop.age)[1] = "recode_age_bucket"

rake_test <- rake(survey_design, 
                  sample.margins = list(~recode_pid_3way, ~recode_age_bucket, ~recode_female),
                  population.margins = list(pop.pid, pop.age, pop.female), 
                  control = list(maxit = 500, epsilon = .000000001))

est_mean("outcome", rake_test)*100 - margin_sim
rake_truth*100 - margin_sim

sum(weights(rake_test))
sum(weights(rake_truth_svyd))
sum(weights(post_stratification_svyd))
sum(weights(kpop_svyd))
compare = cbind(rake = weights(rake_test), calibrate = weights(rake_truth_svyd)*44932,
      ps = weights(post_stratification_svyd)*44932/615, 
      kpop = weights(kpop_svyd)*44932/615)
colSums(compare)
compare[1:10,]
#ok so rake and calibrate are indeed identical that's something but the'yre both failiing


##### do they still fail when we change the outcome to be just pR for example
cces$p_R = p_R
cces$p_D = p_D
cces$diff = p_D - p_R
survey_sim <- cces[sample == 1, ]
survey_design <- suppressWarnings(svydesign(ids = ~1, data = survey_sim))
rake_truth_svyd <- try(calibrate(design = survey_design,
                                 formula = selection_model,
                                 population = targets_demo_truth,
                                 calfun = "raking",
                                 epsilon = .000000001), silent = T)
post_stratification_svyd = svydesign(~1, data = postStrat(survey_sim, 
                                                          cces_counts, "w", 
                                                          strata_pass = "strata"),
                                     weights = ~w)
kpop_svyd <- svydesign(~1, data = survey_sim,
                       weights = kbal_est$w[kbal_data_sampled ==1])




(est_mean("p_R", rake_truth_svyd) - mean(cces$p_R))*100
(est_mean("p_D", rake_truth_svyd) - mean(cces$p_D))*100
(est_mean("diff", rake_truth_svyd) - mean(cces$diff))*100

(est_mean("p_R", post_stratification_svyd) - mean(cces$p_R))*100
(est_mean("p_D", post_stratification_svyd) - mean(cces$p_D))*100
(est_mean("diff", post_stratification_svyd) - mean(cces$diff))*100

(est_mean("p_R", kpop_svyd) - mean(cces$p_R))*100
(est_mean("p_D", kpop_svyd) - mean(cces$p_D))*100
(est_mean("diff", kpop_svyd) - mean(cces$diff))*100

#diff of mean is also mean diff just double checking
est_mean("p_D", kpop_svyd)  - est_mean("p_R", kpop_svyd) 
est_mean("diff", kpop_svyd) 

ht_pR = sum((cces[sample ==1, "p_R"]/p_sample))/nrow(cces) 
ht_pD = sum((cces[sample ==1, "p_D"]/p_sample))/nrow(cces) 
ht_diff = sum((cces[sample ==1, "diff"]/p_sample))/nrow(cces) 
#Hayek
hayek_pR = sum((cces[sample ==1, "p_R"]/p_sample))/sum(1/p_sample)
hayek_pD = sum((cces[sample ==1, "p_D"]/p_sample))/sum(1/p_sample)
hayek_diff = sum((cces[sample ==1, "diff"]/p_sample))/sum(1/p_sample)

(ht_pR - mean(cces$p_R))*100
(ht_pD - mean(cces$p_D))*100
(ht_diff - mean(cces$diff))*100

(hayek_pR - mean(cces$p_R))*100
(hayek_pD - mean(cces$p_D))*100
(hayek_diff - mean(cces$diff))*100

#ok yeah interesting... something about the individual level difference then? and dem?
#look at regrs
s_r = summary(lm(p_R ~ recode_pid_3way + recode_age_bucket + recode_female,data = cces))
s_d = summary(lm(p_D ~ recode_pid_3way + recode_age_bucket + recode_female,data = cces))
s_diff = summary(lm(outcome ~ recode_pid_3way + recode_age_bucket + recode_female,data = cces))
s_r$adj.r.squared
s_d$adj.r.squared
s_diff$adj.r.squared



### so why is post-strat perfect and raking not? bc it gets all the cross sections? but i ahd no interacitons in the model like what
unique(survey_sim$strata)
#so let's try raking on all the interactions to replicate:
rake_strata <- try(calibrate(design = survey_design,
                                 formula = ~strata,
                                 population = create_targets(cces_svy, ~strata),
                                 calfun = "raking",
                                 epsilon = .000000000000001), silent = T)

(est_mean("p_R", rake_strata) - mean(cces$p_R))*100
(est_mean("p_D", rake_strata) - mean(cces$p_D))*100
(est_mean("diff", rake_strata) - mean(cces$diff))*100

###yeah WHAT THE ACTUAL FUCK, the interactions SHOULD NOT MATTER they are not in the model but they clearly do
#doign this the other way:
selection_model_2 = as.formula(~recode_pid_3way + recode_age_bucket + recode_female + 
                                   recode_pid_3way*recode_age_bucket + 
                                   recode_pid_3way*recode_female + 
                                   recode_age_bucket*recode_female+
                                   recode_pid_3way*recode_female*recode_age_bucket)

rake_inter <- try(calibrate(design = survey_design,
                             formula = selection_model_2,
                             population = create_targets(cces_svy, selection_model_2),
                             calfun = "raking",
                             epsilon = .000000001), silent = T)
(est_mean("p_R", rake_inter) - mean(cces$p_R))*100
(est_mean("p_D", rake_inter) - mean(cces$p_D))*100
(est_mean("diff", rake_inter) - mean(cces$diff))*100


#let's run a model on our outcome made WITHOUT interactions and see what those coefs look like
test = lm(p_R ~ recode_pid_3way + recode_age_bucket + recode_female + 
                     recode_pid_3way*recode_age_bucket + 
                     recode_pid_3way*recode_female + 
                     recode_age_bucket*recode_female,
                 cces)
summary(test)

test = lm(outcome ~ recode_pid_3way + recode_age_bucket + recode_female + 
              recode_pid_3way*recode_age_bucket + 
              recode_pid_3way*recode_female + 
              recode_age_bucket*recode_female +
              recode_pid_3way*recode_female*recode_age_bucket,
          cces)
summary(test)

test = lm(p_include ~ recode_pid_3way + recode_age_bucket + recode_female + 
              recode_pid_3way*recode_age_bucket + 
              recode_pid_3way*recode_female + 
              recode_age_bucket*recode_female + 
              ,
          cces)
summary(test)
#fucking SHIT BRO they clearly really matter; OK THEN LET'S EXPLICITLY SET THEM TO BE ZERO



cces_expanded = model.matrix(selection_model_2, data = cces)
#needs to be n x p X p x 1 -> coef matrix is p x 1
coefs = matrix(NA,nrow = ncol(cces_expanded), ncol =1 )
rownames(coefs) = colnames(cces_expanded)
coefs[,1] = c(-10, #intercept
              5, #selection of indep pos
              4, #selection of R pos
              -2, #36-50,
              3, #51-64,
              2, #65+,
              6, #male pos
              rep(0, nrow(coefs) - 7))
xbeta = cces_expanded %*% coefs

test = lm(xbeta ~ recode_pid_3way + recode_age_bucket + recode_female + 
       recode_pid_3way*recode_age_bucket + 
       recode_pid_3way*recode_female + 
       recode_age_bucket*recode_female +
       recode_pid_3way*recode_female*recode_age_bucket,
   cces)
summary(test)

p_include = plogis(xbeta)
sum(p_include)
summary(p_include)
p_include = p_include + rnorm(nrow(cces), 0, sd = sd(p_include))

test = lm(p_include ~ recode_pid_3way + recode_age_bucket + recode_female + 
              recode_pid_3way*recode_age_bucket + 
              recode_pid_3way*recode_female + 
              recode_age_bucket*recode_female +
              recode_pid_3way*recode_female*recode_age_bucket,
           family = binomial(link = "logit"),
          cces)
summary(test)
summary(test)$r.squared
#two options to deal with sample size from 3k or so to 500: 1) rescale manaully 2) adjust coefs
#since we can directly adjust the ceofs i will just do that (making the intercept quite neg)

#out = data.frame(raw = c(summary(p_include), sum(p_include)) )

#rownames(out) = c("Min", "25%", "Median", "Mean", "75%", "Max", "Sum")

# kable(out, format = "latex", booktabs = T, 
#       caption = "Sample Inclusion Probabilities")

#################### DESIGN OUTCOME MODEL ##################
#let's first try modeling pR and pD and subtracting and then doing bernoulli; same vars just diff coefs

coefs_outcome = coefs
#pR
coefs_outcome[,1] = c(-1.3, #intercept
                      .1, #  indep pos
                      3, #  R pos
                      -.1, #36-50,
                      -.3, #51-64,
                      .6, #65+,
                      .2 #male pos
)
xbeta_outcome = cces_expanded %*% coefs_outcome
p_R = plogis(xbeta_outcome)
mean(p_R)
sum(p_R)/nrow(cces)
#pD
coefs_outcome[,1] = c(.65, #intercept
                      .07, #  indep pos
                      -5, #  R pos
                      .4,#36-50,
                      .5, #51-64,
                      -.6, #65+,
                      -.3 #male pos
)
xbeta_outcome = cces_expanded %*% coefs_outcome
p_D = plogis(xbeta_outcome)
mean(p_D)
sum(p_D + p_R >1 )/nrow(cces)
vote_diff = p_D - p_R
mean(vote_diff)
s = summary(lm(vote_diff ~ recode_pid_3way + recode_age_bucket + recode_female,data = cces))
s$adj.r.squared
#oof let's try and add some noise
vote_diff = p_D - p_R
mean(vote_diff)
#vote_diff = vote_diff + rnorm(nrow(cces), mean = 0, sd = sd(vote_diff))
mean(vote_diff)
#svymean(~recode_vote_2016, cces_svy)[1] - svymean(~recode_vote_2016, cces_svy)[3]
# test = vote_diff + rnorm(nrow(cces), 1, sd = sd(vote_diff))
# var(test)
# var(vote_diff)
s = summary(lm(vote_diff ~ recode_pid_3way + recode_age_bucket + recode_female,data = cces))
cat(paste("R^2 outcome is", round(s$adj.r.squared,3), "\n"))
cat(paste("Mean outcome (target) is", round(mean(vote_diff)*100,3)))
cces$outcome = vote_diff