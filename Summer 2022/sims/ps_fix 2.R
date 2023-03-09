


post_stratification_svyd = svydesign(~1, data = postStrat(survey_sim, 
                                                          cces_counts, "w", 
                                                          strata_pass = "strata"),
                                     weights = ~w)

post_stratification <- est_mean("outcome", post_stratification_svyd)

unique(cces$strata)
selection_model
#why is this wrong?
svymean(~recode_age_bucket, post_stratification_svyd)
svymean(~recode_age_bucket, cces_svy)
#ok so it's hitting the strata exactly....
cbind(data.frame(svymean(~strata, post_stratification_svyd))$mean,
data.frame(svymean(~strata, cces_svy))$mean )

#but that misses age... which must mean... wait literally how is this possible 
svymean(~strata, cces_svy)
cces = cces %>% mutate(test_strat = paste(recode_pid_3way,recode_female, recode_age_bucket, sep = "_"))
unique(cces$test_strat)
cces_svy <- svydesign(ids = ~1, data = cces)
sample <- rbinom(nrow(cces), 1, p_include)
survey_sim <- cces[sample == 1, ]

cces_counts <- cces %>%
    group_by(test_strat) %>%
    summarize(n = if(!POPW) {n()} else {sum(commonweight_vv_post, na.rm = TRUE)}) %>%
    ungroup() %>%
    mutate(w = n / sum(n, na.rm = TRUE))

post_stratification_svyd = svydesign(~1, data = postStrat(survey_sim, 
                                                          cces_counts, "w", 
                                                          strata_pass = "test_strat"),
                                     weights = ~w)

post_stratification <- est_mean("outcome", post_stratification_svyd)
post_stratification
#WTF THAT FIXED IT?
out = data.frame(ps= data.frame(svymean(~test_strat, post_stratification_svyd)),
      cces  =data.frame(svymean(~test_strat, cces_svy)))
out[,c(1,3)]
#now how are age buckets
svymean(~recode_age_bucket, post_stratification_svyd)
svymean(~recode_age_bucket, cces_svy)


##### reran sims and now ps is completely unbiased so that's... interesting ill have to think baout why you need to include all cat... it may be that you can drop one, but not all 3 and it was dropping all age 18-36 and missing those interactions i guess?



######### NExt issue is why is rake_truth biased wtf?
#is it the epsilon? is it again the droppoing one cat issue
#it's the epislon ok
rake_truth_svyd <- try(calibrate(design = survey_design,
                                 formula = selection_model,
                                 population = targets_demo_truth,
                                 calfun = "raking",
                                 epsilon = .0000001), silent = T)

rake_truth <- tryCatch(est_mean("outcome", rake_truth_svyd), 
                       error = function(e) NA)
rake_truth



#### test
test = unite(cces, 'test', all.vars(formula_ps), remove = F)
unique(test$test)

paste(all.vars(formula_ps_reduc), sep  = "_")
