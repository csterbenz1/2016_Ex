
count = 0
est_mean <- function(outcome, design) {
    svymean(as.formula(paste0("~", outcome)), design, na.rm = TRUE)[1]
}
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
count = 0

out = replicate(500, {
    
    selection_model = as.formula(~recode_pid_3way + recode_age_bucket + recode_female + recode_race +
                                     recode_relig_6way +
                                     recode_pid_3way:recode_age_bucket
                                 # +recode_age:recode_pid_3way +
                                 #recode_female*recode_pid_3way
    )
    
    cces_expanded = model.matrix(selection_model, data = cces)
    coefs = matrix(NA,nrow = ncol(cces_expanded), ncol =1 )
    rownames(coefs) = colnames(cces_expanded)
    coefs
    coefs[,1] = c(-5, #intercept -5 w race
                  2, #selection of indep pos
                  3, #selection of R pos
                  .15, #36-50,
                  .1, #51-64,
                  .2, #65+,
                  .7, #male
                  .2, #hisp
                  -.2, #other
                  3, #white
                  -3, #relig jewish
                  -.1, #catholic
                  -.4, #not relig
                  -.1, #other race
                  .2,#ind x 36-50
                  1, #rep x 36-50, 
                  .5, #ind x 51-64,
                  .5, #rep x 51-64,
                  -.2, #ind x 65+
                  1 #rep x 65+
    )
    
    xbeta = cces_expanded %*% coefs
    p_include = plogis(xbeta)
    sum(p_include)
    #most likely to select: rep, 65+, male, protestants
    # cces[which(p_include == max(p_include)),
    #      c("recode_pid_3way", "recode_age_bucket", "recode_female", "recode_relig_6way")]
    # #summary(p_include)
    sample <- rbinom(nrow(cces), 1, p_include)
    survey_sim <- cces[sample == 1, ]
    survey_sim %>% group_by(recode_pid_3way, recode_age_bucket) %>% count() %>% mutate(n = n/sum(sample))
    if(length(unique(survey_sim$recode_race)) == 1) { 
        stop("only one cat of race included")
    }
    
    s = survey_sim %>% group_by(recode_pid_3way, recode_age_bucket, recode_race, 
                                recode_relig_6way) %>% count() %>%
        mutate(n_s = round(n/nrow(survey_sim), 3))
    c = cces %>% group_by(recode_pid_3way, recode_age_bucket, recode_race,
                          recode_relig_6way) %>% count() %>%
        mutate(n_c = round(n/nrow(cces), 3))
    nrow(c) - nrow(s)
    
    # s = s %>% mutate(pid_age = paste0(recode_pid_3way, "-", recode_age_bucket))
    # c = c %>% mutate(pid_age = paste0(recode_pid_3way, "-", recode_age_bucket))
    # cat("The Sample loses", nrow(c) - nrow(s), "categories:",
    #     c$pid_age[!(c$pid_age %in% s$pid_age)], "\n" )
    count = (nrow(c) - nrow(s) > 0)
     
     #################### DESIGN OUTCOME MODEL ##################
     coefs_outcome = coefs
     
     coefs_outcome[,1] = c(4.4, #intercept #5 w race
                           -.5, #selection of indep pos
                           -.8, #selection of R pos
                           -1, #36-50,
                           -.1, #51-64,
                           -.2, #65+, 
                           -.4, #male #.4 for male only selection model
                           -1, #hisp
                           -.7, #other
                           -2, #white
                           3, #relig jewish
                           -.1, #catholic
                           -.1, #not relig
                           -.1, #other
                           -2.7,#ind x 36-50
                           -2.8, #rep x 36-50, 
                           -.5, #ind x 51-64,
                           -.5, #rep x 51-64,
                           .5, #ind x 65+
                           -.4 #rep x 65+
     )
     
     coefs_outcome = coefs_outcome/10
     xbeta_outcome = cces_expanded %*% coefs_outcome
     summary(xbeta_outcome)
     if(!bern) {
         #cat(paste("Adding sd(outcome)*",round(noise, 3), "\n")) 
         #set.seed(123123)
         xbeta_outcome = xbeta_outcome + rnorm(nrow(cces), mean = 0, sd = sd(xbeta_outcome)*noise)
         summary(xbeta_outcome)
         
     }
     beyond_support = 0
     if(min(xbeta_outcome) <0 | max(xbeta_outcome) > 1) {
         warning("Outcome beyond prob support for some units when noise is added",
                                                                 immediate. = T)
         beyond_support = 1
     }
     #pD: max should be dem, jewish, female, 18-35
     # cces[which(xbeta_outcome == max(xbeta_outcome)),
     #      c("recode_pid_3way", "recode_age_bucket", "recode_female", "recode_relig_6way")]
     
     #(mean(xbeta_outcome[sample ==1]) - mean(xbeta_outcome))*100
     # cat(paste("Range of outcome w/noise is\n"))
     # cat(paste(summary(xbeta_outcome), "\n"))
     # s = summary(lm(update(selection_model, xbeta_outcome ~ .),data = cces))
     # R2_outcome = s$adj.r.squared
     # cat(paste("R^2 outcome is", round(s$adj.r.squared,3), "\n"))
     # cat(paste("Mean scaled outcome (target) is", round(mean(xbeta_outcome)*100,3)))
     # cat(paste("\nCorr of sampling prob and outcome ", round(cor(xbeta_outcome, p_include),3)))
     #bernoulli draw
     cces$outcome = xbeta_outcome
     margin_sim = mean(cces$outcome)
     cces_svy <- suppressWarnings(svydesign(ids = ~1, data = cces))
     targets_demo_truth <- create_targets(cces_svy, selection_model)
     formula_ps <- selection_model
     cces = cces %>% mutate(strata = paste(recode_pid_3way,
                                           recode_female, 
                                           recode_age_bucket,
                                           recode_race,
                                           recode_relig_6way,
                                           sep = "_"))
     cces_counts <- cces %>%
         group_by(strata) %>%
         summarize(n = if(!POPW) {n()} else {sum(commonweight_vv_post, na.rm = TRUE)}) %>%
         ungroup() %>%
         mutate(w = n / sum(n, na.rm = TRUE))
     
     survey_sim <- cces[sample == 1, ]
     survey_design <- suppressWarnings(svydesign(ids = ~1, data = survey_sim))
     unweighted <- est_mean("outcome", survey_design)
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
    
     rake_truth <- tryCatch(est_mean("outcome", rake_truth_svyd), 
                            error = function(e) NA)
     
     margin_sim = svymean(~outcome, cces_svy)[1]
     
     missing_strata <- unique(cces$strata)[!(unique(cces$strata) %in%
                                                 unique(survey_sim$strata))]
     
     dropped_cells = cces %>% filter(strata %in% missing_strata) %>% group_by(strata) %>% count()
     sum(dropped_cells$n)
     # dropped_cells = data.frame(sum = sum(dropped_cells$n), strata = paste(dropped_cells$strata, collapse = " | "))
     # 
     #note that we no longer have the issue of pew having strata that cces doesn bc
     #we use survey_sim a sample of cces so it will never have diff strata
     
     post_stratification_svyd = svydesign(~1, data = postStrat(survey_sim, 
                                                               cces_counts, "w", 
                                                               strata_pass = "strata", 
                                                               warn = F),
                                          weights = ~w)
     
     post_stratification <- est_mean("outcome", post_stratification_svyd)
     
    c(drop = count, 
      beyond_support = beyond_support,
      drop_ps = sum(dropped_cells$n),
      rake = rake_truth, 
      margin_sim = margin_sim,
      ps = post_stratification, 
      unweighted = unweighted)
     
}, simplify = F)

out
out =out %>% bind_rows()
mean(out$rake.outcome - out$margin_sim.outcome)*100
mean(out$unweighted.outcome- out$margin_sim.outcome)*100
mean(out$ps.outcome- out$margin_sim.outcome)*100
  out$drop
sum(out$beyond_support)
out$drop_ps

out$rake
out$drop





##########################################################
# find what vars covary highly with the interactions and then add them to try
#and force ps to drop cells that are important

colnames(cces)
la = cces %>% select(recode_age_bucket,
                     recode_female,
                     recode_race,
                     recode_region,
                     recode_pid_3way,
                     recode_educ,
                     
                     recode_income_5way,
                     recode_relig_6way,
                     recode_relig,
                     recode_born,
                     recode_attndch_4way) %>% mutate(pid_age = paste0(recode_pid_3way, "-", recode_age_bucket))

plot = la %>% group_by(  recode_pid_3way,recode_income_5way) %>% count()

plot = la %>% mutate(pid_age = paste0(recode_pid_3way, "-", recode_age_bucket),
                     pid_income = paste0(recode_pid_3way, "-", recode_income_5way),
                     pid_educ = paste0(recode_pid_3way, "-", recode_educ),
                     pid_region = paste0(recode_pid_3way, "-", recode_region),
                     pid_attnd = paste0(recode_pid_3way, "-", recode_attndch_4way),
                     pid_born = paste0(recode_pid_3way, "-", recode_born)
                     )
ggplot() +
    geom_histogram(aes(y= pid_income ),stat= "count", data = plot)

ggplot() +
    #geom_histogram(aes(y= recode_income_5way, fill = recode_pid_3way),stat= "count", data = la)
    geom_bar(aes(recode_age_bucket, fill = recode_pid_3way),position="fill", data = la)
