allmargins <- function(pop_weights = TRUE,
                       b = NULL,
                       margins_formula = NULL, 
                       order = NULL,
                       tolerance = 1e-06,
                       pop_weighted_avg = TRUE,
                       maxit = 500) {
    
    pew_nowt <- suppressWarnings(svydesign(ids = ~1, data = pew))
    
    missing_strata <- unique(cces$strata)[!(unique(cces$strata) %in%
                                                unique(pew$strata))]
    missing_strata_reduc <- unique(cces$strata_reduc)[!(unique(cces$strata_reduc) %in%
                                                            unique(pew$strata_reduc))]
    
    missing_strata_all <- unique(cces$strata_all)[!(unique(cces$strata_all) %in%
                                                        unique(pew$strata_all))]
    
    problem <- unique(pew$strata)[!(unique(pew$strata) %in%
                                        unique(cces$strata))]
    problem_reduc <- unique(pew$strata_reduc)[!(unique(pew$strata_reduc)
                                                %in%
                                                    unique(cces$strata_reduc))]
    
    problem_all <- unique(pew$strata_all)[!(unique(pew$strata_all) %in%
                                                unique(cces$strata_all))]
    
    kbal_data_sampled <- c(rep(1, nrow(pew)),
                           rep(0, nrow(cces)))
    
    
    if(pop_weights & maxit == 500 & tolerance == 1e-06) {
        load(paste0(path_weights,"bmaxvar_cat_par_wTRUE_m500_t1e-06_2021-03-30.Rdata"))
        index = 1
        b = out_cat$b
        cat("loading: bmaxvar_cat_par_wTRUE_m500_t1e-06_2021-03-30 \n")
        
    } else if(pop_weights & maxit == 800 & tolerance == 1e-04) {
        #load("./weights/allB_cat_par_wTRUE_m800_t1e-04_2021-03-23.Rdata")
        #cat("loading: allB_cat_par_wTRUE_m800_t1e-04_2021-03-23 \n")
        stop("Not supported currenty")
        b_run <- out_cat$b
        index = which(b_run == b)
        if(length(index) == 0 ) {
            stop("Must enter b in range of b_run")
        }    
    } else if(!pop_weights & tolerance == 1e-06 & maxit == 500) {
        load(paste0(path_weights,"bmaxvar_cat_par_wFALSE_m500_t1e-06_2021-03-30.Rdata"))
        index = 1
        b = out_cat$b
        cat("loading: bmaxvar_cat_par_wFALSE_m500_t1e-06_2021-03-30 \n")
        
    } else if(!pop_weights & tolerance == 1e-04 & maxit == 800) {
        #load("/Users/Ciara/Documents/Cloud Documents/Hazlett:Hartman RA/2016 Election/allB_cat_par_wFALSE_m800_t1e-04_2021-03-23.Rdata")
        #cat("loading: allB_cat_par_wFALSE_m800_t1e-04_2021-03-23.Rdata")
        stop("Not supported currenty")
        b_run <- out_cat$b
        index = which(b_run == b)
        if(length(index) == 0 ) {
            stop("Must enter b in range of b_run")
        }
    }
    
    #updated for the only bmaxvar 
    kpop_w <- weights_cat[[index]]$kpop_w
    kpop_conv_w <- weights_cat[[index]]$kpop_w_conv
    kpop_mf_w <- weights_cat[[index]]$kpop_mf_w
    #kpop_demos_w <- weights_cat[[index]]$kpop_demos_w
    #kpop_demos_wedu_w <- weights_cat[[index]]$kpop_demos_weduc_w
    #kpop_all_w <- weights_cat[[index]]$kpop_all_w
    #temproray until rerun weights
    kpop_demos_w <- weights_cat[[index]]$kpop_contraint_w
    kpop_demos_wedu_w <- weights_cat[[index]]$kpop_contraint_weduc_w
    kpop_all_w <- weights_cat[[index]]$kpop_full_constraint_w
    
    
    b_lab = paste0("b=", round(b, 2))
    
    colnames_bottom = c("Variable", "Orig", "","(Demos)", "(Demos)",
                        "(Demos+Edu)", "(Demos+Edu)", "(All)",  "(All)",
                        "Retro", "Reduc")
    
    if(pop_weights) {
        cces_svy <- svydesign(ids = ~1, weights = ~commonweight_vv_post,
                              data = cces)
    } else {
        cces_svy <- suppressWarnings(svydesign(ids = ~1,
                                               data = cces))
    }
    
    targets_rake_demos_noeduc <- create_targets(cces_svy,
                                                formula_rake_demos_noeduc)
    targets_rake_demos_weduc <- create_targets(cces_svy, formula_rake_demos_weduc)
    targets_rake_all_vars <- create_targets(cces_svy, 
                                            formula_rake_all_vars)
    targets_retrospective <- create_targets(cces_svy, formula_retrospective)
    targets_ps <- svytable(formula = ~strata,
                           design = subset(cces_svy, !(strata %in%
                                                           missing_strata)))
    targets_ps_reduc <- svytable(formula = ~strata_reduc,
                                 design = subset(cces_svy, !(strata_reduc %in%
                                                                 missing_strata_reduc)))
    targets_ps_all <- svytable(formula = ~strata_all,
                               design = subset(cces_svy, !(strata_all %in%
                                                               missing_strata_all)))
    
    rake_demos_noeduc <- calibrate(design = pew_nowt,
                                   formula = formula_rake_demos_noeduc,
                                   population = targets_rake_demos_noeduc,
                                   calfun = "raking")
    rake_demos_noeduc <- svydesign(~1, data = pew, 
                                   weights = weights(rake_demos_noeduc))
    
    rake_demos_weduc <- calibrate(design = pew_nowt,
                                  formula = formula_rake_demos_weduc,
                                  population = targets_rake_demos_weduc,
                                  calfun = "raking")
    rake_demos_weduc <- svydesign(~1, data = pew, 
                                  weights = weights(rake_demos_weduc))
    
    rake_all_vars <- calibrate(design = pew_nowt,
                               formula = formula_rake_all_vars,
                               population = targets_rake_all_vars,
                               calfun = "raking")
    rake_all_vars <- svydesign(~1, data = pew, 
                               weights = weights(rake_all_vars))
    
    #post stratification
    pew_ps <- suppressWarnings(svydesign(ids = ~1, 
                                         data = pew %>% filter(!(strata %in% problem))))
    post_stratification <- postStratify(design = pew_ps,
                                        strata = ~strata,
                                        population = targets_ps)
    post_stratification <- svydesign(~1, data = pew %>% 
                                         filter(!(strata %in% problem)),
                                     weights = weights(post_stratification))
    #reduced
    pew_ps_reduc_xgb <- suppressWarnings(svydesign(ids = ~1, data = pew %>% 
                                                       filter(!(strata_reduc %in% problem_reduc))))
    post_strat_reduc <- postStratify(design = pew_ps_reduc_xgb,
                                     strata = ~strata_reduc,
                                     population = targets_ps_reduc)
    post_strat_reduc <- svydesign(~1, data = pew %>%
                                      filter(!(strata_reduc %in% problem_reduc)),
                                  weights = weights(post_strat_reduc))
    
    pew_ps_all_xgb <- suppressWarnings(svydesign(ids = ~1,  data = pew %>%  
                                                     filter(!(strata_all %in% problem_all))))
    post_strat_all <- postStratify(design = pew_ps_all_xgb,
                                   strata = ~strata_all,
                                   population = targets_ps_all)
    post_strat_all <- svydesign(~1, data = pew %>% 
                                    filter(!(strata_all %in% problem_all)),
                                weights = weights(post_strat_all))
    
    rake_retrospective <- calibrate(design = pew_nowt,
                                    formula = formula_retrospective,
                                    population = targets_retrospective,
                                    calfun = "raking",
                                    force = TRUE)
    rake_retrospective <- svydesign(~1,
                                    data = pew,
                                    weights = weights(rake_retrospective))
    
    # kpop
    kpop <- svydesign(~1, data = pew,
                      weights = kpop_w[kbal_data_sampled ==1])
    
    kpop_conv <- svydesign(~1, data = pew,
                           weights = kpop_conv_w[kbal_data_sampled ==1])
    
    kpop_mf <-  svydesign(~1, data = pew,
                          weights = kpop_mf_w[kbal_data_sampled ==1])
    
    kpop_demos <-  svydesign(~1, data = pew,
                             weights = kpop_demos_w[kbal_data_sampled ==1])
    kpop_demos_wedu <-  svydesign(~1, data = pew,
                                  weights = kpop_demos_wedu_w[kbal_data_sampled ==1])
    
    kpop_all <-  svydesign(~1, data = pew,
                           weights = kpop_all_w[kbal_data_sampled ==1])
    
    
    
    if(is.null(margins_formula)) {
        margins_formula <- ~recode_vote_2016 +
            mod_cces_on_pew_pD + mod_cces_on_pew_pR + mod_cces_on_pew_pO +
            mod_pew_on_pew_pD + mod_pew_on_pew_pR + mod_pew_on_pew_pO + 
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
            mod_cces_on_cces_pD + mod_cces_on_cces_pR + mod_cces_on_cces_pO +
            mod_pew_on_cces_pD + mod_pew_on_cces_pR + mod_pew_on_cces_pO + 
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
    }
    
    margins <- round(cbind(
        pew = svymean(margins_formula, pew_nowt),
        cces_margins = svymean(margins_formula_cces, cces_svy, 
                               na.rm =  TRUE),
        kpop = svymean(margins_formula, kpop),
        kpop_mf = svymean(margins_formula, kpop_mf),
        kpop_conv = svymean(margins_formula, kpop_conv),
        kpop_demos = svymean(margins_formula, kpop_demos),
        kpop_demos_wedu = svymean(margins_formula, kpop_demos_wedu),
        kpop_all = svymean(margins_formula, kpop_all),
        rake_demos_noeduc = svymean(margins_formula, rake_demos_noeduc),
        rake_demos_weduc = svymean(margins_formula, rake_demos_weduc),
        rake_all_vars = svymean(margins_formula, rake_all_vars),
        post_stratification = svymean(margins_formula, post_stratification),
        post_strat_reduc = svymean(margins_formula, post_strat_reduc),
        post_strat_all = svymean(margins_formula, post_strat_all),
        rake_retrospective = svymean(margins_formula, rake_retrospective)) * 100, 5)
    
    rownames(margins) <- sapply(rownames(margins), 
                                function(x) (
                                    if(grepl("recode", x)) {
                                        substr(x, 8,nchar(x) )
                                    } else {
                                        x
                                    }
                                ) )
    #we multiplied the recode_age margins by 100 unnecc above bc it' just a straight mean
    #so let's undo that
    margins["age",] <- margins["age",]/100 
    margins["agesq",] <- margins["agesq",]/100
    margins["agecubed",] <- margins["agecubed",]/100
    
    
    margins_diff <- as.data.frame(margins) %>% 
        mutate(pew = cces_margins - pew,
               kpop = cces_margins - kpop,
               kpop_mf = cces_margins - kpop_mf,
               kpop_conv = cces_margins - kpop_conv,
               kpop_demos = cces_margins - kpop_demos,
               kpop_demos_wedu =cces_margins - kpop_demos_wedu,
               kpop_all = cces_margins - kpop_all,
               rake_demos_noeduc = cces_margins- rake_demos_noeduc,
               rake_demos_weduc = cces_margins- rake_demos_weduc,
               rake_all_vars = cces_margins - rake_all_vars,
               post_stratification = cces_margins-  post_stratification,
               post_strat_reduc = cces_margins-  post_strat_reduc,
               post_strat_all = cces_margins-  post_strat_all,
               rake_retrospective = cces_margins- rake_retrospective)
    rownames(margins_diff) <- rownames(margins)
    
    
    if(pop_weighted_avg) {
        #weighting the difference by the proportion in that bin in cces
        #margins_diff is in percents, so it's the percent difference btwn the sample and the pop on each margin, we want tottake the averge abs value of this differnce, but we dont' just want the straight average, we want it weighted by the percent in that margin in the cces
        #so that means if we just multiply the idfference by the cces margin/100 we only need to sum them in wmabserr
        w_margins_diff <- margins_diff*(margins[, "cces_margins"]/100)
        #dont's square cces:
        w_margins_diff[, "cces_margins"] <- margins_diff[, "cces_margins"]
        #manugally fix the one level margins: then we'll just report the abs diff on these
        w_margins_diff["age",] <- margins_diff["age",]
        w_margins_diff["agesq",] <- margins_diff["agesq",]
        w_margins_diff["agecubed",] <- margins_diff["agecubed",]
        
        
        wmabserr <- w_margins_diff %>% 
            mutate(var = sapply(rownames(margins_diff), 
                                function(y) {
                                    substr(y, 1, 
                                           sapply(rownames(margins_diff), 
                                                  function(x) {
                                                      if(grepl("income",x)) {
                                                          11 #income_5way
                                                          #this is so ridiculous WHY (sunk cost ugh)
                                                      } else if(grepl("xgb",x) |grepl("mod",x) |
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
            group_by(var) %>% summarize(pew_unweighted = sum(abs(pew)), 
                                        kpop = sum(abs(kpop)),
                                        #kpop_conv = sum(abs(kpop_conv)), 
                                        #kpop_mf = sum(abs(kpop_mf)),
                                        kpop_demos = sum(abs(kpop_demos)),
                                        rake_demos_noeduc = sum(abs(rake_demos_noeduc)),
                                        kpop_demos_wedu = sum(abs(kpop_demos_wedu)),
                                        rake_demos_weduc = sum(abs(rake_demos_weduc)),
                                        kpop_all = sum(abs(kpop_all)),
                                        rake_all_vars = sum(abs(rake_all_vars)), 
                                        rake_retrospective = sum(abs(rake_retrospective)),
                                        #post_stratification =
                                        #    sum(abs(post_stratification)),
                                        post_strat_reduc = sum(abs(post_strat_reduc))
                                        #post_strat_all = sum(abs(post_strat_all)),
            )
        wmabserr$var[wmabserr$var == "educ"]<- "educ_6way"
        
    } else {
        #get cces proportions: for weighted L1 dist on age 
        margins_diff[(grep("age_factor", rownames(margins_diff))), ] <-
            (margins_diff[(grep("age_factor", 
                                rownames(margins_diff))),
                          ])* (margins[grep("age_factor", rownames(margins)), "cces_margins"]/100)*length(levels(cces$recode_age_factor))
        
        #HAHAHHAHAHAHA this is horrific I should switch to str_extract now that i belatedly discovered it
        wmabserr <- margins_diff %>% 
            mutate(var = sapply(rownames(margins_diff), 
                                function(y) {
                                    substr(y, 1, 
                                           sapply(rownames(margins_diff), 
                                                  function(x) {
                                                      if(grepl("income",x)) {
                                                          11 #income_5way
                                                          #this is so ridiculous WHY (sunk cost ugh)
                                                      } else if(grepl("xgb",x) |grepl("mod",x) |
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
                                        kpop = mean(abs(kpop)),
                                        #kpop_conv = mean(abs(kpop_conv)), 
                                        #kpop_mf = mean(abs(kpop_mf)),
                                        kpop_demos = mean(abs(kpop_demos)),
                                        rake_demos_noeduc = mean(abs(rake_demos_noeduc)),
                                        kpop_demos_wedu = mean(abs(kpop_demos_wedu)),
                                        rake_demos_weduc = mean(abs(rake_demos_weduc)),
                                        kpop_all = mean(abs(kpop_all)),
                                        rake_all_vars = mean(abs(rake_all_vars)), 
                                        rake_retrospective = mean(abs(rake_retrospective)),
                                        #post_stratification =
                                        #    mean(abs(post_stratification)),
                                        post_strat_reduc = mean(abs(post_strat_reduc))
                                        #post_strat_all = mean(abs(post_strat_all)),
            )
        wmabserr$var[wmabserr$var == "educ"]<- "educ_6way"
    }
    
    
    
    
    #labels
    if(pop_weights) {
        popw_lab = "w/ Pop Weights"
    } else {
        popw_lab = "No pop Weights"
    }
    tol_label = paste0(" tol=",as.character(tolerance))
    maxit_label = paste0( " maxit=", as.character(maxit))
    
    if(pop_weighted_avg) {
        caption_m = paste0("W. Margins (Difference): Unscaled K ", "b=", round(b, 2),
                           popw_lab, " Age Buckets Only (no dnk) ",
                           tol_label, maxit_label)
        caption_o = paste0("Absolute Error on Weighted Margins of Outcomes: Unscaled K ","b=", 
                           round(b, 2),
                           popw_lab, " Age Buckets Only (no dnk) ",
                           tol_label, maxit_label)
        
        caption_err = paste0("Weighted Mean Abs Error: Unscaled K ","b=", 
                             round(b, 2), popw_lab,
                             " Age Buckets Only (no dnk) ", tol_label, maxit_label)   
        
    } else {
        caption_m = paste0("Margins (Difference): Unscaled K ", "b=", round(b, 2),
                           popw_lab, " Age Buckets Only (no dnk) ",
                           tol_label, maxit_label)
        caption_o = paste0("Absolute Error on Margins of Outcomes: Unscaled K ","b=", 
                           round(b, 2),
                           popw_lab, " Age Buckets Only (no dnk) ",
                           tol_label, maxit_label)
        
        caption_err = paste0("Mean Abs Error: Unscaled K ","b=", 
                             round(b, 2), popw_lab,
                             " Age Buckets Only (no dnk) ", tol_label, maxit_label)   
    }
    #var %in% c("agecubed") |
    wmabserr_nooutcomes <- wmabserr %>% 
        filter(!(grepl("xgb", var) | grepl("vote", var))) %>%
        arrange(nchar(var))
    if(is.null(order)) {
        order = c("female", "pid_3way", "age_bucket","race",  "region", "educ_6way",
                  "income_5way", "born",
                  "relig_6way", 
                  "pid_race", "educ_pid",
                  "educ_pid_race", "race_educ_reg", "educ_wh_3way", "midwest_wh_edu",
                  "midwest_edu_race",
                  
                  "agesq", "agecubed", "age_factor")
    }
    
    wmabserr_nooutcomes <- wmabserr_nooutcomes[as.numeric(sapply(order, function(x) which(wmabserr_nooutcomes$var == x))),]
    
    
    
    noedu_aes <- sapply(wmabserr_nooutcomes$var, function(x) {
        ifelse(grepl(x, as.character(formula_rake_demos_noeduc)[2]),
               "gray", "black") })
    wedu_aes <-  sapply(wmabserr_nooutcomes$var, function(x) {
        ifelse(grepl(x, as.character(formula_rake_demos_weduc)[2]) | x == "educ_6way", "gray", "black")
    })
    all_aes <-  sapply(wmabserr_nooutcomes$var, function(x) {
        ifelse(grepl(x, as.character(formula_rake_all_vars)[2]) | x == "educ_6way", "gray", "black")
    })
    
    ps_aes <- sapply(wmabserr_nooutcomes$var, function(x) {
        ifelse(grepl(x, as.character(formula_ps)[2]), "gray", "black") })
    ps_reduc_aes <- sapply(wmabserr_nooutcomes$var, function(x) {
        ifelse(grepl(x, as.character(formula_ps_reduc)[2]), "gray", "black") })
    ps_all_aes <- sapply(wmabserr_nooutcomes$var, function(x) {
        ifelse(grepl(x, as.character(formula_ps_all)[2]), "gray", "black") })
    #interaction picks up race alone so back and fix bc does not include
    retro_aes <- sapply(wmabserr_nooutcomes$var, function(x) {
        ifelse(grepl(x, as.character(formula_retrospective)[2]) & x != "race", "gray", "black") })
    
    
    kpop_demos_aes <- noedu_aes
    kpop_demos_wedu_aes <- wedu_aes
    kpop_all_aes <- all_aes
    
    
    #change name so its clear that recode_educ is 6 way and fixing some labels that are slighly off bc i apparently could not count when i named them/less cleaning in latex names
    
    wmabserr$var[wmabserr$var == "income_5way"]<- "income_6way"
    wmabserr$var[wmabserr$var == "relig_6way"]<- "relig_5way"
    wmabserr$var[wmabserr$var == "educ_wh_3way"]<- "educ_3way_wh"
    wmabserr_nooutcomes$var[wmabserr_nooutcomes$var == "income_5way"]<- "income_6way"
    wmabserr_nooutcomes$var[wmabserr_nooutcomes$var == "relig_6way"]<- "relig_5way"
    wmabserr_nooutcomes$var[wmabserr_nooutcomes$var == "educ_wh_3way"]<- "educ_3way_wh"
    wmabserr$var[wmabserr$var == "race"]<- "race/ethnicity"
    wmabserr_nooutcomes$var[wmabserr_nooutcomes$var == "race"]<- "race/ethnicity"
    
    kable_err <- kable(wmabserr_nooutcomes,
                       format = "latex",
                       caption = caption_err,
                       booktabs = T,
                       digits = 2) %>% kable_styling() %>%
        column_spec(which(colnames(wmabserr_nooutcomes) == "kpop_demos"),
                    color =kpop_demos_aes) %>%
        column_spec(which(colnames(wmabserr_nooutcomes) ==
                              "rake_demos_noeduc"),
                    color = noedu_aes) %>%
        column_spec(which(colnames(wmabserr_nooutcomes) ==
                              "kpop_demos_wedu"),
                    color = kpop_demos_wedu_aes) %>%
        column_spec(which(colnames(wmabserr_nooutcomes) ==
                              "rake_demos_weduc"),
                    color = wedu_aes) %>%
        column_spec(which(colnames(wmabserr_nooutcomes) == "rake_all_vars"),
                    color = all_aes) %>%
        column_spec(which(colnames(wmabserr_nooutcomes) == "kpop_all"),
                    color =kpop_all_aes) %>%
        column_spec(which(colnames(wmabserr_nooutcomes) ==
                              "rake_retrospective"),
                    color =retro_aes) %>% 
        column_spec(which(colnames(wmabserr_nooutcomes) ==
                              "post_strat_reduc"),
                    color =ps_reduc_aes)
    
    kable_om <- kable(wmabserr %>% 
                          filter( grepl("mod", var) | grepl("vote", var)) %>%
                          arrange(nchar(var)),
                      format = "latex",
                      caption = caption_o,
                      booktabs = T,
                      digits = 2)
    
    
    margins_out <- margins_diff[-c(1:9),]
    margins_out <- margins_out[!grepl("age_factor", rownames(margins_out)), ] 
    # margins_out <- margins_out %>% rbind(colMeans(margins_out))
    # rownames(margins_out)[nrow(margins_out)] <- "Column Average"
    
    out <- list(margins = margins,
                w_margins_diff = w_margins_diff,
                margins_diff = margins_diff,
                wmabserr = wmabserr,
                wmabserr_nooutcomes = wmabserr_nooutcomes,
                kable_err = kable_err,
                kable_outmom = kable_om)
    return(out)
}
