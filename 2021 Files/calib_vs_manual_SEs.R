##### Data
path_data= "/Users/Ciara_1/Dropbox/kpop/Updated/application/data/"
cces <- readRDS(paste0(path_data, "cces_lasso_061021.rds"))
pew <- readRDS(paste0(path_data, "pew_lasso_061021.rds"))

#cces_awt <- svydesign(~1, weights = ~commonweight_vv_post, data = cces)
#Not using pop weights 
library(survey)
library(kbal)
library(dplyr)
pew_nowt <- suppressWarnings(svydesign(~1, data = pew))
cces_nowt <- svydesign(ids = ~1, data = cces)



###### Calibrate Object way:
formula_rake_all_vars <- ~recode_age_bucket + recode_female + 
    recode_race + recode_region + recode_pid_3way + recode_educ +
    recode_income_5way + recode_relig_6way + recode_born + recode_attndch_4way


create_targets <- function(target_design, target_formula) {
    target_mf <- model.frame(target_formula, model.frame(target_design))
    target_mm <- model.matrix(target_formula, target_mf)
    wts <- weights(target_design)
    colSums(target_mm * wts) / sum(wts)
}

targets_rake_all_vars <- create_targets(cces_nowt, 
                                        formula_rake_all_vars)

rake_all_vars <- calibrate(design = pew_nowt,
                           formula = formula_rake_all_vars,
                           population = targets_rake_all_vars,
                           calfun = "raking")



######### Ebal Way:
X_all <- bind_rows(pew, cces) %>%
    model.matrix(formula_rake_all_vars, .)
X_all <- X_all[,-1]

target_ebal = c(rep(0, nrow(pew)), rep(1, nrow(cces)) )

ebal_all_vars <- ebalance_custom(Treatment = target_ebal,
                                 X = X_all,
                                 base.weight = NULL,
                                 norm.constant  = NULL,
                                 coefs = NULL ,
                                 max.iterations = 200,
                                 constraint.tolerance = 1e-6,
                                 print.level=0)

########### Compare weights: should be identical
calib_w = weights(rake_all_vars) #weirdly don't see any way to get weights directly from calibrate obj
ebal_w = ebal_all_vars$w/sum(ebal_all_vars$w)
cbind(calib_w[1:5], ebal_w[1:5])
sum(round(ebal_w - calib_w, 5))


########## Compare Svymean SEs:
vote_diff <- quote((recode_vote_2016Democrat - recode_vote_2016Republican) /(recode_vote_2016Democrat + recode_vote_2016Republican + recode_vote_2016Other)) 

calibrate_res  = svycontrast(svymean(~recode_vote_2016, 
                                     rake_all_vars, na.rm = TRUE),
                             vote_diff)

ebal_svy =  svydesign(~1, data = pew,
                      weights = ebal_w)
ebal_res =  svycontrast(svymean(~recode_vote_2016, 
                                ebal_svy, na.rm = TRUE),
                        vote_diff)

#w calibrate object we get smaller SEs?
calibrate_res
ebal_res

#what's actually different between these two: calibrate() provides a $postStrata obj
#if we import that over we see we get the same reduced SEs now
rake_all_vars$postStrata
test_ebal = ebal_svy
test_ebal$postStrata = rake_all_vars$postStrata
svycontrast(svymean(~recode_vote_2016, 
                    test_ebal, na.rm = TRUE),
            vote_diff)

#so what is actually driving this difference? turns out changing stage arg to 1 we get back to the ebal_svy
#what is this arg? according to doc:
#"In a model with two-stage sampling, population totals may be available for the PSUs actually sampled, but not for the whole population. In this situation, calibrating within each PSU reduces with second-stage contribution to variance. This generalizes to multistage sampling. The stage argument specifies which stage of sampling the totals refer to. Stage 0 is full population totals, stage 1 is totals for PSUs, and so on. The default, stage=NULL is interpreted as stage 0 when a single population vector is supplied and stage 1 when a list is supplied. Calibrating to PSU totals will fail (with a message about an exactly singular matrix) for PSUs that have fewer observations than the number of calibration variables."

#ebal way is using the PSU totals, which makes sense as it was never told any information about the population
test_ebal = ebal_svy
test_ebal$postStrata = rake_all_vars$postStrata
test_ebal$postStrata[[1]]$stage = 1
svycontrast(svymean(~recode_vote_2016, 
                    test_ebal, na.rm = TRUE),
            vote_diff)
ebal_res

#calibrate using full pop totals:
test_ebal = ebal_svy
test_ebal$postStrata = rake_all_vars$postStrata
test_ebal$postStrata[[1]]$stage = 0
svycontrast(svymean(~recode_vote_2016, 
                    test_ebal, na.rm = TRUE),
            vote_diff)
calibrate_res


ebal_svy$fpc$sampsize

#difference in SEs starts w survey mean
svymean(~recode_vote_2016, 
        test_ebal, na.rm = TRUE)
svymean(~recode_vote_2016, 
        rake_all_vars, na.rm = TRUE)

rake_all_vars$ps
#digging into the package to see how it's actually coded differently
svyrecvar
#the key difference:
#stage 0 = "full population totals" matches ebal
if (psvar$stage == 0) {
    x <- as.matrix(qr.resid(psvar$qr, x/psvar$w) * 
                       psvar$w)
#stage 1:
} else {
    cal <- c(cal, list(psvar))
}

