#check comparison to svy package



est_mean <- function(outcome, design) {
    svymean(as.formula(paste0("~", outcome)), design, na.rm = TRUE)[1]
}

est_mean("diff_cces_on_pew", kpop_svyd)
svymean(~diff_cces_on_pew, kpop_svyd)[1]
kpop_se = out$SEs["kpop",]
kpop_se



#first do the raking SEs look right for all, no fixed are off
rake_demos_weduc_se
rake_demos_weduc_svyd = rake_demos_weduc
est_mean("diff_cces_on_pew", rake_demos_weduc_svyd)
#weights that sum to 1 work for everythign but not fixed
sum(weights(rake_demos_weduc_svyd))
sum(weights(kpop_svyd))
#test
var_fixed <- function(Y, weights, pop_size) {
    ## note: needs weights that sum to population total
    cat("weights sum to", sum(weights), "\n")
    
    if(round(sum(weights)) == pop_size) { 
        return(Hmisc::wtd.var(Y, weights))
        
    } else if(round(sum(weights)) == 1) {
        warning("weights do not sum to pop_size but to 1, adjusting accordingly", immediate. = T)
        return(Hmisc::wtd.var(Y, weights* pop_size))
        
    } else {
        return("Oops")
    }
    
    # Hmisc::wtd.var(Y, weights)
    # Hmisc::wtd.var(Y, test)
    # Hmisc::wtd.var(Y, test, normwt = T)
    
    # weights = weights(rake_demos_weduc_svyd)
    #ok so wtd.var is internal scaling?
    # sw = sum(weights)
    # x = Y
    # xbar <- sum(weights * x)/sw
    # sum(weights * ((x - xbar)^2))/(sw - 1)
    # Hmisc::wtd.var(Y, weights)
    # Hmisc::wtd.var(Y, weights*pop_size)
    # sum(weights)
    # 
    #return(Hmisc::wtd.var(Y, weights * pop_size))
}


sum(weights(rake_demos_weduc_svyd))
var_fixed(Y = rake_demos_weduc_svyd$variables$diff_cces_on_pew, 
          weights = weights(rake_demos_weduc_svyd),
          pop_size = nrow(cces))

var_fixed(Y = rake_demos_weduc_svyd$variables$diff_cces_on_pew, 
          weights = weights(rake_demos_weduc_svyd)*nrow(cces),
          pop_size = nrow(cces))

out_temp = var_fixed(Y = rake_demos_weduc_svyd$variables$diff_cces_on_pew, 
                weights = weights(rake_demos_weduc_svyd)*nrow(cces),
                pop_size = nrow(cces))
out_temp
sqrt(out_temp)

sqrt(out_temp)/sqrt(nrow(cces))
sqrt(out_temp)/sqrt(length(Y))
rake_demos_weduc_se
length(weights)

#ok so these are the same regardless of teh scaling of the weights, how..
weights = weights(rake_demos_weduc_svyd)*nrow(cces)
sw = sum(weights)
sw
x = Y
xbar <- sum(weights * x)/sw
xbar
#same as estimate
est_mean("diff_cces_on_pew", rake_demos_weduc_svyd)
#var
sum(weights * ((x - xbar)^2))/(sw - 1)

#var is the same weither you use weights that sum to N_pop or to 1 it doesnt matter its just bc the denom sw-1 fucks up otherwise
#so then what 1/sqrt(n) are you supposed to use



#### weighted variance of weights that sum to 1 sqrted
sum(weights(rake_demos_weduc_svyd))
#this should sort of break
Hmisc::wtd.var(Y, weights(rake_demos_weduc_svyd))
#but we can fix it by internally nromalizing to length(Y)
#why are these not identical?
Hmisc::wtd.var(Y, as.numeric(weights(rake_demos_weduc_svyd)), normwt=T)
Hmisc::wtd.var(Y, as.numeric(weights(rake_demos_weduc_svyd))*length(Y))
#bc normwt uses this:
as.numeric(stats::cov.wt(cbind(Y), as.numeric(weights(rake_demos_weduc_svyd))*length(Y))$cov)
#internally it's turtles all the way down so lets call this a bessel correction and say its the same thign

stupid = as.numeric(weights(rake_demos_weduc_svyd))
(stupid* length(Y)/sum(stupid))[1:10]
(stupid*length(Y))[1:10]


#### weighted var of weights that sum to N_pop sqrt/N_pop
sum(weights(rake_demos_weduc_svyd)*nrow(cces))
sqrt(Hmisc::wtd.var(Y, weights(rake_demos_weduc_svyd)*nrow(cces))/nrow(cces) )





#confirm taht fixed weights + weights that sum to pop size /sqrt(n)


### testing
residuals = residuals(lm(update(formula_rake_demos_weduc, diff_cces_on_pew ~ .), 
                         data = rake_demos_weduc_svyd$variables))
rake_demos_weduc_se_2 <- calc_SEs(Y = rake_demos_weduc_svyd$variables$diff_cces_on_pew, 
                                residuals = residuals, 
                                pop_size = nrow(cces), 
                                sample_size = nrow(pew),
                                weights = weights(rake_demos_weduc_svyd))
names(rake_demos_weduc_se_2) = paste0("rake_demos_weduc_", names(rake_demos_weduc_se_2))
rake_demos_weduc_se*100
rake_demos_weduc_se_2*100
old*100
old <- calc_SEs(Y = rake_demos_weduc$variables$diff_cces_on_pew, 
                                residuals = residuals(lm(update(formula_rake_demos_weduc,
                                                                diff_cces_on_pew ~ .), 
                                                         data = rake_demos_weduc$variables)), 
                                pop_size = nrow(cces), 
                                sample_size = nrow(pew),
                                weights = weights(rake_demos_weduc))




#kpop testing
lambdas <- 10^seq(3, -2, by = -.1)
load(paste0(path_weights, "DEBUG_full_kbal_obj_2023-04-06.Rdata"))
x <- model.matrix( ~ ., data = as.data.frame(kbal_est$svdK$v[, 1:kbal_est$numdims]))
fit <- glmnet(x, 
              kpop_svyd$variables$diff_cces_on_pew, alpha = 0, lambda = lambdas)
fit <- glmnet(x, 
              kpop_svyd$variables$diff_cces_on_pew, alpha = 0, lambda = NULL)
cv_fit <- cv.glmnet(x, kpop_svyd$variables$diff_cces_on_pew, alpha = 0, lambda = lambdas)
opt_lambda <- cv_fit$lambda.1se
fit <- cv_fit$glmnet.fit

residuals = kpop_svyd$variables$diff_cces_on_pew - predict(fit, s = opt_lambda, newx = x)

# lm_full = lm(kpop_svyd$variables$diff_cces_on_cces ~ kbal_est$svdK$v[, 1:kbal_est$numdims])
# full_res = residuals(lm_full)
# cbind(residuals, full_res)
kpop_se <- tryCatch(calc_SEs(Y = kpop_svyd$variables$diff_cces_on_pew,
                             residuals = residuals,
                             pop_size = nrow(cces),
                             sample_size = nrow(pew),
                             weights = weights(kpop_svyd)), error = function(e) NA)
kpop_se
out$SEs["kpop",]


rake_demos_weduc_se


residuals = residuals(lm(update(formula_rake_demos_weduc, diff_cces_on_cces ~ .), 
                        data = rake_demos_weduc_svyd$variables))
test <- tryCatch(calc_SEs(Y = rake_demos_weduc_svyd$variables$diff_cces_on_cces,
                             residuals = residuals,
                             pop_size = nrow(cces),
                             weights = weights(rake_demos_weduc_svyd)), error = function(e) NA)
test
weights = as.numeric(weights(rake_demos_weduc_svyd))
sum(weights)
pop_size = 1
sum(weights < 0)

sum((weights^2 - weights / pop_size) * residuals^2)
test
weights = as.numeric(weights(rake_demos_weduc_svyd))*nrow(cces)
pop_size = nrow(cces)
sum((weights^2 - weights / pop_size) * residuals^2)
weights = as.numeric(weights(rake_demos_weduc_svyd))
var_quasi(weights, residuals, nrow(cces))
sqrt(var_quasi(weights, residuals, nrow(cces)))
var_quasi(weights*nrow(cces), residuals, 1)

###################################################
#ran local one shot with new selection model to se SEs and R^2
SEs = data.frame(rake_demos_noeduc_se,
                 rake_demos_weduc_se,
                 rake_all_se,
                 post_stratification_se,
                 post_strat_reduc_se,
                 post_strat_all_se,
                 rake_truth_se)#,
                 #kpop_se)

SE_fixed = SEs[grepl("SE_fixed$", colnames(SEs))]
SE_linear = SEs[grepl("SE_linear$", colnames(SEs))]
SE_quasi = SEs[grepl("SE_quasi$", colnames(SEs))]
SE_chad= SEs[grepl("SE_chad$", colnames(SEs))]


################################
dat = rake_demos_noeduc_svyd$variables


local_outcome = rbinom(nrow(cces), 1, cces$mod_cces_on_cces_pR)
local_outcome = cces$diff_cces_on_cces
#outcome = plogis(var(cces$mod_cces_on_cces_pR) + cces$mod_cces_on_cces_pR)
summary(outcome)

#outcome = cces$diff_cces_on_cces
rake_demos_noeduc_lm = summary(lm(update(formula_rake_demos_noeduc, outcome ~ .), 
                         data = rake_demos_noeduc_svyd$variables))
rake_demos_noeduc_lm$adj.r.squared

rake_demos_weduc_lm = summary(lm(update(formula_rake_demos_weduc, outcome ~ .), 
                         data = rake_demos_weduc_svyd$variables))
rake_all_lm = summary(lm(update(formula_rake_all_vars, outcome ~ .), 
                         data = rake_all_svyd$variables))

R2 = data.frame((rake_demos_noeduc_lm)$adj.r.squared,
                (rake_demos_weduc_lm)$adj.r.squared,
                (rake_all_lm)$adj.r.squared)

R2
outcome = cces$outcome
rake_demos_noeduc_se <- calc_SEs(Y = rake_demos_noeduc_svyd$variables$outcome, 
                                 residuals = residuals(rake_demos_noeduc_lm), 
                                 pop_size = nrow(cces), 
                                 sample_size = sum(sample),
                                 weights = weights(rake_demos_noeduc_svyd))
rake_demos_weduc_se <- calc_SEs(Y = rake_demos_noeduc_svyd$variables$outcome, 
                                residuals = residuals(rake_demos_weduc_lm), 
                                pop_size = nrow(cces),
                                sample_size = sum(sample),
                                weights = weights(rake_demos_weduc_svyd))
rake_all_se <- calc_SEs(Y = rake_demos_noeduc_svyd$variables$outcome,
                        residuals = residuals(rake_all_lm), 
                        pop_size = nrow(cces), 
                        sample_size = sum(sample),
                        weights = weights(rake_all_svyd))
SEs = data.frame(rake_demos_noeduc_se,
                 rake_demos_weduc_se,
                 rake_all_se)



out = matrix(c(SEs[, grep("fixed", colnames(SEs))]/SEs[, grep("chad", colnames(SEs))], 
               R2), 3, 2)
rownames(out) = c("r_demos", "r_demos_edu", "r_all")
colnames(out) = c("ratio", "Adj_R2")
out = out[c(8,7,1:6),]
out = data.frame(out)
out[,1] = as.numeric(out[,1])
out[,2] = as.numeric(out[,2])
out %>% arrange(desc(Adj_R2))
































R2 = data.frame(summary(rake_demos_noeduc_lm)$adj.r.squared,
                summary(rake_demos_weduc_lm)$adj.r.squared,
                summary(rake_all_lm)$adj.r.squared,
                summary(post_stratification_lm)$adj.r.squared,
                summary(post_strat_reduc_lm)$adj.r.squared,
                summary(post_strat_all_lm)$adj.r.squared,
                summary(rake_truth_lm)$adj.r.squared,
                summary(kpop_lm)$adj.r.squared)

SEs
R2
out = matrix(c(SEs[, grep("fixed", colnames(SEs))]/SEs[, grep("chad", colnames(SEs))], 
      R2), 8, 2)
rownames(out) = c("r_demos", "r_demos_edu", "r_all", "ps_orig", "ps_reduc", "ps_all", "r_truth", "kpop")
colnames(out) = c("ratio", "Adj_R2")
out = out[c(8,7,1:6),]
out = data.frame(out)
out[,1] = as.numeric(out[,1])
out[,2] = as.numeric(out[,2])
out %>% arrange(desc(Adj_R2))


comp = matrix(c(SE_fixed, SE_chad), 8,2)
rownames(comp) = c("r_demos", "r_demos_edu", "r_all", "ps_orig", "ps_reduc", "ps_all", "r_truth", "kpop")
colnames(comp) = c("fixed", "chad")
comp = comp[c(8,7,1:6),]
comp = as.data.frame(comp)
comp[,1] = as.numeric(comp[,1])
comp[,2] = as.numeric(comp[,2])
round(comp,4)

comp = cbind(comp, out)
comp %>% arrange(desc(Adj_R2)) %>% round(., 4)
