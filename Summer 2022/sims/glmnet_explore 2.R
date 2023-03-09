x <- matrix(rnorm(100 * 20), 100, 20)
y <- rnorm(100)

fit <- glmnet(x, y )
l1 <- fit$lambda

# default value; this is exactly teh same as the above of course
penalty.factor <- rep(1, ncol(x))  
fit2 <- glmnet(x, y, penalty.factor = penalty.factor)
l2 <- fit2$lambda

# uniformly scale: this does nothing bc it's all relative?
penalty.factor <- rep(2, ncol(x))  
fit3 <- glmnet(x, y, penalty.factor = penalty.factor)
l3 <- fit3$lambda

# not uniformly scale
penalty.factor <- c(rep(1, ncol(x)/2), rep(2, ncol(x)/2))
fit4 <- glmnet(x, y, penalty.factor = penalty.factor)
l4 <- fit4$lambda

# if 0 then mean all entries are equal
sum(l1!=l2)
sum(l1!=l3)
sum(l1!=l4)
cbind(coef(fit, s = 0.1), coef(fit2, s = 0.1), coef(fit3, s = 0.1), coef(fit4, s = 0.1))


#### with user specified lambda; so this rescales the lambdas
fit <- glmnet(x, y, lambda = lambdas )
l1 <- fit$lambda

# default value
penalty.factor <- rep(1, ncol(x))  
fit2 <- glmnet(x, y, penalty.factor = penalty.factor, lambda = lambdas )
l2 <- fit2$lambda

# uniformly scale
penalty.factor <- rep(2, ncol(x))  
fit3 <- glmnet(x, y, penalty.factor = penalty.factor, lambda = lambdas )
l3 <- fit3$lambda

# not uniformly scale
penalty.factor <- c(rep(1, ncol(x)/2), rep(2, ncol(x)/2))
fit4 <- glmnet(x, y, penalty.factor = penalty.factor, lambda = lambdas)
l4 <- fit4$lambda

#look at
#cbind(l1, l2, l3, l4, lambdas)
sum(l1!=l2)
sum(l1!=l3)
sum(l1!=l4)
cbind(coef(fit, s = 0.1), coef(fit2, s = 0.1), coef(fit3, s = 0.1), coef(fit4, s = 0.1))





####################
lambda.seq <- exp(-seq(0.01, 1, length.out = 30)) / 10

fit1 <- glmnet(x, y, lambda = lambda.seq)
l1 <- fit1$lambda

penalty.factor <- rep(1, ncol(x))
fit2 <- glmnet(x, y, penalty.factor = penalty.factor, lambda = lambda.seq)
l2 <- fit2$lambda

penalty.factor <- rep(2, ncol(x))
fit3 <- glmnet(x, y, penalty.factor = penalty.factor, lambda = lambda.seq)
l3 <- fit3$lambda

penalty.factor <- c(rep(1, ncol(x)/2), rep(2, ncol(x)/2))
fit4 <- glmnet(x, y, penalty.factor = penalty.factor, lambda = lambda.seq)
l4 <- fit4$lambda

#same as above bc it rescales internally so all that matters is the relative difference
penalty.factor <- c(rep(2, ncol(x)/2), rep(4, ncol(x)/2)) 
fit5 <- glmnet(x, y, penalty.factor = penalty.factor, lambda = lambda.seq)
l5 <- fit5$lambda


penalty.factor <- c(rep(0, ncol(x)/2), rep(1, ncol(x)/2))
fit6 <- glmnet(x, y, penalty.factor = penalty.factor, lambda = lambda.seq)
l6 <- fit6$lambda

penalty.factor <- c(rep(0, ncol(x)/2), rep(3, ncol(x)/2))
fit7 <- glmnet(x, y, penalty.factor = penalty.factor, lambda = lambda.seq)
l7 <- fit7$lambda

penalty.factor <- c(rep(1, ncol(x)-1),  0)
fit8 <- glmnet(x, y, penalty.factor = penalty.factor, lambda = lambda.seq)
l8 <- fit8$lambda


fit9 <- glmnet(scale(x), scale(y), lambda = lambda.seq)
l9 <- fit9$lambda

penalty.factor <- c(1, rep(0, ncol(x)))
fit10 <- glmnet(scale(x), scale(y), penalty.factor = penalty.factor, lambda = lambda.seq)





# sum(l1!=l2)
# sum(l1!=l3)
# sum(l1!=l4)
# sum(l1!=l5)
# sum(l1 !=l6)
# sum(l1 !=l7)
cbind(coef(fit1, s = 0.1), 
      #coef(fit2, s = 0.1), 
      #coef(fit3, s = 0.1), 
      #coef(fit4, s = 0.1),
      #coef(fit5, s = 0.1), 
      coef(fit6, s = 0.1),
      #coef(fit7, s = 0.1), 
      #coef(fit8, s = 0.1), 
      coef(fit9, s = 0.1), 
      coef(fit10, s= 0.1))
#so it appearst that actually penalty.factor is for WITHOUT the intercept and the intercept is always penalized when it's added by glmnet cool


#############################
V <- data.frame(kbal_dims = kbal_demos_est$svdK$v[, 1:kbal_demos_est$numdims])
#binding constraint cols to
X <- as.matrix(cbind(kbal_demos_est$appended_constraint_cols[kbal_data_sampled==1, ], V))

colnames(X)
class(X)

fit = cv.glmnet(X, kpop_demos_svyd$variables$outcome, alpha = 0, lambda = lambdas)
l1 = fit$lambda

#rep the default

penalty = rep(1, rep(0, ncol(X))-1)
fit2 = cv.glmnet(X, kpop_demos_svyd$variables$outcome, alpha = 0, lambda = lambdas,
                 penalty.factor = penalty)
l2 = fit2$lambda

#now force some to be zero: intercept now not regularized...
penalty = c(rep(0, ncol(kbal_demos_est$appended_constraint_cols)),
            rep(1, 1+kbal_demos_est$numdims))
fit3 = glmnet(X, kpop_demos_svyd$variables$outcome, alpha = 0, lambda = lambdas, 
                     penalty.factor = penalty)
l3 = fit3$lambda


penalty = c(rep(0, 1+ncol(kbal_demos_est$appended_constraint_cols)),
            rep(1, kbal_demos_est$numdims))
fit4 = glmnet(X, kpop_demos_svyd$variables$outcome, alpha = 0, lambda = lambdas, 
              penalty.factor = penalty)
l4 = fit4$lambda


penalty = c(1, rep(0, ncol(kbal_demos_est$appended_constraint_cols)),
            rep(1, kbal_demos_est$numdims))
fit5 = glmnet(X, kpop_demos_svyd$variables$outcome, alpha = 0, lambda = lambdas, 
              penalty.factor = penalty)
l5 = fit4$lambda



###
sum(lambdas != l1)
sum(l1 != l2)
sum(l1 != l3)

#allll that actualyl amtters is the number of 1's bc these penalities are being rescaled internally, so 2-5 all have the same number of 1's even in a different order...
round(cbind(coef(fit, s = 0.1), coef(fit2, s = 0.1), coef(fit3, s = 0.1), coef(fit4, s = 0.1),  coef(fit5, s = 0.1)),3)



###### now check cv.glmenet with penalty works the same way as glmnet with peanlty
penalty = c(rep(0, ncol(kbal_demos_est$appended_constraint_cols)),
            rep(1, kbal_demos_est$numdims))
fit5 = glmnet(X, kpop_demos_svyd$variables$outcome, alpha = 0, lambda = lambdas, 
              penalty.factor = penalty)
cv_fit5 =  cv.glmnet(X, kpop_demos_svyd$variables$outcome, alpha = 0, lambda = lambdas, 
                  penalty.factor = penalty)
#yep!
cbind(coef(fit5, s = 0.1), coef(cv_fit5$glmnet.fit, s = 0.1))

l5 = fit4$lambda



####### does the rescaling of lambda matter in terms of predictions/performance or is it all a relative scaling and who fucking cares
penalty = c(rep(0, ncol(kbal_demos_est$appended_constraint_cols)),
            rep(1, kbal_demos_est$numdims))
fit6 = glmnet(X, kpop_demos_svyd$variables$outcome, alpha = 0,
              penalty.factor = penalty)
#def diff lambdas now
sum(fit6$lambda != fit5$lambda)
#with forced lambda its the same, duh ok good
cbind(coef(fit5, s = 0.1), coef(fit6, s = 0.1))

#now with "best lambda"
#small simulation
j = 0
reps = 500
out = rep(NA, reps)
for(i in 1:reps) {
    seed = round(abs(10000000*rnorm(1)))
    set.seed(seed)
    cv_fit5 =  cv.glmnet(X, kpop_demos_svyd$variables$outcome, alpha = 0, lambda = lambdas, 
                         penalty.factor = penalty)
    set.seed(seed)
    cv_fit6 =  cv.glmnet(X, kpop_demos_svyd$variables$outcome, alpha = 0,
                         penalty.factor = penalty)
    #obviously these lambdas are different?
    cv_fit5$lambda.min
    cv_fit6$lambda.min
    #ok so duh so are the predictions, but which is "better"? they are very close
    cbind(coef(cv_fit5$glmnet.fit, s = cv_fit5$lambda.min), coef(cv_fit6$glmnet.fit, s = cv_fit6$lambda.min))
    #round(cbind(cv_fit5$lambda, cv_fit6$lambda),4)
    #cv error is? ... the same... so yeah i dont think it really matters?
    cv_fit5$cvm[which(cv_fit5$lambda == cv_fit5$lambda.min)]
    cv_fit6$cvm[which(cv_fit6$lambda == cv_fit6$lambda.min)]
    j = j + as.numeric(cv_fit5$cvm[which(cv_fit5$lambda == cv_fit5$lambda.min)] > cv_fit6$cvm[which(cv_fit6$lambda == cv_fit6$lambda.min)])
    out[i] = cv_fit5$cvm[which(cv_fit5$lambda == cv_fit5$lambda.min)] > cv_fit6$cvm[which(cv_fit6$lambda == cv_fit6$lambda.min)]
}
#percent of the time passed in lambdas have lower cv error than internally selected is
#NOTE that this is like way iout in the 7th decimal or something so it really doesnt matter
j/reps
sum(out)/length(out)
#k so it seems fucking fine whatev er
