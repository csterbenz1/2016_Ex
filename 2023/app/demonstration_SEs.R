#demonstration of svy package SE's divergence

#### SE Calc Functions:
########################################
############## variance calc ###########
## Variance functions
var_fixed <- function(Y, weights, pop_size) {
    ## note: needs weights that sum to population total
    if(round(sum(weights)) != pop_size) { weights = weights*pop_size/sum(weights)}
    return(Hmisc::wtd.var(Y, weights))
}

## kott (14) (under poisson)
var_quasi <- function(weights, residuals, pop_size) {
    #moving from kott 14 w sum w =N to weights that sum to 1 + var of total to var of mean:
    #sum^n (w_i^2 - 1/N_pop w_i)e_i^2 
    return(sum((weights^2 - (weights / pop_size))*residuals^2))
}

## kott (15) linearization
var_linear <- function(weights, residuals, sample_size) {
    #moving from kott 14 w sum w =N to weights that sum to 1 + var of total to var of mean:
    # n/(n-1) sum^n (w_i*e_i)^2 - (1/n-1) [sum^n] *using this for now
    # OR approx = sum^n (w_i*e_i)^2 - (1/n) [sum^n]
    n = sample_size
    return((n/(n-1))*sum((weights * residuals)^2) - 1/(n-1) * sum(weights * residuals)^2)
}

## chad
var_chad <- function(weights, residuals) {
    return(sum(weights^2 * residuals^2))
}

## calculate all variances
calc_SEs <- function(Y, residuals, pop_size, weights, sample_size) {
    if(round(sum(weights)) != 1 ) {
        weights = weights/sum(weights)
    }
    return(data.frame(SE_fixed = sqrt(var_fixed(Y, weights, pop_size) / length(Y)),
                      SE_quasi = sqrt(var_quasi(weights, residuals, pop_size)),
                      SE_linear = sqrt(var_linear(weights, residuals, sample_size)),
                      SE_chad = sqrt(var_chad(weights, residuals))))
}


create_targets <- function (target_design, target_formula) {
    target_mf <- model.frame(target_formula, model.frame(target_design))
    target_mm <- model.matrix(target_formula, target_mf)
    wts <- weights(target_design)
    colSums(target_mm * wts) / sum(wts)
}

########## Load Data:
path_data= "/Users/Ciara_1/Dropbox/kpop/Updated/application/data/" 
pew <- readRDS(paste0(path_data, "pew_lasso_061021.rds"))
cces <- readRDS(paste0(path_data, "cces_lasso_061021.rds"))

##################### Start with Rake + Demos + No Edu Ex:

formula_rake_demos_noeduc <- ~recode_age_bucket + recode_female + recode_race +
    recode_region + recode_pid_3way
pop_weights = T
if(pop_weights) {
    cces_svy <- svydesign(ids = ~1, weights = ~commonweight_vv_post,
                          data = cces)
} else {
    cces_svy <- suppressWarnings(svydesign(ids = ~1,
                                           data = cces))
}
targets_rake_demos_noeduc <- create_targets(cces_svy,
                                            formula_rake_demos_noeduc)

pew_nowt <- suppressWarnings(svydesign(ids = ~1, data = pew))

#note that unweighted SEs ~1.5%
svymean(~diff_cces_on_pew, pew_nowt)
#even probabilties assumed, sums to nrow(pew)
sum(pew_nowt$prob)
#argument not discussed in any documentation but I think has to do with multistage sampling
#unweighted sums to nrow(pew)
sum(pew_nowt$allprob)

#using calibrate
rake_calibrate <- calibrate(design = pew_nowt,
                                    formula = formula_rake_demos_noeduc,
                                    population = targets_rake_demos_noeduc,
                                    calfun = "raking")
#resulting object is a svydesign object
rake_calibrate
#now SE is on the order of ~.5%
svymean(~diff_cces_on_pew, rake_calibrate)

#notably this matches the linearization SEs we calculate by hand:
residuals = residuals(lm(update(formula_rake_demos_noeduc, diff_cces_on_pew ~ .), 
                         data = rake_calibrate$variables))
rake_demos_noeduc_se <- calc_SEs(Y = rake_calibrate$variables$diff_cces_on_pew, 
                                 residuals = residuals, 
                                 pop_size = nrow(cces), 
                                 sample_size =nrow(pew),
                                 weights = weights(rake_calibrate))

#we see that our linear SEs exactly match the svymean SE's on the calibrate svydesign
rake_demos_noeduc_se$SE_linear
data.frame(svymean(~diff_cces_on_pew, rake_calibrate))[1,2]

#BUT: if we make a svydesign object directly with the weights from the calibration:
rake_incorr <- svydesign(~1, data = pew, 
                         weights = weights(rake_calibrate))

#we get totally different SEs back to being on the order of ~1.5%
svymean(~diff_cces_on_pew, rake_incorr)

#so what happened?? what is different
#all the same
sum(rake_incorr$cluster != rake_calibrate$cluster)
sum(rake_incorr$strata != rake_calibrate$strata)
rake_incorr$has.strata == rake_calibrate$has.strata
sum(rake_incorr$prob != rake_calibrate$prob)
rake_incorr$fpc
rake_calibrate$fpc
rake_incorr$pps ==rake_calibrate$pps

#only: here is a difference; I can't find any documentation about this parameter but I believe it's involved in multistage sampling perhaps? 
rake_incorr$allprob
rake_calibrate$allprob
#prob and allprob are the same in the svydesign w/weights obj
rake_incorr$prob == rake_incorr$allprob
#but not in the calibrate svydesign obj
rake_calibrate$prob == rake_calibrate$allprob


#is this the issue? can we hack our way from one SE est to the other?
test <- svydesign(ids = ~1, data = pew, 
                  weights = weights(rake_calibrate))
#as before ~1.5%
svymean(~diff_cces_on_pew, test)
#force allprob here to match the calibrate prob to try and hack our way back
test$allprob = rake_calibrate$allprob
#but makes no difference... so I'm at a loss for what's going on,  sometdhing internal that's not an outputted parameter?
svymean(~diff_cces_on_pew, test)
svymean(~diff_cces_on_pew, rake_calibrate)


#######################################################
#ok so this now makes sense w Erin's input:
#svymean() on calibrate object uses residualization and linear SEs
#svymaen() on svydesign object dn know the vars used to create the weights so cn residualize and uses Fixed weight SEs



#one weird thing left is taht unweighted SEs are larger than weighted:
#1) unweighed is not residualized so that is on the order of fixed weight SEs?

#yep on the same order
svymean(~diff_cces_on_pew, pew_nowt)
rake_demos_noeduc_se$SE_fixed

#and then ofc the weighted pew should then be on the same as the linearization SEs as we saw abvoe
svymean(~diff_cces_on_pew, rake_calibrate)
rake_demos_noeduc_se$SE_linear



#2)unweighted + residualization
formula_null = ~1

create_targets <- function (target_design, target_formula) {
    target_mf <- model.frame(target_formula, model.frame(target_design))
    target_mm <- model.matrix(target_formula, target_mf)
    wts <- weights(target_design)
    colSums(target_mm * wts) / sum(wts)
}
targets_null <- create_targets(cces_svy,formula_null)

rake_null<- calibrate(design = pew_nowt,
                            formula = formula_null,
                            population = targets_null,
                            calfun = "raking")
#no change
svymean(~diff_cces_on_pew, rake_null)
svymean(~diff_cces_on_pew, pew_nowt)

#manually:
residuals = residuals(lm(update(formula_null, diff_cces_on_pew ~ .), 
                         data = rake_null$variables))

rake_null_se <- calc_SEs(Y = rake_null$variables$diff_cces_on_pew, 
                                 residuals = residuals, 
                                 pop_size = nrow(cces), 
                                 sample_size =nrow(pew),
                                 weights = weights(rake_null))
#yep so now linearization w just the intercept matches the svymean() SE
rake_null_se$SE_linear
data.frame(svymean(~diff_cces_on_pew, pew_nowt))[1,2]


#################################################################################
# but why are the kpop Fixed SEs so much smaller than the svy package ones?
# does svymean() on svydesign obj match teh fixed SEs or just same order of mag?
svymean(~diff_cces_on_pew, rake_incorr)
rake_demos_noeduc_se$SE_fixed

#ok no so the fixed SE are in general slightly smaller even w the raking, so? is that a bessel corr or pop size thing?
#FPC = sqrt(N-n/N-1): sig = sig/sqrt(n) *(sqrt(N-n/N-1))

#fixed weights are:
weights = weights(rake_demos_weduc_svyd)*nrow(cces)
sw = sum(weights)
#sum to pop N
sw
Y = rake_demos_weduc_svyd$variables$diff_cces_on_pew
x = Y
xbar <- sum(weights * x)/sw
xbar
#same as estimate
data.frame(svymean(~diff_cces_on_pew, rake_demos_weduc_svyd))[1,2]
#var
sum(weights * ((x - xbar)^2))/(sw - 1)
#this is the fixed w est, but it is not hte svymean() on svydesign weights
sqrt(sum(weights * ((x - xbar)^2))/(sw - 1)/(nrow(pew)))
rake_demos_noeduc_se$SE_fixed

#now let's try and adjust the kpop fixed we SEs and se iff they match the svy package:
kpop <- svydesign(~1, data = pew, weights = out$weights$kpop_w)
data.frame(svymean(~diff_cces_on_pew, kpop))[1,2]
kpop_fixed_w = out$SEs["kpop","SE_fixed"]
kpop_fixed_w

fpc_pass = sqrt(nrow(cces) -nrow(pew))/sqrt(nrow(cces) -1)
fpc_pass*kpop_fixed_w

dim(kpop$fpc$sampsize)
kpop_2 <- svydesign(~1, data = pew, fpc = matrix(fpc_pass, nrow = 2052, ncol = 1),
                  weights = out$weights$kpop_w)
kpop_2$fpc$popsize
kpop_2$fpc$sampsize
data.frame(svymean(~diff_cces_on_pew, kpop_2))[1,2]

#so to get the true correct fpc we can see it needs a matrix 2052 x1 and the pop size will be = 2052 / x 
kpop_2 <- svydesign(~1, data = pew, fpc = matrix(nrow(cces), nrow = 2052, ncol = 1),
                    weights = out$weights$kpop_w)
kpop_2$fpc$popsize
kpop_2$fpc$sampsize
data.frame(svymean(~diff_cces_on_pew, kpop_2))[1,2]
kpop_fixed_w


#ok so now can we actually try and match this correct fpc corrected Se number manually:
data.frame(svymean(~diff_cces_on_pew, kpop_2))[1,2]
weights = out$weights$kpop_w*nrow(cces)/nrow(pew)
sw = sum(weights)
#sum to pop N
sw
Y = kpop$variables$diff_cces_on_pew
x = Y
xbar <- sum(weights * x)/sw
xbar
#yep we are replicating the fixed w as is:
sqrt(sum(weights * ((x - xbar)^2))/(sw - 1)/(nrow(pew)))
out$SEs["kpop",]
kpop_fixed_w
#now w fpc:
#but weight yeah this makes the SE SMALLER... so no way its' gonna be the same wtf is svy doign
sqrt(sum(weights * ((x - xbar)^2))/(sw - 1)/(nrow(pew)))*fpc_pass
data.frame(svymean(~diff_cces_on_pew, kpop_2))[1,2]

#fpc saves us like nothing lol ofc bc pew is a tiny percenty of cces
data.frame(svymean(~diff_cces_on_pew, kpop))[1,2] - data.frame(svymean(~diff_cces_on_pew, kpop_2))[1,2]


#### wtf is going on: ok so it is usign teh HT variance calc....
kpop <- svydesign(~1, data = pew, weights = out$weights$kpop_w, variance = "HT")
data.frame(svymean(~diff_cces_on_pew, kpop))[1,2]


data.frame(svymean(~diff_cces_on_pew, kpop, deff = F))

#### from digging thorugh the opackage i can get very close to the same SEs:
#normalize weights to sum to 1
norm_w = kpop_w/sum(kpop_w)
sum(norm_w)
#woutcome vector - weighted mean (est)
y =  matrix(pew$diff_cces_on_pew)
est = sum(matrix(pew$diff_cces_on_pew)*norm_w)
est
w_y = (y - est)*norm_w
var = crossprod(w_y)
#same as:
t(w_y) %*% w_y
var
#very very close, probably som fpc correction here?
sqrt(var)
data.frame(svymean(~diff_cces_on_pew, kpop))

#replicate with unweighted
norm_w = rep(1, nrow(pew))/nrow(pew)
sum(norm_w)
#woutcome vector - weighted mean (est)
y =  matrix(pew$diff_cces_on_pew)
est = sum(matrix(pew$diff_cces_on_pew)*norm_w)
est
w_y = (y - est)*norm_w
var = crossprod(w_y)
#same as:
t(w_y) %*% w_y
var
#very very close, probably som fpc correction here?
sqrt(var)
data.frame(svymean(~diff_cces_on_pew, pew_nowt))
#so this is exactyl chad's se now
rake_null_se 

#vs weighted SEs
#weights here sum to npop
#weights = out$weights$kpop_w*nrow(cces)/nrow(pew)
weights = rep(1, nrow(pew))*nrow(cces)/nrow(pew)
weights
#weights = weights/sum(weights)
sw = sum(weights)
#sum to pop N
sw
Y = pew$diff_cces_on_pew
xbar <- sum(weights * Y)/sw
xbar
#yep we are replicating the fixed w as is:
sqrt(sum(weights * ((Y - xbar)^2))/(sw - 1)/(nrow(pew)))
rake_null_se$SE_fixed
#but chad's SEs are??? interesting... well maybe fine bc yeah theres no residuals here
sqrt(sum(weights * ((Y - xbar)^2))/(sw - 1/nrow(pew))/(nrow(pew)))


#get fixed weights manually w weights that normalize to 1:
weights = rep(1, nrow(pew))*nrow(cces)/nrow(pew)
weights
weights = weights/sum(weights)
sum(weights)
xbar <- sum(weights * Y)
xbar
#drop sw denom
#so slighty differnty but that's the sample correction i think
sqrt(sum(weights * ((Y - xbar)^2)) /(nrow(pew)))
#can i re-intoduce it w weights that sum to 1... no bc the 1/N is in the weights now so they should all really be 
#1/N - 1 (bessel correction but now it has to be divided equally among pew units) im too brain dead to make that work whatever
#does this now match teh fixed weights SE, almost but minus the bessel correction so sure
sqrt(sum(weights * ((Y - xbar)^2))/(nrow(pew)))
#now compare back to Svy package
sqrt(var)
#and it IS THE SAME SO ITS JUST THE BESSEL CORRECTION IN THE PEW_NOWT


#but wtf is goign on w kpop?
out$SEs["kpop",]
data.frame(svymean(~diff_cces_on_pew, kpop))
#yeah so why these don't match the fixed weights w kpop is fucking weird:

#replicate these SEs first
norm_w = kpop_w/sum(kpop_w)
sum(norm_w)
#woutcome vector - weighted mean (est)
y =  matrix(pew$diff_cces_on_pew)
est = sum(matrix(pew$diff_cces_on_pew)*norm_w)
est
w_y = (y - est)*norm_w
var = crossprod(w_y)
sqrt(var)
var
#how is this different from fixed weight ses?
sqrt(sum(norm_w * ((Y - est)^2))/(nrow(pew)) )
#so is this bc the demaeanved vector is squared first rather than the weighted and squared?
y_2 = Y -est
#sum squared vector is ofc the same as the cross prod
sum(y_2^2)
crossprod(y_2)
#buuuut now w the weights:
sum(norm_w*y_2^2)
crossprod(norm_w*y_2)


##################### second look at fixed w SE estimator


#let's use random weights:
weights = runif(100)
Y = rnorm(100, 5,2)
N =  500
#make weights sum to Npop
weights_N = weights/sum(weights)*N
sum(weights)
sw = sum(weights_N)
sw
xbar <- sum(weights_N * Y)/sw
xbar
#same as
weights_1 = weights/sum(weights)
sum(weights_1)
sum(weights_1*Y)
xbar

#Now Fixed W Estimator:
#should want weights_N bc of bessel corr
Hmisc::wtd.var(Y, weights_N)
#does it complain if we use weights_1? yes
Hmisc::wtd.var(Y, weights_1)
#replicate manually
sw
#why is this off by 100?
sum(weights_N * ((Y - xbar)^2))/(sw - 1)/(length(weights))
#ah... it's not divided by n
sum(weights_N * ((Y - xbar)^2))/(sw - 1)


#what is normwt doing
Hmisc::wtd.var(Y, weights_N, normwt = T)
Hmisc::wtd.var(Y, weights_1, normwt = T)
weights_n = weights/sum(weights)*length(weights)
sum(weights_n)
sum(weights_n * ((Y - xbar)^2))/(sum(weights_n) - 1)

stats::cov.wt(cbind(Y), weights_n)$cov



#wait but the documentation says normwt = F corresponds to "frequency weights" which are w_i>= 1 the number of times that unit appears
#in the sample, st that the s


#with wikipedia:
#with FREQ WEIGHTS THAT SUM TO N_POP
sum(weights_N*(Y - xbar)^2)/(sum(weights_N) - 1)
#which yes is this
Hmisc::wtd.var(Y, weights_N)
sqrt(Hmisc::wtd.var(Y, weights_N)/length(weights))

#Reliability weights: unclear what they are supposed to sum to....
#ok so it is the same but WTF DOES IT MEAN
sum(weights_N*(Y - xbar)^2)/(sum(weights_N) - sum(weights_N^2)/sum(weights_N))
Hmisc::wtd.var(Y, weights_N, normwt = T)
sqrt(Hmisc::wtd.var(Y, weights_N, normwt = T)/length(weights))


#and the survey package is doign:sum (1/w)^2 Y^2
test = svydesign(~1, weights = weights_N, data = data.frame(Y))
svymean(~Y, test)
xbar
(data.frame(svymean(~Y, test))[1,2])^2*length(weights)
w_y = (Y - xbar)*weights_1
var = crossprod(w_y)
var

whyme = svymean(~Y, test)
attr(whyme, "var")
sqrt(attr(whyme, "var"))
whyme

var/attr(whyme, "var")
#is this the correction factor?
(sum(weights_N) - 1)/sum(weights_N)



test_2 = svydesign(~1, weights = weights_1, data = data.frame(Y))
again = svymean(~Y, test_2)
attr(again, "var")
attr(whyme, "var")
#this makes no difference 



###### second test:
n = 2429
weights = runif(n)
Y = rnorm(n, 5,2)
N =  1234234
weights_N = weights/sum(weights)*N
sum(weights_N)
test = svydesign(~1, weights = weights_N, data = data.frame(Y))
svmean_obj = svymean(~Y, test)
svmean_obj
attr(svmean_obj, "var")
sqrt(attr(svmean_obj, "var"))

#replicate manually
xbar = sum(Y*weights_N)/sum(weights_N)
xbar
weights_1 = weights/sum(weights)
w_y = (Y - xbar)*weights_1
var = w_y %*% w_y
var
crossprod(w_y)
#now add the bessel correction:
bessel = (sum(weights_N) - 1)/sum(weights_N)
bessel*var
attr(svmean_obj, "var")
var

bessel_2 = sum(weights_N)/(sum(weights_N)-1)
bessel_2
attr(svmean_obj, "var")/var
var*bessel_2
var

bessel_3 = sum(weights_N)/(sum(weights_N)-2)
bessel_3
#bessel_2 = (sum(weights_N) - sum(weights_N^2)/sum(weights_N))
var/(1-mean(1+1))


why = svymean(~Y, test, deff = F)
attr(why, "var")


#### how to replicate w weights that sum to N
w_y = (Y - xbar)*weights_1
var = w_y %*% w_y
var

w_y_2 = (Y - xbar)*weights_N/sum(weights_N)
w_y_2 %*% w_y_2
var
sum( (weights_N/sum(weights_N))^2*(Y-xbar)^2)
sum( weights_1*)


# 
# sqrt(nrow(cces) -1)
# #######sick of this lifting directly from teh pacakge
# pweights<-1/kpop$prob
# psum<-sum(pweights)
# x = matrix(kpop$variables$diff_cces_on_pew)
# dim(x)
# average<-colSums(x*pweights/psum)
# #ok so yes here is how he's getting the mean
# average
# xbar
# x
# #now variance
# #x<-sweep(x,2,average)
# x
# #wtf is sweep
# A <- array(1:24, dim = 4:2)
# A
# ## no warnings in normal use
# sweep(A, 1, 5)
# #ok so sweep is by default subtracting;
# #so this is dthe matrix x of our outcome - the weighted average
# x<-sweep(x,2,average)
# test = kpop$variables$diff_cces_on_pew - xbar
# #exactly teh same... so unnecessary? but ok
# cbind(test,x)
# 
# 
# #now the variance:
# design = kpop
# v<-svyCprod(x*pweights/psum,design$strata,design$cluster[[1]], design$fpc,
#             design$nPSU,design$certainty, design$postStrata)
# 
# design$cluster[[1]]
# strata = design$strata
# design$fpc
# postStrata = design$postStraa
# temp = x*pweights/psum
# 
# 
# x = x*pweights/psum
# strata = design$strata
# psu = design$cluster[[1]]
# fpc = design$fpc
# nPSU = design$nPSU
# certainty = design$certainty
# postStrata = design$postStrata
# kpop$deff
# 
# svyCprod<-function(x, strata, psu, fpc, nPSU, certainty=NULL, postStrata=NULL,
#                    lonely.psu=getOption("survey.lonely.psu")){
#     
#     x<-as.matrix(x)
#     n<-NROW(x)
#     
#     ## Remove post-stratum means, which may cut across PSUs
#     if(!is.null(postStrata)){
#         for (psvar in postStrata){
#             if (inherits(psvar, "greg_calibration") || inherits(psvar, "raking"))
#                 stop("rake() and calibrate() not supported for old-style design objects")
#             psw<-attr(psvar,"weights")
#             psmeans<-rowsum(x/psw,psvar,reorder=TRUE)/as.vector(table(factor(psvar)))
#             x<- x-psmeans[match(psvar,sort(unique(psvar))),]*psw
#         }
#     }
#     
#     ##First collapse over PSUs
#     
#     if (is.null(strata)){
#         strata<-rep("1",n)
#         if (!is.null(nPSU))
#             names(nPSU)<-"1"
#     }else {
#         strata<-as.character(strata) ##can't use factors as indices in for()'
#     }
#     
#     if (is.null(certainty)){
#         certainty<-rep(FALSE,length(strata))
#         names(certainty)<-strata
#     }
#     
#     if (!is.null(psu)){
#         x<-rowsum(x, psu, reorder=FALSE)
#         strata<-strata[!duplicated(psu)]
#         n<-NROW(x)
#     }
#     
#     if (!is.null(nPSU)){
#         obsn<-table(strata)
#         dropped<-nPSU[match(names(obsn),names(nPSU))]-obsn
#         if(sum(dropped)){
#             xtra<-matrix(0,ncol=NCOL(x),nrow=sum(dropped))
#             strata<-c(strata,rep(names(dropped),dropped))
#             if(is.matrix(x))
#                 x<-rbind(x,xtra)
#             else
#                 x<-c(x,xtra)
#             n<-NROW(x)
#         }
#     } else obsn<-table(strata)
#     
#     if(is.null(strata)){
#         x<-t(t(x)-colMeans(x))
#     } else {
#         strata.means<-drop(rowsum(x,strata, reorder=FALSE))/drop(rowsum(rep(1,n),strata, reorder=FALSE))
#         if (!is.matrix(strata.means))
#             strata.means<-matrix(strata.means, ncol=NCOL(x))
#         x<- x- strata.means[ match(strata, unique(strata)),,drop=FALSE]
#     }
#     
#     p<-NCOL(x)
#     v<-matrix(0,p,p)
#     
#     ss<-unique(strata)
#     for(s in ss){
#         this.stratum <- strata %in% s
#         
#         ## original number of PSUs in this stratum 
#         ## before missing data/subsetting
#         this.n <-nPSU[match(s,names(nPSU))]
#         
#         this.df <- this.n/(this.n-1)	
#         
#         if (is.null(fpc))
#             this.fpc <- 1
#         else{
#             this.fpc <- fpc[,2][ fpc[,1]==as.character(s)]
#             this.fpc <- (this.fpc - this.n)/this.fpc
#         }
#         
#         xs<-x[this.stratum,,drop=FALSE]
#         
#         this.certain<-certainty[names(certainty) %in% s]
#         
#         ## stratum with only 1 design cluster leads to undefined variance
#         lonely.psu<-match.arg(lonely.psu, c("remove","adjust","fail",
#                                             "certainty","average"))
#         if (this.n==1 && !this.certain){
#             this.df<-1
#             if (lonely.psu=="fail")
#                 stop("Stratum ",s, " has only one sampling unit.")
#             else if (lonely.psu!="certainty")
#                 warning("Stratum ",s, " has only one sampling unit.")
#             if (lonely.psu=="adjust")
#                 xs<-strata.means[match(s,ss),,drop=FALSE]
#         } else if (obsn[match(s,names(obsn))]==1 && !this.certain){
#             ## stratum with only 1 cluster left after subsetting 
#             warning("Stratum ",s," has only one PSU in this subset.")
#             if (lonely.psu=="adjust")
#                 xs<-strata.means[match(s,ss),,drop=FALSE]
#         }
#         ## add it up
#         if (!this.certain)
#             v<-v+crossprod(xs)*this.df*this.fpc
#     }
#     if (lonely.psu=="average"){
#         v<- v/(1-mean(obsn==1 & !certainty))
#     }
#     v
# }