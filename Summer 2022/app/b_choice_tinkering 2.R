#### testing with mixed kernel
X_1 = c(1, 1, 0)
X_2 = c(2,1,0)
X_3 = rnorm(3)
sd(X_3)
X_3 = c(X_3[1], X_3[1]+sd(X_3), X_3[1]+2*sd(X_3) )
X_3
X_noscale = matrix(c(X_1, X_2, (X_3)), nrow = 3, ncol = 3)
X_noscale
X = matrix(c(X_1, X_2, scale(X_3)), nrow = 3, ncol = 3)
X
raw_diff = abs(X[1,] - X[2,])
raw_diff
sum(raw_diff)^2


#one_hot
cat_cols =  c(1,2)
X_onehot = cbind(one_hot(X[,cat_cols]), X[, -cat_cols])
X_onehot

raw_diff_onehot = abs(X_onehot[1,] - X_onehot[2,])
raw_diff_onehot

sum(raw_diff)
sum(raw_diff_onehot)

#### definitely one option that will work is to multiple the cont col by 2 right>
X_onehot_2 = cbind(one_hot(X[,cat_cols]), 2*X[, -cat_cols])
X_onehot_2
raw_diff_onehot_2 = abs(X_onehot_2[1,] - X_onehot_2[2,])

sum(raw_diff_onehot_2)/2
sum(raw_diff)


########also no *2 but by 2/3 works right? A
#waaait but what if they are the same on cat.... 
sum(raw_diff_onehot)*(2/3)



##### ok so then what best b is being returned. let's make some k's
exp(-sum(raw_diff)/1)
exp(-sum(raw_diff_onehot)*(2/3))
exp(-sum(raw_diff_onehot_2)/2)

K_raw = makeK(allx = X, b = 1, scale = F)
K_raw

K_onehot2 = makeK(allx = X_onehot_2, b = 2, scale = F)
K_onehot2

#but yeah this only correctly scales where they are off by 1 + 1sd, = 3 
K_onehot = makeK(allx = X_onehot, b = (3/2), scale = F)
K_onehot
#to be correct for where they are off by 2 + 1s = 4+1= 5 (units 2,3)

makeK(allx = X_onehot, b = (5/2), scale = F)
K_raw
#for units 1 and 3: 1 + 2s = 2 + 2s = 4?
makeK(allx = X_onehot, b = (4/2), scale = F)

#lol yeah like that doesnt workkkkk
#but why not
K_onehot2 = makeK(allx = X_onehot_2, b = 1, scale = F)
K_onehot2
sum((X_onehot_2[1,] - X_onehot_2[2])^2)
exp(-sum((X_onehot_2[1,] - X_onehot_2[2])^2)/1)
#oh its because it's fucking squared so yeah you need to multiply by...the sqrt of 2???
X_onehot_2sqrt = cbind(one_hot(X[,cat_cols]), sqrt(2)*X[, -cat_cols])
sum((X_onehot_2sqrt[1,] - X_onehot_2sqrt[2])^2)

K_onehot2sqrt = makeK(allx = X_onehot_2sqrt, b = 2, scale = F)
K_onehot2sqrt
K_raw
#FINALLY


#### now check with more levles:
#X1 = cat 2 levels
X_1 = c(1, 1, 0)
#X2 = cat 3 levels
X_2 = c(2, 1,0)
#X3 = cont with differences 1sd for simplicity
X_3 = rnorm(3)
X_3 = c(X_3[1], X_3[1]+sd(X_3), X_3[1]+2*sd(X_3) )
X_3
X = matrix(c(X_1, X_2, scale(X_3)), nrow = 3, ncol = 3)
#final resulting X
X
cat_cols = c(1,2)
X_onehot_2sqrt = cbind(one_hot(X[,cat_cols]), sqrt(2)*X[, -cat_cols])

K_onehot2sqrt = makeK(allx = X_onehot_2sqrt, b = 2, scale = F)
K_onehot2sqrt

#compare with:
K_raw = makeK(allx = X, b = 1, scale = F)
#this is correct and matches the above for unit 1 and 2, and 2 and 3
K_raw[1,2]
K_onehot2sqrt[1,2]
K_raw[3,2]
K_onehot2sqrt[3,2]
#but not unit 1 and 3 since it for categorical column 2 its counting (2-0) = 2 instead of as 1
#as we can see ourselves directly ofc
K_raw[3,1]
K_onehot2sqrt[3,1]
#we go wrong here bc of X2 having 3 levels
X[c(1,3),]
raw_diff = (X[1,] - X[3,])**2
raw_diff
exp(- sum(raw_diff)/1)
K_raw[1,3]
#but we can manually calcualte and see that K_onehot2sqrt[1,3] is correct
raw_diff[2] =1
exp(- sum(raw_diff))
K_onehot2sqrt[1,3]



######################### b max var K ################################
X_1 = c(1, 1, 0)
X_2 = c(0,1,0)
X_3 = rnorm(3)
sd(X_3)
X_3 = c(X_3[1], X_3[1]+sd(X_3), X_3[1]+2*sd(X_3) )
X_3
X = matrix(c(X_1, X_2, scale(X_3)), nrow = 3, ncol = 3)
X
#no double counts:
K_raw = makeK(allx = X, b = 1, useasbases = rep(1,3), scale = F)
K_raw
#for sanity:
raw_diff = (X[1,] - X[2,])**2
exp(- sum(raw_diff))

bout = b_maxvarK(data = X, cat_data = F, useasbases = rep(1,3))
bout
#make best K
K_best = makeK(allx = X, useasbases = rep(1,3), b = bout$b_maxvar, scale = F)
K_best
#double check
raw_diff = (X[2,] - X[3,])**2
exp(- sum(raw_diff)/bout$b_maxvar)
#check var:
index = upper.tri(K_best) | lower.tri(K_best)
var(K_best[index])
bout$var_K


#### now with double counts
cat_cols = c(1,2)
X_onehot = cbind(one_hot(X[, cat_cols]), sqrt(2)*scale(X[, -cat_cols]))
X_onehot
#still can't use cat_data = T bc of Xw
b_double = b_maxvarK(data= X_onehot, cat_data = F, useasbases = rep(1,3))
b_double
bout
#so the max variance is the same, but the b is double for 2x counts
bout$b_maxvar*2
b_double$b_maxvar
#so to make K_best with the X_onehot we need to? use the double or not?
K_double = makeK(allx = X_onehot, useasbases = rep(1,3), scale = F, b =b_double$b_maxvar)
K_double
K_best
#yeahhh fucking duh it's the same b! the be we found is doubled!

#double check with cat data what happens:
X_cat = X[, cat_cols]
b = b_maxvarK(data = X_cat, useasbases = c(1,1), cat_data = F,maxsearch_b = 3)
b
#NOOOO fuck why is this different than cat+data = T???
#is the var calc wrong? like no way:

var_K_man = function(b, data) {
    K <- makeK(data, b = b, useasbases = useasbases, 
               linkernel = FALSE, scale = FALSE)
    diag(K) <- NA
    n = nrow(K) * ncol(K) - sum(useasbases)
    var_k = (1/(n - 1)) * (sum(K^2, na.rm = T) - (1/n) * 
                               sum(K, na.rm = T)^2)
    return(var_k)
}

var_K = function(b, n_d, diag_length) {
    d <- n_d[, 1] %>% pull()
    n_d <- as.vector(n_d[, 2] %>% pull())
    n_d[1] <- n_d[1] - diag_length
    p_d <- n_d/sum(n_d)
    mean_k = sum(exp(-1 * d/b) * p_d)
    var_k = sum((exp(-1 * d/b) - mean_k)^2 * p_d)
    return(var_k)
}
##### onehot way:
#no onehot encoding, so b = 1 here
K <- makeK(allx = X_cat, b = 1, useasbases = rep(1,1,1), linkernel = FALSE, 
           scale = FALSE)
K
raw_counts <- -log(K)
raw_counts
n_d <- data.frame(diff = c(raw_counts)) %>% group_by(diff) %>% 
    summarise(n())
n_d
var_K(b = 1, n_d = n_d, diag_length = 3)
var_K_man(b = 1, data = X_cat)
#FUCK WHY ARE THESE NOT THE SAME? which is right?
k_test = makeK(allx = X_cat, useasbases = rep(1, 3), b = 1, scale = F)
var(k_test[index])
#shit ok so the base way is right?
#why is this wrong? small sample corr?
var_K(b = 1, n_d = n_d, diag_length = 3)
res = k_test[index]
mean(res^2) - mean(res)^2
var(res)
#it is the small n factor then!!
sum(res)/length(res)
mean(res)
#ok so we have shown that the onehot way calculates the correct variance WITHOUT the small sample corr
#if i go change the var_K_man to not have the small sample correction will they then match?
var_K_man_2 = function(b, data) {
    K <- makeK(data, b = b, useasbases = useasbases, 
               linkernel = FALSE, scale = FALSE)
    diag(K) <- NA
    n = nrow(K) * ncol(K) - sum(useasbases)
    var_k = (1/(n)) * (sum(K^2, na.rm = T) - (1/n) * 
                               sum(K, na.rm = T)^2)
    return(var_k)
}
var_K_man_2(b = 1, data = X_cat)
var_K(b = 1, n_d = n_d, diag_length = 2)
#FOR FUCKS SAKE. heart attaack over noooothing




##### ok so back to now to b or 2b for cat data only
b_maxvarK_man<- function (data, useasbases, cat_data = TRUE, maxsearch_b = 2000) 
{
    if (cat_data) {
        K <- makeK(data, b =2, useasbases = useasbases, linkernel = FALSE, 
                   scale = FALSE)
        raw_counts <- -log(K)
        n_d <- data.frame(diff = c(raw_counts)) %>% group_by(diff) %>% 
            summarise(n())
        var_K = function(b, n_d, diag_length) {
            d <- n_d[, 1] %>% pull()
            n_d <- as.vector(n_d[, 2] %>% pull())
            n_d[1] <- n_d[1] - diag_length
            p_d <- n_d/sum(n_d)
            mean_k = sum(exp(-1 * d/b) * p_d)
            var_k = sum((exp(-1 * d/b) - mean_k)^2 * p_d)
            return(var_k)
        }
        res = optimize(var_K, n_d, length(diag(K)), interval = c(0, 
                                                                 maxsearch_b), maximum = TRUE)
    }
    else {
        var_K = function(b, data) {
            K <- makeK(data, b = b, useasbases = useasbases, 
                       linkernel = FALSE, scale = FALSE)
            diag(K) <- NA
            n = nrow(K) * ncol(K) - sum(useasbases)
            var_k = (1/(n)) * (sum(K^2, na.rm = T) - (1/n) * 
                                       sum(K, na.rm = T)^2)
            return(var_k)
        }
        res = optimize(var_K, data, interval = c(0, maxsearch_b), 
                       maximum = TRUE)
    }
    return(list(b_maxvar = res$maximum, var_K = res$objective))
}

b = b_maxvarK_man(data = X_cat, useasbases = c(1,1,1), cat_data = F)
b

X_onehot_cat = one_hot(X_cat)
b2 = b_maxvarK_man(data = X_onehot_cat, useasbases = c(1,1,1), cat_data = T)
b2
#same as
b22 = b_maxvarK(data = X_onehot_cat, useasbases = c(1,1,1), cat_data = T)
b22$b_maxvar*2
#GOD IM SO FUCKING DUMb
#FINALLY THEY ARE IDENTICAL fuck's sake. ughguhguhgughughg all over a small sample correction
#so now if we make K with this b though do we get the right K? seems like not
K_onehot = makeK(allx = X_onehot_cat, b = b2$b_maxvar, useasbases = rep(1,3), scale = F)
K_onehot

K = makeK(allx = X_cat, b = b$b_maxvar, useasbases = rep(1,3), scale = F)
K

K_onehot = makeK(allx = X_onehot_cat, b = 2*b2$b_maxvar, useasbases = rep(1,3), scale = F)
K_onehot
#RIGHT OK SO in the cat_data = T ONLY best b search...the b from cont and cat case are the same so to achieve the actually correct K we need to *2?
#YES BC im dividing the actual counts by 2 so the double counting doesn't fucking matter  bc K in the optimiation always is made with b = 2.... so like 
#one can either do K(b = 2) + b = 2*best_b, orrrr one can just do K(b=1), b = best_b
#i am so fucking stupid
#but for the mixed or cont case... we? dont?
b_nocat = b_maxvarK_man(data = X, useasbases = rep(1,3),cat_data = F)
b_nocat
b_cat = b_maxvarK_man(data = cbind(one_hot(X[,cat_cols]), sqrt(2)*scale(X[,-cat_cols])), 
                                   useasbases = rep(1,3), cat_data = F)
b_cat                      
b_nocat
#RIGHT ofc bc both are using the cat_data = F variance calc so the b that comes out is right for that inputted data



##### Chad's method:
X
X_cat = X[,1:2]
X_onehot_chad = cbind(sqrt(0.5)*one_hot(X[,cat_cols]), X[, -cat_cols] )
X_onehot_chad
makeK(X_onehot_chad, b = 1, scale = F)
K_raw
#many confusing things going on at once:
#first confirm: b=2 + cat requires b*2 with inputted data to be the same as cont
b_oh_b2 = b_maxvarK_man(data = one_hot(X_cat), useasbases = rep(1,3), cat_data = T)
b_cont_dir = b_maxvarK_man(data = X_cat, useasbases = rep(1,3), cat_data = F)
#this produced THE SAME CHOICE OF B
b_oh_b2
b_cont_dir
#so to produce the resulting same K, the OH version needs to be *2
makeK(allx = one_hot(X_cat), b = 2*b_oh_b2$b_maxvar, useasbases = rep(1,3), scale = F)
makeK(allx = X_cat, b = b_cont_dir$b_maxvar, useasbases = rep(1,3), scale = F)
#indeed

#second: make b=1 in cat maxvar and do not multiply by 2:
b_maxvarK_man_b1<- function (data, useasbases, cat_data = TRUE, maxsearch_b = 2000) 
{
    if (cat_data) {
        K <- makeK(data, b =1, useasbases = useasbases, linkernel = FALSE, 
                   scale = FALSE)
        raw_counts <- -log(K)
        n_d <- data.frame(diff = c(raw_counts)) %>% group_by(diff) %>% 
            summarise(n())
        var_K = function(b, n_d, diag_length) {
            d <- n_d[, 1] %>% pull()
            n_d <- as.vector(n_d[, 2] %>% pull())
            n_d[1] <- n_d[1] - diag_length
            p_d <- n_d/sum(n_d)
            mean_k = sum(exp(-1 * d/b) * p_d)
            var_k = sum((exp(-1 * d/b) - mean_k)^2 * p_d)
            return(var_k)
        }
        res = optimize(var_K, n_d, length(diag(K)), interval = c(0, 
                                                                 maxsearch_b), maximum = TRUE)
    }
    else {
        var_K = function(b, data) {
            K <- makeK(data, b = b, useasbases = useasbases, 
                       linkernel = FALSE, scale = FALSE)
            diag(K) <- NA
            n = nrow(K) * ncol(K) - sum(useasbases)
            var_k = (1/(n)) * (sum(K^2, na.rm = T) - (1/n) * 
                                   sum(K, na.rm = T)^2)
            return(var_k)
        }
        res = optimize(var_K, data, interval = c(0, maxsearch_b), 
                       maximum = TRUE)
    }
    return(list(b_maxvar = res$maximum, var_K = res$objective))
}
b_oh_b1 = b_maxvarK_man_b1(data = one_hot(X_cat), useasbases = rep(1,3), cat_data = T)
b_cont_dir = b_maxvarK_man_b1(data = X_cat, useasbases = rep(1,3), cat_data = F)
#now we do not need to multply by 2:
#bc b_oh_b1 is already 2x b_cont_dir
b_oh_b1
b_cont_dir
makeK(allx = one_hot(X_cat), b = b_oh_b1$b_maxvar, useasbases = rep(1,3), scale = F)
makeK(allx = X_cat, b = b_cont_dir$b_maxvar, useasbases = rep(1,3), scale = F)

#third: using chad's method
#onehot is now not double counted, so.... should it not just work as planned?
X_chad = sqrt(0.5)*one_hot(X_cat)
X_chad
b_oh_b1_chad = b_maxvarK_man_b1(data = X_chad, useasbases = rep(1,3), cat_data = T)
b_oh_b1_chad
b_cont_dir
makeK(allx = X_chad, b = b_oh_b1_chad$b_maxvar, useasbases = rep(1,3), scale = F)
makeK(allx = X_cat, b = b_cont_dir$b_maxvar, useasbases = rep(1,3), scale = F)
#this is still 2* but gets the right result in K bc of the 2* in OH X
b_oh_b1
makeK(allx = one_hot(X_cat), b = b_oh_b1$b_maxvar, useasbases = rep(1,3), scale = F)


## just to double check, cehck with more than one level (can no longer compare with cont)
X_cat[1,2] =2
X_cat
X_chad = sqrt(0.5)*one_hot(X_cat)
X_chad
b_oh_b1_chad = b_maxvarK_man_b1(data = X_chad, useasbases = rep(1,3), cat_data = T)
b_oh_b1 = b_maxvarK_man_b1(data = one_hot(X_cat), useasbases = rep(1,3), cat_data = T)
b_oh_b1_chad
#this is twice teh size as expected
b_oh_b1
makeK(allx = X_chad, b = b_oh_b1_chad$b_maxvar, useasbases = rep(1,3), scale = F)
makeK(allx = one_hot(X_cat), b = b_oh_b1$b_maxvar, useasbases = rep(1,3), scale = F)
#wahoooo


#so smartest thing to do is: change to b=1 in maxvarK and *by sqrt(0.5) in onehot
#should we also do that for fully categorical? maybe
#it literally doesnt matter bc the choice of b will compensate but maybe for consistent scale of b for comapriosn? less confusing? or more confusing lol i'll do it for now



########### ok package updated, now test:
devtools::install_github("csterbenz1/KBAL", ref ="cat_kernel", force = T)
detach("package:kbal", unload = TRUE)
library(kbal)
b_maxvarK

#test w cat only but more than 2 levels:
b_maxvarK(data = one_hot(X_cat), useasbases = rep(1,3), cat_data = T)
b_maxvarK(data = X_chad, useasbases = rep(1,3), cat_data = T)
makeK(allx = X_chad, b = b_maxvarK(data = X_chad, useasbases = rep(1,3), cat_data = T)$b_maxvar, useasbases = rep(1,3), scale = F)
makeK(allx = one_hot(X_cat), b = b_maxvarK(data = one_hot(X_cat), useasbases = rep(1,3), cat_data = T)$b_maxvar, useasbases = rep(1,3), scale = F)

#ok so now let's try the full monte:
#first try a case where you can do a cont all vs mixed
X
#shit
test_mixed = kbal(allx = X, sampled = c(1,0,0),
                  cat_columns = c(1,2), 
                  cont_scale = 1, 
                  mixed_data = T)

test_cont = kbal(allx = X, sampled = c(1,0,0),  useasbases = rep(1,3), scale_data = F)
test_cont$K
test_mixed$K
test_cont$b
test_mixed$b
test_mixed$onehot_data



K_man = makeK(allx = X, b= test_cont$b, useasbases = rep(1,3), scale = F)
K_man
test_cont$K
#so it is making the correct kernel matrix with the right choice of b
#variance is also correct
test_cont$maxvar_K
mean(test_cont$K[index]^2)- mean(test_cont$K[index])^2



#next compare mixed and cat only: looks good
test_cat = kbal(allx = X[,c(1,2)], 
                sampled = c(1,0,0),
                cat_data = T)
test_cont_c = kbal(allx = X[,c(1,2)], 
                   sampled = c(1,0,0),
                   scale_data = F)
test_cat$K
test_cont_c$K
test_cat$b
test_cont_c$b
test_cat$maxvar_K
test_cont_c$maxvar_K

#now how is it doing the in house rescaling of the cont data?
test_mixed_2 = kbal(allx = X, sampled = c(1,0,0),
                    cat_columns = c(1,2), 
                    cont_scale = 2, 
                    mixed_data = T, 
                    scale_data = T)
apply(test_mixed_2$onehot_data, 2, sd)
apply(test_mixed_2$onehot_data, 2, var)


test_mixed_2$onehot_data

cat_columns = c(1,2)
allx = X
cont_scale = 2
allx_cont <- t(t(allx[, -cat_columns, drop = F])/(apply(allx[, -cat_columns, drop = F], 2, sd)*(1/cont_scale)))
apply(allx_cont, 2, sd)


apply(scale(test_mixed_2$onehot_data), 2, sd)
#also check where is scale_data being used? only in cont only case?