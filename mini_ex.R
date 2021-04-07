##################### Mini Example in Paper #####################

##########################################
##################### Set up ##############
#1000 unit version
pop <- data.frame(
    republican =  c(rep(1,400), rep(0,400)),
    white = c(rep(1,200), rep(0,200), rep(1,200), rep(0,200)),
    support = c(rep(.8,200), rep(.2,600)))
unique(pop) # white republicans have higher support
table(pop[,1:2])/nrow(pop)
mean(pop$support)  # Target value 0.35

#sample build : correct margins/means, but wrong joint distribution
# 100 unit version
samp <- data.frame( republican = c(rep(1, 40), rep(0,40)),
                    white = c(rep(1,30), rep(0,10), rep(1,10), rep(0,30)),
                    support = c(rep(.8,30), rep(.2,50)))
unique(samp)
table(samp[,1:2])/nrow(samp)
mean(samp$support) #we get .425 because white republicans over represented





##########################################
################ Reweight #################

dat <- as.matrix(rbind(pop,samp))
X <- dat[,1:2]
sampled <- c(rep(0,800), rep(1,80))

#### Raking:
# Run ebal (treatment = population units = 1-sampled)
ebal_out <- kbal::ebalance_custom(Treatment = 1-sampled,
                                  X=X,
                                  constraint.tolerance=1e-6,
                                  print.level=-1)

#already have balanceon margins so have even weights
length(unique(ebal_out$w))
unique(cbind(samp, e_bal_weight = ebal_out$w))
#estimate does not change .425 as we expect
weighted.mean(samp[,3], w = ebal_out$w) 

########## Kbal

#1. Understanding the kernel: compute raw "misses" 
#Kernel matrix is simply: K[i,j] = e^(- sum of raw "misses" between unit i and unit j across all columns of the data) = e^-(euclidean distance between i and j)
#specifically it is the sum of misses with a one-hot encoding of all variables then divided by 2 to compensate for this double counting hence b =2 here
#in this case that is explicitly sum( (rep_i - rep_j)^2 + (notrep_i - notrep_j)^2 + (white_i - white_j)^2 + (notwhite_i - notwhite_j)^2)

#NB: the way to force your df to be factors depends on your version of R, this is for 3.6
onehot_data <- data.frame(apply(X,2, as.factor))
onehot_data <- model.matrix(~ ., onehot_data,
                          contrasts.arg = lapply(onehot_data, contrasts, contrasts=FALSE))
onehot_data <- onehot_data[, -1]
colnames(onehot_data)

K <- kbal::makeK(onehot_data, b=2, useasbases = rep(1, nrow(onehot_data)),
           linkernel = FALSE, scale = FALSE)

raw_counts <- -log(K)
#let see how this looks
#check this: K[i,j]  = exp(- sum of differences )
K[45, 345]
exp(-sum((onehot_data[45, ] - onehot_data[345,])^2)/2)
#same as no double counts:
exp(-sum((X[45, ] - X[345, ])^2))

#taking the log to get back the raw differences
raw_counts[45 ,345]
#we see a double count of one  miss
onehot_data[45, ] - onehot_data[345,]
sum((onehot_data[45, ] - onehot_data[345,])^2)/2
#and similarly without double counts
(X[45, ] - X[345,])^2
sum((X[45, ] - X[345,])^2)


#2. now we find the b that will max the variance of this Kernel
var_K= function(b, n_d){
    p_d <- n_d/ sum(n_d) 
    d = 0:(length(n_d)-1)
    mean_k = sum(exp(-1*d/b)*p_d)
    var_k = sum((exp(-1*d/b)-mean_k)^2 * p_d)
    return(var_k)
}
#run optim:
#we need to know the frequences each number of "misses" and dplyr is much faster than table
n_d <- data.frame(diff = c(raw_counts)) %>% group_by(diff) %>% summarise(n())
n_d <- as.vector(n_d[,2] %>% pull())
n_d
#visually we can see these
hist(as.vector(raw_counts))

res = optimize(var_K, n_d = n_d,
               interval=c(0,2000), maximum=TRUE)
b_maxvar <- res$maximum
b_maxvar



#3. Run kbal with this choice of b
#Now: Kernel balancing for weighting to a population (i.e. kpop) -------
kbalout = kbal::kbal(allx= X,
                     useasbases=rep(1,nrow(X)),
                     sampled = sampled,
                     ebal.convergence = FALSE,
                     b = b_maxvar,
                     fullSVD = TRUE,
                     sampledinpop = FALSE)
# The weights now vary:
plot(kbalout$w[sampled ==1], pch=16)

# And produce correct estimate:
weighted.mean(samp$support, w = kbalout$w[sampled==1])

#down weight non white republicans as we expect
unique(cbind(samp[,-3], k_bal_weight = kbalout$w[sampled==1]))
