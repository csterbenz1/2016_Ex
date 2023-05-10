


############## Simualte 

N = 295939
n = 2534
Y = rnorm(N, 8,2)
y = rnorm(n, 8,2)


#prob-like weights:
#should sum to n
pi = runif(n)
sum(pi)
pi = pi/sum(pi)*n 
sum(pi)
pi = rep(1/n,n)
pi = rep(n/N,n)
sum(pi)
w = 1/pi
sum(w)
N
sum(1/w)



#frequency weights: w >= 1 sum(w) = N
w_N = w/sum(w)*N
sum(w_N)
#normalized weights
w_1 = w/sum(w)

#replicate Hmsic::
sum(w_N*(y - sum(w_N*y)/sum(w_N) )^2)/(sum(w_N) - 1)
Hmisc::wtd.var(y, weights = w_N)

#HT:
(1/N)^2*(n/(n-1))*sum((w*y - sum(w*y/n))^2)
(1/N)^2*(n/(n-1))*sum((pi*y - sum(pi*y/n))^2)
(1/N)^2*(n/(n-1))*sum((w_N*y - sum(w_N*y/n))^2)


mean(Y)
mean(y)
sum(pi*y)
sum(w_N*y)/sum(w_N)
sum(w_1*y)/sum(w_1)
sum(w*y)/sum(w)
sum(w*y)/N
sum(w)


####ok so HT expression is w weights that sum to N
(1/N)^2*(n/(n-1))*sum((w_N*y - sum(w_N*y/n))^2)
(1/N)^2*(n/(n-1))*sum((w_N*y - mean(w_N*y))^2)

test = svydesign(~1, weights = w_N, data = data.frame(y))
attr(svymean(~y, test), "var")
#which is the same! finally:

#before i was trying to replicate the svy pacakge w:
est = sum(w_1*y)
w_y = (y - sum(w_1*y))*w_1
w_y %*% w_y
#whcih was close but not the same
attr(svymean(~y, test), "var")
#bc bessel corr is n/n-1
(w_y %*% w_y)*(n/(n-1))

#just for sanity of translation to w_N
(1/N)^2*(n/(n-1))*sum((w_N*y - mean(w_N*y))^2)
#move N factor back in
(n/(n-1))*sum((w_N*y/N - mean(w_N*y)/N)^2)
#now we have:
(n/(n-1))*sum((w_1*y - mean(w_1*y))^2)
#and this is now
w_y_2= (w_1*y - mean(w_1*y))
cbind(w_y, w_y_2)
cbind(w_y, w_y_2, (y*w_1 - sum(w_1*y)/n) )
#? confused as to why this works but ok...
cbind(w_y, w_y_2, w_1*(y - sum(w_1*y)) )
(y*w_1 - sum(w_1*y))*w_1
sum(w_1/y)


(n/(n-1))*sum((w_1*y - mean(w_1*y))^2)