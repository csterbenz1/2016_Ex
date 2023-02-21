


#the following is a power analysis both by direct calculation and by simulation of difference in means between two bernoulli variables (binary 0-1's); 
#we assume one variable has p(y ==1) of base_reject_prob and a constant treatment effect such that the
#other var has p(y==1) of base_reject_prob +t_effect

power_analysis <- function( n_tot = 100,
                                mu_o = 0,
                                var_t = NULL,
                                var_c = NULL, 
                                alpha = 0.05,
                                n_tasks = 40,
                                base_reject_prob = 0.02,
                                t_effect = 0.05,
                                two_tailed = TRUE, 
                                n_sims = 1000) {
        n_t = n_tot/2
        n_c = n_tot/2
        
        #you can manually specify the variance you expect in yor outcome for the t or c groups or go based on the known
        #variance of a bernoulli variable
        if(is.null(var_t) && is.null(var_c)) {
            #variance of a bernoulli var is p*(1-p)
            var_t = (base_reject_prob+t_effect)*(1-(base_reject_prob+t_effect))
            var_c = base_reject_prob*(1-base_reject_prob)
        }
        #t_stat = treat effect/SE; SE = DIM SE
        t_stat = (t_effect- mu_o)/(sqrt(var_t/n_t + var_c/n_c))
        #critical values
        df = (var_t/n_t + var_c/n_c)^2/(
            (var_t/n_t)^2/(n_t-1) + (var_c/n_c)^2/(n_c-1) )
        #can calclate the power directly:
        if(two_tailed) {
            crit_val = qnorm(1-alpha/2)
            crit_val_tdist = qt(1-alpha/2, df = df)
        } else {
            crit_val = qnorm(1-alpha)
            crit_val_tdist = qt(1-alpha, df = df)
        }
        #power is prob of erroneously rejecting the null: P(reject null | alt = T)
        #this is the ssame as 1 - p(fail to reject null|alt = T)
        # = 1 - P(t_stat < crit_val | alt = T)
        # do some math and some approx for large enought n and this simplifies to pnorm(crit - t_stat)
        calc_power = 1 - pnorm(crit_val - abs(t_stat))
        calc_power_tdist = 1 - pt(crit_val - abs(t_stat), df =df)
        
       
        #double check with a simuation
        simulation = replicate(n_sims, {
            res = data.frame(treated = c(rep(0, n_t),
                                         rep(1, n_c)), 
                             reject = NA)
            for(i in 1:n_tot) {
                if(res$treated[i] == 0) {
                    control = rbinom(n_tasks, 1, prob = base_reject_prob)
                    res$reject[i] = mean(control)
                } else {
                    treated = rbinom(n_tasks, 1, prob = base_reject_prob+t_effect)
                    res$reject[i] = mean(treated)
                }
            }
            t_outcome = res$reject[res$treated ==1]
            c_outcome = res$reject[res$treated ==0]
            t_effect_obs = mean(t_outcome) - mean(c_outcome)
            
            #SE for DIM:
            var_dim = (1/n_t)*var(t_outcome) + (1/n_c)*var(c_outcome)
            t_stat = (t_effect_obs- mu_o)/sqrt(var_dim)
            reject_null = abs(t_stat) >= crit_val
    
            #replicate is returning a T/F value for if we reject the null
            out = data.frame(reject_null = reject_null,
                             t_effect_obs = t_effect_obs, 
                             var_dim = var_dim)
            return(out)
            
        }, simplify = F) %>% bind_rows()
        #power = prob reject the null given that alt is true
        #alt here is specified to be true
        #so it's straight up just the prob we reject the null
        simulate_power = mean(simulation$reject_null)
        avg_t_effect_obs = mean(simulation$t_effect_obs)
        var_t_effect_obs = var(simulation$t_effect_obs)
        avg_var_obs = mean(simulation$var_dim)
        
        return(list(calc_power = calc_power,
                    calc_power_tdist = calc_power_tdist, 
                    sim_power = simulate_power,
                    avg_t_effect_obs = avg_t_effect_obs,
                    #we can also check that the average variance we observe in the simulation is close to our 
                    #varaince estimator for the DIM
                    bootstrap_SE = var_t_effect_obs,
                    avg_var_obs = avg_var_obs, 
                    var_estimator = var_t/n_t + var_c/n_c ))
}



##########testing if we did it right
power_analysis(n_tot = 1000, n_tasks = 1, two_tailed = F, 
               base_reject_prob = 0.03,
               t_effect = 0.02)
#ok so this looks good! our bootstrapped SE marches roughly with our average DIM estimator SE which close to to the full estimator version using the variance of a bernoulli draw 
# we do get an average t effect of aroun 2% as we specified
# we see that as expected we have slightly lower pwoer when we use the t-distirbution but for this large of n it makes little difference 
#does this match canned functions? yes
power.prop.test(n = 500, p1 = 0.03, p2 = 0.05, alternative = "one.sided")

#two tailed check:
power_analysis(n_tot = 1000, n_tasks = 1, two_tailed = T, 
               base_reject_prob = 0.02,
               t_effect = 0.03)
power.prop.test(n = 500, p1 = 0.02, p2 = 0.05, alternative = "two.sided")

#tried to code it so you could also do negative effects let's see if that works, looks like it matches
power_analysis(n_tot = 1000, n_tasks = 1, two_tailed = T, 
               base_reject_prob = 0.05,
               t_effect = -0.02)
power.prop.test(n = 500, p1 = 0.03, p2 = 0.05, alternative = "two.sided")

#####################
#now let's explor hwo the number of tasks each unit responds to effects our power:
power_analysis(n_tot = 1000, n_tasks = 1, two_tailed = F, 
               n_sims = 500,
               base_reject_prob = 0.03,
               t_effect = 0.01)
#now bumping up to 40 as in the original design
power_analysis(n_tot = 1000, n_tasks = 40, two_tailed = F, 
               n_sims = 500,
               base_reject_prob = 0.03,
               t_effect = 0.01)

#ok wow it makes a big difference; this is kind of to be expected; basically what we did was drastically decrease the variance in the outcome by taking the mean of 40 unit-level observations instead of 1
#let's explore more
#vary the number of tasks directly and track the power and variance

test_ntasks = map_dfr(1:40, function(task_no){
                        run = power_analysis(n_tot = 1000, n_tasks = task_no, 
                                       two_tailed = F, 
                                       n_sims = 500,
                                       base_reject_prob = 0.03,
                                       t_effect = 0.01)
                        out = cbind(n_tasks = task_no, t(run))
                        out = as.data.frame(out) %>% bind_rows()
})
#output of this is a bit odd... map is annoying
class(test_ntasks)
class(test_ntasks$n_tasks)
#columns are all lists let's undo that
res = as.data.frame(apply(test_ntasks, 2, function(col) {as.numeric(col)}))

#plot to see how does power in sim seem to change?
ggplot(res) + 
    geom_point(aes(x = n_tasks, y = sim_power) ) +
    theme_bw()
#we seem to get much higher simulated power when we increase the number of tasks each respondent performs
#why is this? because we went from varaince of a single bernoulli draw to variance fo an average of bernoulli draws:
#V( 1/n_tasks sum^n_tasks Y_i) = (1/n_tasks)^2 V(sum^n_tasks Y_i)
#in the simulation all bernoulli draws are random and independent so this simplifies to:
# = (1/n_tasks) V(Y_i)
#as we can see below:

#calc power directly using this variance calculation, using the same parameters and in the above simualtion
n_t = 500
n_c = 500
t_effect = 0.01
base_reject_prob = 0.03
#prior variance whcih doesnt take into account the 1/n-tasks factor
var_t = (base_reject_prob+t_effect)*(1-(base_reject_prob+t_effect))
var_c = base_reject_prob*(1-base_reject_prob)
#now let's add that in:
n_tasks = c(1:40)
var_t_ntasks = var_t*(1/n_tasks)
var_c_ntasks = var_c*(1/n_tasks)
crit_val = qnorm(1 - 0.05)
t_stat_ntasks = (t_effect)/(sqrt(var_t_ntasks/n_t + var_c_ntasks/n_c))
calc_power_n_tasks = 1 - pnorm(crit_val - abs(t_stat_ntasks))

#and we see yes! we can exactly recover the change in simulated variance!
#great so we know what's going on here
ggplot(res) + 
    geom_point(aes(x = n_tasks, y = sim_power, color = "Simulated Var")) + 
    geom_point(aes(x = n_tasks, y = calc_power_n_tasks, color = "Calculated Var")) + 
    theme_bw()

#Do we actually get these gains in power though?
#Next step is to think about whether we realistically gain these precision advantages in our experiment
#unfortunately, the answer is probably not; why?
#because most of of these variance gains are coming from the fact that we assumed all the bernoulli draws/ tasks are INDEPENDENT of each other; in other words all 40 questions a given respondent answers are totally independent
#this is very unlikely to be true in the real world; people's minds do not reset everytime they ahve to evalaute a new pair of signutures, instead waht they said in the previous task or few tasks is very likely to affect what they say in the current and future tasks
#in other words if they just accepted 6 in a row, they might be more likely to reject the 7th simply because they've accepted so many; or maybe since they just rejected a pair of signatures, they might be more likely to reject the next as well or any form of path dependency like this
#next question: how can we deal with this?

#1. in the power analysis, forget the advantage you may slightly gain by having 40 tasks all together; go with the most conservative thing (i think, you should carefully think about this though) of assuming only 1 task, with the most amount of noise
#2. assume some covariance structure. in other words specify what kind of covariance you think the different bernoulli draws/signature match tasks are likely to have, and see how that affects your power; im actually not positive that this covariance structure wont in some way decrease power, i didn't fully work through it

