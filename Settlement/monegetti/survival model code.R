# R code for the survival Weibull-Weibull model fitting for: High-frequency
# sampling and piecewise models reshape dispersal kernels of a common reef coral

# R-3.4.1 for Windows

# Joanne Moneghetti
# College of Science and Engineering
# James Cook University, Australia
# Email: joanne.moneghetti@jcu.edu.au

# Joana Figueiredo
# Halmos College of Natural Sciences and Oceanography
# Nova Southeastern University, USA

# Andrew Baird
# ARC Centre of Excellence for Coral Reef Studies
# James Cook University Australia

# Sean Connolly
# College of Science and Engineering
# ARC Centre of Excellence for Coral Reef Studies
# James Cook University, Australia

################################################################################

rm(list=ls())

# Obtain survival data
setwd("/Volumes/KIMI/COTS/COTS_GBRMPAWorking_group/MFBCRCode/doi_10/") #Location of data
raw.data = read.csv("TenuisGBR2012longtermsurvival.csv", header = TRUE)

# Survival experiment age data (days)
age = unique(raw.data$age..d.)

# Survival experiment survival data (number of larvae alive)
no.surv = raw.data$surv

# Number of time observations
no.obs = length(age)

# Number of intervals
no.intervals = no.obs-1

# Number of replicates per observation time
no.reps = length(no.surv)/no.obs

# Transformed survival data
# surv.rep is a matrix of the survival data, where each column corresponds
# to a replicate and each row corresponds to an observation
surv.rep = matrix(no.surv, nrow=no.obs, ncol=no.reps)

weibweib = function(pars){
  # The weibweib function calculates the negative log-likelihood of a set 
  # of model parameters given the empirical metamorphosis data.
  
  # SURVIVAL LOG LIKELIHOOD
  
  # Larvae either remain alive from time t to time t+1 or die.
  # The likelihood of an individual surviving the interval is described by
  # a piecewise model. The piecewise model is comprised of two functions
  # characterizing the larval survival dynamics: one describing the early
  # portion of the larval life (t<Tcp), and a second characterizing the
  # latter portion (t>Tcp).
  
  # The likelihood that an individual that is alive at time t is still
  # alive at time t+1 follows:
  #     pr.alive
  # The likelihood of an individual that is alive at time t dying before
  # time t+1 is then:
  #     1 - pr.alive
  # Thus the full likelihood is:
  #     (no. of larvae dying between time t and time t+1)*log(1 - pr.alive)
  #           + (no. of larvae alive at time t still alive
  #               at time t+1)*log(pr.alive)
  # summed over all intervals and replicates
  
  # The first mortality parameter for the Weibull model characterizing the
  # early portion of the larval life
  u1 = exp(pars[1])
  
  # The second mortality parameter for the Weibull model characterizing
  # the early portion of the larval life
  v1 = exp(pars[2])
  
  # The first mortality parameter for the Weibull model characterizing the
  # latter portion of the larval life
  u2 = exp(pars[3])
  
  # The second mortality parameter for the Weibull model characterizing
  # latter portion of the larval life
  v2 = exp(pars[4])
  
  # The estimated change point time
  Tcp = exp(pars[5])
  
  # The log-likelihood for intervals
  int.loglikeli = c()
  
  # The log-likelihood for replicates
  rep.loglikeli = c()
  
  for(r in 1:no.reps){
    
    survived = surv.rep[,r]
      # survived is the survival data (number of larvae alive) for
      # replicate r
    
    for(t in 1:no.intervals){
      
      pr.alive = ifelse(age[t+1] < Tcp,
                        exp(-(u1*age[t+1])^v1)/exp(-(u1*age[t])^v1),
                        ifelse(age[t] < Tcp,
                               (exp(-(u1*Tcp)^v1)/exp(-(u1*age[t])^v1))*
                                 (exp(-(u2*age[t+1])^v2)/exp(-(u2*Tcp)^v2)),
                               exp(-(u2*age[t+1])^v2)/exp(-(u2*age[t])^v2)))
        # pr.alive is the likelihood that an individual alive at time t
        # is still alive at time t+1
      
      int.loglikeli[t] = (survived[t]-survived[t+1])*log(1-pr.alive) +
        (survived[t+1]*log(pr.alive))
        # int.loglikeli is a vector of the surivival log-likelihoods,
        # where element t corresponds to the log-likelihood for the tth
        # interval
    }
    
    rep.loglikeli[r] = sum(int.loglikeli)
      # rep.loglikeli is a vector of the sum of the survival
      # log-likelihoods for each replicate, where element r corresponds
      # to the rth replicate
  }
  
  survival.negLL = -sum(rep.loglikeli)
    # survival.negLL is the full survival negative log-likelihood
}

# Fit the Weibull-Weibull model
loginit.pars = c(log(0.400), log(2.892), log(0.019), log(1.716),
                 log(2.583))
  # loginit.pars are initial parameter estimates for the numerical
  # optimization procedure (the values given are the maximum likelihood
  # estimates)
thisfit = optim(loginit.pars, weibweib)

# Parameter estimates
u1 = exp(thisfit$par[1])
v1 = exp(thisfit$par[2])
u2 = exp(thisfit$par[3])
v2 = exp(thisfit$par[4])
Tcp = exp(thisfit$par[5])

# AIC calculation
neg.MLL = thisfit$value
AIC = 2*5+2*neg.MLL

# Plot the survival dynamics
sim.time = seq(0, 94, 1)

pred.alive.ww = c()

for(t in 1:length(age)){
  time = sim.time[t]
  pred.alive.ww[t] = ifelse(time <= Tcp, exp(-(u1*time)^v1),
                            exp(-(u1*Tcp)^v1)*exp(-(u2*(time-Tcp))^v2))
}

windows()
matplot(age, sqrt(surv.rep/50), pch = 20, col = "black",
        ylab = "Proportion alive (square-root scale)", yaxt = "n",
        xlab = "Cohort age (days)", main = "Weibull-Weibull Model")
axis(2, at = c(0, sqrt(0.05), sqrt(0.2), sqrt(0.6), 1),
     labels = c(0, 0.05, 0.2, 0.6, 1))
lines(age, sqrt(pred.alive.ww), col = "red", lwd = 2)

