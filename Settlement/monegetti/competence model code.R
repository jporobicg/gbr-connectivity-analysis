                                        # R code for the competence Weibull-exponential model fitting for:
                                        # High-frequency sampling and piecewise models reshape dispersal kernels of
                                        # a common reef coral

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

###########################################################################

rm(list=ls())

                                        # Obtain competence data
#setwd(("/Volumes/KIMI/COTS/COTS_GBRMPAWorking_group/MFBCRCode/doi_10/")) # Location
setwd('/home/por07g/Documents/Projects/GBR_modeling/Connectivity/R_codes')
raw.data = read.csv("A.tenuisGBR2012metamorphosis.csv", header = TRUE)

                                        # Metamorphosis experiment age data (days)
age = raw.data$age..d.

                                        # Metamorphosis experiment number of metamorphosed larvae data
no.meta = raw.data$meta

                                        # Metamorphosis experiment number of larvae remaining swimming data
no.swimming = raw.data$swimming

weibexp <- function(pars){
                                        # The weibexp function calculates the negative log-likelihood of a set of
                                        # model parameters given the empirical metamorphosis data.

                                        # COMPETENCE LOG LIKELIHOOD

                                        # Larvae either remain swimming or metamorphose.
                                        # The likelihood of an individual metamorphosing is described by a
                                        # piecewise model. The piecewise model is comprised of two functions
                                        # characterizing the larval competence dynamics: one describing the early
                                        # portion of the larval life (tc<t<Tcp, pr.early), and a second
                                        # characterizing the latter portion (t>Tcp, pr.latter).

                                        # Therefore, the likelihood of an individual metamorphosing follows:
                                        #     pr.settle = t<tc, 0; tc<t<Tcp, pr.early; t>Tcp, pr.latter
                                        # The likelihood of an individual swimming is then:
                                        #     pr.swim = 1 - pr.settle
                                        # Thus the full likelihood is:
                                        #     no.meta*log(pr.settle) + no.swimming*log(pr.swim)
                                        # summed over all observations and replicates

                                        # The rate of acquisition of competence, when t>tc
    a = exp(pars[1])

                                        # The first loss of competence parameter for the Weibull model
                                        # characterizing the early portion of the larval life
    b1 = exp(pars[2])

                                        # The second loss of competence parameter for the Weibull model
                                        # characterizing the early portion of the larval life
    v1 = exp(pars[3])

                                        # The loss of competence parameter for the exponential model
                                        # characterizing the latter portion of the larval life
    b2 = exp(pars[4])

                                        # The pre-competent period
    tc = exp(pars[5])

                                        # The estimated change point time
    Tcp = exp(pars[6])

                                        # The log-likelihood for replicates and observations
    loglikeli = c()

                                        # Define the integrals that need to be solved to calculate pr.early
                                        # and pr.latter (according to the Weibull-exponential model)
    int.early = function(tau)
    {a*exp(-a*(tau-tc))*exp(-((b1*(time-tau))^v1))}
    int.latter.1 = function(tau)
    {a*exp(-a*(tau-tc))*exp(-((b1*(Tcp-tau))^v1))*exp(-b2*(time-Tcp))}
    int.latter.2 = function(tau)
    {a*exp(-a*(tau-tc))*exp(-b2*(time-tau))}

    for(t in 1:length(age)){
        time = age[t]

        pr.early = ifelse(time < tc, 0,
                   ifelse(time <= Tcp,
                          integrate(int.early, lower = tc,
                                    upper = time)$value, 0))
                                        # pr.early is the likelihood that an individual alive at time t is
                                        # competent to settle at time t, when t<Tcp

        pr.latter = ifelse(time <= Tcp, 0,
                           integrate(int.latter.1, lower = tc,
                                     upper = Tcp)$value
                           + integrate(int.latter.2, lower = Tcp,
                                       upper = time)$value)
                                        # pr.latter is the likelihood that an individual alive at time t is
                                        # competent to settle at time, when t>Tcp

        pr.settle = ifelse(time <= Tcp, pr.early, pr.latter)
                                        # pr.settle is the likelihood that an individual alive at time t is
                                        # competent to settle at time t

        settle = ifelse(time < tc, 0, no.meta[t]*log(pr.settle))
                                        # settle is the log-likelihood that an individual alive at time t is
                                        # competent to settle at time t

        pr.swim = 1-pr.settle
                                        # pr.swim is the likelihood that an individual alive at time t is
                                        # still swimming at time t

        swim = ifelse(pr.swim == 0, 0, no.swimming[t]*log(pr.swim))
                                        # swim is the log-likelihood that an individual alive at time t is
                                        # still swimming at time t

        loglikeli[t] = settle + swim
                                        # loglikeli is a vector of the competence log-likelihoods,
                                        # where element t corresponds to the log-likelihood for the tth
                                        # observation for a particular replicate (the first n elements
                                        # correspond to the n observations for replicate 1, the next n
                                        # elements correspond to the n observations for replicate 2 and
                                        # so on).
    }

    competence.negLL = -sum(loglikeli)
                                        # competence.negLL is the full competence negative log-likelihood
}


                                        # Fit the Weibull-exponenial model
loginit.pars = c(log(1.292), log(0.001878), log(0.3645), log(0.3969),
                 log(3.332), log(69.91))
                                        # loginit.pars are initial parameter estimates for the numerical
                                        # optimization procedure (the values given are the maximum likelihood
                                        # estimates)
thisfit = optim(par = loginit.pars, weibexp, method = "L-BFGS-B",
                lower = c(-Inf, -Inf, -Inf, -Inf, log(2.604), log(3.729)),
                upper = c(+Inf, +Inf, +Inf, +Inf, log(3.729), log(79.77)),
                control = list(maxit = 500))
                                        # The "L-BFGS-B" methods was used to give the parameters box constraints.
                                        # This was required to constrain:
                                        # - The time delay tc to a minimum of log(day previous to the day on
                                        #   which competent larvae were first observed in all replicates) and a
                                        #   maximum of log(day on which competent larvae were were first observed
                                        #   in all replicates). Required as if tc is greater than the time that
                                        #   settlers were observed then settlement would be occurring when the
                                        #   model predicted competence probability of zero, causing for NaNs to
                                        #   be produced.
                                        # - The estimated change point Tcp to a minimum of log(day on which
                                        #   competent larvae were first observed in all replicates)
                                        #   (i.e., Tcp>tc) and a maximum of log(day of last observation).
                                        #   Required to satisfy the condition that tc<Tcp.

                                        # Parameter estimates
a = exp(thisfit$par[1])
b1 = exp(thisfit$par[2])
v1 = exp(thisfit$par[3])
b2 = exp(thisfit$par[4])
tc = exp(thisfit$par[5])
Tcp = exp(thisfit$par[6])

                                        # AIC calculation
neg.MLL = thisfit$value
AIC = 2*6+2*neg.MLL

                                        # Plot the competence dynamics
sim.time = seq(0, 80, 1)

int.early = function(tau)
{a*exp(-a*(tau-tc))*exp(-((b1*(time-tau))^v1))}
int.latter.1 = function(tau)
{a*exp(-a*(tau-tc))*exp(-((b1*(Tcp-tau))^v1))*exp(-b2*(time-Tcp))}
int.latter.2 = function(tau)
{a*exp(-a*(tau-tc))*exp(-b2*(time-tau))}

pred.settle.we = c()

for(t in 1:length(sim.time)){
    time = sim.time[t]

    pr.early = ifelse(time < tc, 0,
               ifelse(time <= Tcp,
                      integrate(int.early, lower = tc,
                                upper = time)$value, 0))

    pr.latter = ifelse(time <= Tcp, 0,
                       integrate(int.latter.1, lower = tc,
                                 upper = Tcp)$value
                       + integrate(int.latter.2, lower = Tcp,
                                   upper = time)$value)

    pred.settle.we[t] = ifelse(time <= Tcp, pr.early, pr.latter)
}

print(pred.settle.we)
x11()
#windows()
matplot(age, no.meta/20, pch = 16, col = 'black',
        ylab = "Proportion competent (square-root scale)",
        xlab = "Cohort age (days)",
        main = "Weibull-exponential Model")
axis(2, at = c(0, sqrt(0.05), sqrt(0.2), sqrt(0.6), 1),
     labels = c(0, 0.05, 0.2, 0.6, 1))
lines(sim.time, pred.settle.we, col = 'red', lwd = 2)
