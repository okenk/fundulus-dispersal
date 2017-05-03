set.seed(4309832)

# True values
survival <- 0.98
detectability <- 0.05 
sig.disp <- rlnorm(nperiods, 2, 1) 
# Note sig.disp can get smaller over time. This is how it was being estimated. 
# But that makes it very hard to estimate sig.disp when there are few recaptures,
# as there are for the later sampling periods

# Data structure
ncreeks <- 1
nrel <- 200
nsites <- 10
nperiods <- 7
ntraps <- 3
kappa <- 2 # 2 for half-normal, 1 for exponential
# This is one idea that forces sampling sites/events to get farther apart. 
# Could also hard code in values, simulate another way, etc.
distances <- cumsum(c(0, sort(rlnorm(nsites-1, 1.5, 1))))
times <- cumsum(c(1, sort(rlnorm(nperiods-1, 3, 1))))





################### Multinomial model
Derived quantities
hazards <- matrix(0, nrow=nperiods, ncol=nsites, dimnames = list(period=NULL, site=NULL))
for(period in 1:nperiods) {
  hazards[period,] <- detectability * exp(-(distances/sig.disp[period])^kappa)
}
pcapture <- 1-exp(-apply(hazards, 1, sum))

# Random data
nreps <- 1
ncaught <- matrix(0, nrow=nreps, ncol=nperiods, dimnames = list(rep=NULL, period=NULL))
ncaught.bytrap <- array(0, dim = c(nreps, nperiods, nsites), 
                        dimnames = list(rep=NULL, period=NULL, site=NULL))
# Simulate survival explicitly?
for(period in 1:nperiods) {
  ncaught[,period] <- rbinom(nreps, nrel, (survival^times[period]) * pcapture[period])
  ncaught.bytrap[,period,] <- rmultinom(nreps, ncaught[,period], hazards[period,])
}



eval.fun <- function(func, ...) {
  func(...)
}
eval.fun(func=dnorm, x=2, mean=0, sd=1)

################# Poisson model
nreps <- 1
nmods <- 3
count.mat.sim <- array(0, dim = c(nmods, nsites*nperiods, ntraps),
                       dimnames = list(mod=NULL, site.period = NULL, trap=NULL))
fitted.mods <- list()

simulate.counts <- function(distance, times, disp.mod, ntraps, nrel, survival, detectability,
                            half.distn, ...) {
  ## Take dispersal model choice, other parameters, to simulate observed counts at each of
  ## ntraps traps. Returns vector of length ntraps.
  dist.factor <- (1+half.distn) * disp.mod(distance, ...)
  pred.count <- nrel * survival^times * detectability * dist.factor
  rpois(ntraps, pred.count)
}

Parameters <- list(# sig_disp_alpha = 1,
                   # sig_disp_beta = 1,
                   survival = .9,
                   detectability = .5,
                   sig_disp = 5
  )

ind <- 1
for(ii in 1:nreps) {
  # simulate data
  # could resimulate data structure here if desired
  for(period in 1:nperiods) {
    for(site in 1:nsites) {
      # half-normal
      count.mat.sim[1,ind,] <- simulate.counts(distance = distances[site], 
                                                  times=times[period], disp.mod = dnorm,
                                                  ntraps = ntraps, nrel = nrel,
                                                  survival = survival, 
                                                  detectability = detectability, 
                                                  half.distn = TRUE, mean = 0, 
                                                  sd = sig.disp[period])
      # exponential
      count.mat.sim[2,ind,] <- simulate.counts(distance = distances[site], 
                                                  times=times[period], disp.mod = dexp,
                                                  ntraps = ntraps, nrel = nrel,
                                                  survival = survival, 
                                                  detectability = detectability, 
                                                  half.distn = FALSE,
                                                  rate = 1/sig.disp[period])
      # half-cauchy
      count.mat.sim[3,ind,] <- simulate.counts(distance = distances[site], 
                                                  times=times[period], disp.mod = dcauchy,
                                                  ntraps = ntraps, nrel = nrel,
                                                  survival = survival, 
                                                  detectability = detectability, 
                                                  half.distn = TRUE, location = 0, 
                                                  scale = sig.disp[period])
      
      ind <- ind + 1
    }
  }
  # fit data
  for(sim.mod in 1:nmod) {
    for(est.mod in 1:nmod) {
      Data <- list(disp_model = est.mod,
        # ncreeks = ncreeks,
        nrel = nrel,
        nsites = nsites,
        nperiods = nperiods,
        ntraps = ntraps,
        countmat = count.mat.sim[sim.mod,,],
        distances = distances,
        times = times)
      
      model <- MakeADFun(Data, Parameters, DLL="DM_dyn_sig")
      model$env$beSilent()
      opt <- nlminb(start=model$par, objective=model$fn, gradient=model$gr, 
                    lower = c(.001, .000001, 0, .001), #, rep(.001, nperiods)), 
                    upper = c(100000000, 100000, 1, 10000))
      fitted.mods[[sim.mod]][[est.mod]][[ii]] <- sdreport(model)
    }
  }
}
