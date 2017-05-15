require(TMB)
set.seed(4309832)

# Data structure
ncreeks <- 1
nrel <- 200
nsites <- 10
nperiods <- 7
ntraps <- 3
# This is one idea that forces sampling sites/events to get farther apart. 
# Could also hard code in values, simulate another way, etc.
distances <- cumsum(c(0, sort(rlnorm(nsites-1, 1.5, 1))))
times <- cumsum(c(1, sort(rlnorm(nperiods-1, 3, 1))))

# True values
survival <- 0.98
detectability <- 0.7 
sig.disp <- rep(20, nperiods)




################# Poisson model
nreps <- 100
nmods <- 3
count.mat.sim <- array(0, dim = c(nmods, nsites*nperiods, ntraps),
                       dimnames = list(mod=NULL, site.period = NULL, trap=NULL))
fitted.mods <- list()
for(sim.mod in 1:nmods) {
  fitted.mods[[sim.mod]] <- list()
  for(est.mod in 1:nmods) {
    fitted.mods[[sim.mod]][[est.mod]] <- list()
  }
}

simulate.counts <- function(distance, times, disp.mod, ntraps, nrel, survival, detectability,
                            half.distn, ...) {
  ## Take dispersal model choice, other parameters, to simulate observed counts at each of
  ## ntraps traps. Returns vector of length ntraps.
  dist.factor <- (1+half.distn) * disp.mod(distance, ...)
  pred.count <- nrel * survival^times * detectability * dist.factor
  rpois(ntraps, pred.count)
}

Parameters <- list(
  survival = .9,
  detectability = .7,
  sig_disp = 10
)

for(ii in 1:nreps) {
  # simulate data
  # could resimulate data structure here if desired
  ind <- 1
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
  for(sim.mod in 1:nmods) {
    for(est.mod in 1:nmods) {
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
      Opt = nlminb(start=model$par, objective=model$fn, gradient=model$gr, 
                   lower = c(0, 0, .0001), #, rep(.001, nperiods)), 
                   upper = c(1, 1, 100000))
      fitted.mods[[sim.mod]][[est.mod]][[ii]] <- sdreport(model)
    }
  }
}

res <- list()
for(mod in 1:nmods) {
  res[[mod]] <- sapply(fitted.mods[[mod]][[mod]], function(x) x$value) %>%
    apply(1, function(x) c(mn=mean(x), sd=sd(x)))
}
res
# When model is correctly specified, estimates are unbiased. 
# Note: This is really a check to make sure we've coded everything correctly. 