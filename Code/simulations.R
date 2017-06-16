require(TMB)
set.seed(4309832)
dyn.load(dynlib("DM_const_sig"))

# Data structure
dat <- read.data('Data/crk1.dat')
ncreeks <- dat$ncreeks
nrel <- dat$nrel
nsites <- dat$nsites
nperiods <- dat$nperiods
ntraps <- dat$ntraps
distances <- dat$distances
times <- dat$times

# True values
survival <- 0.98
detectability <- 2.6
sig.disp <- rep(20, nperiods)




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
                            overdispersion = NA, half.distn, ...) {
  ## Take dispersal model choice, other parameters, to simulate observed counts at each of
  ## ntraps traps. Returns vector of length ntraps.
  ## If overdispersion parameter is numeric, simulates negative binomials counts. 
  ## Otherwise (NA, NULL, FALSE, etc.) simulates Poison counts.
  dist.factor <- (1+half.distn) * disp.mod(distance, ...)
  pred.count <- nrel * survival^times * detectability * dist.factor
  if(is.numeric(overdispersion)) {
    counts <- rnbinom(ntraps, mu = pred.count, size = overdispersion)
  } else {
    counts <- rpois(ntraps, pred.count)
  }
  return(counts)
}

Parameters <- list(
  survival = survival,
  detectability = detectability,
  sig_disp = sig.disp[1],
  overdispersion = 1
)

disp.mods <- c('normal', 'exponential', 'cauchy')

for(ii in 1:nreps) {
  # simulate data
  ind <- 1
  for(period in 1:nperiods) {
    for(site in 1:nsites) {
      # half-normal
      count.mat.sim[1,ind,] <- simulate.counts(distance = distances[site], 
                                               times=times[period], disp.mod = dnorm,
                                               ntraps = ntraps, nrel = nrel,
                                               survival = survival, 
                                               detectability = detectability, 
                                               #overdispersion = 2,
                                               half.distn = TRUE, mean = 0, 
                                               sd = sig.disp[period])
      # exponential
      count.mat.sim[2,ind,] <- simulate.counts(distance = distances[site], 
                                               times=times[period], disp.mod = dexp,
                                               ntraps = ntraps, nrel = nrel,
                                               survival = survival, 
                                               detectability = detectability, 
                                               #overdispersion = 2,
                                               half.distn = FALSE,
                                               rate = 1/sig.disp[period])
      # half-cauchy
      count.mat.sim[3,ind,] <- simulate.counts(distance = distances[site], 
                                               times=times[period], disp.mod = dcauchy,
                                               ntraps = ntraps, nrel = nrel,
                                               survival = survival, 
                                               detectability = detectability, 
                                               #overdispersion = 2,
                                               half.distn = TRUE, location = 0, 
                                               scale = sig.disp[period])
      
      ind <- ind + 1
    }
  }
  # fit data
  for(sim.mod in 1:nmods) {
    for(est.mod in 1:nmods) {
      input.ls <- make.inputs(dat.ls = list(nrel = nrel,
                                            nsites = nsites,
                                            nperiods = nperiods,
                                            ntraps = ntraps,
                                            count.mat = count.mat.sim[sim.mod,,],
                                            distances = distances,
                                            times = times), 
                              disp.model = disp.mods[est.mod],
                              count.model = 'poisson', dist.cutoff = 50) 
      
      model <- MakeADFun(data = input.ls$Data, parameters = Parameters, 
                         map = input.ls$Map, 
                         DLL="DM_const_sig")
      model$env$beSilent()
      Opt = nlminb(start=model$par, objective=model$fn, gradient=model$gr,
                   lower = c(0, 0, .0001, 0), #, rep(.001, nperiods)), 
                   upper = c(1, 100000, 100000, 100000))
      fitted.mods[[sim.mod]][[est.mod]][[ii]] <- model
    }
  }
}

true.val.mat <- matrix(c(qnorm(.75, 0, sig.disp[1]), 
                         2*(pnorm(50, 0, sig.disp[1]) - 0.5),
                         qexp(.5, 1/sig.disp[1]), 
                         pexp(50, 1/sig.disp[1]),
                         qcauchy(.75, 0, sig.disp[1]), 
                         2*(pcauchy(50, 0, sig.disp[1]) - 0.5)),
                       nrow=3, byrow=TRUE,
                       dimnames = list(model = disp.mods,
                                       val = c('fifty_pct', 'pct_at_dist')))

res <- rmse <- mare <- temp <- list()
aic <- array(0, dim=c(nmods, nmods, nreps), 
             dimnames=list(sim.mod=disp.mods, est.mod=disp.mods, rep=1:nreps))
for(sim.mod in 1:nmods) {
  res[[sim.mod]] <- rmse[[sim.mod]] <- mare[[sim.mod]] <- list()
  for(est.mod in 1:nmods) {
    sdreports <- sapply(fitted.mods[[sim.mod]][[est.mod]], sdreport)
    res[[sim.mod]][[est.mod]] <- apply(sdreports, 2,
             function(x) c(x['par.fixed']$par.fixed['survival'],
                           x['value']$value)) %>% t()
    err <- t(res[[sim.mod]][[est.mod]]) - c(survival, true.val.mat[sim.mod,])
    rel.err <- err/c(survival, true.val.mat[sim.mod,])
    rmse[[sim.mod]][[est.mod]] <- apply(err, 1, function(x) sqrt(mean(x^2)))
    mare[[sim.mod]][[est.mod]] <- apply(rel.err, 1, function(x) median(abs(x)))
    aic[sim.mod, est.mod,] <- sapply(fitted.mods[[sim.mod]][[est.mod]],
                                     function(x) x$fn())
  }
  
  temp[[sim.mod]] <- do.call(rbind, res[[sim.mod]]) %>% data.frame() %>% 
    mutate(est.mod = rep(disp.mods, each=nreps))
}

res.df <- do.call(rbind, temp) %>% 
  mutate(sim.mod = rep(disp.mods, each=nreps*nmods)) %>% rename(surv = survival)


# When model is correctly specified, estimates are unbiased. 
# Note: This is really a check to make sure we've coded everything correctly. 