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

# Derived quantities
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
for(period in 1:nperiods) {
  ncaught[,period] <- rbinom(nreps, nrel, (survival^times[period]) * pcapture[period])
  ncaught.bytrap[,period,] <- rmultinom(nreps, ncaught[,period], hazards[period,])
}
