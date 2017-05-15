require(dplyr)
require(TMB)

# Compile and load model
compile('DM_dyn_sig.cpp')
dyn.load(dynlib("DM_dyn_sig"))

# Read in data
dat <- readLines('crk2.dat')

ncreeks <- gsub('\t', '', dat[2]) %>% as.numeric()
nrel <- gsub('\t', '', dat[4])  %>% as.numeric()
nsites <- gsub('\t', '', dat[6]) %>% as.numeric()
nperiods <- gsub('\t', '', dat[8]) %>% as.numeric()
ntraps <- gsub('\t', '', dat[10]) %>% as.numeric()

count.ls <- strsplit(dat[11 + 1:(ncreeks*nsites*nperiods)], '\t')
count.mat <- do.call(rbind, count.ls)[,1:ntraps] %>% apply(2, as.numeric)

distances <- strsplit(dat[13 + ncreeks*nsites*nperiods], '\t')[[1]] %>% as.numeric()
times <- strsplit(dat[15 + ncreeks*nsites*nperiods], '\t')[[1]] %>% as.numeric()

# Set up model
Data <- list(disp_model = 3,
              # ncreeks = ncreeks,
              nrel = nrel,
              nsites = nsites,
              nperiods = nperiods,
              ntraps = ntraps,
              countmat = count.mat,
              distances = distances,
              times = times)

# Parameters <- list(logit_survival = log(0.99)/log(.01),
#                    logit_detectability = 0,
#                    log_sig_disp = rep(log(5), nperiods))
Parameters <- list(
                   survival = .9,
                   detectability = .7,
                   sig_disp = 10
                   )


model <- MakeADFun(Data, Parameters, DLL="DM_dyn_sig")
model$env$beSilent()

# Fit model
Opt = nlminb(start=model$par, objective=model$fn, gradient=model$gr, 
             lower = c(0, 0, .0001), #, rep(.001, nperiods)), 
             upper = c(1, 1, 100000))

# Examine fitted model
summary(sdreport(model))

report <- model$report()[[1]]
colnames(report) <- c('obscount', 'predcount', 'times', 'distances', 'dist_factor')
plot(report[,1], report[,2], xlab = 'obscount', ylab='predcount')
pearson.resid <- (report[,'obscount'] - report[,'predcount'])/sqrt(report[,'predcount'])
plot(report[,'predcount'], pearson.resid)
# Residuals skew positive, i.e., counts are under-predicted. 
# This is not the case with simulated data.

# Examine data
require(ggplot2)
data.frame(report) %>% group_by(times, distances) %>% 
  summarize(obscount = sum(obscount), check = n()) %>%
  ggplot() + geom_point(aes(x=log(times), y=distances, col=obscount), cex=4)
