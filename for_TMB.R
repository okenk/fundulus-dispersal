## Notes upon leaving: works Olaf's way. Problem with my way is that the integrated probability is numerically 
## equivalent to zero which messes everything up. I don't understand why the integrated probability is zero when
## the density at the midpoint isn't.

require(dplyr)
require(TMB)

compile('DM_dyn_sig.cpp')
dyn.load(dynlib("DM_dyn_sig"))

dat <- readLines('crk1.dat')

ncreeks <- gsub('\t', '', dat[2]) %>% as.numeric()
nrel <- gsub('\t', '', dat[4])  %>% as.numeric()
nsites <- gsub('\t', '', dat[6]) %>% as.numeric()
nperiods <- gsub('\t', '', dat[8]) %>% as.numeric()
ntraps <- gsub('\t', '', dat[10]) %>% as.numeric()

count.ls <- strsplit(dat[11 + 1:(ncreeks*nsites*nperiods)], '\t')
count.mat <- do.call(rbind, count.ls)[,1:ntraps] %>% apply(2, as.numeric)

distances <- strsplit(dat[13 + ncreeks*nsites*nperiods], '\t')[[1]] %>% as.numeric()
times <- strsplit(dat[15 + ncreeks*nsites*nperiods], '\t')[[1]] %>% as.numeric()

Data <- list(# disp_model = 1,
              # ncreeks = ncreeks,
              nrel = nrel,
              nsites = nsites,
              nperiods = nperiods,
              ntraps = ntraps,
              countmat = count.mat,
              distances = distances,
              times = times,
              kappa = 1) # 2 for half normal, 1 for exponential

# Parameters <- list(logit_survival = log(0.99)/log(.01),
#                    logit_detectability = 0,
#                    log_sig_disp = rep(log(5), nperiods))
Parameters <- list(sig_disp_alpha = 1,
                   sig_disp_beta = 1,
                   survival = .9,
                   detectability = .5
                   # sig_disp = rep(5, nperiods)
                   )


model <- MakeADFun(Data, Parameters, DLL="DM_dyn_sig")
# model$report()
model$env$beSilent()
Opt = nlminb(start=model$par, objective=model$fn, gradient=model$gr, 
             lower = c(.001, .000001, 0, .001), #, rep(.001, nperiods)), 
             upper = c(100000000, 100000, 1, 10000))
summary(sdreport(model))

report <- model$report()[[1]]
colnames(report) <- c('obscount', 'predcount', 'distances', 'times')
plot(report[,1], report[,2], xlab = 'obscount', ylab='predcount')
pearson.resid <- (report[,'obscount'] - report[,'predcount'])/sqrt(report[,'predcount'])
plot(report[,'predcount'], pearson.resid)

data.frame(report) %>% group_by(times, distances) %>% 
  summarize(obscount = sum(obscount), check = n()) %>%
  ggplot() + geom_point(aes(x=times, y=distances, col=obscount), cex=3)
