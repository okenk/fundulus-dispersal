## Notes upon leaving: works Olaf's way. Problem with my way is that the integrated probability is numerically 
## equivalent to zero which messes everything up. I don't understand why the integrated probability is zero when
## the density at the midpoint isn't.

require(dplyr)
require(TMB)

compile('DM_dyn_sig.cpp')
dyn.load(dynlib("DM_dyn_sig"))

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

Data <- list(# disp_model = 1,
              # ncreeks = ncreeks,
              nrel = nrel,
              nsites = nsites,
              nperiods = nperiods,
              ntraps = ntraps,
              countmat = count.mat,
              distances = distances,
              times = times,
              kappa = 2) # 2 for half normal, 1 for exponential

# Parameters <- list(logit_survival = log(0.99)/log(.01),
#                    logit_detectability = 0,
#                    log_sig_disp = rep(log(5), nperiods))
Parameters <- list(survival = .5,
                   detectability = .5,
                   sig_disp = rep(5, 1))


model <- MakeADFun(Data, Parameters, DLL="DM_dyn_sig")
model$report()
model$env$beSilent()
Opt = nlminb(start=model$par, objective=model$fn, gradient=model$gr, lower = c(0,0,rep(0.001, 1)), 
             upper = c(1, 1, rep(10000, 1)))
summary(sdreport(model))

report <- model$report()[[1]]
colnames(report) <- c('obscount', 'predcount', 'times', 'distances')
plot(report[,1], report[,2], xlab = 'obscount', ylab='predcount')
pearson.resid <- (report[,'obscount'] - report[,'predcount'])/sqrt(report[,'predcount'])
plot(report[,'predcount'], pearson.resid)
