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

Data <- list(disp_model = 3,
              # ncreeks = ncreeks,
              nrel = nrel,
              nsites = nsites,
              nperiods = nperiods,
              ntraps = ntraps,
              countmat = count.mat,
              distances = distances,
              times = times,
              full_site_width = 5)

Parameters <- list(logit_survival = log(0.99)/log(.01),
                   logit_detectability = 0,
                   log_sig_disp = rep(log(1), nperiods))


model <- MakeADFun(Data, Parameters, DLL="DM_dyn_sig")
model$env$beSilent()
Opt = nlminb(start=model$par, objective=model$fn, gradient=model$gr)
summary(sdreport(model))

report <- model$report()[[1]]
colnames(report) <- c('obscount', 'predcount', 'times', 'distances')
plot(report[,1], report[,2], xlab = 'obscount', ylab='predcount')
pearson.resid <- (report[,'obscount'] - report[,'predcount'])/sqrt(report[,'predcount'])
plot(report[,'predcount'], pearson.resid)


detectability         0.05140536 1.061541e-02
sig_disp             37.22104860 1.082432e+01
sig_disp            149.79671750 4.744857e+01
sig_disp             22.79334599 6.114543e+00
sig_disp             55.92583823 1.556798e+01
sig_disp            118.06119073 4.224284e+01
sig_disp            222.40320588 1.723668e+02
sig_disp             13.00019325 1.042115e+01
sig_disp              2.99624710 7.775782e+00
sig_disp             29.26261636 2.643645e+01
