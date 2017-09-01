aic.emp <- array(0, dim=c(4,3,3), dimnames = list(creek=NULL, 
                                              disp.mod=c('normal', 'exponential', 'cauchy'),
                                              mod = c('constant', 'fixed effect', 'asymptote')))
report.ls <- list()

for(creek in 1:4) {
  filename <- paste('Data/crk', creek, '.dat', sep='')
  dat <- read.data(filename)

  Parameters <- list(
    logit_survival = log(.98/.02),
    log_detectability = log(.7),
    log_sig_disp_mu = log(30),
    log_sig_disp_sig = log(1),
    log_sig_disp_eps = rep(log(1), dat$nperiods), #log(rep(20, 9)),
    log_overdispersion = log(2)
  )
  for(disp.mod in c('normal', 'exponential', 'cauchy')) {
    for(mod in 1:3) {
      to.fit <- make.inputs(dat, disp.mod, 'neg.binom',
                            50, switch(mod, 'constant', 'fixed', 'asymptote'))
      if(mod == 1 | mod == 3) {
        Parameters$log_sig_disp_mu <- log(30)
        Parameters$log_sig_disp_eps <- rep(0, dat$nperiods)
      } else if (mod == 2){
        Parameters$log_sig_disp_mu <- 0
        Parameters$log_sig_disp_eps <- rep(log(30), dat$nperiods)
      }
      aic.emp[creek, disp.mod, mod] <- tryCatch({
        model <- MakeADFun(to.fit$Data, Parameters, DLL="DM_MM_sig", map=to.fit$Map,
                           random = to.fit$Random)
        model$env$beSilent()
        Opt <- nlminb(start=model$par, objective=model$fn, gradient=model$gr)
        npar <- 3 + switch(mod, 1, dat$nperiods, 2)
        2*npar + 2*model$fn()
      }, warning = function(w) {
        return(NA)
      })
      
    }
  }
  report.ls[[creek]] <- model$report()[[1]] %>% data.frame()
  names(report.ls[[creek]]) <- c('obscount', 'predcount', 'times', 'distances', 
                                 'dist_factor')
  report.ls[[creek]]$creek <- paste('Creek', creek)
}
report <- do.call(rbind, report.ls)



## asymptote stuff
filename <- paste('Data/crk', 2, '.dat', sep='')
dat <- read.data(filename)

Parameters <- list(
  logit_survival = log(.98/.02),
  log_detectability = log(.7),
  log_sig_disp_mu = log(30),
  log_sig_disp_sig = log(1),
  log_sig_disp_eps = rep(log(1), dat$nperiods), #log(rep(20, 9)),
  log_overdispersion = log(2)
)
to.fit <- make.inputs(dat, 'exponential', 'neg.binom',
                      50, 'asymptote')
model <- MakeADFun(to.fit$Data, Parameters, DLL="DM_MM_sig", map=to.fit$Map,
                   random = to.fit$Random)
model$env$beSilent()
Opt <- nlminb(start=model$par, objective=model$fn, gradient=model$gr)
2*5 + 2*model$fn()


