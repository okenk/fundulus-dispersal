require(dplyr)
require(TMB)

# Read in data
read.data <- function(filename) {
  dat <- readLines(filename)
  ls.out <- list()
  
  ls.out$ncreeks <- gsub('\t', '', dat[2]) %>% as.numeric()
  ls.out$nrel <- gsub('\t', '', dat[4])  %>% as.numeric()
  ls.out$nsites <- gsub('\t', '', dat[6]) %>% as.numeric()
  ls.out$nperiods <- gsub('\t', '', dat[8]) %>% as.numeric()
  ls.out$ntraps <- gsub('\t', '', dat[10]) %>% as.numeric()
  
  count.ls <- with(ls.out,
                   strsplit(dat[11 + 1:(ncreeks*nsites*nperiods)], '\t'))
  ls.out$count.mat <- do.call(rbind, count.ls)[,1:ls.out$ntraps] %>% apply(2, as.numeric)
  
  ls.out$distances <- with(ls.out,
                           strsplit(dat[13 + ncreeks*nsites*nperiods], '\t')[[1]] %>% 
                             as.numeric())
  ls.out$times <- with(ls.out,
                       strsplit(dat[15 + ncreeks*nsites*nperiods], '\t')[[1]] %>% 
                         as.numeric())
  return(ls.out)
}

# Set up model
make.inputs <- function(dat.ls, disp.model, count.model, dist.cutoff, sigma.type) {
  disp.model.num <- switch(disp.model,
                           normal = 1,
                           exponential = 2,
                           cauchy = 3)
  if(is.null(disp.model.num)) 
    stop('Please choose "normal", "exponential", or "cauchy" for dispersal model')
  
  count.model.num <- switch(count.model,
                            poisson = 1,
                            neg.binom = 2)
  if(is.null(count.model.num))
    stop('Please choose "poisson" or "neg.binom" for count model')
  
  if(!(sigma.type %in% c('random', 'fixed', 'constant')))
    stop('Please choose "random", "fixed", or "constant" for sigma type')
  
  Data <- with(dat.ls, 
               list(disp_model = disp.model.num,
                    count_model = count.model.num,
                    is_random = as.numeric(sigma.type == 'random'),
                    nrel = nrel,
                    nsites = nsites,
                    nperiods = nperiods,
                    ntraps = ntraps,
                    countmat = count.mat,
                    distances = distances,
                    times = times,
                    dist_cutoff = dist.cutoff))
  
  Map <- list()
  Random <- NULL
  
  if(count.model.num == 1) {
    Map$overdispersion <- factor(NA)
  }
  
  if(sigma.type == 'random') {
    Random <- 'log_sig_disp_eps'
  } else {
    Map$log_sig_disp_sig <- factor(NA)
    
    if(sigma.type == 'fixed') 
      Map$log_sig_disp_mu <- factor(NA) 
    
    if(sigma.type == 'constant') 
      Map$log_sig_disp_eps <- rep(factor(NA), dat.ls$nperiods)
  }
  
  return(list(Data = Data, Map = Map, Random = Random))
}

# Compile and load model
compile('DM_const_sig.cpp')
compile('DM_MM_sig.cpp')

dyn.load(dynlib("DM_const_sig"))
dyn.load(dynlib("DM_MM_sig"))

Parameters <- list(
  survival = .9,
  detectability = .7,
  sig_disp = 20,
  overdispersion = 1
)



report.ls <- list()
for(creek in 1:4) {
  filename <- paste('Data/crk', creek, '.dat', sep='')
  dat <- read.data(filename)
  to.fit <- make.inputs(dat, 'normal', 'poisson', 50, 'random')
  # to.fit$Map$log_sig_disp_mu <- to.fit$Map$log_sig_disp_sig <- factor(NA)
  
  Parameters <- list(
    logit_survival = log(.9/.1),
    log_detectability = log(.7),
    log_sig_disp_mu = log(20),
    log_sig_disp_sig = log(1),
    log_sig_disp_eps = rep(log(1), dat$nperiods), #log(rep(20, 9)),
    overdispersion = 1
  )
  
  model <- MakeADFun(to.fit$Data, Parameters, DLL="DM_MM_sig", map=to.fit$Map,
                     random = to.fit$Random)
  # model <- MakeADFun(to.fit$Data, Parameters, DLL="DM_const_sig")
  #                    map = list(overdispersion = factor(NA)))
  model$env$beSilent()

# Fit model
  Opt = nlminb(start=model$par, objective=model$fn, gradient=model$gr) 
  
  # Opt = nlminb(start=model$par, objective=model$fn, gradient=model$gr, 
  #              lower = c(0, 0, .0001, 0), #, rep(.001, nperiods)), 
  #              upper = c(1, 100, 100000, 100000))
  # 
  # Examine fitted model
  # summary(sdreport(model))
  
  report.ls[[creek]] <- model$report()[[1]] %>% data.frame()
  names(report.ls[[creek]]) <- c('obscount', 'predcount', 'times', 'distances', 
                                 'dist_factor')
  report.ls[[creek]]$creek <- paste('Creek', creek)
}
report <- do.call(rbind, report.ls)



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
