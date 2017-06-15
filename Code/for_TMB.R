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
make.inputs <- function(dat.ls, disp.model, count.model, dist.cutoff) {
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
  
  Data <- with(dat.ls, 
               list(disp_model = disp.model.num,
                    count_model = count.model.num,
                    nrel = nrel,
                    nsites = nsites,
                    nperiods = nperiods,
                    ntraps = ntraps,
                    countmat = count.mat,
                    distances = distances,
                    times = times,
                    dist_cutoff = dist.cutoff))
  
  Map <- list()
  if(count.model.num == 1) {
    Map$overdispersion <- factor(NA)
  }
  
  return(list(Data = Data, Map = Map))
}

# Compile and load model
compile('DM_const_sig.cpp')
# compile('DM_MM_sig.cpp')

dyn.load(dynlib("DM_const_sig"))
# dyn.load(dynlib("DM_MM_sig"))

Parameters <- list(
  survival = .9,
  detectability = .7,
  sig_disp = 20,
  overdispersion = 1
)
# Parameters <- list(
#   survival = .9,
#   detectability = .7,
#   sig_disp_alpha = 300,
#   sig_disp_beta = 20,
#   overdispersion = 1
# )

for(creek in 1:4) {
  filename <- paste('Data/crk', creek, '.dat', sep='')
  dat <- read.data(filename)
  to.fit <- make.inputs(dat, 'normal', 'neg.binom', 50)
model <- MakeADFun(to.fit$Data, Parameters, DLL="DM_const_sig", map=to.fit$Map) 
model$env$beSilent()

# Fit model
Opt = nlminb(start=model$par, objective=model$fn, gradient=model$gr, 
             lower = c(0, 0, 0, 0, 0), #, rep(.001, nperiods)), 
             upper = c(1, 100, 10000000, 10000000, 10000000))
 
Opt = nlminb(start=model$par, objective=model$fn, gradient=model$gr, 
             lower = c(0, 0, .0001, 0), #, rep(.001, nperiods)), 
             upper = c(1, 100, 100000, 100000))

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
