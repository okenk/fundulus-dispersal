require(ggplot2)
load('Code/sim-results.RData')

# Simulation estimates ----------------------------------------------------

# plot.sims <- function(rel.err.df, col, offset) {
#   dark.col <- stringr::str_replace(gplots::col2hex(col), 'FF', '99')
#   abline(v=offset, lty=3, col='black')
#   summary.stat.temp <- group_by(rel.err.df, est.mod) %>%
#     do(quants=quantile(.$rel.err, probs = c(.025, .25, .5, .75, .975), na.rm=TRUE))
#   summary.stat <- summary.stat.temp[match(c(disp.mods, 'no dispersal'), summary.stat.temp$est.mod),] %>%
#     filter(!is.na(est.mod))
#     
#   mapply(function(x,y) {
#     segments(x[1:2] + offset, y, x[5:4] + offset, y, lwd=2,
#              col=c('grey80', dark.col))
#     points(x[3] + offset, y, pch=21, col=dark.col, bg=col)}, 
#     summary.stat$quants, length(summary.stat$quants):1)
# }
# 
# max1=.3
# max2=.35
# max3=.75
# 
# # big.df <- res.df
# # res.df <- filter(res.df, est.mod != 'no dispersal')
# 
# png('Figs/poisson-sims.png', height=5.5, width=4.8, units = 'in', res=200)
# par(mfrow=c(n_distinct(res.df$sim.mod),1), mar=rep(.5, 4), oma=c(5,6,4,1))
# plot.max <- n_distinct(res.df$est.mod) + 0.5
# for(mod in disp.mods) {
#   plot(1,1, xlim=c(-max1, max1+max2+2*max3), ylim=c(.5, plot.max), 
#        type='n', axes=FALSE, ann=FALSE)
#   text(-.7, plot.max+0.2, Hmisc::capitalize(mod), xpd=NA, pos=4)
#   text(-.63, (plot.max-0.5):1, Hmisc::capitalize(c(disp.mods, 'no dispersal')), xpd=NA, pos=4)
#   if(mod==disp.mods[1]) {
#     text(c(0, max1+max2, max1+max2+max3), plot.max, xpd=NA, pos=3, offset=2,
#          c('T max', '% w/in 50m', 'Dist @ 50%'),
#          col=stringr::str_replace(gplots::col2hex(c('red', 'blue', 'green')), 'FF', '99'))
#   }
#   
#   temp <- filter(res.df, sim.mod==mod)  %>%
#     mutate(rel.err = (t_max-t.max) / t.max) %>%
#     plot.sims(col='red', offset=0)
#     
#   temp <- filter(res.df, sim.mod==mod)  %>%
#     mutate(rel.err = (pct_at_dist-true.val.mat[sim.mod,'pct_at_dist']) /
#              true.val.mat[sim.mod,'pct_at_dist']) %>%
#    plot.sims(col='blue', offset=max1+max2)
# 
#   temp <- filter(res.df, sim.mod==mod)  %>%
#     mutate(rel.err = (fifty_pct-true.val.mat[sim.mod,'fifty_pct']) /
#              true.val.mat[sim.mod,'fifty_pct']) %>%
#     plot.sims(col='green', offset=max1+max2+max3)
# }
# box('inner')
# axis(1, at=max1+max2+max3+c(0, .25, .5), labels = c('', 0.25, 0.5), line = 1)
# axis(1, at=max1+max2+c(0, .2, .4), labels=c('', .2, .4), line=1)
# axis(1, at=c(0, .2), labels=c('', .2), line=1)
# mtext('Relative error', 1, 3, outer=TRUE)
# dev.off()
# 
# system('open Figs/poisson-sims.png')


# Survival estimates ----------------------------------------------------
t.max.true <- t.max

                                                                       
png('Figs/survival-results.png', height=8, width=8, units = 'in', res=200)
sim.res.df %>%
  mutate(rel.err = (t.max-t.max.true)/t.max.true) %>%
  group_by(sim.mod, est.mod, disp.str) %>%
  summarize(min=quantile(rel.err, .05), low=quantile(rel.err, .25), mid=median(rel.err),
            high=quantile(rel.err, .75), 
            max=ifelse(quantile(rel.err, .95)<.5, quantile(rel.err, .95), .5)) %>%
  ggplot() +
  geom_point(aes(x=mid, y=disp.str), cex=2) +
  geom_segment(aes(x=min, xend=max, y=disp.str, yend=disp.str)) +
  geom_segment(aes(x=low, xend=high, y=disp.str, yend=disp.str), lwd=1.5) +
  geom_vline(xintercept=0, lty=2) +
  facet_grid(est.mod ~ sim.mod, scales = 'free_y', space='free_y', labeller = label_both) +
  ggsidekick::theme_sleek() +
  xlab(expression(paste('Relative error in ', T[max]))) + ylab('Dispersal structure')
dev.off()
system('open Figs/survival-results.png')

sim.res.df %>%
  mutate(rel.err = (t.max-t.max.true)/t.max.true) %>%
  group_by(sim.mod, est.mod, disp.str) %>%
  summarize(min=quantile(rel.err, .05), low=quantile(rel.err, .25), mid=median(rel.err),
            high=quantile(rel.err, .75), 
            max=ifelse(quantile(rel.err, .95)<.5, quantile(rel.err, .95), .5)) %>%
  filter(est.mod != 'no dispersal', disp.str != 'constant') %>%
  arrange(sim.mod, desc(abs(mid))) %>%
  group_by(disp.str) %>% summarize(mean(mid))

# This seems to be sensitive to the choice of tmax function. Propose basing bias 
# conclusions on estimates of M.


# dispersal estimates ----------------------------------------------------
sim.res.df$fifty.pct.true <- sapply(sim.res.df$sim.mod, function(x)
  switch(x,
         normal = qnorm(.75, 0, 32.5),
         exponential = qexp(.5, 1/32.5),
         cauchy = qcauchy(.75, 32.5)))

png('Figs/dispersal-results.png', height=8, width=8, units = 'in', res=200)
sim.res.df %>%
  filter(disp.str != 'fixed', est.mod != 'no dispersal') %>%
  mutate(rel.err = (fifty.pct - fifty.pct.true)/fifty.pct.true) %>%
  group_by(sim.mod, est.mod, disp.str) %>%
  summarize(min=quantile(rel.err, .05), low=quantile(rel.err, .25), mid=median(rel.err),
            high=quantile(rel.err, .75), 
            max=quantile(rel.err, .95)) %>% #ifelse(quantile(rel.err, .95)<.5, quantile(rel.err, .95), .5)) %>%
  ggplot() +
  geom_point(aes(x=mid, y=disp.str), cex=2) +
  geom_segment(aes(x=min, xend=max, y=disp.str, yend=disp.str)) +
  geom_segment(aes(x=low, xend=high, y=disp.str, yend=disp.str), lwd=1.5) +
  geom_vline(xintercept=0, lty=2) +
  facet_grid(est.mod ~ sim.mod, scales = 'free_y', space='free_y', labeller = label_both) +
  ggsidekick::theme_sleek() +
  xlab('Relative error in distance at 50% dispersal, relative to asymptote') +
  ylab('Dispersal structure')
dev.off()
system('open Figs/dispersal-results.png')


# Simulation AIC ----------------------------------------------------------

aic.winner <- apply(aic.sim[,1:3,], c(1,3), which.min) 
winning.pct <- numeric(3)
for(ii in 1:3) {
  winning.pct[ii] <- sum(aic.winner[ii,] == ii)/nreps
}

aic.mns <- apply(aic.sim[,1:3,], 1:2, mean, na.rm=TRUE) 
aic.sds <- apply(aic.sim[,1:3,], 1:2, sd)

paste(round(aic.mns, 1), '±', round(aic.sds, 1)) %>% c(paste0(winning.pct*100, '%')) %>% 
  matrix(nrow=3, dimnames=list(sim.mod=Hmisc::capitalize(disp.mods), 
                               c(Hmisc::capitalize(disp.mods), 
                                 '% correct model chosen'))) %>%
  write.csv(file='Figs/AIC.csv')


# Format results table ----------------------------------------------------

m.ind <- which(names(sdreport(emp.mod.fits[[1]][[1]][[1]])$value)=='annual_mort')
t.max.ind <- which(names(sdreport(emp.mod.fits[[1]][[1]][[1]])$value)=='t_max')
fifty.pct.ind <- which(names(sdreport(emp.mod.fits[[1]][[1]][[1]])$value)=='fifty_pct')
# Those indices should be unaffected by sampling structure for creek or 
# dispersal structure for model

table.width <- 5
print.mat <- matrix('', nrow=22, ncol=2 * table.width + 1)
print.mat[1,1] <- 'Random effect dispersal rate'
print.mat[1,table.width+2] <- 'Asymptotic dispersal rate'
row.marker <- 2
for(disp.mod in c('normal', 'exponential', 'cauchy')) {
  row.marker <- row.marker + 1
  print.mat[row.marker, c(1,table.width+2)] <- Hmisc::capitalize(disp.mod)
  
  row.marker <- row.marker + 1
  print.mat[row.marker, c(1,table.width+2)] <- 'Location'
  print.mat[row.marker, c(2,table.width+3)] <- 'AIC'
  print.mat[row.marker, c(3,table.width+4)] <- 'M (1/yr)'
  print.mat[row.marker, c(4,table.width+5)] <- 'tmax (yrs)'
  print.mat[row.marker, c(5)] <- 'Mean dist @ 50% disp (m)'
  print.mat[row.marker, c(table.width+6)] <- 'Asymptotic dist @ 50% disp (m)'
  
  for(creek in 1:4) {
    row.marker <- row.marker + 1
    print.mat[row.marker, c(1, table.width+2)] <- paste('Creek', creek)
    print.mat[row.marker, c(2, table.width+3)] <- aic.emp[creek, disp.mod, 
                                                          c('random effect', 
                                                            'asymptote')] %>%
      round(1)
    
    print.mat[row.marker, c(3, table.width+4)] <- paste0(
      round(sapply(emp.mod.fits[[creek]][[disp.mod]][3:4], 
                   function(x) sdreport(x)$value[m.ind]), 2),
      '±',
      round(sapply(emp.mod.fits[[creek]][[disp.mod]][3:4], 
                   function(x) sdreport(x)$sd[m.ind]), 2))
    
    print.mat[row.marker, c(4, table.width+5)] <- paste0(
      round(sapply(emp.mod.fits[[creek]][[disp.mod]][3:4], 
                   function(x) sdreport(x)$value[t.max.ind]), 2),
      '±',
      round(sapply(emp.mod.fits[[creek]][[disp.mod]][3:4], 
                   function(x) sdreport(x)$sd[t.max.ind]), 2))
    
    print.mat[row.marker, c(5, table.width+6)] <- paste0(
      round(sapply(emp.mod.fits[[creek]][[disp.mod]][3:4], 
                   function(x) sdreport(x)$value[fifty.pct.ind]), 1),
      '±',
      round(sapply(emp.mod.fits[[creek]][[disp.mod]][3:4], 
                   function(x) sdreport(x)$sd[fifty.pct.ind]), 1))
  }
  row.marker <- row.marker+1
}

# Remove results for models that did not converge
temp <- as.vector(which(is.na(print.mat), arr.ind=TRUE)) + matrix(c(rep(0, table.width-2), 1:(table.width-2)),
                                                                  nrow=2, byrow=TRUE)
print.mat[temp[1,], temp[2,]] <- ''
print.mat[is.na(print.mat)] <- 'Did not converge'

write.table(print.mat, file='Figs/model-fits.csv', row.names = FALSE, col.names = FALSE,
            sep=',')


# Plot data and fits from best model --------------------------------------

require(ggplot2)

best.mods <- apply(aic.emp[,,3:4], 1, function(x) which(x==min(x, na.rm=TRUE), arr.ind=TRUE))
# best.mods[2,2] <- -1
report.ls <- list()

for(creek in 1:4) {
  report.ls[[creek]] <- emp.mod.fits[[creek]][[best.mods[1,creek]]][[best.mods[2,creek]+2]]$report()[[1]] %>% 
    data.frame()
  names(report.ls[[creek]]) <- c('obscount', 'predcount', 'times', 'distances', 
                                 'dist_factor')
  report.ls[[creek]]$creek <- paste('Creek', creek)
}

report <- do.call(rbind, report.ls)

png('Figs/obs-pred-data.png', height=5, width=7, units = 'in', res=200)
report %>%
  group_by(creek, times, distances) %>%
  summarize(obscount=mean(obscount), predcount=first(predcount)) %>%
  ggplot() +
  geom_col(aes(x=factor(times), y=obscount, fill=factor(distances)), position='dodge') +
  geom_point(aes(x=factor(times), y=predcount, bg=factor(distances)), col='black', 
             position=position_dodge(width=0.9), show.legend = FALSE) +
  facet_wrap(~creek, scales = 'free') +
  xlab('Time since release (days)') +
  ylab('Number of recaptures') +
  scale_fill_discrete(name='Distance from\nrelease (m)') +
  ggsidekick::theme_sleek()
dev.off()
system('open Figs/obs-pred-data.png')

pdf('obs-pred-data-talk.pdf')

ggplot() + theme_classic()

report %>%
  filter(creek=='Creek 2', times==2) %>%
  # group_by(times, distances) %>%
  group_by(distances) %>%
  summarize(obscount=mean(obscount), predcount=first(predcount), 
            times=factor(2, levels = unique(report$times[report$creek=='Creek 2']))) %>%
  ggplot() +
  geom_col(aes(x=factor(times), y=obscount, fill=factor(distances)), position='dodge') +
  # geom_point(aes(x=factor(times), y=predcount, bg=factor(distances)), col='black', 
  #            position=position_dodge(width=0.9), show.legend = FALSE) +
  xlab('Time since release (days)') +
  ylab('Number of recaptures') +
  scale_fill_discrete(name='Distance from\nrelease (m)', drop=FALSE) +
  scale_x_discrete(limits=factor(unique(report$times[report$creek=='Creek 2'])), drop=FALSE) +
  ylim(obscount=1.05*c(0,max(report$obscount))) +
  theme_classic(base_size = 18)

report %>%
  filter(creek=='Creek 2') %>%
  group_by(times, distances) %>%
  # group_by(distances) %>%
  summarize(obscount=mean(obscount), predcount=first(predcount)) %>%
  ggplot() +
  geom_col(aes(x=factor(times), y=obscount, fill=factor(distances)), position='dodge') +
  # geom_point(aes(x=factor(times), y=predcount, bg=factor(distances)), col='black', 
  #            position=position_dodge(width=0.9), show.legend = FALSE) +
  xlab('Time since release (days)') +
  ylab('Number of recaptures') +
  scale_fill_discrete(name='Distance from\nrelease (m)', drop=FALSE) +
  ylim(obscount=1.05*c(0,max(report$obscount))) +
  theme_classic(base_size = 18)
dev.off()


