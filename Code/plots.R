
# Simulation estimates ----------------------------------------------------

plot.sims <- function(rel.err.df, col, offset) {
  dark.col <- stringr::str_replace(gplots::col2hex(col), 'FF', '99')
  abline(v=offset, lty=3, col='black')
  summary.stat.temp <- group_by(rel.err.df, est.mod) %>%
    do(quants=quantile(.$rel.err, probs = c(.025, .25, .5, .75, .975), na.rm=TRUE))
  summary.stat <- summary.stat.temp[match(c(disp.mods, 'no dispersal'), summary.stat.temp$est.mod),] %>%
    filter(!is.na(est.mod))
    
  mapply(function(x,y) {
    segments(x[1:2] + offset, y, x[5:4] + offset, y, lwd=2,
             col=c('grey80', dark.col))
    points(x[3] + offset, y, pch=21, col=dark.col, bg=col)}, 
    summary.stat$quants, length(summary.stat$quants):1)
}

max1=.3
max2=.35
max3=.75

# big.df <- res.df
# res.df <- filter(res.df, est.mod != 'no dispersal')

png('Figs/poisson-sims.png', height=5.5, width=4.8, units = 'in', res=200)
par(mfrow=c(n_distinct(res.df$sim.mod),1), mar=rep(.5, 4), oma=c(5,6,4,1))
plot.max <- n_distinct(res.df$est.mod) + 0.5
for(mod in disp.mods) {
  plot(1,1, xlim=c(-max1, max1+max2+2*max3), ylim=c(.5, plot.max), 
       type='n', axes=FALSE, ann=FALSE)
  text(-.7, plot.max+0.2, Hmisc::capitalize(mod), xpd=NA, pos=4)
  text(-.63, (plot.max-0.5):1, Hmisc::capitalize(c(disp.mods, 'no dispersal')), xpd=NA, pos=4)
  if(mod==disp.mods[1]) {
    text(c(0, max1+max2, max1+max2+max3), plot.max, xpd=NA, pos=3, offset=2,
         c('T max', '% w/in 50m', 'Dist @ 50%'),
         col=stringr::str_replace(gplots::col2hex(c('red', 'blue', 'green')), 'FF', '99'))
  }
  
  temp <- filter(res.df, sim.mod==mod)  %>%
    mutate(rel.err = (t_max-t.max) / t.max) %>%
    plot.sims(col='red', offset=0)
    
  temp <- filter(res.df, sim.mod==mod)  %>%
    mutate(rel.err = (pct_at_dist-true.val.mat[sim.mod,'pct_at_dist']) /
             true.val.mat[sim.mod,'pct_at_dist']) %>%
   plot.sims(col='blue', offset=max1+max2)

  temp <- filter(res.df, sim.mod==mod)  %>%
    mutate(rel.err = (fifty_pct-true.val.mat[sim.mod,'fifty_pct']) /
             true.val.mat[sim.mod,'fifty_pct']) %>%
    plot.sims(col='green', offset=max1+max2+max3)
}
box('inner')
axis(1, at=max1+max2+max3+c(0, .25, .5), labels = c('', 0.25, 0.5), line = 1)
axis(1, at=max1+max2+c(0, .2, .4), labels=c('', .2, .4), line=1)
axis(1, at=c(0, .2), labels=c('', .2), line=1)
mtext('Relative error', 1, 3, outer=TRUE)
dev.off()

system('open Figs/poisson-sims.png')


# Dispersal structures ----------------------------------------------------
par(mfrow=c(n_distinct(surv.res.df$sim.mod),1), mar=rep(.5, 4), oma=c(5,6,4,1))
plot.max <- n_distinct(surv.res.df$est.mod) + 0.5
for(mod in disp.mods) {
  plot(1,1, xlim=c(-max1, max1+max2+2*max3), ylim=c(.5, plot.max), 
       type='n', axes=FALSE, ann=FALSE)
  text(-.7, plot.max+0.2, Hmisc::capitalize(mod), xpd=NA, pos=4)
  # text(-.63, (plot.max-0.5):1, Hmisc::capitalize(c(disp.mods, 'no dispersal')), xpd=NA, pos=4)
  # if(mod==disp.mods[1]) {
  #   text(c(0, max1+max2, max1+max2+max3), plot.max, xpd=NA, pos=3, offset=2,
  #        c('T max', '% w/in 50m', 'Dist @ 50%'),
  #        col=stringr::str_replace(gplots::col2hex(c('red', 'blue', 'green')), 'FF', '99'))
  # }
  
  temp <- filter(surv.res.df, sim.mod==mod)  %>%
    mutate(rel.err = (asympBH-t.max) / t.max) %>%
    plot.sims(col='red', offset=0)
  
  temp <- filter(surv.res.df, sim.mod==mod)  %>%
    mutate(rel.err = (asympVB-t.max) / t.max) %>%
   plot.sims(col='blue', offset=max1+max2)

  temp <- filter(surv.res.df, sim.mod==mod)  %>%
    mutate(rel.err = (random-t.max) / t.max) %>%
    plot.sims(col='green', offset=max1+max2+max3)
} 
box('inner')
axis(1, at=max1+max2+max3+c(0, .25, .5), labels = c('', 0.25, 0.5), line = 1)
axis(1, at=max1+max2+c(0, .2, .4), labels=c('', .2, .4), line=1)
axis(1, at=c(0, .2), labels=c('', .2), line=1)
mtext('Relative error', 1, 3, outer=TRUE)



# Simulation AIC ----------------------------------------------------------

aic.winner <- apply(aic.sim, c(1,3), which.min, na.rm=TRUE) 
winning.pct <- numeric(3)
for(ii in 1:3) {
  winning.pct[ii] <- sum(aic.winner[ii,] == ii)/nreps
}

aic.mns <- apply(aic.sim, 1:2, mean, na.rm=TRUE) 
aic.sds <- apply(aic.sim, 1:2, sd)

paste(round(aic.mns, 1), '±', round(aic.sds, 1)) %>% c(winning.pct) %>% 
  matrix(nrow=3, dimnames=list(sim.mod=Hmisc::capitalize(disp.mods), 
                               c(Hmisc::capitalize(disp.mods), 
                                 '% correct model chosen'))) %>%
  write.csv(file='Figs/AIC.csv')


# Data obs  ------------------------------------------------------------
require(ggplot2)

png('Figs/obs-data.png', width=7, height=7, units='in', res=200)
ggplot(report) + facet_wrap(~creek) +
  geom_vline(aes(xintercept=distances), col='gray80') +
  geom_point(aes(x=jitter(distances), y=obscount, col=times), cex=1.5, alpha=.5) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(), legend.title.align = .5) +
  xlab('Distance from release site (m)') +
  ylab('Number of recaptures') + labs(color = 'Days since\nrelease') +
  scale_color_gradient(high='#132B43', low='#56B1F7', trans='log',
                       breaks=c(1, 5, 25, 200))
dev.off()
system('open Figs/obs-data.png')
