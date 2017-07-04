
# Simulation estimates ----------------------------------------------------

plot.sims <- function(rel.err.df, col, offset, sim.mod) {
  abline(v=offset, lty=3, col='black')
  summary.stat <- group_by(rel.err.df, est.mod) %>%
    do(quants=quantile(.$rel.err, probs = c(.025, .25, .5, .75, .975)))
  
  mapply(function(x,y) {
    segments(x[1:2] + offset, y, x[5:4] + offset, y, lwd=2,
             col=c('grey60', 
                   stringr::str_replace(gplots::col2hex(col), 'FF', '80')))
    points(x[3] + offset, y, pch=16, col=col)}, 
    summary.stat$quants, 1:3)
}

max1=.01
max2=.4
max3=.75

png('Figs/poisson-sims.png', height=4.8, width=4.8, units = 'in', res=200)
par(mfrow=c(3,1), mar=rep(.5, 4), oma=c(5,5.5,4,1))
for(mod in disp.mods) {
  plot(1,1, xlim=c(-max1, max1+max2+2*max3), ylim=c(.5, 3.5), 
       type='n', axes=FALSE, ann=FALSE)
  text(-.54, 3.7, Hmisc::capitalize(mod), xpd=NA, pos=4)
  text(-.47, 3:1, Hmisc::capitalize(disp.mods), xpd=NA, pos=4)
  if(mod==disp.mods[1]) {
    text(c(0, max1+max2, max1+max2+max3), 3.5, xpd=NA, pos=3, offset=2,
         c('Survival', '% w/in 50m', 'Dist @ 50%'),
         col=c('red', 'blue', 'green'))
  }
  
  sum.stat <- filter(res.df, sim.mod==mod)  %>%
    mutate(rel.err = (surv-survival) / survival) %>%
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
axis(1, at=c(0, .1), labels=c('', .1), line=1)
mtext('Relative error', 1, 3, outer=TRUE)
dev.off()
shell(shQuote('"C:\\Program Files (x86)\\SumatraPDF\\SumatraPDF.exe" Figs/poisson-sims.pdf')) 

system('"C:/Program Files (x86)/SumatraPDF/SumatraPDF.exe" Figs/poisson-sims.pdf')
system2('open', args='Figs/poisson-sims.png') 

# Simulation AIC ----------------------------------------------------------

aic.winner <- apply(aic, c(1,3), which.min) 
winning.pct <- numeric(3)
for(ii in 1:3) {
  winning.pct[ii] <- sum(aic.winner[ii,] == ii)/nreps
}

aic.mns <- apply(aic, 1:2, mean) 
aic.sds <- apply(aic, 1:2, sd)

paste(round(aic.mns, 1), '±', round(aic.sds, 1)) %>% c(winning.pct) %>% 
  matrix(nrow=3, dimnames=list(sim.mod=Hmisc::capitalize(disp.mods), 
                               c(Hmisc::capitalize(disp.mods), 
                                 '% correct model chosen'))) %>%
  write.csv(file='Figs/AIC.csv')


# Data obs  ------------------------------------------------------------
require(ggplot2)

ggplot(report) + facet_wrap(~creek) +
  geom_vline(aes(xintercept=distances), col='gray80') +
  geom_point(aes(x=jitter(distances), y=obscount, col=times), cex=1.5, alpha=.75) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        strip.background = element_blank(), legend.title.align = .5) +
  xlab('Distance from release site (m)') +
  ylab('Number of recaptures') + labs(color = 'Days since\nrelease') +
  scale_color_gradient(high='#132B43', low='#56B1F7', trans='log',
                       breaks=c(1, 5, 25, 200))
