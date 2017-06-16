plot.sims <- function(rel.err.df, col, offset, sim.mod) {
  abline(v=offset, lty=3, col='black')
  summary.stat <- group_by(rel.err.df, est.mod) %>%
    do(quants=quantile(.$rel.err, probs = c(.025, .25, .5, .75, .975)))
  
  mapply(function(x,y) {
    segments(x[1:2] + offset, y, x[5:4] + offset, y,
             col=c('grey60', col))
    points(x[3] + offset, y, pch=16, col=col)}, 
    summary.stat$quants, 1:3)
}

max1=.01
max2=.4
max3=.75

pdf('Figs/poisson-sims.pdf', height=4.8, width=4.8)
par(mfrow=c(3,1), mar=rep(.5, 4), oma=c(5,5.5,4,1))
for(mod in disp.mods) {
  plot(1,1, xlim=c(-max1, max1+max2+2*max3), ylim=c(.5, 3.5), 
       type='n', axes=FALSE, ann=FALSE)
  text(-.54, 3.7, mod, xpd=NA, pos=4)
  text(-.47, 3:1, disp.mods, xpd=NA, pos=4)
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
mtext('Relative error', 1, 2, outer=TRUE)
dev.off()

aic.winner <- apply(aic, c(1,3), which.min) 
winning.pct <- numeric(3)
for(ii in 1:3) {
  winning.pct[ii] <- sum(aic.winner[ii,] == ii)/nreps
}

aic.mns <- apply(aic, 1:2, mean) 
aic.sds <- apply(aic, 1:2, sd)

pdf('Figs/aic.pdf')
par(mar=c(5,4,7,1))
barplot(t(aic.mns), beside = TRUE, las=1,
        col=c('grey30', 'grey50', 'grey80'))
mtext('AIC', 2, 3)
mtext('Simulation model', 1, 3)
legend(9, 200, disp.mods, fill=c('grey30', 'grey50', 'grey80'), xpd=NA, bty='n', 
       title='Estimation model')
xpos <- .5+c(1:3, 5:7, 9:11)
segments(rep(xpos, 2), rep(as.vector(t(aic.mns)), 2),
       rep(xpos, 2), c(as.vector(t(aic.mns + aic.sds)), 
                       as.vector(t(aic.mns - aic.sds))), xpd=NA)
dev.off()