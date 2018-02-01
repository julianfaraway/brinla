library(brinla)

# Additive Models

library(INLA)
library(mgcv)
set.seed(2) 
dat <- gamSim(1, n = 400, dist = "normal", scale = 2)
str(dat)
formula <- y ~ f(x0, model = 'rw2', scale.model = TRUE) + f(x1, model = 'rw2', scale.model = TRUE) + f(x2, model = 'rw2', scale.model = TRUE) + f(x3, model = 'rw2', scale.model = TRUE)
n.group <- 50
x0.new <- inla.group(dat$x0, n = n.group, method = 'quantile')
x1.new <- inla.group(dat$x1, n = n.group, method = 'quantile')
x2.new <- inla.group(dat$x2, n = n.group, method = 'quantile')
x3.new <- inla.group(dat$x3, n = n.group, method = 'quantile')
dat.inla <- list(y = dat$y, x0 = x0.new, x1 = x1.new, x2 = x2.new, x3 = x3.new)
result <- inla(formula, data = dat.inla)
round(result$summary.fixed, 4)
library(brinla)
bri.band.plot(result, name = 'x0', type = 'random', xlab='', ylab='')
lines(sort(dat$x0), (dat$f0 - mean(dat$f0))[order(dat$x0)], lty = 2)
 round(bri.hyperpar.summary(result), 4)
library(INLA)
data(Munich, package = "brinla")
g <- system.file("demodata/munich.graph", package = "INLA")
formula <- rent ~ f(location, model = "besag", graph = g, scale.model = TRUE) + f(year, model = "rw2", values = seq(1918, 2001), scale.model = TRUE) + f(floor.size, model = "rw2", values = seq(17, 185), scale.model = TRUE) + Gute.Wohnlage + Beste.Wohnlage + Keine.Wwv + Keine.Zh + Kein.Badkach  + Besond.Bad + Gehobene.Kueche + zim1 + zim2 + zim3 + zim4 + zim5 + zim6 - 1
result <- inla(formula, data = Munich, control.predictor = list(compute = TRUE)) 
round(result$summary.fixed, 3)
library(brinla)
bri.band.ggplot(result, name = 'floor.size', type = 'random')
bri.band.ggplot(result, name = 'year', type = 'random')
map.munich(result$summary.random$location$mean)
map.munich(result$summary.random$location$sd)
yhat <- result$summary.fitted.values$mean
residual <- Munich$rent - yhat
plot(yhat, residual, ylab = 'Residual', xlab = 'Fitted value')
abline(0,0)
plot(yhat, Munich$rent, ylab = 'Rent', xlab = 'Fitted value')
abline(0,1)

# Generalized Additive Models

data(kyphosis, package = 'brinla')
str(kyphosis)
library(INLA)
formula1 <- Kyphosis ~ 1 + f(Age, model = 'rw2') + f(StartVert, model = 'rw2') + f(NumVert, model = 'rw2')
result1 <- inla(formula1, family='binomial', data = kyphosis, control.predictor = list(compute = TRUE), control.compute = list(waic = TRUE))	
bri.band.ggplot(result1, name = 'Age', type = 'random') 
bri.band.ggplot(result1, name = 'StartVert', type = 'random')
bri.band.ggplot(result1, name = 'NumVert', type = 'random')
eta <- result1$summary.linear.predictor$mean
phat <- result1$summary.fitted.values$mean
phat.lb <- result1$summary.fitted.values$'0.025quant'
phat.ub <- result1$summary.fitted.values$'0.975quant'
data.plot <- data.frame(eta, phat, phat.lb, phat.ub)
ggplot(data.plot, aes(y = phat, x = eta)) + geom_errorbar(aes(ymin = phat.lb, ymax = phat.ub), width = 0.2, col = 'gray') + geom_point() + theme_bw(base_size = 20) + labs(x = 'Linear predictor', y = 'Probability')
kyphosis$AgeSq <- (kyphosis$Age)^2
formula2 <- Kyphosis ~ 1 + StartVert + NumVert + Age + AgeSq
result2 <- inla(formula2, family='binomial', data = kyphosis, control.predictor = list(compute = TRUE), control.compute = list(waic = TRUE))	
round(result2$summary.fixed, 4)
c(result1$waic$waic, result2$waic$waic)
data(mack, package = 'gamair')
str(mack, vec.len = 2)
loc.obs <- cbind(mack$lon, mack$lat)
plot(loc.obs, cex = 0.2+mack$egg.count/50, cex.axis = 1.5)
plot(loc.obs, cex = 0.2+mack$egg.dens/150, cex.axis = 1.5)
library(INLA)
mesh <- inla.mesh.2d(loc.obs, cutoff = 0.05, max.edge = c(.5,1))
tps <- bri.tps.prior(mesh)
node <- mesh$idx$loc
formula <- egg.count ~ -1 + salinity + c.dist + f(temp.20m, model = 'rw2') + f(node, model = tps, diagonal = 1e-6)
result <- inla(formula, family = 'poisson', data = mack, offset = log(net.area))
round(result$summary.fixed, 4)
bri.band.ggplot(result, name = 'temp.20m', type = 'random')
data(mackp, package = 'gamair') 							
proj <- inla.mesh.projector(mesh, loc = cbind(mackp$lon, mackp$lat))
spa.mean <- inla.mesh.project(proj, result$summary.random$node$mean)
library(fields)
quilt.plot(mackp$lon, mackp$lat, spa.mean, nx = length(unique(mackp$lon)), ny = length(unique(mackp$lat)))
n.pre <- dim(mackp)[1]
y.pre <- c(mack$egg.count, rep(NA, n.pre))
z1.pre <- c(mack$salinity, mackp$salinity)  
z2.pre <- c(mack$c.dist, mackp$c.dist)  
x.pre <- inla.group(c(mack$temp.20m, mackp$temp.20m), n = 100)  
E.pre <- c(mack$net.area, rep(0.25^2, n.pre))  
loc.pre <- rbind(cbind(mack$lon,mack$lat), cbind(mackp$lon,mackp$lat))
mesh2 <- inla.mesh.2d(loc.pre, cutoff = 0.05, max.edge = c(.5, 1))
node2 <- mesh2$idx$loc
tps2 <- bri.tps.prior(mesh2)
mack.pre <- list(egg.count = y.pre, salinity = z1.pre, c.dist = z2.pre, temp.20m = x.pre, node = node2)
formula <- egg.count ~ -1 + salinity + c.dist + f(temp.20m, model = 'rw2') + f(node, model = tps2, diagonal = 1e-6)
link <- rep(NA, length(y.pre))
link[which(is.na(y.pre))] <- 1
result2 <- inla(formula, family = 'poisson', data = mack.pre, offset = log(E.pre), control.predictor = list(link = link, compute = TRUE), control.compute = list(config = TRUE))			 
idx.pre <- which(is.na(y.pre))
res.pre <- result2$summary.fitted.values[idx.pre, ]
library(brinla)
res.exc <- excursions.brinla(result2, name = 'Predictor', ind = idx.pre, u = log(10), type = '>', alpha = 0.05, method = 'NIQC')
quilt.plot(mackp$lon, mackp$lat, res.exc$F, nx = length(unique(mackp$lon)), ny = length(unique(mackp$lat)))
res.exc$E

# Generalized Additive Mixed Models

data(sole, package = 'gamair')
str(sole)
solr <- sole
solr$off <- log(sole$a.1 - sole$a.0)
solr$a <- (sole$a.1 + sole$a.0)/2
solr$t <- solr$t - mean(sole$t)
solr$t <- solr$t/var(sole$t)^0.5
solr$la <- solr$la - mean(sole$la)
solr$lo <- solr$lo - mean(sole$lo)
solr$station <- as.numeric(factor(with(solr, paste(-la, -lo, -t, sep=""))))
solr$eggs <- solr$eggs*1000
solr$off <- solr$off + log(1000)
length(which(solr$eggs == 0))/length(solr$eggs)
library(INLA)
loc <- cbind(solr$lo, solr$la)
mesh <- inla.mesh.2d(loc, max.edge = c(0.1, 0.2))
node <- mesh$idx$loc
library(brinla)
tps <- bri.tps.prior(mesh, constr = TRUE)
formula1 <- eggs ~ a + I(la*t) + I(la*t^2) + I(lo*t) + I(lo*t^2) + f(t, a, model = 'rw2', scale.model = TRUE) + f(node, model = tps, diagonal = 1e-6) + f(station, model = 'iid')
formula0 <- eggs ~ a + I(la*t) + I(la*t^2) + I(lo*t) + I(lo*t^2) + f(t, a, model = 'rw2', scale.model = TRUE) + f(node, model = tps, diagonal = 1e-6)
result0 <- inla(formula0, family = 'zeroinflatednbinomial0', offset = off, data = solr, control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE))
result1 <- inla(formula1, family = 'zeroinflatednbinomial0', offset = off, data = solr, control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE))
result2 <- inla(formula1, family = 'zeroinflatedpoisson0', offset = off, data = solr, control.predictor = list(compute = TRUE), control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE))
c(result0$dic$dic, result1$dic$dic, result2$dic$dic)
c(result0$waic$waic, result1$waic$waic, result2$waic$waic)
result <- result1
round(result$summary.fixed, 4)
bri.band.ggplot(result, name = 't', type = 'random')
proj <- inla.mesh.projector(mesh)
spa.mean <- inla.mesh.project(proj, result$summary.random$node$mean)
library(fields)
image.plot(proj$x, proj$y, spa.mean, xlim=c(-1.5, 2), ylim=c(-1, 1))
points(loc, pch = 19, cex = 0.2)
tmp <- bri.hyperpar.summary(result)
row.names(tmp) <- c("Overdispersion", "Zero-probability", 'SD for t', 'Theta1 for node', 'SD for station')
round(tmp, 4)
length(which(result$cpo$failure > 0))
improved.result <- inla.cpo(result)
idx <- which(solr$eggs != 0)
dat.plot <- data.frame(x = 1:length(idx), y = log(improved.result$cpo$cpo[idx]))
ggplot(dat.plot, aes(x = x, y = y)) + geom_point() 
round(solr[which.min(improved.result$cpo$cpo),], 4)
solr[solr$station==326, c('eggs', 'stage', 'a', 't')]
sessionInfo()
