library(brinla)

# Introduction


# Smoothing Splines

library(INLA); library(brinla)
set.seed(1)
n <- 100
x <- seq(0, 1,, n)
f.true <- (sin(2*pi*x^3))^3
y <- f.true + rnorm(n, sd = 0.2)
data.inla <- list(y = y, x = x)
formula1 <- y ~ -1 + f(x, model = "rw1", constr = FALSE)
result1 <- inla(formula1, data = data.inla)
formula2 <- y ~ -1 + f(x, model = "rw2", constr = FALSE)
result2 <- inla(formula2, data = data.inla)
round(head(result1$summary.random$x), 4)
fhat <- result1$summary.random$x$mean  ## posterior mean
f.lb <- result1$summary.random$x$'0.025quant'  ## 2.5% percentile
f.ub <- result1$summary.random$x$'0.975quant'  ## 97.5% percentile
library(ggplot2)
data.plot <- data.frame(y = y, x = x, f.true = f.true, fhat = fhat, f.lb = f.lb, f.ub = f.ub)
ggplot(data.plot, aes(x = x, y = y)) + geom_line(aes(y = fhat)) + geom_line(aes(y = f.true), linetype = 2) + geom_ribbon(aes(ymin = f.lb, ymax = f.ub), alpha = 0.2) + geom_point(aes(y = y)) + theme_bw(base_size = 20)
result1$summary.hyperpar
result2$summary.hyperpar
round(bri.hyperpar.summary(result1), 4)
round(bri.hyperpar.summary(result2), 4)
fit.ss <- smooth.spline(x, y)
res <- (fit.ss$yin - fit.ss$y)/(1 - fit.ss$lev)  
fhat <- fit.ss$y  ## fitted curve
f.lb <- fhat - 2*sd(res)*sqrt(fit.ss$lev)  ## lower bound
f.ub <- fhat + 2*sd(res)*sqrt(fit.ss$lev)  ## upper bound  
(fit.ss$lambda)
result2$summary.hyperpar$mean[2]/result2$summary.hyperpar$mean[1]
library(mgcv)
fit.gam <- gam(y ~ s(x))
res.gam <- predict(fit.gam, se.fit = TRUE)
fhat <- res.gam$fit  ## fitted curve
f.lb <- res.gam$fit - 2*res.gam$se.fit  ## lower bound
f.ub <- res.gam$fit + 2*res.gam$se.fit  ## upper bound
a1 <- 5e-5
b1 <- 5e-5
lgprior1 <- list(prec = list(param = c(a1, b1)))
a2 <- -0.5
b2 <- 5e-5
lgprior2 <- list(prec = list(param = c(a2, b2)))
formula <- y ~ -1 + f(x, model = "rw2", constr = FALSE, hyper = lgprior2)
result <- inla(formula, data = data.inla, control.family = list(hyper = lgprior1))
set.seed(1)
n <- 100
t <- 0.1
x <- seq(0, t,, n)
f.true <- (sin(2*pi*(x/t)^3))^3
y <- f.true + rnorm(n, sd = 0.2)
data.inla <- list(y = y, x = x)
formula1 <- y ~ -1 + f(x, model = "rw2", constr = FALSE)
result1 <- inla(formula1, data = data.inla)
p <- bri.band.ggplot(result1, name = 'x', type = 'random') 
p + geom_line(aes(y = f.true), linetype = 2)
round(bri.hyperpar.summary(result1), 4)
formula2 <- y ~ -1 + f(x, model = "rw2", constr = FALSE, scale.model = TRUE)
result2 <- inla(formula2, data = data.inla)
round(bri.hyperpar.summary(result2), 4)
formula1 <- rent ~ -1 + f(floor.size, model = 'rw2', constr = FALSE)
formula2 <- rent ~ -1 + f(year, model = 'rw2', constr = FALSE)
data(Munich, package = "brinla")
result1 <- inla(formula1, data = Munich)
result2 <- inla(formula2, data = Munich)
bri.band.plot(result1, name = 'floor.size', alpha = 0.05, xlab = 'Floor size', ylab = 'Rent', type = 'random')
points(Munich$floor.size, Munich$rent, pch = 20, cex = 0.2)
bri.band.plot(result2, name = 'year', alpha = 0.05, xlab = 'Year', ylab = 'Rent', type = 'random')
points(Munich$year, Munich$rent, pch = 20, cex = 0.2)
round(bri.hyperpar.summary(result1), 4)
round(bri.hyperpar.summary(result2), 4)
formula3 <- rent ~ -1 + f(floor.size, model = "rw2", values = seq(17, 185), constr = FALSE)
formula4 <- rent ~ -1 + f(year, model = "rw2", values = seq(1918, 2001), constr = FALSE)
result3 <- inla(formula3, data = Munich)
result4 <- inla(formula4, data = Munich)
x.new <- c(1925, 1938, 1945)
xx <- c(Munich$year, x.new)
yy <- c(Munich$rent, rep(NA, length(x.new)))
data.pred <- list(y = yy, x = xx)
formula5 <- y ~ -1 + f(x, model = 'rw2', constr = FALSE)
result5 <- inla(formula5, data = data.pred, control.predictor = list(compute = TRUE))
ID <- result5$summary.random$x$ID
idx.new <- sapply(x.new, function(x) which(ID==x))
round(result5$summary.random$x[idx.new,], 4)
nsamp <- 10000
pred.marg <- NULL
for(i in 1:length(x.new)){
  ##Sample error precision
  error.prec <- inla.hyperpar.sample(nsamp, result5)[,1] 
  ##Sample new errors
  new.eps <- rnorm(nsamp, mean = 0, sd = 1/sqrt(error.prec))
  ##Sample linear predictor
  pm.new <- result5$marginals.linear.predictor[[which(is.na(data.pred$y))[i]]]
  ##Combine samples
  samp <- inla.rmarginal(nsamp, pm.new) + new.eps 
  pred.marg <- cbind(samp, pred.marg)
}
p.mean <- colMeans(pred.marg)  #mean
p.sd <- apply(pred.marg, 2, sd)  #standard deviation
p.quant <- apply(pred.marg, 2, function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))  #quantiles
data.frame(ID = x.new, mean = p.mean, sd = p.sd, '0.025quant' = p.quant[1,], '0.5quant' = p.quant[2,], '0.975quant' = p.quant[3,], check.names = FALSE)

# Thin-plate Splines\label{sec:tps}

test.fun <- function(x,z,sig.x,sig.z){.75/(pi*sig.x*sig.z)*exp(-(x - .2)^2/sig.x^2-(z - .3)^2/sig.z^2) + .45/(pi*sig.x*sig.z)*exp(-(x - .7)^2/sig.x^2-(z - .8)^2/sig.z^2)}
nrow <- 100  ## number of rows
ncol <- 100  ## number of columns
s.mat <- matrix(NA, nrow = nrow, ncol = ncol)
for(i in 1:nrow){
  for(j in 1:ncol){
    s.mat[i,j] <- test.fun(i/100, j/100, 0.3, 0.4)
  }
}
set.seed(1)
noise.mat <- matrix(rnorm(nrow*ncol, sd = 0.3), nrow, ncol)
y.mat <- s.mat + noise.mat
y <- inla.matrix2vector(y.mat)
formula <- y ~ -1 + f(x, model="rw2d", nrow=nrow, ncol=ncol, constr=F)
data <- data.frame(y = y, x = 1:(nrow*ncol))
result <- inla(formula, data = data)
persp(s.mat, theta = 25, phi = 30, expand = 0.8, xlab='', ylab='', zlab='', ticktype = 'detailed') 
fhat <- result$summary.random$x$mean
fhat.mat <- inla.vector2matrix(fhat, nrow, ncol)
persp((fhat.mat - s.mat)^2, theta = 25, phi = 30, expand = 0.8, xlab = '', ylab = '', zlab = '', ticktype = 'detailed') 
round(bri.hyperpar.summary(result), 4)
data(SPDEtoy)
str(SPDEtoy)
library(fields)
quilt.plot(SPDEtoy$s1, SPDEtoy$s2, SPDEtoy$y)
coords <- as.matrix(SPDEtoy[,1:2])  ## coordinates of data
mesh <- inla.mesh.2d(loc=coords, max.edge=c(0.15, 0.2), cutoff=0.02)
plot(mesh, main='')
tps <- bri.tps.prior(mesh, theta.mean = 0, theta.prec = 0.001)
formula <- y ~ -1 + f(x, model = tps, diagonal = 1e-6)
data.inla <- list(y = SPDEtoy$y, x = mesh$idx$loc)
result <- inla(formula, data = data.inla, control.predictor = list(compute = TRUE))
fhat <- result$summary.random$x$mean
quilt.plot(mesh$loc[,1:2], fhat)
yhat <- result$summary.fitted$mean
plot(SPDEtoy$y, yhat, ylab='Fitted values', xlab='Observed response')
abline(0,1)
round(result$summary.hyperpar[,1:5], 3)
loc.pre <- rbind(c(0.1, 0.1), c(0.5, 0.55), c(0.7, 0.9))
y.pre <- rep(NA, dim(loc.pre)[1])
y2 <- c(y.pre, SPDEtoy$y)
coords2 <- rbind(loc.pre, coords)
mesh2 <- inla.mesh.2d(coords2, max.edge = c(0.15, 0.2), cutoff = 0.02)
tps2 <- bri.tps.prior(mesh2)
formula <- y ~ -1 + f(x, model = tps2)
data2.inla <- list(y = y2, x = mesh2$idx$loc)
result2 <- inla(formula, data = data2.inla, control.predictor = list(compute = TRUE))
(idx.pre <- which(is.na(y2)))
round(result2$summary.fitted[idx.pre,], 3)
pm.samp1 <- result2$marginals.fitted[[idx.pre[1]]]
pm.samp2 <- result2$marginals.fitted[[idx.pre[2]]]
pm.samp3 <- result2$marginals.fitted[[idx.pre[3]]]
inla.hpdmarginal(0.95, pm.samp1)
inla.hpdmarginal(0.95, pm.samp2)
inla.hpdmarginal(0.95, pm.samp3)

# Besag Spatial Model

data(Munich, package = "brinla")
g <- system.file("demodata/munich.graph", package = "INLA")
g.file <- inla.read.graph(g)
str(g.file)
formula <- rent ~ 1 + f(location, model = "besag", graph = g)
result <- inla(formula, data = Munich, control.predictor = list(compute = TRUE))
round(result$summary.fixed, 4)
fhat <- result$summary.random$location$mean
map.munich(fhat)
fhat.sd <- result$summary.random$location$sd
map.munich(fhat.sd)

# Penalized Regression Splines (P-splines)

set.seed(1)
n <- 100
x <- seq(0, 1,, n)
f.true <- (sin(2*pi*x^3))^3
y <- f.true + rnorm(n, sd = 0.2)
library(splines)
p <- 25
B.tmp <- bs(x, df = p, intercept = TRUE)
attributes(B.tmp) <- NULL  ## Remove attributes
Bmat <- as(matrix(B.tmp, n, p), 'sparseMatrix')
data.inla <- list(y = y, x = 1:p)
formula <- y ~ -1 + f(x, model = 'rw1', constr = FALSE)
result <- inla(formula, data=data.inla, control.predictor = list(A = Bmat, compute = TRUE))
round(head(result$summary.linear.predictor), 3)
p <- bri.band.ggplot(result, ind = 1:n, type = 'linear')
p + geom_point(aes(y = y, x = 1:n)) + geom_line(aes(y = f.true, x = 1:n), linetype = 2)
formula <- y ~ -1 + f(x, model = 'rw2', constr = FALSE)

# Adaptive Spline Smoothing\index{adaptive smoothing}

set.seed(1)
n <- 100
x <- seq(0, 1,, n)
f.true <- (sin(2*pi*x^3))^3
y <- f.true + rnorm(n, sd = 0.2)
adapt <- bri.adapt.prior(x, nknot = 5, type = 'spde')
data.inla <- list(y = y, x = adapt$x.ind)
formula <- y ~ -1 + f(x, model = adapt)
result <- inla(formula, data = data.inla)

# Generalized Nonparametric Regression Models

set.seed(2)
n <- 200  #sample size
x <- seq(0, 6,, n)
eta <- sin(x)
Ntrials <- sample(c(1, 5, 10, 15), size = n, replace = TRUE)
prob1 <- exp(eta)/(1 + exp(eta))  ## logit link
prob2 <- pnorm(eta)  ## probit link
prob3 <- 1 - exp(-exp(eta))  ## complementary log-log link
y1 <- rbinom(n, size = Ntrials, prob = prob1)
y2 <- rbinom(n, size = Ntrials, prob = prob2)
y3 <- rbinom(n, size = Ntrials, prob = prob3)
data1 <- list(y = y1, x = x)
data2 <- list(y = y2, x = x)
data3 <- list(y = y3, x = x)
formula <- y ~ -1 + f(x, model = "rw2", constr = FALSE)
result1 <- inla(formula, family = "binomial", data = data1, Ntrials = Ntrials, control.predictor = list(compute = TRUE))
result2 <- inla(formula, family = "binomial", data = data2, Ntrials = Ntrials, control.predictor = list(compute = TRUE), control.family = list(link = 'probit'))
result3 <- inla(formula, family = "binomial", data = data3, Ntrials = Ntrials,  control.predictor = list(compute = TRUE), control.family = list(link = 'cloglog'))
p1 <- bri.band.ggplot(result1, name = 'x', alpha = 0.05, type = 'random')
p1 + geom_line(aes(y = eta), linetype = 2)
set.seed(2)
n <- 200  #sample size
x <- seq(0, 6,, n)
E <- sample(1:10, n, replace = TRUE)
lambda <- E*exp(sin(x))
y4 <- rpois(n, lambda = lambda)
data4 <- list(y = y4, x = x)
formula <- y ~ -1 + f(x, model = "rw2", constr = FALSE)
result4 <- inla(formula, family = "poisson", data = data4, E = E, control.predictor = list(compute = TRUE))
lamb.hat <- E*result4$summary.fitted$mean
yhat.sd <- E*result4$summary.fitted$sd
lamb.lb <- E*result4$summary.fitted$'0.025quant'
lamb.ub <- E*result4$summary.fitted$'0.975quant'
data(Tokyo, package = 'INLA')
str(Tokyo)
formula <- y ~ -1 + f(time, model = "rw2", cyclic = TRUE)
result <- inla(formula, family = "binomial", Ntrials = n, data = Tokyo, control.predictor = list(compute = TRUE))
bri.band.plot(result, name = 'time', alpha = 0.05, type = 'random', xlab = 'Day', ylab = '')
bri.band.plot(result, alpha = 0.05, type = 'fitted', ylim = c(0, 1), xlab = 'Day', ylab = 'Probability')
points(Tokyo$time, Tokyo$y/2, cex = 0.5)
time.new <- seq(length(Tokyo$time) + 1, length.out = 3)
time.pred <- c(Tokyo$time, time.new)
y.pred <- c(Tokyo$y, rep(NA, length(time.new)))
n.pred <- c(Tokyo$n, rep(1, length(time.new)))
Tokyo.pred <- list(y = y.pred, time = time.pred, n = n.pred)
result <- inla(formula, family = "binomial", Ntrials = n, data = Tokyo.pred, control.predictor = list(compute = TRUE, link = 1))
link <- rep(NA, length(y.pred))
link[which(is.na(y.pred))] <- 1
ID <- result$summary.random$time$ID
idx.pred <- sapply(time.new, function(x) which(ID==x))
round(result$summary.random$time[idx.pred,], 4)
round(result$summary.fitted.values[which(is.na(y.pred)),], 4)

# Excursion Set with Uncertainty

data(fossil, package = 'brinla')
str(fossil)
min(diff(sort(fossil$age)))/diff(range(fossil$age))
age.new <- inla.group(fossil$age, n = 100)
(length(unique(age.new)))
inla.data <- list(y = fossil$sr, x = age.new)
formula <- y ~ -1 + f(x, model = 'rw2', constr = FALSE, scale.model = TRUE)
result <- inla(formula, data = inla.data, control.compute = list(config = TRUE))
mar.x <- result$marginals.random$x  #marginal posterior
mar.prob <- 1 - sapply(mar.x, function(x) inla.pmarginal(0.74, x))
result$summary.random$x$ID[mar.prob > 0.95]
res.exc <- excursions.brinla(result, name = 'x', u = 0.74, alpha = 0.05, type = '>', method = 'NI')
res.exc$E
round(head(res.exc$F), 4)
bri.excursions.ggplot(res.exc)
data(Tokyo, package = 'INLA')
formula <- y ~ -1 + f(time, model = "rw2", cyclic = TRUE)
result <- inla(formula, family = "binomial", Ntrials = n, data = Tokyo, control.predictor = list(compute=TRUE), control.compute = list(config = TRUE))
u.fitted <- 0.3  ## threshold
mar.fitted <- result$marginals.fitted.values
mar.prob<- 1-sapply(mar.fitted,function(x) inla.pmarginal(u.fitted,x))
Tokyo$time[mar.prob >= 0.95]
u.pred <- log(u.fitted/(1 - u.fitted))
res.exc <- excursions.brinla(result, name = 'time', u = u.pred, type = '>', method = 'NIQC', alpha = 0.05)
res.exc$E
sessionInfo()
