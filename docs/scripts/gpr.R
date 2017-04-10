
# Introduction

data(fossil, package = "brinla")
fossil$sr <- (fossil$sr - 0.7) * 100
pf <- ggplot(fossil, aes(age, sr)) + geom_point() + xlab("Age") + ylab("Strontium Ratio")
pf + geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"))
nbasis <- 25
library(INLA)
mesh <- inla.mesh.1d(seq(min(fossil$age), max(fossil$age), length.out = nbasis), 
    degree = 2)
alpha <- 2
nu <- alpha - 1/2
sigma0 <- sd(fossil$sr)
rho0 <- 0.25 * (max(fossil$age) - min(fossil$age))
kappa0 <- sqrt(8 * nu)/rho0
tau0 <- 1/(4 * kappa0^3 * sigma0^2)^0.5
spde <- inla.spde2.matern(mesh, alpha = alpha, B.tau = cbind(log(tau0), 
    1, 0), B.kappa = cbind(log(kappa0), 0, 1), theta.prior.prec = 1e-04)
A <- inla.spde.make.A(mesh, loc = fossil$age)
index <- inla.spde.make.index("sinc", n.spde = spde$n.spde)
st.est <- inla.stack(data = list(y = fossil$sr), A = list(A), effects = list(index), 
    tag = "est")
formula <- y ~ -1 + f(sinc, model = spde)
data <- inla.stack.data(st.est)
result <- inla(formula, data = data, family = "gaussian", control.predictor = list(A = inla.stack.A(st.est), 
    compute = TRUE))
ii <- inla.stack.index(st.est, "est")$data
plot(sr ~ age, fossil)
tdx <- fossil$age
lines(tdx, result$summary.fitted.values$mean[ii])
lines(tdx, result$summary.fitted.values$"0.025quant"[ii], lty = 2)
lines(tdx, result$summary.fitted.values$"0.975quant"[ii], lty = 2)
library(brinla)
fg <- bri.gpr(fossil$age, fossil$sr)
plot(sr ~ age, fossil, pch = 20)
lines(fg$xout, fg$mean)
lines(fg$xout, fg$ucb, lty = 2)
lines(fg$xout, fg$lcb, lty = 2)

# Penalized Complexity Priors

spde <- inla.spde2.pcmatern(mesh, alpha = alpha, prior.range = c(5, 0.05), 
    prior.sigma = c(2, 0.05))
formula <- y ~ -1 + f(sinc, model = spde)
resultpc <- inla(formula, data = data, family = "gaussian", control.predictor = list(A = inla.stack.A(st.est), 
    compute = TRUE))
pcmod <- bri.gpr(fossil$age, fossil$sr, pcprior = c(5, 2))

# Credible bands for smoothness

library(brinla)
errorsd <- bri.hyper.sd(result$marginals.hyperpar[[1]])
mres <- inla.spde.result(result, "sinc", spde)
mv <- mres$marginals.variance.nominal[[1]]
sigmad <- as.data.frame(inla.tmarginal(function(x) sqrt(x), mv))
rhod <- mres$marginals.range.nominal[[1]]
plot(y ~ x, errorsd, type = "l", xlab = "sr", ylab = "density")
plot(y ~ x, sigmad, type = "l", xlab = "sr", ylab = "density")
plot(rhod, type = "l", xlab = "age", ylab = "density")
exp(mres$summary.log.kappa[c(4, 6)])
kappa0 <- exp(mres$summary.log.kappa["0.025quant"])[, ]
sigma02 <- exp(mres$summary.log.variance.nominal["0.5quant"])[, ]
tau0 <- 1/(4 * kappa0^3 * sigma02)^0.5
spde <- inla.spde2.matern(mesh, alpha = alpha, constr = FALSE, B.tau = cbind(log(tau0)), 
    B.kappa = cbind(log(kappa0)))
formula <- y ~ -1 + f(sinc, model = spde)
resulta <- inla(formula, data = data, family = "gaussian", control.predictor = list(A = inla.stack.A(st.est), 
    compute = TRUE))
kappa0 <- exp(mres$summary.log.kappa["0.975quant"])[, ]
sigma02 <- exp(mres$summary.log.variance.nominal["0.5quant"])[, ]
tau0 <- 1/(4 * kappa0^3 * sigma02)^0.5
spde <- inla.spde2.matern(mesh, alpha = alpha, constr = FALSE, B.tau = cbind(log(tau0)), 
    B.kappa = cbind(log(kappa0)))
formula <- y ~ -1 + f(sinc, model = spde)
resultb <- inla(formula, data = data, family = "gaussian", control.predictor = list(A = inla.stack.A(st.est), 
    compute = TRUE))
ii <- inla.stack.index(st.est, "est")$data
plot(sr ~ age, fossil, pch = 20)
tdx <- fossil$age
lines(tdx, resulta$summary.fitted.values$mean[ii], lty = 2)
lines(tdx, resultb$summary.fitted.values$mean[ii], lty = 1)
fg <- bri.smoothband(fossil$age, fossil$sr)
plot(sr ~ age, fossil, pch = 20)
lines(fg$xout, fg$rcb, lty = 1)
lines(fg$xout, fg$scb, lty = 2)

# Non-stationary Fields

set.seed(1)
n <- 100
x <- seq(0, 1, length = n)
f.true <- (sin(2 * pi * (x)^3))^3
y <- f.true + rnorm(n, sd = 0.2)
td <- data.frame(y = y, x = x, f.true)
nbasis <- 25
mesh <- inla.mesh.1d(seq(0, 1, length.out = nbasis), degree = 2)
alpha <- 2
nu <- alpha - 1/2
sigma0 <- sd(y)
rho0 <- 0.1
kappa0 <- sqrt(8 * nu)/rho0
tau0 <- 1/(4 * kappa0^3 * sigma0^2)^0.5
spde <- inla.spde2.matern(mesh, alpha = alpha, B.tau = cbind(log(tau0), 
    1, 0), B.kappa = cbind(log(kappa0), 0, 1), theta.prior.prec = 1e-04)
A <- inla.spde.make.A(mesh, loc = td$x)
index <- inla.spde.make.index("sinc", n.spde = spde$n.spde)
st.est <- inla.stack(data = list(y = td$y), A = list(A), effects = list(index), 
    tag = "est")
formula <- y ~ -1 + f(sinc, model = spde)
data <- inla.stack.data(st.est)
result <- inla(formula, data = data, family = "gaussian", control.predictor = list(A = inla.stack.A(st.est), 
    compute = TRUE))
ii <- inla.stack.index(st.est, "est")$data
plot(y ~ x, td, col = gray(0.75))
tdx <- td$x
lines(tdx, result$summary.fitted.values$mean[ii])
lines(tdx, f.true, lty = 2)
basis.T <- as.matrix(inla.mesh.basis(mesh, type = "b.spline", n = 5, degree = 2))
basis.K <- as.matrix(inla.mesh.basis(mesh, type = "b.spline", n = 5, degree = 2))
spde <- inla.spde2.matern(mesh, alpha = alpha, B.tau = cbind(basis.T[-1, 
    ], 0), B.kappa = cbind(0, basis.K[-1, ]/2), theta.prior.prec = 1e-04)
formula <- y ~ -1 + f(sinc, model = spde)
result <- inla(formula, data = data, family = "gaussian", control.predictor = list(A = inla.stack.A(st.est), 
    compute = TRUE))
plot(y ~ x, td, col = gray(0.75))
lines(tdx, result$summary.fitted.values$mean[ii])
lines(tdx, f.true, lty = 2)
fg <- bri.nonstat(td$x, td$y)
plot(y ~ x, td, col = gray(0.75))
lines(f.true ~ x, td, lty = 2)
lines(fg$xout, fg$mean)

# Interpolation with Uncertainty

x <- c(0, 0.1, 0.2, 0.3, 0.7, 1)
y <- c(0, 0.5, 0, -0.5, 1, 0)
td <- data.frame(x, y)
nbasis <- 100
alpha <- 2
mesh <- inla.mesh.1d(seq(0, 1, length.out = nbasis), degree = 2)
spde <- inla.spde2.pcmatern(mesh, alpha = alpha, prior.range = c(0.05, 
    0.1), prior.sigma = c(5, 0.05))
A <- inla.spde.make.A(mesh, loc = td$x)
ngrid <- 101
Ap <- inla.spde.make.A(mesh, loc = seq(0, 1, length.out = ngrid))
index <- inla.spde.make.index("sinc", n.spde = spde$n.spde)
st.est <- inla.stack(data = list(y = td$y), A = list(A), effects = list(index), 
    tag = "est")
st.pred <- inla.stack(data = list(y = NA), A = list(Ap), effects = list(index), 
    tag = "pred")
formula <- y ~ -1 + f(sinc, model = spde)
sestpred <- inla.stack(st.est, st.pred)
result <- inla(formula, data = inla.stack.data(sestpred), family = "gaussian", 
    control.predictor = list(A = inla.stack.A(sestpred), compute = TRUE), 
    control.family(hyper = list(prec = list(fixed = TRUE, initial = 1e+08))))
ii <- inla.stack.index(sestpred, tag = "pred")$data
plot(y ~ x, td, pch = 20, ylim = c(-2, 2))
tdx <- seq(0, 1, length.out = ngrid)
lines(tdx, result$summary.linear.pred$mean[ii])
lines(tdx, result$summary.linear.pred$"0.025quant"[ii], lty = 2)
lines(tdx, result$summary.linear.pred$"0.975quant"[ii], lty = 2)
jj <- which(tdx == 0.5)
margpred <- result$marginals.linear[ii]
plot(margpred[[jj]], type = "l", ylab = "Density")
mres <- inla.spde.result(result, "sinc", spde)
exp(mres$summary.log.range.nominal[c(2, 4, 5, 6, 7)])
sqrt(exp(mres$summary.log.variance.nominal[c(2, 4, 5, 6, 7)]))
spde <- inla.spde2.pcmatern(mesh, alpha = alpha, prior.range = c(0.033269, 
    NA), prior.sigma = c(0.53804, NA))
resultl <- inla(formula, data = inla.stack.data(sestpred), family = "gaussian", 
    control.predictor = list(A = inla.stack.A(sestpred), compute = TRUE), 
    control.family(hyper = list(prec = list(fixed = TRUE, initial = 1e+08))))
spde <- inla.spde2.pcmatern(mesh, alpha = alpha, prior.range = c(0.3462, 
    NA), prior.sigma = c(0.53804, NA))
resulth <- inla(formula, data = inla.stack.data(sestpred), family = "gaussian", 
    control.predictor = list(A = inla.stack.A(sestpred), compute = TRUE), 
    control.family(hyper = list(prec = list(fixed = TRUE, initial = 1e+08))))
plot(y ~ x, td, pch = 20)
tdx <- seq(0, 1, length.out = ngrid)
lines(tdx, result$summary.linear.pred$mean[ii])
lines(tdx, resultl$summary.linear.pred$mean[ii], lty = 2)
lines(tdx, resulth$summary.linear.pred$mean[ii], lty = 3)

# Survival response

data(larynx, package = "brinla")
nbasis <- 25
alpha <- 2
xspat <- larynx$age
mesh <- inla.mesh.1d(seq(min(xspat), max(xspat), length.out = nbasis), 
    degree = 2)
spde <- inla.spde2.pcmatern(mesh, alpha = alpha, prior.range = c(20, 0.1), 
    prior.sigma = c(10, 0.05))
A <- inla.spde.make.A(mesh, loc = xspat)
index <- inla.spde.make.index("sinc", n.spde = spde$n.spde)
st.est <- inla.stack(data = list(time = larynx$time, censor = larynx$delta), 
    A = list(A), effects = list(index), tag = "est")
formula <- inla.surv(time, censor) ~ 0 + f(sinc, model = spde)
data <- inla.stack.data(st.est)
result <- inla(formula, data = data, family = "weibull.surv", control.predictor = list(A = inla.stack.A(st.est), 
    compute = TRUE))
ii <- inla.stack.index(st.est, "est")$data
lcdf <- data.frame(result$summary.linear.predictor[ii, ], larynx)
alpha <- result$summary.hyperpar[1, 1]
lambda <- exp(lcdf$mean)
lcdf$exptime <- lambda^(-1/alpha) * gamma(1/alpha + 1)
lambda <- exp(lcdf$X0.025quant)
lcdf$lcb <- lambda^(-1/alpha) * gamma(1/alpha + 1)
lambda <- exp(lcdf$X0.975quant)
lcdf$ucb <- lambda^(-1/alpha) * gamma(1/alpha + 1)
p <- ggplot(data = lcdf, aes(x = age, y = time)) + geom_point()
p + geom_line(aes(x = age, y = exptime)) + geom_line(aes(x = age, y = ucb), 
    linetype = 2) + geom_line(aes(x = age, y = lcb), linetype = 2)
lambda <- exp(lcdf$mean)
lcdf$hazard <- alpha * lcdf$age^(alpha - 1) * lambda
lambda <- exp(lcdf$X0.025quant)
lcdf$hazlo <- alpha * lcdf$age^(alpha - 1) * lambda
lambda <- exp(lcdf$X0.975quant)
lcdf$hazhi <- alpha * lcdf$age^(alpha - 1) * lambda
ggplot(data = lcdf, aes(x = age, y = hazard)) + geom_line() + geom_line(aes(x = age, 
    y = hazlo), lty = 2) + geom_line(aes(x = age, y = hazhi), lty = 2)
sessionInfo()
