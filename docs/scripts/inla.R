
# Latent Gaussian Models (LGMs)

library(INLA)
names(inla.models()$likelihood)

# Gaussian Markov random fields (GMRFs)


# Laplace Approximation and INLA

set.seed(1)
n <- 50
rho <- 0.8 
prec <- 10 
E <- sample(c(5, 4, 10, 12), size = n, replace = TRUE) 
eta <- arima.sim(list(order=c(1,0,0), ar=rho), n=n, sd=sqrt(1/prec))
y <- rpois(n, E*exp(eta)) 
data <- list(y = y, x = 1:n, E = E) 
formula <- y ~ f(x, model = "ar1")
result.sla <- inla(formula, family = "poisson", data = data, E = E)
result.gau <-  inla(formula, family = "poisson", data = data, E = E, control.inla = list(strategy='gaussian', int.strategy="eb"))
result.la <-  inla(formula, family = "poisson", data = data, E = E, control.inla = list(strategy = 'laplace', int.strategy = "grid", dz = 0.1, diff.logdens = 20))
i <- 50
marg.lst <- list(result.sla$marginals.random$x[[i]], result.gau$marginals.random$x[[i]], result.la$marginals.random$x[[i]])
names(marg.lst) <- as.factor(c(1,2,3))
pp <- data.frame(do.call(rbind, marg.lst))
pp$method <- rep(names(marg.lst), times = sapply(marg.lst, nrow))
library(ggplot2)
ggplot(pp, aes(x = x, y = y, linetype = method)) + geom_line(show.legend = FALSE) + ylab("density") + xlab("") + theme_bw(base_size = 20) 
marg.prec <- prec*(1 - rho^2) 
(marg.sd <- 1/sqrt(marg.prec))
res.hyper <- inla.hyperpar(result.sla)
library(brinla)
p1 <- bri.hyperpar.plot(result.la)
p2 <- bri.hyperpar.plot(res.hyper)
pp <- rbind(p1, p2)
pp$method <- as.factor(c(rep(1, dim(p1)[1]), rep(2, dim(p2)[1])))
ggplot(pp, aes(x = x, y = y, linetype = method)) + geom_line(show.legend = FALSE) + facet_wrap(~parameter, scales = "free") + ylab("density") + xlab("") + theme_bw(base_size = 20) 

# INLA Failure

set.seed(1)
n <- 100
u <- rnorm(n)
eta <- 1 + u
p <- exp(eta)/(1+exp(eta))
y <- rbinom(n, size = 1,  prob = p)
data <- data.frame(y = y, x = 1:n)
result <- inla(y ~ 1 + f(x, model = "iid"), data = data, family = "binomial", Ntrials = 1)
round(bri.hyperpar.summary(result), 4)
round(result$summary.fixed, 4)

# Extensions

result <- inla(..., control.compute = list(dic = TRUE, waic = TRUE))
result <- inla(..., control.compute = list(cpo = TRUE), ...)
inla(..., control.inla = list(int.strategy="grid", diff.logdens=4))
inla(..., control.inla = list(strategy = "laplace", npoints = 21))
inla(..., control.inla = list(lincomb.derived.only = FALSE), ...) 
inla(..., control.inla = list(lincomb.derived.correlation.matrix = TRUE), ...)
sessionInfo()
