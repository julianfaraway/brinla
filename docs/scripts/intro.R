
# Quick Start

data(hubble, package = "brinla")
plot(y ~ x, xlab = "Distance(Mpc)", ylab = "Velocity(km/s)",
    data = hubble)
lmod <- lm(y ~ x - 1, data = hubble)
coef(lmod)
hubtoage <- function(x) 3.09e+19/(x * 60^2 * 24 * 365.25 * 1e+09)
hubtoage(coef(lmod))
(bci <- confint(lmod))
hubtoage(bci)
library(INLA)
imod <- inla(y ~ x - 1, family = "gaussian",
    control.fixed = list(prec = 1e-09), data = hubble)
(ibci <- imod$summary.fixed)
plot(imod$marginals.fixed$x, type = "l", xlab = "beta",
    ylab = "density", xlim = c(60, 100))
abline(v = ibci[c(3, 5)], lty = 2)
hubtoage(ibci[c(1,3,4,5,6)])
ageden <- inla.tmarginal(hubtoage, imod$marginals.fixed$x)
plot(ageden, type = "l", xlab = "Age in billions of years",
    ylab = "density")
abline(v = hubtoage(ibci[c(3, 5)]), lty = 2)
hubtoage(c(10, 15, 20))
imod <- inla(y ~ x - 1, family = "gaussian",
    control.fixed = list(mean = 65, prec = 1/(12^2)), data = hubble)
(ibci <- imod$summary.fixed)
hubtoage(ibci[c(1, 3, 4, 5, 6)])
(uhub <- hubtoage((2016 + 4004 - 1)/1e+09))
imod <- inla(y ~ x - 1, family = "gaussian",
    control.fixed = list(mean = 65, prec = 1/(12^2)), data = hubble)
(ibci <- imod$summary.fixed)
hubtoage(ibci[c(1, 3, 4, 5, 6)])

# Bayes Theory


# Prior and Posterior Distributions


# Model Checking


# Model Selection


# Hypothesis testing


# Bayesian Computation

sessionInfo()
