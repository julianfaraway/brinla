#' ---
#' title: "Bayesian Regression Modeling with INLA"
#' author: "Chapter One: Introduction"
#' output:
#'   github_document
#' ---
#+ echo=FALSE
knitr::opts_chunk$set(comment=NA, 
                      echo = TRUE,
                      fig.path="figs/", 
                      dev = 'svglite',  
                      fig.ext = ".svg",
                      warning=FALSE, 
                      message=FALSE)
library(svglite)
ggplot2::theme_set(ggplot2::theme_bw())
par(mgp=c(1.5,0.5,0), mar=c(3.1,3.1,3.1,0), pch=20)
#' Code from [Bayesian Regression Modeling with INLA](http://julianfaraway.github.io/brinla/)
#' 
#' # Quick Start

data(hubble, package = "brinla")
#+hubble
plot(y ~ x, xlab = "Distance(Mpc)", ylab = "Velocity(km/s)", data = hubble)
#' Use standard linear model
lmod <- lm(y ~ x - 1, data = hubble)
coef(lmod)
hubtoage <- function(x) 3.09e+19/(x * 60^2 * 24 * 365.25 * 1e+09)
hubtoage(coef(lmod))
(bci <- confint(lmod))
hubtoage(bci)
library(INLA)
#' Use a weakly informative prior
imod <- inla(y ~ x - 1, family = "gaussian", control.fixed = list(prec = 1e-09), 
    data = hubble)
(ibci <- imod$summary.fixed)
#' The mean and mode are similar to the published text output but the sd
#' and interval have increased
#' 
#+ hubmarginals
plot(imod$marginals.fixed$x, type = "l", xlab = "beta", ylab = "density", 
    xlim = c(60, 100))
abline(v = ibci[c(3, 5)], lty = 2)
hubtoage(ibci[c(1, 3, 4, 5, 6)])
ageden <- inla.tmarginal(hubtoage, imod$marginals.fixed$x)
#+ hubage
plot(ageden, type = "l", xlab = "Age in billions of years", ylab = "density")
abline(v = hubtoage(ibci[c(3, 5)]), lty = 2)
#' Use a more informative prior
hubtoage(c(10, 15, 20))
imod <- inla(y ~ x - 1, family = "gaussian", control.fixed = list(mean = 65, 
    prec = 1/(12^2)), data = hubble)
(ibci <- imod$summary.fixed)
#' The mean and mode are similar to the published text output but the sd
#' and interval have increased
#' 
hubtoage(ibci[c(1, 3, 4, 5, 6)])
#' Use the Ussher prior which is very wrong
(uhub <- hubtoage((2016 + 4004 - 1)/1e+09))
imod <- inla(y ~ x - 1, family = "gaussian", control.fixed = list(mean = uhub, 
    prec = 1/((0.05 * uhub)^2)), data = hubble)
(ibci <- imod$summary.fixed)
hubtoage(ibci[c(1, 3, 4, 5, 6)])

#' # Version information
sessionInfo()
