# Bayesian Regression with INLA

A book by [Xiaofeng Wang](https://filer.case.edu/xxw17/), [Ryan Yue](http://zicklin.baruch.cuny.edu/faculty/profiles/yu-ryan-yue) and [Julian Faraway](http://people.bath.ac.uk/jjf23/)

INLA stands for Integrated Nested Laplace Approximations. It is used for fitting 
[Latent Gaussian models](https://tgmstat.wordpress.com/2013/10/16/latent-gaussian-models-inla/) (LGM). 
LGMs include a wide range of commonly
used regression models. Unlike MCMC which uses simulation methods, INLA
uses approximation methods for Bayesian model fitting. Within the class
of LGMs, INLA can fit models much faster than MCMC-based methods.

## Chapters

1. Introduction: [intro.R](scripts/intro.R) with [output](scripts/intro.md) and [text online](https://julianfaraway.github.io/brinlabook/)
2. Theory of INLA: [inla.R](scripts/inla.R)
3. Linear Regression: [blr.R](scripts/blr.R)
4. Generalized Linear Models: [glm.R](scripts/glm.R)
5. Generalized Linear Mixed Models [glmm.R](scripts/glmm.R) with [output](scripts/glmm.md) and [text online](https://julianfaraway.github.io/brinlabook/chaglmm.html)
6. Survival Analysis: [surv.R](scripts/surv.R)
7. Random Walk Models for Smoothing: [npr.R](scripts/npr.R)
8. Gaussian Process Regression [gpr.R](scripts/gpr.R) with [output](scripts/gpr.md)
9. Generalized Additive Models: [gam.R](scripts/gam.R) and [text online](https://julianfaraway.github.io/brinlabook/chagpr.html)
10. Errors-in-Variables Regression: [eiv.R](scripts/eiv.R)
11. Miscellaneous Topics: [misc.R](scripts/misc.R) and [extreme values text online](https://julianfaraway.github.io/brinlabook/chamisc.html#sec:extreval)


## R package

The [brinla](https://github.com/julianfaraway/brinla) R package contains data and functions
to support the book.

Visit the  [INLA](http://www.r-inla.org) website to learn much more including
[how to install the INLA R package](https://www.r-inla.org/download-install).

Install the `brinla` package with:

```
devtools::install_github("julianfaraway/brinla")
```

## Errata

Here are the [errata](errata.md). If you find any other errata, please let us know according to the chapter: Ch3, 4, 6, 10 (Xiaofeng Wang),
Ch2, 7 or 9 (Ryan Yue) or Ch1, 5 or 8 (Julian Faraway).

## Examples

Here are some [examples](examples/index.md)

## Purchase

The book can be purchased at the usual online outlets or the [publishers](https://www.crcpress.com/Bayesian-Regression-Modeling-with-INLA/Wang-Ryan-Faraway/p/book/9781498727259)