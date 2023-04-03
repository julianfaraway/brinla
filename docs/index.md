# Bayesian Regression with INLA

A book by [Xiaofeng Wang](https://xfw.github.io/), [Ryan Yue](https://zicklin.baruch.cuny.edu/faculty-profile/yu-yue/) and [Julian Faraway](https://julianfaraway.github.io/)

INLA stands for Integrated Nested Laplace Approximations. It is used for fitting 
[Latent Gaussian models](https://tgmstat.wordpress.com/2013/10/16/latent-gaussian-models-inla/) (LGM). 
LGMs include a wide range of commonly
used regression models. Unlike MCMC which uses simulation methods, INLA
uses approximation methods for Bayesian model fitting. Within the class
of LGMs, INLA can fit models much faster than MCMC-based methods.

## Chapters

1. Introduction: [intro.R](scripts/intro.R) and [text online](https://julianfaraway.github.io/brinlabook/)
2. Theory of INLA: [inla.R](scripts/inla.R) and 
[text online](https://julianfaraway.github.io/brinlabook/theory-of-inla.html)
3. Linear Regression: [blr.R](scripts/blr.R)
4. Generalized Linear Models: [glm.R](scripts/glm.R)
5. Generalized Linear Mixed Models [glmm.R](scripts/glmm.R) and [text online](https://julianfaraway.github.io/brinlabook/chaglmm.html)
6. Survival Analysis: [surv.R](scripts/surv.R)
7. Random Walk Models for Smoothing: [npr.R](scripts/npr.R) and
[text online](https://julianfaraway.github.io/brinlabook/chanpr.html)
8. Gaussian Process Regression [gpr.R](scripts/gpr.R) with [text online](https://julianfaraway.github.io/brinlabook/chagpr.html)
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

The book can be purchased at the usual online outlets or from the publishers [Routledge](https://eur01.safelinks.protection.outlook.com/?url=https%3A%2F%2Fwww.routledge.com%2FBayesian-Regression-Modeling-with-INLA%2FWang-Ryan-Faraway%2Fp%2Fbook%2F9780367572266%3Futm_source%3Dauthor%26utm_medium%3Dshared_link%26utm_campaign%3DB043141_jm1_5ll_6rm_t081_1al_julianfarawayauthorshare&data=05%7C01%7Cjjf23%40bath.ac.uk%7C5229b0dc8d564222c10108db31fa6e3f%7C377e3d224ea1422db0ad8fcc89406b9e%7C0%7C0%7C638158723883874621%7CUnknown%7CTWFpbGZsb3d8eyJWIjoiMC4wLjAwMDAiLCJQIjoiV2luMzIiLCJBTiI6Ik1haWwiLCJXVCI6Mn0%3D%7C3000%7C%7C%7C&sdata=vhy4wINWPGtr6s8DMhrgD8ybz8usUz4B7Fuy9g0eG28%3D&reserved=0)