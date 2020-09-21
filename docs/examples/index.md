---
title: INLA for Bayesian Regression Models
author: "[Julian Faraway](http://people.bath.ac.uk/jjf23/)"
output:
 html_document
---

INLA stands for Integrated Nested Laplace Approximations. It is used for fitting 
[Latent Gaussian models](https://tgmstat.wordpress.com/2013/10/16/latent-gaussian-models-inla/) (LGM). 
LGMs include a wide range of commonly
used regression models. Unlike MCMC which uses simulation methods, INLA
uses approximation methods for Bayesian model fitting. Within the class
of LGMs, INLA can fit models much faster than MCMC-based methods.

Visit the  [INLA](http://www.r-inla.org) website to learn much more including
[how to install the INLA R package](http://www.r-inla.org/download).

I am the author of the book entitled [Bayesian Regression Models](http://julianfaraway.github.io/brinla/) with
[Xiaofeng Wang](https://filer.case.edu/xxw17/) and [Ryan Yue](http://zicklin.baruch.cuny.edu/faculty/profiles/yu-ryan-yue)

You will need to [install the brinla R package](https://github.com/julianfaraway/brinla) to run
many of the examples described here.

# Worked Examples

- [Introduction](intro.md) - a simple example concerning Hubble's law.
- [Linear regression](chicago.html) - Chicago Insurance data
- [One-way ANOVA](reeds.html) - just one random effect
- [Ridge regression](ridge.html) - Ridge regression with meat spectroscopy data
- [Gaussian Process regression](gpreg.html) - GP regression on fossil data
- [Confidence bands for smoothness](smoothband.html) - GP regression to determine uncertainty regarding smoothness
- [Non-stationary smoothing](nonstat.html) - GP regression with variable smoothing
- [Generalized Extreme Values](gev.html) - fitting maximum annual river flows
- [Define your own prior](prior.html) - using a half Cauchy prior for the SD of a random effect.
- [Linear regression with bounded parameters](frenchpres.html) - French presidential election example with slope parameters bounded in [0,1]

See also [linear mixed models examples](../inla/index.html) from my
[Extending Linear Models with R](http://people.bath.ac.uk/jjf23/ELM/index.html) book.

