# Errata in Bayesian Regression with INLA

INLA is not strongly dependent on random number generation (unlike MCMC!) but it does sometimes
use simulation as part of the calculation. You will find that if you repeat the same command,
you may get slightly different output. Hence, you may find that if you repeat the same commands
as seen in the book, you might not get exactly the same answer. This is usually not a serious concern bearing in mind that the last letter of INLA stands for approximation and exact answers
are not to be expected.

- p7 The code for the fit based on the Ussher prior does not match the text. The code should read:
```
imod <- inla(y ~ x -1, family="gaussian", control.fixed=list(mean=uhub, prec=1/((0.05*uhub)^2)), data=hubble)
```
As it happens, this makes no important difference to the numerical outcome.

- p104: the `inla.contrib.sd` function is no longer available in the INLA package. You can compute this information as `bri.hyperpar.summary` function as demonstrated on the following page.

- p111: 3rd line of R code should now read `mean(lvsamp[,'site:3'] > lvsamp[,'site:1'])`. At some point, the parameter labeling in INLA changed which broke the previous version of this line. 

- p118: Both `inla()` calls on this page now also require the option `control.compute=list(return.marginals.predictor=TRUE)` because these quantities are expensive
to compute and are not always needed.

- p139: Section 5.7.1 should be ignored. This correction is not available in the current version of 
INLA because it is not regarded as sufficiently stable to be relied upon.

- p218: The specification for `inla.spde2.matern` has changed arguments requiring

```
spde <- inla.spde2.matern(mesh, alpha=alpha, constr = FALSE,
  prior.tau = tau0,
  prior.kappa = kappa0,
  theta.prior.prec = 1e5)
```

The value of `theta.prior.prec` is intentionally large to ensure prior is respected.

- p273: should be 340 cases, not visits.