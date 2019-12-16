# Errata in Bayesian Regression with INLA

- p7 The code for the fit based on the Ussher prior does not match the text. The code should read:
```
imod <- inla(y ~ x -1, family="gaussian", control.fixed=list(mean=uhub, prec=1/((0.05*uhub)^2)), data=hubble)
```
As it happens, this makes no important difference to the numerical outcome.

- p104: the `inla.contrib.sd` function is no longer available in the INLA package. You can compute this information as `bri.hyperpar.summary` function as demonstrated on the following page.

- p111: 3rd line of R code should now read `mean(lvsamp[,'site:3'] > lvsamp[,'site:1'])`. At some point, the parameter labeling in INLA changed which broke the previous version of this line. 

- p273: should be 340 cases, not visits.