# Errata in Bayesian Regression with INLA

- p104: the `inla.contrib.sd` function is no longer available in the INLA package. You can compute this information as `brinla.hyperpar.summary` function as demonstrated on the following page.

- p111: 3rd line of R code should now read `mean(lvsamp[,'site:3'] > lvsamp[,'site:1'])`. At some point, the parameter labeling in INLA changed which broke the previous version of this line. 

- p273: should be 340 cases, not visits.