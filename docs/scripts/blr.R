library(brinla)

# Introduction\label{intro}


# Bayesian Inference for Linear Regression

library(inla)
library(ggplot2, GGally)
data(usair, package = "brinla")
pairs.chart <- ggpairs(usair, lower = list(continuous = "cor"), upper = list(continuous = "points", combo = "dot")) + ggplot2::theme(axis.text = element_text(size = 6)) + ggplot2::theme_bw() 
print(pairs.chart)
usair.formula1 <-  SO2 ~ negtemp + manuf + wind + precip + days
usair.lm1 <- lm(usair.formula1, data = usair)
round(coef(summary(usair.lm1)), 4)
round(summary(usair.lm1)$sigma, 4)
usair.inla1 <- inla(usair.formula1, data = usair, control.compute = list(dic = TRUE, cpo = TRUE))
round(usair.inla1$summary.fixed, 4)
round(usair.inla1$summary.hyperpar, 4)
inla.emarginal(fun = function(x) 1/sqrt(x), marg = usair.inla1$marginals.hyperpar$`Precision for the Gaussian observations`)
library(brinla)
round(bri.hyperpar.summary(usair.inla1),4)
bri.hyperpar.plot(usair.inla1)
usair.inla2 <- inla(usair.formula1, data = usair, control.compute = list(dic = TRUE, cpo =TRUE), control.fixed = list(mean.intercept = 100, prec.intercept = 10^(-2), mean = list(negtemp = 2, wind = -3, default =0), prec = 1), control.family = list(hyper = list(prec = list(prior="gaussian", param =c(0,1)))))

# Prediction

new.data <- data.frame(negtemp = c(-50, -60, -40), manuf = c(150, 100, 400), pop = c(200, 100, 300), wind = c(6, 7, 8), precip = c(10, 30, 20), days = c(20, 100, 40))
predict(usair.lm1, new.data, se.fit = TRUE)
usair.combined <- rbind(usair, data.frame(SO2 = c(NA, NA, NA), new.data))
usair.link <- c(rep(NA, nrow(usair)), rep(1, nrow(new.data)))
usair.inla1.pred <- inla(usair.formula1, data = usair.combined, control.predictor = list(link = usair.link))
usair.inla1.pred$summary.fitted.values[(nrow(usair)+1):nrow(usair.combined),]

# Model Selection and Checking

library(MASS)
usair.step <- stepAIC(usair.lm1, trace = FALSE)
usair.step$anova
usair.formula2 <- SO2 ~ negtemp + manuf + wind + precip
usair.lm2 <- lm(usair.formula2, data = usair)
round(coef(summary(usair.lm2)), 4)
usair.inla3 <- inla(usair.formula2, data = usair, control.compute = list(dic = TRUE, cpo = TRUE))
round(usair.inla3$summary.fixed, 4)
c(usair.inla1$dic$dic, usair.inla3$dic$dic)
usair.inla3.pred <- inla(usair.formula2, data = usair, control.predictor = list(link = 1, compute = TRUE))
post.predicted.pval <- vector(mode = "numeric", length = nrow(usair))
for(i in (1:nrow(usair))) {
  post.predicted.pval[i] <- inla.pmarginal(q=usair$SO2[i], marginal = usair.inla3.pred$marginals.fitted.values[[i]])
}
hist(post.predicted.pval, main="", breaks = 10, xlab="Posterior predictive p-value")
sum(usair.inla2$cpo$failure)
hist(usair.inla3$cpo$pit, main="", breaks = 10, xlab = "PIT")
qqplot(qunif(ppoints(length(usair.inla3$cpo$pit))), usair.inla3$cpo$pit, main = "Q-Q plot for Unif(0,1)", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(usair.inla3$cpo$pit, distribution = function(p) qunif(p), prob = c(0.1, 0.9))
LPML1 <- sum(log(usair.inla1$cpo$cpo)); LPML3 <- sum(log(usair.inla3$cpo$cpo))
c(LPML1, LPML3)
plot(usair.inla1$cpo$cpo, usair.inla3$cpo$cpo, xlab="CPO for the full model", ylab="CPO for the reduced model")
abline(0,1)
bri.lmresid.plot(usair.inla2)
bri.lmresid.plot(usair.inla2, usair1$negtemp, xlab = "negtemp", smooth =TRUE)

# Robust Regression

usair.inla4 <- inla(usair.formula2, family = "T", data = usair, control.compute = list(dic = TRUE, cpo = TRUE))
round(usair.inla4$summary.fixed, 4)

# Analysis of Variance

data(painrelief, , package = "brinla")
painrelief$PainLevel <- as.factor(painrelief$PainLevel)
painrelief$Codeine <- as.factor(painrelief$Codeine)
painrelief$Acupuncture <- as.factor(painrelief$Acupuncture)
painrelief.inla <- inla(Relief ~ PainLevel + Codeine*Acupuncture, data = painrelief)
round(painrelief.inla$summary.fixed, 4)

# Ridge Regression for Multicollinearity

data(frencheconomy, package = "brinla")
round(cor(frencheconomy[,-1]),4)
n <- nrow(frencheconomy)
frencheconomy$beta1 <- rep(1,n)
frencheconomy$beta2 <- rep(2,n)
frencheconomy$beta3 <- rep(3,n)
# Set the priors
param.beta = list(prec = list(param = c(1.0e-3, 1.0e-3)))
formula.ridge  = IMPORT ~ f(beta1, DOPROD, model="iid", values = c(1,2,3), hyper = param.beta) + f(beta2, STOCK, copy="beta1", fixed=T) + f(beta3, CONSUM, copy="beta1", fixed=T)
frencheconomy.ridge <- inla(formula.ridge, data=frencheconomy)
formula <- IMPORT ~  DOPROD + STOCK + CONSUM
frencheconomy.inla <- inla(formula, data = frencheconomy, control.fixed = list(prec = 1.0e-3), control.family = list(hyper = param.beta))
# compare results
naive.est = frencheconomy.inla$summary.fixed
ridge.est = rbind(frencheconomy.ridge$summary.fixed, frencheconomy.ridge$summary.random$beta1[,-1])
round(naive.est,4)
round(ridge.est,4)

# Regression with Autoregressive Errors

data(nzunemploy, package = "brinla")
nzunemploy$time <- 1:nrow(nzunemploy)
qplot(time, value, data = gather(nzunemploy[,c(2,3,5)], variable, value, -time), geom = "line") + geom_vline(xintercept = 90) + facet_grid(variable ~ ., scale = "free") + ylab("Unemployment rate") + theme_bw()
nzunemploy$centeredadult = with(nzunemploy, adult - mean(adult))
nzunemploy$centeredadult = with(nzunemploy, adult - mean(adult))
formula1 <- youth ~ centeredadult*policy
nzunemploy.inla1 <- inla(formula1, data= nzunemploy)
round(nzunemploy.inla1$summary.fixed, 4)
nzunemploy.res1 <- bri.lmresid.plot(nzunemploy.inla1, type="o")
acf(nzunemploy.res1$resid, main = "")
acf(nzunemploy.res1$resid, type = "partial", main="")
formula2 <- youth ~ centeredadult*policy + f(time, model = "ar1")
nzunemploy.inla2 <- inla(formula2, data = nzunemploy, control.family = list(hyper = list(prec = list(initial = 15, fixed = TRUE))))
summary(nzunemploy.inla2)
library(nlme)
nzunemploy.gls = gls(youth ~ centeredadult*policy, correlation = corAR1(form=~1), data = nzunemploy)
summary(nzunemploy.gls)
sessionInfo()
