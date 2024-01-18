#' This file estimates the peer effect model using the data set that is built from 
#' the file 1.smoking.data.R

rm(list = ls())
library(haven)
library(dplyr)
library(PartialNetwork)
library(AER)
library(Matrix)
library(MASS)
library(ggplot2)
library(fastDummies)
set.seed(1234)

## Add your working directory
proot <- c("~/AHdata/FastFood")
root  <- sapply(proot, dir.exists)
root  <- proot[root][1]
setwd(root)

## data
load(file = "DataFastFood.2Stage.rda")
G         <- lapply(1:nsc, function(x) GM[[x]] + GF[[x]])
Gnorm     <- norm.network(G)
M         <- nsc
nvec      <- sapply(G, nrow)
n         <- sum(nvec)
ncum      <- c(0, cumsum(nvec))
k         <- 1e4

## Explanatory variables
exp.var   <- c("Female", "Age", "Hispanic", "Grade7_8", "Grade9_10", "Grade11_12", "Black", "Asian", 
               "Otherrace", "WithBoth",  "Wkyallowance", "MoEduLHigh", "MoEduSomCO", "MoEduUniv", 
               "MoEduMiss",  "FaEduLHigh", "FaEduSomCO", "FaEduUniv", "FaEduMiss", "MJobProf", 
               "MJobOhter", "MJobMiss", "FJobProf", "FJobOhter", "FJobMiss")

X         <- cbind("Intercept" = 1, as.matrix(datafull[,exp.var]))
y         <- datafull$Fastfood
Gy        <- peer.avg(Gnorm, y)
GX        <- peer.avg(Gnorm, X[,-1]); colnames(GX) <- paste0("G:", colnames(X)[-1])

## Functions
# Instruments construction
fZ        <- function(k, GX){
  nc      <- ncol(GX)
  out     <- matrix(0, nrow(GX), nc*(k - 1))
  out[,1:nc] <- peer.avg(Gnorm, GX)
  for (x in 3:k) {
    out[,(nc*(x - 2) + 1):(nc*(x - 1))] <- peer.avg(Gnorm, out[,(nc*(x - 3) + 1):(nc*(x - 2))])
  }
  out
  apply(out, 2, function(x) x/sd(x)); colnames(out) <- paste0("V", 1:ncol(out))
  out
}

##################### OLS
(ols      <- summary(lm(y ~ -1 + Gy + X + GX)))
write.csv(ols$coefficients, file = "ols.estim.csv")

##################### Classical method with GGX as instrument
GGX       <- peer.avg(Gnorm, GX); colnames(GGX) <- paste0("GG:", colnames(X)[-1])
(estim    <- summary(ivreg(y ~ -1 + Gy + X + GX|-1 + X + GX + GGX), diagnostic = TRUE))
write.csv(estim$coefficients, file = "gmm.estim.csv")
write.csv(estim$diagnostics, file = "gmm.diagnostics.csv")

##################### Optimal GMM
Ey        <- peer.avg(lapply(1:M, function(m) solve(diag(nvec[m]) - estim$coefficients["Gy","Estimate"]*Gnorm[[m]])),
                      cbind(X, GX)%*%estim$coefficients[-1,"Estimate"])
EGy       <- peer.avg(Gnorm, Ey)   
(estopt   <- summary(ivreg(y ~ -1 + Gy + X + GX|-1 + X + GX + EGy), diagnostic = TRUE))
write.csv(estopt$coefficients, file = "ogmm.estim.csv")
write.csv(estopt$diagnostics, file = "ogmm.diagnostics.csv")

##################### Many instruments
GpX       <- fZ(11, GX[,])
dim(GpX)
R         <- cbind("Gy" = Gy, X, GX)
U         <- cbind(X, GX, GpX)
summary(ivreg(y ~ -1 + Gy + X + GX|-1 + X + GX + GpX), diagnostic = TRUE)

########## First stage
step11   <- summary(lm(y ~ -1 + U))
step12   <- summary(lm(Gy ~ -1 + U))
gamh1    <- step11$coefficients[,"Estimate"]
gamh2    <- step12$coefficients[,"Estimate"]
gamh     <- c(gamh1, gamh2)
covgamh  <- bdiag(solve(crossprod(U)), solve(crossprod(U))) 
S0       <- crossprod(cbind(U*step11$residuals, U*step12$residuals))
covgamh  <- covgamh%*%S0%*%covgamh
yhat     <- c(U %*% gamh1)
Gyhat    <- c(U %*% gamh2)

########## Second stage
# Estimation
Rh      <- cbind(Gyhat, X, GX)
step2   <- summary(lm(yhat ~ - 1 + Rh))
(thetah <- step2$coefficients[,"Estimate"])

# Naive standard errors
cothe   <- solve(crossprod(Rh)) %*% crossprod(Rh*step2$residuals) %*% solve(crossprod(Rh))
snai1   <- sqrt(diag(cothe))*sqrt(n)
(snai1/sqrt(n))

# psi
gamhs   <- mvrnorm(n = k, mu = gamh, Sigma = covgamh)
yhats   <- sapply(1:k, function(x) U %*% gamhs[x, 1:ncol(U)])
Gyhats  <- sapply(1:k, function(x) U %*% gamhs[x, (1 + ncol(U)):(2*ncol(U))])
An      <- 2*crossprod(Rh)/n
En      <- sapply(1:k, function(x) 2*crossprod(cbind(Gyhats[,x], Rh[,-1]), yhats[,x] - cbind(Gyhats[,x], Rh[,-1]) %*% thetah))/sqrt(n)
psi     <- sapply(1:k, function(x) solve(An, En[,x]))

# Standard errors
snai2   <- sqrt(diag(solve(An) %*% (var(t(En))*k/(k - 1)) %*% solve(An)))
snai2/sqrt(n)

# Bias reduction
(thetas <- thetah - solve(An, rowMeans(En))/sqrt(n))
En      <- sapply(1:k, function(x) 2*crossprod(cbind(Gyhats[,x], Rh[,-1]), yhats[,x] - cbind(Gyhats[,x], Rh[,-1]) %*% thetas))/sqrt(n)
Omega   <- rowMeans(En)
psis    <- sapply(1:k, function(x) solve(An, En[,x] - Omega)) 

snai2s  <- sqrt(diag(solve(An) %*% (var(t(En))*k/(k - 1)) %*% solve(An)))
snai2s/sqrt(n)

estMI   <- data.frame(coef   = thetah,
                      sder1  = snai1/sqrt(n),
                      sder2  = snai2/sqrt(n),
                      ICInf  = apply(thetah - psi/sqrt(n), 1, function(x) quantile(x, prob = 0.025)),
                      ICSup  = apply(thetah - psi/sqrt(n), 1, function(x) quantile(x, prob = 0.975)),
                      coefs  = thetas,
                      sder2  = snai2s/sqrt(n),
                      ICInfs = apply(thetas - psis/sqrt(n), 1, function(x) quantile(x, prob = 0.025)),
                      ICSups = apply(thetas - psis/sqrt(n), 1, function(x) quantile(x, prob = 0.975)))

write.csv(estMI, file = "gmmlarge.estim.csv")

datapl  <- data.frame(Type     = 1,
                      Model    = 5:1,
                      Estimate = c(ols$coefficients[1,"Estimate"],
                                   estim$coefficients[1,"Estimate"],
                                   estopt$coefficients[1,"Estimate"],
                                   thetah[1], thetas[1]),
                      Sderr    = c(ols$coefficients[1,"Std. Error"],
                                   estim$coefficients[1,"Std. Error"],
                                   estopt$coefficients[1,"Std. Error"],
                                   snai2[1]/sqrt(n), snai2s[1]/sqrt(n)),
                      ICInf    = c(ols$coefficients[1,"Estimate"] + ols$coefficients[1,"Std. Error"]*qnorm(0.025),
                                   estim$coefficients[1,"Estimate"] + estim$coefficients[1,"Std. Error"]*qnorm(0.025),
                                   estopt$coefficients[1,"Estimate"] + estopt$coefficients[1,"Std. Error"]*qnorm(0.025),
                                   estMI$ICInf[1], estMI$ICInfs[1]),
                      ICSup    = c(ols$coefficients[1,"Estimate"] + ols$coefficients[1,"Std. Error"]*qnorm(0.975),
                                   estim$coefficients[1,"Estimate"] + estim$coefficients[1,"Std. Error"]*qnorm(0.975),
                                   estopt$coefficients[1,"Estimate"] + estopt$coefficients[1,"Std. Error"]*qnorm(0.975),
                                   estMI$ICSup[1], estMI$ICSups[1])) %>%
  mutate(ICInfNorm = Estimate + Sderr*qnorm(0.025), ICISupNorm = Estimate + Sderr*qnorm(0.975))

##################### ADDING SCHOOL FIXE EFFECTS AS DUMMIES
X         <- cbind(as.matrix(dummy_cols(datafull["SSCID2"])[,-1]), as.matrix(datafull[,exp.var]))
GX        <- peer.avg(Gnorm, X[,-(1:16)]); colnames(GX) <- paste0("G:", colnames(X)[-(1:16)])

##################### OLS
(ols      <- summary(lm(y ~ -1 + Gy + X + GX)))
write.csv(ols$coefficients, file = "ols.fe.estim.csv")

##################### Classical method with GGX as instrument
GGX       <- peer.avg(Gnorm, GX[,]); colnames(GGX) <- paste0("GG:", colnames(X)[-(1:16)])
(estim    <- summary(ivreg(y ~ -1 + Gy + X + GX|-1 + X + GX + GGX), diagnostic = TRUE))
write.csv(estim$coefficients, file = "gmm.fe.estim.csv")
write.csv(estim$diagnostics, file = "gmm.fe.diagnostics.csv")

##################### Optimal GMM
Ey        <- peer.avg(lapply(1:M, function(m) solve(diag(nvec[m]) - estim$coefficients["Gy","Estimate"]*Gnorm[[m]])),
                      cbind(X, GX)%*%estim$coefficients[-1,"Estimate"])
EGy       <- peer.avg(Gnorm, Ey)   
(estopt   <- summary(ivreg(y ~ -1 + Gy + X + GX|-1 + X + GX + EGy), diagnostic = TRUE))
write.csv(estopt$coefficients, file = "ogmm.fe.estim.csv")
write.csv(estopt$diagnostics, file = "ogmm.fe.diagnostics.csv")

##################### Many instruments
########## Function to construct the number of instrument
GpX       <- fZ(11, GX)
dim(GpX)
R         <- cbind("Gy" = Gy, X, GX)
U         <- cbind(X, GX, GpX)

summary(ivreg(y ~ -1 + Gy + X + GX|-1 + X + GX + GpX), diagnostic = TRUE)
########## First stage
step11   <- summary(lm(y ~ -1 + U))
step12   <- summary(lm(Gy ~ -1 + U))
gamh1    <- step11$coefficients[,"Estimate"]
gamh2    <- step12$coefficients[,"Estimate"]
gamh     <- c(gamh1, gamh2)
covgamh  <- bdiag(solve(crossprod(U)), solve(crossprod(U))) 
S0       <- crossprod(cbind(U*step11$residuals, U*step12$residuals))
covgamh  <- covgamh%*%S0%*%covgamh
yhat     <- c(U %*% gamh1)
Gyhat    <- c(U %*% gamh2)

########## Second stage
# Estimation
Rh       <- cbind(Gyhat, X, GX)
step2    <- summary(lm(yhat ~ - 1 + Rh))
(thetah  <- step2$coefficients[,"Estimate"])

# Naive standard errors
cothe   <- solve(crossprod(Rh)) %*% crossprod(Rh*step2$residuals) %*% solve(crossprod(Rh))
snai1   <- sqrt(diag(cothe))*sqrt(n)
(snai1/sqrt(n))

# psi
gamhs   <- mvrnorm(n = k, mu = gamh, Sigma = covgamh)
yhats   <- sapply(1:k, function(x) U %*% gamhs[x, 1:ncol(U)])
Gyhats  <- sapply(1:k, function(x) U %*% gamhs[x, (1 + ncol(U)):(2*ncol(U))])
An      <- 2*crossprod(Rh)/n
En      <- sapply(1:k, function(x) 2*crossprod(cbind(Gyhats[,x], Rh[,-1]), yhats[,x] - cbind(Gyhats[,x], Rh[,-1]) %*% thetah))/sqrt(n)
psi     <- sapply(1:k, function(x) solve(An, En[,x]))

# Standard errors
snai2   <- sqrt(diag(solve(An) %*% (var(t(En))*k/(k - 1)) %*% solve(An)))
snai2/sqrt(n)

# Bias reduction
(thetas <- thetah - solve(An, rowMeans(En))/sqrt(n))
En      <- sapply(1:k, function(x) 2*crossprod(cbind(Gyhats[,x], Rh[,-1]), yhats[,x] - cbind(Gyhats[,x], Rh[,-1]) %*% thetas))/sqrt(n)
Omega   <- rowMeans(En)
psis    <- sapply(1:k, function(x) solve(An, En[,x] - Omega)) 

snai2s  <- sqrt(diag(solve(An) %*% (var(t(En))*k/(k - 1)) %*% solve(An)))
snai2s/sqrt(n)

estMI   <- data.frame(coef   = thetah,
                      sder1  = snai1/sqrt(n),
                      sder2  = snai2/sqrt(n),
                      ICInf  = apply(thetah - psi/sqrt(n), 1, function(x) quantile(x, prob = 0.025)),
                      ICSup  = apply(thetah - psi/sqrt(n), 1, function(x) quantile(x, prob = 0.975)),
                      coefs  = thetas,
                      sder2  = snai2s/sqrt(n),
                      ICInfs = apply(thetas - psis/sqrt(n), 1, function(x) quantile(x, prob = 0.025)),
                      ICSups = apply(thetas - psis/sqrt(n), 1, function(x) quantile(x, prob = 0.975)))

write.csv(estMI, file = "gmmlarge.fe.estim.csv")

datapl.f <- data.frame(Type     = 2,
                       Model    = 5:1,
                       Estimate = c(ols$coefficients[1,"Estimate"],
                                    estim$coefficients[1,"Estimate"],
                                    estopt$coefficients[1,"Estimate"],
                                    thetah[1], thetas[1]),
                       Sderr    = c(ols$coefficients[1,"Std. Error"],
                                    estim$coefficients[1,"Std. Error"],
                                    estopt$coefficients[1,"Std. Error"],
                                    snai2[1]/sqrt(n), snai2s[1]/sqrt(n)),
                       ICInf    = c(ols$coefficients[1,"Estimate"] + ols$coefficients[1,"Std. Error"]*qnorm(0.025),
                                    estim$coefficients[1,"Estimate"] + estim$coefficients[1,"Std. Error"]*qnorm(0.025),
                                    estopt$coefficients[1,"Estimate"] + estopt$coefficients[1,"Std. Error"]*qnorm(0.025),
                                    estMI$ICInf[1], estMI$ICInfs[1]),
                       ICSup    = c(ols$coefficients[1,"Estimate"] + ols$coefficients[1,"Std. Error"]*qnorm(0.975),
                                    estim$coefficients[1,"Estimate"] + estim$coefficients[1,"Std. Error"]*qnorm(0.975),
                                    estopt$coefficients[1,"Estimate"] + estopt$coefficients[1,"Std. Error"]*qnorm(0.975),
                                    estMI$ICSup[1], estMI$ICSups[1])) %>%
  mutate(ICInfNorm = Estimate + Sderr*qnorm(0.025), ICISupNorm = Estimate + Sderr*qnorm(0.975))

##################### PLOT CIs
dataplot <- rbind(datapl, datapl.f) %>% 
  mutate(Type  = factor(Type, labels = c(expression(paste("Without ", "Fixed ", "Effects")),
                                        expression(paste("With ", "Fixed ", "Effects")))),
         Model = factor(Model, labels = c("DIV-MI", "IV-MI", "OIV", "CIV", "OLS")))

(graph   <- ggplot(dataplot, aes(y = Model, x = Estimate)) + theme_bw() +
    geom_errorbarh(aes(xmin = ICInf, xmax = ICSup),
                   height = 0.1, linetype = 1, linewidth = 0.6,               
                   position = position_dodge(.1), colour = "#555555") +
    geom_point(shape = 8, colour = "#EE2090", size = 1.5) +
    geom_text(aes(x = Estimate, y = Model, label = format(round(Estimate, 3), nsmall = 2)), nudge_y = 0.25, size = 3, colour = "#EE2090") + 
    geom_text(aes(x = ICInf, y = Model, label = format(round(ICInf, 3), nsmall = 2)), nudge_x = -0.12, size = 3, colour = "#555555") +
    geom_text(aes(x = ICSup, y = Model, label = format(round(ICSup, 3), nsmall = 2)), nudge_x = 0.12, size = 3, colour = "#555555") +
    xlim(x = c(-0.8, 0.65)) +
    facet_wrap(. ~ Type, labeller = label_parsed, dir="v", ncol = 2))

ggsave("app:plot.pdf", plot = graph, device = "pdf", width = 7, height = 3)
