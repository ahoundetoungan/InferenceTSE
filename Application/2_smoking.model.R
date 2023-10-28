#' This file estimates the peer effect model using the dataset that is built from 
#' the file 1_smoking.data.R

##################### Libraries
rm(list = ls())
library(dplyr)
library(ggplot2)
library(MASS)
library(Matrix)
library(PartialNetwork)

##################### load data
load(file = "AHdata.rda")

exp.var  <- va.names[-length(va.names)]
Gnorm    <- norm.network(G)
X        <- cbind(Intercept = 1, as.matrix(mydata[,exp.var]))
y        <- mydata$smoke
Gy       <- peer.avg(Gnorm, y)
GX       <- peer.avg(Gnorm, X)
GGX      <- peer.avg(Gnorm, GX)
var.net  <- colnames(Xlogit)
var.net  <- var.net[!(var.net %in% c("female.i", "female.j", "ID.i", "ID.j", "sschlcde.j"))]
M        <- nsch
nvec     <- sapply(G, nrow)
n        <- sum(nvec)

######################### Taking the network as observe
J          <- lapply(nvec, function(x) diag(x) - matrix(1/x, x, x))
JX         <- peer.avg(J, X[,-1])
JGX        <- peer.avg(J, GX[,-1])
JGGX       <- peer.avg(J, GGX[,-1])
W          <- solve(crossprod(as.matrix(cbind(JX, JGX, JGGX))))
gmm        <- smmSAR(formula      = y ~ -1 + X[,-1],
                     contextual    = TRUE,
                     dnetwork      = G,
                     fixed.effects = TRUE,
                     W             = W,
                     smm.ctr       = list(R = 1L, iv.power = 2, opt.tol = 1e-3, print = TRUE))
saveRDS(gmm, file = "gmm.RDS")
sgmm       <- summary(gmm)
write.csv(cbind(coef = sgmm$estimates, se = sqrt(diag(sgmm$cov))), file = "gmm.csv")

######################### Controlling for unmatched and non-nominated friends
# We first estimate the distribution of the network.
# We use a censored Poisson regression and a censored normal distribution. 
# The dependent variable is the degree. 
# As the number of friends is censored for female and male friends, 
# we run a regression for the number of declared female friends and another
# regression on the number of declared male friends.
# We include school fixed effects

########### Number of declared female friends
sel     <- (Xlogit[,"female.j"] == 1)
Nicum   <- c(0, cumsum(table(Xlogit[sel, "ID.i"])))
####### Poisson distribution
# Starting point
NetP1f  <- ENetPoi(rho = rep(0, M + length(var.net)), d = mydata$friendfall,
                   x = Xlogit[sel, var.net], bound = 5, Nicum = Nicum,
                   group = as.numeric(factor(mydata$sschlcde, labels = 0:(M - 1))) - 1,
                   n = n, M = M, maxit = 1e6, eps_f = 1e-11, eps_g = 1e-11)
saveRDS(NetP1f, file = "NetP1f.RDS")
# Optimization using Newton Raphson
NetP2f  <- NRNetPoi(rho = readRDS("NetP1f.RDS")$estimate, d = mydata$friendfall,
                    x = Xlogit[sel, var.net], bound = 5, Nicum = Nicum,
                    group = as.numeric(factor(mydata$sschlcde, labels = 0:(M - 1))) - 1,
                    n = n, M = M, maxit = 100, tol = 1e-4)
saveRDS(NetP2f, file = "NetP2f.RDS")

####### Normal regression
# Starting point 1
rho0     <- StartNetNor(readRDS("NetP1f.RDS")$estimate, d = mydata$friendfall,
                        x = Xlogit[sel, var.net], bound = 5, Nicum = Nicum,
                        group = as.numeric(factor(mydata$sschlcde, labels = 0:(M - 1))) - 1,
                        n = n, M = M)
# Starting point 2 (Better than the first one)
NetN1f  <- ENetNor(rho = rho0, d = mydata$friendfall,
                   x = Xlogit[sel, var.net], bound = 5, Nicum = Nicum,
                   group = as.numeric(factor(mydata$sschlcde, labels = 0:(M - 1))) - 1,
                   n = n, M = M, maxit = 1e6, eps_f = 1e-11, eps_g = 1e-11)
saveRDS(NetN1f, file = "NetN1f.RDS")
# Optimization using Newton Raphson
NetN2f  <- NRNetNor(rho = readRDS("NetN1f.RDS")$estimate, d = mydata$friendfall,
                    x = Xlogit[sel, var.net], bound = 5, Nicum = Nicum,
                    group = as.numeric(factor(mydata$sschlcde, labels = 0:(M - 1))) - 1,
                    n = n, M = M, maxit = 100, tol = 1e-4)
saveRDS(NetN2f, file = "NetN2f.RDS")

########### Number of declared male friends
sel     <- (Xlogit[,"female.j"] == 0)
Nicum   <- c(0, cumsum(table(Xlogit[sel, "ID.i"])))
####### Poisson distribution
# Starting point
NetP1m  <- ENetPoi(rho = readRDS("NetP1m.RDS")$estimate, d = mydata$friendmall,
                   x = Xlogit[sel, var.net], bound = 5, Nicum = Nicum,
                   group = as.numeric(factor(mydata$sschlcde, labels = 0:(M - 1))) - 1,
                   n = n, M = M, maxit = 1e6, eps_f = 1e-11, eps_g = 1e-11)
saveRDS(NetP1m, file = "NetP1m.RDS")
# Optimization using Newton Raphson
NetP2m  <- NRNetPoi(rho = readRDS("NetP1m.RDS")$estimate, d = mydata$friendmall,
                    x = Xlogit[sel, var.net], bound = 5, Nicum = Nicum,
                    group = as.numeric(factor(mydata$sschlcde, labels = 0:(M - 1))) - 1,
                    n = n, M = M, maxit = 100, tol = 1e-4)
saveRDS(NetP2m, file = "NetP2m.RDS")

####### Normal distribution
# Starting point 1
rho0    <- StartNetNor(readRDS("NetP1m.RDS")$estimate, d = mydata$friendmall,
                       x = Xlogit[sel, var.net], bound = 5, Nicum = Nicum,
                       group = as.numeric(factor(mydata$sschlcde, labels = 0:(M - 1))) - 1,
                       n = n, M = M)
# Starting point 2 (Better than the first one)
NetN1m  <- ENetNor(rho = rho0, d = mydata$friendmall,
                   x = Xlogit[sel, var.net], bound = 5, Nicum = Nicum,
                   group = as.numeric(factor(mydata$sschlcde, labels = 0:(M - 1))) - 1,
                   n = n, M = M, maxit = 1e6, eps_f = 1e-11, eps_g = 1e-11)
saveRDS(NetN1m, file = "NetN1m.RDS")
# Optimization using Newton Raphson
NetN2m  <- NRNetNor(rho = readRDS("NetN1m.RDS")$estimate, d = mydata$friendmall,
                    x = Xlogit[sel, var.net], bound = 5, Nicum = Nicum,
                    group = as.numeric(factor(mydata$sschlcde, labels = 0:(M - 1))) - 1,
                    n = n, M = M, maxit = 100, tol = 1e-4)
saveRDS(NetN2m, file = "NetN2m.RDS")
sumNet      <- tibble("Variable"     = c(paste0("school", 1:M), var.net, "log(se2)")) %>%
  mutate(`PCoef Female` = c(NetP2f$estimate, NA),
         `PSter Female` = c(sqrt(diag(NetP2f$covmat)), NA),
         `PStat Female` = `PCoef Female`/`PSter Female`,
         `PPval Female` = 2*(1 - pnorm(abs(`PStat Female`))),
         `NCoef Female` = NetN2f$estimate,
         `NSter Female` = sqrt(diag(NetN2f$covmat)),
         `NStat Female` = `NCoef Female`/`NSter Female`,
         `NPval Female` = 2*(1 - pnorm(abs(`NStat Female`))),
         `PCoef Male`   = c(NetP2m$estimate, NA),
         `PSter Male`   = c(sqrt(diag(NetP2m$covmat)), NA),
         `PStat Male`   = `PCoef Male`/`PSter Male`,
         `PPval Male`   = 2*(1 - pnorm(abs(`PStat Male`))),
         `NCoef Male`   = NetN2m$estimate,
         `NSter Male`   = sqrt(diag(NetN2m$covmat)),
         `NStat Male`   = `NCoef Male`/`NSter Male`,
         `NPval Male`   = 2*(1 - pnorm(abs(`NStat Male`))))
write.csv(sumNet, file = "summary.network.csv", row.names = FALSE)

########### Network distribution
IDi     <- Xlogit[,"ID.i"]
groupj  <- as.numeric(factor(Xlogit[,"sschlcde.j"], labels = 0:(M - 1))) - 1
femalej <- Xlogit[,"female.j"]
NX      <- Xlogit[,var.net]
G       <- mat.to.vec(G)
Gobtmis <- mat.to.vec(Gobtmis)
gc()

# Coefficient rho and joint covariance because the estimates of 
# rhofemale and rhomale are linked to the same students.
NetP2f  <- readRDS("NetP2f.RDS")
NetP2m  <- readRDS("NetP2m.RDS")
NetN2f  <- readRDS("NetN2f.RDS")
NetN2m  <- readRDS("NetN2m.RDS")

rhofemP <- NetP2f$estimate
rhomalP <- NetP2m$estimate
rhofemN <- NetN2f$estimate
rhomalN <- NetN2m$estimate

inHessP <- as.matrix(bdiag(solve(NetP2f$Hess), solve(NetP2m$Hess)))
cholCoP <- chol(tcrossprod(rbind(NetP2f$grad, NetP2m$grad)))
cholCoP <- inHessP %*% t(cholCoP)  #CholCoP %*% t(CholCoP) is the covariance matrix

inHessN <- as.matrix(bdiag(solve(NetN2f$Hess), solve(NetN2m$Hess)))
cholCoN <- chol(tcrossprod(rbind(NetN2f$grad, NetN2m$grad)))
cholCoN <- inHessN %*% t(cholCoN)  #CholCoN %*% t(CholCoN) is the covariance matrix

rm(list = c("NetP2f", "NetP2m", "NetN2f", "NetN2m"))

Kx      <- length(var.net)
nparmsP <- M + Kx
nparmsN <- M + Kx + 1

# A function that generates the network distribution
fNet    <- function(rhofemale, rhomale, nparms, NX, groupj, femalej, G, Gobtmis, 
                    Kx, M, nvec, cholcovrho = NULL){
  rhof  <- rhofemale
  rhom  <- rhomale
  if(!is.null(cholcovrho)){
    return (vec.to.mat(estimPsim(rho = c(rhofemale, rhofemale), cholcovrho = cholcovrho, X = NX, groupj = groupj, 
                                 femalej = femalej, G = G, Gobs = Gobtmis, Kx = Kx, M = M, nparms = nparms), 
                       nvec))
  }
  vec.to.mat(estimP(rhofemale = rhof, rhomale = rhom, X = NX, groupj = groupj, 
                    femalej = femalej, G = G, Gobs = Gobtmis, 
                    Kx = Kx, M = M), nvec)
}
dnetPoi <- fNet(rhofemale = rhofemP, rhomale = rhomalP, nparms = nparmsP, NX = NX, 
                groupj = groupj, femalej = femalej, G = G, Gobtmis = Gobtmis, 
                Kx = Kx, M = M, nvec = nvec)
dnetNor <- fNet(rhofemale = rhofemN, rhomale = rhomalN, nparms = nparmsM, NX = NX, 
                groupj = groupj, femalej = femalej, G = G, Gobtmis = Gobtmis, 
                Kx = Kx, M = M, nvec = nvec)

########### SMM model With the Poisson distribution
# constructing W = (Z'Z)^{-1} taking into account school fixed effects
# Centering matrix
set.seed(123)
J          <- lapply(nvec, function(x) diag(x) - matrix(1/x, x, x))
JX         <- peer.avg(J, X[,-1])
WPoi       <- lapply(1:1e3, function(x){
  cat("Iteration: ", x, "/1000", sep = "", "\n")
  Gsimnorm <- norm.network(sim.network(dnetPoi))
  GXsim    <- peer.avg(Gsimnorm, X)
  GGXsim   <- peer.avg(Gsimnorm, GXsim)
  JGXsim   <- peer.avg(J, GXsim[,-1])
  JGGXsim  <- peer.avg(J, GGXsim[,-1])
  crossprod(as.matrix(cbind(JX, JGXsim, JGGXsim)))})
WPoi       <- solve(Reduce('+', WPoi)/(N*1e3))
saveRDS(WPoi, file = "WPoi.RDS")

WPoi   <- readRDS("WPoi.RDS")
set.seed(123)
R      <- 50
sgmmP  <- smmSAR(formula       = y ~ -1 + X[,-1],
                 contextual    = TRUE,
                 dnetwork      = dnetPoi,
                 fixed.effects = TRUE,
                 W             = WPoi,
                 smm.ctr       = list(R = R, iv.power = 2, opt.tol = 1e-3, print = TRUE))

saveRDS(sgmmP, file = paste0("sgmmP.R=", R, ".RDS"))

###### Inference
R        <- 50
sgmmP    <- readRDS(paste0("sgmmP.R=", R, ".RDS"))
set.seed(123)
# Without taking into account the first stage
summary(sgmmP) 
# Our method
fNetargs <- list(rhofemale = rhofemP, rhomale = rhomalP, nparms = nparmsP, NX = NX, 
                 groupj = groupj, femalej = femalej, G = G, Gobtmis = Gobtmis, Kx = Kx, 
                 M = M, nvec = nvec, cholcovrho = cholCoP)  
ssgmmP   <- simudist(sgmmP, dnetwork = dnetPoi, .fun = fNet, .args = fNetargs, sim = 1000L)
saveRDS(ssgmmP, file = paste0("ssgmmP.R=", R, ".RDS"))
ssgmmP   <- readRDS(paste0("ssgmmP.R=", R, ".RDS"))
outsgmmP <- t(apply(ssgmmP, 1, function(x) c(SdErr2 = sd(x), quantile(x, probs = c(0.025, 0.975)))))
write.csv(cbind(coef = sgmmP$estimates, SdErr1 = sqrt(diag(summary(sgmmP)$cov)), outsgmmP), 
          file = "sgmmP.csv")

########### SMM model With the Normal distribution
# constructing W = (Z'Z)^{-1} taking into account school fixed effects
# Centering matrix
set.seed(123)
J          <- lapply(nvec, function(x) diag(x) - matrix(1/x, x, x))
JX         <- peer.avg(J, X[,-1])
WNor       <- lapply(1:1e3, function(x){
  cat("Iteration: ", x, "/1000", sep = "", "\n")
  Gsimnorm <- norm.network(sim.network(dnetNor))
  GXsim    <- peer.avg(Gsimnorm, X)
  GGXsim   <- peer.avg(Gsimnorm, GXsim)
  JGXsim   <- peer.avg(J, GXsim[,-1])
  JGGXsim  <- peer.avg(J, GGXsim[,-1])
  crossprod(as.matrix(cbind(JX, JGXsim, JGGXsim)))})
WNor       <- solve(Reduce('+', WNor)/(N*1e3))
saveRDS(WNor, file = "WNor.RDS")

WNor  <- readRDS("WNor.RDS")
set.seed(123)
R     <- 50
sgmmN <- smmSAR(formula       = y ~ -1 + X[,-1],
                contextual    = TRUE,
                dnetwork      = dnetNor,
                fixed.effects = TRUE,
                W             = WNor,
                smm.ctr       = list(R = R, iv.power = 2, opt.tol = 1e-3, print = TRUE))

saveRDS(sgmmN, file = paste0("sgmmN.R=", R, ".RDS"))

###### Inference
R        <- 50
sgmmN    <- readRDS(paste0("sgmmN.R=", R, ".RDS"))
set.seed(123)
# Without taking into account the first stage
summary(sgmmN)
# Our method
fNetargs <- list(rhofemale = rhofemN, rhomale = rhomalN, nparms = nparmsN, NX = NX,
                 groupj = groupj, femalej = femalej, G = G, Gobtmis = Gobtmis, Kx = Kx,
                 M = M, nvec = nvec, cholcovrho = cholCoN)
ssgmmN   <- simudist(sgmmN, dnetwork = dnetNor, .fun = fNet, .args = fNetargs, sim = 1000)
saveRDS(ssgmmN, file = paste0("ssgmmN.R=", R, ".RDS"))
ssgmmN   <- readRDS(paste0("ssgmmN.R=", R, ".RDS"))
outsgmmN <- t(apply(ssgmmN, 1, function(x) c(SdErr2 = sd(x), quantile(x, probs = c(0.025, 0.975)))))
write.csv(cbind(coef = sgmmN$estimates, SdErr1 = sqrt(diag(summary(sgmmN)$cov)), outsgmmN), 
          file = "sgmmN.csv")
