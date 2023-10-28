# This script replicates simulation results for Data Generating Process (DGP) A.
# Estimations are performed in parallel to enhance efficiency.

## Libraries
rm(list = ls())
library(dplyr)
library(MASS)
library(doParallel)
library(ggplot2)

## Loading R and C++ functions
source("tsesimu.rfunctions.R")
Rcpp::sourceCpp("tsesimu.cppfunctions.cpp")

## True values of the parameters
theta0    <- c(0, 1, -0.5)
MTE0      <- theta0[2] + 2*theta0[3]*0.5

## The values at which the CDF are evaluated
t         <- seq(-15, 15, 0.1)
lt        <- length(t)

## Sample sizes
nvec <- c(250, 500, 1e3, 2e3)

## This function performs a single iteration of the Monte Carlo simulation.
fsim   <- function(N, k = 1e3){# k is kappa in the paper
  ## Data
  Z    <- cbind(1, runif(N)) #Z,i ∼ Uniform[0,0.8]
  V    <- runif(N) #Vi|Zi ∼ Uniform[0,1]
  U0   <- runif(N, -1, 1) #U0i|Zi ,Vi ∼ Uniform[−1,1]
  U1   <- runif(N, -0.5, 1.5 - 2*V) #U1i |Zi,Vi ∼ Uniform[−0.5,1.5−2Vi ].
  Y0   <- U0 #Yi(0) = U0i
  Y1   <- 0.5 + U1 #Yi (1) = 0.5+U1i
  Ts   <- 1*(c(Z %*% c(0.1, 0.7)) >= V) #Selection equation
  Y    <- Ts*Y1 + (1 - Ts)*Y0 #Yi = Ti Yi (1)+(1−Ti )Yi (0)
  
  ## First stage
  step1   <- summary(lm(Ts ~ -1 + Z))
  bhat    <- step1$coefficients[,"Estimate"]
  covbhat <- solve(crossprod(Z)) %*% crossprod(Z*step1$residuals) %*% solve(crossprod(Z))
  phat    <- c(Z %*% bhat)
  
  ## Second stage
  # Estimation
  step2   <- summary(lm(Y ~ phat + I(phat^2)))
  thetah  <- step2$coefficients[,"Estimate"]
  MTEh    <- c(thetah[2] + thetah[3]); names(MTEh) <- "MTE"
  
  # Naive standard errors
  Xh      <- cbind(1, phat, phat^2)
  cothe   <- solve(crossprod(Xh)) %*% crossprod(Xh*step2$residuals) %*% solve(crossprod(Xh))
  snai1T  <- sqrt(diag(cothe))*sqrt(N)
  snai1M  <- c(t(c(0, 1, 1)) %*% cothe %*% c(0, 1, 1)*sqrt(N)) #Using Delta method
  
  # psi
  phats   <- Z %*% t(mvrnorm(n = k, mu = bhat, Sigma = covbhat))
  mhat    <- cbind(1, phat, phat^2)
  mhats   <- lapply(1:k, function(x) cbind(1, phats[,x], phats[,x]^2))
  An      <- -2*crossprod(mhat)/N
  En      <- sapply(1:k, function(x) 2*crossprod(mhats[[x]], (mhat - mhats[[x]]) %*% thetah)/sqrt(N))
  Vn      <- lapply(mhats, function(x) 4*crossprod(x*step2$residuals)/N)
  psiT    <- sapply(1:k, function(x) solve(An, t(chol(Vn[[x]])) %*% rnorm(3) + En[,x]))
  psiM    <- psiT[2,] + 2*psiT[3,]*0.5
  
  # Standard errors
  snai2T  <- apply(psiT, 1, sd)
  snai2M  <- sd(psiM)
  
  # Quantiles using each inference method
  quant1  <- c(thetah, MTEh) + t(sapply(c(snai1T, MTE = snai1M)/sqrt(N), function(x){
    x*qnorm(c(0.005, 0.025, 0.05, 0.950, 0.975,  0.995))})) 
  quant2  <- c(thetah, MTEh) + t(sapply(c(snai2T, MTE = snai2M)/sqrt(N), function(x){
    x*qnorm(c(0.005, 0.025, 0.05, 0.950, 0.975,  0.995))})) 
  quant   <- c(thetah, MTEh) - t(apply(rbind(psiT, MTE = psiM), 1, function(x){
    quantile(x/sqrt(N), probs = c(0.005, 0.025, 0.05, 0.950, 0.975,  0.995))})) 
  rownames(quant1)  <- c(names(thetah), "MTE")
  rownames(quant2)  <- c(names(thetah), "MTE")
  rownames(quant)   <- c(names(thetah), "MTE")
  colnames(quant1)  <- c("0.5%", "2.5%", "5%", "95%", "97.5%",  "99.5%")
  colnames(quant2)  <- c("0.5%", "2.5%", "5%", "95%", "97.5%",  "99.5%")
  colnames(quant)   <- c("0.5%", "2.5%", "5%", "95%", "97.5%",  "99.5%")
  
  list(delta  = sqrt(N)*c(thetah - theta0, MTE = MTEh - MTE0),
       psi    = apply(rbind(psiT, psiM), 1, function(x) fcdf(x, t, k, lt)),
       snai1  = c(snai1T, MTE = snai1M),
       snai2  = c(snai2T, MTE = snai2M),
       quant1 = quant1,
       quant2 = quant2,
       quant  = quant)
}

# Simulations
set.seed(1234)
sim  <- 1e4
k    <- 1e3
out1 <- mclapply(1:sim, function(...) fsim(nvec[1], k), mc.cores = 15L)
out2 <- mclapply(1:sim, function(...) fsim(nvec[2], k), mc.cores = 15L)
out3 <- mclapply(1:sim, function(...) fsim(nvec[3], k), mc.cores = 12L)
out4 <- mclapply(1:sim, function(...) fsim(nvec[4], k), mc.cores = 10L)
out  <- list(out1, out2, out3, out4)

# Tableau of proportions of CI contenting the true value
coverage <- lapply(c("N = 250" = 1, "N = 500" = 2, "N = 1000" = 3, "N = 2000" = 4), function(s){
  f.coverage.vec(out[[s]], c(theta0, MTE0))})
print(coverage)

# Plot CDFs
fdata    <- function(s){
  Ns     <- nvec[s]
  datas  <- out[[s]]
  Nsl    <- ifelse(Ns < 1000, Ns, paste0(Ns/1000, ",000"))
  F0     <- apply(sapply(datas, function(x) x$delta), 1, function(x) fcdf(x, t, sim, lt))
  Fh     <- sapply(1:4, function(x1) rowMeans(sapply(datas, function(x2) x2$psi[,x1])))
  Norn1  <- sapply(1:4, function(x1) pnorm(t, 0, mean(sapply(datas, function(x2) x2$snai1[x1]))))
  Norn2  <- sapply(1:4, function(x1) pnorm(t, 0, median(sapply(datas, function(x2) x2$snai2[x1]))))
  
  data.frame(t    = rep(t, 4),
             F0   = c(F0[,1], F0[,2], F0[,3], F0[,4]),
             Fh   = c(Fh[,1], Fh[,2], Fh[,3], Fh[,4]),
             Nn1  = c(Norn1[,1], Norn1[,2], Norn1[,3], Norn1[,4]),
             Nn2  = c(Norn2[,1], Norn2[,2], Norn2[,3], Norn2[,4]),
             N    = rep(s, 4*lt),
             type = rep((4*s - 3):(4*s), each = lt),
             parm = rep(1:4, each = lt))  %>% 
    mutate(type   = factor(type, labels = c(expression(paste(sqrt(n)*(hat(theta)[n][","][1] - theta[0][","][1]), " with n = ", !!Nsl)),
                                            expression(paste(sqrt(n)*(hat(theta)[n][","][2] - theta[0][","][2]), " with n = ", !!Nsl)), 
                                            expression(paste(sqrt(n)*(hat(theta)[n][","][3] - theta[0][","][3]), " with n = ", !!Nsl)),
                                            expression(paste(sqrt(n)*(hat(tau)["MTE"] - tau["MTE"]), " with n = ", !!Nsl))))) 
}

dataplot <- bind_rows(lapply(1:4, fdata))

(graph   <- ggplot(dataplot %>% filter(N %in% 1:4, parm == 4), aes(x = t, y = F0)) + theme_bw() + 
    geom_line(aes(col = "A", lty = "A")) + 
    geom_line(aes(y = Fh, col = "B", lty = "B")) + 
    geom_line(aes(y = Nn1, col = "C", lty = "C")) + 
    geom_line(aes(y = Nn2, col = "D", lty = "D")) + 
    facet_wrap(. ~ type, scales = "free", labeller = label_parsed) +
    # xlab(expression(paste(sqrt(n),"(",hat(theta)[n], '-', theta[0],")"))) +
    xlab("") +
    ylab("Probability") + 
    scale_colour_manual("", values = c("black", "#e01b89", "#1b1be0", "#1b1be0"), 
                        labels = c("A" = "True distribution", "B" = "Estimation", "C" = 'Normal 1', "D" = 'Normal 2')) +
    scale_linetype_manual("", values = c("A" = 1, "B" = 2, "C" = 3, "D" = 4),
                          labels = c("A" = "True distribution", "B" = "Estimation", "C" = 'Normal 1', "D" = 'Normal 2')) +
    theme(legend.position = c(0.06, 0.75), legend.text.align = 0,
          legend.background = element_rect(fill='transparent'),
          legend.text = element_text(size = 8)) +
    facet_wrap(. ~ type, scales = "free", labeller = label_parsed, dir="v", ncol = 4))

# export figures
ggsave("simu:iv.pdf", plot = graph, device = "pdf", path = "Simulations", width = 10, height = 2.5)
