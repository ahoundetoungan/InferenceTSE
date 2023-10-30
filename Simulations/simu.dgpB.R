# This script replicates simulation results for Data Generating Process (DGP) B.
# Estimations are performed in parallel to enhance efficiency.

## Libraries
rm(list = ls())
library(MASS)
library(splines)
library(dplyr)
library(doParallel)
library(ggplot2)

## Loading R and C++ functions
source("tsesimu.rfunctions.R")
Rcpp::sourceCpp("tsesimu.cppfunctions.cpp")

## True values of the parameters
theta0    <- c(0.5, 2)

## The values at which the CDF are evaluated
t         <- seq(-15, 15, 0.01)
lt        <- length(t)

## Sample sizes
nvec      <- c(250, 500, 1e3, 2e3)

## This function performs a single iteration of the Monte Carlo simulation.
fsim      <- function(N, k){# k is kappa in the paper
  ## Data
  nz      <- c(1, 0.95, 0.90, 0.85)[nvec == N]
  n       <- round(N^nz)
  Z       <- sort(runif(N, 0, 10))
  p       <- sin(Z*acos(-1))^2
  d       <- ifelse(p > runif(N), 1, 0)
  lambda  <- exp(theta0[1] + theta0[2]*p)
  y       <- rpois(N, lambda)
  
  ## First stage
  # Spline regression
  knots   <- seq(0, 9.5, 0.5)
  nknots  <- length(knots)
  nZ      <- nknots + 3
  Zbs     <- bs(Z, degree = 3L, knots = knots); colnames(Zbs) <- paste0("Z", 1:nZ)
  data    <- data.frame(d = d, Zbs, select = 0); data$select[sample(1:N, n)] <- 1
  step1   <- lm(paste0("d ~ -1 + ", paste0("Z", 1:nZ, collapse = "+")), 
                data = data %>% filter(select == 1))
  Zbs1    <- as.matrix(data %>% filter(select == 1) %>% dplyr::select(!!paste0("Z", 1:nZ)))
  phat    <- c(Zbs %*% step1$coefficients)
  # # check the fit
  # plot(Z, p, type = "line", col = "blue")
  # lines(Z, phat, col = "red")
  # cor(phat, p)
  
  ## Second stage
  # Estimation
  data    <- data.frame(y = y, phat = phat)
  step2   <- summary(glm(y ~ phat, data = data, family = poisson(link = "log")))
  thetah  <- step2$coefficients[,"Estimate"]
  
  # Naive standard errors
  snai1   <- step2$coefficients[,"Std. Error"]*sqrt(N)
  
  # psi
  tmp     <- solve(crossprod(Zbs1))
  varBs   <- tmp %*% crossprod(Zbs1*step1$residuals) %*% tmp
  phats   <- Zbs %*% t(mvrnorm(n = k, mu = step1$coefficients, Sigma = varBs))
  mhat    <- cbind(1, phat)
  mhats   <- lapply(1:k, function(x) cbind(1, phats[,x]))
  lhat    <- c(exp(mhat %*% thetah)) #lambdahat
  lhats   <- sapply(mhats, function(x) c(exp(x %*% thetah))) #lambdahat_s
  An      <- -crossprod(mhat, mhat*lhat)/N
  En      <- sapply(1:k, function(x) crossprod(mhats[[x]], lhats[,x] - lhat)/sqrt(N))
  Vn      <- lapply(mhats, function(x) crossprod(x, x*lhat)/N)
  psi     <- sapply(1:k, function(x) solve(An, t(chol(Vn[[x]])) %*% rnorm(2) + En[,x]))
  
  # Standard errors
  snai2   <- apply(psi, 1, sd)
  
  # Quantiles using each inference method
  quant1  <- thetah + t(sapply(snai1/sqrt(N), function(x){
    x*qnorm(c(0.005, 0.025, 0.05, 0.950, 0.975,  0.995))})) 
  quant2  <- thetah + t(sapply(snai2/sqrt(N), function(x){
    x*qnorm(c(0.005, 0.025, 0.05, 0.950, 0.975,  0.995))})) 
  quant   <- thetah - t(apply(psi, 1, function(x){
    quantile(x/sqrt(N), probs = c(0.005, 0.025, 0.05, 0.950, 0.975,  0.995))})) 
  rownames(quant1)  <- names(thetah)
  rownames(quant2)  <- names(thetah)
  rownames(quant)   <- names(thetah)
  colnames(quant1)  <- c("0.5%", "2.5%", "5%", "95%", "97.5%",  "99.5%")
  colnames(quant2)  <- c("0.5%", "2.5%", "5%", "95%", "97.5%",  "99.5%")
  colnames(quant)   <- c("0.5%", "2.5%", "5%", "95%", "97.5%",  "99.5%")
  
  list(delta  = sqrt(N)*(thetah - theta0),
       psi    = apply(psi, 1, function(x) fcdf(x, t, k, lt)),
       snai1  = snai1,
       snai2  = snai2,
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
out3 <- mclapply(1:sim, function(...) fsim(nvec[3], k), mc.cores = 15L)
out4 <- mclapply(1:sim, function(...) fsim(nvec[4], k), mc.cores = 15L)
out  <- list(out1, out2, out3, out4)

# Tableau of proportions of CI contenting the true value
coverage <- lapply(c("N = 250" = 1, "N = 500" = 2, "N = 1000" = 3, "N = 2000" = 4), function(s){
  f.coverage.vec(out[[s]], theta0)})
print(coverage)

# Plot CDFs
fdata    <- function(s){
  Ns     <- nvec[s]
  datas  <- out[[s]]
  Nsl    <- ifelse(Ns < 1000, Ns, paste0(Ns/1000, ",000"))
  F0     <- apply(sapply(datas, function(x) x$delta), 1, function(x) fcdf(x, t, sim, lt))
  Fh     <- sapply(1:2, function(x1) rowMeans(sapply(datas, function(x2) x2$psi[,x1])))
  Norn1  <- sapply(1:2, function(x1) pnorm(t, 0, mean(sapply(datas, function(x2) x2$snai1[x1]))))
  Norn2  <- sapply(1:2, function(x1) pnorm(t, 0, median(sapply(datas, function(x2) x2$snai2[x1]))))
  
  data.frame(t    = rep(t, 2),
             F0   = c(F0[,1], F0[,2]),
             Fh   = c(Fh[,1], Fh[,2]),
             Nn1  = c(Norn1[,1], Norn1[,2]),
             Nn2  = c(Norn2[,1], Norn2[,2]),
             N    = rep(s, 2*lt),
             type = rep((2*s -1):(2*s), each = lt))  %>% 
    mutate(type   = factor(type, labels = c(expression(paste(sqrt(n)*(hat(theta)[n][","][1] - theta[0][","][1]), " with n = ", !!Nsl)), 
                                            expression(paste(sqrt(n)*(hat(theta)[n][","][2] - theta[0][","][2]), " with n = ", !!Nsl))))) 
}

dataplot <- bind_rows(lapply(1:4, fdata))

(graph   <- ggplot(dataplot %>% filter(N %in% 1:4), aes(x = t, y = F0)) + theme_bw() + 
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
    # theme(legend.position = c(0.06, 0.9), legend.text.align = 0,
    #       legend.background = element_rect(fill='transparent'),
    #       legend.text = element_text(size = 8)) +
    theme(legend.position = "none") + 
    facet_wrap(. ~ type, scales = "free", labeller = label_parsed, dir="v", ncol = 4))

# export figures
ggsave("simu:poisson.pdf", plot = graph, device = "pdf", width = 10, height = 5)

