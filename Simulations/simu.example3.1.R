# This script replicates simulation results for Data Generating Process (DGP) A.
# Estimations are performed in parallel to enhance efficiency.

## Add your working directory
proot <- c("~/Dropbox/Papers - In progress/2-steps Estimation/2-steps M estimator", 
           "~/2StageInference")
root  <- sapply(proot, dir.exists)
root  <- proot[root][1]
setwd(root)

## Libraries
rm(list = ls())
library(dplyr)
library(MASS)
library(doParallel)
library(ggplot2)
library(patchwork)

## Loading C++ functions
Rcpp::sourceCpp("Simulations/Codes/tsesimu.cppfunctions.cpp")

## This function performs a single iteration of the Monte Carlo simulation.
fsim   <- function(N, k = 1e3){# k is kappa in the paper
  ## Data
  Z    <- cbind(1, runif(N)) #instrument
  b0   <- Z%*%c(0.1, 0.8) #Prob(d = 1)
  Y    <- theta0*b0 + runif(N, -1, 1) #outcome
  d    <- 1*(b0 > runif(N))
  
  ## First stage
  step1   <- summary(lm(d ~ -1 + Z))
  gamahat <- step1$coefficients[,"Estimate"]
  covgamh <- solve(crossprod(Z)) %*% crossprod(Z*step1$residuals) %*% solve(crossprod(Z))
  bhat    <- c(Z %*% gamahat)
  
  ## Second stage
  # Estimation
  step2   <- summary(lm(Y ~ -1 + bhat))
  thetah  <- step2$coefficients[,"Estimate"]
  
  # Naive standard errors
  cothe   <- solve(sum(bhat*bhat))^2 * crossprod(bhat*step2$residuals) 
  sdnai   <- sqrt(diag(cothe))*sqrt(N)
  
  # psi
  bhats   <- Z %*% t(mvrnorm(n = k, mu = gamahat, Sigma = covgamh))
  An      <- 2*sum(bhat^2)/N
  En      <- sapply(1:k, function(x) 2*sum(bhats[,x]*(bhat - bhats[,x])*thetah)/sqrt(N))
  Vh      <- 4*(step2$sigma^2)*sum(bhat^2)/N
  psi     <- sapply(1:k, function(x) (sqrt(Vh)*rnorm(1) + En[x])/An)
  
  # Standard errors
  sderr   <- sqrt(Vh + var(En)*k/(k - 1))/abs(An)
  
  list(delta  = sqrt(N)*(thetah - theta0),
       psi    = psi,
       sdnai  = sdnai,
       sderr  = sderr)
}

# Simulations
RNGkind("L'Ecuyer-CMRG")
set.seed(1234)

theta0 <- 1
nvec   <- c(1000, 2000)
sim    <- 1e4
k      <- 5e3 

out1   <- mclapply(1:sim, function(...) fsim(nvec[1], k), mc.cores = 15L)
out2   <- mclapply(1:sim, function(...) fsim(nvec[2], k), mc.cores = 15L)
out    <- list(out1, out2)

# Plot CDFs
# The values at which the CDF are evaluated
t      <- seq(-8, 8, 0.001) 
lt     <- length(t)

fdata   <- function(s){
  Ns    <- nvec[s]
  datas <- out[[s]]
  Nsl   <- ifelse(Ns < 1000, Ns, paste0(Ns/1000, ",000"))
  F0    <- fcdf(sapply(datas, function(x) x$delta), t, sim, lt)
  Fh    <- rowMeans(do.call(cbind, mclapply(datas, function(x) fcdf(x$psi, t, k, lt), mc.cores = 15L)))
  H1    <- rowMeans(do.call(cbind, mclapply(datas, function(x) pnorm(t, 0, x$sdnai), mc.cores = 11L)))
  H2    <- rowMeans(do.call(cbind, mclapply(datas, function(x) pnorm(t, 0, x$sderr), mc.cores = 11L)))
  
  data.frame(s    = s,
             t    = t,
             F0   = F0,
             Fh   = Fh,
             H1   = H1,
             H2   = H2,
             N    = Ns,
             parm = "theta0",
             type = rep(s, lt))  %>% 
    mutate(type   = factor(type, labels = expression(paste(sqrt(n)*(hat(theta)[n] - theta[0]), " with ", n, " = ", !!Nsl)))) 
}

dataplot <- bind_rows(lapply(1:2, fdata))

# Distances
(outdist <- dataplot %>% group_by(parm, N) %>% 
    summarise(across(c("H1", "H2", "Fh"), ~ sum(abs(.x - F0))*(max(t) - min(t))/length(t))))

# CDF of Delta
graph         <- list()
for (k in 1:length(nvec)) {
  distF       <- format(round(outdist$Fh[k], 3), nsmall = 3)
  distH1      <- format(round(outdist$H1[k], 3), nsmall = 3)
  distH2      <- format(round(outdist$H2[k], 3), nsmall = 3)
  graph[[k]]  <- ggplot(dataplot %>% filter(t >= -4, t <= 4, s == k), aes(x = t, y = F0)) + theme_bw() + 
    geom_line(aes(col = "A", lty = "A")) + 
    geom_line(aes(y = Fh, col = "B", lty = "B")) + 
    geom_line(aes(y = H1, col = "C", lty = "C")) + 
    geom_line(aes(y = H2, col = "D", lty = "D")) + 
    facet_wrap(. ~ type, scales = "free", labeller = label_parsed) +
    xlab("") + 
    ylab(ifelse(k >= 2, "", "Probability")) + 
    scale_colour_manual("", values = c("black", "#e01b89", "#1b1be0", "#1b1be0"), 
                        labels = c("A" = expr(F[0]), 
                                   "B" = expr(paste(hat(F)[n], " (", !!distF, ")")), 
                                   "C" = expr(paste(hat(H)[1], " (", !!distH1, ")")),
                                   "D" = expr(paste(hat(H)[2], " (", !!distH2, ")")))) + 
    scale_linetype_manual("", values = c("A" = 1, "B" = 2, "C" = 3, "D" = 4),
                          labels = c("A" = expr(F[0]), 
                                     "B" = expr(paste(hat(F)[n], " (", !!distF, ")")), 
                                     "C" = expr(paste(hat(H)[1], " (", !!distH1, ")")),
                                     "D" = expr(paste(hat(H)[2], " (", !!distH2, ")")))) +
    theme(legend.position = c(0.22, 0.8), legend.text.align = 0,
          legend.background = element_rect(fill='transparent'),
          legend.text = element_text(size = 8),
          plot.margin = margin(0, 0, 0, 0, "cm")) 
  if(k >= 2){
    graph[[k]]  <- graph[[k]] +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
  }
}

# export figures
ggsave("simu:example3.1.pdf", path = "Simulations", plot = wrap_plots(graph, ncol = length(nvec)), 
       device = "pdf", width = 5, height = 2.8)
