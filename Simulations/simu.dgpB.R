# This script replicates simulation results for Data Generating Process (DGP) B.
# Estimations are performed in parallel to enhance efficiency.

## Add your working directory
proot <- c("~/Dropbox/2-steps M estimator", 
           "~/Dropbox/2StageInference")
root  <- sapply(proot, dir.exists)
root  <- proot[root][1]
setwd(root)

## Libraries
rm(list = ls())
library(dplyr)
library(MASS)
library(Matrix)
library(doParallel)
library(ggplot2)
library(patchwork)

## Loading C++ functions
Rcpp::sourceCpp("Simulations/Codes/tsesimu.cppfunctions.cpp")

## This function performs a single iteration of the Monte Carlo simulation.
fsim   <- function(N, Kz, k = 1e3){# k is kappa in the paper
  ## Data
  Z    <- cbind(1, matrix(runif(N*Kz, 0, 0.2), N))
  e    <- runif(N, -1, 1) 
  d    <- 1*((0.2*Z[,1] + Z[,2] + Z[,3] + Z[,4] + Z[,5]) > (0.5*(e + 1.2))) 
  y    <- theta0*d + e 
  
  ## First stage
  Zint    <- Z
  step11  <- summary(lm(y ~ -1 + Zint))
  step12  <- summary(lm(d ~ -1 + Zint))
  gamh1   <- step11$coefficients[,"Estimate"]
  gamh2   <- step12$coefficients[,"Estimate"]
  gamh    <- c(gamh1, gamh2)
  covgamh <- bdiag(solve(crossprod(Zint)), solve(crossprod(Zint))) 
  S0      <- crossprod(cbind(Zint*step11$residuals, Zint*step12$residuals))
  covgamh <- covgamh%*%S0%*%covgamh
  yhat     <- c(Zint %*% gamh1)
  dhat     <- c(Zint %*% gamh2)
  
  ## Second stage
  # Estimation
  step2   <- summary(lm(yhat ~ - 1 + dhat))
  thetah  <- step2$coefficients[,"Estimate"]
  
  # Naive standard errors
  cothe   <- solve(crossprod(dhat)) %*% crossprod(dhat*step2$residuals) %*% solve(crossprod(dhat))
  sdnai   <- sqrt(diag(cothe))*sqrt(N)
  
  # psi
  gamhs   <- mvrnorm(n = k, mu = gamh, Sigma = covgamh)
  yhats   <- sapply(1:k, function(x) c(Zint %*% gamhs[x, 1:ncol(Zint)]))
  dhats   <- sapply(1:k, function(x) c(Zint %*% gamhs[x, (1 + ncol(Zint)):(2*ncol(Zint))]))
  An      <- 2*sum(dhat^2)/N
  En      <- sapply(1:k, function(x) 2*sum(dhats[,x]*(yhats[,x] - dhats[,x] * thetah))/sqrt(N))
  psi     <- sapply(1:k, function(x) En[x]/An)
  
  # Standard errors
  sderr   <- sqrt((var(En)*k/(k - 1))/An^2)
  
  # Bias reduction
  thetahs  <- thetah - mean(En)/(An*sqrt(N))
  En      <- sapply(1:k, function(x) 2*sum(dhats[,x]*(yhats[,x] - dhats[,x] * thetahs))/sqrt(N))
  Omega   <- mean(En)
  psis    <- sapply(1:k, function(x) (En[x] - Omega)/An)
  
  list(delta   = sqrt(N)*(thetah - theta0),
       deltas  = sqrt(N)*(thetahs - theta0),
       psi     = psi,
       psis    = psis,
       sdnai   = sdnai,
       sderr   = sderr,
       thetah  = thetah,
       thetahs  = thetahs)
}

# Simulations
RNGkind("L'Ecuyer-CMRG")
set.seed(1234)

theta0 <- 1
nvec   <- c(250, 500, 1e3, 2e3)
knvec  <- lapply(nvec, function(x) round(c(2, 4)*sqrt(x)))
sim    <- 1e4
k      <- 1e3

out1   <- mclapply(1:sim, function(...) fsim(nvec[1], knvec[[1]][1], k), mc.cores = 15L)
out2   <- mclapply(1:sim, function(...) fsim(nvec[1], knvec[[1]][2], k), mc.cores = 15L)
out3   <- mclapply(1:sim, function(...) fsim(nvec[2], knvec[[2]][1], k), mc.cores = 15L)
out4   <- mclapply(1:sim, function(...) fsim(nvec[2], knvec[[2]][2], k), mc.cores = 15L)
out5   <- mclapply(1:sim, function(...) fsim(nvec[3], knvec[[3]][1], k), mc.cores = 12L)
out6   <- mclapply(1:sim, function(...) fsim(nvec[3], knvec[[3]][2], k), mc.cores = 12L)
out7   <- mclapply(1:sim, function(...) fsim(nvec[4], knvec[[4]][1], k), mc.cores = 10L)
out8   <- mclapply(1:sim, function(...) fsim(nvec[4], knvec[[4]][2], k), mc.cores = 10L)
out    <- list(out1, out2, out3, out4, out5, out6, out7, out8)
saveRDS(out, file = "Simulations/simu:ivlarge.RDS")

# Plot CDFs
# The values at which the CDF are evaluated
t       <- seq(-20, 20, 0.001)
lt      <- length(t)

out     <- readRDS("Simulations/simu:ivlarge.RDS")
fdata   <- function(s){
  s1    <- ceiling(s/2)
  Ns    <- nvec[s1]
  Kzs   <- knvec[[s1]][ifelse((2*s1) == s, 2, 1)]
  datas <- out[[s]]
  Nsl   <- ifelse(Ns < 1000, Ns, paste0(Ns/1000, ",000"))
  F0    <- fcdf(sapply(datas, function(x) x$delta), t, sim, lt)
  F0s   <- fcdf(sapply(datas, function(x) x$deltas), t, sim, lt)
  Fh    <- rowMeans(do.call(cbind, mclapply(datas, function(x) fcdf(x$psi, t, k, lt), mc.cores = 15L)))
  Fhs   <- rowMeans(do.call(cbind, mclapply(datas, function(x) fcdf(x$psis, t, k, lt), mc.cores = 15L)))
  H1    <- rowMeans(do.call(cbind, mclapply(datas, function(x) pnorm(t, 0, x$sdnai), mc.cores = 11L)))
  H2    <- rowMeans(do.call(cbind, mclapply(datas, function(x) pnorm(t, 0, x$sderr), mc.cores = 11L)))
  
  data.frame(s    = s,
             parm = "theta0",
             type = rep(s1, lt),
             N    = Ns,
             Kz   = Kzs,
             t    = t,
             F0   = F0,
             F0s  = F0s,
             Fh   = Fh,
             Fhs  = Fhs,
             H1   = H1,
             H2   = H2)  %>% 
    mutate(type   = factor(type, labels = expression(paste(sqrt(n)*(hat(theta)[n] - theta[0]), 
                                                           ", ", n, " = ", !!Nsl, ", ", k[n], " = ", !!Kzs))),
           typec  = factor(type, labels = expression(paste(sqrt(n)*(hat(theta)[n]^"*" - theta[0]), 
                                                           ", ", n, " = ", !!Nsl, ", ", k[n], " = ", !!Kzs))))
}

dataplot <- bind_rows(lapply(1:8, fdata))

# Distance
(outdist <- dataplot %>% group_by(s, parm, N, Kz) %>% 
  summarise(across(c("H1", "H2", "Fh"), ~ sum(abs(.x - F0))*(max(t) - min(t))/length(t))) %>%
  inner_join(dataplot %>% group_by(s, parm, N, Kz) %>% 
               summarise(Fhs = sum(abs(Fhs - F0s))*(max(t) - min(t))/length(t)), by  = c("s", "parm", "N", "Kz")))

# CDF of Delta
graph1        <- list()
for (k in 1:(2*length(nvec))) {
  distF       <- format(round(outdist$Fh[k], 3), nsmall = 3)
  distH1      <- format(round(outdist$H1[k], 3), nsmall = 3)
  distH2      <- format(round(outdist$H2[k], 3), nsmall = 3)
  graph1[[k]] <- ggplot(dataplot %>% filter(t >= -12, t <= 10, s == k), aes(x = t, y = F0)) + theme_bw() + 
    geom_line(aes(col = "A", lty = "A")) + 
    geom_line(aes(y = Fh, col = "B", lty = "B")) + 
    geom_line(aes(y = H1, col = "C", lty = "C")) + 
    geom_line(aes(y = H2, col = "D", lty = "D")) + 
    facet_wrap(. ~ type, scales = "free", labeller = label_parsed) +
    xlab("") + 
    ylab(ifelse(k >= 3, "", "Probability")) + 
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
    theme(legend.position = c(0.2, 0.75), legend.text.align = 0,
          legend.background = element_rect(fill='transparent'),
          legend.text = element_text(size = 8),
          plot.margin = margin(0, 0, 0, 0, "cm")) 
  if(k >= 3){
    graph1[[k]] <- graph1[[k]] +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
  }
}

# CDF of Deltastart
graph2        <- list()
for (k in 1:(2*length(nvec))) {
  distF       <- format(round(outdist$Fhs[k], 3), nsmall = 3)
  graph2[[k]] <- ggplot(dataplot %>% filter(t >= -10, t <= 10, s == k), aes(x = t, y = F0s)) + theme_bw() + 
    geom_line(aes(col = "A", lty = "A")) + 
    geom_line(aes(y = Fhs, col = "B", lty = "B")) + 
    facet_wrap(. ~ type, scales = "free", labeller = label_parsed) +
    xlab("") + 
    ylab(ifelse(k >= 3, "", "Probability")) + 
    scale_colour_manual("", values = c("black", "#e01b89", "#1b1be0", "#1b1be0"), 
                        labels = c("A" = expr(F[0]^"*"), 
                                   "B" = expr(paste(hat(F)[n]^"*", " (", !!distF, ")")))) + 
    scale_linetype_manual("", values = c("A" = 1, "B" = 2, "C" = 3, "D" = 4),
                          labels = c("A" = expr(F[0]^"*"), 
                                     "B" = expr(paste(hat(F)[n]^"*", " (", !!distF, ")")))) +
    theme(legend.position = c(0.2, 0.9), legend.text.align = 0,
          legend.background = element_rect(fill='transparent'),
          legend.text = element_text(size = 8),
          plot.margin = margin(0, 0, 0, 0, "cm")) 
  if(k >= 3){
    graph2[[k]] <- graph2[[k]] +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
  }
}

# export figures
ggsave("simu:ivlarge.pdf", path = "Simulations", plot = wrap_plots(graph1, ncol = length(nvec), byrow = FALSE), 
       device = "pdf", width = 10, height = 5)
ggsave("simu:ivlargeR.pdf", path = "Simulations", plot = wrap_plots(graph2, ncol = length(nvec), byrow = FALSE), 
       device = "pdf", width = 10, height = 5)

# Bias
fbias <- function(s){
  s1    <- ceiling(s/2)
  Ns    <- nvec[s1]
  Kzs   <- knvec[[s1]][ifelse((2*s1) == s, 2, 1)]
  datas <- out[[s]]
  theh  <- sapply(datas, function(x) x$thetah)
  thes  <- sapply(datas, function(x) x$thetahs)
  data.frame(parm = "theta0",
             N    = Ns,
             Kz   = Kzs,
             the0 = theta0,
             theh = theh,
             thes = thes)
}

(databias <- bind_rows(lapply(1:8, fbias)) %>% group_by(parm, N, Kz) %>%
  summarise(Mean1 = mean(theh), BIAS1 = mean(theh - the0), SD1 = sd(theh), 
            RMSE1 = sqrt(mean((theh - the0)^2)), MAE1 = mean(abs(theh - the0)),
            Mean2 = mean(thes), BIAS2 = mean(thes - the0), SD2 = sd(thes), 
            RMSE2 = sqrt(mean((thes - the0)^2)), MAE2 = mean(abs(thes - the0))))
write.csv(databias, file = "Simulations/simu:ivlarge.bias.csv", row.names = FALSE)