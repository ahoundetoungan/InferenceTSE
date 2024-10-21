# This script replicates simulation results for Data Generating Process (DGP) A.
# Estimations are performed in parallel to enhance efficiency.

## Add your working directory
proot <- c("~/Dropbox/Academy/1.Papers/2-steps Estimation/2-steps M estimator",
           "~/2StageInference")
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
fsim   <- function(N, k = 1e3){# k is kappa in the paper
  ## Data
  Z    <- runif(N) 
  e    <- runif(N, -1, 1) 
  d    <- 1*((2*Z) > (e + 1.2)) 
  y    <- theta0*d + e 
  
  ## First stage
  Zint    <- cbind(1, Z)
  step11  <- summary(lm(y ~ -1 + Zint))
  step12  <- summary(lm(d ~ -1 + Zint))
  gamh1   <- step11$coefficients[,"Estimate"]
  gamh2   <- step12$coefficients[,"Estimate"]
  gamh    <- c(gamh1, gamh2)
  covgamh <- bdiag(solve(crossprod(Zint)), solve(crossprod(Zint))) 
  S0      <- crossprod(cbind(Zint*step11$residuals, Zint*step12$residuals))
  covgamh <- covgamh%*%S0%*%covgamh
  yhat    <- c(Zint %*% gamh1)
  dhat    <- c(Zint %*% gamh2)
  
  ## Second stage
  # Estimation
  step2   <- summary(lm(yhat ~ - 1 + dhat))
  thetah  <- step2$coefficients[,"Estimate"]
  
  # # Naive standard errors
  # cothe   <- solve(crossprod(dhat)) %*% crossprod(dhat*step2$residuals) %*% solve(crossprod(dhat))
  # sdnai   <- sqrt(diag(cothe))*sqrt(N)
  
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
  thetahs <- thetah - mean(En)/(An*sqrt(N))
  En      <- sapply(1:k, function(x) 2*sum(dhats[,x]*(yhats[,x] - dhats[,x] * thetahs))/sqrt(N))
  Omega   <- mean(En)
  psis    <- sapply(1:k, function(x) (En[x] - Omega)/An)
  
  list(delta   = sqrt(N)*(thetah - theta0),
       deltas  = sqrt(N)*(thetahs - theta0),
       psi     = psi,
       psis    = psis,
       sderr   = sderr,
       thetah  = thetah,
       thetahs = thetahs)
}

# Simulations
RNGkind("L'Ecuyer-CMRG")
set.seed(1234)

theta0 <- 1
nvec   <- c(250, 500, 1e3, 2e3)
sim    <- 1e4
k      <- 1e3

out1   <- mclapply(1:sim, function(...) fsim(nvec[1], k), mc.cores = 15L)
out2   <- mclapply(1:sim, function(...) fsim(nvec[2], k), mc.cores = 15L)
out3   <- mclapply(1:sim, function(...) fsim(nvec[3], k), mc.cores = 12L)
out4   <- mclapply(1:sim, function(...) fsim(nvec[4], k), mc.cores = 10L)
out    <- list(out1, out2, out3, out4)
saveRDS(out, file = "Simulations/simu.iv.RDS")

# Plot CDFs
# The values at which the CDF are evaluated
t       <- seq(-10, 10, 0.001)
lt      <- length(t)

out     <- readRDS("Simulations/simu.iv.RDS")

fdata   <- function(s){
  Ns    <- nvec[s]
  datas <- out[[s]]
  Nsl   <- ifelse(Ns < 1000, Ns, paste0(Ns/1000, ",000"))
  F0    <- fcdf(sapply(datas, function(x) x$delta), t, sim, lt)
  F0s   <- fcdf(sapply(datas, function(x) x$deltas), t, sim, lt)
  Fh    <- rowMeans(do.call(cbind, mclapply(datas, function(x) fcdf(x$psi, t, k, lt), mc.cores = 15L)))
  Fhs   <- rowMeans(do.call(cbind, mclapply(datas, function(x) fcdf(x$psis, t, k, lt), mc.cores = 15L)))
  H     <- rowMeans(do.call(cbind, mclapply(datas, function(x) pnorm(t, 0, x$sderr), mc.cores = 11L)))
  
  
  data.frame(s    = s,
             parm = "theta0",
             type = rep(s, lt),
             N    = Ns,
             t    = t,
             F0   = F0,
             F0s  = F0s,
             Fh   = Fh,
             Fhs  = Fhs,
             H    = H)  %>% 
    mutate(type   = factor(type, labels = expression(paste(sqrt(n)*(hat(theta)[n] - theta[0]), ", ", n, " = ", !!Nsl))),
           types  = factor(type, labels = expression(paste(sqrt(n)*(theta[paste(n, ",", kappa)]^"*" - theta[0]), ", ", n, " = ", !!Nsl))))
}

dataplot <- bind_rows(lapply(1:4, fdata))


# Distances
(outdist <- dataplot %>% group_by(s, parm, N) %>% 
    summarise(across(c("H", "Fh"), ~ sum(abs(.x - F0))*(max(t) - min(t))/length(t)))  %>%
    inner_join(dataplot %>% group_by(s, parm, N) %>% 
                 summarise(Fhs = sum(abs(Fhs - F0s))*(max(t) - min(t))/length(t)), by  = c("s", "parm", "N")))

# CDF of Delta
graph1        <- list()
for (k in 1:length(nvec)) {
  distF       <- format(round(outdist$Fh[k], 3), nsmall = 3)
  distH       <- format(round(outdist$H[k], 3), nsmall = 3)
  graph1[[k]] <- ggplot(dataplot %>% filter(t >= -5, t <= 5, s == k), aes(x = t, y = F0)) + theme_bw() + 
    geom_line(aes(col = "A", lty = "A")) + 
    geom_line(aes(y = Fh, col = "B", lty = "B")) + 
    geom_line(aes(y = H, col = "C", lty = "C")) + 
    facet_wrap(. ~ type, scales = "free", labeller = label_parsed) +
    xlab("") + 
    ylab(ifelse(k >= 2, "", "Probability")) + 
    scale_colour_manual("", values = c("black", "#e01b89", "#1b1be0"), 
                        labels = c("A" = expr(F[0]), 
                                   "B" = expr(paste(hat(F)[n], " (", !!distF, ")")), 
                                   "C" = expr(paste(hat(H)[n], " (", !!distH, ")")))) + 
    scale_linetype_manual("", values = c("A" = 1, "B" = 2, "C" = 3),
                          labels = c("A" = expr(F[0]), 
                                     "B" = expr(paste(hat(F)[n], " (", !!distF, ")")), 
                                     "C" = expr(paste(hat(H)[n], " (", !!distH, ")")))) +
    theme(legend.position = c(0.2, 0.75), legend.text.align = 0,
          legend.background = element_rect(fill='transparent'),
          legend.text = element_text(size = 8),
          plot.margin = margin(0, 0, 0, 0, "cm")) 
  if(k >= 2){
    graph1[[k]] <- graph1[[k]] +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
  }
}

# CDF of Deltastart
graph2        <- list()
for (k in 1:length(nvec)) {
  distF       <- format(round(outdist$Fhs[k], 3), nsmall = 3)
  graph2[[k]] <- ggplot(dataplot %>% filter(t >= -5, t <= 5, s == k), aes(x = t, y = F0s)) + theme_bw() + 
    geom_line(aes(col = "A", lty = "A")) + 
    geom_line(aes(y = Fhs, col = "B", lty = "B")) + 
    facet_wrap(. ~ types, scales = "free", labeller = label_parsed) +
    xlab("") + 
    ylab(ifelse(k >= 2, "", "Probability")) + 
    scale_colour_manual("", values = c("black", "#e01b89"), 
                        labels = c("A" = expr(F[0]^"*"), 
                                   "B" = expr(paste(hat(F)[n]^"*", " (", !!distF, ")")))) + 
    scale_linetype_manual("", values = c("A" = 1, "B" = 2),
                          labels = c("A" = expr(F[0]^"*"), 
                                     "B" = expr(paste(hat(F)[n]^"*", " (", !!distF, ")")))) +
    theme(legend.position = c(0.2, 0.9), legend.text.align = 0,
          legend.background = element_rect(fill='transparent'),
          legend.text = element_text(size = 8),
          plot.margin = margin(0, 0, 0, 0, "cm")) 
  if(k >= 2){
    graph2[[k]] <- graph2[[k]] +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
  }
}

# export figures
ggsave("simu.iv.pdf", path = "Simulations", plot = wrap_plots(graph1, ncol = length(nvec)), 
       device = "pdf", width = 10, height = 2.5)
ggsave("simu.ivR.pdf", path = "Simulations", plot = wrap_plots(graph2, ncol = length(nvec)), 
       device = "pdf", width = 10, height = 2.5)

# Coverage rates
fcrate  <- function(s){
  Ns    <- nvec[s]
  datas <- out[[s]]
  IC1   <- sapply(datas, function(x) x$thetah + qnorm(c(0.025, 0.05, 0.95, 0.975), 0, x$sderr/sqrt(Ns)))
  IC    <- sapply(datas, function(x) quantile(x$thetah - x$psi/sqrt(Ns), probs = c(0.025, 0.05, 0.95, 0.975), na.rm = TRUE))
  ICs   <- sapply(datas, function(x) quantile(x$thetahs - x$psis/sqrt(Ns), probs = c(0.025, 0.05, 0.95, 0.975), na.rm = TRUE))
  
  crate <- c(mean(theta0 >= IC1[1,] & theta0 <= IC1[4,]), mean(theta0 <= IC1[3,]), mean(theta0 >= IC1[2,]),
             mean(theta0 >= IC[1,] & theta0 <= IC[4,]), mean(theta0 <= IC[3,]), mean(theta0 >= IC[2,]),
             mean(theta0 >= ICs[1,] & theta0 <= ICs[4,]), mean(theta0 <= ICs[3,]), mean(theta0 >= ICs[2,]))
  crate <- as.data.frame(crate)
  colnames(crate) <- paste0("N=", Ns)
  rownames(crate) <- paste(rep(c("Naive", "Simu", "Simu-Debiased"), each = 3), 
                           rep(c("2-Sides", "Lower 1-Side", "Upper 1 Side"), 3))
  as.data.frame(crate)
}

(crate  <- do.call(cbind, lapply(1:4, fcrate)))
write.csv(crate, file = "Simulations/simu.iv.crate.csv")

# Bias
fbias <- function(s){
  Ns    <- nvec[s]
  datas <- out[[s]]
  theh  <- sapply(datas, function(x) x$thetah)
  thes  <- sapply(datas, function(x) x$thetahs)
  data.frame(parm = "theta0",
             N    = Ns,
             the0 = theta0,
             theh = theh,
             thes = thes)
}

(databias <- bind_rows(lapply(1:4, fbias)) %>% group_by(parm, N) %>%
    summarise(Mean1 = mean(theh), BIAS1 = mean(theh - the0), SD1 = sd(theh), 
              RMSE1 = sqrt(mean((theh - the0)^2)), MAE1 = mean(abs(theh - the0)),
              Mean2 = mean(thes), BIAS2 = mean(thes - the0), SD2 = sd(thes), 
              RMSE2 = sqrt(mean((thes - the0)^2)), MAE2 = mean(abs(thes - the0))))
write.csv(databias, file = "Simulations/simu.iv.bias.csv", row.names = FALSE)