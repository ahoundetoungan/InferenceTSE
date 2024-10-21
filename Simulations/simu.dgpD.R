# This script replicates simulation results for Data Generating Process (DGP) D.
# Estimations are performed in parallel to enhance efficiency.

## Add your working directory
proot <- c("~/Dropbox/Academy/1.Papers/2-steps Estimation/2-steps M estimator",
           "~/2StageInference")
root  <- sapply(proot, dir.exists)
root  <- proot[root][1]
setwd(root)

## Libraries
rm(list = ls())
library(copula)
library(MASS)
library(numDeriv)
library(Matrix)
library(dplyr)
library(doParallel)
library(ggplot2)
library(patchwork)

## Loading C++ functions
Rcpp::sourceCpp("Simulations/Codes/tsesimu.cppfunctions.cpp")

## This function performs a single iteration of the Monte Carlo simulation.
fsim   <- function(N, S = 2, k = 1e3){#S is the number of returns and k is kappa in the paper
  out         <- NULL
  tryCatch(
    expr = {
      ## Data
      # simulate the innovations of the garch
      csim    <- rCopula(N, claytonCopula(param = theta0, dim = S))
      z       <- apply(csim, 2, function(x) qt(x, df = nu)*sqrt(1 - 2/nu))
      
      # simulate the returns
      dt      <- filterdata(N, S, beta0, phi0, z)
      
      ## First stage
      th0     <- c(fbetat(c(phi0, beta0, nu)))
      st1     <- lapply(1:S, function(s){
        optim(par = th0, fn = fllhgarcht, N = N, dt = dt[,s],
              control = list(maxit = 1e6, abstol = 1e-13, reltol = 1e-13))})
      coeh    <- sapply(st1, function(s) fbeta(s$par))
      nuh     <- coeh[6,]
      zh      <- sapply(1:S, function(s){filterz(N, coeh[3:5, s], coeh[1:2, s], dt[,s])})
      pzh     <- sapply(1:S, function(s) pt(zh[,s]/sqrt(1 - 2/nuh[s]), df = nuh[s]))
      
      ## Joint covariance matrix of the first-stage estimator
      # Jacobian of the likelihood
      JAC     <- lapply(1:S, function(s){
        jacobian(func = function(x) c(fllhigarcht(x, N, dt[,s])), x = st1[[s]]$par,
                 method.args = list(eps = 1e-13, zero.tol = 1e-13, r = 4, v = 2))})
      # Hessian of the likelihood
      HES     <- lapply(1:S, function(s){
        hessian(func = function(x) c(-fllhgarcht(x, N, dt[,s])), x = st1[[s]]$par,
                method.args = list(eps = 1e-13, zero.tol = 1e-13, r = 4, v = 2))})
      bHES    <- solve(bdiag(HES))
      ## Covariance matrix
      MCOt    <- N * bHES %*% fHAC(do.call(cbind, JAC), N, r = 6*S) %*% t(bHES)
      
      ## Second stage
      # Estimation
      st2     <- optimize(fllhclayton, interval = c(-10, 3), u = pzh, dim = S, tol = 1e-11)
      lth     <- st2$minimum
      
      # An
      An      <- -sum(d2lclaytonpdf(u = pzh, ltheta = lth, dim = S))/N
      # Naive standard errors
      # sdnai   <- c(sqrt(fHAC(d1lclaytonpdf(u = pzh, ltheta = lth, dim = S), N, r = 1))/An)
      
      # psi
      coet    <- mvrnorm(n = k, mu = c(apply(coeh, 2, fbetat)), Sigma = MCOt)
      coes    <- lapply(1:S, function(s) t(apply(coet[,((s - 1)*6 + 1):(s*6)], 1, fbeta)))
      nus     <- sapply(coes, function(s) s[,6])
      zs      <- lapply(1:k, function(x){sapply(1:S, function(s){
        filterz(N, coes[[s]][x, 3:5], coes[[s]][x, 1:2], dt[,s])})})
      pzs     <- lapply(1:k, function(x) sapply(1:S, function(s) pt(zs[[x]][,s]/sqrt(1 - 2/nus[x, s]), df = nus[x, s])))
      En      <- sapply(1:k, function(x) sum(d1lclaytonpdf(u = pzs[[x]], ltheta = lth, dim = S)))/sqrt(N)
      # Vh      <- fHAC(d1lclaytonpdf(u = pzh, ltheta = lth, dim = S), N, r = 6*S)
      # psi     <- sapply(1:k, function(x) (sqrt(Vh)*rnorm(1) + En[x])/An)
      
      # Standard errors
      JACcdf  <- lapply(1:S, function(s){
        jacobian(func = function(x) c(fGigarcht(x, N, dt[,s])), x = st1[[s]]$par,
                 method.args = list(eps = 1e-13, zero.tol = 1e-13, r = 4, v = 2))})
      dudClay <- dud1lclaytonpdf(u = pzh, ltheta = lth, dim = S)
      AAn     <- rbind(cbind(-bdiag(HES)/N, 0),
                       c(apply(do.call(cbind, lapply(1:S, function(s) JACcdf[[s]]*dudClay[,s])), 2, mean), An))
      sderr   <- solve(AAn, fHAC(cbind(do.call(cbind, JAC), d1lclaytonpdf(u = pzh, ltheta = lth, dim = S)), N, r = 6*S + 1, st = 0)) %*% t(solve(AAn))
      sderr   <- sqrt(diag(sderr)[6*S + 1])
      psi     <- sderr*rnorm(k) + mean(En)/An
      
      # Bias reduction
      lths1   <- lth - mean(En)/(An*sqrt(N))
      lths2   <- lth - median(En)/(An*sqrt(N))
      
      # Inference: mean correction
      An      <- -sum(d2lclaytonpdf(u = pzh, ltheta = lths1, dim = S))/N
      dudClay <- dud1lclaytonpdf(u = pzh, ltheta = lths1, dim = S)
      AAn[6*S + 1,] <- c(-colMeans(do.call(cbind, lapply(1:S, function(s) JACcdf[[s]]*dudClay[,s]))), An)
      sderr1  <- solve(AAn, fHAC(cbind(do.call(cbind, JAC), d1lclaytonpdf(u = pzh, ltheta = lths1, dim = S)), N, r = 6*S + 1)) %*% t(solve(AAn))
      sderr1  <- sqrt(diag(sderr1)[6*S + 1])
      psis1   <- sderr1*rnorm(k)
      
      # # Inference: mean correction
      An      <- -sum(d2lclaytonpdf(u = pzh, ltheta = lths2, dim = S))/N
      dudClay <- dud1lclaytonpdf(u = pzh, ltheta = lths2, dim = S)
      AAn[6*S + 1,] <- c(-colMeans(do.call(cbind, lapply(1:S, function(s) JACcdf[[s]]*dudClay[,s]))), An)
      sderr2  <- solve(AAn, fHAC(cbind(do.call(cbind, JAC), d1lclaytonpdf(u = pzh, ltheta = lths2, dim = S)), N, r = 6*S + 1)) %*% t(solve(AAn))
      sderr2  <- sqrt(diag(sderr2)[6*S + 1])
      psis2   <- sderr*rnorm(k)
      
      out     <- list(delta    = sqrt(N)*(lth - log(theta0)),
                      deltas1  = sqrt(N)*(lths1 - log(theta0)),
                      deltas2  = sqrt(N)*(lths2 - log(theta0)),
                      psi      = psi,
                      psis1    = psis1,
                      psis2    = psis2,
                      sderr    = sderr,
                      thetah   = exp(lth),
                      thetahs1 = exp(lths1),
                      thetahs2 = exp(lths2))
      out
    },
    error = function(N){ 
      NULL
    }
  )
  out
}

# Simulations
RNGkind("L'Ecuyer-CMRG")
set.seed(1234)

phi0   <- c(0, 0.4) #AR parameters
beta0  <- c(0.05, 0.05, 0.9) #GARCH parameters
nu     <- 6 #Student df
theta0 <- 4 #Copula parameter
nvec   <- c(250, 500, 1e3, 2e3)
sim    <- 1e4
k      <- 1e3
S      <- c(2, 3, 5, 8)

out1   <- mclapply(1:sim, function(...) fsim(nvec[1], S[1], k), mc.cores = 4L)
out2   <- mclapply(1:sim, function(...) fsim(nvec[2], S[2], k), mc.cores = 4L)
out3   <- mclapply(1:sim, function(...) fsim(nvec[3], S[3], k), mc.cores = 4L)
out4   <- mclapply(1:sim, function(...) fsim(nvec[4], S[4], k), mc.cores = 4L)
out    <- list(out1, out2, out3, out4)
saveRDS(out, file = "Simulations/simu.garch.RDS")

# Plot CDFs
# The values at which the CDF are evaluated
t      <- seq(-40, 40, 0.001)
lt     <- length(t)

out      <- readRDS("Simulations/simu.garch.RDS")
fdata    <- function(s){
  Ns    <- nvec[s]
  Ss    <- S[s]
  datas <- out[[s]]
  datas <- datas[!sapply(datas, is.null)]
  sim   <- length(datas)
  Nsl   <- ifelse(Ns < 1000, Ns, paste0(substr(Ns, 1, nchar(Ns) - 3), ",", substr(Ns, nchar(Ns) - 2, nchar(Ns))))
  F0    <- fcdf(sapply(datas, function(x) x$delta), t, sim, lt)
  F0s1  <- fcdf(sapply(datas, function(x) x$deltas1), t, sim, lt)
  F0s2  <- fcdf(sapply(datas, function(x) x$deltas2), t, sim, lt)
  Fh    <- rowMeans(do.call(cbind, mclapply(datas, function(x) fcdf(x$psi, t, k, lt), mc.cores = 15L)))
  Fhs1  <- rowMeans(do.call(cbind, mclapply(datas, function(x) fcdf(x$psis1, t, k, lt), mc.cores = 15L)))
  Fhs2  <- rowMeans(do.call(cbind, mclapply(datas, function(x) fcdf(x$psis2, t, k, lt), mc.cores = 15L)))
  H     <- rowMeans(do.call(cbind, mclapply(datas, function(x) pnorm(t, 0, x$sderr), mc.cores = 11L)), na.rm = TRUE)
  
  
  data.frame(s    = s,
             parm = "theta0",
             type = rep(s, lt),
             N    = Ns,
             t    = t,
             F0   = F0,
             F0s1 = F0s1,
             F0s2 = F0s2,
             Fh   = Fh,
             Fhs1 = Fhs1,
             Fhs2 = Fhs2,
             H    = H)  %>% 
    mutate(type   = factor(type, labels = expression(paste(sqrt(n)*(log(hat(theta)[n]) - log(theta[0])), 
                                                           ", ", n, " = ", !!Nsl, ", ", k[n], " = ", !!Ss))),
           types   = factor(type, labels = expression(paste(sqrt(n)*(log(theta[paste(n, ",", kappa)]^"*") - log(theta[0])), 
                                                            ", ", n, " = ", !!Nsl, ", ", k[n], " = ", !!Ss)))) 
}

dataplot <- bind_rows(lapply(1:4, fdata))

# Distances
(outdist <- dataplot %>% group_by(s, parm, N) %>% 
    summarise(across(c("H", "Fh"), ~ sum(abs(.x - F0))*(max(t) - min(t))/length(t))) %>%
    inner_join(dataplot %>% group_by(s, parm, N) %>% 
                 summarise(Fhs1 = sum(abs(Fhs1 - F0s1))*(max(t) - min(t))/length(t)), by  = c("s", "parm", "N")) %>%
    inner_join(dataplot %>% group_by(s, parm, N) %>% 
                 summarise(Fhs2 = sum(abs(Fhs2 - F0s2))*(max(t) - min(t))/length(t)), by  = c("s", "parm", "N")))

# CDF of Delta
graph1        <- list()
for (k in 1:length(nvec)) {
  distF       <- format(round(outdist$Fh[k], 3), nsmall = 3)
  distH       <- format(round(outdist$H[k], 3), nsmall = 3)
  graph1[[k]] <- ggplot(dataplot %>% filter(t >= -12, t <= 8, s == k), aes(x = t, y = F0)) + theme_bw() + 
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

# CDF of Deltastart: mean correction
graph2        <- list()
for (k in 1:length(nvec)) {
  distF       <- format(round(outdist$Fhs1[k], 3), nsmall = 3)
  graph2[[k]] <- ggplot(dataplot %>% filter(t >= -8, t <= 8, s == k), aes(x = t, y = F0s1)) + theme_bw() + 
    geom_line(aes(col = "A", lty = "A")) + 
    geom_line(aes(y = Fhs1, col = "B", lty = "B")) + 
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

# CDF of Deltastart: mean correction
graph3        <- list()
for (k in 1:length(nvec)) {
  distF       <- format(round(outdist$Fhs2[k], 3), nsmall = 3)
  graph3[[k]] <- ggplot(dataplot %>% filter(t >= -8, t <= 8, s == k), aes(x = t, y = F0s2)) + theme_bw() + 
    geom_line(aes(col = "A", lty = "A")) + 
    geom_line(aes(y = Fhs2, col = "B", lty = "B")) + 
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
    graph3[[k]] <- graph3[[k]] +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
  }
}

# export figures
ggsave("simu.garch.pdf", path = "Simulations", plot = wrap_plots(graph1, ncol = length(nvec)), 
       device = "pdf", width = 10, height = 2.5)
ggsave("simu.garchR.pdf", path = "Simulations", plot = wrap_plots(graph2, ncol = length(nvec)), 
       device = "pdf", width = 10, height = 2.5)
ggsave("simu.garch.median.pdf", path = "Simulations", plot = wrap_plots(graph3, ncol = length(nvec)), 
       device = "pdf", width = 10, height = 2.5)

# Coverage rates
fcrate  <- function(s){
  Ns    <- nvec[s]
  datas <- out[[s]]
  datas <- datas[!sapply(datas, is.null)]
  datas <- datas[!sapply(datas, function(x) is.na(x$thetahs1))]
  lth   <- log(theta0)
  IC1   <- sapply(datas, function(x) log(x$thetah) + qnorm(c(0.025, 0.05, 0.95, 0.975), 0, x$sderr/sqrt(Ns)))
  IC    <- sapply(datas, function(x) quantile(log(x$thetah) - x$psi/sqrt(Ns), probs = c(0.025, 0.05, 0.95, 0.975), na.rm = TRUE))
  ICs1  <- sapply(datas, function(x) quantile(log(x$thetahs1) - x$psis1/sqrt(Ns), probs = c(0.025, 0.05, 0.95, 0.975), na.rm = TRUE))
  ICs2  <- sapply(datas, function(x) quantile(log(x$thetahs2) - x$psis2/sqrt(Ns), probs = c(0.025, 0.05, 0.95, 0.975), na.rm = TRUE))
  
  crate <- c(mean(lth >= IC1[1,] & lth <= IC1[4,], na.rm = TRUE), mean(lth <= IC1[3,], na.rm = TRUE), mean(lth >= IC1[2,], na.rm = TRUE),
             mean(lth >= IC[1,] & lth <= IC[4,], na.rm = TRUE), mean(lth <= IC[3,], na.rm = TRUE), mean(lth >= IC[2,], na.rm = TRUE),
             mean(lth >= ICs1[1,] & lth <= ICs1[4,], na.rm = TRUE), mean(lth <= ICs1[3,], na.rm = TRUE), mean(lth >= ICs1[2,], na.rm = TRUE),
             mean(lth >= ICs2[1,] & lth <= ICs2[4,], na.rm = TRUE), mean(lth <= ICs2[3,], na.rm = TRUE), mean(lth >= ICs2[2,], na.rm = TRUE))
  crate <- as.data.frame(crate)
  colnames(crate) <- paste0("N=", Ns)
  rownames(crate) <- paste(rep(c("Naive", "Simu", "Simu-Debiased-Mean", "Simu-Debiased-Median"), each = 3), 
                           rep(c("2-Sides", "Lower 1-Side", "Upper 1 Side"), 4))
  as.data.frame(crate)
}

(crate  <- do.call(cbind, lapply(1:4, fcrate)))
write.csv(crate, file = "Simulations/simu.garch.crate.csv")

# Bias
fbias <- function(s){
  Ns    <- nvec[s]
  datas <- out[[s]]
  datas <- datas[!sapply(datas, is.null)]
  datas <- datas[!sapply(datas, function(x) is.na(x$thetahs1))]
  theh  <- sapply(datas, function(x) x$thetah)
  thes1 <- sapply(datas, function(x) x$thetahs1)
  thes2 <- sapply(datas, function(x) x$thetahs2)
  data.frame(parm  = "theta0",
             N     = Ns,
             the0  = theta0,
             theh  = theh,
             thes1 = thes1,
             thes2 = thes2)
}

(databias <- bind_rows(lapply(1:4, fbias)) %>% group_by(parm, N) %>%
    summarise(Mean1 = mean(theh), BIAS1 = mean(theh - the0), SD1 = sd(theh), 
              RMSE1 = sqrt(mean((theh - the0)^2)), MAE1 = mean(abs(theh - the0)),
              Mean2 = mean(thes1), BIAS2 = mean(thes1 - the0), SD2 = sd(thes1), 
              RMSE2 = sqrt(mean((thes1 - the0)^2)), MAE2 = mean(abs(thes1 - the0)),
              Mean3 = mean(thes2), BIAS3 = mean(thes2 - the0), SD3 = sd(thes2), 
              RMSE3 = sqrt(mean((thes2 - the0)^2)), MAE3 = mean(abs(thes2 - the0))))
write.csv(databias, file = "Simulations/simu.garch.bias.csv", row.names = FALSE)