# This script replicates simulation results for Data Generating Process (DGP) C.
# Estimations are performed in parallel to enhance efficiency.

## Add your working directory
proot <- c("~/Dropbox/2-steps M estimator", 
           "~/Dropbox/2StageInference")
root  <- sapply(proot, dir.exists)
root  <- proot[root][1]
setwd(root)

## Libraries
rm(list = ls())
library(MASS)
library(splines)
library(dplyr)
library(doParallel)
library(ggplot2)
library(patchwork)

## Loading C++ functions
Rcpp::sourceCpp("Simulations/Codes/tsesimu.cppfunctions.cpp")

## This function performs a single iteration of the Monte Carlo simulation.
fsim      <- function(N, k){# k is kappa in the paper
  tryCatch(
    expr = {
      ## Data
      nz      <- c(1, 0.985, 0.945, 0.91)[nvec == N]
      n       <- round(N^nz)
      Z       <- runif(N, 0, 10)
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
      sdnai   <- step2$coefficients[,"Std. Error"]*sqrt(N)
      
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
      Vh      <- crossprod(mhat, mhat*lhat)/N
      psi     <- sapply(1:k, function(x) solve(An, t(chol(Vh)) %*% rnorm(2) + En[,x]))
      
      # Standard errors
      sderr   <- solve(An) %*% (Vh + var(t(En))*k/(k - 1)) %*% solve(An)
      sderr   <- sqrt(diag(sderr))
      
      # Bias reduction
      thetahs  <- thetah - solve(An, apply(En, 1, mean))/sqrt(N)
      lhat    <- c(exp(mhat %*% thetahs)) #lambdahat
      lhats   <- sapply(mhats, function(x) c(exp(x %*% thetahs))) #lambdahat_s
      An      <- -crossprod(mhat, mhat*lhat)/N
      En      <- sapply(1:k, function(x) crossprod(mhats[[x]], lhats[,x] - lhat)/sqrt(N))
      Omega   <- apply(En, 1, mean)
      Vh      <- crossprod(mhat, mhat*lhat)/N
      psis    <- sapply(1:k, function(x) solve(An, t(chol(Vh)) %*% rnorm(2) + En[,x] - Omega))
      
      list(delta   = sqrt(N)*(thetah - theta0),
                      deltas  = sqrt(N)*(thetahs - theta0),
                      psi     = psi,
                      psis    = psis,
                      sdnai   = sdnai,
                      sderr   = sderr,
                      thetah  = thetah,
                      thetahs = thetahs)
    },
    error = function(N){ 
      NULL
    }
  )
}

# Simulations
RNGkind("L'Ecuyer-CMRG")
set.seed(1234)

theta0 <- c(-0.8, 2)
nvec   <- c(250, 500, 1e3, 2e3)
sim    <- 1e4
k      <- 1e3

out1   <- mclapply(1:sim, function(...) fsim(nvec[1], k), mc.cores = 15L)
out2   <- mclapply(1:sim, function(...) fsim(nvec[2], k), mc.cores = 15L)
out3   <- mclapply(1:sim, function(...) fsim(nvec[3], k), mc.cores = 15L)
out4   <- mclapply(1:sim, function(...) fsim(nvec[4], k), mc.cores = 15L)
out    <- list(out1, out2, out3, out4)
saveRDS(out, file = "Simulations/simu:poisson.RDS")

# Plot CDFs
# The values at which the CDF are evaluated
t        <- seq(-30, 24, 0.001)
lt       <- length(t)

out      <- readRDS("Simulations/simu:poisson.RDS")

fdata    <- function(s){
  Ns     <- nvec[s]
  nz     <- round(Ns^(c(1, 0.985, 0.945, 0.91)[nvec == Ns]))
  datas  <- out[[s]]
  datas  <- datas[!sapply(datas, is.null)]
  sim   <- length(datas)
  Nsl    <- ifelse(Ns < 1000, Ns, paste0(Ns/1000, ",000"))
  nzl    <- ifelse(nz < 1000, nz, paste0(floor(nz/1000), ",", substr(nz/1000, 3, 5)))
  F0     <- apply(sapply(datas, function(x) x$delta), 1, function(x) fcdf(x, t, sim, lt))
  F0s    <- apply(sapply(datas, function(x) x$deltas), 1, function(x) fcdf(x, t, sim, lt))
  Fh     <- sapply(1:2, function(x1) rowMeans(do.call(cbind, mclapply(datas, function(x2) fcdf(x2$psi[x1,], t, k, lt), mc.cores = 15L))))
  Fhs    <- sapply(1:2, function(x1) rowMeans(do.call(cbind, mclapply(datas, function(x2) fcdf(x2$psis[x1,], t, k, lt), mc.cores = 15L))))
  H1     <- sapply(1:2, function(x1) rowMeans(do.call(cbind, mclapply(datas, function(x2) pnorm(t, 0, x2$sdnai[x1]), mc.cores = 11L))))
  H2     <- sapply(1:2, function(x1) rowMeans(do.call(cbind, mclapply(datas, function(x2) pnorm(t, 0, x2$sderr[x1]), mc.cores = 11L))))
  
  
  data.frame(s    = rep((2*s -1):(2*s), each = lt),
             parm = c(rep("theta01", lt), rep("theta02", lt)),
             type = rep((2*s -1):(2*s), each = lt),
             N    = Ns,
             nz   = nz,
             t    = rep(t, 2),
             F0   = c(F0[,1], F0[,2]),
             F0s  = c(F0s[,1], F0s[,2]),
             Fh   = c(Fh[,1], Fh[,2]),
             Fhs  = c(Fhs[,1], Fhs[,2]),
             H1   = c(H1[,1], H1[,2]),
             H2   = c(H2[,1], H2[,2]))  %>% 
    mutate(type   = factor(type, labels = c(expression(paste(sqrt(n)*(hat(theta)["n,1"] - theta["0,1"]), 
                                                             ", ", n, " = ", !!Nsl, ", ", n, "*", " = ", !!nzl)), 
                                            expression(paste(sqrt(n)*(hat(theta)["n,2"] - theta["0,2"]), 
                                                             ", ", n, " = ", !!Nsl, ", ", n, "*", " = ", !!nzl)))),
           typec  = factor(type, labels = c(expression(paste(sqrt(n)*(hat(theta)["n,1"]^{"*"} - theta["0,1"]), 
                                                             ", ", n, " = ", !!Nsl, ", ", n, "*", " = ", !!nzl)), 
                                            expression(paste(sqrt(n)*(hat(theta)["n,2"]^{"*"} - theta["0,2"]), 
                                                             ", ", n, " = ", !!Nsl, ", ", n, "*", " = ", !!nzl))))) 
}

dataplot <- bind_rows(lapply(1:4, fdata))

# Distances
(outdist <- dataplot %>% group_by(s, parm, N, nz) %>% 
    summarise(across(c("H1", "H2", "Fh"), ~ sum(abs(.x - F0))*(max(t) - min(t))/length(t))) %>% 
    inner_join(dataplot %>% group_by(s, parm, N, nz) %>% 
                 summarise(Fhs = sum(abs(Fhs - F0s))*(max(t) - min(t))/length(t)), by  = c("s", "parm", "N", "nz")))

# CDF of Delta
graph1        <- list()
for (k in 1:(2*length(nvec))) {
  distF       <- format(round(outdist$Fh[k], 3), nsmall = 3)
  distH1      <- format(round(outdist$H1[k], 3), nsmall = 3)
  distH2      <- format(round(outdist$H2[k], 3), nsmall = 3)
  graph1[[k]] <- ggplot(dataplot %>% filter(t >= -15, t <= 12, s == k), aes(x = t, y = F0)) + theme_bw() + 
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
  graph2[[k]] <- ggplot(dataplot %>% filter(t >= -12, t <= 12, s == k), aes(x = t, y = F0s)) + theme_bw() + 
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
ggsave("simu:poisson.pdf", path = "Simulations", plot = wrap_plots(graph1, ncol = length(nvec), byrow = FALSE), 
       device = "pdf", width = 10, height = 5)
ggsave("simu:poissonR.pdf", path = "Simulations", plot = wrap_plots(graph2, ncol = length(nvec), byrow = FALSE), 
       device = "pdf", width = 10, height = 5)

# Bias
fbias <- function(s){
  Ns    <- nvec[s]
  nz    <- round(Ns^(c(1, 0.985, 0.945, 0.91)[nvec == Ns]))
  datas <- out[[s]]
  datas <- datas[!sapply(datas, is.null)]
  sim   <- length(datas)
  theh  <- c(t(sapply(datas, function(x) x$thetah)))
  thes  <- c(t(sapply(datas, function(x) x$thetahs)))
  data.frame(parm = rep(c("theta01", "theta02"), each = sim),
             N    = Ns,
             nz   = nz,
             the0 = rep(theta0, each = sim),
             theh = theh,
             thes = thes)
}

(databias <- bind_rows(lapply(1:4, fbias)) %>% group_by(parm, N, nz) %>%
    summarise(Mean1 = mean(theh), BIAS1 = mean(theh - the0), SD1 = sd(theh), 
              RMSE1 = sqrt(mean((theh - the0)^2)), MAE1 = mean(abs(theh - the0)),
              Mean2 = mean(thes), BIAS2 = mean(thes - the0), SD2 = sd(thes), 
              RMSE2 = sqrt(mean((thes - the0)^2)), MAE2 = mean(abs(thes - the0))))
write.csv(databias, file = "Simulations/simu:poisson.bias.csv", row.names = FALSE)