# This script replicates simulation results for Data Generating Process (DGP) C.
# Estimations are performed in parallel to enhance efficiency.

## Libraries
rm(list = ls())
library(copula)
library(MASS)
library(numDeriv)
library(Matrix)
library(dplyr)
library(doParallel)
library(ggplot2)

## Loading R and C++ functions
source("tsesimu.rfunctions.R")
Rcpp::sourceCpp("tsesimu.cppfunctions.cpp")

## True values of the parameters
phi0   <- c(0, 0.4) #AR parameters
beta0  <- c(0.05, 0.05, 0.9) #GARCH parameters
nu     <- 6 #Student df
theta0 <- 4 #Copula parameter

## The values at which the CDF are evaluated
t      <- seq(-12, 8, 0.05)
lt     <- length(t)

## Sample sizes
nvec <- c(250, 500, 1e3, 2e3)

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
      zh      <- sapply(1:S, function(s){
        filterz(N, coeh[3:5, s], coeh[1:2, s], dt[,s])})
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
      
      # Naive standard errors
      An      <- -sum(d2lclaytonpdf(u = pzh, ltheta = lth, dim = S))/N
      snai1   <- c(sqrt(fHAC(d1lclaytonpdf(u = pzh, ltheta = lth, dim = S), N, r = 1))/An)
      
      # psi
      coet    <- mvrnorm(n = k, mu = c(apply(coeh, 2, fbetat)), Sigma = MCOt)
      coes    <- lapply(1:S, function(s) t(apply(coet[,((s - 1)*6 + 1):(s*6)], 1, fbeta)))
      nus     <- sapply(coes, function(s) s[,6])
      zs      <- lapply(1:k, function(x){sapply(1:S, function(s){
        filterz(N, coes[[s]][x, 3:5], coes[[s]][x, 1:2], dt[,s])})})
      pzs     <- lapply(1:k, function(x) sapply(1:S, function(s) pt(zs[[x]][,s]/sqrt(1 - 2/nus[x, s]), df = nus[x, s])))
      En      <- sapply(1:k, function(x) sum(d1lclaytonpdf(u = pzs[[x]], ltheta = lth, dim = S)))/sqrt(N)
      Vn      <- sapply(1:k, function(x) fHAC(d1lclaytonpdf(u = pzs[[x]], ltheta = lth, dim = S), N, r = 6*S))
      psi     <- sapply(1:k, function(x) (sqrt(Vn[x])*rnorm(1) + En[x])/An)
      
      # Standard errors
      snai2   <- sd(psi)
      
      # Quantiles using each inference method
      quant1  <- lth + snai1/sqrt(N)*qnorm(c(0.005, 0.025, 0.05, 0.950, 0.975,  0.995))
      quant2  <- lth + snai2/sqrt(N)*qnorm(c(0.005, 0.025, 0.05, 0.950, 0.975,  0.995))
      quant   <- lth - quantile(psi/sqrt(N), probs = c(0.005, 0.025, 0.05, 0.950, 0.975,  0.995))
      names(quant1)  <- c("0.5%", "2.5%", "5%", "95%", "97.5%",  "99.5%")
      names(quant2)  <- c("0.5%", "2.5%", "5%", "95%", "97.5%",  "99.5%")
      names(quant)   <- c("0.5%", "2.5%", "5%", "95%", "97.5%",  "99.5%")
      
      out     <- list(delta  = sqrt(N)*(lth - log(theta0)),
                      psi    = fcdf(psi, t, k, lt),
                      snai1  = snai1,
                      snai2  = snai2,
                      quant1 = quant1,
                      quant2 = quant2,
                      quant  = quant)
      out
    },
    error = function(N){ 
      NULL
    }
  )
  out
}

# Simulations
set.seed(1234)
sim  <- 1e4
k    <- 1e3
S    <- c(2, 3, 5, 8)
out1 <- mclapply(1:sim, function(...) fsim(nvec[1], S[1], k), mc.cores = 4L)
out2 <- mclapply(1:sim, function(...) fsim(nvec[2], S[2], k), mc.cores = 4L)
out3 <- mclapply(1:sim, function(...) fsim(nvec[3], S[3], k), mc.cores = 4L)
out4 <- mclapply(1:sim, function(...) fsim(nvec[4], S[4], k), mc.cores = 4L)
out  <- list(out1, out2, out3, out4)

# Tableau of proportions of CI contenting the true value
coverage <- lapply(c("N = 250" = 1, "N = 500" = 2, "N = 1000" = 3, "N = 1000" = 4), function(s){
  f.coverage.sca(out[[s]], log(theta0))})
print(coverage)

# Plot CDFs
fdata    <- function(s){
  tryCatch(
    expr = {
      Ns     <- nvec[s]
      datas  <- out[[s]]
      datas  <- datas[!sapply(datas, is.null)]
      sim    <- length(datas)
      Nsl    <- ifelse(Ns < 1000, Ns, paste0(substr(Ns, 1, nchar(Ns) - 3), ",", substr(Ns, nchar(Ns) - 2, nchar(Ns))))
      F0     <- fcdf(sapply(datas, function(x) x$delta), t, sim, lt)
      Fh     <- rowMeans(sapply(datas, function(x) x$psi))
      Norn1  <- pnorm(t, 0, mean(sapply(datas, function(x) x$snai1), na.rm = TRUE))
      Norn2  <- pnorm(t, 0, median(sapply(datas, function(x) x$snai2), na.rm = TRUE))
      
      data.frame(t    = t,
                 F0   = F0,
                 Fh   = Fh,
                 Nn1  = Norn1,
                 Nn2  = Norn2,
                 N    = rep(s, lt),
                 type = rep(s, lt))  %>% 
        mutate(type   = factor(type, labels = expression(paste(sqrt(n)*(log(hat(theta)[n]) - log(theta[0])), " with n = ", !!Nsl)))) 
    },
    error = function(N){ 
      NULL
    })
}

dataplot <- bind_rows(lapply(1:4, fdata))

(graph   <- ggplot(dataplot %>% filter(N %in% 1:4), aes(x = t, y = F0)) + theme_bw() + 
    geom_line(aes(col = "A", lty = "A")) + 
    geom_line(aes(y = Fh, col = "B", lty = "B")) + 
    geom_line(aes(y = Nn1, col = "C", lty = "C")) + 
    geom_line(aes(y = Nn2, col = "D", lty = "D")) + 
    facet_wrap(. ~ type, scales = "free", labeller = label_parsed) +
    # xlab(expression(paste(sqrt(n)*(log(hat(theta)[n]) - log(theta[0]))))) +
    xlab("") + 
    ylab("Probability") + 
    scale_colour_manual("", values = c("black", "#e01b89", "#1b1be0", "#1b1be0"), 
                        labels = c("A" = "True distribution", "B" = "Estimation", "C" = 'Normal 1', "D" = 'Normal 2')) +
    scale_linetype_manual("", values = c("A" = 1, "B" = 2, "C" = 3, "D" = 4),
                          labels = c("A" = "True distribution", "B" = "Estimation", "C" = 'Normal 1', "D" = 'Normal 2')) +
    # theme(legend.position = c(0.06, 0.75), legend.text.align = 0,
    #       legend.background = element_rect(fill='transparent'),
    #       legend.text = element_text(size = 8)) +
    theme(legend.position = "none") + 
    facet_wrap(. ~ type, scales = "free", labeller = label_parsed, dir="v", ncol = 4))

# export figures
ggsave("simu:garch.pdf", plot = graph, device = "pdf", width = 10, height = 2.5)
