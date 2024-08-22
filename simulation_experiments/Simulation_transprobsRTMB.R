## Simulation experiments

## packages
# install.packages("RTMB")
library(RTMB)
# devtools::install_github("janoleko/LaMa") # development version
library(LaMa)
# install.packages("scales")
library(scales) # for transparent colors


# Simulation study functions ----------------------------------------------

# functions for transition probabilities (on logit scale)
sf1 = function(x) -2 + sin(3*pi*x) + exp(1.5*x)
sf2 = function(x) 2 + cos(4*pi*x) - 2*exp(x)

# sim = simHMM(z[1:2000])
# ord = order(z[1:2000])
# Gamma = sim$Gamma
# 
# par(mfrow = c(1,2))
# plot(z[1:2000][ord], Gamma[1,2,ord], type = "l", ylim = c(0,1))
# plot(z[1:2000][ord], Gamma[2,1,ord], type = "l", ylim = c(0,1))


# function to simulate data from the inhomogeneous HMM
simHMM = function(Gamma, delta = rep(0.5,2), mu = c(1,5), sigma = c(1,3)){
  n = dim(Gamma)[3]
  # simulating Markov chain
  s = rep(NA, n)
  s[1] = sample(1:2, 1, prob = delta)
  for(t in 2:n){
    s[t] = sample(1:2, 1, prob = Gamma_sim[s[t-1],,t])
  }
  # simulating observations
  x = rnorm(n, mean = mu[s], sd = sigma[s])
  
  return(list(x = x, s = s))
}

# penalized likelihood function for non-parametric fit
pnll_sim = function(par){
  getAll(par, dat)
  sigma = exp(logsigma) # transformation for unconstraint optimization
  delta = c(1, exp(logitdelta))
  delta = delta / sum(delta)  # initial distribution
  Gamma = tpm_g(Z, cbind(beta0, betaSpline), ad = T) # building tpm with LaMa
  # state-dependent probabilties
  allprobs = cbind(dnorm(x, mu[1], sigma[1]), dnorm(x, mu[2], sigma[2]))
  # forward algorithm + penalty
  -forward_g(delta, Gamma, allprobs, ad = T) + penalty(betaSpline, S, lambda)
}


# Simulation --------------------------------------------------------------

# number of simulation runs
simiter = 250

# drawing covariate values
set.seed(1)
z = runif(1e4)

# computing tpm based on covariate
Gamma_sim = array(dim = c(2,2,1e4))
for(t in 1:1e4){
  G = diag(2)
  G[!G] = c(exp(sf2(z[t])), exp(sf1(z[t])))
  Gamma_sim[,,t] = G / rowSums(G)
}

# creating model matrices
nb = 15 # number of basis functions
modmat = make_matrices(~s(z, bs = "ps", k = nb), 
                       data = data.frame(z = z))
Zfull = modmat$Z
S = modmat$S

# initial parameter list
par = list(mu = c(1,5), 
           logsigma = log(c(1,3)), 
           logitdelta = 0, 
           beta0 = c(-2,2),
           betaSpline = matrix(0, 2, ncol = nb-1))


## three simulation scenarios

## n = 1000
n = 1000
Z = Zfull[1:n,]
dat = list(Z = Z, S = S, lambda = rep(1000, 2))

# simulation
modsSim1000 = list()
for(l in 1:simiter) {
  cat("Model", l, "\n")
  
  sim = simHMM(Gamma_sim[,,1:n])
  dat$x = sim$x
  
  thismod = tryCatch(pql(pnll_sim, par, dat, random = "betaSpline", 
                         silent = 2, maxiter = 100),
                     error = function(e) "An error occurred")
  
  if(!is.list(thismod)){
    cat("Error\n")
    next
  } else{
    modsSim1000[[l]] = thismod
  }
}
# modsSim1000 = modsSim1000[!sapply(modsSim1000, is.null)]
# saveRDS(modsSim1000[1:200], "./simulation_experiments/mods/n1000.rds")


# n = 2000
n = 2000
Z = Zfull[1:n,]
dat = list(Z = Z, S = S, lambda = rep(1000, 2))

# simulation
modsSim2000 = list()
for(l in 1:simiter) {
  cat("Model", l, "\n")
  
  sim = simHMM(Gamma_sim[,,1:n])
  dat$x = sim$x
  
  thismod = tryCatch(pql(pnll_sim, par, dat, random = "betaSpline", 
                         silent = 2, maxiter = 100),
                     error = function(e) "An error occurred")
  
  if(!is.list(thismod)){
    cat("Error\n")
    next
  } else{
    modsSim2000[[l]] = thismod
  }
}
# modsSim2000 = modsSim2000[!sapply(modsSim2000, is.null)]
# saveRDS(modsSim2000[1:200], "./simulation_experiments/mods/n2000.rds")


# n = 5000
n = 5000
Z = Zfull[1:n,]
dat = list(Z = Z, S = S, lambda = rep(1000, 2))

# simulation
modsSim5000 = list()
for(l in 1:simiter) {
  cat("Model", l, "\n")
  
  sim = simHMM(Gamma_sim[,,1:n])
  dat$x = sim$x
  
  thismod = tryCatch(pql(pnll_sim, par, dat, random = "betaSpline", 
                         silent = 1, maxiter = 100),
                     error = function(e) "An error occurred")
  
  if(!is.list(thismod)){
    cat("Error\n")
    next
  } else{
    modsSim5000[[l]] = thismod
  }
}
# modsSim5000 = modsSim5000[!sapply(modsSim5000, is.null)]
# saveRDS(modsSim5000[1:200], "./simulation_experiments/mods/n5000.rds")




# Visualizing results -----------------------------------------------------

zseq = seq(0, 1, length = 200)
Z_plot = pred_matrix(modmat, data.frame(z = zseq))

Gamma_plot = array(dim = c(2,2,200))
for(t in 1:200){
  G = diag(2)
  G[!G] = c(exp(sf2(zseq[t])), exp(sf1(zseq[t])))
  Gamma_plot[,,t] = G / rowSums(G)
}

# modsSim1000 = readRDS("./simulation_experiments/mods/n1000.rds")
# modsSim2000 = readRDS("./simulation_experiments/mods/n2000.rds")
# modsSim5000 = readRDS("./simulation_experiments/mods/n5000.rds")

allmods = list(modsSim1000, modsSim2000, modsSim5000)
Ts = c(1000, 2000, 5000)


## Results for gamma_12

# pdf("./simulation_experiments/figs/transprobs_simulation_new.pdf", width = 7, height = 4.5)

par(mfrow = c(2,3), mar = c(5,4.4,2,1)+0.1)
i = 1; j = 2
for(m in 1:3){
  modsSim = allmods[[m]]

  plot(zseq, Gamma_plot[i,j,], type = "l", lwd = 1, col = "black", main = paste0("T = ", Ts[m]),
       ylab = expression(gamma[12]^(t)), xlab = expression(z[t]), bty = "n", ylim = c(0, 1))
  for(k in 1:length(modsSim)){
    # beta_k = matrix(modsSim[[k]]$mod$estimate[5 + 1:(2*(nb+1))], ncol = 2)
    beta_k = modsSim[[k]]$beta
    Gamma_k = tpm_g(Z_plot, beta_k)

    lines(zseq, Gamma_k[i,j,], col = alpha("orange", 0.1), lwd = 1)
  }
  lines(zseq, Gamma_plot[i,j,], lwd = 1.5, col = "black")

}
legend(-0.05, 1, legend = c("true function", "estimates"), col = c("black", "orange"), lwd = c(1, 1), bty = "n")

# dev.off()


# Trace plots

for(m in 1:3){
  modsSim = allmods[[m]]

  plot(NA, bty = "n", xlim = c(0,30), ylim = c(0, 1000),
       ylab = "penalty strength", xlab = "iteration")

  for(k in 1:length(modsSim)){
    lambdas = matrix(unlist(modsSim[[k]]$all_lambda), ncol = 2, byrow = TRUE)[,2]
    lines(lambdas, col = alpha("orange", 0.1))
  }

  lastlambdas = rep(NA, length(modsSim))
  for(k in 1:length(modsSim)){
    lastlambdas[k] = modsSim[[k]]$lambda[2]
  }

  iters = sapply(1:k, function(k) length(modsSim[[k]]$all_lambda))

  abline(v = median(iters), col = "black", lwd = 1, lty = 3)

  abline(h = median(lastlambdas), col = "black", lwd = 1, lty = 3)

}
# legend("topright", legend = c("mean final penalty strength"), col = c("black"), lwd = c(2), lty = c(2), bty = "n")

# dev.off()

# we see that convergence is much faster and much better for more data, as expected
# more stable, because Laplace approximation gets better with more data


## Results for gamma_21

# pdf("./simulation_experiments/figs/transprobs_simulation2_new.pdf", width = 7, height = 4.5)

par(mfrow = c(2,3), mar = c(5,4.4,2,1)+0.1)
i = 2; j = 1
for(m in 1:3){
  modsSim = allmods[[m]]
  
  plot(zseq, Gamma_plot[i,j,], type = "l", lwd = 1, col = "black", main = paste0("T = ", Ts[m]),
       ylab = expression(gamma[21]^(t)), xlab = expression(z[t]), bty = "n", ylim = c(0, 1))
  for(k in 1:length(modsSim)){
    # beta_k = matrix(modsSim[[k]]$mod$estimate[5 + 1:(2*(nb+1))], ncol = 2)
    beta_k = modsSim[[k]]$beta
    Gamma_k = tpm_g(Z_plot, beta_k)
    
    lines(zseq, Gamma_k[i,j,], col = alpha("deepskyblue", 0.1), lwd = 1)
  }
  lines(zseq, Gamma_plot[i,j,], lwd = 1.5, col = "black")
  
}
legend("topright", legend = c("true function", "estimates"), col = c("black", "deepskyblue"), lwd = c(1, 1), bty = "n")


# Trace plots

for(m in 1:3){
  modsSim = allmods[[m]]
  
  plot(NA, bty = "n", xlim = c(0,30), ylim = c(0, 1000),
       ylab = "penalty strength", xlab = "iteration")
  
  for(k in 1:length(modsSim)){
    lambdas = matrix(unlist(modsSim[[k]]$all_lambda), ncol = 2, byrow = TRUE)[,1]
    lines(lambdas, col = alpha("deepskyblue", 0.1))
  }
  
  lastlambdas = rep(NA, length(modsSim))
  for(k in 1:length(modsSim)){
    lastlambdas[k] = modsSim[[k]]$lambda[1]
  }
  
  iters = sapply(1:k, function(k) length(modsSim[[k]]$all_lambda))
  
  abline(v = median(iters), col = "black", lwd = 1, lty = 3)
  
  abline(h = median(lastlambdas), col = "black", lwd = 1, lty = 3)
  
}
# legend("topright", legend = c("mean final penalty strength"), col = c("black"), lwd = c(2), lty = c(2), bty = "n")

# dev.off()

