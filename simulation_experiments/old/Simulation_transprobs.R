

# Simulation study functions ----------------------------------------------

# functions for transprobs
sf1 = function(x) -2 + sin(3*pi*x) + exp(1.5*x)
sf2 = function(x) 2 + cos(4*pi*x) - 2*exp(x)

# sim = simHMM(z[1:2000])
# ord = order(z[1:2000])
# Gamma = sim$Gamma
# 
# par(mfrow = c(1,2))
# plot(z[1:2000][ord], Gamma[1,2,ord], type = "l", ylim = c(0,1))
# plot(z[1:2000][ord], Gamma[2,1,ord], type = "l", ylim = c(0,1))

# function to simulate data
simHMM = function(z, delta = rep(0.5,2), mu = c(1,5), sigma = c(1,3)){
  # defining functions for transprobs
  # sf1 = function(x) 0.05 * (x + sin(x)*x^2 + 5*exp(0.3*x))
  # sf2 = function(x) 0.05 * (x + 0.4*sin(-x)*x^2 + 10*exp(-0.3*x))
  n = length(z) # length of simulated time series
  # building tpm array with specified functions
  Gamma_sim = array(dim = c(2,2,n))
  for(t in 1:n){
    G = diag(2)
    G[!G] = c(exp(sf2(z[t])),
              exp(sf1(z[t])))
    Gamma_sim[,,t] = G / rowSums(G)
  }
  # simulating Markov chain
  s = rep(NA, n)
  s[1] = sample(1:2, 1, prob = delta)
  for(t in 2:n){
    s[t] = sample(1:2, 1, prob = Gamma_sim[s[t-1],,t])
  }
  # simulating observations
  x = rnorm(n, mean = mu[s], sd = sigma[s])
  
  return(list(x = x, s = s, Gamma = Gamma_sim))
}

# likelihood function for non-parametric fit
mllk_sim = function(theta.star, x, Z, S, lambda){
  # state-dependent parameters
  mu = theta.star[1:2]
  sigma = exp(theta.star[3:4])
  # initial distribution
  delta = c(1, exp(theta.star[5])); delta = delta / sum(delta)
  # spline coefficients
  nb = ncol(Z) - 1
  beta = matrix(theta.star[5 + 1:(2*(nb+1))], ncol = 2)
  # building tpm with LaMa
  Gamma = LaMa::tpm_g(Z[,-1], t(beta))
  # state-dependent probabilties
  allprobs = matrix(1, nrow = length(x), ncol = 2)
  for(j in 1:2) allprobs[,j] = dnorm(x, mu[j], sigma[j])
  # calculating the penalty
  pen = t(as.numeric(beta[-1,])) %*% kronecker(diag(lambda), S) %*% as.numeric(beta[-1,])
  # forward algorithm with LaMa
  -LaMa::forward_g(delta, Gamma, allprobs) + 0.5 * pen
}

# function to fit the model via marginal ML
fitMML = function(x, Z, S, theta0, mllk, REind = matrix(5 + 1:(2*ncol(Z)), nrow = 2, byrow = TRUE)[,-1],
                  lambda0 = c(1000, 1000), maxiter = 30, tol = 0.05, gradtol = 1e-6, print.level = 0) {
  Lambdas = matrix(NA, maxiter, 2)
  Lambdas[1,] = lambda0
  mods = list()
  # initial parameter values
  theta.star = theta0
  for(k in 1:maxiter){
    t1 = Sys.time()
    mods[[k]] = nlm(mllk, theta.star, x = x, Z = Z, lambda = Lambdas[k,], S = S,
                    iterlim = 1000, print.level = print.level, gradtol = gradtol, hessian = TRUE)
    mods[[k]]$time = Sys.time()-t1

    theta.star = mods[[k]]$estimate # saves theta.star as starting value for next iteration
    J_inv = MASS::ginv(mods[[k]]$hessian) # inverse penalized hessian
      
    # updating all penalty strengths state-dependent process
    for(i in 1:2){
      edoF = sum(diag(diag(nrow(S)) - Lambdas[k, i] * J_inv[REind[i,], REind[i,]] %*% S))
      penalty = t(theta.star[REind[i,]]) %*% S %*% theta.star[REind[i,]]
      Lambdas[k+1, i] = as.numeric((edoF-1) / penalty)
    }
    Lambdas[k+1, which(Lambdas[k+1,] < 0)] = 0
      
    if(mean(abs(Lambdas[k+1,] - Lambdas[k,]) / Lambdas[k,]) < tol){
      t1 = Sys.time()
      mods[[k+1]] = nlm(mllk, theta.star, x = x, Z = Z, lambda = Lambdas[k+1,], S = S,
                        iterlim = 1000, print.level = print.level, hessian = TRUE)
      mods[[k+1]]$time = Sys.time()-t1
      break
    }
  }
  Lambdas = as.matrix(na.omit(Lambdas))
    
  return(list(
    mod = mods[[length(mods)]], 
    lambda = Lambdas[nrow(Lambdas),], 
    estTime = time,
    allLambdas = Lambdas)
  )
}

# one repetition
oneRep = function(z, Z, S, theta0, mllk, simHMM, fitMML, REind = matrix(5 + 1:(2*ncol(Z)), nrow = 2, byrow = TRUE)[,-1],
                  lambda0 = c(5000, 5000), maxiter = 30, tol = 0.05, gradtol = 1e-6, print.level = 0){
  sim = simHMM(z)
  fitMML(sim$x, Z, S, theta0, mllk, REind, lambda0, maxiter, tol, gradtol, print.level)
}



# Simulation --------------------------------------------------------------

nb = 15 # number of basis functions

mu0 = c(1,5)
sigma0 = c(1,3)
beta0 = c(2,-2)
theta0 = c(mu0, log(sigma0), 0, beta0[1], rep(0, nb), beta0[2], rep(0, nb))

set.seed(1)
z = runif(1e4)


# n = 1000
n = 1000

gam_pre = mgcv::gam(y ~ s(z, bs = "ps", k = nb+1), 
                    data = data.frame(y = 1, dummy = 1, z = z[1:n]),
                    knots = list(seq(0, 1, length = nb)),
                    fit = FALSE)
Z = gam_pre$X
S = gam_pre$S[[1]]

# simulation

simiter = 250
modsSim = list()
for(l in 1:simiter) {
  cat("\nModel", l)
  modsSim[[l]] = tryCatch(
    oneRep(z[1:n], Z, S, theta0, mllk_sim, simHMM, fitMML, tol = 0.01),
    error = function(e) "An error occurred")
}
saveRDS(modsSim, "./simulation_experiments/mods/mods_1000_2.rds")


# n = 2000
n = 2000

gam_pre = mgcv::gam(y ~ s(z, bs = "ps", k = nb+1), 
                    data = data.frame(y = 1, dummy = 1, z = z[1:n]),
                    knots = list(seq(0, 1, length = nb)),
                    fit = FALSE)
Z = gam_pre$X
S = gam_pre$S[[1]]

# simulation

simiter = 230
modsSim = list()
for(l in 1:simiter) {
  cat("\nModel", l)
  modsSim[[l]] = tryCatch(
    oneRep(z[1:n], Z, S, theta0, mllk_sim, simHMM, fitMML, tol = 0.01),
    error = function(e) "An error occurred")
}
saveRDS(modsSim, "./simulation_experiments/mods/mods_2000_2.rds")


# n = 5000
n = 5000

gam_pre = mgcv::gam(y ~ s(z, bs = "ps", k = nb+1), 
                    data = data.frame(y = 1, dummy = 1, z = z[1:n]),
                    knots = list(seq(0, 1, length = nb)),
                    fit = FALSE)
Z = gam_pre$X
S = gam_pre$S[[1]]

# simulation

simiter = 230
modsSim = list()
for(l in 1:simiter) {
  cat("\nModel", l)
  modsSim[[l]] = tryCatch(
    oneRep(z[1:n], Z, S, theta0, mllk_sim, simHMM, fitMML, tol = 0.01),
    error = function(e) "An error occurred")
}
saveRDS(modsSim, "./simulation_experiments/mods/mods_5000_2.rds")






# Visualization -----------------------------------------------------------

zseq = seq(0, 1, length = 200)
Z_plot = gam_pre = mgcv::gam(y ~ s(z, bs = "ps", k = nb+1), 
                             data = data.frame(y = 1, dummy = 1, z = zseq),
                             knots = list(seq(0, 1, length = nb)),
                             fit = FALSE)$X
Gamma = simHMM(zseq)$Gamma

modsSim1000 = readRDS("./simulation_experiments/mods/mods_1000_2.rds")
ind = which(sapply(modsSim1000, is.list))
250-length(ind) # 20 models yielded an error
modsSim1000 = modsSim1000[ind[1:200]]

modsSim2000 = readRDS("./simulation_experiments/mods/mods_2000_2.rds")
ind2 = which(sapply(modsSim2000, is.list))
230-length(ind2) # 4 models yielded an error
modsSim2000 = modsSim2000[ind2[1:200]]

modsSim5000 = readRDS("./simulation_experiments/mods/mods_5000_2.rds")
ind3 = which(sapply(modsSim5000, is.list))
230-length(ind3) # 2 models yielded an error
modsSim5000 = modsSim5000[ind3[1:200]]

allmods = list(modsSim1000, modsSim2000, modsSim5000)
Ts = c(1000, 2000, 5000)



# Results for gamma_12 ----------------------------------------------------

pdf("./simulation_experiments/figs/transprobs_simulation.pdf", width = 7, height = 4.5)

par(mfrow = c(2,3), mar = c(5,4.4,2,1)+0.1)
i = 1; j = 2
for(m in 1:3){
  modsSim = allmods[[m]]
  
  plot(zseq, Gamma[i,j,], type = "l", lwd = 1, col = "black", main = paste0("T = ", Ts[m]),
       ylab = expression(gamma[12]^(t)), xlab = expression(z[t]), bty = "n", ylim = c(0, 1))
  for(k in 1:length(modsSim)){
    beta_k = matrix(modsSim[[k]]$mod$estimate[5 + 1:(2*(nb+1))], ncol = 2)
    Gamma_k = LaMa::tpm_g(Z_plot[,-1], t(beta_k))
    
    lines(zseq, Gamma_k[i,j,], col = scales::alpha("orange", 0.1), lwd = 1)
  }
  lines(zseq, Gamma[i,j,], lwd = 1.5, col = "black")
  
}
legend(-0.05, 1, legend = c("true function", "estimates"), col = c("black", "orange"), lwd = c(1, 1), bty = "n")
# dev.off()


# Trace plots

for(m in 1:3){
  modsSim = allmods[[m]]
  
  lambdas = modsSim[[1]]$allLambdas[,2]
  plot(lambdas, type = "l", col = scales::alpha("orange", 0.1), bty = "n", xlim = c(0,30), ylim = c(0, 2000),
       ylab = "penalty strength", xlab = "iteration")
  
  for(k in 2:length(modsSim)){
    lambdas = modsSim[[k]]$allLambdas[,2]
    lines(lambdas, col = scales::alpha("orange", 0.1))
  }
  
  lastlambdas = rep(NA, length(modsSim))
  for(k in 1:length(modsSim)){
    lastlambdas[k] = modsSim[[k]]$lambda[2]
  }
  
  meaniter = 0
  for(k in 1:length(modsSim)) meaniter = meaniter + nrow(modsSim[[k]]$allLambdas)
  meaniter = meaniter / length(modsSim)
  
  abline(v = meaniter, col = "black", lwd = 1, lty = 3)
  
  abline(h = mean(lastlambdas), col = "black", lwd = 1, lty = 2)
  
}
# legend("topright", legend = c("mean final penalty strength"), col = c("black"), lwd = c(2), lty = c(2), bty = "n")


dev.off()

# we see that convergence is much faster and much better for more data, as expected
# more stable, because Laplace approximation gets better with more data


# Results for gamma_21 ----------------------------------------------------

pdf("./simulation_experiments/figs/transprobs_simulation2.pdf", width = 7, height = 4.5)

par(mfrow = c(2,3), mar = c(5,4.4,2,1)+0.1)
i = 2; j = 1
for(m in 1:3){
  modsSim = allmods[[m]]
  
  plot(zseq, Gamma[i,j,], type = "l", lwd = 1, col = "black", main = paste0("T = ", Ts[m]),
       ylab = expression(gamma[21]^(t)), xlab = expression(z[t]), bty = "n", ylim = c(0, 1))
  for(k in 1:length(modsSim)){
    beta_k = matrix(modsSim[[k]]$mod$estimate[5 + 1:(2*(nb+1))], ncol = 2)
    Gamma_k = LaMa::tpm_g(Z_plot[,-1], t(beta_k))
    
    lines(zseq, Gamma_k[i,j,], col = scales::alpha("deepskyblue", 0.1), lwd = 1)
  }
  lines(zseq, Gamma[i,j,], lwd = 1.5, col = "black")
  
}
legend("topright", legend = c("true function", "estimates"), col = c("black", "deepskyblue"), lwd = c(1, 1), bty = "n")


# Trace plots

for(m in 1:3){
  modsSim = allmods[[m]]
  
  lambdas = modsSim[[1]]$allLambdas[,1]
  plot(lambdas, type = "l", col = scales::alpha("deepskyblue", 0.1), bty = "n", xlim = c(0,30), ylim = c(0, 500),
       ylab = "penalty strength", xlab = "iteration")
  
  for(k in 2:length(modsSim)){
    lambdas = modsSim[[k]]$allLambdas[,1]
    lines(lambdas, col = scales::alpha("deepskyblue", 0.1))
  }
  
  lastlambdas = rep(NA, length(modsSim))
  for(k in 1:length(modsSim)){
    lastlambdas[k] = modsSim[[k]]$lambda[1]
  }
  meaniter = 0
  for(k in 1:length(modsSim)) meaniter = meaniter + nrow(modsSim[[k]]$allLambdas)
  meaniter = meaniter / length(modsSim)
  
  abline(v = meaniter, col = "black", lwd = 1, lty = 3)
  
  abline(h = mean(lastlambdas), col = "black", lwd = 1, lty = 2)
  
}
# legend("topright", legend = c("mean final penalty strength"), col = c("black"), lwd = c(2), lty = c(2), bty = "n")

dev.off()

