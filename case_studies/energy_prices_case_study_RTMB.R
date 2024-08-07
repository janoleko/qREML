
# Case study: Spanish energy prices ---------------------------------------

## packages

# install.packages("RTMB")
library(RTMB) # automatic differentiation
# install.packages("MSwM")
library(MSwM) # data
# install.packages("LaMa")
library(LaMa) # HMM functions
# install.packages("mgcv")
library(mgcv) # spline design and penalty matrices


# Loading the data --------------------------------------------------------

data(energy, package = "MSwM")


# EDA ---------------------------------------------------------------------

nrow(energy)
head(energy)


# Model fitting -----------------------------------------------------------

mllk = function(par) {
  getAll(par, dat)
  
  Gamma = tpm(eta) # tpm
  delta = stationary(Gamma) # stationary distribution
  beta = matrix(beta, ncol = 2) # spline coefs mean
  alpha = matrix(alpha, ncol = 2) # spline coefs sd
  Mu = cbind(1, Z) %*% beta # mean
  Sigma = exp(cbind(1, Z) %*% alpha) # sd
  allprobs = cbind(dnorm(price, Mu[,1], Sigma[,1]), dnorm(price, Mu[,2], Sigma[,2])) # state-dependent densities
  
  pen = 0 # penalty
  for(i in 1:2) pen = pen + lambda[i] * (t(as.numeric(beta[-1,i])) %*% S %*% as.numeric(beta[-1,i])) # quadform
  for(i in 1:2) pen = pen + lambda[2+i] * (t(as.numeric(alpha[-1,i])) %*% S %*% as.numeric(alpha[-1,i])) # quadform
  
  REPORT(beta)
  REPORT(alpha)
  
  - forward(delta, Gamma, allprobs, ad = T) + 0.5 * pen
}


# design and penalty matrix
nb = 12
xseq = seq(min(energy$Oil), max(energy$Oil), length = 200) # sequence for prediction
gam_setup = gam(dummy ~ s(Oil, k = nb, bs = "ps"), 
                data = data.frame(dummy = 1, Oil = c(xseq, energy$Oil)), fit = FALSE)
Z_plot = gam_setup$X[1:200,] # design matrix for prediction
Z = gam_setup$X[-(1:200),-1] # design matrix for model fit
S = gam_setup$S[[1]]

# checking the rank deficiency for correction in algorithm
(m = nrow(S) - Matrix::rankMatrix(S)[1])

S_sparse = as(S, "sparseMatrix")

# Optimisation via marginal ML --------------------------------------------

maxiter = 50 # maximum number of iterations
tol = 0.01 # relative tolerance for convergence, sufficient 
Lambdas = matrix(NA, maxiter+1, 4)
Lambdas[1,] = c(rep(10000, 4)) # initial penalty strengths
mods = list()
# object that contains the indices for each random effect/ spline coef vector. If these are of differnt lengths, use list and change below
REind = matrix(2 + 1:(4*(nb)), nrow = 4, byrow = TRUE)[,-1] # each row is the index of one RE
alpha_sm = 1 # exponential smoothing parameter

# initial parameter
par = list(eta = rep(-4, 2), 
           # beta = c(4.5, seq(-2, 2, length = nb-1), 6.5, seq(-1, 1, length = nb-1)),
           beta = c(2, rep(0, nb-1), 5, rep(0, nb-1)),
           alpha = rep(0, 2*nb))
dat = list(price = energy$Price, Z = Z, S = S, lambda = Lambdas[1,])

# making AD loss function
obj = MakeADFun(mllk, par, silent = TRUE) 

# updating algorithm
T1 = Sys.time()
for(k in 1:maxiter){
  cat("\n\n- Iteration", k, "-")
  t1 = Sys.time()
  obj = MakeADFun(mllk, par, silent = TRUE) 
  opt = nlminb(obj$par, obj$fn, obj$gr, control = list(rel.tol = 1e-10))
  cat("\nEstimation time:", Sys.time()-t1)
  mod = obj$report()
  mods[[k]] = mod
  coefs = cbind(mod$beta, mod$alpha)
  J = obj$he()
  J_inv = solve(J)
  for(i in 1:4){ # updating all penalty strengths
    edoF = nrow(S) - sum(diag(Lambdas[k,i] * J_inv[REind[i,], REind[i,]] %*% S))
    penalty = t(as.numeric(coefs[-1,i])) %*% S %*% as.numeric(coefs[-1,i])
    lambda_new = as.numeric((edoF - m) / penalty) # m is correction if S does not have full rank
    Lambdas[k+1, i] = alpha_sm * lambda_new + (1-alpha_sm) * Lambdas[k, i]
  }
  dat$lambda = Lambdas[k+1,]
  par$beta = mod$beta; par$alpha = mod$alpha
  
  cat("\nSmoothing strengths:", round(Lambdas[k+1, ], 4))
  if(max(abs(Lambdas[k+1,] - Lambdas[k,]) / Lambdas[k,]) < tol){
    cat("\n\n")
    break
  }
}
Sys.time()-T1

Lambdas = as.matrix(na.omit(Lambdas))
lambda_hat = Lambdas[nrow(Lambdas),]
par(mfrow = c(1,4))
# for(j in 1:4){
#   plot(Lambdas[,j], type = "l", lwd = 2, main = paste("lambda", j))
# }

mod_final = mods[[length(mods)]]

beta = mod_final$beta
alpha = mod_final$alpha
states = viterbi(mod_final$delta, mod_final$Gamma, mod_final$allprobs) # decoding most probable state sequence

# xseq = seq(min(energy$Oil), max(energy$Oil), length = 200)
# Z_plot = cbind(1, splines::splineDesign(knots, xseq, outer.ok = TRUE)[,-1])

# Generate the new design matrix using the PredictMat function
# Z_plot = cbind(1, PredictMat(gam_setup$smooth[[1]], data.frame(dummy = 1, Price = 1, Oil = xseq)))

Mu_plot = Z_plot %*% t(beta)
Sigma_plot = exp(Z_plot %*% alpha)

color = c("orange", "deepskyblue")

# pdf("./case_studies/figs/energy_oil.pdf", width = 8, height = 4.5)

par(mfrow = c(1,2), mar = c(5,4,3,1))
plot(energy$Oil, energy$Price, pch = 20, bty = "n", col = scales::alpha(color[states], 0.1),
     xlab = "oil price", ylab = "energy price")
for(j in 1:2) lines(xseq, Mu_plot[,j], col = color[j], lwd = 3)

qseq = qnorm(seq(0.5, 0.95, length = 4))
for(i in qseq){
  for(j in 1:2){
    lines(xseq, Mu_plot[,j] + i*Sigma_plot[,j], col = scales::alpha(color[j], 0.7), lwd = 1, lty = 2)
    lines(xseq, Mu_plot[,j] - i*Sigma_plot[,j], col = scales::alpha(color[j], 0.7), lwd = 1, lty = 2)
  }
}

legend("topright", bty = "n", legend = paste("state", 1:2), col = color, lwd = 3)

plot(NA, xlim = c(0, nrow(energy)), ylim = c(1,10), bty = "n",
     xlab = "time", ylab = "energy price")
segments(x0 = 1:(nrow(energy)-1), x1 = 2:nrow(energy),
         y0 = energy$Price[-nrow(energy)], y1 = energy$Price[-1], col = color[states[-1]], lwd = 0.5)

# dev.off()

# par(mfrow = c(1,1))
# plot(Sigma_plot[,1], ylim = c(0,1.5), type = "l", col = color[1], lwd = 2, bty = "n",
#      xlab = "z", ylab = "standard deviation")
# lines(Sigma_plot[,2], col = color[2], lwd = 2)



# Plotting the model sequence ---------------------------------------------

length(mods)
plotind = c(1,2,3,5,7,length(mods))

# xseq = seq(min(energy$Oil), max(energy$Oil), length = 200)
# Z_plot = cbind(1, splines::splineDesign(knots, xseq, outer.ok = TRUE)[,-1])

# pdf("./case_studies/figs/energy_oil_modseq.pdf", width = 8, height = 5.5)

par(mfrow = c(2,3), mar = c(5,4,3,1))
for(m in plotind){
  delta = mods[[m]]$delta
  Gamma = mods[[m]]$Gamma
  beta = mods[[m]]$beta
  alpha = mods[[m]]$alpha
  
  Mu = cbind(1, Z) %*% beta
  Sigma = exp(cbind(1, Z) %*% alpha)
  allprobs = cbind(dnorm(energy$Price, Mu[,1], Sigma[,1]),dnorm(energy$Price, Mu[,2], Sigma[,2]))
  states = viterbi(c(delta), Gamma, allprobs)
  
  Mu_plot = Z_plot %*% beta
  Sigma_plot = exp(Z_plot %*% alpha)
  
  plot(energy$Oil, energy$Price, pch = 20, bty = "n", col = scales::alpha(color[states], 0.1),
       xlab = "oil price", ylab = "energy price", main = paste("Iteration", m))
  for(j in 1:2) lines(xseq, Mu_plot[,j], col = color[j], lwd = 2)
}

#dev.off()


mllk2 = function(par) {
  getAll(par, dat)
  
  Gamma = tpm(eta) # tpm
  delta = stationary(Gamma) # stationary distribution
  beta = matrix(c(beta0, betavec), ncol = 2, byrow = T) # spline coefs mean
  alpha = matrix(c(alpha0, alphavec), ncol = 2, byrow = T) # spline coefs sd
  Mu = cbind(1, Z) %*% beta # mean
  Sigma = exp(cbind(1, Z) %*% alpha) # sd
  allprobs = cbind(dnorm(price, Mu[,1], Sigma[,1]), dnorm(price, Mu[,2], Sigma[,2])) # state-dependent densities
  
  REPORT(beta)
  REPORT(alpha)
  
  l = forward(delta, Gamma, allprobs, ad = TRUE)
  l = l + sum(dgmrf(t(beta[-1,]), 0, S, log = TRUE))
  l = l + sum(dgmrf(t(alpha[-1,]), 0, S, log = TRUE))
  -l
}

# initial parameter
par = list(eta = rep(-4, 2), 
           beta0 = c(2, 5), betavec = rep(0, 2*(nb-1)),
           alpha0 = c(0,0), alphavec = rep(0, 2*(nb-1)))
dat = list(price = energy$Price, Z = Z, S = S_sparse, lambda = Lambdas[1,])

# making AD loss function
obj2 = MakeADFun(mllk2, par, random = c("betavec", "alphavec")) 

opt2 = nlminb(obj2$par, obj2$fn, obj2$gr) 
mod2 = obj2$report()

beta = mod2$beta
alpha = mod2$alpha
states = viterbi(mod2$delta, mod2$Gamma, mod2$allprobs) # decoding most probable state sequence

Mu_plot = Z_plot %*% beta
Sigma_plot = exp(Z_plot %*% alpha)

color = c("orange", "deepskyblue")

# pdf("./case_studies/figs/energy_oil.pdf", width = 8, height = 4.5)

par(mfrow = c(1,2), mar = c(5,4,3,1))
plot(energy$Oil, energy$Price, pch = 20, bty = "n", col = scales::alpha(color[states], 0.1),
     xlab = "oil price", ylab = "energy price")
for(j in 1:2) lines(xseq, Mu_plot[,j], col = color[j], lwd = 3)

qseq = qnorm(seq(0.5, 0.95, length = 4))
for(i in qseq){
  for(j in 1:2){
    lines(xseq, Mu_plot[,j] + i*Sigma_plot[,j], col = scales::alpha(color[j], 0.7), lwd = 1, lty = 2)
    lines(xseq, Mu_plot[,j] - i*Sigma_plot[,j], col = scales::alpha(color[j], 0.7), lwd = 1, lty = 2)
  }
}

legend("topright", bty = "n", legend = paste("state", 1:2), col = color, lwd = 3)






###########################################################################
# With pql function -------------------------------------------------------

## packages
# install.packages("RTMB")
library(RTMB) # automatic differentiation
# install.packages("MSwM")
library(MSwM) # data set
# install.packages("LaMa")
library(LaMa) # HMM functions
# install.packages("mgcv")
library(mgcv) # spline design and penalty matrices


# Loading the data --------------------------------------------------------

data(energy, package = "MSwM")


# Writing the penalized negative log-likelihood function ------------------

mllk = function(par) {
  getAll(par, dat)
  
  Gamma = tpm(eta) # computing the tpm
  delta = stationary(Gamma) # stationary distribution
  
  beta = cbind(beta0, betaspline) # mean parameter matrix
  alpha = cbind(alpha0, alphaspline) # sd parameter matrix

  Mu = Z %*% t(beta) # mean
  Sigma = exp(Z %*% t(alpha)) # sd
  
  allprobs = cbind(dnorm(price, Mu[,1], Sigma[,1]), 
                   dnorm(price, Mu[,2], Sigma[,2])) # state-dependent densities

  REPORT(beta) # reporting parameter matrices for later use
  REPORT(alpha)
  
  - forward(delta, Gamma, allprobs, ad = T) + 
    penalty(list(betaspline, alphaspline), S, lambda) # penalty function does the heavy lifting
}



# Building model matrices with mgcv ---------------------------------------

nb = 12 # number of basis functions
modmat = make_matrices(~ s(Oil, k = nb, bs = "ps"), energy) # model matrices Z and S
Z = modmat$Z
S = modmat$S

# sequence for prediction
xseq = seq(min(energy$Oil), max(energy$Oil), length = 200) # sequence for prediction
Z_pred = pred_matrix(modmat, newdata = data.frame(Oil = xseq))



# Model fit ---------------------------------------------------------------

# initial parameter list
par = list(eta = rep(-4, 2),
           beta0 = c(2, 5),
           betaspline = matrix(0, nrow = 2, ncol = nb-1),
           alpha0 = c(0, 0),
           alphaspline = matrix(0, nrow = 2, ncol = nb-1))

# data, model matrices and initial penalty strength
dat = list(price = energy$Price, Z = Z, S = S, lambda = list(rep(100, 2), rep(100, 2)))

mod = pql(mllk, par, dat, random = c("betaspline", "alphaspline"))


# Visualizing results -----------------------------------------------------

beta = mod$beta
alpha = mod$alpha
states = viterbi(mod$delta, mod$Gamma, mod$allprobs) # decoding most probable state sequence

Mu_plot = Z_pred %*% t(beta)
Sigma_plot = exp(Z_pred %*% t(alpha))

color = c("orange", "deepskyblue")

# pdf("./case_studies/figs/energy_oil.pdf", width = 8, height = 4.5)

par(mfrow = c(1,2), mar = c(5,4,3,1))
plot(energy$Oil, energy$Price, pch = 20, bty = "n", col = scales::alpha(color[states], 0.1),
     xlab = "oil price", ylab = "energy price")
for(j in 1:2) lines(xseq, Mu_plot[,j], col = color[j], lwd = 3)

qseq = qnorm(seq(0.5, 0.95, length = 4))
for(i in qseq){
  for(j in 1:2){
    lines(xseq, Mu_plot[,j] + i*Sigma_plot[,j], col = scales::alpha(color[j], 0.7), lwd = 1, lty = 2)
    lines(xseq, Mu_plot[,j] - i*Sigma_plot[,j], col = scales::alpha(color[j], 0.7), lwd = 1, lty = 2)
  }
}

legend("topright", bty = "n", legend = paste("state", 1:2), col = color, lwd = 3)

plot(NA, xlim = c(0, nrow(energy)), ylim = c(1,10), bty = "n",
     xlab = "time", ylab = "energy price")
segments(x0 = 1:(nrow(energy)-1), x1 = 2:nrow(energy),
         y0 = energy$Price[-nrow(energy)], y1 = energy$Price[-1], col = color[states[-1]], lwd = 0.5)


