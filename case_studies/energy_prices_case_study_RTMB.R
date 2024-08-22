## Spanish energy prices case study

# install.packages("RTMB")
library(RTMB) # automatic differentiation
# install.packages("MSwM")
library(MSwM) # data
# devtools::install_github("janoleko/LaMa") # development version
library(LaMa) # HMM functions and pql
# install.packages("scales")
library(scales) # just for transparency in plotting


## data
data(energy, package = "MSwM")
nrow(energy)
head(energy)

## colors
color = c("orange", "deepskyblue")


# Fitting MS-GAMLSS -------------------------------------------------------

## penalized likelihood
pnll = function(par) {
  getAll(par, dat)

  Gamma = tpm(eta) # computing the tpm
  delta = stationary(Gamma) # stationary distribution

  beta = cbind(beta0, betaspline) # mean parameter matrix
  alpha = cbind(alpha0, alphaspline) # sd parameter matrix

  Mu = Z %*% t(beta) # mean
  Sigma = exp(Z %*% t(alpha)) # sd

  allprobs = cbind(dnorm(price, Mu[,1], Sigma[,1]),
                   dnorm(price, Mu[,2], Sigma[,2])) # state-dependent densities

  REPORT(beta) # reporting for later use
  REPORT(alpha) # reporting for later use

  - forward(delta, Gamma, allprobs, ad = T) +
    penalty(list(betaspline, alphaspline), S, lambda) # penalty does the heavy lifting and reports to pql
}

## prepwork: model matrices (uses mgcv under the hood)
nb = 12 # number of basis functions
modmat = make_matrices(~ s(Oil, k = nb, bs = "ps"), energy) # model matrices Z and S
Z = modmat$Z # design matrix
S = modmat$S # penalty matrix (list)

## prediction matrix for visualizing state-dependent distributions later
xseq = seq(min(energy$Oil), max(energy$Oil), length = 200) # sequence for prediction
Z_pred = pred_matrix(modmat, newdata = data.frame(Oil = xseq))

## initial parameter list
par = list(eta = rep(-4, 2),
           beta0 = c(2, 5),
           betaspline = matrix(0, nrow = 2, ncol = nb-1),
           alpha0 = c(0, 0),
           alphaspline = matrix(0, nrow = 2, ncol = nb-1))

## data, model matrices and initial penalty strength
dat = list(price = energy$Price, Z = Z, S = S, 
           lambda = rep(1e5, 4))

## model fit
system.time(
  mod <- pql(pnll, par, dat, random = c("betaspline", "alphaspline"), 
             silent = 1, saveall = TRUE)
)

## estimated penalty strength
round(mod$lambda, 2)


## extracting parameters
beta = mod$beta # mean parameter matrix
alpha = mod$alpha # sd parameter matrix
(Gamma = mod$Gamma) # t.p.m.
round(c(Gamma[1,2], Gamma[2,1])^-1) # mean dwell time
(delta = mod$delta) # stationary distribution
mod$states = viterbi(mod$delta, mod$Gamma, mod$allprobs) # decoding most probable state sequence

Mu_plot = Z_pred %*% t(beta)
Sigma_plot = exp(Z_pred %*% t(alpha))


## visualizing results
# pdf("./case_studies/figs/energy_oil.pdf", width = 8, height = 4.5)

par(mfrow = c(1,2), mar = c(5,4,3,1))
plot(energy$Oil, energy$Price, pch = 20, bty = "n", col = alpha(color[mod$states], 0.1),
     xlab = "oil price", ylab = "energy price")
for(j in 1:2) lines(xseq, Mu_plot[,j], col = color[j], lwd = 3)

qseq = qnorm(seq(0.5, 0.95, length = 4))
for(i in qseq){
  for(j in 1:2){
    lines(xseq, Mu_plot[,j] + i*Sigma_plot[,j], col = alpha(color[j], 0.7), lwd = 1, lty = 2)
    lines(xseq, Mu_plot[,j] - i*Sigma_plot[,j], col = alpha(color[j], 0.7), lwd = 1, lty = 2)
  }
}

legend("topright", bty = "n", legend = paste("state", 1:2), col = color, lwd = 3)

plot(NA, xlim = c(0, nrow(energy)), ylim = c(1,10), bty = "n",
     xlab = "time", ylab = "energy price")
segments(x0 = 1:(nrow(energy)-1), x1 = 2:nrow(energy),
         y0 = energy$Price[-nrow(energy)], y1 = energy$Price[-1], col = color[mod$states[-1]], lwd = 0.5)

# dev.off()


## plotting the model sequence
mods = mod$allmods
length(mods)
plotind = c(1,2,3,5,10,length(mods))

# pdf("./case_studies/figs/energy_oil_modseq.pdf", width = 8, height = 5.5)

par(mfrow = c(2,3), mar = c(5,4,3,1))
for(m in plotind){
  delta = mods[[m]]$delta
  beta = mods[[m]]$beta
  alpha = mods[[m]]$alpha
  
  Mu = Z %*% t(beta)
  Sigma = exp(Z %*% t(alpha))
  allprobs = cbind(dnorm(energy$Price, Mu[,1], Sigma[,1]),dnorm(energy$Price, Mu[,2], Sigma[,2]))
  states = viterbi(delta, mods[[m]]$Gamma, mods[[m]]$allprobs)
  
  Mu_plot = Z_pred %*% t(beta)
  Sigma_plot = exp(Z_pred %*% t(alpha))
  
  plot(energy$Oil, energy$Price, pch = 20, bty = "n", col = alpha(color[states], 0.1),
       xlab = "oil price", ylab = "energy price", main = paste("Iteration", m))
  for(j in 1:2) lines(xseq, Mu_plot[,j], col = color[j], lwd = 2)
}

# dev.off()




# Using full Laplace method -----------------------------------------------

## joint likelihood (has complete normal density for random effects)
jnll = function(par) {
  getAll(par, dat)
  
  Gamma = tpm(eta) # computing the tpm
  delta = stationary(Gamma) # stationary distribution
  
  beta = cbind(beta0, betaspline) # mean parameter matrix
  alpha = cbind(alpha0, alphaspline) # sd parameter matrix
  
  Mu = Z %*% t(beta) # mean
  Sigma = exp(Z %*% t(alpha)) # sd
  
  allprobs = cbind(dnorm(price, Mu[,1], Sigma[,1]), 
                   dnorm(price, Mu[,2], Sigma[,2])) # state-dependent densities
  
  REPORT(beta) # reporting for later use
  REPORT(alpha) # reporting for later use
  
  lambda = exp(loglambda)
  # rellk = 0
  # for(i in 1:2) rellk = rellk + dgmrf(betaspline[i,], 0, lambda[i]*S, log = T)
  # for(i in 1:2) rellk = rellk + dgmrf(alphaspline[i,], 0, lambda[2+i]*S, log = T)
  
  - forward(delta, Gamma, allprobs, ad = T) - 
    sum(dgmrf2(rbind(betaspline, alphaspline), 0, S, lambda, log = TRUE))
}


## prepwork: model matrices (uses mgcv under the hood)
nb = 12 # number of basis functions
modmat = make_matrices(~ s(Oil, k = nb, bs = "ps"), energy) # model matrices Z and S
Z = modmat$Z # design matrix
S = modmat$S # penalty matrix (list)

## prediction matrix for visualizing state-dependent distributions later
xseq = seq(min(energy$Oil), max(energy$Oil), length = 200) # sequence for prediction
Z_pred = pred_matrix(modmat, newdata = data.frame(Oil = xseq))

## initial parameter list
par = list(eta = rep(-3, 2),
           beta0 = c(2, 5),
           betaspline = matrix(0, nrow = 2, ncol = nb-1),
           alpha0 = c(0, 0),
           alphaspline = matrix(0, nrow = 2, ncol = nb-1),
           loglambda = rep(log(1000), 4))

## data, model matrices and initial penalty strength
dat = list(price = energy$Price, Z = Z, S = S[[1]])

## model fitting
obj2 = MakeADFun(jnll, par, random = c("betaspline", "alphaspline"))

system.time(
  opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr) 
)

sdr = sdreport(obj2)
round(exp(as.list(sdr, "Estimate")$loglambda), 2) # estimated lambdas 

mod2 = obj2$report()

beta = mod2$beta
alpha = mod2$alpha
mod2$states = viterbi(mod2$delta, mod2$Gamma, mod2$allprobs) # decoding most probable state sequence

Mu_plot = Z_pred %*% t(beta)
Sigma_plot = exp(Z_pred %*% t(alpha))


## visualizing results
# pdf("./case_studies/figs/energy_oil.pdf", width = 8, height = 4.5)

par(mfrow = c(1,2), mar = c(5,4,3,1))
plot(energy$Oil, energy$Price, pch = 20, bty = "n", col = alpha(color[mod2$states], 0.1),
     xlab = "oil price", ylab = "energy price")
for(j in 1:2) lines(xseq, Mu_plot[,j], col = color[j], lwd = 3)

qseq = qnorm(seq(0.5, 0.95, length = 4))
for(i in qseq){
  for(j in 1:2){
    lines(xseq, Mu_plot[,j] + i*Sigma_plot[,j], col = alpha(color[j], 0.7), lwd = 1, lty = 2)
    lines(xseq, Mu_plot[,j] - i*Sigma_plot[,j], col = alpha(color[j], 0.7), lwd = 1, lty = 2)
  }
}

legend("topright", bty = "n", legend = paste("state", 1:2), col = color, lwd = 3)

plot(NA, xlim = c(0, nrow(energy)), ylim = c(1,10), bty = "n",
     xlab = "time", ylab = "energy price")
segments(x0 = 1:(nrow(energy)-1), x1 = 2:nrow(energy),
         y0 = energy$Price[-nrow(energy)], y1 = energy$Price[-1], col = color[mod2$states[-1]], lwd = 0.5)
