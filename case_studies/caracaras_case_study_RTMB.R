## Caracaras case study

## packages
# install.packages("RTMB")
library(RTMB)
# devtools::install_github("janoleko/LaMa") # development version
library(LaMa)
# install.packages("scales")
library(scales) # for transparent colors

## data
data = read.csv("./data/caracara_downsampled.csv")
head(data)
nrow(data)

## color palette
color = c("orange", "deepskyblue", "seagreen3")


# Simple parametric HMM ---------------------------------------------------

## likelihood function for simple homogeneous HMM
mllk_sim = function(par){
  getAll(par, dat)
  sigma = exp(logsigma)
  Gamma = tpm(eta)
  delta = stationary(Gamma)
  
  REPORT(mu)
  REPORT(sigma)
  
  allprobs = matrix(1, length(x), N)
  ind = which(!is.na(x))
  for(j in 1:N) allprobs[ind, j] = dnorm(x[ind], mu[j], sigma[j])
  
  -forward(delta, Gamma, allprobs, ad = T)
}

## first 3 state model
N = 3
par = list(mu = c(-5, -4, -2), # state-dep mean
           logsigma = rep(log(0.3), 3), # state-dep sd
           eta = rep(-2, 6)) # tpm params

dat = list(x = data$logVDBA, N = 3)

## fitting model
obj3 = MakeADFun(mllk_sim, par, silent = TRUE)
opt3 = nlminb(obj3$par, obj3$fn, obj3$gr)

## extracting parameters
mod3 = obj3$report()
delta3 = mod3$delta
mu3 = mod3$mu
sigma3 = mod3$sigma

## plotting estimated state-dependent distributions

#pdf("./case_studies/figs/caracaras_simple.pdf", width = 8.5, height = 4)

par(mfrow = c(1,2))

hist(data$logVDBA, breaks = 50, prob = T, bor = "white", main = "", xlab = "log(VeDBA)", ylab = "density")
for(j in 1:N) curve(delta3[j] * dnorm(x, mu3[j], sigma3[j]), add = T, col = color[j], lwd = 2, n = 200)
curve(delta3[1] * dnorm(x, mu3[1], sigma3[1]) + 
        delta3[2] * dnorm(x, mu3[2], sigma3[2]) + 
        delta3[3] * dnorm(x, mu3[3], sigma3[3]), add = T, col = "black", lty = 2, lwd = 2, n = 200)
legend("topright", legend = c("resting", "low activity", "high activity"), 
       col = color, lwd = 2, bty = "n")


mod3$AIC = 2 * obj3$fn() + 2 * length(obj3$par)
mod3$BIC = 2 * obj3$fn() + log(nrow(data)) * length(obj3$par)


## then 4 state model
N = 4
par = list(mu = c(-5, -4.5, -3, -1.5), # state-dep mean
           logsigma = rep(log(0.3), 4), # state-dep sd
           eta = rep(-2, 12)) # tpm params

dat = list(x = data$logVDBA, N = 4)

## fitting model
obj4 = MakeADFun(mllk_sim, par, silent = TRUE)
opt4 = nlminb(obj4$par, obj4$fn, obj4$gr)

## extracting parameters
mod4 = obj4$report()
delta4 = mod4$delta
mu4 = mod4$mu
sigma4 = mod4$sigma

## plotting estimated state-dependent distributions
color2 = c("#FFEF00", "orange", "deepskyblue", "seagreen4")

hist(data$logVDBA, breaks = 50, prob = T, bor = "white", main = "", xlab = "log(VeDBA)", ylab = "density")
for(j in 2:N) curve(delta4[j] * dnorm(x, mu4[j], sigma4[j]), add = T, col = color2[j], lwd = 2, n = 200)
curve(delta4[1] * dnorm(x, mu4[1], sigma4[1]), add = T, col = color2[1], lwd = 2, n = 200)
curve(delta4[1] * dnorm(x, mu4[1], sigma4[1]) + 
        delta4[2] * dnorm(x, mu4[2], sigma4[2]) + 
        delta4[3] * dnorm(x, mu4[3], sigma4[3])+
        delta4[4] * dnorm(x, mu4[4], sigma4[4]), add = T, col = "black", lty = 2, lwd = 2, n = 200)
legend("topright", legend = c("resting 1", "resting 2", "low activity", "high activity"), 
       col = color2, lwd = 2, bty = "n")

mod4$AIC = 2 * obj4$fn() + 2 * length(obj4$par)
mod4$BIC = 2 * obj4$fn() + log(nrow(data)) * length(obj4$par)

#dev.off()




# Nonparametric emission distributions ------------------------------------

## penalized likelihood
pnll = function(par) {
  getAll(par, dat)
  
  Gamma = tpm(eta) # transition probability matrix
  delta = stationary(Gamma) # stationary distribution
  
  ## computing mixture weights
  A = exp(cbind(rep(0, N), beta)) # first column fixed at zero
  A = A / rowSums(A) # multinomial logit link
  
  ## computing state-dependent distributions
  allprobs = matrix(1, nrow = nrow(Z), ncol = N)
  ind = which(!is.na(x))
  allprobs[ind,] = Z[ind,] %*% t(A) # linear: sum over weighted basis functions
  
  REPORT(A) # report A for convenience later
  
  -forward(delta, Gamma, allprobs, ad = T) + 
    penalty(beta, S, lambda) # penalty function reports to pql
}


## prepwork: model matrices
modmat = make_matrices_dens(data$logVDBA, 25)
Z = modmat$Z # design matrix
S = modmat$S # penalty matrix
# knots = modmat$knots # knots
basis_pos = modmat$basis_pos # basis positions for parameter starting values

## prediction matrix for visualizing state-dependent distributions later
xseq = seq(min(data$logVDBA, na.rm=T), max(data$logVDBA, na.rm=T), length = 300)
Z_pred = pred_matrix_dens(modmat, xseq)

## starting values
N = 3
beta = log(rbind(dnorm(basis_pos, -5, 1.5),
                 dnorm(basis_pos, -4, 1.5),
                 dnorm(basis_pos, -2.5, 1.5)))
beta = beta - beta[,1] # subtracting first column for smooth transition to fixed zero coefs 

par = list(eta = rep(-3, N*(N-1)), beta = beta)
dat = list(x = data$logVDBA, Z = Z, N = 3, S = S, lambda = rep(30, 3))

## model fitting via pql
system.time(
  mod <- pql(pnll, par, dat, random = "beta")
)

## estimated penalty strength
round(mod$lambda, 2)

## extracting parameters
A = mod$A # weight matrix
delta = mod$delta # stationary distribution
mod$states = viterbi(mod$delta, mod$Gamma, mod$allprobs) # state decoding

## AIC and BIC
npar = mod$n_fixpar + sum(unlist(mod$edf))
(mod$AIC = -2*mod$llk + 2*npar)
(mod$BIC = -2*mod$llk + log(nrow(data))*npar)

## visualizing results

# pdf("./case_studies/figs/caracara_log.pdf", width = 8.5, height = 4)

par(mfrow = c(1,2), mar = c(5,4,3.5,1) + 0.1)

hist(data$logVDBA, prob = T, breaks = 50, bor = "white",
     main = "", xlab = "log(VeDBA)", ylab = "density")

dens = matrix(NA, nrow = length(xseq), ncol = N)
for(j in 1:N){
  dens[,j] = delta[j] * Z_pred %*% A[j,]
  lines(xseq, dens[,j], col = color[j], lwd = 2)
}
lines(xseq, rowSums(dens), col = 1, lty = 2, lwd = 2)

legend("topright", legend = c("resting", "low activity", "high activity"), 
       col = color, lwd = 2, bty = "n")


# plot(data$VDBA[plotind], type = "h", col = color[mod$states[plotind]], 
#      xlab = "time", ylab = "VeDBA", bty = "n", lwd = 0.3)

plotind = 4000:8000
plot(data$logVDBA[plotind], pch = 20, col = alpha(color[mod$states[plotind]], 0.3), 
     xlab = "time", ylab = "logVeDBA", bty = "n")

# dev.off()


















# Full Laplace with REML --------------------------------------------------

jnll = function(par) {
  getAll(par, dat)
  
  Gamma = tpm(eta) # transition probability matrix
  delta = stationary(Gamma) # stationary distribution
  
  ## computing mixture weights
  A = exp(cbind(rep(0, N), beta)) # first column fixed at zero
  A = A / rowSums(A) # multinomial logit link
  
  ## computing state-dependent distributions
  allprobs = matrix(1, nrow = nrow(Z), ncol = N)
  ind = which(!is.na(x))
  allprobs[ind,] = Z[ind,] %*% t(A) # linear: sum over weighted basis functions
  
  REPORT(A) # report A for convenience later
  
  lambda = exp(loglambda)
  
  -forward(delta, Gamma, allprobs, ad = T) - 
    sum(dgmrf2(beta, 0, S, lambda, log = TRUE))
}

## Good starting values via guessing normal distributions
N = 3
beta = log(rbind(dnorm(basis_pos, -5, 1),
                 dnorm(basis_pos, -4, 1),
                 dnorm(basis_pos, -2.5, 1)))
beta = beta - beta[,1] # subtracting first column for smooth transition to fixed zero coefs 

par = list(eta = rep(-3, N*(N-1)), 
           beta = beta, 
           loglambda = rep(log(15), 3))
dat = list(x = data$logVDBA, Z = Z, N = 3, S = S)

obj = MakeADFun(jnll, par, random = "beta")
opt = nlminb(obj$par, obj$fn, obj$gr)
obj$fn()
