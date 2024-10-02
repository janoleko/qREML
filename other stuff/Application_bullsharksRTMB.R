
## packages
library(LaMa)
library(RTMB)

## data
load("~/Downloads/bull_sharks_workspace_for_plots.RData")
sharks = x
IDs = unique(sharks$SharkID) # unique IDs

hist(sharks$log_ODBA, breaks = 100, prob = TRUE)

## color
color = c("orange", "deepskyblue")


### simple model
mllk_sim = function(par){
  getAll(par, dat)
  sigma = exp(logsigma)
  REPORT(mu); REPORT(sigma)
  
  Gamma = tpm(eta)
  delta = stationary(Gamma)
  
  allprobs = matrix(1, nrow = length(log_ODBA), ncol = N)
  ind = which(!is.na(log_ODBA))
  for (j in 1:N){
    allprobs[ind, j] = dnorm(log_ODBA[ind], mu[j], sigma[j])
  }
  -forward(delta, Gamma, allprobs, trackID, ad = T)
}

N = 2
par = list(mu = c(-0.7, -0.2), logsigma = c(0.25, 0.25), eta = rep(-3, 2))
dat = list(log_ODBA = sharks$log_ODBA, trackID = sharks$SharkID, N = 2)

## model fitting
obj_sim = MakeADFun(mllk_sim, par, silent = T)
opt_sim = nlminb(obj_sim$par, obj_sim$fn, obj_sim$gr)

## extracting parameters
mod_sim = obj_sim$report()
mu = mod_sim$mu
sigma = mod_sim$sigma
Gamma = mod_sim$Gamma
delta = mod_sim$delta

## visualizing
hist(sharks$log_ODBA, breaks = 100, prob = TRUE)
for (j in 1:N){
  curve(delta[j]*dnorm(x, mu[j], sigma[j]), add = TRUE, col = color[j], lwd = 2)
}
curve(delta[1]*dnorm(x, mu[1], sigma[1]) + delta[2]*dnorm(x, mu[2], sigma[2]), 
      add = TRUE, lwd = 2, lty = 2)


### Flexible big model as in the paper
# GAM effects over time for each individual, 
# non-parametric effects of time of day on the transition probabilities (pooled)

## penalized likelihood function
pnll = function(par){
  getAll(par, data)
  
  ## state process
  beta = cbind(beta0, betaSpline1, betaSpline2)
  Gamma = tpm_g(Z_s, beta, ad = T)
  delta = c(1, exp(logitdelta))
  delta = delta / sum(delta); REPORT(delta)
  
  ## state-dependent process (GAM)
  smooth = rowSums(Z_sd[,-1] * A[IDnum,]) # GAM means for state 1 and each individual (on log scale)
  Mu = exp(cbind(alpha0[1] + smooth, # mean state 1: intercept + smooth
                 alpha0[2] + smooth + mu2[IDnum])) # mean state 2: intercept + smooth + individual-specific deviatio
  Sigma = exp(logSigma)
  sigma = Sigma[IDnum,]
  REPORT(Mu); REPORT(Sigma); REPORT(alpha0); REPORT(A); REPORT(mu2)
  
  allprobs = matrix(1, nrow = length(ODBA), ncol = N)
  ind = which(!is.na(ODBA))
  allprobs[ind,] = cbind(dgamma2(ODBA[ind], Mu[ind,1], sigma[ind,1]),
                         dgamma2(ODBA[ind], Mu[ind,2], sigma[ind,2]))

  -forward_g(t(delta), Gamma, allprobs, trackID, ad = T) + 
    penalty(list(betaSpline1, betaSpline2, A), S, lambda)
}

## prepwork: design and penalty matrices
k_s = 10
modmat_s = make_matrices(~ s(Hour, bs = "cp", k = k_s) +
                           s(Hour, bs = "cp", k = k_s, by = AvgTemp),
                         data = sharks, knots = list(Hour = c(0, 24)))

Z_s = modmat_s$Z
S_s = modmat_s$S # 1. s(time), 2. AvgTemp * s(time)

k_sd = 15
modmat_sd = make_matrices(~ s(IDindex, bs = "ps", k = k_sd), data = sharks)
Z_sd = modmat_sd$Z
S_sd = modmat_sd$S # 1. s(IDindex)

S = c(S_s, S_sd)


## initial parameters and data
par = list(beta0 = rep(-2, 2),                           # state process intercepts
           betaSpline1 = matrix(0, 2, k_s-1),            # state process spline coef s(time)
           betaSpline2 = matrix(0, 2, k_s),              # state process spline coef AvgTemp * s(time)
           logitdelta = 0,                               # initial distribution
           alpha0 = rep(-0.5, 2),                        # intercept state-dependent mean
           A = matrix(0, length(IDs), k_sd-1),           # state-dep process spline coef for state one but each individual
           mu2 = rep(0.3, length(IDs)),                  # individual-specific deviations from state 1 mean for state 2 mean
           logSigma = matrix(log(0.16), length(IDs), 2)) # individual-specific standard deviations

data = list(ODBA = sharks$ODBA, IDnum = sharks$IDnum, trackID = sharks$SharkID, N = 2,
           Z_s = Z_s, Z_sd = Z_sd, 
           S = S,
           lambda = c(rep(1e3, 4), rep(1e5, length(IDs))))

## model fitting
system.time(
  mod <- qreml(pnll, par, data, random = c("betaSpline1", "betaSpline2", "A"),
             maxiter = 500, silent = 1)
)

length(mod$par_vec)


## extracting parameters
# state process coefficients
beta = mod$beta
# state-dependent means based on smooth for each shark
Mu = mod$Mu

## decode states
mod$states = viterbi_g(mod$delta, mod$Gamma, mod$allprobs, mod$trackID)

## visualizing results
# state-dependent process
par(mfrow = c(2,3))
for(j in 1:6){
  sharkInd = which(sharks$IDnum == j)
  plot(sharks$ODBA[sharkInd], pch = 16, col = scales::alpha(color[mod$states[sharkInd]], 0.2), 
       xlab = "time", ylab = "ODBA", main = IDs[j], bty = "n")
  lines(Mu[sharkInd, 1], col = "white", lwd = 3) # white border
  lines(Mu[sharkInd, 2], col = "white", lwd = 3) # white border
  lines(Mu[sharkInd, 1], col = color[1], lwd = 2)
  lines(Mu[sharkInd, 2], col = color[2], lwd = 2)
}

pdf("~/Desktop/Laplace approx plots/bullsharks_state_dep.pdf", width = 8, height = 3.5)

lim = c(2, 1.4, 1.6)
par(mfrow = c(1,3))
for(i in 1:3){
  j = c(1,3,6)[i]
  sharkInd = which(sharks$IDnum == j)
  plot(sharks$ODBA[sharkInd], pch = 16, col = scales::alpha(color[mod$states[sharkInd]], 0.2), 
       xlab = "time", ylab = "ODBA", main = IDs[j], bty = "n", xaxt = "n", ylim = c(0.2, lim[i]))
  lines(Mu[sharkInd, 1], col = "white", lwd = 4) # white border
  lines(Mu[sharkInd, 2], col = "white", lwd = 4) # white border
  lines(Mu[sharkInd, 1], col = color[1], lwd = 2)
  lines(Mu[sharkInd, 2], col = color[2], lwd = 2)
  
  axis(1, at = seq(0, round(max(sharks$IDindex[sharkInd]), -2), by = 600), labels = seq(0, round(max(sharks$IDindex[sharkInd]), -2), by = 600))
}

dev.off()

# state process
tod_seq = seq(0, 24, length = 250)
nlevels = 12
temp_seq = seq(20, 33, length = nlevels)
tempcolor = colorRampPalette(c("blue", "red"))(nlevels)
par(mfrow = c(1,1))

plot(NA, type = "l", xlim = c(0, 24), ylim = c(0,1), xaxt = "n",
     xlab = "time of day",ylab = "foraging probability", bty = "n")
axis(1, at = seq(0, 24, by = 6), labels = seq(0, 24, by = 6))
polygon(x = c(0, 6, 6, 0), y = c(0, 0, 1, 1), col = "gray95", border = "white")
polygon(x = c(18, 24, 24, 18), y = c(0, 0, 1, 1), col = "gray95", border = "white")

for(i in 1:length(temp_seq)){
  Z_pred_s = pred_matrix(modmat_s, data.frame(Hour = tod_seq, AvgTemp = temp_seq[i]))
  Gamma = tpm_g(Z_pred_s, beta)
  Delta = matrix(NA, dim(Gamma)[3], ncol = N)
  for(k in 1:dim(Gamma)[3]){
    Delta[k,] = LaMa::stationary(Gamma[,,k])
  }
  lines(tod_seq, Delta[,2], lwd = 2, col = scales::alpha(tempcolor[i], 0.8))
}
legend(x = 12, y = 1, border = NA, lwd = 2,
       legend = paste0(round(rev(temp_seq[c(1,3,6,9,12)])),"°"), 
       col = rev(tempcolor[c(1,3,6,9,12)]), bty = "n")

# pdf("~/Desktop/Laplace approx plots/bullsharks_state.pdf", width = 8, height = 5)

# true periodically stationary distribution
plot(NA, type = "l", xlim = c(0, 24), ylim = c(0,1), xaxt = "n",
     xlab = "time of day",ylab = "foraging probability", bty = "n")
axis(1, at = seq(0, 24, by = 6), labels = seq(0, 24, by = 6))
polygon(x = c(0, 6, 6, 0), y = c(0, 0, 1, 1), col = "gray95", border = "white")
polygon(x = c(18, 24, 24, 18), y = c(0, 0, 1, 1), col = "gray95", border = "white")

for(i in 1:length(temp_seq)){
  Z_pred_s = pred_matrix(modmat_s, data.frame(Hour = 0:23, AvgTemp = temp_seq[i]))
  Gamma = tpm_g(Z_pred_s, beta)
  Delta = stationary_p(Gamma)
  Delta = rbind(Delta, Delta[1,])
  lines(0:24, type = "b", Delta[,2], pch = 16, lwd = 1, col = scales::alpha(tempcolor[i], 0.8))
}
legend(x = 18, y = 1, border = NA, pch = 16,
       legend = paste0(round(rev(temp_seq[c(1,3,6,9,12)])),"°"), 
       col = scales::alpha(rev(tempcolor[c(1,3,6,9,12)]), 0.8), bty = "n")

# dev.off()


# Adding random intercept -------------------------------------------------

## penalized likelihood function
pnll2 = function(par){
  getAll(par, data)
  
  ## state process
  beta = cbind(beta0, betaRI, betaSpline1, betaSpline2)
  Gamma = tpm_g(Z_s, beta, ad = T)
  delta = c(1, exp(logitdelta))
  delta = delta / sum(delta); REPORT(delta)
  
  ## state-dependent process (GAM)
  smooth = rowSums(Z_sd[,-1] * A[IDnum,]) # GAM means for state 1 and each individual (on log scale)
  Mu = exp(cbind(alpha0[1] + smooth, # mean state 1: intercept + smooth
                 alpha0[2] + smooth + mu2[IDnum])) # mean state 2: intercept + smooth + individual-specific deviatio
  Sigma = exp(logSigma)
  sigma = Sigma[IDnum,]
  REPORT(Mu); REPORT(Sigma); REPORT(alpha0); REPORT(A); REPORT(mu2)
  
  allprobs = matrix(1, nrow = length(ODBA), ncol = N)
  ind = which(!is.na(ODBA))
  allprobs[ind,] = cbind(dgamma2(ODBA[ind], Mu[ind,1], sigma[ind,1]),
                         dgamma2(ODBA[ind], Mu[ind,2], sigma[ind,2]))
  
  -forward_g(t(delta), Gamma, allprobs, trackID, ad = T) + 
    penalty(list(betaRI, betaSpline1, betaSpline2, A), S, lambda)
}

## prepwork: design and penalty matrices
sharks$SharkID = as.factor(sharks$SharkID)

k_s = 10
modmat_s = make_matrices(~ s(SharkID, bs = "re") +
                           s(Hour, bs = "cp", k = k_s) +
                           s(Hour, bs = "cp", k = k_s, by = AvgTemp),
                         data = sharks, knots = list(Hour = c(0, 24)))

Z_s = modmat_s$Z
S_s = modmat_s$S # 1. SharkID, 2. s(time), 3. AvgTemp * s(time)

k_sd = 15
modmat_sd = make_matrices(~ s(IDindex, bs = "ps", k = k_sd), data = sharks)
Z_sd = modmat_sd$Z
S_sd = modmat_sd$S # 1. s(IDindex)

S = c(S_s, S_sd)


## initial parameters and data
par = list(beta0 = rep(-2, 2),                           # state process intercepts
           betaRI = matrix(0, 2, length(IDs)),           # random intercepts
           betaSpline1 = matrix(0, 2, k_s-1),            # state process spline coef s(time)
           betaSpline2 = matrix(0, 2, k_s),              # state process spline coef AvgTemp * s(time)
           logitdelta = 0,                               # initial distribution
           alpha0 = rep(-0.5, 2),                        # intercept state-dependent mean
           A = matrix(0, length(IDs), k_sd-1),           # state-dep process spline coef for state one but each individual
           mu2 = rep(0.3, length(IDs)),                  # individual-specific deviations from state 1 mean for state 2 mean
           logSigma = matrix(log(0.16), length(IDs), 2)) # individual-specific standard deviations

data = list(ODBA = sharks$ODBA, IDnum = sharks$IDnum, trackID = sharks$SharkID, N = 2,
            Z_s = Z_s, Z_sd = Z_sd, 
            S = S,
            lambda = c(rep(100, 2), rep(1e3, 4), rep(1e5, length(IDs))))

## model fitting
system.time(
  mod2 <- qreml(pnll2, par, data, random = c("betaRI", "betaSpline1", "betaSpline2", "A"),
             maxiter = 500, silent = 1, alpha = 0.1)
)
