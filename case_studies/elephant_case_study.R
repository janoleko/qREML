## Ivory coast elephant case study

## packages
library(mvtnorm)
# install.packages("RTMB")
library(RTMB)
# devtools::install_github("janoleko/LaMa") # development version
library(LaMa)
# install.packages("scales") # for muted colors
library(scales)


## data
data = read.csv("./data/elephant_data.csv")
head(data)
nrow(data)


## colors
color = c("orange", "deepskyblue")


# Fitting homogeneous and parametric HMMs ---------------------------------

mllk_hom = function(par){
  getAll(par, dat)
  ## transforming parameters
  mu = exp(logmu)
  sigma = exp(logsigma)
  kappa = exp(logkappa)
  Gamma = tpm(eta)
  delta = stationary(Gamma)
  ## reporting quantities of interest later
  REPORT(mu)
  REPORT(sigma)
  REPORT(kappa)
  ## computing state-dependent densities
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j]) *
      dvm(angle[ind], 0, kappa[j])
  }
  ## forward algorithm
  -forward(delta, Gamma, allprobs, ad = T)
}

N = 2
par = list(logmu = log(c(0.2, 2)),
           logsigma = log(c(0.2, 2)),
           logkappa = log(c(0.2, 1)),
           eta = rep(-2, 2))

dat = list(step = data$step, angle = data$angle, N = 2)

obj_hom = MakeADFun(mllk_hom, par, silent = TRUE)
opt_hom = nlminb(obj_hom$par, obj_hom$fn, obj_hom$gr)

mod_hom = obj_hom$report()

mu = mod_hom$mu
sigma = mod_hom$sigma
kappa = mod_hom$kappa
delta = mod_hom$delta

mod_hom$states = viterbi(mod_hom$delta, mod_hom$Gamma, mod_hom$allprobs)

mod_hom$AIC = 2 * obj_hom$fn() + 2 * length(obj_hom$par)
mod_hom$BIC = 2 * obj_hom$fn() + log(nrow(data)) * length(obj_hom$par)


# pdf("./case_studies/figs/elephant_marginal.pdf", width = 8, height = 4)

par(mfrow = c(1,2))
hist(data$step, breaks = 100, prob = T, bor = "white", xlim = c(0,5), main = "", xlab = "step length", ylab = "density")
curve(delta[1]*dgamma2(x, mu[1], sigma[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta[2]*dgamma2(x, mu[2], sigma[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta[1]*dgamma2(x, mu[1], sigma[1])+delta[2]*dgamma2(x, mu[2], sigma[2]),
        add = T, lwd = 2, lty = 2, n = 500)
legend("topright", legend = c("encamped", "exploratory", "marginal"), col = c(color[1], color[2], "black"), 
       lty = c(1,1,2), bty = "n")

hist(data$angle, breaks = 20, prob = T, bor = "white", main = "", xlab = "turning angle", ylab = "density")
curve(delta[1]*dvm(x, 0, kappa[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta[2]*dvm(x, 0, kappa[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta[1]*dvm(x, 0, kappa[1])+delta[2]*dvm(x, 0, kappa[2]), 
      add = T, lwd = 2, lty = 2, n = 500)

# dev.off()


## periodic variation in state process
mllk_par = function(par){
  getAll(par, dat)
  ## transforming parameters
  mu = exp(logmu)
  sigma = exp(logsigma)
  kappa = exp(logkappa)
  Gamma = tpm_g(Z, beta, ad = T)
  delta = stationary_p(Gamma, t = tod[1], ad = T)
  ## reporting quantities of interest later
  REPORT(mu)
  REPORT(sigma)
  REPORT(kappa)
  ## computing state-dependent densities
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j]) *
      dvm(angle[ind], 0, kappa[j])
  }
  ## forward algorithm
  -forward_g(delta, Gamma[,,tod], allprobs, ad = T)
}

## prepwork: building trigonometric basis design matrix
Z = trigBasisExp((1:12)*2-1, degree = 2)

par = list(logmu = log(c(0.2, 2)),
           logsigma = log(c(0.2, 2)),
           logkappa = log(c(0.2, 1)),
           beta = matrix(c(rep(-2,2), rep(0, 8)), nrow = 2))

dat = list(step = data$step, angle = data$angle, N = 2, tod = data$tod)

obj_par = MakeADFun(mllk_par, par, silent = TRUE)
opt_par = nlminb(obj_par$par, obj_par$fn, obj_par$gr)

mod_par = obj_par$report()
mod_par$Gamma = tpm_g(Z, mod_par$beta)
mod_par$Delta = stationary_p(mod_par$Gamma)

mod_par$AIC = 2 * obj_par$fn() + 2 * length(obj_par$par)
mod_par$BIC = 2 * obj_par$fn() + log(nrow(data)) * length(obj_par$par)


## confidence intervals by sampling from MLE distribution
B = 10000 # number of samples
mod_par$H = obj_par$he() # saving hessian
sdr = sdreport(obj_par)

thetas_par = rmvnorm(B, sdr$par.fixed, sigma = solve(mod_par$H))

n = 200 # plotting resolution
tod_seq = seq(0,24, length = n) # sequence of time of day values
Z_pred = trigBasisExp(tod_seq, degree = 2) # prediction design matrix

Gamma_boot_par = array(dim = c(2, 2, n, B))
for(b in 1:B){
  beta_boot = matrix(thetas_par[b, 7:ncol(thetas_par)], nrow = N*(N-1))
  Gamma_boot_par[,,,b] = tpm_g(Z_pred, beta_boot)
}
GammaCI_par = apply(Gamma_boot_par, c(1,2,3), quantile, probs = c(0.025, 0.975))
Gamma_plot_par = tpm_p(tod = tod_seq, L = 24, mod_par$beta, degree = 2)



# Transition probabilities as smooth functions ----------------------------

## penalized likelihood function
pnll = function(par) {
  getAll(par, dat)
  
  mu = exp(logmu); REPORT(mu)
  sigma = exp(logsigma); REPORT(sigma)
  kappa = exp(logkappa); REPORT(kappa)
  
  Gamma = tpm_g(Z, cbind(beta0, betaspline), ad = T)
  delta = stationary_p(Gamma, t = tod[1], ad = T)
  
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j]) * 
      dvm(angle[ind],0,kappa[j])
  }
  
  -forward_g(delta, Gamma[,,tod], allprobs, ad = T) + 
    penalty(betaspline, S, lambda) # computes 0.5 * lambda t(b) S b
}

modmat = make_matrices(~ s(tod, bs = "cp"), # default number of basis functions
                       data = data.frame(tod = (1:12)*2-1),
                       knots = list(tod = c(0, 24))) # telling mgcv where to wrap the basis
Z = modmat$Z
S = modmat$S


# initial parameter values
N = 2 
par = list(logmu = log(c(0.35, 1.1)),
           logsigma = log(c(0.25, 0.75)),
           logkappa = log(c(0.2, 0.7)),
           beta0 = c(-2,-2),
           betaspline = matrix(0, nrow = N*(N-1), ncol = ncol(Z)-1))

dat = list(step = data$step, angle = data$angle,
           tod = data$tod,
           N = 2, Z = Z, S = S, lambda = rep(1e5, 2))

system.time(
  mod <- qreml(pnll, par, dat, random = "betaspline", saveall = T, silent = 1)
)

## extracting parameters
beta = mod$beta

# penalty strength
round(mod$lambda, 3)

# edf
round(mod$edf[[1]], 2)

# sim reporting -> Monte Carlo sampling reported quantities from the MLE distribution
simreport = sdreportMC(mod$obj_joint, what = c("beta", "mu", "sigma"), n = 5000) # joint uncertainty
# simreport_cond = sdreportMC(mod$obj, what = c("beta", "mu", "sigma")) # conditinoal uncertainty
betaSim = simreport$beta
# betaSim = simreport_cond$beta

# calculating transition probabilities and stationary distribution for plotting
Z_pred = pred_matrix(modmat, newdata = data.frame(tod = tod_seq)) # prediction design matrix
Gamma_plot = tpm_g(Z_pred, beta)

# calculating periodically stationary distribution
Delta_cont = matrix(NA, length(tod_seq), N)
for(t in 1:length(tod_seq)){
  t_seq = (tod_seq[t] + (1:12)*2 - 1) %% 24
  Z_cont = pred_matrix(modmat, newdata = data.frame(tod = t_seq))
  G = tpm_g(Z_cont, beta)
  Delta_cont[t,] = stationary_p(G, t = 1)
}


# computing confidence bands from Monte Carlo samples
# getting the mle out of the final model
mle = mod$par_vec
B = dim(betaSim)[3]; N = 2
Gamma_boot = array(dim = c(2, 2, length(tod_seq), B))
Delta_boot = array(dim = c(length(tod_seq), 2, B))
for(b in 1:B){
  beta = betaSim[,,b]
  Gamma_boot[,,,b] = tpm_g(Z_pred, beta)
}

# periodically stationary distribution
for(t in 1:length(tod_seq)) {
  if(t %% 10 == 0) {
    print(paste0(round(100*t/length(tod_seq)), "%") )
  }
  t_seq = (tod_seq[t] + (1:12)*2 - 1) %% 24
  Z_cont = pred_matrix(modmat, newdata = data.frame(tod = t_seq))
  
  for(b in 1:B){
    beta = betaSim[,,b]
    Delta_boot[t,,b] = stationary_p(tpm_g(Z_cont, beta), t = 1)
  }
}

# pointwise confidence intervals
GammaCI = apply(Gamma_boot, c(1,2,3), quantile, probs = c(0.025, 0.975), na.rm = T)
DeltaCI = apply(Delta_boot, c(1,2), quantile, probs = c(0.025, 0.975))




## visualizing results: comparing parametric and nonparametric

# pdf("./case_studies/figs/elephant_transprobs.pdf", width = 8, height = 4)

par(mfrow = c(1,2))

# parametric fit
plot(tod_seq, Gamma_plot_par[1,2,], type = "l", lwd = 1, bty = "n", main = "parametric",
     xlab = "time of day", ylab ="transition probability", xaxt = "n", ylim = c(0,1))
polygon(c(tod_seq, rev(tod_seq)), c(GammaCI_par[1,1,2,], rev(GammaCI_par[2,1,2,])), 
        col = alpha("black", 0.1), border = F)
lines(tod_seq, Gamma_plot_par[2,1,], lwd = 1, lty = 5)
polygon(c(tod_seq, rev(tod_seq)), c(GammaCI_par[1,2,1,], rev(GammaCI_par[2,2,1,])), 
        col = alpha("black", 0.1), border = F)
legend(x = -0.5, y = 1.02, lty = c(1,5), y.intersp = 1.3,
       legend = c(expression(gamma[12]^(t)), expression(gamma[21]^(t))), bty = "n")
axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))

# non-parametric fit
plot(tod_seq, Gamma_plot[1,2,], type = "l", lwd = 1, bty = "n", main = "non-parametric",
     xlab = "time of day", ylab ="transition probability", xaxt = "n", ylim = c(0,1))
polygon(c(tod_seq, rev(tod_seq)), c(GammaCI[1,1,2,], rev(GammaCI[2,1,2,])), 
        col = alpha("black", 0.1), border = F)
lines(tod_seq, Gamma_plot[2,1,], lwd = 1, lty = 5)
polygon(c(tod_seq, rev(tod_seq)), c(GammaCI[1,2,1,], rev(GammaCI[2,2,1,])), 
        col = alpha("black", 0.1), border = F)
# legend(x = -0.5, y = 1.02, lty = c(1,2), y.intersp = 1.3,
#        legend = c(expression(gamma[12]^(t)), expression(gamma[21]^(t))), bty = "n")
axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))

# dev.off()


## nonparametric periodically stationary distribution

sun_cycle_colors <- c(
  "#0b0d3e",  # 00:00 - Midnight (Night)
  "#0d1046",  # 00:30
  "#0f124e",  # 01:00
  "#121557",  # 01:30
  "#14185f",  # 02:00
  "#161a67",  # 02:30
  "#191d6f",  # 03:00
  "#1b1f78",  # 03:30
  "#1d227f",  # 04:00
  "#1f2587",  # 04:30
  "#22288f",  # 05:00
  "#564d7b",  # 05:30
  "#877060",  # 06:00
  "#e67e22",  # 06:30 - Sunrise (Vibrant Orange)
  "#f28c32",  # 07:00
  "#f89b42",  # 07:30
  "#fba554",  # 08:00 - Morning
  "#fcaf67",  # 08:30
  "#fcc978",  # 09:00
  "#fcd68d",  # 09:30
  "#fce3a2",  # 10:00
  "#fde1b7",  # 10:30
  "#feeacb",  # 11:00
  "#d8eef5",  # 11:30 - Gradual shift from morning to day
  "#b7e9f9",  # 12:00 - Noon (Light Blue)
  "#a2e2f9",  # 12:30
  "#8cdcf9",  # 13:00
  "#75d5fa",  # 13:30
  "#5ecefa",  # 14:00
  "#47c7fa",  # 14:30
  "#30c0fb",  # 15:00
  "#2ca9e2",  # 15:30
  "#2593cb",  # 16:00 - Afternoon (Deeper Blue)
  "#1f7db4",  # 16:30
  "#19679d",  # 17:00
  "#7f6d63",  # 17:30
  "#e57328",  # 18:00 - Sunset (Vibrant Orange)
  "#cc6925",  # 18:30
  "#b36022",  # 19:00
  "#9a571f",  # 19:30
  "#7c4550",  # 20:00 - Nightfall
  "#603a59",  # 20:30
  "#4c2f55",  # 21:00
  "#38264b",  # 21:30
  "#291f3f",  # 22:00 - Nighttime
  "#1c1732",  # 22:30
  "#140f28",  # 23:00 - Late Night
  "#0e0a1e"   # 23:30
)



# pdf("./case_studies/figs/elephant_stationary.pdf", width = 7, height = 4.5)

par(mfrow = c(1,1))
plot(NA, bty = "n", ylim = c(0,1), xlim = c(0,24), ylab = "Pr(exploratory)", xlab = "time of day", xaxt = "n")
polygon(x = c(0, 6.5, 6.5, 0), y = c(0, 0, 1, 1), col = "gray95", border = "white")
polygon(x = c(18.5, 24, 24, 18.5), y = c(0, 0, 1, 1), col = "gray95", border = "white")
# polygon(x = c(0, 24, 24, 0), y = c(-0.05, -0.05, -0.01, -0.01), col = "black", border = "black")
for(t in 0:47){
  polygon(x = c(t/2, (t+1)/2, (t+1)/2, t/2), y = c(-0.04, -0.04, -0.01, -0.01), col = sun_cycle_colors[t+1], border = sun_cycle_colors[t+1])
}
# for(t in 0:95){
#   polygon(x = c(t/4, (t+1)/4, (t+1)/4, t/4), y = c(-0.04, -0.04, -0.01, -0.01), col = sun_cycle_colors[t+1], border = sun_cycle_colors[t+1])
# }
lines(tod_seq, Delta_cont[,2], lwd = 1.5, col = "black",
     bty = "n", ylim = c(0,1), ylab = "Pr(exploratory)", xlab = "time of day", xaxt = "n")
polygon(c(tod_seq, rev(tod_seq)), c(DeltaCI[1,,2], rev(DeltaCI[2,,2])), 
        col = alpha("black", 0.1), border = F)
axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))

# dev.off()



# Plotting the model sequence ---------------------------------------------

## plotting the model sequence
mods = mod$allmods
length(mods)
plotind = c(1,2,3,5,7,length(mods))

# pdf("./case_studies/figs/elephant_modseq.pdf", width = 8.5, height = 5.5)

par(mfrow = c(2,3), mar = c(5,4,3,1))
for(m in plotind){
  beta = mods[[m]]$beta
  
  Gamma_plot = tpm_g(Z_pred, beta)
  
  plot(tod_seq, Gamma_plot[1,2,], type = "l", lwd = 1, bty = "n", main = paste("Iteration", m),
       xlab = "time of day", ylab ="transition probability", xaxt = "n", ylim = c(0,1))
  lines(tod_seq, Gamma_plot[2,1,], lwd = 1, lty = 5)
  legend(x = -0.5, y = 1.02, lty = c(1,2), y.intersp = 1.3,
         legend = c(expression(gamma[12]^(t)), expression(gamma[21]^(t))), bty = "n")
  axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))

}

#dev.off()




# Full REML and marginal ML -----------------------------------------------

## joint likelihood (has complete normal density for random effects)
jnll = function(par) {
  getAll(par, dat)
  
  mu = exp(logmu); REPORT(mu)
  sigma = exp(logsigma); REPORT(sigma)
  kappa = exp(logkappa); REPORT(kappa)
  
  Gamma = tpm_g(Z, cbind(beta0, betaspline), ad = T)
  delta = stationary_p(Gamma, t = tod[1], ad = T)
  
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j]) * 
      dvm(angle[ind],0,kappa[j])
  }
  lambda = exp(loglambda)
  
  -forward_g(delta, Gamma[,,tod], allprobs, ad = T) - 
    dgmrf(betaspline[1,], 0, lambda[1]*S, log = TRUE) -
    dgmrf(betaspline[2,], 0, lambda[2]*S, log = TRUE)
}

## initial parameters
par = list(logmu = log(c(0.35, 1.1)),
           logsigma = log(c(0.25, 0.75)),
           logkappa = log(c(0.2, 0.7)),
           beta0 = c(-2,-2),
           betaspline = matrix(0, nrow = N*(N-1), ncol = ncol(Z)-1),
           loglambda = rep(log(10),2))

dat = list(step = data$step, angle = data$angle, 
           tod = data$tod, 
           N = 2, Z = Z, S = as(S[[1]], "sparseMatrix"))

## creating objective function
## full REML
obj1 = MakeADFun(jnll, par, random = names(par)[names(par)!="loglambda"])

# model fitting
system.time(
  opt1 <- nlminb(obj1$par, obj1$fn, obj1$gr, control = list(iter.max = 500))
)

## marginal ML
obj2 = MakeADFun(jnll, par, random = "betaspline")
# model fitting
system.time(
  opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr, control = list(iter.max = 500))
)

mod2 = obj1$report()
beta2 = mod2$beta

Gamma_plot2 = tpm_g(Z_pred, beta2)

par(mfrow = c(1,1))
plot(tod_seq, Gamma_plot2[1,2,], type = "l", lwd = 2, bty = "n",
     xlab = "time of day", ylab = "transition probability", xaxt = "n")
lines(tod_seq, Gamma_plot2[2,1,], lwd = 2, lty = 3)
axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))
