
# Packages ----------------------------------------------------------------

library(dplyr)
library(LaMa)
library(RTMB)
library(mgcv)
# library(CircStats)


# Reading in the data -----------------------------------------------------

data = read.csv("./data/elephant_ivory_coast.csv")
data$date = strptime(data$timestamp, format = "%Y-%m-%d %H:%M:%S")
data$tod = as.numeric(format(data$date, format = "%H"))
data$tod = data$tod/2+1
data$doy = as.numeric(format(data$date, format = "%j"))
data = data %>% dplyr::select(location.long, location.lat, doy, tod, timestamp)
data = moveHMM::prepData(data, coordNames = c("location.long", "location.lat"))


# EDA ---------------------------------------------------------------------

nrow(data)
data$timestamp[1]
data$timestamp[nrow(data)]


# Defining color ----------------------------------------------------------

color = c("orange", "deepskyblue")

# Fitting models ----------------------------------------------------------

## homogeneos HMM

mllk_hom = function(theta.star, N, X){
  mu = exp(theta.star[1:N])
  sigma = exp(theta.star[N + 1:N])
  kappa = exp(theta.star[2*N + 1:N])
  Gamma = LaMa::tpm(theta.star[3*N + 1:(N*(N-1))])
  delta = LaMa::stationary(Gamma)
  allprobs = matrix(1, nrow = nrow(X), ncol = N)
  ind = which(!is.na(X$step) & !is.na(X$angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma(X$step[ind],shape=mu[j]^2/sigma[j]^2,scale=sigma[j]^2/mu[j])*
      CircStats::dvm(X$angle[ind], mu = 0, kappa = kappa[j])
  }
  # returning the negative log likelihood
  -LaMa::forward(delta, Gamma, allprobs)
}

theta.star = c(log(c(0.2, 2, # mu
                     0.2, 2, # sigma
                     0.2, 1)), # kappa
               rep(-2, 2)) # Gamma
N = 2
mod_hom = nlm(mllk_hom, theta.star, N = 2, X = data,
              iterlim = 1000, print.level = 2)
theta.star = mod_hom$estimate 
mu0 = exp(theta.star[1:N])
sigma0 = exp(theta.star[N + 1:N])
kappa0 = exp(theta.star[2*N + 1:N])
Gamma = LaMa::tpm(theta.star[3*N + 1:(N*(N-1))])
delta = LaMa::stationary(Gamma)
beta0 = mod_hom$estimate[3*N + 1:(N*(N-1))]

AIC_hom = 2*mod_hom$minimum + 2*length(mod_hom$estimate)
BIC_hom = 2*mod_hom$minimum + log(nrow(data))*length(mod_hom$estimate)

# pdf("./case_studies/figs/elephant_marginal.pdf", width = 8, height = 4)

par(mfrow = c(1,2))
hist(data$step, breaks = 100, prob = T, bor = "white", xlim = c(0,5), main = "", xlab = "step length", ylab = "density")
curve(delta[1]*dgamma(x, shape = mu0[1]^2/sigma0[1]^2, scale = sigma0[1]^2/mu0[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta[2]*dgamma(x, shape = mu0[2]^2/sigma0[2]^2, scale = sigma0[2]^2/mu0[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta[1]*dgamma(x, shape = mu0[1]^2/sigma0[1]^2, scale = sigma0[1]^2/mu0[1])+
        delta[2]*dgamma(x, shape = mu0[2]^2/sigma0[2]^2, scale = sigma0[2]^2/mu0[2]), add = T, lwd = 2, lty = 2, n = 500)
legend("topright", legend = c("encamped", "exploratory", "marginal"), col = c(color[1], color[2], "black"), 
       lty = c(1,1,2), bty = "n")

hist(data$angle, breaks = 20, prob = T, bor = "white", main = "", xlab = "turning angle", ylab = "density")
curve(delta[1]*CircStats::dvm(x, 0, kappa0[1]), add = T, lwd = 2, col = color[1], n = 500)
curve(delta[2]*CircStats::dvm(x, 0, kappa0[2]), add = T, lwd = 2, col = color[2], n = 500)
curve(delta[1]*CircStats::dvm(x, 0, kappa0[1])+
        delta[2]*CircStats::dvm(x, 0, kappa0[2]), add = T, lwd = 2, lty = 2, n = 500)

# dev.off()


# Parametric fit ----------------------------------------------------------

mllk_par = function(theta.star, X, N=2, K=2){
  beta = matrix(theta.star[1:((N*(N-1))*(1+2*K))], nrow = N*(N-1), ncol = 1+2*K)
  Gamma = LaMa::tpm_p(tod = 1:12*2-1, L = 24, beta = beta, degree = K)
  delta = LaMa::stationary_p(Gamma, t = X$tod[1])
  mu = exp(theta.star[(N*(N-1))*(1+2*K)+1:N])
  sigma = exp(theta.star[(N*(N-1))*(1+2*K)+N+1:N])
  kappa = exp(theta.star[(N*(N-1))*(1+2*K)+2*N+1:N])
  
  allprobs = matrix(1, nrow = nrow(X), ncol = N)
  ind = which(!is.na(X$step) & !is.na(X$angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(X$step[ind], mu[j], sigma[j])*
      CircStats::dvm(X$angle[ind], mu = 0, kappa = kappa[j])
  }
  
  -LaMa::forward_g(delta, Gamma[,,X$tod[-1]], allprobs)
}

theta.star_par = c(rep(-2,2), rep(0,8), 
                   log(mu0), log(sigma0), log(kappa0))
N = 2; K = 2

mod_par = nlm(mllk_par, theta.star_par, X = data, iterlim = 1000, print.level = 2, hessian = T)

## AIC and BIC

AIC_par = 2*mod_par$minimum + 2 * length(mod_par$estimate)
BIC_par = 2*mod_par$minimum + log(nrow(data)) * length(mod_par$estimate)

beta2 = matrix(mod_par$estimate[1:((N*(N-1))*(1+2*K))], nrow = N*(N-1), ncol = 1+2*K)
Gamma2 = tpm_p(tod = 1:12*2, L = 24, beta = beta2, degree = K)
Delta2 = stationary_p(Gamma2)
mu2 = exp(mod_par$estimate[(N*(N-1))*(1+2*K)+1:N])
sigma2 = exp(mod_par$estimate[(N*(N-1))*(1+2*K)+N+1:N])
kappa2 = exp(mod_par$estimate[(N*(N-1))*(1+2*K)+2*N+1:N])

B = 10000
thetas_par = mvtnorm::rmvnorm(B, mod_par$estimate, sigma = solve(mod_par$hessian))

n = 300
tod_seq = seq(0,24,length=n)

Gamma_boot_par = array(dim = c(2, 2, n, B))
for(b in 1:B){
  beta_par = matrix(thetas_par[b,1:((N*(N-1))*(1+2*K))], nrow = N*(N-1), ncol = 1+2*K)
  Gamma_boot_par[,,,b] = LaMa::tpm_p(tod = tod_seq, L = 24, beta_par, degree = K)
}
GammaCI_par = apply(Gamma_boot_par, c(1,2,3), quantile, probs = c(0.025, 0.975))
Gamma_plot_par = LaMa::tpm_p(tod = tod_seq, L = 24, beta2, degree = K)




# Modelling the transition probabilities as smooth functions --------------

mllk_mgcv = function(par) {
  getAll(par, dat)
  mu = exp(logmu)
  sigma = exp(logsigma)
  kappa = exp(logkappa)
  nb = ncol(Z)
  beta = matrix(betavec, nrow = nb, ncol = N*(N-1))
  Gamma = tpm_g(Z, t(beta), ad = T)
  delta = stationary_p(Gamma, t = tod[1], ad = T)
  
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j]) * dvm(angle[ind],0,kappa[j])
  }
  
  pen = 0
  for(i in 1:(N*(N-1))) pen = pen + lambda[i] * (t(beta[-1,i]) %*% S %*% beta[-1,i])
  
  l = forward_g(delta, Gamma[,,tod], allprobs, ad = T)
  
  REPORT(Gamma)
  REPORT(mu)
  REPORT(sigma)
  REPORT(kappa)
  
  -l + 0.5 * pen
}

nb = 10
knots = seq(0, 24, length.out = nb+1)
# todseq = seq(0, 24, length = 200)
gam_prefit = gam(y ~ s(tod, bs = "cp", k = nb), 
                       data = data.frame(dummy = 1, tod = (1:12)*2-1, y = 1), 
                       knots = list(tod = knots), fit = FALSE)
Z = gam_prefit$X
S = gam_prefit$S[[1]]
m = nrow(S) - Matrix::rankMatrix(S)

# hyperparameters for outer maximization
maxiter = 50 # maximum number of iterations
tol = 0.01 # relative tolerance for convergence
gradtol = 1e-6 # relative gradient tolerance for nlm
alpha_sm = 1

Lambdas = matrix(NA, maxiter, 2)
Lambdas[1,] = c(1e6, 1e6)
mods = list()

N = 2
# defining the indices of the spline coefficients
REind = matrix(3*N + 1:(ncol(Z)*N*(N-1)), nrow = 2, byrow = TRUE)[,-1]

# initial parameter values
par = list(logmu = log(c(0.35, 1.1)),
           logsigma = log(c(0.25, 0.75)),
           logkappa = log(c(0.2, 0.7)),
           betavec = c(-2, rep(0, ncol(Z)-1), -2, rep(0, ncol(Z)-1)))
dat = list(step = data$step, angle = data$angle, tod = data$tod, N = 2, Z = Z, S = S, lambda = Lambdas[1,])

# updating algorithm
T1 = Sys.time()
for(k in 1:maxiter){
  cat("\n\n- Iteration", k, "-")
  t1 = Sys.time()
  obj = MakeADFun(mllk_mgcv, par, silent = TRUE) 
  opt = nlminb(obj$par, obj$fn, obj$gr, control = list(rel.tol = 1e-15))
  mod = obj$report()
  coefs = t(mod$beta)
  J = obj$he()
  J_inv = MASS::ginv(J)
  mod$I = J_inv
  mods[[k]] = mod
  cat("\nEstimation time:", Sys.time()-t1)
  for(i in 1:(N*(N-1))){ # updating all penalty strengths
    edoF = nrow(S) - sum(diag(Lambdas[k,i] * J_inv[REind[i,], REind[i,]] %*% S))
    penalty = t(as.numeric(coefs[-1,i])) %*% S %*% as.numeric(coefs[-1,i])
    lambda_new = as.numeric((edoF - m) / penalty) # m is correction if S does not have full rank
    Lambdas[k+1, i] = alpha_sm * lambda_new + (1-alpha_sm) * Lambdas[k, i]
  }
  dat$lambda = Lambdas[k+1,]
  par$betavec = as.numeric(mod$beta)
  cat("\nSmoothing strengths:", round(Lambdas[k+1, ], 4))
  if(max(abs(Lambdas[k+1,] - Lambdas[k,]) / Lambdas[k,]) < tol){
    cat("\n\n")
    break
  }
}
Sys.time()-T1

Lambdas = as.matrix(na.omit(Lambdas))
par(mfrow = c(1,2))
for(i in 1:ncol(Lambdas)){
  plot(Lambdas[,i], type = "l", lwd = 2, main = paste("lambda", i), 
       ylim = c(0,10), bty = "n", xlab = "iteration", ylab = "penalty strength")
}
lambda_hat = Lambdas[nrow(Lambdas),]
round(lambda_hat, 3)

## effective degrees of freedom
edoF1 = sum(diag(diag(rep(1, nrow(S))) - lambda_hat[1] * J_inv[REind[1,], REind[1,]] %*% S))
edoF2 = sum(diag(diag(rep(1, nrow(S))) - lambda_hat[2] * J_inv[REind[2,], REind[2,]] %*% S))


## AIC and BIC

nll = obj$fn()
AIC = 2*nll + 2 * (4*N + edoF1 + edoF2)
BIC = 2*nll + log(nrow(data)) * (4*N + edoF1 + edoF2)

# AIC BIC table
AICBIC = cbind(c(AIC_hom, BIC_hom), c(AIC_par, BIC_par), c(AIC, BIC))
colnames(AICBIC) = c("homogeneous", "parametric", "non-parametric")
rownames(AICBIC) = c("AIC", "BIC")

round(AICBIC, 1)


# extracting parameters
mod_final = mods[[length(mods)]]
# mod_final$hessian = obj$he()
delta = mod_final$delta
Gamma = mod_final$Gamma
beta = mod_final$beta
allprobs = mod_final$allprobs

# calculating transition probabilities and stationary distribution for plotting
tod_seq = seq(0, 24, length = 300)
Z_plot = gam(y ~ s(tod, bs = "cp", k = nb), 
                   data = data.frame(dummy = 1, tod = tod_seq, y = 1), 
                   knots = list(tod = knots), fit = FALSE)$X

Gamma_plot = tpm_g(Z_plot, beta)
par(mfrow = c(N,N))
for(i in 1:N){
  for(j in 1:N){
    plot(tod_seq, Gamma_plot[i,j,], type = "l", lwd = 2, bty = "n", ylim = c(0,1))
  }
}

tod_seq = seq(0, 24, length = 300)
Delta_cont = matrix(NA, length(tod_seq), 2)
for(t in 1:length(tod_seq)){
  t_seq = (tod_seq[t] + (1:12)*2 - 1) %% 24
  Z_cont = gam(dummy ~ s(tod, bs = "cp", k = nb), 
                     data = data.frame(dummy = 1, tod = t_seq), 
                     knots = list(tod = knots), fit = FALSE)$X
  G = tpm_g(Z_cont, beta)
  Delta_cont[t,] = stationary_p(G, t = 1)
}

# computing confidence bands
# getting the mle out of the final model
sdrl = as.list(sdreport(obj), "Estimate")
mle = unlist(sdrl)

B = 1000; N = 2
set.seed(123)
theta_boot = mvtnorm::rmvnorm(10*B, mle, mod_final$I)
Delta_boot = array(dim = c(300, 2, B))
Gamma_boot = array(dim = c(2, 2, 300, 10*B))
for(b in 1:(10*B)){
  p = 3*N
  betavec = theta_boot[b, p + 1:(ncol(Z_plot)*N*(N-1))]
  beta = matrix(betavec, nrow = ncol(Z_plot), ncol = N*(N-1))
  Gamma_boot[,,,b] = tpm_g(Z_plot, t(beta))
}
# for(b in 1:B){
#   if(b %% 10 == 0) print(b)
#   for(t in 1:length(tod_seq)){
#     p = 3*N
#     beta = matrix(theta_boot[b, p + 1:(ncol(Z)*N*(N-1))], nrow = ncol(Z), ncol = N*(N-1))
# 
#     t_seq = (tod_seq[t] + (1:12)*2 - 1) %% 24
#     Z_cont = mgcv::gam(y ~ s(tod, bs = "cp", k = nb),
#                        data = data.frame(dummy = 1, tod = t_seq, y = 1),
#                        knots = list(tod = knots), fit = FALSE)$X
#     G = LaMa::tpm_g(Z_cont[,-1], t(beta))
#     Delta_boot[t,,b] = LaMa::stationary_p(G, t = 1)
#   }
# }
# saveRDS(Delta_boot, "./case_studies/objects/delta_boot2.rds")
delta_boot = readRDS("./case_studies/objects/delta_boot2.rds")
GammaCI = apply(Gamma_boot, c(1,2,3), quantile, probs = c(0.025, 0.975))
DeltaCI = apply(delta_boot, c(1,2), quantile, probs = c(0.025, 0.975))



# Plots -------------------------------------------------------------------

# pdf("./case_studies/figs/elephant_transprobs.pdf", width = 8, height = 4)

par(mfrow = c(1,2))

# parametric fit
plot(tod_seq, Gamma_plot_par[1,2,], type = "l", lwd = 1, bty = "n", main = "parametric",
     xlab = "time of day", ylab ="transition probability", xaxt = "n", ylim = c(0,1))
polygon(c(tod_seq, rev(tod_seq)), c(GammaCI_par[1,1,2,], rev(GammaCI_par[2,1,2,])), 
        col = scales::alpha("black", 0.1), border = F)
lines(tod_seq, Gamma_plot_par[2,1,], lwd = 1, lty = 5)
polygon(c(tod_seq, rev(tod_seq)), c(GammaCI_par[1,2,1,], rev(GammaCI_par[2,2,1,])), 
        col = scales::alpha("black", 0.1), border = F)
legend(x = -0.5, y = 1.02, lty = c(1,5), y.intersp = 1.3,
       legend = c(expression(gamma[12]^(t)), expression(gamma[21]^(t))), bty = "n")
axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))

# non-parametric fit
plot(tod_seq, Gamma_plot[1,2,], type = "l", lwd = 1, bty = "n", main = "non-parametric",
     xlab = "time of day", ylab ="transition probability", xaxt = "n")
polygon(c(tod_seq, rev(tod_seq)), c(GammaCI[1,1,2,], rev(GammaCI[2,1,2,])), 
        col = scales::alpha("black", 0.1), border = F)
lines(tod_seq, Gamma_plot[2,1,], lwd = 1, lty = 5)
polygon(c(tod_seq, rev(tod_seq)), c(GammaCI[1,2,1,], rev(GammaCI[2,2,1,])), 
        col = scales::alpha("black", 0.1), border = F)
# legend(x = -0.5, y = 1.02, lty = c(1,2), y.intersp = 1.3,
#        legend = c(expression(gamma[12]^(t)), expression(gamma[21]^(t))), bty = "n")
axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))

# dev.off()

# non-parametric delta
# pdf("./case_studies/figs/elephant_stationary.pdf", width = 5, height = 4)
par(mfrow = c(1,1))
plot(tod_seq, Delta_cont[,2], type = "l", lwd = 1, col = "deepskyblue",
     bty = "n", ylim = c(0,1), ylab = "Pr(exploratory)", xlab = "time of day", xaxt = "n")
polygon(c(tod_seq, rev(tod_seq)), c(DeltaCI[1,,2], rev(DeltaCI[2,,2])), 
        col = scales::alpha("deepskyblue", 0.2), border = F)
axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))
# dev.off()


# Plotting the model sequence ---------------------------------------------

itershow = c(1,2,3,4,5,10)
# tod_seq = seq(0,24, length = 300)

# pdf("./case_studies/figs/elephant_model_sequence.pdf", width = 7.5, height = 5)

par(mfrow = c(2,3))
for(m in itershow){
  beta = mods[[m]]$beta
  Gamma_plot = tpm_g(Z_plot, beta)
  
  plot(tod_seq, Gamma_plot[1,2,], type = "l", lwd = 1, main = paste("Iteration", m),
       bty = "n", ylim = c(0,1), xlab = "time of day", ylab = "Transition probability", xaxt = "n")
  lines(tod_seq, Gamma_plot[2,1,], lty = 2)
  axis(1, at = seq(0,24,by=4), labels = seq(0,24,by=4))
  
  if(m == 1){
    legend(x = -0.5, y = 1.02, lty = c(1,2), y.intersp = 1.3,
           legend = c(expression(gamma[12]^(t)), expression(gamma[21]^(t))), bty = "n")
  }
}

# dev.off()





# Writing a function for pql ----------------------------------------------

library(RTMB)
library(LaMa)
library(mgcv)

mllk = function(par) {
  getAll(par, dat)
  
  mu = exp(logmu)
  sigma = exp(logsigma)
  kappa = exp(logkappa)
  
  Gamma = tpm_g(Z, cbind(beta0, betaspline), ad = T)
  delta = stationary_p(Gamma, t = tod[1], ad = T)
  
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j]) * dvm(angle[ind],0,kappa[j])
  }
  
  l = forward_g(delta, Gamma[,,tod], allprobs, ad = T)
  
  REPORT(mu)
  REPORT(sigma)
  REPORT(kappa)
  
  -l + penalty(betaspline, S, lambda) # computes 0.5 * lambda t(b) S b
}



nb = 10

modmat = make_matrices(~s(tod, bs = "cp", k = nb), data = data.frame(tod = (1:12)*2-1))
Z = modmat$Z
S = modmat$S

# initial parameter values

N = 2
par = list(logmu = log(c(0.35, 1.1)),
           logsigma = log(c(0.25, 0.75)),
           logkappa = log(c(0.2, 0.7)),
           beta0 = c(-2,2),
           betaspline = matrix(0, nrow = N*(N-1), ncol = ncol(Z)-1))

dat = list(step = data$step, angle = data$angle, 
           tod = data$tod, 
           N = 2, Z = Z, S = S, lambda = rep(10000, 2))

mod = pql(mllk, par, dat, random = "betaspline")


## extracting parameters
beta = mod$beta

# calculating transition probabilities and stationary distribution for plotting
tod_seq = seq(0, 24, length = 300)

Z_pred = pred_matrix(modmat, newdata = data.frame(tod = tod_seq))

Gamma_plot = tpm_g(Z_pred, beta)

par(mfrow = c(N,N), mar = c(5,4.5,4,2)+0.1)
for(i in 1:N){
  for(j in 1:N){
    plot(tod_seq, Gamma_plot[i,j,], type = "l", lwd = 2, bty = "n", 
         ylim = c(0,1), ylab = bquote(gamma[.(i)*.(j)]^t),
         xlab = "time of day", xaxt = "n")
    axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))
  }
}


Delta_cont = matrix(NA, length(tod_seq), 2)

for(t in 1:length(tod_seq)){
  t_seq = (tod_seq[t] + (1:12)*2 - 1) %% 24
  Z_cont = pred_matrix(modmat, newdata = data.frame(tod = t_seq))
  G = tpm_g(Z_cont, beta)
  Delta_cont[t,] = stationary_p(G, t = 1)
}












# With random effects -----------------------------------------------------

mllk2 = function(par) {
  getAll(par, dat)
  mu = exp(logmu)
  sigma = exp(logsigma)
  kappa = exp(logkappa)
  nb = ncol(Z)
  beta = matrix(c(beta0, betavec), nrow = nb, ncol = N*(N-1), byrow = TRUE)
  Gamma = tpm_g(Z, t(beta), ad = T)
  delta = stationary_p(Gamma, t = tod[1], ad = T)
  
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j]) * dvm(angle[ind],0,kappa[j])
  }
  
  l = forward_g(delta, Gamma[,,tod], allprobs, ad = T)
  
  REPORT(Gamma)
  REPORT(mu)
  REPORT(sigma)
  REPORT(kappa)
  
  l = l + dgmrf(t(beta[-1,1]), 0, exp(loglambda[1]) * S, log = TRUE)
  l = l + dgmrf(t(beta[-1,2]), 0, exp(loglambda[2]) * S, log = TRUE)
  -l
}

S_sparse = as(S, "sparseMatrix")

# initial parameter values
par = list(logmu = log(c(0.35, 1.1)),
           logsigma = log(c(0.25, 0.75)),
           logkappa = log(c(0.2, 0.7)),
           beta0 = c(-2,-2),
           loglambda = log(c(1000, 1000)),
           betavec = rep(0, 2*(ncol(Z)-1)))

dat = list(step = data$step, angle = data$angle, 
           tod = data$tod, N = 2, Z = Z, S = S_sparse)

obj2 = MakeADFun(mllk2, par, random = c("betavec"))

t1 = Sys.time()
opt2 = nlminb(obj2$par, obj2$fn, obj2$gr)
Sys.time()-t1

mod2 = obj2$report()

delta = mod2$delta
Gamma = mod2$Gamma
beta = mod2$beta
allprobs = mod2$allprobs

# calculating transition probabilities and stationary distribution for plotting
tod_seq = seq(0, 24, length = 300)
Z_plot = gam(dummy ~ s(tod, bs = "cp", k = nb), 
             data = data.frame(dummy = 1, tod = tod_seq), 
             knots = list(tod = knots), fit = FALSE)$X

Gamma_plot = tpm_g(Z_plot, beta)
par(mfrow = c(N,N))
for(i in 1:N){
  for(j in 1:N){
    plot(tod_seq, Gamma_plot[i,j,], type = "l", lwd = 2, bty = "n", ylim = c(0,1))
  }
}



# Simulate data from the model --------------------------------------------
library(mgcv)
mod = mod_final
mod$beta
mod$delta
mod$mu

nb = 10
knots = seq(0, 24, length.out = nb+1)
# todseq = seq(0, 24, length = 200)
gam_prefit = gam(dummy ~ s(tod, bs = "cp", k = nb), 
                       data = data.frame(dummy = 1, tod = 1:24), 
                       knots = list(tod = knots), fit = FALSE)
Z = gam_prefit$X

Gamma = tpm_g(Z, mod$beta)

n = 10000
s = rep(NA, n)

tod = rep(1:24, n/20)
tod = tod[(length(tod)-n+1):length(tod)]

elephant = data.frame(tod = tod)

delta = stationary_p(Gamma, t = tod[1])

s[1] = sample(1:2, 1, prob = delta)
for(t in 2:n) {
  s[t] = sample(1:2, 1, prob = Gamma[s[t-1],,tod[t]])
}

elephant$state = s  

step = rgamma(n, shape = mod$mu[s]^2/mod$sigma[s]^2, 
              scale = mod$sigma[s]^2/mod$mu[s])
angle = rep(NA, n)
for(t in 1:n) angle[t] = rvm(1, pi, mod$kappa[s[t]]) - pi

elephant$step = step
elephant$angle = angle

elephant = elephant[,c(1,3,4,2)]
elephant$angle[c(1,n)] = NA
elephant$step[n] = NA

head(elephant)

saveRDS(elephant, "~/Desktop/elephant_sim.rds")
