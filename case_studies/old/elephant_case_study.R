
# Packages ----------------------------------------------------------------

library(dplyr)
library(LaMa)


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
    allprobs[ind,j] = dgamma(X$step[ind],shape=mu[j]^2/sigma[j]^2,scale=sigma[j]^2/mu[j])*
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

mllk_mgcv = function(theta.star, X, N=2, Z, lambda, S){
  mu = exp(theta.star[1:N])
  sigma = exp(theta.star[N + 1:N])
  kappa = exp(theta.star[2*N + 1:N]); p = 3*N
  nb = ncol(Z)
  beta = matrix(theta.star[p + 1:(nb*N*(N-1))], nrow = nb, ncol = N*(N-1))
  
  Gamma = LaMa::tpm_g(Z[,-1], t(beta))
  delta = LaMa::stationary_p(Gamma, t = X$tod[1])
  
  allprobs = matrix(1, nrow = nrow(X), ncol = N)
  ind = which(!is.na(X$step) & !is.na(X$angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma(X$step[ind],shape=mu[j]^2/sigma[j]^2,scale=sigma[j]^2/mu[j])*
      CircStats::dvm(X$angle[ind], mu = 0, kappa = kappa[j])
  }
  
  pen = t(as.numeric(beta[-1,])) %*% kronecker(diag(lambda), S) %*% as.numeric(beta[-1,])
  
  -LaMa::forward_g(delta, Gamma[,,X$tod], allprobs) + 0.5 * pen
}

nb = 10
knots = seq(0, 24, length.out = nb+1)
# todseq = seq(0, 24, length = 200)
gam_prefit = mgcv::gam(y ~ s(tod, bs = "cp", k = nb), 
                       data = data.frame(dummy = 1, tod = (1:12)*2-1, y = 1), 
                       knots = list(tod = knots), fit = FALSE)
Z = gam_prefit$X
S = gam_prefit$S[[1]]

# hyperparameters for outer maximization
maxiter = 50 # maximum number of iterations
tol = 0.01 # relative tolerance for convergence
gradtol = 1e-6 # relative gradient tolerance for nlm
alpha = 0.99

Lambdas = matrix(NA, maxiter, 2)
Lambdas[1,] = c(1e6, 1e6)
mods = list()

# defining the indices of the spline coefficients
REind = matrix(3*N + 1:(ncol(Z)*N*(N-1)), nrow = 2, byrow = TRUE)[,-1]

N = 2
# initial parameter values
theta.star = c(log(c(0.35, 1.1, 0.25, 0.75, 0.2, 0.7)), # state-dependent process
               -2, rep(0, ncol(Z)-1), -2, rep(0, ncol(Z)-1)) # state process

T1 = Sys.time()
for(k in 1:maxiter){
  cat("\n\n- Iteration", k, "-")
  if(k == 1){
    print.level = 2
  } else print.level = 0
  
  t1 = Sys.time()
  mods[[k]] = nlm(mllk_mgcv, theta.star, X = data, N = 2, Z = Z, lambda = Lambdas[k,], S = S,
                  iterlim = 1000, print.level = print.level, gradtol = gradtol, hessian = TRUE, stepmax = 1000)
  Sys.time()-t1; cat("\nEstimation time:", Sys.time()-t1); cat("\nIterations:", mods[[k]]$iterations)
  
  theta.star = mods[[k]]$estimate # saves theta.star as starting value for next iteration
  J_inv = MASS::ginv(mods[[k]]$hessian ) # inverse penalized hessian
  
  # updating all penalty strengths state-dependent process
  for(i in 1:(N*(N-1))){
    edoF = sum(diag(diag(nrow(S)) - Lambdas[k, i] * J_inv[REind[i,], REind[i,]] %*% S))
    penalty = t(theta.star[REind[i,]]) %*% S %*% theta.star[REind[i,]]
    lambda_new = as.numeric(edoF / penalty) # no -1 because D_cyclic has full rank
    Lambdas[k+1, i] = alpha * lambda_new + (1-alpha) * Lambdas[k, i] # exponential smoothing
  }
  
  cat("\nSmoothing strengths:", round(Lambdas[k+1,], 4))
  if(mean(abs(Lambdas[k+1,] - Lambdas[k,]) / Lambdas[k,]) < tol){
    cat("\n\n- Final model fit")
    t1 = Sys.time()
    mods[[k+1]] = nlm(mllk_mgcv, theta.star, X = data, N = 2, Z = Z, lambda = Lambdas[k+1,], S = S,
                      iterlim = 1000, print.level = print.level, hessian = TRUE, stepmax = 1000)
    Sys.time()-t1; cat("\nEstimation time:", Sys.time()-t1); cat("\nIterations:", mods[[k+1]]$iterations)
    
    cat("\n\n")
    break
  }
}
cat("\nTotal estimation time:", Sys.time()-T1, "\n")

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

# extracting parameters
theta.star = mods[[length(mods)]]$estimate


## AIC and BIC

AIC = 2*mods[[length(mods)]]$minimum + 2 * (4*N + edoF1 + edoF2)
BIC = 2*mods[[length(mods)]]$minimum + log(nrow(data)) * (4*N + edoF1 + edoF2)

# AIC BIC table
AICBIC = cbind(c(AIC_hom, BIC_hom), c(AIC_par, BIC_par), c(AIC, BIC))
colnames(AICBIC) = c("homogeneous", "parametric", "non-parametric")
rownames(AICBIC) = c("AIC", "BIC")

round(AICBIC, 1)


# extracting parameters
N = 2
mu = exp(theta.star[1:N])
sigma = exp(theta.star[N + 1:N])
kappa = exp(theta.star[2*N + 1:N]); p = 3*N
beta = matrix(theta.star[p + 1:(ncol(Z)*N*(N-1))], nrow = ncol(Z), ncol = N*(N-1))
Gamma = LaMa::tpm_g(Z[,-1], t(beta))
delta = LaMa::stationary_p(Gamma, t = data$tod[1])
allprobs = matrix(1, nrow = nrow(data), ncol = N)
ind = which(!is.na(data$step) & !is.na(data$angle))
for(j in 1:N){
  allprobs[ind,j] = dgamma(data$step[ind],shape=mu[j]^2/sigma[j]^2,scale=sigma[j]^2/mu[j])*
    CircStats::dvm(data$angle[ind], mu = 0, kappa = kappa[j])
}

# calculating transition probabilities and stationary distribution for plotting
tod_seq = seq(0, 24, length = 300)
Z_plot = mgcv::gam(y ~ s(tod, bs = "cp", k = nb), 
                   data = data.frame(dummy = 1, tod = tod_seq, y = 1), 
                   knots = list(tod = knots), fit = FALSE)$X

Gamma_plot = LaMa::tpm_g(Z_plot[,-1], t(beta))
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
  Z_cont = mgcv::gam(y ~ s(tod, bs = "cp", k = nb), 
                     data = data.frame(dummy = 1, tod = t_seq, y = 1), 
                     knots = list(tod = knots), fit = FALSE)$X
  G = LaMa::tpm_g(Z_cont[,-1], t(beta))
  Delta_cont[t,] = LaMa::stationary_p(G, t = 1)
}

# computing confidence bands
B = 1000; N = 2
set.seed(123)
theta_boot = mvtnorm::rmvnorm(10*B, theta.star, MASS::ginv(mods[[length(mods)]]$hessian))
# Delta_boot = array(dim = c(300, 2, B))
Gamma_boot = array(dim = c(2, 2, 300, 10*B))
for(b in 1:(10*B)){
  p = 3*N
  beta = matrix(theta_boot[b, p + 1:(ncol(Z)*N*(N-1))], nrow = ncol(Z), ncol = N*(N-1))
  Gamma_boot[,,,b] = LaMa::tpm_g(Z_plot[,-1], t(beta))
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
# saveRDS(Delta_boot, "./case_studies/objects/delta_boot.rds")
delta_boot = readRDS("./case_studies/objects/delta_boot.rds")
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

itershow = c(1,2,3,4,5,12)
# tod_seq = seq(0,24, length = 300)

# pdf("./case_studies/figs/elephant_model_sequence.pdf", width = 7.5, height = 5)

par(mfrow = c(2,3))
for(m in itershow){
  theta.star = mods[[m]]$estimate
  N = 2; p = 3*N
  beta = matrix(theta.star[p + 1:(ncol(Z)*N*(N-1))], nrow = ncol(Z), ncol = N*(N-1))
  
  Gamma_plot = LaMa::tpm_g(Z_plot[,-1], t(beta))
  
  plot(tod_seq, Gamma_plot[1,2,], type = "l", lwd = 1, main = paste("Iteration", m),
       bty = "n", ylim = c(0,1), xlab = "time of day", ylab = "Transition probability", xaxt = "n")
  lines(tod_seq, Gamma_plot[2,1,], lty = 1)
  axis(1, at = seq(0,24,by=4), labels = seq(0,24,by=4))
  
  if(m == 1){
    legend(x = -0.5, y = 1.02, lty = c(1,2), y.intersp = 1.3,
           legend = c(expression(gamma[12]^(t)), expression(gamma[21]^(t))), bty = "n")
  }
}

# dev.off()
