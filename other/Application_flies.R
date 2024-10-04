load("./data/fruitflies.RData")

mllk_flies_sim = function(theta.star, X, N=2, trackInd){
  Gamma = LaMa::tpm(theta.star[1:(N*(N-1))]); p = N*(N-1)
  delta = LaMa::stationary(Gamma)
  mu = exp(theta.star[p + 1:N])
  phi = exp(theta.star[p + N + 1:N])
  
  allprobs = matrix(1, nrow = nrow(X), ncol = N)
  ind = which(!is.na(X$activity))
  for (j in 1:N){
    allprobs[ind, j] = dnbinom(X$activity[ind], mu=mu[j], size=1/phi[j])
  }
  -LaMa::forward(delta, Gamma, allprobs, trackInd)
}

par(mfrow = c(1,1))
hist(data$activity)
N=2
theta.star = c(rep(-2,2), log(c(4, 55, 10, 0.5)))

trackInd = LaMa::calc_trackInd(as.vector(data$ID))
mod_sim_flies = nlm(mllk_flies_sim, theta.star, X = data, N = 2, trackInd = trackInd,
                    print.level=2, iterlim = 1000)
theta.star = mod_sim_flies$estimate
Gamma = LaMa::tpm(theta.star[1:(N*(N-1))]); p = N*(N-1)
delta = LaMa::stationary(Gamma)
mu = exp(theta.star[p + 1:N])
phi = exp(theta.star[p + N + 1:N])

color = c("orange", "deepskyblue")
hist(data$activity, prob = TRUE, breaks = 100, bor = "white")
for(j in 1:N){
  curve(delta[j]*dnbinom(x, mu=mu[j], size=1/phi[j]), add = TRUE, col = color[j], lwd = 2)
}


# splines

mllk_flies = function(theta.star, X, N=2, Z, lambda, S, trackInd){
  # assigning state-dependent parameters
  mu = exp(theta.star[1:N])
  phi = exp(theta.star[N + 1:N]); p = 2*N
  
  # assigning state parameters, beta0 intercept, bLD and bDD matrix
  nb = ncol(Z)-1
  beta0 = matrix(theta.star[p + 1:(2*N*(N-1))], nrow = 2); p = p + 2*N*(N-1)
  bLD = theta.star[p + 1:((N*(N-1))*nb)]; p = p + (N*(N-1)) * nb
  bDD = theta.star[p + 1:((N*(N-1))*nb)]
  b = c(bLD, bDD)
  beta = array(0, dim = c(N, nb+1, 2))
  beta[,-1,1] = matrix(bLD, nrow = N, ncol = nb, byrow = TRUE) #1st slice LD
  beta[,-1,2] = matrix(bDD, nrow = N, ncol = nb, byrow = TRUE) #2nd slice DD
  beta[,1,1] = beta0[,1]
  beta[,1,2] = beta0[,2]
  
  # only computing Gamma for the 48 distinct time points
  Gamma = array(dim=c(N,N,48,2))
  Gamma[,,,1] = LaMa::tpm_g(Z[,-1], beta[,,1])
  Gamma[,,,2] = LaMa::tpm_g(Z[,-1], beta[,,2])
  
  # cleverly indexing Gamma to get in in the right order for all tracks
  ind = which(X$condition == "LD")
  Gamma_track = array(dim = c(N,N,nrow(X)))
  Gamma_track[,,ind] = Gamma[,,X$tod[ind],1]
  Gamma_track[,,-ind] = Gamma[,,X$tod[-ind],2]

  # computing intial, periodically stationary distribution, same for all flies
  delta = LaMa::stationary_p(Gamma[,,,1], t=X$tod[1])
  
  # computing state-dependent probs
  allprobs = matrix(1, nrow = nrow(X), ncol = N)
  ind = which(!is.na(X$activity))
  for (j in 1:N){
    allprobs[ind, j] = dnbinom(X$activity[ind], mu=mu[j], size=1/phi[j])
  }
  
  # computing the penalty
  pen = b %*% kronecker(diag(lambda), S) %*% b
  
  # forward algorithm + penalty
  -LaMa::forward_g(delta, Gamma_track, allprobs, trackInd) + 0.5 * pen # 0.5 is important!
}

nb = 10 # number of basis functions
knots = seq(0, 24, length.out = nb+2) 
# fixing equidistant knots as mgcv otherwise selects them data driven (which is annoying when predicting later)
tod = sort(unique(data$tod)) / 2 # unique time points
# building desing and penalty matrix with mgcv
gam_pre = mgcv::gam(dummy ~ s(tod, bs = "cp", k = nb+1), 
                       data = data.frame(dummy = 1, tod = tod), 
                       knots = list(tod = knots), fit = FALSE)
Z = gam_pre$X
S = gam_pre$S[[1]]

# initial parameter vector
theta.star = c(log(c(mu, phi)), # state-dependent
               rep(-2,2), # intercept LD
               rep(-3,2), # intercept DD
               rep(0, 2*N*(N-1)*nb)) # spline coefs

trackInd = LaMa::calc_trackInd(as.vector(data$ID))

# iterative procedure

maxiter = 100 # maximum number of iterations
tol = 0.01 # relative tolerance for convergence, sufficient
gradtol = 1e-6 # relative gradient tolerance for nlm (1e-6 is the default)
print.level = 0 # print level for nlm
alpha = 0.9

Lambdas = matrix(NA, maxiter+1, N*(N-1)*2)
Lambdas[1,] = rep(1000,4)
mods = list()

p = 2*N*(N-1)+2*N
# matrix containing all the indices of the random effects
REind = p + matrix(1:(2*nb*N*(N-1)), nrow = N*(N-1)*2, byrow = TRUE) # each row is the index of one RE

# updating algorithm
for(k in 1:maxiter){
  cat("\n\n- Iteration", k, "-")
  if(k == 1) print.level = 2
  else print.level = 0
  
  # model fitting
  t1 = Sys.time()
  mod = nlm(mllk_flies, theta.star, N=2, X=data, Z=Z, lambda = Lambdas[k,], S=S, trackInd=trackInd,
            iterlim = 1000, print.level = print.level, hessian = TRUE, gradtol = gradtol)
  Sys.time()-t1; cat("\nEstimation time:", Sys.time()-t1); cat("\nIterations:", mod$iterations)
  
  theta.star = mod$estimate # saves theta.star as starting value for next iteration
  mods[[k]] = mod # saving model
  
  J_inv = MASS::ginv(mod$hessian ) # inverse penalized hessian
  
  # updating all penalty strengths
  for(i in 1:(N*(N-1)*2)){
    edoF = sum(diag(diag(rep(1, nrow(S))) - Lambdas[k, i] * J_inv[REind[i,], REind[i,]] %*% S))
    penalty = t(theta.star[REind[i,]]) %*% S %*% theta.star[REind[i,]]
    lambda_new = as.numeric(edoF / penalty)
    Lambdas[k+1, i] = alpha*lambda_new + (1-alpha)*Lambdas[k, i]
  }
  cat("\nSmoothing strengths:", round(Lambdas[k+1, ], 4))
  if(mean(abs(Lambdas[k+1,] - Lambdas[k,]) / Lambdas[k,]) < tol){
    cat("\n\n")
    break
  }
}

# trace plots of penalty strenghts
Lambdas = as.matrix(na.omit(Lambdas))
par(mfrow = c(1,ncol(Lambdas)))
for(i in 1:ncol(Lambdas)){
  plot(Lambdas[,i], type = "l", lwd = 2, main = paste("lambda", i))
}

# assinging estimated parameter:
mod_final = mods[[length(mods)]]
theta.star = mod_final$estimate

# assigning state-dependent parameters
mu = exp(theta.star[1:N])
phi = exp(theta.star[N + 1:N]); p = 2*N

# assigning state parameters, beta0 intercept, bLD and bDD matrix
nb = ncol(Z)-1
beta0 = matrix(theta.star[p + 1:(2*N*(N-1))], nrow = 2); p = p + 2*N*(N-1)
bLD = theta.star[p + 1:((N*(N-1))*nb)]; p=p+(N*(N-1))*nb
bDD = theta.star[p + 1:((N*(N-1))*nb)]
b = c(bLD, bDD)
beta = array(0, dim = c(N, nb+1, 2))
beta[,-1,1] = matrix(bLD, nrow = N, ncol = nb, byrow = TRUE) #1st slice LD
beta[,-1,2] = matrix(bDD, nrow = N, ncol = nb, byrow = TRUE) #2nd slice DD
beta[,1,1] = beta0[,1]
beta[,1,2] = beta0[,2]

Gamma = array(dim=c(N,N,48,2))
Gamma[,,,1] = LaMa::tpm_g(Z[,-1], beta[,,1])
Gamma[,,,2] = LaMa::tpm_g(Z[,-1], beta[,,2])

# par(mfrow = c(1,1))
# plot(Gamma[1,2,,2], type = "l", col = color[1], lwd = 2)
# 
# DeltaLD = LaMa::stationary_p(Gamma[,,,1]) # stationary for LD
# DeltaDD = LaMa::stationary_p(Gamma[,,,2]) # stationary for DD
# plot(tod, DeltaLD[,2], type = "l", col = color[2], lwd = 2,
#      bty = "n", xlab = "time of day", ylab = "Pr(active)", xaxt = "n")
# axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))

# plotting smoothly
n = 500
knots = seq(0, 24, length.out = nb+2)
Delta = matrix(NA, n, 2)
todseq = seq(0, 24, length.out = n)
for(t in 1:length(todseq)){
  time = (todseq[t] + 1:24) %% 24
  Z = mgcv::gam(y ~ s(tod, bs = "cp", k = nb+1), 
            data = data.frame(dummy = 1, tod = time, y = 1), 
            knots = list(tod = knots), fit = FALSE)$X
  
  Gamma_t = LaMa::tpm_g(Z[,-1], beta[,,1]) # for LD because it is the first slice of beta
  Delta[t,] = LaMa::stationary_p(Gamma_t, 1)
}
par(mfrow = c(1,1))
plot(todseq, Delta[,2], type = "l", lwd = 2, col = color[2], ylim = c(0,1), 
     bty = "n", xlab = "time of day", ylab = "Pr(active)", xaxt = "n")
axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))

Z_plot = gam_pre = mgcv::gam(y ~ s(tod, bs = "cp", k = nb+1), 
                    data = data.frame(dummy = 1, tod = seq(0,24,length=300), y = 1), 
                    knots = list(tod = knots), fit = FALSE)$X
Gamma_plot = LaMa::tpm_g(Z_plot[,-1], beta[,,1])

par(mfrow = c(2,2))
for(i in 1:2){for(j in 1:2){
  plot(seq(0,24,length=300), Gamma_plot[i,j,], type = "l", lwd = 2, ylim = c(0,1))
}}


# Adding a random intercept -----------------------------------------------

mllk_fliesRE = function(theta.star, X, N=2, Z, lambda, S, trackInd){
  # assigning state-dependent parameters
  mu = exp(theta.star[1:N])
  phi = exp(theta.star[N + 1:N]); p = 2*N
  
  # assigning state parameters, beta0 intercept, bLD and bDD matrix
  nb = ncol(Z) - 1
  beta0 = matrix(theta.star[p + 1:(2*N*(N-1))], nrow = 2); p = p + 2*N*(N-1)
  bLD = theta.star[p + 1:((N*(N-1))*nb)]; p = p + (N*(N-1))*nb
  bDD = theta.star[p + 1:((N*(N-1))*nb)]
  beta = array(0, dim = c(nb+1, N*(N-1), 2))
  beta[-1,,1] = matrix(bLD, nrow = nb, ncol = N*(N-1)) # 1st slice LD
  beta[-1,,2] = matrix(bDD, nrow = nb, ncol = N*(N-1)) # 2nd slice DD
  beta[1,,] = beta0

  
  # building the t.p.m. (for LD and DD separately)
  ind = which(X$condition == "LD")
  Gamma = array(dim = c(N,N,nrow(X)))
  Gamma[,,ind] = LaMa::tpm_g(Z[ind,-1], t(beta[,,1])) # LD
  Gamma[,,-ind] = LaMa::tpm_g(Z[-ind,-1], t(beta[,,2])) # DD
  
  # computing intial, periodically stationary distribution for all flies
  nAnimals = length(trackInd)
  Delta = matrix(NA, nAnimals, 2)
  for(a in 1:nAnimals) Delta[a,] = LaMa::stationary_p(Gamma[,,trackInd[a] + 0:47], t = 1)
  
  # computing state-dependent probs
  allprobs = matrix(1, nrow = nrow(X), ncol = N)
  ind = which(!is.na(X$activity))
  for (j in 1:N){
    allprobs[ind, j] = dnbinom(X$activity[ind], mu = mu[j], size = 1/phi[j])
  }
  
  # computing the penalty (\sum_{i=1}^8 b_i^t S_i b_i)
  pen = t(as.numeric(beta[1 + 1:nAnimals,, 1])) %*% kronecker(diag(lambda[1:2]), S[[1]]) %*% as.numeric(beta[1 + 1:nAnimals,, 1]) + # RE LD
    t(as.numeric(beta[1 + 1:nAnimals,, 2])) %*% kronecker(diag(lambda[3:4]), S[[1]]) %*% as.numeric(beta[1 + 1:nAnimals,, 2]) + # RE DD
    t(as.numeric(beta[1 + nAnimals + 1:nrow(S[[2]]),, 1])) %*% kronecker(diag(lambda[5:6]), S[[2]]) %*% as.numeric(beta[1 + nAnimals + 1:nrow(S[[2]]),, 1]) + # spline LD
    t(as.numeric(beta[1 + nAnimals + 1:nrow(S[[2]]),, 2])) %*% kronecker(diag(lambda[7:8]), S[[2]]) %*% as.numeric(beta[1 + nAnimals + 1:nrow(S[[2]]),, 2]) # spline DD
    
  # forward algorithm + penalty
  -LaMa::forward_g(Delta, Gamma, allprobs, trackInd) + 0.5 * pen # 0.5 is important!
}

nb = 12 # number of basis functions
knots = seq(0, 24, length.out = nb+2) 
# fixing equidistant knots as mgcv otherwise selects them data driven (which is annoying when predicting later)
# building desing and penalty matrix with mgcv
gam_pre = mgcv::gam(y ~ s(ID, bs = "re") + s(tod, bs = "cp", k = nb+1), 
                    data = data.frame(dummy = 1, tod = data$tod/2, ID = data$ID, y = 1), 
                    knots = list(tod = knots), fit = FALSE)
Z = gam_pre$X
S = gam_pre$S

# initial parameter vector
theta.star = c(log(c(mu, phi)), # state-dependent
               rep(-2,2), # intercept LD
               rep(-2,2), # intercept DD
               rep(0, 2 * N*(N-1) * (ncol(Z)-1))) # RE and spline coefs

trackInd = LaMa::calc_trackInd(as.vector(data$ID))

# iterative procedure

maxiter = 100 # maximum number of iterations
tol = 0.001 # relative tolerance for convergence
alpha = 0.75 # exponential smoothing parameter for penalty strengths

Lambdas = matrix(NA, maxiter+1, N*(N-1)*4)
# Lambdas[1,] = c(rep(5, 4), rep(0, 4))
Lambdas[1,] = c(rep(50, 4), rep(200, 4))
mods = list()

# building a list that contains the indices for each random effect and spline vector
p = 2*N*(N-1)+2*N
nAnimals = length(unique(data$ID))
REind = list()
for(i in 1:4) REind[[i]] = p + (i-1)*nAnimals + 1:nAnimals
for(i in 5:8) REind[[i]] = p + 4 * nAnimals + (i-5)*nb + 1:nb


# updating algorithm
T1 = Sys.time()
for(k in 1:maxiter){
  cat("\n\n- Iteration", k, "-")
  if(k == 1) print.level = 2
  else print.level = 0
  
  # model fitting
  t1 = Sys.time()
  mod = nlm(mllk_fliesRE, theta.star, N=2, X=data, Z=Z, lambda = Lambdas[k,], S=S, trackInd=trackInd,
            iterlim = 1000, print.level = print.level, hessian = TRUE)
  Sys.time()-t1; cat("\nEstimation time:", Sys.time()-t1); cat("\nIterations:", mod$iterations)
  
  theta.star = mod$estimate # saves theta.star as starting value for next iteration
  mods[[k]] = mod # saving model
  J_inv = MASS::ginv(mod$hessian) # inverse penalized hessian
  
  # updating all penalty strengths
  for(i in 1:8){
    if(i %in% 1:4){
      S_i = S[[1]]
    } else{ S_i = S[[2]] }
    edoF = sum(diag(diag(rep(1, nrow(S_i))) - Lambdas[k, i] * J_inv[REind[[i]], REind[[i]]] %*% S_i))
    penalty = t(theta.star[REind[[i]]]) %*% S_i %*% theta.star[REind[[i]]]
    lambda_new = as.numeric(edoF / penalty)
    Lambdas[k+1, i] = alpha*lambda_new + (1-alpha)*Lambdas[k, i]
  }
  Lambdas[k+1, which(Lambdas[k+1,] < 0)] = 0 # constraining lambda to be non-negative
  
  cat("\nSmoothing strengths:", round(Lambdas[k+1, ], 4))
  if(mean(abs(Lambdas[k+1,] - Lambdas[k,]) / Lambdas[k,]) < tol){
    cat("\n\n")
    break
  }
}
cat("\nTotal estimation time:", Sys.time()-T1, "\n")

# trace plots of penalty strenghts
Lambdas = as.matrix(na.omit(Lambdas))
lambda_hat = Lambdas[nrow(Lambdas),]
par(mfrow = c(2,4))
for(i in 1:ncol(Lambdas)){
  plot(Lambdas[,i], type = "l", lwd = 2, main = paste("lambda", i), bty = "n", ylim = c(0,10))
  abline(h = lambda_hat[i], col = "blue")
}

# standard deviations of the random intercepts
# gamma_21_LD, gamma_12_LD, gamma_21_DD, gamma_12_DD
sqrt(1 / lambda_hat[1:4])

# assinging estimated parameter:
mod_final = mods[[length(mods)]]
theta.star = mod_final$estimate

# state-dependent parameters
mu = exp(theta.star[1:N])
phi = exp(theta.star[N + 1:N]); p = 2*N

# assigning state parameters, beta0 intercept, bLD and bDD matrix
nb = ncol(Z)-1
beta0 = matrix(theta.star[p + 1:(2*N*(N-1))], nrow = 2); p = p + 2*N*(N-1)
bLD = theta.star[p + 1:((N*(N-1))*nb)]; p=p+(N*(N-1))*nb
bDD = theta.star[p + 1:((N*(N-1))*nb)]
beta = array(0, dim = c(N, nb+1, 2))
beta[,-1,1] = matrix(bLD, nrow = N, ncol = nb, byrow = TRUE) #1st slice LD
beta[,-1,2] = matrix(bDD, nrow = N, ncol = nb, byrow = TRUE) #2nd slice DD
beta[,1,1] = beta0[,1]
beta[,1,2] = beta0[,2]


# Plotting the stationary distribution(s) ---------------------------------

n = 200 # plotting smoothness
nb = 12
todseq = seq(0, 24, length.out = n)
Delta_mean = array(dim = c(n, 2, 2))
Delta = array(dim = c(n, 2, 2, nAnimals))

# computing all stationary distributions
for(cond in 1:2){
  for(t in 1:n){
    time = (todseq[t] + 1:48 / 2) %% 24
    Z_splines = mgcv::gam(y ~ s(tod, bs = "cp", k = nb+1), 
                          data = data.frame(dummy = 1, tod = time, y = 1), 
                          knots = list(tod = knots), fit = FALSE)$X
    Z_t = cbind(matrix(0, nrow = 48, ncol = nAnimals), Z_splines[,-1])
    Gamma_t = LaMa::tpm_g(Z_t, beta[,,cond])
    Delta_mean[t,,cond] = LaMa::stationary_p(Gamma_t, 1)
  }
}
for(a in 1:nAnimals){
  print(a)
  for(cond in 1:2){
    for(t in 1:n){
      time = (todseq[t] + 1:48 / 2) %% 24
      Z_splines = mgcv::gam(y ~ s(tod, bs = "cp", k = nb+1), 
                            data = data.frame(dummy = 1, tod = time, y = 1), 
                            knots = list(tod = knots), fit = FALSE)$X
      Z_t = cbind(matrix(0, nrow = 48, ncol = nAnimals), Z_splines[,-1])
      Z_t[,a] = 1
      Gamma_t = LaMa::tpm_g(Z_t, beta[,,cond])
      Delta[t,,cond,a] = LaMa::stationary_p(Gamma_t, 1)
    }
  }
}


# plotting
conditions = c("LD", "DD")
state = 2

par(mfrow = c(1,2))
for(cond in 1:2) {
  plot(NA, xlim = c(0,24), ylim = c(0,1), xlab = "time of day", ylab = "Pr(active)", 
       bty = "n", main = conditions[cond], xaxt = "n")
  axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))
  
  for(a in 1:nAnimals){
    lines(seq(0,24,length=n), Delta[,state,cond,a], col = scales::alpha(a, 0.3), lwd = 2)
  }
  lines(seq(0,24,length=n), Delta_mean[,state,cond], col = "black", lwd = 4)
}





# Parametric with sine and cosine, but with random slope ------------------

library(LaMa)

mllk_par = function(theta.star, X, N=2, Z, trackInd){
  # assigning state-dependent parameters
  mu = exp(theta.star[1:N])
  phi = exp(theta.star[N + 1:N]); p = 2*N
  
  # assigning state parameters, beta0 intercept, bLD and bDD matrix
  nb = ncol(Z)
  beta = array(theta.star[p + 1:(2*(nb+1)*N*(N-1))], dim = c(N*(N-1), nb+1, 2))

  # building the t.p.m. (for LD and DD separately)
  Gamma = array(dim = c(N,N,nrow(X)))
  ind = which(X$condition == "LD")
  Gamma[,,ind] = LaMa::tpm_g(Z[ind, ], beta[,,1]) # LD
  Gamma[,,-ind] = LaMa::tpm_g(Z[-ind, ], beta[,,2]) # DD
  
  # computing intial, periodically stationary distribution for all flies
  nAnimals = length(trackInd)
  Delta = matrix(NA, nAnimals, 2)
  for(a in 1:nAnimals) Delta[a,] = LaMa::stationary_p(Gamma[,,trackInd[a] + 0:47], t = 1)
  
  # computing state-dependent probs
  allprobs = matrix(1, nrow = nrow(X), ncol = N)
  ind = which(!is.na(X$activity))
  for (j in 1:N){
    allprobs[ind, j] = dnbinom(X$activity[ind], mu=mu[j], size=1/phi[j])
  }

  # forward algorithm + penalty
  -LaMa::forward_g(Delta, Gamma, allprobs, trackInd)
}

# building sine-cosine design matrix
Z_par = trigBasisExp(data$tod/2, degree = 4)

# initial parameter vector
theta.star = c(log(c(mu, phi)), # state-dependent
               rep(-2,2), # intercept LD
               rep(0, N*(N-1) * ncol(Z_par)), # coefficients LD
               rep(-2,2), # intercept DD
               rep(0, N*(N-1) * ncol(Z_par))) # coefficients DD

mod_par = nlm(mllk_par, theta.star, X = data, Z = Z_par, trackInd = trackInd,
              iterlim = 1000, print.level = 2, hessian = TRUE)

theta.star = mod_par$estimate
N = 2

# assigning state-dependent parameters
mu = exp(theta.star[1:N])
phi = exp(theta.star[N + 1:N]); p = 2*N

# assigning state parameters, beta0 intercept, bLD and bDD matrix
nb = ncol(Z_par)
beta = array(theta.star[p + 1:(2*(nb+1)*N*(N-1))], dim = c(N*(N-1), nb+1, 2))

# building the t.p.m. (for LD and DD separately)
Z_plot = trigBasisExp(1:48/2, degree = 4)

Gamma_plot = array(dim = c(N,N,48,2))
Gamma_plot[,,,1] = tpm_g(Z_plot, beta[,,1]) # LD
Gamma_plot[,,,2] = tpm_g(Z_plot, beta[,,2]) # LD

DeltaLD = stationary_p(Gamma_plot[,,,1])
DeltaDD = stationary_p(Gamma_plot[,,,2])

plot(DeltaLD[,2], type = "l", col = color[2], lwd = 2, ylim = c(0,1), bty = "n")
lines(DeltaDD[,2], col = color[1], lwd = 2)



# Parametric model with random slope --------------------------------------

library(mgcv)


mllk_re = function(theta.star, X, N=2, Z, lambda, trackInd, degree){
  # assigning state-dependent parameters
  mu = exp(theta.star[1:N])
  phi = exp(theta.star[N + 1:N]); p = 2*N
  
  # assigning state process parameters
  nAnimals = length(trackInd)
  beta_fe = array(theta.star[p + 1 : (2*(2*degree+1)*N*(N-1))], dim = c(N*(N-1), (2*degree + 1), 2)); p = p + 2*(2*degree+1)*N*(N-1)
  beta_rintercept = array(theta.star[p + 1:(2*N*(N-1)*nAnimals)], dim = c(N*(N-1), nAnimals, 2)); p = p + 2*N*(N-1)*nAnimals
  beta_rslope = array(theta.star[p + 1:(2*N*(N-1)*nAnimals*degree*2)], dim = c(N*(N-1), nAnimals*degree*2, 2))
  
  beta = array(dim = c(N*(N-1), ncol(beta_fe) + ncol(beta_rintercept) + ncol(beta_rslope), 2))
  beta[,,1] = cbind(beta_fe[,1,1], beta_rintercept[,,1], beta_fe[,-1,1], beta_rslope[,,1])
  beta[,,2] = cbind(beta_fe[,1,2], beta_rintercept[,,2], beta_fe[,-1,2], beta_rslope[,,2])

  # building the t.p.m. (for LD and DD separately)
  Gamma = array(dim = c(N,N,nrow(X)))
  ind = which(X$condition == "LD")
  Gamma[,,ind] = LaMa::tpm_g(Z[ind, -1], beta[,,1]) # LD
  Gamma[,,-ind] = LaMa::tpm_g(Z[-ind, -1], beta[,,2]) # DD
  
  # computing intial, periodically stationary distribution for all flies
  nAnimals = length(trackInd)
  Delta = matrix(NA, nAnimals, 2)
  for(a in 1:nAnimals) Delta[a,] = LaMa::stationary_p(Gamma[,,trackInd[a] + 0:47], t = 1)
  
  # computing state-dependent probs
  allprobs = matrix(1, nrow = nrow(X), ncol = N)
  ind = which(!is.na(X$activity))
  for (j in 1:N){
    allprobs[ind, j] = dnbinom(X$activity[ind], mu=mu[j], size=1/phi[j])
  }
  
  # computing the penalty
  pen = lambda[1] * sum(beta_rintercept[,,1]^2) + lambda[2] * sum(beta_rintercept[,,2]^2) +
    lambda[3] * sum(beta_rslope[,,1]^2) + lambda[4] * sum(beta_rslope[,,2]^2)
  
  # forward algorithm + penalty
  -LaMa::forward_g(Delta, Gamma, allprobs, trackInd) + 0.5 * pen # 0.5 is important!
}

## building the design matrix
degree = 3 # number of sines and cosines

gam_pre = gam(y ~ s(ID, bs = "re"),
              data = data.frame(dummy = 1, ID = data$ID, y = 1),
              fit = FALSE)
Z = gam_pre$X
Z = cbind(Z, trigBasisExp(data$tod/2, degree = degree))
for(i in 1:length(unique(data$ID))) Z = cbind(Z, Z[,i+1] * trigBasisExp(data$tod/2, degree = degree))



# MMLE algorithm ----------------------------------------------------------

# separate tracks
trackInd = calc_trackInd(as.vector(data$ID))
nAnimals = length(trackInd)

# initial parameter vector
theta.star = c(log(c(mu, phi)), # state-dependent
               rep(-2,2), rep(0, degree*2*2), # fixed effects LD
               rep(-2,2), rep(0, degree*2*2), # fixed effects DD
               rep(0, 2*2*nAnimals), # random intercept
               rep(0, 2*2*nAnimals*(degree*2))) # random slope

# hyperparameters
maxiter = 100 # maximum number of iterations
tol = 0.01 # relative tolerance for convergence, sufficient
print.level = 0 # print level for nlm
alpha = 0.9 # exponential smoothing parameter for penalty strengths

Lambdas = matrix(NA, maxiter+1, 4)
Lambdas[1,] = rep(100, 4)
mods = list()

# building list that contains the indices of the RE vectors
REind = list()
n_fe_par = 2*N + 2*(2*degree+1)*N*(N-1)
REind[[1]] = n_fe_par + 1:(N*(N-1)*nAnimals)
REind[[2]] = n_fe_par + N*(N-1)*nAnimals + 1:(N*(N-1)*nAnimals)
REind[[3]] = n_fe_par + 2*N*(N-1)*nAnimals + 1:(N*(N-1)*nAnimals*degree*2)
REind[[4]] = n_fe_par + 2*N*(N-1)*nAnimals + N*(N-1)*nAnimals*degree*2 + 1:(N*(N-1)*nAnimals*degree*2)

# updating algorithm
T1 = Sys.time()
for(k in 1:maxiter){
  cat("\n\n- Iteration", k, "-")
  if(k == 1) print.level = 2
  else print.level = 0
  
  # model fitting
  t1 = Sys.time()
  mod = nlm(mllk_re, theta.star, X = data, Z = Z, lambda = Lambdas[k,], trackInd = trackInd, degree = degree,
            iterlim = 1000, print.level = print.level, hessian = TRUE)
  Sys.time()-t1; cat("\nEstimation time:", Sys.time()-t1); cat("\nIterations:", mod$iterations)
  
  theta.star = mod$estimate # saves theta.star as starting value for next iteration
  mods[[k]] = mod # saving model
  J_inv = MASS::ginv(mod$hessian) # inverse penalized hessian
  
  # updating all penalty strengths
  for(i in 1:4){
    edoF = sum(diag(diag(rep(1, length(REind[[i]]))) - Lambdas[k,i] * J_inv[REind[[i]], REind[[i]]]))
    penalty = sum(theta.star[REind[[i]]]^2)
    lambda_new = as.numeric(edoF / penalty)
    Lambdas[k+1, i] = alpha * lambda_new + (1-alpha) * Lambdas[k, i]
  }
  Lambdas[k+1, which(Lambdas[k+1,] < 0)] = 0 # constraining lambda to be non-negative
  
  cat("\nSmoothing strengths:", round(Lambdas[k+1, ], 4))
  if(mean(abs(Lambdas[k+1,] - Lambdas[k,]) / Lambdas[k,]) < tol){
    cat("\n\n")
    break
  }
}
cat("\nTotal estimation time:", Sys.time()-T1, "\n")

# trace plots of penalty strenghts
Lambdas = as.matrix(na.omit(Lambdas))
lambda_hat = Lambdas[nrow(Lambdas),]
par(mfrow = c(1,4))
for(i in 1:ncol(Lambdas)){
  plot(Lambdas[,i], type = "l", lwd = 2, main = paste("lambda", i), bty = "n", ylim = c(0, 100))
  abline(h = lambda_hat[i], col = "deepskyblue")
}
sqrt(1/lambda_hat) # standard deviations of the random intercepts and slopes

mod_final = mods[[length(mods)]]
theta.star = mod_final$estimate
N = 2

# assigning state-dependent parameters
mu = exp(theta.star[1:N])
phi = exp(theta.star[N + 1:N]); p = 2*N

# assigning state process parameters
nAnimals = length(trackInd)
beta_fe = array(theta.star[p + 1 : (2*(2*degree+1)*N*(N-1))], dim = c(N*(N-1), (2*degree + 1), 2)); p = p + 2*(2*degree+1)*N*(N-1)
beta_rintercept = array(theta.star[p + 1:(2*N*(N-1)*nAnimals)], dim = c(N*(N-1), nAnimals, 2)); p = p + 2*N*(N-1)*nAnimals
beta_rslope = array(theta.star[p + 1:(2*N*(N-1)*nAnimals*degree*2)], dim = c(N*(N-1), nAnimals*degree*2, 2))

beta = array(dim = c(N*(N-1), ncol(beta_fe) + ncol(beta_rintercept) + ncol(beta_rslope), 2))
beta[,,1] = cbind(beta_fe[,1,1], beta_rintercept[,,1], beta_fe[,-1,1], beta_rslope[,,1])
beta[,,2] = cbind(beta_fe[,1,2], beta_rintercept[,,2], beta_fe[,-1,2], beta_rslope[,,2])


# Plotting the stationary distribution(s) ---------------------------------

n = 300 # plotting smoothness
Delta_mean = array(dim = c(n, 2, 2))
Delta = array(dim = c(n, 2, 2, nAnimals))
todseq = seq(0, 24, length = n)

for(cond in 1:2){
  for(t in 1:n){
    time = (todseq[t] + (1:48 / 2)) %% 24
    Z_trig = trigBasisExp(time, degree = degree)
    Z = matrix(0, nrow = 48, ncol = dim(beta)[2]-1)
    Z[, nAnimals + 1:ncol(Z_trig)] = Z_trig
    
    Gamma_t = LaMa::tpm_g(Z, beta[,,cond])
    Delta_mean[t,,cond] = LaMa::stationary_p(Gamma_t, 1)
  }
}
for(a in 1:nAnimals){
  print(a)
  for(cond in 1:2){
    for(t in 1:n){
      time = (todseq[t] + (1:48 / 2)) %% 24
      Z_trig = trigBasisExp(time, degree = degree)
      Z = matrix(0, nrow = 48, ncol = dim(beta)[2]-1)
      Z[, nAnimals + 1:ncol(Z_trig)] = Z_trig
      Z[, a] = 1
      Z[, nAnimals + 2*degree + (a-1)*2*degree + 1:(2*degree)] = Z_trig
      
      Gamma_t = LaMa::tpm_g(Z, beta[,,cond])
      Delta[t,,cond,a] = LaMa::stationary_p(Gamma_t, 1)
    }
  }
}


par(mfrow = c(1,1))
conditions = c("LD", "DD")
cond = 2
plot(NA, xlim = c(0,24), ylim = c(0,1), xlab = "time of day", ylab = "Pr(active)", 
     bty = "n", main = conditions[cond], xaxt = "n")
axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))

for(a in 1:nAnimals) lines(seq(0,24,length=n), Delta[,2,cond,a], col = scales::alpha(a+1, 0.4))
lines(todseq, Delta_mean[,2,cond], type = "l", lwd = 4)
