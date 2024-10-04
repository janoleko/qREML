load("~/Downloads/bull_sharks_workspace_for_plots.RData")

sharks = x
IDs = unique(sharks$SharkID)

hist(sharks$log_ODBA, breaks = 100, prob = TRUE)

mllk_sim = function(theta.star, X, N, trackInd){
  mu = theta.star[1:N]
  sigma = exp(theta.star[N + 1:N])
  
  Gamma = LaMa::tpm(theta.star[2*N + 1:(N*(N-1))]);
  delta = LaMa::stationary(Gamma)
  
  allprobs = matrix(1, nrow = nrow(X), ncol = N)
  ind = which(!is.na(X$log_ODBA))
  for (j in 1:N){
    allprobs[ind, j] = dnorm(X$log_ODBA[ind], mu[j], sigma[j])
  }
  -LaMa::forward(delta, Gamma, allprobs, trackInd)
}

library(LaMa)

trackInd = calc_trackInd(sharks$SharkID)

N = 2

mu0 = c(-0.7, -0.2)
sigma0 = c(0.25, 0.25)

theta.star = c(mu0, log(sigma0), rep(-3, N*(N-1)))

mod0 = nlm(mllk_sim, theta.star, X = sharks, N = N, trackInd = trackInd,
           print.level = 2, iterlim = 1000)

theta.star = mod0$estimate
mu = theta.star[1:N]
sigma = exp(theta.star[N + 1:N])

Gamma = LaMa::tpm(theta.star[2*N + 1:(N*(N-1))]);
delta = LaMa::stationary(Gamma)

hist(sharks$log_ODBA, breaks = 100, prob = TRUE)
color = c("orange", "deepskyblue")
for (j in 1:N){
  curve(delta[j]*dnorm(x, mu[j], sigma[j]), add = TRUE, col = color[j], lwd = 2)
}
curve(delta[1]*dnorm(x, mu[1], sigma[1]) + delta[2]*dnorm(x, mu[2], sigma[2]), 
      add = TRUE, lwd = 2, lty = 2)



# spline model

mllk = function(theta.star, X, N, trackInd, Z, lambda, D){
  mu = theta.star[1:N]
  sigma = exp(theta.star[N + 1:N]); p = 2*N
  
  beta0 = theta.star[p + 1:(N*(N-1))]; p = p + N*(N-1)
  b = theta.star[p + 1:((ncol(Z)-1) * N*(N-1))]
  
  Beta = matrix(0, nrow = N*(N-1), ncol = ncol(Z)+1)
  Beta[,1] = beta0
  Beta[,-(1:2)] = matrix(b, nrow = N*(N-1), ncol = ncol(Z)-1, byrow = TRUE)
  
  Gamma = LaMa::tpm_g(Z, Beta)
  # Delta = matrix(NA, nrow = length(unique(X$SharkID)), ncol = N)
  # for(i in 1:length(unique(X$SharkID))){
  #   Delta[i,] = LaMa::stationary_p(Gamma, X$Hour[trackInd[i]]+1)
  # }
  delta = rep(1/N, N)
  
  allprobs = matrix(1, nrow = nrow(X), ncol = N)
  ind = which(!is.na(X$log_ODBA))
  for (j in 1:N){
    allprobs[ind, j] = dnorm(X$log_ODBA[ind], mu[j], sigma[j])
  }
  
  pen = t(b) %*% kronecker(diag(lambda), D) %*% b
  
  -LaMa::forward_g(Delta, Gamma[,,X$Hour+1], allprobs, trackInd) + 0.5 * pen
}

nb = 15 # number of basis functions
knots = 24 * 0:nb / nb # knots
Z = mgcv::cSplineDes(0:23, knots) ## cyclic spline design matrix
L = WH:::build_D_mat(nb-1, 2); D = t(L)%*%L # difference penalty matrix

# iterative procedure

maxiter = 100 # maximum number of iterations
tol = 0.01 # relative tolerance for convergence, sufficient
gradtol = 1e-6 # relative gradient tolerance for nlm (1e-6 is the default)
print.level = 0 # print level for nlm

Lambdas = matrix(NA, maxiter+1, N*(N-1))
Lambdas[1,] = rep(100, N*(N-1))
mods = list()

p = N*(N-1) + 2*N
REind = matrix(1:((nb-1)*N*(N-1)), nrow = N*(N-1), byrow = TRUE) # each row is the index of one RE

theta.star = c(mu, log(sigma), # state-dependent pars
               rep(-3, N*(N-1)), # beta0
               rep(0, (ncol(Z)-1) * N*(N-1))) # spline coefs

# updating algorithm
for(k in 1:maxiter){
  cat("\n\n- Iteration", k, "-")
  if(k == 1) print.level = 2
  else print.level = 0
  
  t1 = Sys.time()
  mod = nlm(mllk, theta.star, X = sharks, N = N, trackInd = trackInd, 
            Z = Z, lambda = Lambdas[k,], D = D, 
            print.level = print.level, iterlim = 1000, gradtol = gradtol, hessian=T)
  Sys.time()-t1
  cat("\nEstimation time:", Sys.time()-t1); cat("\nIterations:", mod$iterations)
  
  theta.star = mod$estimate # saves theta.star as starting value for next iteration
  mods[[k]] = mod # saving model
  J_p = mod$hessian # assigning hessian
  b = theta.star[p + 1:(N*(N-1)*(nb-1))]
  
  # updating all penalty strengths
  for(i in 1:(N*(N-1))){
    edoF = sum(diag(diag(rep(1, nrow(D))) - Lambdas[k, i] * MASS::ginv(J_p)[p+REind[i,], p+REind[i,]] %*% D))
    penalty = t(b[REind[i,]]) %*% D %*% b[REind[i,]]
    Lambdas[k+1, i] = as.numeric((edoF - 1) / penalty)
  }
  cat("\nSmoothing strengths:", round(Lambdas[k+1, ], 4))
  if(mean(abs(Lambdas[k+1,] - Lambdas[k,]) / Lambdas[k,]) < tol){
    cat("\n\n")
    break
  }
}

Lambdas = as.matrix(na.omit(Lambdas))
par(mfrow = c(1,ncol(Lambdas)))
for(i in 1:ncol(Lambdas)){
  plot(Lambdas[,i], type = "l", lwd = 2, main = paste("lambda", i))
}

mod_final = nlm(mllk, theta.star, X = sharks, N = N, trackInd = trackInd, 
          Z = Z, lambda = Lambdas[k,], D = D, 
          print.level = 2, iterlim = 1000, hessian=T)

theta.star = mod$estimate

mu = theta.star[1:N]
sigma = exp(theta.star[N + 1:N]); p = 2*N

beta0 = theta.star[p + 1:(N*(N-1))]; p = p + N*(N-1)
b = theta.star[p + 1:((ncol(Z)-1) * N*(N-1))]

Beta = matrix(0, nrow = N*(N-1), ncol = ncol(Z)+1)
Beta[,1] = beta0
Beta[,-(1:2)] = matrix(b, nrow = N*(N-1), ncol = ncol(Z)-1, byrow = TRUE)

n = 300
tod_seq = seq(0, 24, length = n)
Z_plot = mgcv::cSplineDes(tod_seq, knots) ## cyclic spline design matrix
Gamma = LaMa::tpm_g(Z_plot, Beta)

Delta = stationary_p(tpm_g(Z, Beta))
Delta2 = matrix(NA, nrow = n, ncol = N)
for(i in 1:n){
  Delta2[i,] = LaMa::stationary(Gamma[,,i])
}

par(mfrow = c(1,1))
plot(tod_seq, Delta2[,2], type = "l", col = color[2], ylim = c(0,1), lwd = 2)

plot(0:23, Delta[,2], type = "l", col = color[2], ylim = c(0,1), lwd = 2, 
     bty = "n", xlab = "time of day", ylab = "probability of being active")
points(0:23, Delta[,2], pch = 19, col = color[2])


# decoding states
allprobs = matrix(1, nrow = nrow(sharks), ncol = N)
ind = which(!is.na(sharks$log_ODBA))
for (j in 1:N){
  allprobs[ind, j] = dnorm(sharks$log_ODBA[ind], mu[j], sigma[j])
}
states = viterbi_g(c(0.5,0.5), Gamma[,,sharks$Hour+1], allprobs)

ind = 1:1000
plot(sharks$log_ODBA[ind], type = "h", col = color[states[ind]])





## parametric model
Z2 = trigBasisExp(0:23, degree = 2)

mllk2 = function(theta.star, X, N, trackInd, Z2){
  mu = theta.star[1:N]
  sigma = exp(theta.star[N + 1:N]); p = 2*N
  
  Beta = matrix(theta.star[p + 1:((N*(N-1)) * (ncol(Z2) + 1))], 
                nrow = N*(N-1), ncol = ncol(Z2)+1)
  
  Gamma = LaMa::tpm_g(Z2, Beta)
  Delta = matrix(NA, nrow = length(unique(X$SharkID)), ncol = N)
  for(i in 1:length(unique(X$SharkID))){
    Delta[i,] = LaMa::stationary_p(Gamma, X$Hour[trackInd[i]]+1)
  }
  # delta = rep(1/N, N)
  
  allprobs = matrix(1, nrow = nrow(X), ncol = N)
  ind = which(!is.na(X$log_ODBA))
  for (j in 1:N){
    allprobs[ind, j] = dnorm(X$log_ODBA[ind], mu[j], sigma[j])
  }
  
  -LaMa::forward_g(Delta, Gamma[,,X$Hour+1], allprobs, trackInd)
}

theta.star = c(mu, log(sigma), # state-dependent pars
               rep(-3, N*(N-1)), # beta0
               rep(0, (N*(N-1)) * ncol(Z2))) # spline coefs

mod2 = nlm(mllk2, theta.star, X = sharks, N = N, trackInd = trackInd, Z2 = Z2,
           print.level = 2, iterlim = 1000)

theta.star = mod2$estimate

mu = theta.star[1:N]
sigma = exp(theta.star[N + 1:N]); p = 2*N

Beta = matrix(theta.star[p + 1:((N*(N-1)) * (ncol(Z2) + 1))], 
              nrow = N*(N-1), ncol = ncol(Z2)+1)
Gamma = LaMa::tpm_g(Z2, Beta)

Delta = stationary_p(Gamma)

plot(0:23 + 0.5, Delta[,2], type = "l", col = color[2], ylim = c(0,1), lwd = 2, bty = "n")
points(0:23 + 0.5, Delta[,2], pch = 19, col = color[2])







# New model as in the paper -----------------------------------------------

## GAM effects over time for each individual, 
## non-parametric effects of time of day on the transition probabilities (pooled)

mllk = function(theta.star, X, trackInd, Z_sdd, Z_s, lambda, D_sdd, D_s){
  k = length(trackInd) # number of individuals
  N=2 # fix N at 2 for this application

  # parameters for the GAM in the state-dependent process
  b1 = theta.star[1:(k * ncol(Z_sdd))]; p = k * ncol(Z_sdd) # RE vector 1
  B1 = matrix(b1, nrow = k, ncol = ncol(Z_sdd), byrow = TRUE) # arrange in matrix
  mu2 = theta.star[p + 1:k]; p = p + k # individual-specific deviations from state 1 mean
  Sigma = matrix(exp(theta.star[p + 1:(N*k)]), nrow = k, ncol = N); p = p+N*k 
  
  # parameters for the GAM in the transition probabilities
  beta0 = theta.star[p + 1:(N*(N-1))]; p = p + N*(N-1) # intercepts
  beta1 = theta.star[p + 1:(N*(N-1))]; p = p + N*(N-1) # coefs for AvgTemp
  nb = (ncol(Z_s)-1)/2 # number of basis functions in state process
  b2 = b3 = rep(0, nb*N*(N-1))
  b2[-c(1,nb+1)] = theta.star[p + 1:((nb-1) * N*(N-1))]; p = p + (nb-1)*N*(N-1) # RE vector 2
  b3[-c(1,nb+1)] = theta.star[p + 1:((nb-1) * N*(N-1))]; p = p + (nb-1)*N*(N-1) # RE vector 3
  Beta = matrix(NA, nrow = N*(N-1), ncol = 2*nb+2)
  Beta[,1:2] = cbind(beta0, beta1)
  Beta[,2 + 1:nb] = matrix(b2, nrow = N*(N-1), ncol = nb, byrow = TRUE)
  Beta[,2 + nb + 1:nb] = matrix(b3, nrow = N*(N-1), ncol = nb, byrow = TRUE)
  
  # initial distribution
  delta = c(1, exp(theta.star[p + 1:(N-1)]))
  delta = delta / sum(delta)
  
  # calculating the transition probability matrix
  Gamma = LaMa::tpm_g(Z_s, Beta)

  # calculating the state-dependent probabilities
  smooth = rowSums(Z_sdd * B1[X$IDnum,]) # GAM means for state 1 and each individual (on log scale)
  Mu1 = exp(smooth)
  Mu2 = exp(smooth  + mu2[X$IDnum]) # means for state two by individual-specific deviations
  sigma1 = Sigma[X$IDnum, 1]
  sigma2 = Sigma[X$IDnum, 2]
  
  allprobs = matrix(1, nrow = nrow(X), ncol = N)
  ind = which(!is.na(X$ODBA))
  allprobs[ind,1] = dgamma(X$ODBA[ind], shape = Mu1[ind]^2/sigma1[ind]^2, scale = sigma1[ind]^2/Mu1[ind])
  allprobs[ind,2] = dgamma(X$ODBA[ind], shape = Mu2[ind]^2/sigma2[ind]^2, scale = sigma2[ind]^2/Mu2[ind])
  
  # calculating the penalty term
  pen = t(b1) %*% kronecker(diag(lambda[1:k]), D_sdd) %*% b1 + 
    t(b2) %*% kronecker(diag(lambda[k + 1:(N*(N-1))]), D_s) %*% b2 +
    t(b3) %*% kronecker(diag(lambda[k + N*(N-1) + 1:(N*(N-1))]), D_s) %*% b3
  
  # forward algorithm for separate tracks + penalty term
  -LaMa::forward_g(delta, Gamma, allprobs, trackInd) + 0.5 * pen
}

## design matrices
# state process design matrix
nb1 = 15 # number of basis functions for cyclic splines in state process
knots = 24 * 0:nb1 / nb1 # knots
Z_tod = mgcv::cSplineDes(sharks$Hour, knots) ## cyclic spline design matrix
Z_s = cbind(sharks$AvgTemp, Z_tod, matrix(sharks$AvgTemp, nrow=nrow(Z_tod), ncol=ncol(Z_tod)) * Z_tod) # design matrix with interaction block
L = WH:::build_D_mat(nb1, 2); D_s = t(L)%*%L # difference penalty matrix

# state-dependent process design matrix
nb2 = 15
ord = 4
degree = ord-1
nrknots = nb2 - (ord-2) 
knots = seq(-600, 3000, length = nrknots + 2 * degree)
Z_sdd = spline.des(knots, sharks$IDindex, ord, outer.ok=T)$design
L = WH:::build_D_mat(nb2, 2); D_sdd = t(L)%*%L # difference penalty matrix

trackInd = LaMa::calc_trackInd(sharks$SharkID)
K = length(trackInd)
mu = c(rep(-0.5, nb2 * K), rep(0.3, K))
sigma = rep(0.16, 2*K)
beta0 = c(-2, -2)
beta1 = c(-0.2, 0.2)
b2 = rep(0, (nb1-1) * N*(N-1)) # spline coefs
b3 = rep(0, (nb1-1) * N*(N-1)) # spline coefs
delta = 0.9
theta.star = c(mu, log(sigma), beta0, beta1, b2, b3, delta)

lambda = c(c(10,50,50,30,30,10,10)/2, # Rolands smoothing strengths
           rep(0, 2), rep(0, 2))

t1 = Sys.time()
mod2 = nlm(mllk, theta.star, X = sharks, trackInd = trackInd, Z_sdd = Z_sdd, Z_s = Z_s, 
           lambda = lambda, D_sdd = D_sdd, D_s = D_s,
           print.level = 2, iterlim = 1000)
Sys.time()-t1

theta.star = mod2$estimate
k = length(trackInd); N=2

# parameters for the GAM in the state-dependent process
b1 = theta.star[1:(k * ncol(Z_sdd))]; p = k * ncol(Z_sdd) # RE vector 1
B1 = matrix(b1, nrow = k, ncol = ncol(Z_sdd), byrow = TRUE) # arrange in matrix
mu2 = theta.star[p + 1:k]; p = p + k # individual-specific deviations from state 1 mean
Sigma = matrix(exp(theta.star[p + 1:(N*k)]), nrow = k, ncol = N); p = p+N*k 

# parameters for the GAM in the transition probabilities
beta0 = theta.star[p + 1:(N*(N-1))]; p = p + N*(N-1) # intercepts
beta1 = theta.star[p + 1:(N*(N-1))]; p = p + N*(N-1) # coefs for AvgTemp
nb = (ncol(Z_s)-1)/2 # number of basis functions in state process
b2 = b3 = rep(0, nb*N*(N-1))
b2[-c(1,nb+1)] = theta.star[p + 1:((nb-1) * N*(N-1))]; p = p + (nb-1)*N*(N-1) # RE vector 2
b3[-c(1,nb+1)] = theta.star[p + 1:((nb-1) * N*(N-1))]; p = p + (nb-1)*N*(N-1) # RE vector 3
Beta = matrix(NA, nrow = N*(N-1), ncol = 2*nb+2)
Beta[,1:2] = cbind(beta0, beta1)
Beta[,2 + 1:nb] = matrix(b2, nrow = N*(N-1), ncol = nb, byrow = TRUE)
Beta[,2 + nb + 1:nb] = matrix(b3, nrow = N*(N-1), ncol = nb, byrow = TRUE)

# initial distribution
delta = c(1, exp(theta.star[p + 1:(N-1)]))
delta = delta / sum(delta)


## visualize results  
color = c("orange", "deepskyblue")

# state-dependent process design matrix for plotting
nb = 15; ord = 4; degree = ord-1; nrknots = nb - (ord-2) 
knots = seq(-600, 3000,length = nrknots + 2 * degree)
Z_sdd_plot = spline.des(knots, 1:max(sharks$IDindex), ord, outer.ok=T)$design

par(mfrow = c(4,2))
for(j in 1:7){
  smooth = Z_sdd_plot %*% B1[j,]
  Mu1 = exp(smooth)
  Mu2 = exp(smooth  + mu2[j])
  
  plot(sharks$ODBA[which(sharks$IDnum == j)], pch = 20, xlab = "time", ylab = "ODBA", main = paste("ID", j), bty = "n")
  lines(Mu1, col = color[1], lwd = 2)
  lines(Mu2, col = color[2], lwd = 2)
}

# state process GAM
tod_seq = seq(0, 24, length = 300)
nb = 15; knots = 24 * 0:nb / nb # knots
Z_tod_plot = mgcv::cSplineDes(tod_seq, knots) ## cyclic spline design matrix

temp_seq = seq(17, 35, by = 4)
par(mfrow = c(3,2))

for(i in 1:length(temp_seq)){
  Z_s_plot = cbind(rep(temp_seq[i], nrow(Z_tod_plot)), Z_tod_plot, temp_seq[i] * Z_tod_plot) # design matrix with interaction block
  Gamma = LaMa::tpm_g(Z_s_plot, Beta)
  
  Delta = matrix(NA, dim(Gamma)[3], ncol = N)
  for(k in 1:dim(Gamma)[3]){
    Delta[k,] = LaMa::stationary(Gamma[,,k])
  }
  
  plot(tod_seq, Delta[,2], type = "l", lwd = 2, col = color[2], xlab = "time of day",
       ylab = "probability of being active", main = paste("Temp:", temp_seq[i]), bty = "n")
}





# Marginal ML -------------------------------------------------------------

## design matrices
# state process design matrix
nb = 15 # number of basis functions
knots = 24 * 0:nb / nb # knots
Z_tod = mgcv::cSplineDes(sharks$Hour, knots) ## cyclic spline design matrix
Z_s = cbind(sharks$AvgTemp, Z_tod, matrix(sharks$AvgTemp,nrow=nrow(Z_tod), ncol=ncol(Z_tod)) * Z_tod) # design matrix with interaction block
L = WH:::build_D_mat(nb-1, 2); D_s = t(L)%*%L # difference penalty matrix

# state-dependent process design matrix
nb = 15
ord = 4
degree = ord-1
nrknots = nb - (ord-2) 
knots = seq(-600, 3000,length = nrknots + 2 * degree)
Z_sdd = spline.des(knots, sharks$IDindex, ord, outer.ok=T)$design
L = WH:::build_D_mat(nb, 2); D_sdd = t(L)%*%L # difference penalty matrix

# iterative procedure

maxiter = 100 # maximum number of iterations
tol = 0.05 # relative tolerance for convergence, sufficient
gradtol = 1e-6 # relative gradient tolerance for nlm (1e-6 is the default)
print.level = 0 # print level for nlm

Lambdas = matrix(NA, maxiter+1, 11)
Lambdas[1,] = c(rep(200, 7), # state-dependent GAM
                rep(10, 2), rep(1000, 2)) # state GAM
mods = list()

K = length(trackInd)
REind_sdd = matrix(1:(K*15), nrow = K, ncol = 15, byrow = TRUE)
REind_s = matrix(K*15 + N*K + K + N*(N-1) + 1:(2*14*(N*(N-1))), nrow = N*(N-1)*2, ncol = 14, byrow = TRUE)

mu = c(rep(-0.5, nb*K), rep(0.3, 7))
sigma = rep(0.16, 14)
beta0 = c(-1, -1.5)
b2 = b3 = rep(0, (ncol(Z_s)/2-1) * N*(N-1)) # spline coefs
delta = 0.9
theta.star = c(mu, log(sigma), beta0, b2, b3, delta)

# updating algorithm
for(k in 1:maxiter){
  cat("\n\n- Iteration", k, "-")
  if(k == 1){
    print.level = 2
  } else print.level = 0
  
  t1 = Sys.time()
  mod = nlm(mllk, theta.star, X = sharks, trackInd = trackInd, Z_sdd = Z_sdd, Z_s = Z_s, 
             lambda = Lambdas[k,], D_sdd = D_sdd, D_s = D_s,
             print.level = print.level, iterlim = 1000, gradtol = gradtol, hessian = TRUE)
  Sys.time()-t1
  cat("\nEstimation time:", Sys.time()-t1); cat("\nIterations:", mod$iterations)
  
  theta.star = mod$estimate # saves theta.star as starting value for next iteration
  mods[[k]] = mod # saving model
  J = mod$hessian # assigning hessian
  J_inv = MASS::ginv(J) # inverse penalized hessian
  
  b1 = theta.star[1:(k * ncol(Z_sdd))] # RE vector 1
  p = K*15 + N*K + K + N*(N-1); nb = ncol(Z_s)/2
  
  # updating all penalty strengths state-dependent process
  
  for(i in 1:K){
    edoF = sum(diag(diag(rep(1, nrow(D_sdd))) - Lambdas[k, i] * J_inv[REind_sdd[i,], REind_sdd[i,]] %*% D_sdd))
    penalty = t(theta.star[REind_sdd[i,]]) %*% D_sdd %*% theta.star[REind_sdd[i,]]
    Lambdas[k+1, i] = as.numeric((edoF - 1) / penalty)
  }
  for(i in 1:(N*(N-1)*2)){
    edoF = sum(diag(diag(rep(1, nrow(D_s))) - Lambdas[k, K + i] * J_inv[REind_s[i,], REind_s[i,]] %*% D_s))
    penalty = t(theta.star[REind_s[i,]]) %*% D_s %*% theta.star[REind_s[i,]]
    Lambdas[k+1, K+i] = as.numeric((edoF - 1) / penalty)
  }
  
  cat("\nSmoothing strengths:", round(Lambdas[k+1,], 4))
  if(mean(abs(Lambdas[k+1,] - Lambdas[k,]) / Lambdas[k,]) < tol){
    cat("\n\n")
    break
  }
}

par(mfrow = c(4,3))
Lambdas = as.matrix(na.omit(Lambdas))
for(i in 1:ncol(Lambdas)){
  plot(Lambdas[,i], type = "l", lwd = 2, main = paste("lambda", i))
}











# Simpler model: GAM in sdp - but sin cos in state process ----------------

mllk2 = function(theta.star, X, trackInd, Z_sd, Z_s, lambda, D_sd){
  K = length(trackInd); N = 2
  
  # parameters for the GAM in the state-dependent process
  nb_sd = ncol(Z_sd)
  b1 = theta.star[1:(K * nb_sd)]; p = K * nb_sd # RE vector
  B1 = matrix(b1, nrow = K, ncol = nb_sd, byrow = TRUE)
  mu2 = theta.star[p + 1:K]; p = p + K # individual-specific deviations from state 1 mean
  Sigma = matrix(exp(theta.star[p + 1:(N*K)]), nrow = K, ncol = N); p = p+N*K 
  
  # parameters for the transition probabilities
  nb_s = ncol(Z_s)
  beta = matrix(theta.star[p + 1:(N*(N-1)*(nb_s + 1))], nrow = N*(N-1), ncol = nb_s + 1)
  p = p + N*(N-1)*(nb_s + 1)
  
  # initial distribution
  delta = c(1, exp(theta.star[p + 1:(N-1)]))
  delta = delta / sum(delta)
  
  # calculating the transition probability matrix
  Gamma = LaMa::tpm_g(Z_s, beta)
  
  # calculating the state-dependent probabilities
  smooth = rowSums(Z_sd * B1[X$IDnum,]) # GAM means for state 1 and each individual (on log scale)
  Mu1 = exp(smooth)
  Mu2 = exp(smooth  + mu2[X$IDnum]) # means for state two by individual-specific deviations
  sigma1 = Sigma[X$IDnum, 1]
  sigma2 = Sigma[X$IDnum, 2]
  
  allprobs = matrix(1, nrow = nrow(X), ncol = N)
  ind = which(!is.na(X$ODBA))
  allprobs[ind,1] = dgamma(X$ODBA[ind], shape = Mu1[ind]^2/sigma1[ind]^2, scale = sigma1[ind]^2/Mu1[ind])
  allprobs[ind,2] = dgamma(X$ODBA[ind], shape = Mu2[ind]^2/sigma2[ind]^2, scale = sigma2[ind]^2/Mu2[ind])
  
  # calculating the penalty term
  pen = t(b1) %*% kronecker(diag(lambda), D_sd) %*% b1
  
  # forward algorithm for separate tracks + penalty term
  -LaMa::forward_g(delta, Gamma, allprobs, trackInd) + 0.5 * pen
}

# state process design matrix
Z_trig = LaMa::trigBasisExp(sharks$Hour, L = 24, degree = 2)
Z_s = cbind(sharks$AvgTemp, Z_trig, matrix(sharks$AvgTemp, nrow=nrow(Z_trig), ncol=ncol(Z_trig)) * Z_trig) # design matrix with interaction block

# state-dependent process design matrix
ord = 4; nb = 15
degree = ord-1
nrknots = nb - (ord-2) 
knots = seq(-600, 3000, length = nrknots + 2 * degree)
Z_sd = spline.des(knots, sharks$IDindex, ord, outer.ok=T)$design
L = WH:::build_D_mat(nb, 2)
D_sd = t(L)%*%L # difference penalty matrix

trackInd = LaMa::calc_trackInd(sharks$SharkID)
K = length(trackInd); N = 2

mu = c(rep(-0.5, nb*K), rep(0.3, 7))
sigma = rep(0.16, 14)
beta0 = c(-2, -2)
beta = rep(0, ncol(Z_s) * N*(N-1)) # spline coefs
delta = 0.9
theta.star = c(mu, log(sigma), beta0, beta, delta)

lambda = c(10,50,50,30,30,10,10)*2

t1 = Sys.time()
mod_trig = nlm(mllk2, theta.star, X = sharks, trackInd = trackInd, Z_sd = Z_sd, Z_s = Z_s, 
               lambda = lambda, D_sd = D_sd,
               print.level = 2, iterlim = 5000)
Sys.time()-t1

theta.star = mod_trig$estimate

K = length(trackInd); N = 2

# parameters for the GAM in the state-dependent process
nb_sd = ncol(Z_sd)
b1 = theta.star[1:(K * nb_sd)]; p = K * nb_sd # RE vector
B1 = matrix(b1, nrow = K, ncol = nb_sd, byrow = TRUE)
mu2 = theta.star[p + 1:K]; p = p + K # individual-specific deviations from state 1 mean
Sigma = matrix(exp(theta.star[p + 1:(N*K)]), nrow = K, ncol = N); p = p+N*K 

# parameters for the transition probabilities
nb_s = ncol(Z_s)
beta = matrix(theta.star[p + 1:(N*(N-1)*(nb_s + 1))], nrow = N*(N-1), ncol = nb_s + 1)
p = p + N*(N-1)*(nb_s + 1)

# initial distribution
delta = c(1, exp(theta.star[p + 1:(N-1)]))
delta = delta / sum(delta)



## visualizing results  

### state-dependent distribution
color = c("orange", "deepskyblue")

# state-dependent process design matrix for plotting
nb = 15; ord = 4; degree = ord-1; nrknots = nb - (ord-2) 
knots = seq(-600, 3000,length = nrknots + 2 * degree)
Z_sdd_plot = spline.des(knots, 1:max(sharks$IDindex), ord, outer.ok=T)$design

par(mfrow = c(4,2))
for(j in 1:7){
  smooth = Z_sdd_plot %*% B1[j,]
  Mu1 = exp(smooth)
  Mu2 = exp(smooth  + mu2[j])
  
  plot(sharks$ODBA[which(sharks$IDnum == j)], pch = 20, xlab = "time", ylab = "ODBA", main = paste("ID", j), bty = "n")
  lines(Mu1, col = color[1], lwd = 2)
  lines(Mu2, col = color[2], lwd = 2)
}
# looking good


### state process
tod_seq = seq(0, 24, length = 300)
Z_trig_plot = LaMa::trigBasisExp(tod_seq, L = 24, degree = 2)

temp_seq = seq(17, 35, by = 3)

par(mfrow = c(2,4))
for(i in 1:length(temp_seq)){
  Z_s_plot = cbind(rep(temp_seq[i], nrow(Z_trig_plot)), Z_trig_plot, 
                   matrix(temp_seq[i], nrow=nrow(Z_trig_plot), ncol=ncol(Z_trig_plot)) * Z_trig_plot)
  Gamma_plot = LaMa::tpm_g(Z_s_plot, beta)
  
  Delta = matrix(NA, dim(Gamma_plot)[3], ncol = N)
  for(k in 1:dim(Gamma_plot)[3]){
    Delta[k,] = LaMa::stationary(Gamma_plot[,,k])
  }
  plot(tod_seq, Delta[,2], type = "l", lwd = 2, col = color[2], xlab = "time of day",
       ylab = "probability of being active", main = paste("Temp:", temp_seq[i]), bty = "n")
}
# also looking as in the paper





# Now again, non-parametric modelling of the state process ----------------

mllk3 = function(theta.star, X, trackInd, Z_sd, Z_s, lambda, D_sd, D_s){
  K = length(trackInd); N = 2
  
  # parameters for the GAM in the state-dependent process
  nb_sd = ncol(Z_sd)
  b1 = theta.star[1:(K * nb_sd)]; p = K * nb_sd # RE vector
  B1 = matrix(b1, nrow = K, ncol = nb_sd, byrow = TRUE)
  mu2 = theta.star[p + 1:K]; p = p + K # individual-specific deviations from state 1 mean
  Sigma = matrix(exp(theta.star[p + 1:(N*K)]), nrow = K, ncol = N); p = p + N*K 
  
  # parameters for the transition probabilities
  nb_s = (ncol(Z_s)-1)/2
  beta_fix = matrix(theta.star[p + 1:(N*(N-1)*2)], nrow = N*(N-1), ncol = 2); p = p + N*(N-1)*2
  b2 = theta.star[p + 1:((nb_s-1) * N*(N-1) * 2)]; p = p + (nb_s-1)*N*(N-1)*2 # RE vector 2
  Beta = matrix(0, nrow = N*(N-1), ncol = 2*nb_s+2)
  Beta[,1:2] = beta_fix
  Beta[,-c(1:2, 2+c(1, nb_s+1))] = matrix(b2, nrow = N*(N-1), ncol = (nb_s-1)*2, byrow = TRUE)

  # initial distribution
  delta = c(1, exp(theta.star[p + 1:(N-1)]))
  delta = delta / sum(delta)
  
  # calculating the transition probability matrix
  Gamma = LaMa::tpm_g(Z_s, Beta)
  
  # calculating the state-dependent probabilities
  smooth = rowSums(Z_sd * B1[X$IDnum,]) # GAM means for state 1 and each individual (on log scale)
  Mu1 = exp(smooth)
  Mu2 = exp(smooth  + mu2[X$IDnum]) # means for state two by individual-specific deviations
  sigma1 = Sigma[X$IDnum, 1]
  sigma2 = Sigma[X$IDnum, 2]
  
  allprobs = matrix(1, nrow = nrow(X), ncol = N)
  ind = which(!is.na(X$ODBA))
  allprobs[ind,1] = dgamma(X$ODBA[ind], shape = Mu1[ind]^2/sigma1[ind]^2, scale = sigma1[ind]^2/Mu1[ind])
  allprobs[ind,2] = dgamma(X$ODBA[ind], shape = Mu2[ind]^2/sigma2[ind]^2, scale = sigma2[ind]^2/Mu2[ind])
  
  # calculating the penalty term
  pen = t(b1) %*% kronecker(diag(lambda[1:K]), D_sd) %*% b1 + 
    # t(c(Beta[1,-(1:2)], Beta[2,-(1:2)])) %*% kronecker(diag(lambda[K + 1:(N*(N-1)*2)]), D_s) %*% c(Beta[1,-(1:2)], Beta[2,-(1:2)])
    t(b2) %*% kronecker(diag(lambda[K + 1:(N*(N-1)*2)]), D_s) %*% b2
  
  # forward algorithm for separate tracks + penalty term
  -LaMa::forward_g(delta, Gamma, allprobs, trackInd) + 0.5 * pen
}

# state process design matrix
nb_s = 10 # number of basis functions
knots = 24 * 0:nb_s / nb_s # knots
Z_tod = mgcv::cSplineDes(sharks$Hour, knots) ## cyclic spline design matrix
Z_s = cbind(sharks$AvgTemp, Z_tod, 
            matrix(sharks$AvgTemp, nrow=nrow(Z_tod), ncol=ncol(Z_tod)) * Z_tod)
L_bar = WH:::build_D_mat(nb_s, 2)
# custom modification for 2 more (circular) differences 
# it seems to be okay to make this nb_1-1 times nb_1-1 (see Schellhase and Kauermann)
L = matrix(0, nb_s, nb_s)
L[1:(nb_s-2),1:(nb_s)] = L_bar
L[nb_s-1, 1] = 1; L[nb_s-1, (nb_s-1):(nb_s)] = c(1, -2)
L[nb_s,1:2] = c(-2, 1); L[nb_s, nb_s] = 1
D_s = t(L[,-1])%*%L[,-1] # difference penalty matrix
# this has full rank and is invertible

# state-dependent process design matrix
ord = 4; nb_sd = 15
degree = ord-1
nrknots = nb_sd - (ord-2) 
knots_sd = seq(-600, 3000, length = nrknots + 2 * degree)
Z_sd = spline.des(knots_sd, sharks$IDindex, ord, outer.ok = T)$design
L = WH:::build_D_mat(nb_sd, 2)
D_sd = t(L)%*%L # difference penalty matrix

trackInd = LaMa::calc_trackInd(sharks$SharkID)
K = length(trackInd); N = 2

mu = c(rep(-0.5, nb_sd * K), rep(0.3, 7))
sigma = rep(0.16, 14)
beta0 = c(-1.5, -1.5)
beta1 = c(0, 0)
b2 = rep(0, (nb_s-1) * N*(N-1) *2) # spline coefs
delta = 0.9
theta.star = c(mu, log(sigma), beta0, beta1, b2, delta)

lambda = c(c(10,50,50,30,30,10,10)*2,
           rep(30, 2), rep(5000, 2))

t1 = Sys.time()
mod3 = nlm(mllk3, theta.star, X = sharks, trackInd = trackInd, Z_sd = Z_sd, Z_s = Z_s, 
           lambda = lambda, D_sd = D_sd, D_s = D_s,
           print.level = 2, iterlim = 2000, hessian = TRUE)
Sys.time()-t1

theta.star = mod3$estimate

K = length(trackInd); N = 2

# parameters for the GAM in the state-dependent process
nb_sd = ncol(Z_sd)
b1 = theta.star[1:(K * nb_sd)]; p = K * nb_sd # RE vector
B1 = matrix(b1, nrow = K, ncol = nb_sd, byrow = TRUE)
mu2 = theta.star[p + 1:K]; p = p + K # individual-specific deviations from state 1 mean
Sigma = matrix(exp(theta.star[p + 1:(N*K)]), nrow = K, ncol = N); p = p + N*K 

# parameters for the transition probabilities
nb_s = (ncol(Z_s)-1)/2
beta_fix = matrix(theta.star[p + 1:(N*(N-1)*2)], nrow = N*(N-1), ncol = 2); p = p + N*(N-1)*2
b2 = theta.star[p + 1:((nb_s-1) * N*(N-1) * 2)]; p = p + (nb_s-1)*N*(N-1)*2 # RE vector 2
Beta = matrix(0, nrow = N*(N-1), ncol = 2*nb_s+2)
Beta[,1:2] = beta_fix
Beta[,-c(1:2, 2+c(1, nb_s+1))] = matrix(b2, nrow = N*(N-1), ncol = (nb_s-1)*2, byrow = TRUE)
beta = Beta[,-(1:2)]

par(mfrow = c(2,2))
plot(beta[1,1:nb_s], type = "l")
plot(beta[1,nb_s+1:nb_s], type = "l")
plot(beta[2,1:nb_s], type = "l")
plot(beta[2,nb_s+1:nb_s], type = "l")

# initial distribution
delta = c(1, exp(theta.star[p + 1:(N-1)]))
delta = delta / sum(delta)


## visualizing results
### state-dependent distribution
color = c("orange", "deepskyblue")

# state-dependent process design matrix for plotting
nb = 15; ord = 4; degree = ord-1; nrknots = nb - (ord-2) 
knots = seq(-600, 3000,length = nrknots + 2 * degree)
Z_sdd_plot = spline.des(knots, 1:max(sharks$IDindex), ord, outer.ok=T)$design

par(mfrow = c(3,3))
for(j in 1:7){
  smooth = Z_sdd_plot %*% B1[j,]
  Mu1 = exp(smooth)
  Mu2 = exp(smooth  + mu2[j])
  
  plot(sharks$ODBA[which(sharks$IDnum == j)], pch = 20, xlab = "time", ylab = "ODBA", main = paste("ID", j), bty = "n")
  lines(Mu1, col = color[1], lwd = 2)
  lines(Mu2, col = color[2], lwd = 2)
}
# looking good

## state process
tod_seq = seq(0, 24, length = 300)
nb_s = 10 # number of basis functions
knots = 24 * 0:nb_s / nb_s # knots
Z_tod_plot = mgcv::cSplineDes(tod_seq, knots) ## cyclic spline design matrix

temp_seq = seq(17, 35, by = 3)

par(mfrow = c(2,4))
for(i in 1:length(temp_seq)){
  Z_s_plot = cbind(rep(temp_seq[i], nrow(Z_tod_plot)), Z_tod_plot, 
                   matrix(temp_seq[i], nrow=nrow(Z_tod_plot), ncol=ncol(Z_tod_plot)) * Z_tod_plot)
  Gamma_plot = LaMa::tpm_g(Z_s_plot, Beta)
  
  Delta = matrix(NA, dim(Gamma_plot)[3], ncol = N)
  for(k in 1:dim(Gamma_plot)[3]){
    Delta[k,] = LaMa::stationary(Gamma_plot[,,k])
  }
  plot(tod_seq, Delta[,2], type = "l", lwd = 2, col = color[2], xlab = "time of day",
       ylab = "probability of being active", main = paste("Temp:", temp_seq[i]), bty = "n")
}


temp = 20
Z_s_plot = cbind(rep(temp, nrow(Z_tod_plot)), Z_tod_plot, 
                 matrix(temp, nrow=nrow(Z_tod_plot), ncol=ncol(Z_tod_plot)) * Z_tod_plot)
Gamma_plot = LaMa::tpm_g(Z_s_plot, Beta)

par(mfrow = c(2,2))
for(i in 1:N){
  for(j in 1:N){
    plot(Gamma_plot[i,j,], type = "l", lwd = 2)
  }
}

par(mfrow = c(2,2))
for(i in 1:N){
    Eta = cbind(1, Z_s_plot) %*% Beta[i,]
    plot(Eta, type = "l", lwd = 2)
    plot(Beta[i,-(1:2)])
}






# Marginal ML -------------------------------------------------------------

# state process design matrix
nb_s = 10 # number of basis functions
knots = 24 * 0:nb_s / nb_s # knots
Z_tod = mgcv::cSplineDes(sharks$Hour, knots) ## cyclic spline design matrix
Z_s = cbind(sharks$AvgTemp, Z_tod, 
            matrix(sharks$AvgTemp, nrow=nrow(Z_tod), ncol=ncol(Z_tod)) * Z_tod)
L_bar = WH:::build_D_mat(nb_s, 2)
# custom modification for 2 more (circular) differences 
# it seems to be okay to make this nb_1-1 times nb_1-1 (see Schellhase and Kauermann)
L = matrix(0, nb_s, nb_s)
L[1:(nb_s-2),1:(nb_s)] = L_bar
L[nb_s-1, 1] = 1; L[nb_s-1, (nb_s-1):(nb_s)] = c(1, -2)
L[nb_s,1:2] = c(-2, 1); L[nb_s, nb_s] = 1
D_s = t(L[,-1])%*%L[,-1] # difference penalty matrix

# state-dependent process design matrix
ord = 4; nb_sd = 15
degree = ord-1
nrknots = nb_sd - (ord-2) 
knots_sd = seq(-600, 3000, length = nrknots + 2 * degree)
Z_sd = spline.des(knots_sd, sharks$IDindex, ord, outer.ok = T)$design
L = WH:::build_D_mat(nb_sd, 2)
D_sd = t(L)%*%L # difference penalty matrix

trackInd = LaMa::calc_trackInd(sharks$SharkID)
K = length(trackInd); N = 2

# iterative procedure

maxiter = 100 # maximum number of iterations
tol = 0.05 # relative tolerance for convergence, sufficient
gradtol = 1e-6 # relative gradient tolerance for nlm (1e-6 is the default)
print.level = 0 # print level for nlm

Lambdas = matrix(NA, maxiter+1, 11)
Lambdas[1,] = c(c(10,50,50,30,30,10,10)*3,
                rep(30, 2), rep(5000, 2))
mods = list()

K = length(trackInd)
REind_sd = matrix(1:(K*nb_sd), nrow = K, ncol = nb_sd, byrow = TRUE)
REind_s = matrix(K*nb_sd + K + K*N + N*(N-1)*2 + 1:(2*(nb_s-1)*(N*(N-1))), nrow = N*(N-1)*2, ncol = nb_s-1, byrow = TRUE)

mu = c(rep(-0.5, nb_sd * K), rep(0.3, 7))
sigma = rep(0.16, 14)
beta0 = c(-1.5, -1.5)
beta1 = c(0, 0)
b2 = rep(0, (nb_s-1) * N*(N-1) *2) # spline coefs
delta = 0.9
theta.star = c(mu, log(sigma), beta0, beta1, b2, delta)

# updating algorithm
for(k in 1:maxiter){
  cat("\n\n- Iteration", k, "-")
  if(k == 1){
    print.level = 2
  } else print.level = 0
  
  t1 = Sys.time()
  mod = nlm(mllk3, theta.star, X = sharks, trackInd = trackInd, Z_sd = Z_sd, Z_s = Z_s, 
            lambda = Lambdas[k,], D_sd = D_sd, D_s = D_s,
            print.level = print.level, iterlim = 2000, gradtol = gradtol, hessian = TRUE)
  Sys.time()-t1
  cat("\nEstimation time:", Sys.time()-t1); cat("\nIterations:", mod$iterations)
  
  theta.star = mod$estimate # saves theta.star as starting value for next iteration
  mods[[k]] = mod # saving model
  J = mod$hessian # assigning hessian
  J_inv = MASS::ginv(J) # inverse penalized hessian
  
  # updating all penalty strengths state-dependent process
  for(i in 1:K){
    edoF = sum(diag(diag(rep(1, nrow(D_sd))) - Lambdas[k, i] * J_inv[REind_sd[i,], REind_sd[i,]] %*% D_sd))
    penalty = t(theta.star[REind_sd[i,]]) %*% D_sd %*% theta.star[REind_sd[i,]]
    Lambdas[k+1, i] = as.numeric((edoF - 1) / penalty)
  }
  for(i in 1:(N*(N-1)*2)){
    edoF = sum(diag(diag(rep(1, nrow(D_s))) - Lambdas[k, K + i] * J_inv[REind_s[i,], REind_s[i,]] %*% D_s))
    penalty = t(theta.star[REind_s[i,]]) %*% D_s %*% theta.star[REind_s[i,]]
    Lambdas[k+1, K+i] = as.numeric(edoF / penalty) # no -1 because cyclic penalty matrix has full rank
  }
  
  cat("\nSmoothing strengths:", round(Lambdas[k+1,], 4))
  if(mean(abs(Lambdas[k+1,] - Lambdas[k,]) / Lambdas[k,]) < tol){
    cat("\n\n")
    break
  }
}

par(mfrow = c(4,3))
Lambdas = as.matrix(na.omit(Lambdas))
for(i in 1:ncol(Lambdas)){
  plot(Lambdas[,i], type = "l", lwd = 2, main = paste("lambda", i))
}