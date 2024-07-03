
# Loading the data --------------------------------------------------------

load("./data/fruitflies.RData")


# Starting with simple homogeneous model ----------------------------------

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



# Spline model with random intercept --------------------------------------

mllk_fliesRE = function(theta.star, X, N=2, Z, lambda, S, trackInd){
  # assigning state-dependent parameters
  mu = exp(theta.star[1:N])
  phi = exp(theta.star[N + 1:N]); p = 2*N
  
  # assigning state parameters, beta0 intercept, bLD and bDD matrix
  nb = ncol(Z) - 1
  beta0 = matrix(theta.star[p + 1:(2*N*(N-1))], nrow = 2); p = p + 2*N*(N-1) # intercepts
  bLD = theta.star[p + 1:((N*(N-1))*nb)]; p = p + (N*(N-1))*nb # both random intercept and spline coefs LD
  bDD = theta.star[p + 1:((N*(N-1))*nb)] # both, random intercept and spline coefs LD
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
  # I know this is terrible
  pen = t(as.numeric(beta[1 + 1:nAnimals,, 1])) %*% kronecker(diag(lambda[1:2]), S[[1]]) %*% as.numeric(beta[1 + 1:nAnimals,, 1]) + # random intercepts LD
    t(as.numeric(beta[1 + 1:nAnimals,, 2])) %*% kronecker(diag(lambda[3:4]), S[[1]]) %*% as.numeric(beta[1 + 1:nAnimals,, 2]) + # random intercepts DD
    t(as.numeric(beta[1 + nAnimals + 1:nrow(S[[2]]),, 1])) %*% kronecker(diag(lambda[5:6]), S[[2]]) %*% as.numeric(beta[1 + nAnimals + 1:nrow(S[[2]]),, 1]) + # splines LD
    t(as.numeric(beta[1 + nAnimals + 1:nrow(S[[2]]),, 2])) %*% kronecker(diag(lambda[7:8]), S[[2]]) %*% as.numeric(beta[1 + nAnimals + 1:nrow(S[[2]]),, 2]) # splines DD
  
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
Z = gam_pre$X # design matrix containing dummies and spline design
S = gam_pre$S # list of 2 penalty matrices, first for random intercept, second for spline

# initial parameter vector
theta.star = c(log(c(mu, phi)), # state-dependent parameters
               rep(-2,2), # intercept LD
               rep(-2,2), # intercept DD
               rep(0, 2 * N*(N-1) * (ncol(Z)-1))) # RE and spline coefs

# calculating trackInd for LaMa
trackInd = LaMa::calc_trackInd(as.vector(data$ID))


# Model fitting via marginal ML --------

# hyperparamters
maxiter = 100 # maximum number of iterations
tol = 0.01 # relative tolerance for convergence, sufficient
alpha = 0.85 # exponential smoothing parameter for penalty strengths

# matrix of penalty strengths/ inverse variances
Lambdas = matrix(NA, maxiter+1, N*(N-1)*4)
Lambdas[1,] = c(rep(50, 4), rep(100, 4))
mods = list() # modlist 

# building a list that contains the indices for each random effect and spline vector
p = 2*N*(N-1)+2*N
nAnimals = length(unique(data$ID))
REind = list()
for(i in 1:4) REind[[i]] = p + (i-1)*nAnimals + 1:nAnimals # random intercepts
for(i in 5:8) REind[[i]] = p + 4 * nAnimals + (i-5)*nb + 1:nb # splines

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
    
    edoF = sum(diag(diag(rep(1, nrow(S_i))) - Lambdas[k, i] * J_inv[REind[[i]], REind[[i]]] %*% S_i)) # effective # of parameters for this random effect/ spline
    penalty = t(theta.star[REind[[i]]]) %*% S_i %*% theta.star[REind[[i]]] # penalty for this random effect/ spline
    lambda_new = as.numeric(edoF / penalty) # update
    Lambdas[k+1, i] = alpha*lambda_new + (1-alpha)*Lambdas[k, i] # exponential smoothing of update
  }
  Lambdas[k+1, which(Lambdas[k+1,] < 0)] = 0 # constraining lambda to be non-negative (numerical reasons)
  
  cat("\nSmoothing strengths:", round(Lambdas[k+1, ], 4))
  if(mean(abs(Lambdas[k+1,] - Lambdas[k,]) / Lambdas[k,]) < tol){
    cat("\n\n")
    break
  }
}
cat("\nTotal estimation time:", Sys.time()-T1, "\n")

# trace plots of penalty strenghts
Lambdas = as.matrix(na.omit(Lambdas))
lambda_hat = Lambdas[nrow(Lambdas),] # final lambda
par(mfrow = c(2,4))
for(i in 1:ncol(Lambdas)){
  plot(Lambdas[,i], type = "l", lwd = 2, main = paste("lambda", i), bty = "n")
}

# standard deviations of the state-dependent intercepts
sqrt(1 / lambda_hat[1:4])

# assinging estimated parameter:
mod_final = mods[[length(mods)]]
theta.star = mod_final$estimate

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

n = 300 # plotting smoothness
nb = 12
todseq = seq(0, 24, length.out = n)
Delta_mean = array(dim = c(n, 2, 2))
Delta = array(dim = c(n, 2, 2, nAnimals))

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
# this takes some time
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
cond = 1 # 1 = LD, 2 = DD
state = 2

par(mfrow = c(1,1))
plot(NA, xlim = c(0,24), ylim = c(0,1), xlab = "time of day", ylab = "Pr(active)", 
     bty = "n", main = conditions[cond], xaxt = "n")
axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))

for(a in 1:nAnimals){
  lines(seq(0,24,length=n), Delta[,state,cond,a], col = scales::alpha(a, 0.4))
}
lines(seq(0,24,length=n), Delta_mean[,state,cond], col = "black", lwd = 4)

