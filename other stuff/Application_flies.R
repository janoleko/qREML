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

mllk_flies = function(theta.star, X, N=2, Z, lambda, D, trackInd){
  mu = exp(theta.star[1:N])
  phi = exp(theta.star[N + 1:N]); p = 2*N
  
  beta0 = matrix(theta.star[p + 1:(2*N*(N-1))], nrow = 2); p = p + 2*N*(N-1)
  nb = ncol(Z)
  bLD = theta.star[p + 1:((N*(N-1))*(nb-1))]; p=p+(N*(N-1))*(nb-1)
  bDD = theta.star[p + 1:((N*(N-1))*(nb-1))]
  b = c(bLD, bDD)
  beta = array(0, dim = c(N, nb+1, 2))
  beta[,-(1:2),1] = matrix(bLD, nrow = N, ncol = nb-1, byrow = TRUE) #1st slice LD
  beta[,-(1:2),2] = matrix(bDD, nrow = N, ncol = nb-1, byrow = TRUE) #2nd slice DD
  beta[,1,1] = beta0[,1]
  beta[,1,2] = beta0[,2]
  
  Gamma = array(dim=c(N,N,48,2))
  Gamma[,,,1] = LaMa::tpm_g(Z, beta[,,1])
  Gamma[,,,2] = LaMa::tpm_g(Z, beta[,,2])
  
  ind = which(X$condition == "LD")
  Gamma_track = array(dim = c(N,N,nrow(X)))
  Gamma_track[,,ind] = Gamma[,,X$tod[ind],1]
  Gamma_track[,,-ind] = Gamma[,,X$tod[-ind],2]

  delta = LaMa::stationary_p(Gamma[,,,1], t=X$tod[1])
  
  allprobs = matrix(1, nrow = nrow(X), ncol = N)
  ind = which(!is.na(X$activity))
  for (j in 1:N){
    allprobs[ind, j] = dnbinom(X$activity[ind], mu=mu[j], size=1/phi[j])
  }
  
  pen = b %*% kronecker(diag(lambda), D) %*% b
  
  -LaMa::forward_g(delta, Gamma_track, allprobs, trackInd) + 0.5 * pen
}


nb = 20 # number of basis functions
x = 1:48 / 2
k = 24 * 0:nb / nb # knots
Z = mgcv::cSplineDes(x, k) ## cyclic spline design matrix

L = WH:::build_D_mat(nb-1, 2)
D = t(L)%*%L


theta.star = c(log(c(mu, phi)), 
               rep(-2,2), rep(-3,2), rep(0, 2*N*(N-1)*(nb-1)))

trackInd = LaMa::calc_trackInd(as.vector(data$ID))

# iterative procedure

maxiter = 100 # maximum number of iterations
tol = 0.01 # relative tolerance for convergence, sufficient
gradtol = 1e-5 # relative gradient tolerance for nlm (1e-6 is the default)
print.level = 0 # print level for nlm

Lambdas = matrix(NA, maxiter+1, N*(N-1)*2)
Lambdas[1,] = rep(60,4)
mods = list()

p = 2*N*(N-1)+2*N
REind = matrix(1:(2*(nb-1)*N*(N-1)), nrow = N*(N-1)*2, byrow = TRUE) # each row is the index of one RE

# updating algorithm
for(k in 1:maxiter){
  cat("\n\n- Iteration", k, "-")
  if(k == 1) print.level = 2
  else print.level = 0
  
  t1 = Sys.time()
  mod = nlm(mllk_flies, theta.star, N=2, X=data, Z=Z, lambda = Lambdas[k,], D=D, trackInd=trackInd,
            iterlim = 1000, print.level = print.level, hessian = TRUE, gradtol = gradtol)
  Sys.time()-t1
  cat("\nEstimation time:", Sys.time()-t1)
  cat("\nIterations:", mod$iterations)
  
  theta.star = mod$estimate # saves theta.star as starting value for next iteration
  mods[[k]] = mod # saving model
  
  J_p = mod$hessian # assigning hessian
  b = theta.star[p + 1:(2*N*(N-1)*(nb-1))]
  
  # updating all penalty strengths
  for(i in 1:(N*(N-1)*2)){
    edoF = sum(diag( # trace
      diag(rep(1, nrow(D))) - Lambdas[k, i] * MASS::ginv(J_p)[p+REind[i,], p+REind[i,]] %*% D
    ))
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

mod_final = nlm(mllk_flies, theta.star, N=2, X=data, Z=Z, lambda = Lambdas[nrow(Lambdas),], D=D, trackInd=trackInd,
          iterlim = 1000, print.level = 2, hessian = TRUE)
theta.star = mod_final$estimate

# plot
mu = exp(theta.star[1:N])
phi = exp(theta.star[N + 1:N]); p = 2*N

beta0 = matrix(theta.star[p + 1:(2*N*(N-1))], nrow = 2); p = p + 2*N*(N-1)
nb = ncol(Z)
bLD = theta.star[p + 1:((N*(N-1))*(nb-1))]; p=p+(N*(N-1))*(nb-1)
bDD = theta.star[p + 1:((N*(N-1))*(nb-1))]
b = c(bLD, bDD)
beta = array(0, dim = c(N, nb+1, 2))
beta[,-(1:2),1] = matrix(bLD, nrow = N, ncol = nb-1, byrow = TRUE) #1st slice LD
beta[,-(1:2),2] = matrix(bDD, nrow = N, ncol = nb-1, byrow = TRUE) #2nd slice DD
beta[,1,1] = beta0[,1]
beta[,1,2] = beta0[,2]

Gamma = array(dim=c(N,N,48,2))
Gamma[,,,1] = LaMa::tpm_g(Z, beta[,,1])
Gamma[,,,2] = LaMa::tpm_g(Z, beta[,,2])

par(mfrow = c(1,1))
plot(Gamma[1,2,,2], type = "l", col = color[1], lwd = 2)

Delta1 = LaMa::stationary_p(Gamma[,,,1])
Delta2 = LaMa::stationary_p(Gamma[,,,2])
plot(Delta1[,2], type = "l", col = color[2], lwd = 2)

n = 500
k = 24 * 0:nb / nb # knots
Delta = matrix(NA, n, 2)
todseq = seq(0, 24, length.out = n)
for(t in 1:length(todseq)){
  time = (todseq[t] + 1:24) %% 24
  Gamma_t = LaMa::tpm_g(mgcv::cSplineDes(time, k), beta[,,1])
  Delta[t,] = LaMa::stationary_p(Gamma_t, 1)
}
par(mfrow = c(1,1))
plot(todseq, Delta[,2], type = "l", lwd = 2, col = color[2])


length(mods)
par(mfrow = c(3,3))
for(m in 1:9){
  theta.star = mods[[m]]$estimate
  p = 2*N
  beta0 = matrix(theta.star[p + 1:(2*N*(N-1))], nrow = 2); p = p + 2*N*(N-1)
  nb = ncol(Z)
  bLD = theta.star[p + 1:((N*(N-1))*(nb-1))]; p=p+(N*(N-1))*(nb-1)
  bDD = theta.star[p + 1:((N*(N-1))*(nb-1))]
  b = c(bLD, bDD)
  beta = array(0, dim = c(N, nb+1, 2))
  beta[,-(1:2),1] = matrix(bLD, nrow = N, ncol = nb-1, byrow = TRUE) #1st slice LD
  beta[,-(1:2),2] = matrix(bDD, nrow = N, ncol = nb-1, byrow = TRUE) #2nd slice DD
  beta[,1,1] = beta0[,1]
  beta[,1,2] = beta0[,2]
  
  Gamma = array(dim=c(N,N,48,2))
  Gamma[,,,1] = LaMa::tpm_g(Z, beta[,,1])
  Gamma[,,,2] = LaMa::tpm_g(Z, beta[,,2])
  
  Delta1 = LaMa::stationary_p(Gamma[,,,1])
  Delta2 = LaMa::stationary_p(Gamma[,,,2])
  plot(Delta1[,2], type = "l", col = color[2], lwd = 2)
}

