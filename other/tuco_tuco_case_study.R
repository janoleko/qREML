load("./data/vedba_tucos.RData")

tuco = vedba_tucos[[2]]
tuco[which(tuco==0)] = NA
ltuco = log(tuco)
length(ltuco)

3456000 / length(ltuco)

# thinning the data
ltuco = ltuco[seq(1, length(ltuco), by = 10)]

hist(ltuco, breaks = 100, prob = T)

mllk_simple = function(theta.star, x, N){
  Gamma = LaMa::tpm(theta.star[1:(N*(N-1))]); p = N*(N-1)
  delta = LaMa::stationary(Gamma)
  mu = theta.star[p + 1:N]
  sigma = exp(theta.star[p + N + 1:N])
  
  allprobs = matrix(1, nrow = length(x), ncol = N)
  ind = which(!is.na(x))
  for (j in 1:N){
    allprobs[ind, j] = dnorm(x[ind], mu[j], sigma[j])
  }

  -LaMa::forward(delta, Gamma, allprobs)
}

# first fit simple model

mu0 = c(-4.3, -2.4, -1)
sigma0 = c(0.5,1,0.5)
hist(ltuco, breaks = 100, prob = T)
N = 3
for(j in 1:N){
  curve(0.3*dnorm(x, mu0[j], sigma0[j]), add = T, col = j)
}

theta.star = c(rep(-3, 6),
               mu0, log(sigma0))

mod_sim_tuco = nlm(mllk_simple, theta.star, x = ltuco, N = 3, print.level = 2)

theta.star = mod_sim_tuco$estimate
N=3
Gamma = LaMa::tpm(theta.star[1:(N*(N-1))]); p = N*(N-1)
delta = LaMa::stationary(Gamma)
mu = theta.star[p + 1:N]
sigma = exp(theta.star[p + N + 1:N])

par(mfrow = c(1,1))
color = c("orange", "deepskyblue", "seagreen2")
hist(ltuco, breaks = 100, prob = T)
for(j in 1:N){
  curve(delta[j]*dnorm(x, mu[j], sigma[j]), add = T, col = color[j], lwd = 2)
}


# flexible model via marginal ML

mllk_tuco = function(theta.star, x, N, B, lambda, D){
  Gamma = LaMa::tpm(theta.star[1:(N*(N-1))]); p = N*(N-1)
  delta = LaMa::stationary(Gamma)
  
  nb = ncol(B) # number of basis functions
  b = theta.star[p + 1:((nb-1)*N)]
  Beta = matrix(0, nrow = N, ncol = nb) # spline coefficent matrix
  Beta[,-1] = matrix(b, nrow = N, ncol = nb-1, byrow = TRUE)
  A = exp(Beta)
  A = A / rowSums(A)
  
  allprobs = matrix(1, nrow = length(x), ncol = N)
  ind = which(!is.na(x))
  for (j in 1:N){
    allprobs[ind, j] = rowSums(B[ind,] * matrix(A[j,], nrow=length(ind), ncol=nb, byrow=T))
  }
  
  pen = t(b) %*% kronecker(diag(lambda), D) %*% b
  
  -LaMa::forward(delta, Gamma, allprobs) + 0.5 * pen
}

# building B-spline design matrix
build_desmat = function(x, K, ord = 4){
  nb = K+1
  minX = min(x, na.rm = T)
  maxX = max(x, na.rm = T)
  
  n = length(x)
  degree = ord-1
  nrknots = nb - (degree-1) 
  d = (maxX-minX) / nrknots
  bm = c(minX - degree*d, maxX + degree*d)
  knots = seq(bm[1], bm[2], length = nrknots+2*degree)
  
  # numerical integration for normalizing the B-spline basis functions
  npoints = 10000
  xseq =  seq(bm[1], bm[2], length=npoints)
  B0 = splines::spline.des(knots, xseq, degree+1, outer.ok=T)$design # unnormalized
  w = rep(NA, nb)
  h = diff(c(knots[1], knots[length(knots)])) / npoints
  for (k in 1:nb){
    w[k] = (h* sum(B0[,k]))^(-1) 
    # this computes the integrals of the B-spline basis functions (which are then standardized below)
  } 
  # actual data design matrix
  
  B = matrix(NA, nrow = n, ncol = nb)
  ind = which(!is.na(x))
  B[ind,] = t(t(splines::spline.des(knots, x[ind], degree+1, outer.ok = TRUE)$design) * w) 
  
  list(B=B, knots=knots, w=w)
}

K = 40
ord = 4
degree = ord-1

desmat = build_desmat(ltuco, K=K, ord = ord)
B=desmat$B
knots=desmat$knots
w=desmat$w

L = WH:::build_D_mat(K+1, 2)
D = t(L[,-1])%*%L[,-1]

b.pos = knots[(degree+1):(length(knots)-degree+1)]
b0 = log(cbind(dnorm(b.pos+0.5, mu[1], sigma[1]),
               dnorm(b.pos+0.5, mu[2], sigma[2]),
               dnorm(b.pos+0.5, mu[3], sigma[3])))
b0 = b0 - matrix(apply(b0, 2, min), nrow = length(b.pos), ncol = N, byrow = TRUE)

theta.star = c(mod_sim_tuco$estimate[1:(N*(N-1))], as.numeric(b0))
               

# iterative procedure

maxiter = 50 # maximum number of iterations
tol = 0.05 # relative tolerance for convergence, sufficient
gradtol = 1e-6 # relative gradient tolerance for nlm (1e-6 is the default)

Lambdas = matrix(NA, maxiter+1, N)
Lambdas[1,] = c(300,300,200)
mods = list()

p = N*(N-1)
REind = p+matrix(1:(K*N), nrow = N, byrow = TRUE) # each row is the index of one RE

# updating algorithm
for(k in 1:maxiter){
  cat("\n\n- Iteration", k, "-")
  if(k == 1){
    print.level = 2
  } else print.level = 0
  
  t1 = Sys.time()
  mods[[k]] = nlm(mllk_tuco, theta.star, N=3, x=ltuco, B=B, lambda = Lambdas[k,], D=D, 
            iterlim = 1000, print.level = print.level, hessian = TRUE, gradtol = gradtol)
  Sys.time()-t1
  cat("\nEstimation time:", Sys.time()-t1)
  cat("\nIterations:", mods[[k]]$iterations)
  
  theta.star = mods[[k]]$estimate # saves theta.star as starting value for next iteration

  J_p = mods[[k]]$hessian # assigning hessian
  J_inv = MASS::ginv(J_p)
  
  # updating all penalty strengths
  for(i in 1:N){
    edoF = sum(diag( # trace
      diag(rep(1, nrow(D))) - Lambdas[k, i] * J_inv[REind[i,], REind[i,]] %*% D
    ))
    penalty = t(theta.star[REind[i,]]) %*% D %*% theta.star[REind[i,]]
    
    Lambdas[k+1, i] = as.numeric((edoF - 1) / penalty)
  }
  
  cat("\nSmoothing strengths:", round(Lambdas[k+1, ], 4))
  
  if(mean(abs(Lambdas[k+1,] - Lambdas[k,]) / Lambdas[k,]) < tol){
    cat("\n\n")
    break
  }
}

Lambdas = as.matrix(na.omit(Lambdas))
par(mfrow = c(1,3))
for(j in 1:N){
  plot(Lambdas[,j], type = "l", lwd = 2, main = paste("lambda", j), ylim = c(0,20))
}


t1 = Sys.time()
mod_final = nlm(mllk_tuco, theta.star, N=3, x=ltuco, B=B, lambda = Lambdas[nrow(Lambdas),], D=D, 
                iterlim = 1000, print.level = 2, hessian = TRUE)
Sys.time()-t1
cat("\nEstimation time:", Sys.time()-t1)

theta.star = mod_final$estimate

# length(mods)
# theta.star = mods[[2]]$estimate

# plot results

Gamma = LaMa::tpm(theta.star[1:(N*(N-1))]); p = N*(N-1)
delta = LaMa::stationary(Gamma)

nb = ncol(B) # number of basis functions
b = theta.star[p + 1:((nb-1)*N)]
Beta = matrix(0, nrow = N, ncol = nb) # spline coefficent matrix
Beta[,-1] = matrix(b, nrow = N, ncol = nb-1, byrow = TRUE)
A = exp(Beta)
A = A / rowSums(A)

# starting values should be chosen such that model resembles simple Gaussian model
color = c("orange", "deepskyblue", "seagreen2")
plotseq = seq(min(ltuco, na.rm=T), max(ltuco, na.rm=T), length=500)
Bplot = splines::spline.des(knots, plotseq, degree+1, outer.ok=T)$design * w[1]

pdf("./case_studies/figs/tuco_marginal.pdf", width = 5.5, height = 4)
par(mfrow = c(1,1))
hist(ltuco, breaks = 50, prob = TRUE, bor = "white", 
     xlab = "logVEDBA", ylab = "density", main = "")
dens = matrix(NA, nrow = length(plotseq), ncol = N)
for(j in 1:N){
  dens[,j] = delta[j] * rowSums(Bplot * matrix(A[j,], nrow=length(plotseq), ncol=nb, byrow=T))
  lines(plotseq,
        dens[,j], col = color[j], lwd = 2)
  # for(k in 1:nb){
  #   lines(plotseq, Bplot[,k] * A[j,k], col = color[j], lwd = 1)
  # }
}
lines(plotseq, rowSums(dens), col = 1, lty = 2, lwd = 2)
legend("topleft", legend = c("state 1", "state 2", "state 3", "marginal"), 
       col = c(color, 1), lty = c(1,1,1,2), lwd = c(2,2,2,2), bty = "n")
dev.off()


# plot all models
length(mods)
par(mfrow = c(2,2))
for(m in 1:length(mods)){
  theta.star = mods[[m]]$estimate
  Gamma = LaMa::tpm(theta.star[1:(N*(N-1))]); p = N*(N-1)
  delta = LaMa::stationary(Gamma)
  nb = ncol(B) # number of basis functions
  b = theta.star[p + 1:((nb-1)*N)]
  Beta = matrix(0, nrow = N, ncol = nb) # spline coefficent matrix
  Beta[,-1] = matrix(b, nrow = N, ncol = nb-1, byrow = TRUE)
  A = exp(Beta)
  A = A / rowSums(A)
  
  color = c("orange", "deepskyblue", "seagreen2")
  plotseq = seq(min(ltuco, na.rm=T), max(ltuco, na.rm=T), length=500)
  Bplot = splines::spline.des(knots, plotseq, degree+1, outer.ok=T)$design * w[1]
  
  hist(ltuco, breaks = 50, prob = TRUE, bor = "white", xlab = "logVEDBA", main = "")
  dens = matrix(NA, nrow = length(plotseq), ncol = N)
  for(j in 1:N){
    dens[,j] = delta[j] * rowSums(Bplot * matrix(A[j,], nrow=length(plotseq), ncol=nb, byrow=T))
    lines(plotseq,
          dens[,j], col = color[j], lwd = 2)
  }
  lines(plotseq, rowSums(dens), col = 1, lty = 2, lwd = 2)
}


## state decoding
allprobs = matrix(1, nrow = length(ltuco), ncol = N)
ind = which(!is.na(ltuco))
for (j in 1:N){
  allprobs[ind, j] = rowSums(B[ind,] * matrix(A[j,], nrow=length(ind), ncol=nb, byrow=T))
}

states = LaMa::viterbi(delta, Gamma, allprobs)

par(mfrow = c(1,1))
plotind = 1:2000
plot(ltuco[plotind], pch = 20, col = color[states[plotind]], bty = "n")
