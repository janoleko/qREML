## probably acceleration data best

data = read.csv("http://www.rolandlangrock.com/whitetip_data.csv")

data$lODBA = log(data$ODBA)

hist(data$lODBA, breaks = 100, prob = TRUE)


mllk_sim = function(theta.star, X, N){
  Gamma = LaMa::tpm(theta.star[1:(N*(N-1))])
  delta = LaMa::stationary(Gamma)
  p = N*(N-1)
  mu = theta.star[p + 1:N]
  sigma = exp(theta.star[p + N + 1:N])
  
  allprobs = matrix(1, nrow = nrow(X), ncol = N)
  ind = which(!is.na(X$lODBA))
  for (j in 1:N){
    allprobs[ind, j] = dnorm(X$lODBA[ind], mu[j], sigma[j])
  }
  -LaMa::forward(delta, Gamma, allprobs)
}

# first fit simple model
theta.star = c(rep(-2, 6),
               c(-4, -3.5, -2.5),
               log(c(1,1,1)))

mod_sim = nlm(mllk_sim, theta.star, X = data, N = 3)
theta.star = mod_sim$estimate
N=3

Gamma = LaMa::tpm(theta.star[1:(N*(N-1))]); p = N*(N-1)
delta = LaMa::stationary(Gamma)
mu = theta.star[p + 1:N]
sigma = exp(theta.star[p + N + 1:N])


# flexible model via marginal ML

# building B-spline design matrix
ord = 4
K = 20
nb = 2*K+1
minX = min(data$lODBA, na.rm = T)
maxX = max(data$lODBA, na.rm = T)
bm = c(minX - 0.05 * (maxX - minX), maxX + 0.05 * (maxX - minX))

n = nrow(data)
degree = ord-1
nrknots = nb - (ord-2) 
knots = seq(bm[1], bm[2], length = nrknots+2*degree)
npoints = 10000
xseq =  seq(bm[1], bm[2], length=npoints)
B0 = splines::spline.des(knots, xseq, degree+1, outer.ok=T)$design # unnormalized
w = rep(NA, nb)
h = diff(c(knots[1], knots[length(knots)]))/npoints
for (k in 1:nb){
  w[k] = (h* sum(B0[,k]))^(-1) 
  # this computes the integrals of the B-spline basis functions (which are then standardized below)
} 
# actual data design matrix
B = t(t(splines::spline.des(knots, data$lODBA, degree+1, outer.ok=T)$design) * w) 


mllk = function(theta.star, X, N, B, lambda, D){
  Gamma = LaMa::tpm(theta.star[1:(N*(N-1))]); p = N*(N-1)
  delta = LaMa::stationary(Gamma)
  
  nb = ncol(B)-1
  b = numeric(N*(nb+1))
  b[-(0:(N-1)*(nb+1) + 1)] = theta.star[p + 1:(nb*N)]
  A = matrix(exp(b), nrow = N, ncol = nb+1, byrow = TRUE)
  A = A / rowSums(A)
  
  allprobs = matrix(1, nrow = nrow(X), ncol = N)
  ind = which(!is.na(X$lODBA))
  for (j in 1:N){
    allprobs[ind, j] = rowSums(B[ind,] * matrix(A[j,], nrow=length(ind), ncol=nb+1, byrow=T))
  }
  
  pen = t(b[-(0:(N-1)*(nb+1) + 1)]) %*% kronecker(diag(lambda), D) %*% b[-(0:(N-1)*(nb+1) + 1)]
  
  -LaMa::forward(delta, Gamma, allprobs) + 0.5 * pen
}


diff_order = 2
L = WH:::build_D_mat(2*K, diff_order)
D = t(L)%*%L

b.pos = knots[(degree+1):(length(knots)-degree+1)]
theta.star = c(mod_sim$estimate[1:(N*(N-1))],
               log(100*c(dnorm(b.pos, mu[1], sigma[1]),
                        dnorm(b.pos, mu[2], sigma[2]),
                        dnorm(b.pos, mu[3], sigma[3]))))

# iterative procedure

maxiter = 50 # maximum number of iterations
tol = 0.02 # relative tolerance for convergence, sufficient
gradtol = 1e-5 # relative gradient tolerance for nlm (1e-6 is the default)
print.level = 0 # print level for nlm

Lambdas = matrix(NA, maxiter+1, N)
Lambdas[1,] = c(150, 20, 80)
mods = list()

p = N*(N-1)
REind = matrix(1:(2*K*N), nrow = N, byrow = TRUE) # each row is the index of one RE

# updating algorithm
for(k in 1:maxiter){
  cat("\n\n- Iteration", k, "-")
  if(k == 1) print.level = 2
  else print.level = 0
  
  # cat("\nFitting model...")
  t1 = Sys.time()
  mod = nlm(mllk, theta.star, N=3, X=data, B=B, lambda = Lambdas[k,], D=D, 
            iterlim = 1000, print.level = print.level, hessian = TRUE, gradtol = gradtol)
  Sys.time()-t1
  cat("\nEstimation time:", Sys.time()-t1)
  cat("\nIterations:", mod$iterations)
  
  theta.star = mod$estimate # saves theta.star as starting value for next iteration
  mods[[k]] = mod # saving model
  
  J_p = mod$hessian # assigning hessian
  b = theta.star[p + 1:(2*K*(N))]
  
  # updating all penalty strengths
  for(i in 1:N){
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
par(mfrow = c(1,3))
plot(Lambdas[,1], type = "l", lwd = 2, main = "lambda1")
plot(Lambdas[,2], type = "l", lwd = 2, main = "lambda2")
plot(Lambdas[,3], type = "l", lwd = 2, main = "lambda3")


t1 = Sys.time()
mod_final = nlm(mllk, theta.star, N=3, X=data, B=B, lambda = Lambdas[nrow(Lambdas),], D=D, 
          iterlim = 1000, print.level = 2, hessian = TRUE)
Sys.time()-t1
cat("\nEstimation time:", Sys.time()-t1)

theta.star = mod_final$estimate


# plot results

Gamma = LaMa::tpm(theta.star[1:(N*(N-1))]); p = N*(N-1)
delta = LaMa::stationary(Gamma)

nb = 2*K
b = numeric(N*(nb+1))
b[-(0:(N-1)*(nb+1) + 1)] = theta.star[p + 1:(nb*N)]
A = matrix(exp(b), nrow = N, ncol = nb+1, byrow = TRUE)
A = A / rowSums(A)

# starting values should be chosen such that model resembles simple Gaussian model
color = c("orange", "deepskyblue", "seagreen2")
plotseq = seq(bm[1], bm[2], length=500)
Bplot = splines::spline.des(knots, plotseq, degree+1, outer.ok=T)$design # unnormalized

par(mfrow = c(1,1))
hist(data$lODBA, breaks = 80, prob = TRUE, bor = "white", xlim = c(-5, -1))
dens = matrix(NA, nrow = length(plotseq), ncol = N)
for(j in 1:N){
  dens[,j] = delta[j]*rowSums(Bplot * w[1] * matrix(A[j,], nrow=length(plotseq), ncol=nb+1, byrow=T))
  lines(plotseq,
        dens[,j], col = color[j], lwd = 2)
  for(k in 1:(nb+1)){
    lines(plotseq, Bplot[,k] * A[j,k], col = color[j], lwd = 1)
  }
}
lines(xseq, rowSums(dens), col = 1, lty = 2, lwd = 2)


