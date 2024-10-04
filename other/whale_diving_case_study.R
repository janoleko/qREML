# Movebank Data Repository: https://doi.org/10.5441/001/1.47h576f2
data = read.csv("~/Downloads/blue_dives.csv")

colnames(data)

unique(data$animal.id)

data1 = data[which(data$animal.id==unique(data$animal.id)[1]),]
data2 = data[which(data$animal.id==unique(data$animal.id)[2]),]

datai = data[which(data$animal.id==unique(data$animal.id)[6]),]
hist(log(datai$MaxDiveDepth), breaks = 20)

plot(datai$MaxDiveDepth, type = "h")

data2 = datai


# potential candidates: 3, 6

nrow(data1)
nrow(data2)


plot(data2$MaxDiveDepth, type = "h")

hist(data2$MaxDiveDepth, prob = T, breaks = 25)

hist(log(data2$MaxDiveDepth), prob = T, breaks = 25)
curve(0.5*dnorm(x, 2.8, 0.3), add = T)
curve(0.4*dnorm(x, 4.1, 0.4), add = T)
curve(0.1*dnorm(x, 5.4, 0.4), add = T)


mllk = function(theta.star, N, x) {
  mu = theta.star[1:N]
  sigma = exp(theta.star[N + 1:N])
  
  Gamma = LaMa::tpm(theta.star[2*N + 1:(N*(N-1))])
  delta = LaMa::stationary(Gamma)
  
  allprobs = matrix(1, length(x), N)
  ind = which(!is.na(x))
  for(j in 1:N) allprobs[ind, j] = dnorm(x[ind], mu[j], sigma[j])
  
  -LaMa::forward(delta, Gamma, allprobs)
}
x = log(data2$MaxDiveDepth)

theta.star = c(2.8, 4.1, 5.4, log(c(0.3, 0.4, 0.4)),
               rep(-3, 6))
mod1 = nlm(mllk, theta.star, N = 3, x = log(data2$MaxDiveDepth), print.level = 2)

theta.star = mod1$estimate

N = 3
mu = theta.star[1:N]
sigma = exp(theta.star[N + 1:N])

Gamma = LaMa::tpm(theta.star[2*N + 1:(N*(N-1))])
delta = LaMa::stationary(Gamma)

allprobs = matrix(1, nrow(data2), N)
ind = which(!is.na(log(data2$MaxDiveDepth)))
for(j in 1:N) allprobs[ind, j] = dnorm(log(data2$MaxDiveDepth)[ind], mu[j], sigma[j])

color = c("orange", "deepskyblue", "seagreen2", "plum")
hist(log(data2$MaxDiveDepth), prob = T, breaks = 30, bor = "white")
for(j in 1:N) curve(delta[j] * dnorm(x, mu[j], sigma[j]), add = T, col = color[j], lwd = 2)
curve(delta[1]*dnorm(x, mu[1], sigma[1])+
        delta[2]*dnorm(x, mu[2], sigma[2])+
        delta[3]*dnorm(x, mu[3], sigma[3]), add = T, lwd = 2, lty = 2)

states = LaMa::viterbi(delta, Gamma, allprobs)

plot(log(data2$MaxDiveDepth), type = "h", col = color[states], lwd = 1)


## 4 states
theta.star = c(2.2, 3.5, 4, 5, log(c(0.25, 0.5,0.5, 1)),
               rep(-4, 12))
mod2 = nlm(mllk, theta.star, N = 4, x = log(data2$MaxDiveDepth), 
           print.level = 2, iterlim = 1000)

theta.star = mod2$estimate

N = 4

mu = theta.star[1:N]
sigma = exp(theta.star[N + 1:N])

Gamma = LaMa::tpm(theta.star[2*N + 1:(N*(N-1))])
delta = LaMa::stationary(Gamma)

allprobs = matrix(1, nrow(data2), N)
ind = which(!is.na(log(data2$MaxDiveDepth)))
for(j in 1:N) allprobs[ind, j] = dnorm(log(data2$MaxDiveDepth)[ind], mu[j], sigma[j])

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
lines(plotseq, rowSums(dens), col = 1, lty = 2, lwd = 2)for(j in 1:N) curve(delta[j] * dnorm(x, mu[j], sigma[j]), add = T, col = color[j], lwd = 2)
curve(delta[1]*dnorm(x, mu[1], sigma[1])+
        delta[2]*dnorm(x, mu[2], sigma[2])+
        delta[3]*dnorm(x, mu[3], sigma[3])+
        delta[4]*dnorm(x, mu[4], sigma[4]), add = T, lwd = 2, lty = 2)

states2 = LaMa::viterbi(delta, Gamma, allprobs)

plot(log(data2$MaxDiveDepth), type = "h", col = color[states2], lwd = 1)



# Non-parametric model ----------------------------------------------------

mllk_np = function(theta.star, x, N, B, lambda, S){
  Gamma = LaMa::tpm(theta.star[1:(N*(N-1))])
  delta = LaMa::stationary(Gamma)
  
  nb = ncol(B) # number of basis functions
  b = theta.star[N*(N-1) + 1:((nb-1)*N)]
  Beta = matrix(0, nrow = N, ncol = nb) # spline coefficent matrix
  Beta[,-1] = matrix(b, nrow = N, ncol = nb-1, byrow = TRUE)
  A = exp(Beta)
  A = A / rowSums(A)
  
  allprobs = matrix(1, nrow = length(x), ncol = N)
  ind = which(!is.na(x))
  for (j in 1:N){
    allprobs[ind, j] = rowSums(B[ind,] * matrix(A[j,], nrow=length(ind), ncol=nb, byrow=T))
  }
  
  pen = t(b) %*% kronecker(diag(lambda), S) %*% b
  
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

K = 20
ord = 4
degree = ord-1

desmat = build_desmat(log(data2$MaxDiveDepth), K=K, ord = ord)
B=desmat$B
knots=desmat$knots
w=desmat$w

L = WH:::build_D_mat(K+1, 2)
S = t(L[,-1])%*%L[,-1]

m = nrow(S) - as.numeric(Matrix::rankMatrix(S))

b.pos = knots[(degree+1):(length(knots)-degree+1)]
b0 = log(cbind(dnorm(b.pos, mu[1], sigma[1]),
               dnorm(b.pos, mu[2], sigma[2]),
               dnorm(b.pos, mu[3], sigma[3])))
b0 = b0 - matrix(apply(b0, 2, min), nrow = length(b.pos), ncol = N, byrow = TRUE)

theta.star = c(rep(-2, N*(N-1)), as.numeric(b0))


# iterative procedure

maxiter = 50 # maximum number of iterations
tol = 0.01 # relative tolerance for convergence, sufficient
gradtol = 1e-6 # relative gradient tolerance for nlm (1e-6 is the default)

Lambdas = matrix(NA, maxiter+1, N)
Lambdas[1,] = rep(100, 3)
mods = list()

REind = N*(N-1) + matrix(1:(K*N), nrow = N, byrow = TRUE) # each row is the index of one RE

# updating algorithm
for(k in 1:maxiter){
  cat("\n\n- Iteration", k, "-")
  if(k == 1){
    print.level = 2
  } else print.level = 0
  
  t1 = Sys.time()
  mods[[k]] = nlm(mllk_np, theta.star, N=3, x=log(data2$MaxDiveDepth), B=B, lambda = Lambdas[k,], S=S, 
                  iterlim = 1000, print.level = print.level, hessian = TRUE, gradtol = gradtol, stepmax = 80)
  Sys.time()-t1
  cat("\nEstimation time:", Sys.time()-t1)
  cat("\nIterations:", mods[[k]]$iterations)
  
  theta.star = mods[[k]]$estimate # saves theta.star as starting value for next iteration
  J_inv = MASS::ginv(mods[[k]]$hessian)
  
  # updating all penalty strengths
  for(i in 1:N){
    edoF = sum(diag( # trace
      diag(rep(1, nrow(S))) - Lambdas[k, i] * J_inv[REind[i,], REind[i,]] %*% S
    ))
    penalty = t(theta.star[REind[i,]]) %*% S %*% theta.star[REind[i,]]
    
    Lambdas[k+1, i] = as.numeric((edoF - m) / penalty)
  }
  Lambdas[k+1, which(Lambdas[k+1,] < 0)] = 0
  
  cat("\nSmoothing strengths:", round(Lambdas[k+1, ], 4))
  
  if(mean(abs(Lambdas[k+1,] - Lambdas[k,]) / (Lambdas[k,] + 1e-10)) < tol){
    cat("\n\n")
    break
  }
}


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
ldepth = log(data2$MaxDiveDepth)
plotseq = seq(min(ldepth, na.rm=T), max(ldepth, na.rm=T), length=500)
Bplot = splines::spline.des(knots, plotseq, degree+1, outer.ok=T)$design * w[1]

par(mfrow = c(1,1))
hist(log(data2$MaxDiveDepth), prob = T, breaks = 50, bor = "white")
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

# decoding
allprobs = matrix(1, nrow = length(ldepth), ncol = N)
ind = which(!is.na(ldepth))
for (j in 1:N){
  allprobs[ind, j] = rowSums(B[ind,] * matrix(A[j,], nrow=length(ind), ncol=nb, byrow=T))
}

states = LaMa::viterbi(delta, Gamma, allprobs)

plot(exp(ldepth), type = "h", col = color[states], lwd = 1)
