# 
# # loading the data of first bird
# data = read.csv("./data/caracara.csv")
# 
# 
# # downsampling to 0.1 Hz
# thin_ind = (1:floor(nrow(data) / 100)) * 100
# data = data[thin_ind,] 
# 
# # calculating log (base e) of VDBA
# data$logVDBA = log(data$VDBA)
# 
# write.csv(data, "./data/caracara_downsampled.csv", row.names = F)

data = read.csv("./data/caracara_downsampled.csv")


# EDA ---------------------------------------------------------------------

head(data)



# Colors ------------------------------------------------------------------

color = c("orange", "deepskyblue", "seagreen4", "plum")


# Simple parametric HMM ---------------------------------------------------

# likelihood function
mllk_simple = function(theta.star, N, x) {
  mu = theta.star[1:N]
  sigma = exp(theta.star[N + 1:N])
  Gamma = LaMa::tpm(theta.star[2*N + 1:(N*(N-1))])
  delta = LaMa::stationary(Gamma)
  allprobs = matrix(1, length(x), N)
  ind = which(!is.na(x))
  for(j in 1:N) allprobs[ind, j] = dnorm(x[ind], mu[j], sigma[j])
  -LaMa::forward(delta, Gamma, allprobs)
}

mu0 = c(-5, -4, -2)
sigma0 = c(0.3, 0.3, 0.3)

theta.star0 = c(mu0, log(sigma0), rep(-2, 6))

mod_simple3 = nlm(mllk_simple, theta.star0, N = 3, x = data$logVDBA, print.level = 2, stepmax = 10, iterlim = 1000)

theta.star = mod_simple3$estimate
N = 3
mu3 = theta.star[1:N]
sigma3 = exp(theta.star[N + 1:N])
Gamma3 = LaMa::tpm(theta.star[2*N + 1:(N*(N-1))])
delta3 = LaMa::stationary(Gamma3)

hist(data$logVDBA, breaks = 50, prob = T, bor = "white", main = "", xlab = "log(VeDBA)", ylab = "density")
for(j in 1:N) curve(delta3[j] * dnorm(x, mu3[j], sigma3[j]), add = T, col = color[j], lwd = 2, n = 200)
curve(delta3[1] * dnorm(x, mu3[1], sigma3[1]) + 
        delta3[2] * dnorm(x, mu3[2], sigma3[2]) + 
        delta3[3] * dnorm(x, mu3[3], sigma3[3]), add = T, col = "black", lty = 2, lwd = 2, n = 200)

AIC3 = 2 * mod_simple3$minimum + 2 * length(mod_simple3$estimate)
BIC3 = 2 * mod_simple3$minimum + log(nrow(data)) * length(mod_simple3$estimate)

# 4 states ----------------------------------------------------------------

mu0 = c(-5, -4.5, -3, -1.5)
sigma0 = c(0.3, 0.3, 0.3, 0.3)
theta.star0 = c(mu0, log(sigma0), rep(-4, 12))

mod_simple4 = nlm(mllk_simple, theta.star0, N = 4, x = data$logVDBA, print.level = 2, stepmax = 10, iterlim = 1000)

theta.star = mod_simple4$estimate
N = 4
mu4 = theta.star[1:N]
sigma4 = exp(theta.star[N + 1:N])
Gamma4 = LaMa::tpm(theta.star[2*N + 1:(N*(N-1))])
delta4 = LaMa::stationary(Gamma4)

AIC4 = 2 * mod_simple4$minimum + 2 * length(mod_simple4$estimate)
BIC4 = 2 * mod_simple4$minimum + log(nrow(data)) * length(mod_simple4$estimate)


# Visualising both simple models together ---------------------------------

# pdf("./case_studies/figs/caracaras_simple.pdf", width = 8.5, height = 4)

par(mfrow = c(1,2))

# 3 states
N = 3
hist(data$logVDBA, breaks = 50, prob = T, bor = "white", main = "", xlab = "log(VeDBA)", ylab = "density")
for(j in 1:N) curve(delta3[j] * dnorm(x, mu3[j], sigma3[j]), add = T, col = color[j], lwd = 2, n = 200)
curve(delta3[1] * dnorm(x, mu3[1], sigma3[1]) + 
        delta3[2] * dnorm(x, mu3[2], sigma3[2]) + 
        delta3[3] * dnorm(x, mu3[3], sigma3[3]), add = T, col = "black", lty = 2, lwd = 2, n = 200)

legend("topright", legend = c("resting", "low activity", "high activity"), 
       col = color, lwd = 2, bty = "n")

# 4 states
N = 4
color2 = c("#FFEF00", "orange", "deepskyblue", "seagreen4")
hist(data$logVDBA, breaks = 50, prob = T, bor = "white", main = "", xlab = "log(VeDBA)", ylab = "density")
for(j in 2:N) curve(delta4[j] * dnorm(x, mu4[j], sigma4[j]), add = T, col = color2[j], lwd = 2, n = 200)
curve(delta4[1] * dnorm(x, mu4[1], sigma4[1]), add = T, col = color2[1], lwd = 2, n = 200)
curve(delta4[1] * dnorm(x, mu4[1], sigma4[1]) + 
        delta4[2] * dnorm(x, mu4[2], sigma4[2]) + 
        delta4[3] * dnorm(x, mu4[3], sigma4[3])+
        delta4[4] * dnorm(x, mu4[4], sigma4[4]), add = T, col = "black", lty = 2, lwd = 2, n = 200)

legend("topright", legend = c("resting 1", "resting 2", "low activity", "high activity"), 
       col = color2, lwd = 2, bty = "n")

# dev.off()

## corresponding lognormal distributions - for better comparison to original paper
par(mfrow = c(1,1))
color2 = c("#FFEF00", "orange", "deepskyblue", "seagreen4")
hist(data$VDBA, breaks = 150, prob = T, bor = "white", main = "", xlab = "log(VeDBA)", ylab = "density", xlim = c(0, 0.2), ylim = c(0,30))
for(j in 1:N) curve(delta4[j] * dlnorm(x, mu4[j], sigma4[j]), add = T, col = color2[j], lwd = 2, n = 700)


# state decoding
allprobs = matrix(1, length(data$logVDBA), N)
ind = which(!is.na(data$logVDBA))
for(j in 1:N) allprobs[ind, j] = dnorm(data$logVDBA[ind], mu[j], sigma[j])

states4 = LaMa::viterbi(delta, Gamma, allprobs)

plot(data$VDBA[1:1000], type = "h", col = color[states4], xlab = "time", ylab = "log(VeDBA")





# Gamma model -------------------------------------------------------------

# likelihood function
mllk_gamma = function(theta.star, N, x) {
  mu = exp(theta.star[1:N])
  sigma = exp(theta.star[N + 1:N])
  Gamma = LaMa::tpm(theta.star[2*N + 1:(N*(N-1))])
  delta = LaMa::stationary(Gamma)
  allprobs = matrix(1, length(x), N)
  ind = which(!is.na(x))
  for(j in 1:N) allprobs[ind, j] = dgamma(x[ind], shape=mu[j]^2/sigma[j]^2, scale=sigma[j]^2/mu[j])
  -LaMa::forward(delta, Gamma, allprobs)
}

# 3 states

mu0 = c(0.01, 0.05, 0.15)
sigma0 = c(0.01, 0.05, 0.15)

theta.star0 = c(log(mu0), log(sigma0), rep(-2, 6))

mod_gamma3 = nlm(mllk_gamma, theta.star0, N = 3, x = data$VDBA, print.level = 2, stepmax = 10, iterlim = 1000)

theta.star = mod_gamma3$estimate
N = 3
mu3_g = exp(theta.star[1:N])
sigma3_g = exp(theta.star[N + 1:N])
Gamma3_g = LaMa::tpm(theta.star[2*N + 1:(N*(N-1))])
delta3_g = LaMa::stationary(Gamma3_g)

N = 3
hist(data$VDBA, breaks = 150, prob = T, bor = "white", main = "", xlab = "log(VeDBA)", ylab = "density", xlim = c(0, 0.2), ylim = c(0,30))
for(j in 1:N) curve(delta3_g[j] * dgamma(x, shape=mu3_g[j]^2/sigma3_g[j]^2, scale=sigma3_g[j]^2/mu3_g[j]), add = T, col = color[j], lwd = 2, n = 500)
curve(delta3_g[1] * dgamma(x, shape=mu3_g[1]^2/sigma3_g[1]^2, scale=sigma3_g[1]^2/mu3_g[1])+ 
        delta3_g[2] * dgamma(x, shape=mu3_g[2]^2/sigma3_g[2]^2, scale=sigma3_g[2]^2/mu3_g[2])+ 
        delta3_g[3] * dgamma(x, shape=mu3_g[3]^2/sigma3_g[3]^2, scale=sigma3_g[3]^2/mu3_g[3]),
      add = T, col = "black", lty = 2, lwd = 2, n = 200)

legend("topright", legend = c("resting", "low activity", "high activity"), 
       col = color, lwd = 2, bty = "n")

# 4 states

mu0 = c(0.005, 0.01, 0.04, 0.3)
sigma0 = c(0.02, 0.01, 0.04, 0.2)

theta.star0 = c(log(mu0), log(sigma0), rep(-4, 12))

mod_gamma4 = nlm(mllk_gamma, theta.star0, N = 4, x = data$VDBA, print.level = 2, stepmax = 10, iterlim = 1000)

theta.star = mod_gamma4$estimate
N = 4
mu4_g = exp(theta.star[1:N])
sigma4_g = exp(theta.star[N + 1:N])
Gamma4_g = LaMa::tpm(theta.star[2*N + 1:(N*(N-1))])
delta4_g = LaMa::stationary(Gamma4_g)

# 3 states
N = 4
hist(data$VDBA, breaks = 150, prob = T, bor = "white", main = "", xlab = "log(VeDBA)", ylab = "density", xlim = c(0, 0.2), ylim = c(0,30))
for(j in 1:N) curve(delta4_g[j] * dgamma(x, shape=mu4_g[j]^2/sigma4_g[j]^2, scale=sigma4_g[j]^2/mu4_g[j]), add = T, col = color[j], lwd = 2, n = 500)
curve(delta4_g[1] * dgamma(x, shape=mu4_g[1]^2/sigma4_g[1]^2, scale=sigma4_g[1]^2/mu4_g[1])+ 
        delta4_g[2] * dgamma(x, shape=mu4_g[2]^2/sigma4_g[2]^2, scale=sigma4_g[2]^2/mu4_g[2])+ 
        delta4_g[3] * dgamma(x, shape=mu4_g[3]^2/sigma4_g[3]^2, scale=sigma4_g[3]^2/mu4_g[3]),
      add = T, col = "black", lty = 2, lwd = 2, n = 200)

legend("topright", legend = c("resting", "low activity", "high activity"), 
       col = color, lwd = 2, bty = "n")


###########################################################################
# Non-parametric model ----------------------------------------------------

# likelihood function
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
  allprobs[ind, ] = B[ind,] %*% t(A)
  
  pen = t(b) %*% kronecker(diag(lambda), S) %*% b
  
  -LaMa::forward(delta, Gamma, allprobs) + 0.5 * pen
}

# function that builds B-spline design matrix with basis functions integrating to one
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

# building the design matrix
K = 25
ord = 4
degree = ord-1
desmat = build_desmat(data$logVDBA, K = K, ord = ord)
B = desmat$B
knots = desmat$knots
w = desmat$w

# building the penalty matrix (penalizing squared second order differences)
L = WH:::build_D_mat(K+1, 2)
S = t(L[,-1])%*%L[,-1]
m = nrow(S) - as.numeric(Matrix::rankMatrix(S)) # rank deficiency (for correction)

# initial values for spline coefficients (from simple model)
N = 3
b.pos = knots[(degree+1):(length(knots)-degree+1)]
b0 = log(cbind(dnorm(b.pos, mu3[1], sigma3[1] * 1.3), # initially more overlap
               dnorm(b.pos, mu3[2], sigma3[2] * 1.3),
               dnorm(b.pos, mu3[3], sigma3[3] * 1.3)))
b0 = b0 - matrix(apply(b0, 2, min), nrow = length(b.pos), ncol = N, byrow = TRUE)

# initial parameter vector
# theta.star = c(mod_simple$estimate[1:(N*(N-1))], as.numeric(b0))
theta.star = c(rep(-3, N*(N-1)), as.numeric(b0))



# Model fitting via marginal ML -------------------------------------------

maxiter = 50 # maximum number of iterations
tol = 0.01 # relative tolerance for convergence, sufficient

Lambdas = matrix(NA, maxiter + 1, N)
Lambdas[1,] = c(10, 100, 100)
mods = list()

REind = N*(N-1) + matrix(1:(K*N), nrow = N, byrow = TRUE) # each row is the index of one RE

# updating algorithm
for(k in 1:maxiter){
  cat("\n\n- Iteration", k, "-")
  if(k == 1){
    print.level = 2
  } else print.level = 0
  
  t1 = Sys.time()
  mods[[k]] = nlm(mllk_np, theta.star, N=3, x=data$log10VDBA, B=B, lambda = Lambdas[k,], S=S, 
                  iterlim = 1000, print.level = print.level, hessian = TRUE, stepmax = 100)
  Sys.time()-t1
  cat("\nEstimation time:", Sys.time()-t1)
  cat("\nIterations:", mods[[k]]$iterations)
  
  theta.star = mods[[k]]$estimate # saves theta.star as starting value for next iteration
  J_inv = MASS::ginv(mods[[k]]$hessian)
  
  # updating all penalty strengths
  for(i in 1:N){
    edoF = sum(diag(diag(rep(1, nrow(S))) - Lambdas[k, i] * J_inv[REind[i,], REind[i,]] %*% S))
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

Lambdas = as.matrix(na.omit(Lambdas))
par(mfrow = c(1,N))
for(j in 1:N){
  plot(Lambdas[,j], type = "l", lwd = 2, main = paste("lambda", j))
}
lambda_hat = Lambdas[nrow(Lambdas),]


# theta.star = mods[[1]]$estimate
theta.star = mods[[length(mods)]]$estimate

# assigning parameters
N=3
Gamma = LaMa::tpm(theta.star[1:(N*(N-1))])
delta = LaMa::stationary(Gamma)
nb = ncol(B) # number of basis functions
b = theta.star[N*(N-1) + 1:((nb-1)*N)]
Beta = matrix(0, nrow = N, ncol = nb) # spline coefficent matrix
Beta[,-1] = matrix(b, nrow = N, ncol = nb-1, byrow = TRUE)
A = exp(Beta)
A = A / rowSums(A)

# AIC and BIC
edof = rep(NA, N)
for(i in 1:N){
  edof[i] = sum(diag(diag(rep(1, nrow(S))) - lambda_hat[i] * J_inv[REind[i,], REind[i,]] %*% S))
}

nllk = mllk_np(theta.star, x = data$logVDBA, N = 3, B = B, S = S, lambda = rep(0, N))
AIC = 2 * nllk + 2 * (N*(N-1) + sum(edof))
BIC = 2 * nllk + log(nrow(data)) * (N*(N-1) + sum(edof))

## comparing all models

AICBIC = cbind(m3 = c(AIC3, BIC3), m3 = c(AIC4, BIC4), m_np = c(AIC, BIC))
rownames(AICBIC) = c("AIC", "BIC")
round(AICBIC,1)

# decoding states
allprobs = matrix(1, nrow = length(data$logVDBA), ncol = N)
ind = which(!is.na(data$logVDBA))
allprobs[ind,] = B[ind,] %*% t(A)

states = LaMa::viterbi(delta, Gamma, allprobs)

# plotting state-dependent distributions
plotseq = seq(min(data$logVDBA, na.rm=T), max(data$logVDBA, na.rm=T), length = 500)
Bplot = splines::spline.des(knots, plotseq, degree+1, outer.ok=T)$design * w[1]

# pdf("./case_studies/figs/caracara.pdf", width = 8.5, height = 4)

par(mfrow = c(1,2), mar = c(5,4,3.5,1)+0.1)
hist(data$logVDBA, prob = T, breaks = 50, bor = "white",
     main = "", xlab = "log(VeDBA)", ylab = "density")
dens = matrix(NA, nrow = length(plotseq), ncol = N)
for(j in 1:N){
  dens[,j] = delta[j] * rowSums(Bplot * matrix(A[j,], nrow=length(plotseq), ncol=nb, byrow=T))
  lines(plotseq, dens[,j], col = color[j], lwd = 2)
}
lines(plotseq, rowSums(dens), col = 1, lty = 2, lwd = 2)

legend("topright", legend = c("resting", "low activity", "high activity"), 
       col = color, lwd = 2, bty = "n")

plotind = 4000:6000
plot(data$VDBA[plotind], type = "h", col = color[states[plotind]], 
     xlab = "time", ylab = "VeDBA", bty = "n", lwd = 0.3)

# dev.off()



# Non-parametric with 4 states --------------------------------------------

N = 4
# initial values for spline coefficients (from simple model)
b.pos = knots[(degree+1):(length(knots)-degree+1)]
b0 = log(cbind(dnorm(b.pos, mu[1], sigma[1] * 1.5), # initially more overlap
               dnorm(b.pos, mu[2], sigma[2] * 1.5),
               dnorm(b.pos, mu[3], sigma[3] * 1.5),
               dnorm(b.pos, mu[4], sigma[4] * 1.5)))
b0 = b0 - matrix(apply(b0, 2, min), nrow = length(b.pos), ncol = N, byrow = TRUE)

# initial parameter vector
# theta.star = c(mod_simple$estimate[1:(N*(N-1))], as.numeric(b0))
theta.star = c(rep(-4, N*(N-1)), as.numeric(b0))


# Model fitting via marginal ML

maxiter = 50 # maximum number of iterations
tol = 0.01 # relative tolerance for convergence, sufficient

Lambdas4 = matrix(NA, maxiter + 1, N)
Lambdas4[1,] = c(15, 15, 100, 100)
mods4 = list()

REind = N*(N-1) + matrix(1:(K*N), nrow = N, byrow = TRUE) # each row is the index of one RE

# updating algorithm
for(k in 1:maxiter){
  cat("\n\n- Iteration", k, "-")
  if(k == 1){
    print.level = 2
  } else print.level = 0
  
  t1 = Sys.time()
  mods4[[k]] = nlm(mllk_np, theta.star, N=4, x=data$log10VDBA, B=B, lambda = Lambdas4[k,], S=S, 
                  iterlim = 1000, print.level = print.level, hessian = TRUE, stepmax = 100)
  Sys.time()-t1
  cat("\nEstimation time:", Sys.time()-t1)
  cat("\nIterations:", mods4[[k]]$iterations)
  
  theta.star = mods4[[k]]$estimate # saves theta.star as starting value for next iteration
  J_inv = MASS::ginv(mods4[[k]]$hessian)
  
  # updating all penalty strengths
  for(i in 1:N){
    edoF = sum(diag(diag(rep(1, nrow(S))) - Lambdas4[k, i] * J_inv[REind[i,], REind[i,]] %*% S))
    penalty = t(theta.star[REind[i,]]) %*% S %*% theta.star[REind[i,]]
    Lambdas4[k+1, i] = as.numeric((edoF - m) / penalty)
  }
  Lambdas4[k+1, which(Lambdas4[k+1,] < 0)] = 0
  
  cat("\nSmoothing strengths:", round(Lambdas4[k+1, ], 4))
  if(mean(abs(Lambdas4[k+1,] - Lambdas4[k,]) / (Lambdas4[k,] + 1e-10)) < tol){
    cat("\n\n")
    break
  }
}

Lambdas4 = as.matrix(na.omit(Lambdas4))
par(mfrow = c(1,N))
for(j in 1:N){
  plot(Lambdas4[,j], type = "l", lwd = 2, main = paste("lambda", j))
}


## visualising results
N=4
theta.star = mods4[[length(mods)]]$estimate

# assigning parameters
Gamma = LaMa::tpm(theta.star[1:(N*(N-1))])
delta = LaMa::stationary(Gamma)
nb = ncol(B) # number of basis functions
b = theta.star[N*(N-1) + 1:((nb-1)*N)]
Beta = matrix(0, nrow = N, ncol = nb) # spline coefficent matrix
Beta[,-1] = matrix(b, nrow = N, ncol = nb-1, byrow = TRUE)
A = exp(Beta)
A = A / rowSums(A)

# decoding states
allprobs = matrix(1, nrow = length(data$logVDBA), ncol = N)
ind = which(!is.na(data$logVDBA))
for (j in 1:N) allprobs[ind, j] = B[ind,] %*% A[j,]

states4 = LaMa::viterbi(delta, Gamma, allprobs)

# plotting state-dependent distributions
plotseq = seq(min(data$logVDBA, na.rm=T), max(data$logVDBA, na.rm=T), length = 500)
Bplot = splines::spline.des(knots, plotseq, degree+1, outer.ok=T)$design * w[1]

# pdf("./case_studies/figs/caracaras_marginal.pdf", width = 6.5, height = 4.5)

par(mfrow = c(1,1))
hist(data$logVDBA, prob = T, breaks = 50, bor = "white",
     main = "", xlab = "log(VeDBA)", ylab = "density")
dens = matrix(NA, nrow = length(plotseq), ncol = N)
for(j in 1:N){
  dens[,j] = delta[j] * rowSums(Bplot * matrix(A[j,], nrow=length(plotseq), ncol=nb, byrow=T))
  lines(plotseq, dens[,j], col = color[j], lwd = 2)
}
lines(plotseq, rowSums(dens), col = 1, lty = 2, lwd = 2)

legend("topright", legend = c("resting", "low activity", "high activity", "high activity 2"), 
       col = color, lwd = 2, bty = "n")

# dev.off()

par(mfrow = c(1,1))

plot(data$logVDBA, pch = 20, col = color[states], xlab = "time", ylab = "log(VeDBA)")
plot(data$VDBA, type = "h", col = color[states4], xlab = "time", ylab = "VeDBA", bty = "n", ylim = c(0, 0.4))




L1 = cbind(diag(5),0)
diag(L1[,-1]) = -1

L2 = L1[-1,-1] %*% L1

S = t(L2)%*%L2

determinant(S, logarithm = TRUE)
