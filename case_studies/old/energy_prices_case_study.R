
# Case study: Spanish energy prices ---------------------------------------

## packages

# install.packages("MSwM")
library(MSwM)
# install.packages("LaMa")
library(LaMa)
# install.packages("mgcv")
library(mgcv)


# Loading the data --------------------------------------------------------

data(energy, package = "MSwM")


# EDA ---------------------------------------------------------------------

nrow(energy)
head(energy)


# Model fitting -----------------------------------------------------------

# likelihood function
mllk_mgcv = function(theta.star, X, Z, lambda, S){
  Gamma = LaMa::tpm(theta.star[1:2])
  delta = c(1, exp(theta.star[3]))
  delta = delta / sum(delta)
  
  nb = ncol(Z)
  beta = matrix(theta.star[3 + 1:(2*(nb+1))], ncol = 2) # parameter matrix for mean
  alpha = matrix(theta.star[3 + 2*(nb+1) + 1:(2*(nb+1))], ncol = 2) # parameter matrix for sd
  
  Mu = cbind(1, Z) %*% beta
  Sigma = exp(cbind(1, Z) %*% alpha)
  
  allprobs = cbind(dnorm(X$Price, Mu[,1], Sigma[,1]), dnorm(X$Price, Mu[,2], Sigma[,2]))
  
  # excluding the intercept, stacking into one vector, and penalising each RE separately by using kronecker product
  pen = t(as.numeric(beta[-1,])) %*% kronecker(diag(lambda[1:2]), S) %*% as.numeric(beta[-1,]) + # sum of penalties for mean
    t(as.numeric(alpha[-1,])) %*% kronecker(diag(lambda[3:4]), S) %*% as.numeric(alpha[-1,]) # sum of penalties for sd

  -LaMa::forward(delta, Gamma, allprobs) + 0.5 * pen # the 0.5 is important for the math to work!
}

# design and penalty matrix
# nb = 10
# gam_setup = gam(Price ~ s(Oil, k = nb, bs = "cs"), data = cbind(dummy = 1, energy), fit = FALSE)
# Z = gam_setup$X[,-1]
# S = gam_setup$S[[1]]

# make my own xdesign and penalty matrix -> works better
nb = 10 # number of basis functions, one will be excluded because b_0 = 0
d = diff(seq(min(energy$Oil), max(energy$Oil), length = nb-2))[1]
knots = seq(min(energy$Oil)-3*d, max(energy$Oil)+3*d, by = d)
Z = splines::splineDesign(knots, energy$Oil, outer.ok = TRUE)[,-1] # dropping first column because b_0 = 0

# penalty matrix
L = diff(diag(nb), differences = 2) # second-order difference matrix
S = t(L[,-1])%*%L[,-1] # dropping first column because b_0 = 0
m = nrow(S) - Matrix::rankMatrix(S)[1] # checking the rank deficiency for correction in algorithm (in this case 1)

# Eigen decomposition
# eigenS = eigen(S)
# U = eigenS$vectors
# La = diag(eigenS$values)

# Optimisation via marginal ML --------------------------------------------

maxiter = 50 # maximum number of iterations
tol = 0.01 # relative tolerance for convergence, sufficient 
Lambdas = matrix(NA, maxiter+1, 4)
Lambdas[1,] = c(rep(100000, 4)) # initial penalty strengths
mods = list()
# object that contains the indices for each random effect/ spline coef vector. If these are of differnt lengths, use list and change below
REind = matrix(3 + 1:(4*(nb)), nrow = 4, byrow = TRUE)[,-1] # each row is the index of one RE
alpha_sm = 1

# initial parameter
theta.star = c(rep(-4, 2), 1,
               # 3, rep(0, nb-1), 6, rep(0, nb-1),
               4.5, seq(-2, 2, length = nb-1), 6.5, seq(-1, 1, length = nb-1),
               rep(0, 2*nb))

# updating algorithm
T1 = Sys.time()
for(k in 1:maxiter){
  cat("\n\n- Iteration", k, "-")
  t1 = Sys.time()
  mods[[k]] = nlm(mllk_mgcv, theta.star, X = energy, Z=Z, lambda = Lambdas[k,], S=S, 
                  iterlim = 1000, print.level = 0, hessian = TRUE) # fitting the model with current penalty strengths
  Sys.time()-t1
  cat("\nEstimation time:", Sys.time()-t1)
  cat("\nIterations:", mods[[k]]$iterations)
  theta.star = mods[[k]]$estimate # saves theta.star as starting value for next iteration
  J_inv = MASS::ginv(mods[[k]]$hessian) # inverting hessian
  # updating all penalty strengths
  for(i in 1:4){
    # computing the effective degrees of freedom
    edoF = nrow(S) - sum(diag(Lambdas[k, i] * J_inv[REind[i,], REind[i,]] %*% S))
    penalty = t(theta.star[REind[i,]]) %*% S %*% theta.star[REind[i,]]
    # updating the penalty strength
    lambda_new = as.numeric((edoF - m) / penalty) # m is correction if S does not have full rank
    Lambdas[k+1, i] = alpha_sm * lambda_new + (1-alpha_sm) * Lambdas[k, i]
  }
  cat("\nSmoothing strengths:", round(Lambdas[k+1, ], 4))
  if(max(abs(Lambdas[k+1,] - Lambdas[k,]) / Lambdas[k,]) < tol){
    cat("\n\n")
    break
  }
}
Sys.time()-T1

Lambdas = as.matrix(na.omit(Lambdas))
lambda_hat = Lambdas[nrow(Lambdas),]
par(mfrow = c(1,4))
for(j in 1:4){
  plot(Lambdas[,j], type = "l", lwd = 2, main = paste("lambda", j), ylim = c(0, 10))
}
round(lambda_hat, 3)

theta.star = mods[[length(mods)]]$estimate

Gamma = tpm(theta.star[1:2])
delta = c(1, exp(theta.star[3]))
delta = delta / sum(delta)
beta = matrix(theta.star[3 + 1:(2*nb)], ncol = 2)
alpha = matrix(theta.star[3 + 2*nb + 1:(2*nb)], ncol = 2)
Mu = cbind(1, Z) %*% beta
Sigma = exp(cbind(1, Z) %*% alpha)
allprobs = cbind(dnorm(energy$Price, Mu[,1], Sigma[,1]),
                 dnorm(energy$Price, Mu[,2], Sigma[,2]))

states = viterbi(delta, Gamma, allprobs) # decoding most probable state sequence

xseq = seq(min(energy$Oil), max(energy$Oil), length = 200)
Z_plot = cbind(1, splines::splineDesign(knots, xseq, outer.ok = TRUE)[,-1])

# Generate the new design matrix using the PredictMat function
# Z_plot = cbind(1, PredictMat(gam_setup$smooth[[1]], data.frame(dummy = 1, Price = 1, Oil = xseq)))

Mu_plot = Z_plot %*% beta
Sigma_plot = exp(Z_plot %*% alpha)

color = c("orange", "deepskyblue")

# pdf("./case_studies/figs/energy_oil.pdf", width = 8, height = 4.5)

par(mfrow = c(1,2), mar = c(5,4,3,1))
plot(energy$Oil, energy$Price, pch = 20, bty = "n", col = scales::alpha(color[states], 0.1),
     xlab = "oil price", ylab = "energy price")
for(j in 1:2) lines(xseq, Mu_plot[,j], col = color[j], lwd = 3)

qseq = qnorm(seq(0.5, 0.95, length = 4))
for(i in qseq){
  for(j in 1:2){
    lines(xseq, Mu_plot[,j] + i*Sigma_plot[,j], col = scales::alpha(color[j], 0.7), lwd = 1, lty = 2)
    lines(xseq, Mu_plot[,j] - i*Sigma_plot[,j], col = scales::alpha(color[j], 0.7), lwd = 1, lty = 2)
  }
}

legend("topright", bty = "n", legend = paste("state", 1:2), col = color, lwd = 3)

plot(NA, xlim = c(0, nrow(energy)), ylim = c(1,10), bty = "n",
     xlab = "time", ylab = "energy price")
segments(x0 = 1:(nrow(energy)-1), x1 = 2:nrow(energy),
         y0 = energy$Price[-nrow(energy)], y1 = energy$Price[-1], col = color[states[-1]], lwd = 0.5)

# dev.off()

# par(mfrow = c(1,1))
# plot(Sigma_plot[,1], ylim = c(0,1.5), type = "l", col = color[1], lwd = 2, bty = "n",
#      xlab = "z", ylab = "standard deviation")
# lines(Sigma_plot[,2], col = color[2], lwd = 2)



# Plotting the model sequence ---------------------------------------------

length(mods)

plotind = c(1,2,3,5,10,length(mods))

xseq = seq(min(energy$Oil), max(energy$Oil), length = 200)
Z_plot = cbind(1, splines::splineDesign(knots, xseq, outer.ok = TRUE)[,-1])

# pdf("./case_studies/figs/energy_oil_modseq.pdf", width = 8, height = 5.5)

par(mfrow = c(2,3), mar = c(5,4,3,1))

for(m in plotind){
  theta.star = mods[[m]]$estimate
  
  beta = matrix(theta.star[3 + 1:(2*(nb))], ncol = 2)
  alpha = matrix(theta.star[3 + 2*(nb) + 1:(2*(nb))], ncol = 2)
  Mu = cbind(1, Z) %*% beta
  Sigma = exp(cbind(1, Z) %*% alpha)
  allprobs = cbind(dnorm(energy$Price, Mu[,1], Sigma[,1]),
                   dnorm(energy$Price, Mu[,2], Sigma[,2]))
  states = LaMa::viterbi(delta, Gamma, allprobs)
  
  Mu_plot = Z_plot %*% beta
  Sigma_plot = exp(Z_plot %*% alpha)
  
  plot(energy$Oil, energy$Price, pch = 20, bty = "n", col = scales::alpha(color[states], 0.1),
       xlab = "oil price", ylab = "energy price", main = paste("Iteration", m))
  for(j in 1:2) lines(xseq, Mu_plot[,j], col = color[j], lwd = 2)
}

#dev.off()

