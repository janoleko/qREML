
# Case study: Spanish energy prices ---------------------------------------

## packages

# install.packages("MSwM")
library(MSwM)
# install.packages("LaMa")
library(LaMa)


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
  beta = matrix(theta.star[3 + 1:(2*(nb))], ncol = 2) # parameter matrix for mean
  alpha = matrix(theta.star[3 + 2*(nb) + 1:(2*(nb))], ncol = 2) # parameter matrix for sd
  
  Mu = Z %*% beta
  Sigma = exp(Z %*% alpha)
  
  allprobs = cbind(dnorm(X$Price, Mu[,1], Sigma[,1]), dnorm(X$Price, Mu[,2], Sigma[,2]))
  
  # excluding the intercept, stacking into one vector, and penalising each RE separately by using kronecker product
  pen = t(as.numeric(beta[-1,])) %*% kronecker(diag(lambda[1:2]), S) %*% as.numeric(beta[-1,]) + # sum of penalties for mean
    t(as.numeric(alpha[-1,])) %*% kronecker(diag(lambda[3:4]), S) %*% as.numeric(alpha[-1,]) # sum of penalties for sd

  -LaMa::forward(delta, Gamma, allprobs) + 0.5 * pen # the 0.5 is important for the math to work!
}

# design and penalty matrix
# nb = 12
# gam_setup = mgcv::gam(Price ~ s(Oil, k = nb, bs = "ps"), data = cbind(dummy = 1, energy), fit = FALSE)
# Z = gam_setup$X
# S = gam_setup$S[[1]]

# make my own design and penalty matrix -> works better
nb = 10
d = diff(seq(min(energy$Oil), max(energy$Oil), length = (nb-1)-2))[1]
knots = seq(min(energy$Oil)-3*d, max(energy$Oil)+3*d, by = d)
# x = seq(min(energy$Oil), max(energy$Oil), length = 200)
Z = splines::splineDesign(knots, energy$Oil, outer.ok = TRUE)
Z = cbind(1, Z)

# penalty matrix
L = WH:::build_D_mat(nb-1, 2) # builds second-order difference matrix
# S = t(L[,-1])%*%L[,-1]
S = t(L)%*%L # penalty matrix is L^t L
m = nrow(S) - Matrix::rankMatrix(S)[1] # checking the rank deficiency for correction in algorithm

# m = 0


# Optimisation via marginal ML --------------------------------------------

maxiter = 50 # maximum number of iterations
tol = 0.01 # relative tolerance for convergence, sufficient 
gradtol = 1e-6 # relative gradient tolerance for nlm (1e-6 is the default)
Lambdas = matrix(NA, maxiter+1, 4)
Lambdas[1,] = c(rep(1000, 4)) # initial penalty strengths
mods = list()
# object that contains the indices for each random effect/ spline coef vector. If these are of differnt lengths, use list and change below
REind = matrix(3 + 1:(4*(nb)), nrow = 4, byrow = TRUE)[,-1] # each row is the index of one RE

# initial parameter
theta.star = c(rep(-4, 2), 1,
               # 3, rep(0, nb-1), 6, rep(0, nb-1),
               4.5, seq(-2,2, length = nb-1), 6.5, seq(-1,1, length = nb-1),
               rep(1, 2*nb))

# updating algorithm
T1 = Sys.time()
for(k in 1:maxiter){
  cat("\n\n- Iteration", k, "-")
  t1 = Sys.time()
  mods[[k]] = nlm(mllk_mgcv, theta.star, X = energy, Z=Z, lambda = Lambdas[k,], S=S, 
                  iterlim = 1000, print.level = 0, hessian = TRUE, gradtol = gradtol) # fitting the model with current penalty strengths
  Sys.time()-t1
  cat("\nEstimation time:", Sys.time()-t1)
  cat("\nIterations:", mods[[k]]$iterations)
  theta.star = mods[[k]]$estimate # saves theta.star as starting value for next iteration
  J_inv = MASS::ginv(mods[[k]]$hessian) # inverting hessian
  # updating all penalty strengths
  for(i in 1:4){
    # computing the effective degrees of freedom
    edoF = sum(diag( # trace
      diag(rep(1, nrow(S))) - Lambdas[k, i] * J_inv[REind[i,], REind[i,]] %*% S
    ))
    # computing each penalty separately via the quadratic form
    penalty = t(theta.star[REind[i,]]) %*% S %*% theta.star[REind[i,]]
    # updating the penalty strength
    Lambdas[k+1, i] = as.numeric((edoF - m) / penalty) # m is correction if S does not have full rank
  }
  cat("\nSmoothing strengths:", round(Lambdas[k+1, ], 4))
  if(max(abs(Lambdas[k+1,] - Lambdas[k,]) / Lambdas[k,]) < tol){
    cat("\n\n")
    break
  }
}
Sys.time()-T1

Lambdas = as.matrix(na.omit(Lambdas))
par(mfrow = c(1,4))
for(j in 1:4){
  plot(Lambdas[,j], type = "l", lwd = 2, main = paste("lambda", j))
}
lambda_hat = Lambdas[nrow(Lambdas),]
round(lambda_hat, 3)

theta.star = mods[[length(mods)]]$estimate

Gamma = LaMa::tpm(theta.star[1:2])
delta = c(1, exp(theta.star[3]))
delta = delta / sum(delta)
beta = matrix(theta.star[3 + 1:(2*(nb))], ncol = 2)
alpha = matrix(theta.star[3 + 2*(nb) + 1:(2*(nb))], ncol = 2)
Mu = Z %*% beta
Sigma = exp(Z %*% alpha)
allprobs = cbind(dnorm(energy$Price, Mu[,1], Sigma[,1]),
                 dnorm(energy$Price, Mu[,2], Sigma[,2]))

states = LaMa::viterbi(delta, Gamma, allprobs)

xseq = seq(min(energy$Oil), max(energy$Oil), length = 200)
# Z_plot = mgcv::gam(Price ~ s(Oil, k = nb, bs = "ps"), data = data.frame(dummy = 1, Price=1, Oil = xseq), fit = FALSE)$X
Z_plot = cbind(1, splines::splineDesign(knots, xseq, outer.ok = TRUE))

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

plotind = c(1,2,3,5,10,20)

xseq = seq(min(energy$Oil), max(energy$Oil), length = 200)
Z_plot = cbind(1, splines::splineDesign(knots, xseq, outer.ok = TRUE))

pdf("./case_studies/figs/energy_oil_modseq.pdf", width = 8, height = 5.5)

par(mfrow = c(2,3), mar = c(5,4,3,1))

for(m in plotind){
  theta.star = mods[[m]]$estimate
  
  beta = matrix(theta.star[3 + 1:(2*(nb))], ncol = 2)
  alpha = matrix(theta.star[3 + 2*(nb) + 1:(2*(nb))], ncol = 2)
  Mu = Z %*% beta
  Sigma = exp(Z %*% alpha)
  allprobs = cbind(dnorm(energy$Price, Mu[,1], Sigma[,1]),
                   dnorm(energy$Price, Mu[,2], Sigma[,2]))
  states = LaMa::viterbi(delta, Gamma, allprobs)
  
  Mu_plot = Z_plot %*% beta
  Sigma_plot = exp(Z_plot %*% alpha)
  
  plot(energy$Oil, energy$Price, pch = 20, bty = "n", col = scales::alpha(color[states], 0.1),
       xlab = "oil price", ylab = "energy price", main = paste("Iteration", m))
  for(j in 1:2) lines(xseq, Mu_plot[,j], col = color[j], lwd = 2)
}

dev.off()











# Testing with Gas --------------------------------------------------------

# make my own design and penalty matrix
nb = 15
z = energy$Gas
d = diff(seq(min(z), max(z), length = (nb-1)-2))[1]
knots = seq(min(z)-3*d, max(z)+3*d, by = d)
x = seq(min(z), max(z), length = 200)
Z = splines::splineDesign(knots, z, outer.ok = TRUE)
# Z_smooth = splines::splineDesign(knots, x, outer.ok = TRUE)
Z = cbind(1, Z)

# penalty matrix
L = WH:::build_D_mat(nb-1, 2)
# S = t(L[,-1])%*%L[,-1]
S = t(L)%*%L
m = nrow(S) - Matrix::rankMatrix(S)[1]

# IPML
maxiter = 50 # maximum number of iterations
tol = 0.01 # relative tolerance for convergence, sufficient
gradtol = 1e-6 # relative gradient tolerance for nlm (1e-6 is the default)
Lambdas = matrix(NA, maxiter+1, 4)
Lambdas[1,] = c(rep(100, 2), rep(100, 2))
mods = list()
REind = matrix(3 + 1:(4*(nb)), nrow = 4, byrow = TRUE)[,-1] # each row is the index of one RE

# initial parameter
theta.star = c(rep(-4, 2), 1,
               # 3, rep(0, nb-1), 6, rep(0, nb-1),
               4.5, seq(-2,2, length = nb-1), 6.5, seq(-1,1, length = nb-1),
               rep(0, 2*nb))

# updating algorithm
T1 = Sys.time()
for(k in 1:maxiter){
  cat("\n\n- Iteration", k, "-")
  # if(k == 1){
  #   print.level = 2
  # } else 
  print.level = 0
  t1 = Sys.time()
  mods[[k]] = nlm(mllk_mgcv, theta.star, X = energy, Z=Z, lambda = Lambdas[k,], S=S, 
                  iterlim = 1000, print.level = print.level, hessian = TRUE, gradtol = gradtol)
  Sys.time()-t1
  cat("\nEstimation time:", Sys.time()-t1)
  cat("\nIterations:", mods[[k]]$iterations)
  theta.star = mods[[k]]$estimate # saves theta.star as starting value for next iteration
  J_p = mods[[k]]$hessian # assigning hessian
  J_inv = MASS::ginv(J_p)
  # updating all penalty strengths
  for(i in 1:4){
    edoF = sum(diag( # trace
      diag(rep(1, nrow(S))) - Lambdas[k, i] * J_inv[REind[i,], REind[i,]] %*% S
    ))
    penalty = t(theta.star[REind[i,]]) %*% S %*% theta.star[REind[i,]]
    Lambdas[k+1, i] = as.numeric((edoF - 1) / penalty)
  }
  cat("\nSmoothing strengths:", round(Lambdas[k+1, ], 4))
  if(mean(abs(Lambdas[k+1,] - Lambdas[k,]) / Lambdas[k,]) < tol){
    cat("\n\n")
    break
  }
}
Sys.time()-T1

Lambdas = as.matrix(na.omit(Lambdas))
par(mfrow = c(1,4))
for(j in 1:4){
  plot(Lambdas[,j], type = "l", lwd = 2, main = paste("lambda", j))
}

theta.star = mods[[length(mods)]]$estimate

Gamma = LaMa::tpm(theta.star[1:2])
delta = c(1, exp(theta.star[3]))
delta = delta / sum(delta)
beta = matrix(theta.star[3 + 1:(2*(nb))], ncol = 2)
alpha = matrix(theta.star[3 + 2*(nb) + 1:(2*(nb))], ncol = 2)
Mu = Z %*% beta
Sigma = exp(Z %*% alpha)
allprobs = cbind(dnorm(energy$Price, Mu[,1], Sigma[,1]),
                 dnorm(energy$Price, Mu[,2], Sigma[,2]))

states = LaMa::viterbi(delta, Gamma, allprobs)

zseq = seq(min(z), max(z), length = 200)
Z_plot = cbind(1, splines::splineDesign(knots, zseq, outer.ok = TRUE))

Mu_plot = Z_plot %*% beta
Sigma_plot = exp(Z_plot %*% alpha)

color = c("orange", "deepskyblue")

pdf("./case_studies/figs/energy_gas.pdf", width = 8, height = 4)

par(mfrow = c(1,2), mar = c(5,4,3,1))
plot(z, energy$Price, pch = 20, bty = "n", col = scales::alpha(color[states], 0.2),
     xlab = "gas price", ylab = "energy price", xlim = c(0, 60))
# for(j in 1:2) lines(zseq, Mu_plot[,j], lwd = 3)
for(j in 1:2) lines(zseq, Mu_plot[,j], col = color[j], lwd = 2.5)

# qseq = qnorm(seq(0.5, 0.9, length = 2))
qseq = 1
for(i in qseq){
  for(j in 1:2){
    lines(zseq, Mu_plot[,j] + i*Sigma_plot[,j], col = scales::alpha(color[j], 0.7), lwd = 1, lty = 2)
    lines(zseq, Mu_plot[,j] - i*Sigma_plot[,j], col = scales::alpha(color[j], 0.7), lwd = 1, lty = 2)
  }
}

plot(NA, xlim = c(0, nrow(energy)), ylim = c(1,10), bty = "n",
     xlab = "time", ylab = "energy price")
segments(x0 = 1:(nrow(energy)-1), x1 = 2:nrow(energy),
         y0 = energy$Price[-nrow(energy)], y1 = energy$Price[-1], col = color[states[-1]], lwd = 0.5)

dev.off()

par(mfrow = c(1,1))
plot(Sigma_plot[,1], ylim = c(0,1.5), type = "l", col = color[1], lwd = 2, bty = "n",
     xlab = "z", ylab = "standard deviation")
lines(Sigma_plot[,2], col = color[2], lwd = 2)








# New application ---------------------------------------------------------

# France energy prices

data = read.csv("./data/energy_france.csv")
data = data[-(1:2),]

plot(as.Date(data$date), data$price, type = "l")
plot(data$eurdol, data$price)

data$price[which(data$price<=0)] = NA

plot(data$price)
data = data[1:1500,]


mllk_mgcv2 = function(theta.star, x, Z, lambda, S){
  Gamma = LaMa::tpm(theta.star[1:2])
  delta = c(1, exp(theta.star[3]))
  delta = delta / sum(delta)
  
  nb = ncol(Z)
  beta = matrix(theta.star[3 + 1:(2*(nb))], ncol = 2)
  alpha = matrix(theta.star[3 + 2*(nb) + 1:(2*(nb))], ncol = 2)
  
  Mu = Z %*% beta
  Sigma = exp(Z %*% alpha)
  
  allprobs = matrix(1, nrow = length(x), ncol = 2)
  ind = which(!is.na(x))
  allprobs[ind,] = cbind(dnorm(x[ind], Mu[ind,1], Sigma[ind,1]),
                   dnorm(x[ind], Mu[ind,2], Sigma[ind,2]))
  
  pen = t(as.numeric(beta[-1,])) %*% kronecker(diag(lambda[1:2]), S) %*% as.numeric(beta[-1,]) +
    t(as.numeric(alpha[-1,])) %*% kronecker(diag(lambda[3:4]), S) %*% as.numeric(alpha[-1,])
  
  -LaMa::forward(delta, Gamma, allprobs) + 0.5 * pen
}

# make my own design and penalty matrix
nb = 15
# z = data$eurdol
z = data$gas
d = diff(seq(min(z), max(z), length = (nb-1)-2))[1]
knots = seq(min(z)-3*d, max(z)+3*d, by = d)
Z = cbind(1, splines::splineDesign(knots, z, outer.ok = TRUE))

# penalty matrix
L = WH:::build_D_mat(nb-1, 2)
S = t(L)%*%L
m = nrow(S) - Matrix::rankMatrix(S)[1]

# IPML
maxiter = 50 # maximum number of iterations
tol = 0.05 # relative tolerance for convergence, sufficient
gradtol = 1e-6 # relative gradient tolerance for nlm (1e-6 is the default)
Lambdas = matrix(NA, maxiter+1, 4)
Lambdas[1,] = c(rep(1, 2), 50, 150)
mods = list()
REind = matrix(3 + 1:(4*(nb)), nrow = 4, byrow = TRUE)[,-1] # each row is the index of one RE

# initial parameter
theta.star = c(rep(-4, 2), 1, # Gamma and delta
               30, seq(-5,5,length=nb-1), 60, seq(-15,15,length=nb-1), # coefficients for Mu
               log(10), rep(0, nb-1), log(10), rep(0, nb-1)) # coefficients for Sigma

# updating algorithm
T1 = Sys.time()
for(k in 1:maxiter){
  cat("\n\n- Iteration", k, "-")
  print.level = 0
  t1 = Sys.time()
  mods[[k]] = nlm(mllk_mgcv2, theta.star, x = data$price, Z=Z, lambda = Lambdas[k,], S=S, 
                  iterlim = 1000, print.level = 0, hessian = TRUE, gradtol = gradtol)
  Sys.time()-t1
  cat("\nEstimation time:", Sys.time()-t1)
  cat("\nIterations:", mods[[k]]$iterations)
  theta.star = mods[[k]]$estimate # saves theta.star as starting value for next iteration
  J_p = mods[[k]]$hessian # assigning hessian
  J_inv = MASS::ginv(J_p)
  # updating all penalty strengths
  for(i in 1:4){
    edoF = sum(diag( # trace
      diag(rep(1, nrow(S))) - Lambdas[k, i] * J_inv[REind[i,], REind[i,]] %*% S
    ))
    penalty = t(theta.star[REind[i,]]) %*% S %*% theta.star[REind[i,]]
    Lambdas[k+1, i] = as.numeric((edoF - m) / penalty)
  }
  cat("\nSmoothing strengths:", round(Lambdas[k+1, ], 4))
  if(mean(abs(Lambdas[k+1,] - Lambdas[k,]) / Lambdas[k,]) < tol){
    cat("\n\n")
    break
  }
}
Sys.time()-T1

Lambdas = as.matrix(na.omit(Lambdas))
par(mfrow = c(1,4))
for(j in 1:4){
  plot(Lambdas[,j], type = "l", lwd = 2, main = paste("lambda", j))
}

theta.star = mods[[length(mods)]]$estimate

Gamma = LaMa::tpm(theta.star[1:2])
delta = c(1, exp(theta.star[3]))
delta = delta / sum(delta)
beta = matrix(theta.star[3 + 1:(2*(nb))], ncol = 2)
alpha = matrix(theta.star[3 + 2*(nb) + 1:(2*(nb))], ncol = 2)
Mu = Z %*% beta
Sigma = exp(Z %*% alpha)

allprobs = matrix(1, nrow = nrow(data), ncol = 2)
ind = which(!is.na(data$price))
allprobs[ind,] = cbind(dnorm(data$price[ind], Mu[ind,1], Sigma[ind,1]),
                 dnorm(data$price[ind], Mu[ind,2], Sigma[ind,2]))

states = LaMa::viterbi(delta, Gamma, allprobs)

zseq = seq(min(z), max(z), length = 200)
Z_plot = cbind(1, splines::splineDesign(knots, zseq, outer.ok = TRUE))

Mu_plot = Z_plot %*% beta
Sigma_plot = exp(Z_plot %*% alpha)

color = c("orange", "deepskyblue")

# pdf("./case_studies/figs/energy_france.pdf", width = 8, height = 4)

par(mfrow = c(1,2), mar = c(5,4,3,1))
plot(z, data$price, pch = 20, bty = "n", col = scales::alpha(color[states], 0.2),
     xlab = "euro-dollar", ylab = "french energy price")
# for(j in 1:2) lines(zseq, Mu_plot[,j], lwd = 3)
for(j in 1:2) lines(zseq, Mu_plot[,j], col = color[j], lwd = 2.5)

qseq = qnorm(seq(0.5, 0.975, length = 3))
for(i in qseq){
  for(j in 1:2){
    lines(zseq, Mu_plot[,j] + i*Sigma_plot[,j], col = scales::alpha(color[j], 0.5), lwd = 1, lty = 2)
    lines(zseq, Mu_plot[,j] - i*Sigma_plot[,j], col = scales::alpha(color[j], 0.5), lwd = 1, lty = 2)
  }
}

plot(NA, xlim = c(0, nrow(data)), ylim = c(0,180), bty = "n",
     xlab = "time", ylab = "french energy price")
segments(x0 = 1:(nrow(data)-1), x1 = 2:nrow(data),
         y0 = data$price[-nrow(data)], y1 = data$price[-1], col = color[states[-1]], lwd = 0.5)

# dev.off()

# Create a data frame for ggplot
library(ggplot2)
library(ggridges)

plot_data <- data.frame(z = rep(zseq, 2), 
                        Mu = c(Mu_plot[,1], Mu_plot[,2]), 
                        Sigma = c(Sigma_plot[,1], Sigma_plot[,2]), 
                        state = factor(rep(1:2, each = length(zseq))))


p <- ggplot() +
  geom_point(data = data.frame(z = z, price = data$price, state = factor(states)), 
             aes(x = z, y = price, color = state), alpha = 0.3, size = 2) +
  scale_color_manual(values = color) +
  geom_line(data = plot_data, aes(x = z, y = Mu, color = state), size = 1) +
  geom_ribbon(data = plot_data, aes(x = z, ymin = Mu - qnorm(0.7) * Sigma, ymax = Mu + qnorm(0.7) * Sigma, fill = state), alpha = 0.2) +
  geom_ribbon(data = plot_data, aes(x = z, ymin = Mu - qnorm(0.9) * Sigma, ymax = Mu + qnorm(0.9) * Sigma, fill = state), alpha = 0.2) +
  labs(x = "euro-dollar", y = "french energy price") +
  theme_minimal()

p


par(mfrow = c(1,1))
plot(Sigma_plot[,1], type = "l", col = color[1], lwd = 2, bty = "n",
     xlab = "z", ylab = "standard deviation", ylim = c(0,30))
lines(Sigma_plot[,2], col = color[2], lwd = 2)




# Oil as covariate --------------------------------------------------------

# make my own design and penalty matrix
nb = 20
z = data$oil
d = diff(seq(min(z), max(z), length = (nb-1)-2))[1]
knots = seq(min(z)-3*d, max(z)+3*d, by = d)
x = seq(min(z), max(z), length = 200)
Z = splines::splineDesign(knots, z, outer.ok = TRUE)
# Z_smooth = splines::splineDesign(knots, x, outer.ok = TRUE)
Z = cbind(1, Z)

# penalty matrix
L = WH:::build_D_mat(nb-1, 2)
# S = t(L[,-1])%*%L[,-1]
S = t(L)%*%L

# IPML
maxiter = 50 # maximum number of iterations
tol = 0.1 # relative tolerance for convergence, sufficient
gradtol = 1e-5 # relative gradient tolerance for nlm (1e-6 is the default)
Lambdas = matrix(NA, maxiter+1, 4)
Lambdas[1,] = c(rep(0.05, 2), 1, )
mods = list()
REind = matrix(3 + 1:(4*(nb)), nrow = 4, byrow = TRUE)[,-1] # each row is the index of one RE

# initial parameter
theta.star = c(rep(-4, 2), 1, # Gamma and delta
               40, rep(0, nb-1), 70, rep(0, nb-1), # coefficients for Mu
               log(10), rep(0, nb-1), log(15), rep(0, nb-1)) # coefficients for Sigma

# updating algorithm
T1 = Sys.time()
for(k in 1:maxiter){
  cat("\n\n- Iteration", k, "-")
  print.level = 0
  t1 = Sys.time()
  mods[[k]] = nlm(mllk_mgcv2, theta.star, x = data$price, Z=Z, lambda = Lambdas[k,], S=S, 
                  iterlim = 1000, print.level = 0, hessian = TRUE, gradtol = gradtol)
  Sys.time()-t1
  cat("\nEstimation time:", Sys.time()-t1)
  cat("\nIterations:", mods[[k]]$iterations)
  theta.star = mods[[k]]$estimate # saves theta.star as starting value for next iteration
  J_p = mods[[k]]$hessian # assigning hessian
  J_inv = MASS::ginv(J_p)
  # updating all penalty strengths
  for(i in 1:4){
    edoF = sum(diag( # trace
      diag(rep(1, nrow(S))) - Lambdas[k, i] * J_inv[REind[i,], REind[i,]] %*% S
    ))
    penalty = t(theta.star[REind[i,]]) %*% S %*% theta.star[REind[i,]]
    Lambdas[k+1, i] = as.numeric((edoF - 1) / penalty)
  }
  cat("\nSmoothing strengths:", round(Lambdas[k+1, ], 4))
  if(mean(abs(Lambdas[k+1,] - Lambdas[k,]) / Lambdas[k,]) < tol){
    cat("\n\n")
    break
  }
}
Sys.time()-T1

Lambdas = as.matrix(na.omit(Lambdas))
par(mfrow = c(1,4))
for(j in 1:4){
  plot(Lambdas[,j], type = "l", lwd = 2, main = paste("lambda", j))
}

theta.star = mods[[length(mods)]]$estimate

Gamma = LaMa::tpm(theta.star[1:2])
delta = c(1, exp(theta.star[3]))
delta = delta / sum(delta)
beta = matrix(theta.star[3 + 1:(2*(nb))], ncol = 2)
alpha = matrix(theta.star[3 + 2*(nb) + 1:(2*(nb))], ncol = 2)
Mu = Z %*% beta
Sigma = exp(Z %*% alpha)

allprobs = matrix(1, nrow = nrow(data), ncol = 2)
ind = which(!is.na(data$price))
allprobs[ind,] = cbind(dnorm(data$price[ind], Mu[ind,1], Sigma[ind,1]),
                       dnorm(data$price[ind], Mu[ind,2], Sigma[ind,2]))

states = LaMa::viterbi(delta, Gamma, allprobs)

zseq = seq(min(z), max(z), length = 200)
Z_plot = cbind(1, splines::splineDesign(knots, zseq, outer.ok = TRUE))

Mu_plot = Z_plot %*% beta
Sigma_plot = exp(Z_plot %*% alpha)

color = c("orange", "deepskyblue")

# pdf("./case_studies/figs/energy_france.pdf", width = 8, height = 4)

par(mfrow = c(1,2), mar = c(5,4,3,1))
plot(z, data$price, pch = 20, bty = "n", col = scales::alpha(color[states], 0.1),
     xlab = "euro-dollar", ylab = "french energy price")
# for(j in 1:2) lines(zseq, Mu_plot[,j], lwd = 3)
for(j in 1:2) lines(zseq, Mu_plot[,j], col = color[j], lwd = 2.5)

qseq = qnorm(seq(0.5, 0.975, length = 3))
for(i in qseq){
  for(j in 1:2){
    lines(zseq, Mu_plot[,j] + i*Sigma_plot[,j], col = scales::alpha(color[j], 0.5), lwd = 1, lty = 2)
    lines(zseq, Mu_plot[,j] - i*Sigma_plot[,j], col = scales::alpha(color[j], 0.5), lwd = 1, lty = 2)
  }
}

plot(NA, xlim = c(0, nrow(data)), ylim = c(0,180), bty = "n",
     xlab = "time", ylab = "french energy price")
segments(x0 = 1:(nrow(data)-1), x1 = 2:nrow(data),
         y0 = data$price[-nrow(data)], y1 = data$price[-1], col = color[states[-1]], lwd = 0.5)

# dev.off()

