
# Load package and data
library(hmmTMB)

data <- read.csv("./data/petrels.csv")

data = data[,-(5:6)]

data$ID = as.factor(data$ID)

# Transform ID to factor for random effects data$ID <- factor(data$ID)
data$time <- as.POSIXct(data$time)

head(data)

# Initial parameters
step_mean0 <- c(1, 6, 20)
step_sd0 <- c(1, 5, 10)
angle_mean0 <- c(0, 0, 0)
angle_rho0 <- c(0.8, 0.8, 0.9)
par0 <- list(step = list(mean = step_mean0, sd = step_sd0),
             angle = list(mu = angle_mean0, rho = angle_rho0)) 

# Observation distributions
dists <- list(step = "gamma2", angle = "wrpcauchy")

# Create Observation object
obs <- Observation$new(data = data,
                       dists = dists,
                       n_states = 3,
                       par = par0)

# Model formulas
f <- "~ s(ID, bs = 're') + s(d2c, k = 10, bs = 'cs')"
tpm_structure <- matrix(c(".",   f, "~1",
                          f, ".", f,
                          "~1",  f,  "."),
                        ncol = 3, byrow = TRUE)

# Initial transition probabilities
tpm0 <- matrix(c(0.9, 0.1, 0,
                 0.1, 0.8, 0.1,
                 0, 0.1, 0.9),
               ncol = 3, byrow = TRUE)

# Create MarkovChain object
hid <- MarkovChain$new(n_states = 3,
                       formula = tpm_structure,
                       data = data,
                       tpm = tpm0,
                       initial_state = "stationary")

# List of fixed parameters
fixpar <- list(obs = c("angle.mu.state1.(Intercept)" = NA,
                       "angle.mu.state2.(Intercept)" = NA,
                       "angle.mu.state3.(Intercept)" = NA),
               hid = c("S1>S3.(Intercept)" = NA,
                       "S3>S1.(Intercept)" = NA))

# Create HMM object
hmm <- HMM$new(obs = obs, hid = hid, fixpar = fixpar)

# Fitting the model
t1 = Sys.time()
hmm$fit()
Sys.time()-t1

# Summary of the model
hmm$plot_ts("lon", "lat") +
  coord_map("mercator") +
  geom_point(size = 0.3) +
  labs(x = "longitude", y = "latitude")

hmm$plot_dist("step") +
  coord_cartesian(ylim = c(0, 0.25)) +
  theme(legend.position = c(0.8, 0.8)) +
  labs(x = "step (km)")

hmm$plot_dist("angle") +
  theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(-pi, pi, by = pi/2),
                     labels = expression(-pi, -pi/2, 0, pi/2, pi))

# Transition prob Pr(3 -> 2)
hmm$plot(what = "tpm", var = "d2c", i = 3, j = 2) +
  labs(x = "distance to centre (km)")

# Stationary state probabilities
hmm$plot(what = "delta", var = "d2c") +
  theme(legend.position = "top", legend.margin = margin(c(0, 0, -10, 0))) +
  labs(title = NULL, x = "distance to centre (km)")

hmm$plot(what = "delta", var = "ID", covs = list(d2c = 1500))



# Using my own code and MMLE ----------------------------------------------

library(LaMa)

mllk = function(theta.star, X, Z, S, lambda, trackInd) {
  N=3
  mu = exp(theta.star[1:N])
  sigma = exp(theta.star[N + 1:N])
  kappa = exp(theta.star[2*N + 1:N])
  p = 3*N
  
  # Transition probabilities
  n_re = nrow(S[[1]])
  nb = nrow(S[[2]])
  
  # coefficients
  beta0 = matrix(theta.star[p + 1:4])
  beta_re = matrix(theta.star[p + 4 + 1:(n_re*4)], ncol = 4)
  beta_sm = matrix(theta.star[p + 4 + n_re*4 + 1:(nb*4)], ncol = 4)
  
  Beta = matrix(0, nrow = ncol(Z), ncol = N*(N-1))
  Beta[1, -c(2,5)] = beta0
  Beta[1, c(2,5)] = -(10^2)
  Beta[2:(n_re+1), -c(2,5)] = beta_re
  Beta[(n_re+2):nrow(Beta), -c(2,5)] = beta_sm
  
  Gamma = LaMa::tpm_g(Z[,-1], t(Beta))
  
  # Initial distribution
  Delta = matrix(NA, n_re, N)
  for(i in 1:n_re){
    Delta[i,] = LaMa::stationary(Gamma[,,trackInd[i]])
  }
  
  # Observation probabilities
  allprobs = matrix(1, nrow = nrow(X), ncol = N)
  ind = which(!is.na(X$step) & !is.na(X$angle))
  for(j in 1:N){
    allprobs[ind, j] = dgamma(X$step[ind],shape=mu[j]^2/sigma[j]^2,scale=sigma[j]^2/mu[j]) * 
      CircStats::dvm(X$angle[ind], 0, kappa[j])
  }
  
  pen = t(as.numeric(beta_re)) %*% kronecker(diag(lambda[1:4]), S[[1]]) %*% as.numeric(beta_re) + 
    t(as.numeric(beta_sm)) %*% kronecker(diag(lambda[5:8]), S[[2]]) %*% as.numeric(beta_sm)
  
  - LaMa::forward_g(Delta, Gamma, allprobs, trackInd = trackInd) + 0.5 * pen
}

# Design and penalty matrix
gam_pre = mgcv::gam(y ~ s(ID, bs = 're') + s(d2c, k = 10, bs = 'cs'), 
                    data = data.frame(dummy = 1, ID = data$ID, d2c = data$d2c, y = 1), fit = FALSE)
Z = gam_pre$X
S = gam_pre$S

m = nrow(S[[2]]) - Matrix::rankMatrix(S[[2]])
# S has full rank


# trackInd for independent tracks
trackInd = calc_trackInd(as.character(data$ID))


### Optimisation

# hyperparameters for outer maximization
maxiter = 100 # maximum number of iterations
tol = 0.01 # relative tolerance for convergence
gradtol = 1e-6 # relative gradient tolerance for nlm

Lambdas = matrix(NA, maxiter, 8)
Lambdas[1,] = c(rep(200, 4), rep(20000, 4))
mods = list()

# defining the indices of the spline coefficients
REind = c(
  as.list(data.frame(t(matrix(3*N + 4 + 1:(nrow(S[[1]])*4), nrow = 4, byrow = TRUE)))),
  as.list(data.frame(t(matrix(3*N + 4 + nrow(S[[1]])*4 + 1:(nrow(S[[2]])*4), nrow = 4, byrow = TRUE))))
)

# Initial values
theta.star = c(
  log(step_mean0), log(step_sd0), log(angle_rho0), # observation parameters
  rep(-2, 4), # beta0
  rep(0, 4 * nrow(S[[1]])), # beta_re
  rep(0, 4 * nrow(S[[2]])) # beta_sm
)


# algorithm
T1 = Sys.time()
for(k in 1:maxiter){
  cat("\n\n- Iteration", k, "-")
  if(k == 1){
    print.level = 2
  } else print.level = 0
  
  t1 = Sys.time()
  mods[[k]] = nlm(mllk, theta.star, X = data, Z = Z, S = S, lambda = Lambdas[k,], trackInd = trackInd,
                  iterlim = 2000, print.level = print.level, gradtol = gradtol, hessian = TRUE, stepmax = 100)
  Sys.time()-t1; cat("\nEstimation time:", Sys.time()-t1); cat("\nIterations:", mods[[k]]$iterations)
  
  theta.star = mods[[k]]$estimate # saves theta.star as starting value for next iteration
  J_inv = MASS::ginv(mods[[k]]$hessian) # inverse penalized hessian
  
  # updating all penalty strengths state-dependent process
  for(i in 1:8){
    if(i %in% 1:4){
      S_i = S[[1]]
    } else{ S_i = S[[2]] }
    edoF = sum(diag(diag(rep(1, nrow(S_i))) - Lambdas[k, i] * J_inv[REind[[i]], REind[[i]]] %*% S_i))
    penalty = t(theta.star[REind[[i]]]) %*% S_i %*% theta.star[REind[[i]]]
    Lambdas[k+1, i] = as.numeric(edoF / penalty)
  }
  Lambdas[k+1, which(Lambdas[k+1] < 0)] = 0 # ensures that the penalty strengths are non-negative
  
  cat("\nSmoothing strengths:", round(Lambdas[k+1,], 4))
  if(mean(abs(Lambdas[k+1,] - Lambdas[k,]) / Lambdas[k,]) < tol){
    # cat("\n\n- Final model fit")
    # t1 = Sys.time()
    # mods[[k+1]] = nlm(mllk_mgcv, theta.star, X = data, N = 2, Z = Z, lambda = Lambdas[k+1,], S = S,
    #                   iterlim = 1000, print.level = print.level, hessian = TRUE, stepmax = 1000)
    # Sys.time()-t1; cat("\nEstimation time:", Sys.time()-t1); cat("\nIterations:", mods[[k+1]]$iterations)
    # 
    # cat("\n\n")
    break
  }
}
cat("\nTotal estimation time:", Sys.time()-T1, "\n")

# Trace plots of lambdas
Lambdas = as.matrix(na.omit(Lambdas))
par(mfrow = c(2,4))
for(i in 1:ncol(Lambdas)){
  plot(Lambdas[,i], type = "l", lwd = 2, main = paste("lambda", i))
}
lambda_hat = Lambdas[nrow(Lambdas),]             

sqrt(1 / lambda_hat[1:4])          


# Extracting the results
N=3
mu = exp(theta.star[1:N])
sigma = exp(theta.star[N + 1:N])
kappa = exp(theta.star[2*N + 1:N])
p = 3*N

# Transition probabilities
n_re = nrow(S[[1]])
nb = nrow(S[[2]])

# coefficients
beta0 = matrix(theta.star[p + 1:4])
beta_re = matrix(theta.star[p + 4 + 1:(n_re*4)], ncol = 4)
beta_sm = matrix(theta.star[p + 4 + n_re*4 + 1:(nb*4)], ncol = 4)

Beta = matrix(0, nrow = ncol(Z), ncol = N*(N-1))
Beta[1, -c(2,5)] = beta0
Beta[1, c(2,5)] = -(10^2)
Beta[2:(n_re+1), -c(2,5)] = beta_re
Beta[(n_re+2):nrow(Beta), -c(2,5)] = beta_sm

Gamma = LaMa::tpm_g(Z[,-1], t(Beta))   

# Initial distribution
Delta = matrix(NA, n_re, N)
for(i in 1:n_re){
  Delta[i,] = LaMa::stationary(Gamma[,,trackInd[i]])
}
delta = colMeans(Delta)

# Observation probabilities
allprobs = matrix(1, nrow = nrow(data), ncol = N)
ind = which(!is.na(data$step) & !is.na(data$angle))
for(j in 1:N){
  allprobs[ind, j] = dgamma(data$step[ind],shape=mu[j]^2/sigma[j]^2,scale=sigma[j]^2/mu[j]) * 
    CircStats::dvm(data$angle[ind], 0, kappa[j])
}

# decoding states
states = c()
for(i in 1:n_re){
  ind_i = which(data$ID == unique(data$ID)[i])
  states = c(states, LaMa::viterbi(Delta[i,], Gamma[,,ind_i[-1]], allprobs[ind_i,]))
}

# state frequencies
delta_bar = prop.table(table(states))


# Visualize state-dependent distributions
color = c("orange", "deepskyblue", "seagreen2")

par(mfrow = c(1,2))
hist(data$step, prob = TRUE, bor = "white", xlim = c(0, 40), breaks = 250)
curve(delta_bar[1]*dgamma(x, shape = mu[1]^2/sigma[1]^2, scale = sigma[1]^2/mu[1]), add = TRUE, col = color[1], lwd = 2)
curve(delta_bar[2]*dgamma(x, shape = mu[2]^2/sigma[2]^2, scale = sigma[2]^2/mu[2]), add = TRUE, col = color[2], lwd = 2)
curve(delta_bar[3]*dgamma(x, shape = mu[3]^2/sigma[3]^2, scale = sigma[3]^2/mu[3]), add = TRUE, col = color[3], lwd = 2)

hist(data$angle, prob = TRUE, bor = "white", ylim = c(0, 1.5), breaks = 30)
curve(delta_bar[1]*CircStats::dvm(x, 0, kappa[1]), add = TRUE, col = color[1], lwd = 2)
curve(delta_bar[2]*CircStats::dvm(x, 0, kappa[2]), add = TRUE, col = color[2], lwd = 2)
curve(delta_bar[3]*CircStats::dvm(x, 0, kappa[3]), add = TRUE, col = color[3], lwd = 2)


# Visualize d2c
i = 3; j = 2
ord = order(data$d2c)

# xseq = seq(min(data$d2c), max(data$d2c), length.out = 200)

Z_plot = mgcv::gam(y ~ s(ID, bs = 're') + s(d2c, k = 10, bs = 'cs'), 
                    data = data.frame(dummy = 1, ID = data$ID, d2c = data$d2c, y = 1), fit = FALSE)$X

# Z_plot = mgcv::gam(y ~ s(d2c, k = 10, bs = 'cs'), 
#                    data = data.frame(dummy = 1, d2c = xseq, y = 1), fit = FALSE)$X

Z_plot = cbind(matrix(0, nrow(data), 10), Z_plot[ord,12:20])
Gamma_plot = LaMa::tpm_g(Z_plot, t(Beta))  

par(mfrow = c(1,1))
plot(data$d2c[ord], Gamma_plot[i,j,], type = "l", col = "deepskyblue", lwd = 2, ylim = c(0,1), 
     xlab = "distance to centre (km)", ylab = "Pr(3 -> 2)", bty = "n")


# Visualizing random effects (for d2c ~ 1500 km)

Z_plot_re = Z_plot[rep(which.min(abs(data$d2c-1500)), 10),]
Z_plot_re[,1:n_re] = diag(rep(1,n_re))

Gamma_re = LaMa::tpm_g(Z_plot_re, t(Beta))
Delta_re = matrix(NA, n_re, N)
for(i in 1:n_re) Delta_re[i,] = LaMa::stationary(Gamma_re[,,i])

plot(1:10, Delta_re[,1], ylim = c(0,1), col = color[1], pch = 16, bty = "n",
     xlab = "ID", ylab = "Pr(state)")
points(1:10, Delta_re[,2], col = color[2], pch = 16)
points(1:10, Delta_re[,3], col = color[3], pch = 16)
