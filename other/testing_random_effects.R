
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
system.time(
  hmm$fit()
)

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

## packages
library(LaMa)
library(RTMB)

## data
data = read.csv("./data/petrels.csv")
data = data[,-(5:6)]
data$ID = as.factor(data$ID)

## penalized likelihood
pnll = function(par){
  getAll(par, dat)
  
  mu = exp(logmu); REPORT(mu)
  sigma = exp(logsigma); REPORT(sigma)
  kappa = exp(logkappa); REPORT(kappa)
  
  beta = matrix(0, 6, ncol(Z))
  beta[c(2,5), 1] = -1e3
  beta[-c(2,5),] = cbind(beta0, betaRe, betaSpline)
  
  Gamma = tpm_g(Z, beta, ad = T, byrow = T)
  
  uIDs = unique(ID)
  Delta = matrix(NaN, length(uIDs), 3)
  for(i in 1:length(uIDs)) Delta[i,] = stationary(Gamma[,,which(ID == uIDs[i])[1]])

  allprobs = matrix(1, length(step), 3)
  ind = which(!is.na(step) & !is.na(angle))
  for(j in 1:3) allprobs[ind,j] = dgamma2(step[ind], mu[j], sigma[j]) * dvm(angle[ind], 0, kappa[j])
  
  -forward_g(Delta, Gamma, allprobs, trackID = ID, ad = T) + 
    penalty(list(betaRe, betaSpline), S, lambda)
}

## model matrices
k = 10
modmat = make_matrices( ~s(ID, bs = 're') + s(d2c, k = k, bs = 'cs'), data = data)
Z = modmat$Z
S = modmat$S

## initial parameter
nAnimals = length(unique(data$ID))
par = list(logmu = log(c(1, 6, 20)),
           logsigma = log(c(1, 5, 10)),
           logkappa = log(c(0.8, 0.8, 0.9)),
           beta0 = rep(-2,4),
           betaRe = matrix(0, 4, nAnimals),
           betaSpline = matrix(0, 4, k-1))

dat = list(step = data$step, angle = data$angle, ID = data$ID,
           Z = Z, S = S, lambda = c(rep(1e2, 4), rep(1e3, 4)), alpha = 0.1)

system.time(
  mod <- pql(pnll, par, dat, random = c("betaRe", "betaSpline"),
             alpha = 0.2, inner_tol = 1e-10, tol = 1e-6)
)

# extracting parameters
mu = mod$mu
sigma = mod$sigma
kappa = mod$kappa
beta = mod$beta
Delta = mod$delta

# RE standard deviations
sds = 1/sqrt(mod$lambda[1:4])
names(sds) = paste0("beta", c("12", "21", "23", "32"))
sds

# decoding states
states = viterbi_g(Delta, mod$Gamma, mod$allprobs, mod$trackID)
# state frequencies
delta = prop.table(table(states))


# Visualize state-dependent distributions
color = c("orange", "deepskyblue", "seagreen2")

par(mfrow = c(1,2))
hist(data$step, prob = TRUE, bor = "white", xlim = c(0, 40), breaks = 250, xlab = "step", main = "")
for(j in 1:3) curve(delta[j] * dgamma2(x, mu[j], sigma[j]), add = TRUE, col = color[j], lwd = 2)

hist(data$angle, prob = TRUE, bor = "white", ylim = c(0, 1.5), breaks = 30, xlab = "angle", main = "")
for(j in 1:3) curve(delta[j] * dvm(x, 0, kappa[j]), add = TRUE, col = color[j], lwd = 2)



# Visualize d2c
i = 3; j = 2
d2cseq = seq(min(data$d2c), max(data$d2c), length.out = 200)

uIDs = unique(data$ID)

Z_pred = pred_matrix(modmat, data.frame(d2c = d2cseq, ID = uIDs[1]))
Gamma_pred = tpm_g(Z_pred, beta, byrow = T)  

par(mfrow = c(1,1))
plot(d2cseq, Gamma_pred[i,j,], type = "l", col = "deepskyblue", lwd = 2, ylim = c(0,1), 
     xlab = "distance to centre (km)", ylab = "Pr(3 -> 2)", bty = "n")


# Visualizing random effects (for d2c ~ 1500 km)

Z_pred2 = pred_matrix(modmat, data.frame(d2c = 1500, ID = uIDs))
GammaID = tpm_g(Z_pred2, beta, byrow = T)
DeltaRe = t(sapply(1:length(uIDs), function(i) stationary(GammaID[,,i])))

plot(NA, xlim = c(0,10), ylim = c(0,1), bty = "n", xlab = "ID", ylab = "Pr(state)")
for(i in 1:3) points(1:10, DeltaRe[,i], col = color[i], pch = 16)

