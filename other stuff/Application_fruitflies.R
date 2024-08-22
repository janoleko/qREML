
## packages
# install.packages("RTMB")
library(RTMB)
# devtools::install_github("janoleko/LaMa")
library(LaMa)
library(Matrix)


## data
load("./data/fruitflies.RData")
data$time = data$tod / 2

## color
color = c("orange", "deepskyblue")


## penalized log-likelihood function
pnll = function(par){
  "[<-" <- ADoverload("[<-")
  getAll(par, dat) # extract everything from lists
  
  # state-dependent process parameters (mean and dispersion)
  mu = exp(logmu)
  phi = exp(logphi)
  
  beta = cbind(beta0, betaRI[1:2,], betaRI[3:4,], betaSpline[1:2,], betaSpline[3:4,])
  # intercept LD and DD, random intercepts LD and DD, spline coefficients LD and DD
  # looks a little complicated because interaction with light schedule inflates design matrix
  
  Gamma = tpm_g(Z, beta, ad = T)
  
  # periodically stationary initial distribution for each fly
  nAnimals = length(trackInd)
  Delta = matrix(0, nAnimals, 2)
  for(a in 1:nAnimals) Delta[a,] = stationary_p(Gamma[,,trackInd[a] + 0:47], t = 1, ad = T)

  # state-dependent probabilities
  allprobs = matrix(1, nrow = length(activity), ncol = 2)
  ind = which(!is.na(activity))
  size = 1/phi; prob = size/(size + mu) # reparametrization of negbinom
  allprobs[ind,] = cbind(dnbinom(activity[ind], prob[1], size = size[1]),
                         dnbinom(activity[ind], prob[2], size = size[2]))
  
  # forward algorithm + penalty()                      
  - forward_g(Delta, Gamma, allprobs, trackID, ad = T) + 
    penalty(list(betaRI, betaSpline), S, lambda)
}


## model prep
nb = 10 # number of basis functions
modmat = make_matrices(~ condition + s(ID, bs = "re", by = condition) + 
                         s(time, bs = "cp", k = nb, by = condition),
                       data = data,
                       knots = list(time = c(0,24)))

Z = modmat$Z # design matrix
S = modmat$S # list of 4 penalty matrices
S = S[c(1,3)] # only 2 unique ones -> one for RE, one for spline

nAnimals = length(unique(data$ID))

# initial parameters
par = list(logmu = log(c(4, 55)),             # state-dependent mean
           logphi = log(c(10, 0.5)),          # state-dependent dispersion
           beta0 = matrix(-2, 2, 2),          # state process intercepts
           betaRI = matrix(0, 4, nAnimals),   # state process random intercepts
           betaSpline = matrix(0, 4, (nb-1))) # state process spline coefficients

# data list
dat = list(activity = data$activity,
           condition = data$condition,
           trackID = data$ID,
           trackInd = calc_trackInd(as.vector(data$ID)),
           Z = Z, 
           S = S,
           lambda = rep(100, 8)) # initial lambda: lenght equals total number of random effects


## model fitting
system.time(
  mod <- pql(pnll,                               # passing penalized log-likelihood function
             par,                                # passing initial parameter list
             dat,                                # passing dat list with lambda!
             random = c("betaRI", "betaSpline")) # specifying random effects
)

## extracting parameters
beta = mod$beta

## plotting periodically stationary distribution
L = 48 # number of unique times of day
Delta_mean = array(dim = c(L, 2, 2)) # mean array
Delta = array(dim = c(L, 2, 2, nAnimals)) # array for all flies
IDs = unique(data$ID) # unique IDs
conds = unique(data$condition) # unique conditions

# estimated stationary distribution for each fly
for(cond in 1:2){
  for(a in 1:nAnimals){
    Z_pred = pred_matrix(modmat, 
                         data.frame(time = 1:L/2, 
                                    ID = IDs[a],
                                    condition = conds[cond]))
    Gamma_here = tpm_g(Z_pred, beta)
    Delta[,,cond,a] = stationary_p(Gamma_here)
  }
}
for(cond in 1:2){
  Z_pred = pred_matrix(modmat, 
                       data.frame(time = 1:L/2, 
                                  ID = IDs[1],
                                  condition = conds[cond]))
  Z_pred[,2 + (1-(cond-1))*nAnimals + 1] = 0 # pretty hacky but works
  Gamma_here = tpm_g(Z_pred, beta)
  Delta_mean[,,cond] = stationary_p(Gamma_here)
}


# plotting
conditions = c("LD", "DD")
cond = 2 # 1 = LD, 2 = DD
state = 2

par(mfrow = c(1,2))
for(cond in 1:2){
  plot(NA, xlim = c(0,24), ylim = c(0,1), xlab = "time of day", ylab = "Pr(active)", 
       bty = "n", main = conditions[cond], xaxt = "n")
  axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))
  # individual flies
  for(a in 1:nAnimals){
    lines(1:L/2, Delta[,state,cond,a], col = scales::alpha(a, 0.3), type = "b", pch = 20)
  }
  # mean
  points(1:L/2, Delta_mean[,state,cond], type = "b", lwd = 2.5, pch = 16)
}


## plotting periodically stationary distribution
n = 200 # plotting smoothness
todseq = seq(0, 24, length.out = n) # time of day sequence
Delta_mean = array(dim = c(n, 2, 2)) # mean array
Delta = array(dim = c(n, 2, 2, nAnimals)) # array for all flies
IDs = unique(data$ID) # unique IDs
conds = unique(data$condition) # unique conditions

# estimated stationary distribution for each fly (can take some time)
for(t in 1:n){
  print(t)
  time = (todseq[t] + 1:48 / 2) %% 24
  for(cond in 1:2){
    for(a in 1:nAnimals){
      Z_pred = pred_matrix(modmat, 
                           data.frame(time = time, 
                                      ID = IDs[a],
                                      condition = conds[cond]))
      Gamma_t = tpm_g(Z_pred, beta)
      Delta[t,,cond,a] = stationary_p(Gamma_t, t = 1)
    }
  }
}
# mean stationary distribution
for(t in 1:n){
  print(t)
  time = (todseq[t] + 1:48 / 2) %% 24
  for(cond in 1:2){
      Z_pred = pred_matrix(modmat, 
                           data.frame(time = time, 
                                      ID = IDs[1],
                                      condition = conds[cond]))
      Z_pred[,2 + (1-(cond-1))*nAnimals + 1] = 0 # pretty hacky but works
      Gamma_t = tpm_g(Z_pred, beta)
      Delta_mean[t,,cond] = stationary_p(Gamma_t, t = 1)
  }
}

# plotting
conditions = c("LD", "DD")
cond = 2 # 1 = LD, 2 = DD
state = 2

par(mfrow = c(1,2))
for(cond in 1:2){
  plot(NA, xlim = c(0,24), ylim = c(0,1), xlab = "time of day", ylab = "Pr(active)", 
       bty = "n", main = conditions[cond], xaxt = "n")
  axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))
  # individual flies
  for(a in 1:nAnimals){
    lines(seq(0,24,length=n), Delta[,state,cond,a], col = scales::alpha(a, 0.5))
  }
  # mean
  lines(seq(0,24,length=n), Delta_mean[,state,cond], col = "black", lwd = 3)
}





# Just for fun: adding state-dependent random mean AND dispersion ---------

# for efficient indexing in likelihood
data$IDnum = as.numeric(data$ID)

## penalized log-likelihood function
pnll2 = function(par){
  #"[<-" <- ADoverload("[<-")
  getAll(par, dat) # extract everything from lists
  
  # state-dependent process parameters (mean and dispersion)
  mu = exp(logmu + logmuRI); REPORT(mu) # with random means
  phi = exp(logphi + logphiRI); REPORT(phi) # and random dispersions
  
  # state process parameters
  beta = cbind(beta0, betaRI[1:2,], betaRI[3:4,], betaSpline[1:2,], betaSpline[3:4,])
  betaLD = cbind(beta0[,1], betaRI[1:2,], betaSpline[1:2,]); REPORT(betaLD) # reporting
  betaDD = cbind(beta0[,2], betaRI[3:4,], betaSpline[3:4,]); REPORT(betaDD) # reporting
  
  # tpm calculation
  Gamma = tpm_g(Z, beta, ad = T)
  
  # periodically stationary initial distribution for each fly
  nAnimals = length(trackInd)
  Delta = matrix(0, nAnimals, 2)
  for(a in 1:nAnimals) Delta[a,] = stationary_p(Gamma[,,trackInd[a] + 0:47], t = 1, ad = T)
  
  # state-dependent probabilities
  allprobs = matrix(1, nrow = length(activity), ncol = 2)
  ind = which(!is.na(activity))
  size = 1/phi; prob = size/(size + mu) # reparametrization of negbinom
  allprobs[ind,] = cbind(dnbinom(activity[ind], size[1, IDnum[ind]], prob[1, IDnum[ind]]),
                         dnbinom(activity[ind], size[2, IDnum[ind]], prob[2, IDnum[ind]]))
  
  # forward algorithm + penalty()                      
  - forward_g(Delta, Gamma, allprobs, trackID, ad = T) + 
    penalty(list(logmuRI, logphiRI, betaRI, betaSpline), S, lambda)
}


## model prep
nb = 10 # number of basis functions
modmat = make_matrices(~ condition + s(ID, bs = "re", by = condition) + 
                         s(time, bs = "cp", k = nb, by = condition),
                       data = data,
                       knots = list(time = c(0,24)))

Z = modmat$Z # design matrix
S = modmat$S # list of 4 penalty matrices
S = S[c(1,1,1,3)] # 4 different REs: 1:3 iid RE, 4 spline

nAnimals = length(unique(data$ID))

# initial parameters
par = list(logmu = log(c(4, 55)),             # state-dependent mean
           logphi = log(c(10, 0.5)),          # state-dependent dispersion
           logmuRI = matrix(0, 2, nAnimals),  # state-dependent random mean
           logphiRI = matrix(0, 2, nAnimals), # state-dependent random dispersion
           beta0 = matrix(-2, 2, 2),          # state process intercepts
           betaRI = matrix(0, 4, nAnimals),   # state process random intercepts
           betaSpline = matrix(0, 4, (nb-1))) # state process spline coefficients

# data list
dat = list(activity = data$activity,
           condition = data$condition,
           trackID = data$ID,
           trackInd = calc_trackInd(as.vector(data$ID)),
           IDnum = data$IDnum,
           Z = Z, 
           S = S,
           lambda = rep(100, 12)) # initial lambda: length equals total number of random effects


## model fitting
system.time(
  mod2 <- pql(pnll2, par, dat,
              random = c("logmuRI", "logphiRI", "betaRI", "betaSpline")) # specifying random effects
)

## number of parameters and AIC BIC
npar = length(mod2$par_vec)
npar_eff = mod2$n_fixpar + sum(unlist(mod2$edf))
llk = mod2$llk
(AIC = -2 * llk + 2 * npar_eff)
(BIC = -2 * llk + log(nrow(data)) * npar_eff)


## extracting parameters
beta = mod2$beta
mu = mod2$mu
phi = mod2$phi
mu_mean = exp(mod2$par$logmu)
phi_mean = exp(mod2$par$logphi)

## RE distribution
lambda = mod$lambda
# sds of random mean
1 / sqrt(lambda[1:2])
# sds of random dispersion
1 / sqrt(lambda[3:4])
# sds of random intercepts
matrix(1 / sqrt(lambda[5:8]), nrow = 2)

## plot all state-dependent distributions
par(mfrow = c(1,1))
plot(NA, xlim = c(0, 130), ylim = c(0,0.05), ylab = "prob", xlab = "activity", bty = "n")
for(state in 1:2){
  for(k in 1:nAnimals){
    lines(stats::dnbinom(0:130, mu = mu[state,k], size = 1/phi[state,k]), 
          col = scales::alpha(color[state], 0.2), lwd = 1)
  }
  points(stats::dnbinom(0:130, mu = mu_mean[state], size = 1/phi_mean[state]), 
        col = scales::alpha(color[state], 0.9), pch = 20)
}


## plotting periodically stationary distribution
L = 48 # number of unique times of day
Delta_mean = array(dim = c(L, 2, 2)) # mean array
Delta = array(dim = c(L, 2, 2, nAnimals)) # array for all flies
IDs = unique(data$ID) # unique IDs
conds = unique(data$condition) # unique conditions

# estimated stationary distribution for each fly (can take some time)
for(cond in 1:2){
  for(a in 1:nAnimals){
    Z_pred = pred_matrix(modmat, 
                         data.frame(time = 1:L/2, 
                                    ID = IDs[a],
                                    condition = conds[cond]))
    Gamma_here = tpm_g(Z_pred, beta)
    Delta[,,cond,a] = stationary_p(Gamma_here)
  }
}
for(cond in 1:2){
  Z_pred = pred_matrix(modmat, 
                       data.frame(time = 1:L/2, 
                                  ID = IDs[1],
                                  condition = conds[cond]))
  Z_pred[,2 + (1-(cond-1))*nAnimals + 1] = 0 # pretty hacky but works
  Gamma_here = tpm_g(Z_pred, beta)
  Delta_mean[,,cond] = stationary_p(Gamma_here)
}


# plotting
conditions = c("LD", "DD")
cond = 2 # 1 = LD, 2 = DD
state = 2

par(mfrow = c(1,2))
for(cond in 1:2){
  plot(NA, xlim = c(0,24), ylim = c(0,1), xlab = "time of day", ylab = "Pr(active)", 
       bty = "n", main = conditions[cond], xaxt = "n")
  axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))
  # individual flies
  for(a in 1:nAnimals){
    lines(1:L/2, Delta[,state,cond,a], col = scales::alpha(a, 0.3), type = "b", pch = 20)
  }
  # mean
  points(1:L/2, Delta_mean[,state,cond], type = "b", lwd = 2.5, pch = 16)
}


## smooth plotting
n = 200 # plotting smoothness
todseq = seq(0, 24, length.out = n) # time of day sequence
Delta_mean = array(dim = c(n, 2, 2)) # mean array
Delta = array(dim = c(n, 2, 2, nAnimals)) # array for all flies
IDs = unique(data$ID) # unique IDs
conds = unique(data$condition) # unique conditions

# estimated stationary distribution for each fly (can take some time)
for(t in 1:n){
  print(t)
  time = (todseq[t] + 1:48 / 2) %% 24
  for(cond in 1:2){
    for(a in 1:nAnimals){
      Z_pred = pred_matrix(modmat, 
                           data.frame(time = time, 
                                      ID = IDs[a],
                                      condition = conds[cond]))
      Gamma_t = tpm_g(Z_pred, beta)
      Delta[t,,cond,a] = stationary_p(Gamma_t, t = 1)
    }
  }
}
for(t in 1:n){
  print(t)
  time = (todseq[t] + 1:48 / 2) %% 24
  for(cond in 1:2){
    Z_pred = pred_matrix(modmat, 
                         data.frame(time = time, 
                                    ID = IDs[1],
                                    condition = conds[cond]))
    Z_pred[,2 + (1-(cond-1))*nAnimals + 1] = 0 # pretty hacky but works
    Gamma_t = tpm_g(Z_pred, beta)
    Delta_mean[t,,cond] = stationary_p(Gamma_t, t = 1)
  }
}

# plotting
conditions = c("LD", "DD")
cond = 2 # 1 = LD, 2 = DD
state = 2

par(mfrow = c(1,2))
for(cond in 1:2){
  plot(NA, xlim = c(0,24), ylim = c(0,1), xlab = "time of day", ylab = "Pr(active)", 
       bty = "n", main = conditions[cond], xaxt = "n")
  axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))
  # individual flies
  for(a in 1:nAnimals){
    lines(seq(0,24,length=n), Delta[,state,cond,a], col = scales::alpha(a, 0.5))
  }
  # mean
  lines(seq(0,24,length=n), Delta_mean[,state,cond], col = "black", lwd = 3)
}

par = mod2$par
parvec = unlist(mod2$par)
skeleton = as.relistable(par)
relist(parvec, skeleton)

obj = mod$obj

obj$par = obj$par
mo = obj$report(obj$par+2)
mo$beta[,1]



# Same with hmmTMB --------------------------------------------------------

library(hmmTMB)


# Initial parameters for the state-dependent process
size0 <- c(0.05, 1.5)
prob0 <- c(0.02, 0.02)

# Initial parameter is a named list, with a list for each data stream
par0 <- list(activity = list(size = size0, prob = prob0)) 

# We need to specify the observation distributions (distribution and initial parameter list need to match).
dists <- list(activity = "nbinom")

# Now, we can create Observation object
obs <- Observation$new(data = data, # data frame that contains the observation variable, potentially an ID column and covariates
                       dists = dists, # list of distributions
                       n_states = 2, # number of hidden states
                       par = par0) # initial parameter list

# Then, we create a MarkovChain object
hid <- MarkovChain$new(formula = ~ condition + s(ID, bs = "re", by = condition) + 
                         s(time, bs = "cp", k = 10, by = condition),
                       n_states = 2, # number of states, obviously needs to match obs
                       data = data, # data frame
                       initial_state = "stationary") # stationary initial distribution

# Now, create HMM object
hmm <- HMM$new(obs = obs, hid = hid)

# and fit the model
system.time(
  hmm$fit()
)

