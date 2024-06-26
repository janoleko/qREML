library(tidyverse)
library(moveHMM)

dat = read.csv("~/Downloads/Adlie penguins, Sphenisciforms, foraging activity during chick-rearing period in Adlie Land.csv")
  
dat = dat %>%  
  select(ID=tag.local.identifier, x=location.long, y=location.lat, time=timestamp)


# regularizing by imputing NAs

dat$time = ymd_hms(dat$time)

datlist = split(dat, dat$ID)
for(i in 1:length(datlist)){
  timeseries = tibble(time = seq(datlist[[i]]$time[1], datlist[[i]]$time[nrow(datlist[[i]])], by = "min"))
  datlist[[i]] = timeseries %>% left_join(datlist[[i]], by = "time")
  datlist[[i]]$ID = as.character(i)
}

dat = bind_rows(datlist)

# calculating step lengths and turning angles
dat = prepData(dat)

datlist = split(dat, dat$ID)
sapply(datlist, nrow)

hist(dat$step, breaks = 600, prob = TRUE, xlim = c(0,0.4))

## get hours from time
dat$hour = hour(dat$time)
dat$minute = minute(dat$time)
dat$tod = dat$hour*60 + dat$minute
dat$doy = yday(dat$time)

## get day of year
boxplot(dat$step ~ as.factor(dat$tod), ylim = c(0, 0.2))

plot(dat$step[1:200])

dat$step

View(dat)




# Storcs ------------------------------------------------------------------

library(tidyverse)
library(moveHMM)
library(lubridate)

storcs = read_csv("~/Downloads/30min_49_breeding_birds2024.csv")
nrow(storcs)

colnames(storcs)
unique(storcs$individual.local.identifier)

storcs_sub = storcs %>% filter(individual.local.identifier == "Adebar...DER.AX142..eobs.4360.")

nrow(storcs_sub)

colnames(storcs_sub)

plot(storcs_sub$coords.x1, storcs_sub$coords.x2, type = "l")
unique(storcs_sub$trackId)

storcs_sub$timestamp

dat = storcs_sub %>% dplyr::select(x = long, y = lat, time = timestamp)

dat = prepData(dat)
head(dat)

hist(dat$step, breaks = 800, prob = TRUE, xlim = c(0, 10))

# format time
dat$time = ymd_hms(dat$time)
dat$tod = hour(dat$time) + minute(dat$time)/60
dat$doy = yday(dat$time)

boxplot(dat$step ~ as.factor(dat$doy))


hist(dat$angle, breaks = 100)


# fitting HMMs

mllk = function(theta.star, X, N){
  mu = exp(theta.star[1:N])
  sigma = exp(theta.star[N+1:N])
  mu.turn = theta.star[2*N+1:N]
  kappa = exp(theta.star[3*N+1:N])
  
  Gamma = LaMa::tpm(theta.star[4*N + 1:(N*(N-1))])
  delta = LaMa::stationary(Gamma)
  
  allprobs = matrix(1, nrow(X), N)
  ind = which(!is.na(X$step) & !is.na(X$angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma(X$step[ind], shape = mu[j]^2/sigma[j]^2, scale = sigma[j]^2/mu[j])*
      CircStats::dvm(X$angle[ind], mu.turn[j], kappa[j])
  }
  
  -LaMa::forward(delta, Gamma, allprobs)
}

mu0 = c(0.05, 1, 2)
sigma0 = c(0.1, 1, 2)
mu.turn0 = c(0, 0, 0)
kappa0 = c(0.1, 1, 2)

theta.star0 = c(log(c(mu0, sigma0)), mu.turn0, log(kappa0), rep(-2, 6))

mod0 = nlm(mllk, theta.star0, X = dat, N=3, print.level = 2, iterlim = 1000)

theta.star = mod0$estimate; N=3

mu = exp(theta.star[1:N])
sigma = exp(theta.star[N+1:N])
mu.turn = theta.star[2*N+1:N]
kappa = exp(theta.star[3*N+1:N])
Gamma = LaMa::tpm(theta.star[4*N + 1:(N*(N-1))])
delta = LaMa::stationary(Gamma)

allprobs = matrix(1, nrow(dat), N)
ind = which(!is.na(dat$step) & !is.na(dat$angle))
for(j in 1:N){
  allprobs[ind,j] = dgamma(dat$step[ind], shape = mu[j]^2/sigma[j]^2, scale = sigma[j]^2/mu[j])*
    CircStats::dvm(dat$angle[ind], mu.turn[j], kappa[j])
}

color = c("orange", "deepskyblue", "seagreen2")
hist(dat$step, breaks = 2000, prob = TRUE, xlim = c(0, 8), ylim = c(0, 3), bor = "white")
for(j in 1:N){
  curve(delta[j]*dgamma(x, shape = mu[j]^2/sigma[j]^2, scale = sigma[j]^2/mu[j]), add = TRUE, col = color[j], lwd = 2, n = 500)
}

states = LaMa::viterbi(delta, Gamma, allprobs)

plot(dat$x, dat$y, pch = 20, col = color[states], lwd = 2)

library(maps)

# Plot the map
map("world", fill = FALSE, col = "lightgray", bg = "white", xlim = c(-7, 11), ylim = c(35, 50))

# Add your points on top of the map
points(dat$x, dat$y, pch = 20, col = color[states], lwd = 2)


plotind = 2000:6000
plot(dat$step[plotind], pch = 20, col = color[states[plotind]])
nrow(dat)

nrow(dat)/(48*366) # years


## model with inhomogeneity for migratory phase

migInd = 2000:6000
datsub = dat[migInd,]

mllk = function(theta.star, X, N, Z, lambda, D){
  mu = exp(theta.star[1:N])
  sigma = exp(theta.star[N+1:N])
  mu.turn = theta.star[2*N+1:N]
  kappa = exp(theta.star[3*N+1:N]); p = 4*N
  
  beta0 = theta.star[p + 1:(N*(N-1))]; p = p + N*(N-1)
  nb = ncol(Z)
  b = theta.star[p + 1:(N*(N-1)*(nb-1))];
  beta = matrix(0, N*(N-1), nb+1)
  beta[,1] = beta0
  beta[,-(1:2)] = matrix(b, N*(N-1), nb-1, byrow = TRUE)
  
  Gamma = LaMa::tpm_g(Z, beta)
  delta = c(1, exp(theta.star[p + N*(N-1)*(nb-1) + 1:(N-1)]))
  delta = delta / sum(delta)
  
  allprobs = matrix(1, nrow(X), N)
  ind = which(!is.na(X$step) & !is.na(X$angle))
  for(j in 1:N){
    allprobs[ind,j] = dgamma(X$step[ind], shape = mu[j]^2/sigma[j]^2, scale = sigma[j]^2/mu[j])*
      CircStats::dvm(X$angle[ind], mu.turn[j], kappa[j])
  }
  
  pen = b %*% kronecker(diag(lambda), D) %*% b
  
  -LaMa::forward_g(delta, Gamma[,,-1], allprobs) + 0.5 * pen
}

nk = 15 # number of basis functions
k = 24 * 0:nk / nk # knots
Z = mgcv::cSplineDes(dat$tod[migInd], k) ## cyclic spline design matrix

L = WH:::build_D_mat(nb-1, 2)
D = t(L)%*%L

theta.star = c(log(mu), log(sigma), mu.turn, log(kappa), 
               mod0$estimate[(4*N+1):length(mod0$estimate)],
               rep(0, (ncol(Z)-1)*N*(N-1)),
               rep(0, N-1))

mod1 = nlm(mllk, theta.star, X = datsub, N=3, Z = Z, lambda = rep(10, N*(N-1)), D=D,
            print.level = 2, iterlim = 1000)

theta.star = mod1$estimate; N=3

mu = exp(theta.star[1:N])
sigma = exp(theta.star[N+1:N])
mu.turn = theta.star[2*N+1:N]
kappa = exp(theta.star[3*N+1:N]); p = 4*N

beta0 = theta.star[p + 1:(N*(N-1))]; p = p + N*(N-1)
nb = ncol(Z)
b = theta.star[p + 1:(N*(N-1)*(nb-1))];
beta = matrix(0, N*(N-1), nb+1)
beta[,1] = beta0
beta[,-(1:2)] = matrix(b, N*(N-1), nb-1, byrow = TRUE)

todseq = seq(0,24, length = 200)
Zplot = mgcv::cSplineDes(todseq, k)

Gamma = LaMa::tpm_g(Zplot, beta)

plot(Gamma[1,1,], type = "l", col = "orange", lwd = 2)