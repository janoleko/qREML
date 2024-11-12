reshape_lambda <- function(num_elements, lambda) {
  start <- 1
  result <- lapply(num_elements, function(len) {
    # Extract sub-vector from lambda based on the number of elements
    sub_vector <- lambda[start:(start + len - 1)]
    start <<- start + len
    return(sub_vector)
  })
  return(result)
}

penalty_tp = function(tp_coef, S_tp, omega, lambda_tp){
  
  pen_tp = 0
  
  # Convert tp_coef to a list of matrices (even if originally a vector)
  if (!is.list(tp_coef)) {
    tp_coef = list(tp_coef)
  }
  
  tp_coef = lapply(tp_coef, function(x) {
    if (is.vector(x)) {
      matrix(x, nrow = 1)  # Convert vectors to 1-row matrices
    } else {
      x  # Leave matrices unchanged
    }
  })
  
  # Get number of distinct random effects (of the same structure)
  n_tp = length(tp_coef)
  
  # Get the number of similar random effects for each distinct random effect
  tp_lengths = sapply(tp_coef, nrow)  # All elements are matrices now
  
  # Precompute start and end indices for lambda
  end = cumsum(tp_lengths)
  start = c(1, end[-length(end)] + 1)
  
  # S_tp needs to be a list
  ## either a list containing 2 matrices, i.e. only one tp
  ## or a nested list of lists containing 2 matrices, i.e. multiple tps
  if (!is.list(S_tp[[1]])) { # if not nested list, convert to nested list
    S_tp = list(S_tp)
  }
  if (length(S_tp) == 1) { # if nested list but only one tp, replicate it
    S_tp = rep(S_tp, n_tp)
  }
  
  RTMB::REPORT(S_tp) # Report penalty matrix list
  
  # Initialize penalty variables
  Pen_lambda = Pen_omega = vector("list", n_tp)
  
  # Loop over distinct random effects - each now a matrix
  for (i in 1:n_tp) {
    this_tp = tp_coef[[i]]  # current_re is always a matrix now
    
    this_omega = omega[start[i]:end[i]] # get sub-vector of omegas
    this_lambda = lambda_tp[start[i]:end[i]] # get sub-vector of lambdas
    
    quadform_lambda = rep(NaN, nrow(this_tp))
    
    # Loop over similar random effects
    for(j in 1:nrow(this_tp)){
      
      S1 = kronecker(S_tp[[i]][[1]], diag(nrow(S_tp[[i]][[2]])))
      S2 = kronecker(diag(nrow(S_tp[[i]][[1]])), S_tp[[i]][[2]])
      
      thisS = this_omega[j] * S1 + (1 - this_omega[j]) * S2
      
      k = nrow(thisS)
      
      quadform_lambda[j] = crossprod(this_tp[j,], thisS[-k, -k] %*% this_tp[j,])
      # quadform_omega[j] = this_lambda[j] * crossprod(this_tp[j,], (S1-S2)[-k,-k] %*% this_tp[j,])
    }
    
    Pen_lambda[[i]] = quadform_lambda
    # Pen_omega[[i]] = this_lambda * quadform_omega
    
    pen_tp = pen_tp + sum(this_lambda * quadform_lambda)
  }
  
  RTMB::REPORT(Pen_lambda)
  # RTMB::REPORT(Pen_omega)
  
  0.5 * pen_tp
}

library(RTMB)
library(LaMa)
library(mgcv)
library(Matrix)

data = read.csv("/Users/jan-ole/R/Packages_on_Git/qREML/data/elephant_data.csv")
head(data)


## set up design and penalty matrices
k = 10
gam0 = gam(dummy ~ ti(tod, doy, bs = c("cp","cp"), k = k),
                data = data.frame(dummy = 1, tod = data$tod*2-1, doy = data$doy),
                knots = list(tod = seq(0, 24, length = k+1),
                             doy = seq(0, 366, length = k+1)), fit = FALSE)
gam1 = gam(dummy ~ ti(tod, doy, bs = c("cp","cp"), k = k),
           data = data.frame(dummy = 1, tod = data$tod*2-1, doy = data$doy),
           knots = list(tod = seq(0, 24, length = k+1),
                        doy = seq(0, 366, length = k+1)))
Z = gam0$X
S_tp = list(gam1$smooth[[1]]$margin[[1]]$S[[1]], gam1$smooth[[1]]$margin[[2]]$S[[1]])

tp_coef = matrix(rnorm((ncol(Z)-2)*2), nrow = 2)

omega = c(0.5, 0.5)
lambda_tp = c(10,10)

penalty_tp(tp_coef, S_tp, omega, lambda_tp)

# in each iteration, for tp, find uniroot of this guy for each omega
lprime = function(omega, b, lambda, S, values1, values2, J_inv){
  S1 = kronecker(S[[1]], diag(nrow(S[[2]])))
  S2 = kronecker(diag(nrow(S[[1]])), S[[2]])
  
  Sdiff = (S1 - S2)
  k = nrow(Sdiff)
  
  - 0.5 * lambda * crossprod(b, Sdiff[-k,-k] %*% b) - # lambda b^t (S1-S2) b
    0.5 * lambda * sum(rowSums(J_inv * Sdiff[-k,-k])) + # lambda tr(J^-1 (S1-S2))
    0.5 * sum((values1 - values2) / (omega * values1 + (1-omega) * values2)) # tr(S^-1 (S1-S2))
}


xseq = seq(0,1,length = 200)

b1 = seq(1,4,length=9)
b2 = seq(-4,2, length=9)
b = as.numeric(outer(b1,b2))[-81]
lambda = 0.001
S = S
Omega1 = eigen(S[[1]])$values
Omega2 = eigen(S[[2]])$values
J_inv = matrix(10, 80, 80)

yseq = sapply(xseq, function(x) lprime(x, b, lambda, S, Omega1, Omega2, diag(80)))

plot(xseq, yseq, type = "l")

m = mean(c(lprime(1, b, lambda, S, Omega1, Omega2, J_inv), lprime(0, b, lambda, S, Omega1, Omega2, J_inv)))

lprime2 = function(omega) lprime(omega, b, lambda, S, Omega1, Omega2, J_inv) - m

uniroot(lprime2, c(0,1))


yseq = sapply(xseq, lprime2)
plot(xseq, yseq, type = "l")

