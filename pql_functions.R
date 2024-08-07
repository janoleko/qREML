penalty = function(re_coef, S, lambda) {

  # convert paramter matrix to list of length 1 if necessary  
  if(is.matrix(re_coef)) {
    re_coef = list(re_coef)
  }
  
  n_re = length(re_coef) # number of different random effects (matrices)
  
  # convert penalty matrix to list of length n_re if necessary
  if(is.matrix(S)){
    S = list(S)
    S = rep(S, n_re)
  }
  
  # if S is a list of length 1, then repeat it n_re times
  if(is.list(S) & length(S) == 1){
    S = rep(S, n_re)
  }
  
  # if lambda is only vector, convert to list of length 1
  if(!is.list(lambda)){
    lambda = list(lambda)
  }
  
  RTMB::REPORT(re_coef) # parameter list that is reported
  RTMB::REPORT(lambda) # lambda is reported
  RTMB::REPORT(S) # penalty matrix list is reported
  
  Pen = list() # penalty list that is reported and used for update in pql
  pen = 0 # initializing penalty that will be returned
  
  # loop over different random effects (matrices)
  for(i in 1:n_re){
    re = as.matrix(re_coef[[i]])
    
    Pen[[i]] = numeric(nrow(re))
    
    # for each, loop over rows and compute penalty
    for(j in 1:nrow(re)){
      Pen[[i]][j] = t(re[j,]) %*% S[[i]] %*% re[j,]
      
      pen = pen + lambda[[i]][j] * Pen[[i]][j]
    }
  }
  
  RTMB::REPORT(Pen) # reporting the penalty list for pql update
  
  0.5 * pen
}


pql = function(pnll, par, dat, random,
               alpha_sm = 0.98, maxiter = 50, tol = 1e-2, 
               silent = TRUE, saveall = FALSE) {
  
  # loading RTMB and LaMa if necessary
  if(!("RTMB" %in% loadedNamespaces())) library(RTMB)
  if(!("LaMa" %in% loadedNamespaces())) library(LaMa)
  
  # setting the environment for mllk to be the local environment such that it pull the right lambda
  environment(pnll) = environment() 
  
  # number of random effects, each one can be a matrix where each row is a random effect, but then they have the same penalty structure
  n_re = length(random) 
  
  allmods = list() # creating a list to save all model objects
  
  ## finding the indices of the random effect
  obj = MakeADFun(pnll, par)
  mod0 = obj$report() # getting all necessary information from penalty report
  S = mod0$S # penalty matrix/ matrices
  
  # finding the indices of the random effects to later index Hessian
  re_inds = list() 
  for(i in 1:n_re){
    re_dim = dim(as.matrix(par[[random[i]]]))
    re_inds[[i]] = matrix(which(names(obj$par) == random[i]), nrow = re_dim[1], ncol = re_dim[2])
    # for some weird reason the indices are shuffled around by MakeADFun, so not byrow = TRUE
  }
  
  # initializing list of penalty strength parameters with the one contained in dat
  Lambdas = list()
  if(!is.list(dat$lambda)) {
    lambda0 = list(dat$lambda)
  } else{
    lambda0 = dat$lambda
  }
  Lambdas[[1]] = lambda0
  
  cat("\nInitializing with lambda0:", round(unlist(Lambdas[[1]]), 4))
  
  # computing rank deficiency for each penalty matrix to use in correction term
  m = numeric(length(S)) 
  for(i in 1:length(m)) {
    m[i] = nrow(S[[i]]) - Matrix::rankMatrix(S[[i]])
  } 
  
  ## updating algorithm
  Start = Sys.time()
  
  # loop over outer iterations until maxiter or convergence
  for(k in 1:maxiter){
    cat("\nouter", k)
    cat("\n inner fit...")
    
    # creating the objective function
    obj = MakeADFun(func=pnll, parameters=par, silent=silent)
    
    # fitting the model conditional on lambda
    opt = stats::nlminb(start=obj$par, objective=obj$fn, gradient=obj$gr)
    
    mod = obj$report() # reporting to extract random effects
    
    J = obj$he() # saving current Hessian
    J_inv = MASS::ginv(J) # computing Fisher information
    mod$Fisher = J_inv # saving fisher information in model object for convenience
    allmods[[k]] = mod # saving entire model object
    
    ## updating all lambdas
    lambdas_k = list() # temporary lambda list
    
    # looping over random effects (matrices)
    for(i in 1:n_re){
      
      lambdas_k[[i]] = numeric(nrow(re_inds[[i]])) # initializing lambda vector for i-th random effect
      
      # looping over similar random effects
      for(j in 1:nrow(re_inds[[i]])){
        idx = re_inds[[i]][j,] # indices of this random effect
        
        # effective degrees of freedom for this random effect: J^-1_p J
        edoF = nrow(S[[i]]) - sum(diag(Lambdas[[k]][[i]][j] * J_inv[idx, idx] %*% S[[i]]))
        
        # calculating new lambda based on updating rule
        lambda_new = as.numeric((edoF - m[i]) / mod$Pen[[i]][j]) # m is correction if S_i does not have full rank
        
        # potentially smooting new lambda
        lambdas_k[[i]][j] = alpha_sm * lambda_new + (1-alpha_sm) * Lambdas[[k]][[i]][j]
      }
      
      # minimum of zero for penalty strengths
      lambdas_k[[i]][which(lambdas_k[[i]] < 0)] = 0
    
    }
    Lambdas[[k+1]] = lambdas_k
    
    if(!is.list(dat$lambda)) {
      dat$lambda = unlist(Lambdas[[k+1]])
    } else{
      dat$lambda = Lambdas[[k+1]]
    }
    
    # sdreport to get estimate in list form for good initialization of RE in next iteration
    sdr = sdreport(obj) 
    parlist = as.list(sdr, "Estimate")
    
    for(i in 1:n_re) {
      par[[random[i]]] = parlist[[random[i]]]
    }
    
    cat("\n lambda:", round(unlist(Lambdas[[k+1]]), 4))
    
    # convergence check
    if(max(abs(unlist(Lambdas[[k+1]]) - unlist(Lambdas[[k]])) / unlist(Lambdas[[k]])) < tol){
      cat("\nConverged\n")
      break
    }
    
    if(k == maxiter) cat("\nNo convergence\n")
  }
  Sys.time() - Start
  
  mod$obj = obj
  
  if(saveall) {
    mod$allmods = allmods
  }
  
  return(mod)
}


# Regression matrix functions ---------------------------------------------

make_matrices = function(formula, data){
  gam_setup = gam(formula = update(formula, dummy ~ .),
                  data = cbind(dummy = 1, data), 
                  fit = FALSE)
  Z = gam_setup$X
  S = gam_setup$S
  formula = gam_setup$formula
  return(list(Z = Z, S = S, formula = formula, data = data))
}
pred_matrix = function(model_matrices, newdata) {
  gam_setup0 = gam(model_matrices$form, 
                   data = cbind(dummy = 1, model_matrices$data))
  predict(gam_setup0, newdata = cbind(dummy = 1, newdata), type = "lpmatrix")
}


# Density matrix functions ------------------------------------------------

make_matrices_dens = function(x, K, degree = 3, npoints = 1e4, diff_order = 2){
  ## building the design matrix
  rangex = range(x, na.rm = TRUE)
  nObs = length(x)
  ord = degree + 1
  nrknots = K - (degree-1) 
  d = (rangex[2] - rangex[1]) / nrknots
  bm = c(rangex[1] - degree*d, rangex[2] + degree*d)
  knots = seq(bm[1], bm[2], length = nrknots + 2*degree)
  # numerical integration for normalizing the B-spline basis functions
  xseq =  seq(bm[1], bm[2], length = npoints)
  B0 = splines::spline.des(knots, xseq, degree+1, outer.ok=T)$design # unnormalized
  w = rep(NA, K)
  h = diff(c(knots[1], knots[length(knots)])) / npoints
  for (k in 1:K){
    w[k] = (h* sum(B0[,k]))^(-1) 
    # this computes the integrals of the B-spline basis functions (which are then standardized below)
  } 
  # actual data design matrix
  B = matrix(NA, nrow = nObs, ncol = K)
  ind = which(!is.na(x))
  B[ind,] = t(t(splines::spline.des(knots, x[ind], degree+1, outer.ok = TRUE)$design) * w) 
  
  # basis positions
  basis_pos = knots[(degree+1):(length(knots)-degree+1)]
  
  ## building the penalty matrix
  L = diff(diag(K), differences = diff_order) # second-order difference matrix
  S = t(L[,-1])%*%L[,-1] # leaving out first column
  cat("Leaving out first column of S, fix first column of parameter matrix at zero!")
  
  list(Z=B, S = S, knots=knots, w=w, degree = degree, basis_pos = basis_pos)
}

pred_matrix_dens = function(model_matrices, x_pred){
  knots = model_matrices$knots
  degree = model_matrices$degree
  w = model_matrices$w
  
  B = splines::spline.des(knots, x_pred, degree+1, outer.ok=T)$design
  sweep(B, 2, w, FUN = "*")
}