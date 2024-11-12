# helper function for penalty and qreml
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

penalty2 = function(re_coef, S, lambda, omega = NULL) {

  # Convert re_coef to a list of matrices (even if originally a vector)
  if (!is.list(re_coef)) {
    re_coef = list(re_coef)
  }
  
  re_coef = lapply(re_coef, function(x) {
    if (is.vector(x)) {
      matrix(x, nrow = 1)  # Convert vectors to 1-row matrices
    } else {
      x  # Leave matrices unchanged
    }
  })
  
  # Get number of distinct random effects (of the same structure)
  n_re = length(re_coef)
  
  # Get the number of similar random effects for each distinct random effect
  re_lengths = sapply(re_coef, nrow)  # All elements are matrices now
  
  if ()
  
  
  
  
  # Ensure S is a list of length n_re, replicating it if necessary
  if (length(S_tp) == 1) {
    S_tp = rep(S_tp, n_re)
  }

  RTMB::REPORT(S_tp) # separately report the tensor product penalty matrices

  # Precompute start and end indices for lambda
  end = cumsum(re_lengths)
  start = c(1, end[-length(end)] + 1)
  
  # Initialize penalty variables
  Pen_lambda = Pen_omega = vector("list", n_re)
  pen = 0
  
  # Loop over distinct random effects - each now a matrix
  for (i in 1:n_re) {
    current_re = re_coef[[i]]  # current_re is always a matrix now
    
    thislambda = lambda[start[i]:end[i]]
    thisomega = omega[start[i]:end[i]]
    
    quadform_lambda = quadform_omega = rep(NaN, nrow(current_re))
    
    # Loop over similar random effects
    for(j in 1:nrow(current_re)){
      S1 = kronecker(S_tp[[i]][[1]], diag(nrow(S_tp[[i]][[2]])))
      S2 = kronecker(diag(nrow(S_tp[[i]][[1]])), S_tp[[i]][[2]])
      
      thisS = thisomega[j] * S1 + (1 - thisomega[j]) * S2
      
      k = nrow(thisS)
      
      quadform_lambda[j] = crossprod(current_re[j,], thisS[-k, -k] %*% current_re[j,])
      quadform_omega[j] = thislambda[j] * crossprod(current_re[j,], (S1-S2)[-k,-k] %*% current_re[j,])
    }
    
    Pen_lambda[[i]] = quadform_lambda
    Pen_omega[[i]] = thislambda[j]
    
    pen = pen + sum(lambda[start[i]:end[i]] * quadform_lambda)
  }
  
  RTMB::REPORT(Pen_lambda)
  RTMB::REPORT(Pen_omega)
  
  0.5 * pen
}

qreml2 = function(pnll, # penalized negative log-likelihood function
                  par, # initial parameter list
                  dat, # initial dat object, currently needs to be called dat!
                  random, # names of parameters in par that are random effects/ penalized
                  psnames = c("lambda", "omega"), # name given to the penalty parameter in dat
                  alpha = 0, # exponential smoothing parameter
                  maxiter = 100, # maximum number of iterations
                  tol = 1e-5, # tolerance for convergence
                  control = list(reltol = 1e-10, maxit = 1000), # control list for inner optimization
                  silent = 1, # print level
                  saveall = FALSE) # save all intermediate models?
{
  
  # setting the argument name for par because later updated par is returned
  argname_par = as.character(substitute(par))
  argname_dat = as.character(substitute(dat))
  
  # setting the environment for mllk to be the local environment such that it pull the right lambda
  environment(pnll) = environment() 
  
  # number of random effects, each one can be a matrix where each row is a random effect, but then they have the same penalty structure
  n_re = length(random) 
  
  # list to save all model objects
  allmods = list() 
  
  # initial lambda and potential omega locally
  psnames_lambda = psnames[1]
  psnames_omega = psnames[2]
  lambda = dat[[psnames_lambda]]
  omega = dat[[psnames_omega]]
  
  # experimentally, changing the name of the data object in pnll to dat
  if(argname_dat != "dat"){
    body(pnll) <- parse(text=gsub(argname_dat, "dat", deparse(body(pnll))))
  }
  
  # creating the objective function as wrapper around pnll to pull lambda from local
  f = function(par){
    environment(pnll) = environment()
    
    "[<-" <- ADoverload("[<-") # overloading assignment operators, currently necessary
    "c" <- ADoverload("c")
    "diag<-" <- ADoverload("diag<-")
    
    getLambda = function(x) lambda
    
    dat[[psnames_lambda]] = DataEval(getLambda, rep(advector(1), 0))
    
    pnll(par)
  }
  
  # creating the RTMB objective function
  if(silent %in% 0:1) cat("Creating AD function\n")
  
  obj = MakeADFun(func = f, parameters = par, silent = TRUE) # silent and replacing with own prints
  newpar = obj$par # saving initial parameter value as vector to initialize optimization in loop
  
  # own printing of maximum gradient component if silent = 0
  if(silent == 0){
    newgrad = function(par){
      gr = obj$gr(par)
      cat(" inner mgc:", max(abs(gr)), "\n")
      gr
    }
  } else{
    newgrad = obj$gr
  }
  
  # prepwork
  mod0 = obj$report() # getting all necessary information from penalty report
  S = mod0$S # penalty matrix/ matrices in list format
  # S_dims = sapply(S, nrow)
  
  # finding the indices of the random effects to later index Hessian
  re_inds = list() 
  for(i in 1:n_re){
    re_dim = dim(as.matrix(par[[random[i]]]))
    # if(re_dim[2] == S_dims[i]){
    #   byrow = FALSE
    # } else{
    #   byrow = TRUE
    # }
    re_inds[[i]] = matrix(which(names(obj$par) == random[i]), nrow = re_dim[1], ncol = re_dim[2])# , byrow = byrow)
    if(dim(re_inds[[i]])[2] == 1) re_inds[[i]] = t(re_inds[[i]]) # if only one column, then transpose
  }
  
  # get number of similar random effects for each distinct random effect (of same structure)
  re_lengths = sapply(re_inds, function(x) if (is.vector(x)) 1 else nrow(x))
  
  # initialize list of penalty strength parameters
  Lambdas = list()
  Lambdas[[1]] = reshape_lambda(re_lengths, lambda) # reshaping to match structure of random effects
  
  if(silent < 2){
    cat("Initializing with", paste0(penalty, ":"), round(lambda, 3), "\n")
  }
  
  # computing rank deficiency for each penalty matrix to use in correction term
  m = numeric(length(S)) 
  for(i in 1:length(m)) {
    m[i] = nrow(S[[i]]) - Matrix::rankMatrix(S[[i]])
  } 
  
  ### updating algorithm
  # loop over outer iterations until convergence or maxiter
  for(k in 1:maxiter){
    
    # fitting the model conditional on lambda: current local lambda will be pulled by f
    opt = stats::optim(newpar, obj$fn, newgrad, 
                       method = "BFGS", hessian = TRUE, # return hessian in the end
                       control = control)
    
    
    # setting new optimum par for next iteration
    newpar = opt$par 
    
    # reporting to extract penalties
    mod = obj$report() 
    
    # evaluating current Hessian
    # J = obj$he()
    J = opt$hessian
    
    # computing inverse Hessian
    J_inv = MASS::ginv(J) 
    
    # saving entire model object
    if(saveall){
      allmods[[k]] = mod
    }
    
    ## updating all lambdas
    lambdas_k = list() # temporary lambda list
    
    # looping over distinct random effects (matrices)
    for(i in 1:n_re){
      # initializing lambda vector for i-th random effect
      lambdas_k[[i]] = numeric(nrow(re_inds[[i]]))
      
      # looping over similar random effects
      for(j in 1:nrow(re_inds[[i]])){
        idx = re_inds[[i]][j,] # indices of this random effect
        
        # effective degrees of freedom for this random effect: J^-1_p J
        edoF = nrow(S[[i]]) - Lambdas[[k]][[i]][j] * sum(rowSums(J_inv[idx, idx] * S[[i]])) # trace(J^-1 \lambda S)
        
        # calculating new lambda based on updating rule
        lambda_new = as.numeric((edoF - m[i]) / mod$Pen[[i]][j]) # m is correction if S_i does not have full rank
        
        # potentially smoothing new lambda
        lambdas_k[[i]][j] = (1-alpha) * lambda_new + alpha * Lambdas[[k]][[i]][j]
      }
      
      # minimum of zero for penalty strengths
      lambdas_k[[i]][which(lambdas_k[[i]] < 0)] = 0
    }
    
    # assigning new lambda to global list
    Lambdas[[k+1]] = lambdas_k
    
    # updating lambda vector locally for next iteration
    lambda = unlist(lambdas_k) 
    
    if(silent < 2){
      cat("outer", k, "-", paste0(penalty, ":"), round(lambda, 3), "\n")
    }
    
    # convergence check
    # if(all(abs(lambda - unlist(Lambdas[[k]])) / unlist(Lambdas[[k]])) < tol)){
    if(max(abs(lambda - unlist(Lambdas[[k]])) / unlist(Lambdas[[k]])) < tol){
      if(silent < 2){
        cat("Converged\n")
      }
      break
    }
    
    if(k == maxiter) warning("No convergence\n")
  }
  
  # assign RTMB obj to return object
  mod$obj <- obj
  
  # if all intermediate models should be returned, assign
  if(saveall) {
    mod$allmods = allmods
  }
  
  # assign final lambda to return object
  mod[[penalty]] = lambda
  
  # assigning all lambdas to return object
  mod[[paste0("all_", penalty)]] = Lambdas
  
  # calculating unpenalized log-likelihood at final parameter values
  lambda = rep(0, length(lambda))
  dat[[penalty]] = lambda
  
  # format parameter to list
  skeleton = utils::as.relistable(par)
  parlist = utils::relist(opt$par, skeleton)
  mod[[argname_par]] = parlist # and assing to return object
  
  # assign estimated parameter as vector
  mod[[paste0(argname_par, "_vec")]] = opt$par
  
  # assign log-likelihood at optimum to return object
  mod$llk = -pnll(parlist)
  
  ## calculating effective degrees of freedom for final model
  mod$edf = list()
  for(i in 1:n_re){
    edoF_i = numeric(nrow(re_inds[[i]]))
    for(j in 1:nrow(re_inds[[i]])){
      idx = re_inds[[i]][j,]
      edoF_i[j] = edoF = nrow(S[[i]]) - Lambdas[[k]][[i]][j] * sum(rowSums(J_inv[idx, idx] * S[[i]]))
    }
    mod$edf[[i]] = edoF_i
  }
  
  # number of fixed parameters
  mod$n_fixpar = length(unlist(par[!(names(par) %in% random)]))
  
  # assing conditinoal Hessian
  mod$Hessian_conditional = J
  
  # removing penalty list from model object
  mod = mod[names(mod) != "Pen"] 
  
  #############################
  ### constructing joint object
  parlist$loglambda = log(mod[[penalty]])
  
  # finding the number of similar random effects for each random effect
  # indvec = rep(1:n_re, times = re_lengths)
  
  # computing log determinants
  logdetS = numeric(length(S))
  for(i in 1:length(S)){
    logdetS[i] = determinant(S[[i]])$modulus
  }
  
  ## defining joint negative log-likelihood
  jnll = function(par) {
    
    environment(pnll) = environment()
    
    "[<-" <- ADoverload("[<-") # overloading assignment operators, currently necessary
    "c" <- ADoverload("c")
    "diag<-" <- ADoverload("diag<-")
    
    dat[[penalty]] = exp(par$loglambda)
    
    l_p = -pnll(par[names(par) != "loglambda"])
    
    ## computing additive constants (missing from only penalized likelihood)
    const = 0
    for(i in 1:n_re){
      for(j in 1:nrow(re_inds[[i]])){
        k = length(re_inds[[i]][j,])
        
        if(i == 1){
          loglam = par$loglambda[j]
        } else{
          loglam = par$loglambda[re_lengths[i-1] + j]
        }
        
        const = const - k * log(2*pi) + k * loglam + logdetS[i]
      }
    }
    
    l_joint = l_p + 0.5 * const
    -l_joint
  }
  
  # creating joint AD object
  obj_joint = MakeADFun(jnll, parlist,
                        random = names(par)[names(par) != "loglambda"]) # REML, everything random except lambda
  
  # assigning object to return object
  mod$obj_joint = obj_joint
  
  return(mod)
}

