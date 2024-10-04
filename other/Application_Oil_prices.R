# data(energy, package = "MSwM")
# 
# plot(energy$EurDol, energy$Price, pch = 20, bty = "n")
# 
# # simple HMM
# mllk_sim = function(theta.star, X){
#   Gamma = LaMa::tpm(theta.star[1:2])
#   mu = theta.star[3:4]
#   sigma = exp(theta.star[5:6])
#   delta = c(1, exp(theta.star[7]))
#   delta = delta/sum(delta)
#   allprobs = cbind(dnorm(X$Price, mu[1], sigma[1]),
#                    dnorm(X$Price, mu[2], sigma[2]))
#   -LaMa::forward(delta, Gamma, allprobs)
# }
# 
# theta.star = c(rep(-3,2), 3, 6, log(1), log(1), 0)
# mod_sim = nlm(mllk_sim, theta.star, X = energy, iterlim = 1000, print.level = 2)
# theta.star = mod_sim$estimate
# Gamma = LaMa::tpm(theta.star[1:2])
# mu = theta.star[3:4]
# sigma = exp(theta.star[5:6])
# delta = c(1, exp(theta.star[7]))
# delta = delta/sum(delta)
# allprobs = cbind(dnorm(energy$Price, mu[1], sigma[1]),
#                  dnorm(energy$Price, mu[2], sigma[2]))
# 
# states_sim = LaMa::viterbi(delta, Gamma, allprobs)
# 
# color = c("orange", "deepskyblue")
# plot(energy$EurDol, energy$Price, pch = 16, bty = "n", 
#      col = scales::alpha(color[states_sim], 0.5), xlab = "EurDol", ylab = "Price")
# 
# 
# 
# # MS-GAM
# 
# mllkGAM = function(theta.star, X, Z, lambda, D){
#   Gamma = LaMa::tpm(theta.star[1:2])
#   beta0 = theta.star[3:4]
#   delta = c(1, exp(theta.star[5]))
#   delta = delta / sum(delta)
#   
#   nb = ncol(Z)
#   b = theta.star[5 + 1:(2*(nb-1))]
#   beta = cbind(beta0, c(0,0),
#                matrix(b, nrow = 2, ncol = nb-1, byrow = TRUE))
#   # alpha = matrix(theta.star[5 + 2*(nb-1) + 1:8], nrow = 2)
#   sigma = exp(theta.star[5 + 2*(nb-1) + 1:2])
#   Mu = cbind(1, Z) %*% t(beta)
#   # Sigma = exp(cbind(1, poly(X$EurDol, 3)) %*% t(alpha))
#   allprobs = cbind(dnorm(X$Price, Mu[,1], sigma[1]),
#                    dnorm(X$Price, Mu[,2], sigma[2]))
#   
#   pen = t(b) %*% kronecker(diag(lambda), D) %*% b
#   -LaMa::forward(delta, Gamma, allprobs) + 0.5 * pen
# }
# 
# # design matrix
# # nb = 10
# # gam_setup = mgcv::gam(Price ~ s(EurDol, k = nb+1, bs = "ps"), data = cbind(dummy = 1, energy), fit = FALSE)
# # Z = gam_setup$X
# # S = gam_setup$S[[1]]
# 
# nb = 10
# d = diff(seq(min(energy$EurDol), max(energy$EurDol), length = nb-2))[1]
# knots = seq(min(energy$EurDol)-3*d, max(energy$EurDol)+3*d, by = d)
# x = seq(min(energy$EurDol), max(energy$EurDol), length = 200)
# # Z = splines::splineDesign(knots, x, outer.ok = TRUE)
# Z = splines::splineDesign(knots, energy$EurDol, outer.ok = TRUE)
# 
# # plot(x, Z[,1], type = "l", ylim = c(0,1))
# # for(i in 2:nb) lines(x, Z[,i], col = i)
#   
# # Z = splines::bs(energy$EurDol, df = nb+1)[,-(nb+1)]
# # Z = splines::bs(energy$EurDol, df = nb)
# 
# # ord = 4
# # degree = ord-1
# # nrknots = nb - (degree-1) 
# # knots = seq(0.4, 1.4, length = nrknots+2*degree)
# # Z = splines::spline.des(knots, energy$EurDol, degree+1, outer.ok=T)$design
# 
# # penalty matrix
# L = WH:::build_D_mat(nb, 2)
# D = t(L[,-1])%*%L[,-1]
# 
# theta.star = c(rep(-3,2), 3, 6, 1,
#                rep(0, 2*(nb-1)),
#                log(1), log(1), 
#                rep(0, 2))
# 
# 
# # IPML
# maxiter = 50 # maximum number of iterations
# tol = 0.02 # relative tolerance for convergence, sufficient
# gradtol = 1e-6 # relative gradient tolerance for nlm (1e-6 is the default)
# Lambdas = matrix(NA, maxiter+1, 2)
# Lambdas[1,] = c(100, 100)
# mods = list()
# REind = matrix(5 + 1:(2*(nb-1)), nrow = 2, byrow = TRUE) # each row is the index of one RE
# 
# # updating algorithm
# for(k in 1:maxiter){
#   cat("\n\n- Iteration", k, "-")
#   if(k == 1){
#     print.level = 2
#   } else print.level = 0
#   
#   t1 = Sys.time()
#   mods[[k]] = nlm(mllkGAM, theta.star, X = energy, Z=Z, lambda = Lambdas[k,], D=D, 
#                   iterlim = 1000, print.level = print.level, hessian = TRUE, gradtol = gradtol)
#   Sys.time()-t1
#   cat("\nEstimation time:", Sys.time()-t1)
#   cat("\nIterations:", mods[[k]]$iterations)
#   
#   theta.star = mods[[k]]$estimate # saves theta.star as starting value for next iteration
#   
#   J_p = mods[[k]]$hessian # assigning hessian
#   J_inv = MASS::ginv(J_p)
#   
#   # updating all penalty strengths
#   for(i in 1:2){
#     edoF = sum(diag( # trace
#       diag(rep(1, nrow(D))) - Lambdas[k, i] * J_inv[REind[i,], REind[i,]] %*% D
#     ))
#     penalty = t(theta.star[REind[i,]]) %*% D %*% theta.star[REind[i,]]
#     
#     Lambdas[k+1, i] = as.numeric((edoF - 1) / penalty)
#   }
#   
#   cat("\nSmoothing strengths:", round(Lambdas[k+1, ], 4))
#   
#   if(mean(abs(Lambdas[k+1,] - Lambdas[k,]) / Lambdas[k,]) < tol){
#     cat("\n\n")
#     break
#   }
# }
# 
# Lambdas = as.matrix(na.omit(Lambdas))
# par(mfrow = c(1,2))
# for(j in 1:2){
#   plot(Lambdas[,j], type = "l", lwd = 2, main = paste("lambda", j), ylim = c(0,20))
# }
# 
# 
# theta.star = mods[[length(mods)]]$estimate
# 
# Gamma = LaMa::tpm(theta.star[1:2])
# beta0 = theta.star[3:4]
# delta = c(1, exp(theta.star[5]))
# delta = delta / sum(delta)
# nb = ncol(Z)
# b = theta.star[5 + 1:(2*(nb-1))]
# beta = cbind(beta0, c(0,0),
#              matrix(b, nrow = 2, ncol = nb-1, byrow = TRUE))
# # alpha = matrix(theta.star[5 + 2*(nb-1) + 1:8], nrow = 2)
# sigma = exp(theta.star[5 + 2*(nb-1) + 1:2])
# Mu = cbind(1, Z) %*% t(beta)
# # Sigma = exp(cbind(1, poly(X$EurDol, 3)) %*% t(alpha))
# allprobs = cbind(dnorm(energy$Price, Mu[,1], sigma[1]),
#                  dnorm(energy$Price, Mu[,2], sigma[2]))
# 
# states = LaMa::viterbi(delta, Gamma, allprobs)
# 
# # xseq = seq(min(energy$EurDol), max(energy$EurDol), length = nrow(energy))
# # Z_plot = gam(Price ~ s(EurDol, k = nb+1, bs = "cr"), 
# #              data = data.frame(dummy = 1, EurDol = xseq, Price = energy$Price), fit = FALSE)$X
# 
# 
# Z_plot = splines::splineDesign(knots, xseq, outer.ok = TRUE)
# # Z_plot = splines::bs(xseq, df = nb+1)[,-(nb+1)]
# # Z_plot = splines::bs(xseq, df = nb)
# # Z_plot = splines::spline.des(knots, xseq, degree+1, outer.ok=T)$design
# 
# Mu_plot = cbind(1,Z_plot) %*% t(beta)
# # Sigma_plot = exp(cbind(1, poly(xseq, 3)) %*% t(alpha))
# 
# color = c("orange", "deepskyblue")
# par(mfrow = c(1,1))
# plot(energy$EurDol, energy$Price, pch = 16, bty = "n", col = scales::alpha(color[states], 0.5))
# lines(xseq, Mu_plot[,1], col = color[1], lwd = 3)
# lines(xseq, Mu_plot[,1]+sigma[1], col = color[1], lty = 2, lwd = 2)
# lines(xseq, Mu_plot[,1]-sigma[1], col = color[1], lty = 2, lwd = 2)
# 
# lines(xseq, Mu_plot[,2], col = color[2], lwd = 3)
# lines(xseq, Mu_plot[,2]+sigma[2], col = color[2], lty = 2, lwd = 2)
# lines(xseq, Mu_plot[,2]-sigma[2], col = color[2], lty = 2, lwd = 2)
# 
# #lines(xseq, Mu_plot[,1] + 1*Sigma_plot[,1], col = color[1], lty = 2, lwd = 2)
# #lines(xseq, Mu_plot[,1] - 1*Sigma_plot[,1], col = color[1], lty = 2, lwd = 2)
# 
# #lines(xseq, Mu_plot[,2] + 1*Sigma_plot[,2], col = color[2], lty = 2, lwd = 2)
# #lines(xseq, Mu_plot[,2] - 1*Sigma_plot[,2], col = color[2], lty = 2, lwd = 2)
# 
# plot(xseq, exp(cbind(1, poly(xseq, 3)) %*% alpha[1,]),
#      type = "l", lwd = 2, col = color[1], bty = "n", ylim = c(0, 4))
# lines(xseq, exp(cbind(1, poly(xseq, 3)) %*% alpha[2,]), col = color[2], lwd = 2)
# 
# 
# 
# 
# # MS GAMLSS ---------------------------------------------------------------
# 
# mllkGAMLSS = function(theta.star, X, Z, Z_sigma, lambda, D){
#   Gamma = LaMa::tpm(theta.star[1:2])
#   beta0 = theta.star[3:4]
#   delta = c(1, exp(theta.star[5]))
#   delta = delta / sum(delta)
#   
#   nb = ncol(Z)
#   b = theta.star[5 + 1:(2*(nb-1))]
#   beta = cbind(beta0, c(0,0),
#                matrix(b, nrow = 2, ncol = nb-1, byrow = TRUE))
#   alpha = matrix(theta.star[5 + 2*(nb-1) + 1:8], nrow = 2)
#   # sigma = exp(theta.star[5 + 2*(nb-1) + 1:2])
#   Mu = cbind(1, Z) %*% t(beta)
#   Sigma = exp(cbind(1, Z_sigma) %*% t(alpha))
#   allprobs = cbind(dnorm(X$Price, Mu[,1], Sigma[,1]),
#                    dnorm(X$Price, Mu[,2], Sigma[,2]))
#   
#   pen = t(b) %*% kronecker(diag(lambda), D) %*% b
#   -LaMa::forward(delta, Gamma, allprobs) + 0.5 * pen
# }
# 
# # design matrices
# 
# nb = 20
# d = diff(seq(min(energy$EurDol), max(energy$EurDol), length = nb-2))[1]
# knots = seq(min(energy$EurDol)-3*d, max(energy$EurDol)+3*d, by = d)
# x = seq(min(energy$EurDol), max(energy$EurDol), length = 200)
# Z = splines::splineDesign(knots, energy$EurDol, outer.ok = TRUE)
# 
# Z_sigma = poly(energy$EurDol, 3)
# 
# # penalty matrix
# L = WH:::build_D_mat(nb, 2)
# D = t(L[,-1])%*%L[,-1]
# 
# # initial parameter
# theta.star = c(rep(-3,2), 3, 6, 1,
#                rep(0, 2*(nb-1)),
#                log(1), log(1), 
#                rep(0, 8))
# 
# 
# # IPML
# maxiter = 50 # maximum number of iterations
# tol = 0.1 # relative tolerance for convergence, sufficient
# gradtol = 1e-5 # relative gradient tolerance for nlm (1e-6 is the default)
# Lambdas = matrix(NA, maxiter+1, 2)
# Lambdas[1,] = c(100, 100)
# mods = list()
# REind = matrix(5 + 1:(2*(nb-1)), nrow = 2, byrow = TRUE) # each row is the index of one RE
# 
# # updating algorithm
# for(k in 1:maxiter){
#   cat("\n\n- Iteration", k, "-")
#   if(k == 1){
#     print.level = 2
#   } else print.level = 0
#   
#   t1 = Sys.time()
#   mods[[k]] = nlm(mllkGAMLSS, theta.star, X = energy, Z=Z, Z_sigma=Z_sigma, lambda = Lambdas[k,], D=D, 
#                   iterlim = 1000, print.level = print.level, hessian = TRUE, gradtol = gradtol)
#   Sys.time()-t1
#   cat("\nEstimation time:", Sys.time()-t1)
#   cat("\nIterations:", mods[[k]]$iterations)
#   
#   theta.star = mods[[k]]$estimate # saves theta.star as starting value for next iteration
#   
#   J_p = mods[[k]]$hessian # assigning hessian
#   J_inv = MASS::ginv(J_p)
#   
#   # updating all penalty strengths
#   for(i in 1:2){
#     edoF = sum(diag( # trace
#       diag(rep(1, nrow(D))) - Lambdas[k, i] * J_inv[REind[i,], REind[i,]] %*% D
#     ))
#     penalty = t(theta.star[REind[i,]]) %*% D %*% theta.star[REind[i,]]
#     
#     Lambdas[k+1, i] = as.numeric((edoF - 1) / penalty)
#   }
#   
#   cat("\nSmoothing strengths:", round(Lambdas[k+1, ], 4))
#   
#   if(mean(abs(Lambdas[k+1,] - Lambdas[k,]) / Lambdas[k,]) < tol){
#     cat("\n\n")
#     break
#   }
# }
# 
# Lambdas = as.matrix(na.omit(Lambdas))
# par(mfrow = c(1,2))
# for(j in 1:2){
#   plot(Lambdas[,j], type = "l", lwd = 2, main = paste("lambda", j), ylim = c(0,20))
# }
# 
# 
# Gamma = LaMa::tpm(theta.star[1:2])
# beta0 = theta.star[3:4]
# delta = c(1, exp(theta.star[5]))
# delta = delta / sum(delta)
# 
# nb = ncol(Z)
# b = theta.star[5 + 1:(2*(nb-1))]
# beta = cbind(beta0, c(0,0),
#              matrix(b, nrow = 2, ncol = nb-1, byrow = TRUE))
# alpha = matrix(theta.star[5 + 2*(nb-1) + 1:8], nrow = 2)
# # sigma = exp(theta.star[5 + 2*(nb-1) + 1:2])
# Mu = cbind(1, Z) %*% t(beta)
# Sigma = exp(cbind(1, Z_sigma) %*% t(alpha))
# allprobs = cbind(dnorm(X$Price, Mu[,1], Sigma[,1]),
#                  dnorm(X$Price, Mu[,2], Sigma[,2]))
# 
# states = LaMa::viterbi(delta, Gamma, allprobs)
# 
# xseq = seq(min(energy$EurDol), max(energy$EurDol), length = 200)
# Z_plot = splines::splineDesign(knots, xseq, outer.ok = TRUE)
# 
# Mu_plot = cbind(1,Z_plot) %*% t(beta)
# Sigma_plot = exp(cbind(1, poly(xseq, 3)) %*% t(alpha))
# 
# color = c("orange", "deepskyblue")
# par(mfrow = c(1,1))
# plot(energy$EurDol, energy$Price, pch = 16, bty = "n", col = scales::alpha(color[states], 0.5))
# lines(xseq, Mu_plot[,1], col = color[1], lwd = 3)
# lines(xseq, Mu_plot[,1]+ Sigma_plot[,1], col = color[1], lty = 2, lwd = 1)
# lines(xseq, Mu_plot[,1]- Sigma_plot[,1], col = color[1], lty = 2, lwd = 1)
# 
# lines(xseq, Mu_plot[,2], col = color[2], lwd = 3)
# lines(xseq, Mu_plot[,2]+ Sigma_plot[,2], col = color[2], lty = 2, lwd = 1)
# lines(xseq, Mu_plot[,2]- Sigma_plot[,2], col = color[2], lty = 2, lwd = 1)
# 
# 
# 
# 
# # flexible modelling of the sd
# 
# mllkGAMLSS2 = function(theta.star, X, Z, lambda, D){
#   Gamma = LaMa::tpm(theta.star[1:2])
#   delta = c(1, exp(theta.star[3]))
#   delta = delta / sum(delta)
#   
#   nb = ncol(Z)
#   beta = matrix(theta.star[3 + 1:(2*(nb))], ncol = 2)
#   alpha = matrix(theta.star[3 + 2*(nb) + 1:(2*(nb))], ncol = 2)
#   
#   Mu = Z %*% beta
#   Sigma = exp(Z %*% alpha)
#   
#   allprobs = cbind(dnorm(X$Price, Mu[,1], Sigma[,1]),
#                    dnorm(X$Price, Mu[,2], Sigma[,2]))
#   
#   pen = t(as.numeric(beta)) %*% kronecker(diag(lambda[1:2]), D) %*% as.numeric(beta) +
#     t(as.numeric(alpha)) %*% kronecker(diag(lambda[3:4]), D) %*% as.numeric(alpha)
#   
#   # pen = t(as.numeric(beta[-1,])) %*% kronecker(diag(lambda[1:2]), D) %*% as.numeric(beta[-1,]) + 
#   #   t(as.numeric(alpha[-1,])) %*% kronecker(diag(lambda[3:4]), D) %*% as.numeric(alpha[-1,])
#   
#   -LaMa::forward(delta, Gamma, allprobs) + 0.5 * pen
# }
# 
# # design matrix
# nb = 15
# d = diff(seq(min(energy$EurDol), max(energy$EurDol), length = nb-2))[1]
# knots = seq(min(energy$EurDol)-3*d, max(energy$EurDol)+3*d, by = d)
# x = seq(min(energy$EurDol), max(energy$EurDol), length = 200)
# Z = splines::splineDesign(knots, energy$EurDol, outer.ok = TRUE)
# 
# # nb = 15
# # gam_setup = gam(Price ~ s(EurDol, k = nb, bs = "cs"), data = cbind(dummy = 1, energy), fit = FALSE)
# # Z = gam_setup$X
# # D = gam_setup$S[[1]]
# 
# # penalty matrix
# L = WH:::build_D_mat(nb, 2)
# D = t(L)%*%L
# 
# # IPML
# maxiter = 50 # maximum number of iterations
# tol = 0.1 # relative tolerance for convergence, sufficient
# gradtol = 1e-5 # relative gradient tolerance for nlm (1e-6 is the default)
# Lambdas = matrix(NA, maxiter+1, 4)
# Lambdas[1,] = rep(100, 4)
# mods = list()
# REind = matrix(3 + 1:(4*(nb)), nrow = 4, byrow = TRUE) # each row is the index of one RE
# # REind = matrix(c(4+1:(nb-1),
# #                  19+1:(nb-1),
# #                  34+1:(nb-1),
# #                  49+1:(nb-1)), nrow = 4, byrow = TRUE) # each row is the index of one RE
# 
# # initial parameter
# theta.star = c(rep(-2, 2), 1,
#                rep(3, nb), rep(5, nb),
#                rep(0, 2*nb))
# # theta.star = c(rep(-2, 2), 1,
# #                3, rep(0, (nb-1)), 5, rep(0, (nb-1)),
# #                rep(0, 2*nb))
# 
# # updating algorithm
# T1 = Sys.time()
# for(k in 1:maxiter){
#   cat("\n\n- Iteration", k, "-")
#   if(k == 1){
#     print.level = 2
#   } else print.level = 0
#   
#   t1 = Sys.time()
#   mods[[k]] = nlm(mllkGAMLSS2, theta.star, X = energy, Z=Z, lambda = Lambdas[k,], D=D, 
#                   iterlim = 1000, print.level = print.level, hessian = TRUE, gradtol = gradtol)
#   Sys.time()-t1
#   cat("\nEstimation time:", Sys.time()-t1)
#   cat("\nIterations:", mods[[k]]$iterations)
#   
#   theta.star = mods[[k]]$estimate # saves theta.star as starting value for next iteration
#   
#   J_p = mods[[k]]$hessian # assigning hessian
#   J_inv = MASS::ginv(J_p)
#   
#   # updating all penalty strengths
#   for(i in 1:4){
#     edoF = sum(diag( # trace
#       diag(rep(1, nrow(D))) - Lambdas[k, i] * J_inv[REind[i,], REind[i,]] %*% D
#     ))
#     penalty = t(theta.star[REind[i,]]) %*% D %*% theta.star[REind[i,]]
#     Lambdas[k+1, i] = as.numeric((edoF - 1) / penalty)
#   }
#   
#   cat("\nSmoothing strengths:", round(Lambdas[k+1, ], 4))
#   
#   if(mean(abs(Lambdas[k+1,] - Lambdas[k,]) / Lambdas[k,]) < tol){
#     cat("\n\n")
#     break
#   }
# }
# Sys.time()-T1
# 
# Lambdas = as.matrix(na.omit(Lambdas))
# par(mfrow = c(1,4))
# for(j in 1:4){
#   plot(Lambdas[,j], type = "l", lwd = 2, main = paste("lambda", j))
# }
# 
# Gamma = LaMa::tpm(theta.star[1:2])
# delta = c(1, exp(theta.star[3]))
# delta = delta / sum(delta)
# 
# nb = ncol(Z)
# beta = matrix(theta.star[3 + 1:(2*(nb))], ncol = 2)
# alpha = matrix(theta.star[3 + 2*(nb) + 1:(2*(nb))], ncol = 2)
# 
# Mu = Z %*% beta
# Sigma = exp(Z %*% alpha)
# 
# allprobs = cbind(dnorm(energy$Price, Mu[,1], Sigma[,1]),
#                  dnorm(energy$Price, Mu[,2], Sigma[,2]))
# 
# states = LaMa::viterbi(delta, Gamma, allprobs)
# 
# xseq = seq(min(energy$EurDol), max(energy$EurDol), length = 200)
# Z_plot = splines::splineDesign(knots, xseq, outer.ok = TRUE)
# # Z_plot = gam(Price ~ s(EurDol, k = nb, bs = "cs"), data = data.frame(dummy = 1, Price = 1, EurDol = xseq), fit = FALSE)$X
# Mu_plot = Z_plot %*% beta
# Sigma_plot = exp(Z_plot %*% alpha)
# 
# 
# color = c("orange", "deepskyblue")
# par(mfrow = c(1,1))
# plot(energy$EurDol, energy$Price, pch = 16, bty = "n", col = scales::alpha(color[states], 0.5),
#      xlab = "EurDol", ylab = "Price")
# for(j in 1:2) lines(xseq, Mu_plot[,j], col = color[j], lwd = 3)
# 
# qseq = qnorm(seq(0.5, 0.975, length = 10))
# for(i in qseq){
#   for(j in 1:2){
#     lines(xseq, Mu_plot[,j]+ i*Sigma_plot[,j], col = scales::alpha(color[j], 0.5), lwd = 1, lty = 2)
#     lines(xseq, Mu_plot[,j]- i*Sigma_plot[,j], col = scales::alpha(color[j], 0.5), lwd = 1, lty = 2)
#   }
# }
# 
# splot(xseq, Sigma_plot[,1], type = "l", lwd = 3, col = color[1], bty = "n", ylim = c(0,3))
# lines(xseq, Sigma_plot[,2], type = "l", lwd = 3, col = color[2])
# 
# 
# plot(energy$Price, type = "l", bty = "n", col = "gray")
# points(energy$Price, pch = 20, col = scales::alpha(color[states], 0.5))
# 
# 
# 



# With mgcv building blocks -----------------------------------------------

data(energy, package = "MSwM")

mllk_mgcv = function(theta.star, X, Z, lambda, S){
  Gamma = LaMa::tpm(theta.star[1:2])
  delta = c(1, exp(theta.star[3]))
  delta = delta / sum(delta)
  
  nb = ncol(Z)
  beta = matrix(theta.star[3 + 1:(2*(nb))], ncol = 2)
  alpha = matrix(theta.star[3 + 2*(nb) + 1:(2*(nb))], ncol = 2)
  
  Mu = Z %*% beta
  Sigma = exp(Z %*% alpha)
  
  allprobs = cbind(dnorm(X$Price, Mu[,1], Sigma[,1]),
                   dnorm(X$Price, Mu[,2], Sigma[,2]))
  
  pen = t(as.numeric(beta[-1,])) %*% kronecker(diag(lambda[1:2]), S) %*% as.numeric(beta[-1,]) +
    t(as.numeric(alpha[-1,])) %*% kronecker(diag(lambda[3:4]), S) %*% as.numeric(alpha[-1,])
  
  -LaMa::forward(delta, Gamma, allprobs) + 0.5 * pen
}

# design and penalty matrix
nb = 15
gam_setup = mgcv::gam(Price ~ s(EurDol, k = nb, bs = "ps"), data = cbind(dummy = 1, energy), fit = FALSE)
Z = gam_setup$X
S = gam_setup$S[[1]]

# IPML
maxiter = 50 # maximum number of iterations
tol = 0.05 # relative tolerance for convergence, sufficient
gradtol = 1e-5 # relative gradient tolerance for nlm (1e-6 is the default)
Lambdas = matrix(NA, maxiter+1, 4)
Lambdas[1,] = c(rep(100, 2), rep(1000, 2))
mods = list()
REind = matrix(3 + 1:(4*(nb)), nrow = 4, byrow = TRUE)[,-1] # each row is the index of one RE

# initial parameter
theta.star = c(rep(-2, 2), 1,
               3, rep(0, nb-1), 5, rep(0, nb-1),
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

xseq = seq(min(energy$EurDol), max(energy$EurDol), length = 200)
Z_plot = mgcv::gam(Price ~ s(EurDol, k = nb, bs = "ps"), data = data.frame(dummy = 1, Price=1, EurDol = xseq), fit = FALSE)$X

Mu_plot = Z_plot %*% beta
Sigma_plot = exp(Z_plot %*% alpha)

color = c("orange", "deepskyblue")
par(mfrow = c(1,1))
plot(energy$EurDol, energy$Price, pch = 16, bty = "n", col = scales::alpha(color[states], 0.4),
     xlab = "EurDol", ylab = "Price")
for(j in 1:2) lines(xseq, Mu_plot[,j], col = color[j], lwd = 3)

qseq = qnorm(seq(0.5, 0.95, length = 10))
for(i in qseq){
  for(j in 1:2){
    lines(xseq, Mu_plot[,j]+ i*Sigma_plot[,j], col = scales::alpha(color[j], 0.5), lwd = 1, lty = 2)
    lines(xseq, Mu_plot[,j]- i*Sigma_plot[,j], col = scales::alpha(color[j], 0.5), lwd = 1, lty = 2)
  }
}

plot(energy$Price, col = color[states])





# hmmTMB ------------------------------------------------------------------

library(hmmTMB)

data(energy, package = "MSwM")

hid <- MarkovChain$new(data = energy, n_states = 2)
# List of observation distributions
dists <- list(Price = "norm")
# List of initial parameters
par0 <- list(Price = list(mean = c(3, 6), sd = c(1, 1)))
# List of formulas
# f <- list(Price = list(mean = ~ s(EurDol, k = 15, bs = "ps"),
#                        # sd = ~ poly(EurDol, 3)))
#                        sd = ~ s(EurDol, k = 15, bs = "ps")))
f <- list(Price = list(mean = ~ s(Oil, k = 10, bs = "ps"),
                       sd = ~ s(Oil, k = 10, bs = "ps")))
# Create observation model
obs <- Observation$new(data = energy,
                       n_states = 2,
                       dists = dists,
                       par = par0,
                       formulas = f)

hmm <- HMM$new(hid = hid, obs = obs)

t1 = Sys.time()
hmm$fit()
Sys.time()-t1

energy$viterbi <- factor(paste0("State ", hmm$viterbi()))

hmm$plot(what = "obspar", var = "Oil", i = "Price.mean") +
  geom_point(aes(x = Oil, y = Price, fill = viterbi, col = viterbi),
             data = energy, alpha = 0.3) +
  theme(legend.position = "none")

hmm$plot(what = "obspar", var = "EurDol", i = "Price.sd") +
  theme(legend.position = c(0.3, 0.7))




# Replicating the Application in Timos paper ------------------------------

# energy prices explained by oil prices

data(energy, package = "MSwM")

mllk_mgcv = function(theta.star, X, Z, lambda, S){
  Gamma = LaMa::tpm(theta.star[1:2])
  delta = c(1, exp(theta.star[3]))
  delta = delta / sum(delta)
  
  nb = ncol(Z)
  beta = matrix(theta.star[3 + 1:(2*(nb))], ncol = 2)
  alpha = matrix(theta.star[3 + 2*(nb) + 1:(2*(nb))], ncol = 2)
  
  Mu = Z %*% beta
  Sigma = exp(Z %*% alpha)
  
  allprobs = cbind(dnorm(X$Price, Mu[,1], Sigma[,1]),
                   dnorm(X$Price, Mu[,2], Sigma[,2]))
  
  pen = t(as.numeric(beta[-1,])) %*% kronecker(diag(lambda[1:2]), S) %*% as.numeric(beta[-1,]) +
    t(as.numeric(alpha[-1,])) %*% kronecker(diag(lambda[3:4]), S) %*% as.numeric(alpha[-1,])

  -LaMa::forward(delta, Gamma, allprobs) + 0.5 * pen
}

# design and penalty matrix
# nb = 15
# gam_setup = mgcv::gam(Price ~ s(Oil, k = nb, bs = "ps"), data = cbind(dummy = 1, energy), fit = FALSE)
# Z = gam_setup$X
# S = gam_setup$S[[1]]

# make my own design and penalty matrix
nb = 20
d = diff(seq(min(energy$Oil), max(energy$Oil), length = (nb-1)-2))[1]
knots = seq(min(energy$Oil)-3*d, max(energy$Oil)+3*d, by = d)
x = seq(min(energy$Oil), max(energy$Oil), length = 200)
Z = splines::splineDesign(knots, energy$Oil, outer.ok = TRUE)
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
Lambdas[1,] = c(rep(100, 2), rep(10000, 2))
mods = list()
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

xseq = seq(min(energy$Oil), max(energy$Oil), length = 200)
# Z_plot = mgcv::gam(Price ~ s(Oil, k = nb, bs = "ps"), data = data.frame(dummy = 1, Price=1, Oil = xseq), fit = FALSE)$X
Z_plot = cbind(1, splines::splineDesign(knots, xseq, outer.ok = TRUE))

Mu_plot = Z_plot %*% beta
Sigma_plot = exp(Z_plot %*% alpha)

color = c("orange", "deepskyblue")
par(mfrow = c(1,1))
plot(energy$Oil, energy$Price, pch = 16, bty = "n", col = scales::alpha(color[states], 0.25),
     xlab = "Oil price", ylab = "Price")
for(j in 1:2) lines(xseq, Mu_plot[,j], col = color[j], lwd = 3)

qseq = qnorm(seq(0.5, 0.95, length = 8))
for(i in qseq){
  for(j in 1:2){
    lines(xseq, Mu_plot[,j] + i*Sigma_plot[,j], col = scales::alpha(color[j], 0.7), lwd = 1, lty = 2)
    lines(xseq, Mu_plot[,j] - i*Sigma_plot[,j], col = scales::alpha(color[j], 0.7), lwd = 1, lty = 2)
  }
}

plot(energy$Price, col = color[states])



# Testing with Gas --------------------------------------------------------

# make my own design and penalty matrix
nb = 25
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

# IPML
maxiter = 50 # maximum number of iterations
tol = 0.01 # relative tolerance for convergence, sufficient
gradtol = 1e-6 # relative gradient tolerance for nlm (1e-6 is the default)
Lambdas = matrix(NA, maxiter+1, 4)
Lambdas[1,] = c(rep(100, 2), rep(1000, 2))
mods = list()
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

qseq = qnorm(seq(0.5, 0.95, length = 5))
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








# New application ---------------------------------------------------------
library(fHMM)

bit = fHMM::download_data("BTC-USD", from = "2014-09-17")
eurdol = fHMM::download_data("EURUSD=X", from = "2007-02-01")
gold = fHMM::download_data("GC=F", from = "2000-08-30")
usd = fHMM::download_data("DX-Y.NYB", from = "2000-08-30")

data = data.frame(Date = gold$Date, Gold = gold$Close)

library(tidyverse)

data = data %>% left_join(usd[c("Date", "Close")], by = "Date")

colnames(data)[3] = "USD"

## impute NAs in usd column
data$EurDol[is.na(data$EurDol)] = data$EurDol[which(is.na(data$EurDol)) - 1]


nvda = fHMM::download_data("NVDA", from = "2004-08-19")
goog = fHMM::download_data("GOOG", from = "2004-08-19")

nvda$return = c(NA, diff(log(nvda$Close)))

data = data.frame(Date = nvda$Date, Price = nvda$Close, Google = goog$Close)

plot(data$Google, data$Price, pch = 20, bty = "n")

par(mfrow = c(2,1))
plot(data$Gold, type = "l")
plot(data$USD, type = "l")

par(mfrow = c(1,1))
plot(data$USD, data$Gold, pch = 20)






