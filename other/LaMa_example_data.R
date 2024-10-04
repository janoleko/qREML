library(LaMa)


# Elephant example --------------------------------------------------------

## elephant data
# steps and turns
# complex relationship between time of day and transition probabilities


## parameters for the state-dependent process
mu = c(0.3, 2.5)
sigma = c(0.2, 1.5)
kappa = c(0.1, 1.5)

## parameter matrix for state process
beta = matrix(c(-2,-2,0.2,0.4,2,2,3,1,-2,4.5,7.5,-6,5,-3,0.5,2,-5,0.5,-5,10), 
              nrow = 2)

## design matrix for state process
knots = seq(0, 24, length.out = 11)
modmat = make_matrices(~ s(tod, bs = "cp", k = 10),
                       data = data.frame(tod = 1:24),
                       knots = list(tod = knots))
## tpm
Gamma = tpm_g(modmat$Z, beta)

## simulate data
n = 1e4
s = rep(NA, n)

tod = rep(1:24, n/20)
tod = tod[(length(tod)-n+1):length(tod)]

## initial distribution
delta = stationary_p(Gamma, t = tod[1])

## simulating state process
s[1] = sample(1:2, 1, prob = delta)
for(t in 2:n) {
  s[t] = sample(1:2, 1, prob = Gamma[s[t-1],,tod[t]])
}

## simulating state-dependent process
step = rgamma2(n, mu[s], sigma[s])

angle = rep(NA, n)
for(t in 1:n) angle[t] = CircStats::rvm(1, pi, kappa[s[t]]) - pi

elephant = data.frame(tod = tod, step = step, angle = angle, state = s)
elephant$step[n] = NA
elephant$angle[c(1,n)] = NA

saveRDS(elephant, "~/Desktop/elephant_sim.rds")



f <- function(...) {
  # Gather all arguments into a list
  args_list <- list(...)
  
  # Capture the names of the arguments
  arg_names <- sapply(substitute(list(...))[-1], deparse)
  
  # Combine names and the list of arguments into a named list
  names(args_list) <- arg_names
  
  return(args_list)
}
x = 1
y = 1

f(x,y)

?penalty


re_coef = list(x,y)
re_coef
