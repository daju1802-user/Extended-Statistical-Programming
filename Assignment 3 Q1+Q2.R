##This function implements the initial part of the Covid-19 deconvolution model.
##The goal is to construct the design matrix required for estimation
##the design matrices needed to estimate the daily infection curve f(t) 
## Observed deaths t_death and infection-to-death distribution π(d) are given.
## The model assumes a Poisson likelihood
## The infection curve f(t) is assumed to be smooth, and represented using a B-spline.
## The function performs the following:
## 1. Defines the infection-to-death probability distribution π(d)
##    for durations d = 1:80 days using a lognormal distribution
## 2. Constructs the B-spline basis matrix tilde_X for f(t):
##    K + 4 evenly spaced knots are used.
##    The middle K-2 knots cover the interval over which f(t) is evaluated.
## 3. Constructs the design matrix X for mu_i:
##    Each row corresponds to one observed death day t_i.
##    Each column corresponds to a B-spline coefficient beta_k.
##    Each mu_i is computed as a weighted sum of the basis functions
##    evaluated at the infection times t_i - j, weighted by p(j).
## 4. Computes the smoothing penalty matrix S:
##    S = t(D) %*% D, where D is the second-difference matrix.
##    Penalizes roughness in beta, encouraging smooth f(t).

data <- read.table("engcov.txt", header = TRUE)
t_death <- data$julian
y <- data$nhs

create_design_matrices <- function(t_death, K = 80) {
  ## Define infection-to-death probability distribution
  d <- 1:80
  edur <- 3.151
  sdur <- 0.469
  pd <- dlnorm(d, edur, sdur)
  pd <- pd / sum(pd)  # Normalize to probability distribution
  
  n <- length(t_death)  # Number of death observations
  ## Determine time range for f(t)
  ## f(t) needs to cover infections that could lead to observed deaths
  t_min <- min(t_death) - 80  # Go back maximum duration before first death
  t_max <- max(t_death)       # Up to last death day
  t_f <- t_min:t_max          # All time points for infection function
  
  ## Create B-spline basis matrix tilde_X
  ## Create K+4 evenly spaced knots covering the range of f(t)
  knots <- seq(t_min, t_max, length.out = K + 4)
  tilde_X <- splineDesign(knots = knots, x = t_f, outer.ok = TRUE)
  
  ## Create design matrix X for death model
  # Create matrix of all possible (i,j) combinations
  i_vec <- rep(1:n, times = pmin(29 + 1:n, 80))
  j_vec <- unlist(lapply(1:n, function(i) 1:pmin(29 + i, 80))) 
  
  # Corresponding infection times
  infection_times <- t_death[i_vec] - j_vec
  t_f_idx <- match(infection_times, t_f)  # Find the row index corresponding to these infection times in t_f
  # Weighted tilde_X rows by probability
  weighted_rows <- tilde_X[t_f_idx, , drop = FALSE] * pd[j_vec]
  
  # Add the weighted basis function rows belonging to the same death date i
  X <- rowsum(weighted_rows, group = i_vec)
  
  ## Penalty matrix S 
  D <- diff(diag(K), differences = 2)
  S <- crossprod(D)
  
  ## Return all matrices
  return(list(
    tilde_X = tilde_X,
    X = X,
    S = S,
    t_f = t_f,
    pd = pd
  ))
}


##This function computes the penalized negative log-likelihood and its gradient 
##for the Poisson deconvolution model used to estimate the daily infection curve f(t).
##The model assumes that observed daily deaths y_i follow a Poisson distribution
##A smoothness penalty on beta is included to control overfitting.
##The functions are designed for use with 'optim' for numerical optimization.

penalized_nll <- function(gamma, y, X, S, lambda) {
  beta <- exp(gamma)  # ensure positivity of beta
  mu <- as.vector(X %*% beta)   
  mu[mu <= 1e-12] <- 1e-12   # avoid log(0)
  
  # negative log-likelihood
  nll <- sum(mu - y * log(mu))
  
  # add smoothness penalty term
  penalty <- 0.5 * lambda * crossprod(beta, S %*% beta)
  return(nll + penalty)
}


penalized_grad <- function(gamma, y, X, S, lambda) {
  beta <- exp(gamma)
  mu <- as.vector(X %*% beta)
  mu[mu <= 1e-12] <- 1e-12
  
  # gradient of log-likelihood part
  grad_ll <- t(X) %*% (1 - y / mu)
  
  # gradient of penalty part
  grad_pen <- lambda * (S %*% beta)
  
  # chain rule
  grad <- as.vector(beta * (grad_ll + grad_pen))
  return(grad)
}

# ---- Finite differencing test ----
# Choose initial gamma
gamma0 <- rep(0, ncol(mats$X))
lambda <- 5e-5

# Analytical gradient
g_analytical <- penalized_grad(gamma0, y, mats$X, mats$S, lambda)

# Numerical gradient (finite differencing)
eps <- 1e-6
g_numeric <- sapply(1:length(gamma0), function(k) {
  gp <- gamma0; gp[k] <- gp[k] + eps
  gm <- gamma0; gm[k] <- gm[k] - eps
  (penalized_nll(gp, y, mats$X, mats$S, lambda) -
      penalized_nll(gm, y, mats$X, mats$S, lambda)) / (2 * eps)
})

# Compare
diff <- max(abs(g_analytical - g_numeric))
cat("Maximum absolute difference between analytic and numeric gradients:", diff, "\n")

