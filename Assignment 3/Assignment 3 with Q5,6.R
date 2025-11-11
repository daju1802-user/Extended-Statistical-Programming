library(splines)
data <- read.table("engcov.txt", header = TRUE)
t <- data$julian  # Day of year
y <- data$nhs
n <- length(y)  # Number of observations


#============================================
# Step 1: Create model matrices X_tilde, X, and penalty matrix S
#============================================
d <- 1:80
edur <- 3.151
sdur <- 0.469
pd <- dlnorm(d, edur, sdur)
pd <- pd / sum(pd)  # Normalize to probability distribution
K <- 80  # Number of basis functions

build_matrices <- function(t,k,pd){
  ## Construct design matrices for the Poisson deconvolution model.
  ## This function prepares the matrices required to estimate the infection curve f(t),
  ## which is related to observed deaths through a convolution with the infection-to-death distribution.
  lower <- min(t) - 30
  upper <- max(t)
  
  # create a sequence of interior knots evenly spaced over the time range
  mid_squence <- seq(lower, upper, length.out = k-2)
  gap <- mid_squence[2] - mid_squence[1]
  
  # extend the knots 3 steps beyond both ends to handle cubic splines
  knot <- c(
    seq(from = lower - 3*gap, to = lower - gap, length.out = 3),
    mid_squence,
    seq(from = upper + gap, to = upper + 3*gap, length.out = 3)
  )
  
  # construct the B-spline basis matrix between lower and upper
  X_tilde <- splineDesign(knot,lower:upper,outer.ok = TRUE) #outer.ok = TRUE
  r <- nrow(X_tilde) # number of time points in f(t) grid
  
  X <- matrix(0, n, ncol(X_tilde))
  # build X by summing contributions of past infections to each day's deaths
  for (i in 1:n) {
    max_j <- min(29 + i, length(pd))
    # compute indices of infection times corresponding to death day
    idx <- t[i] - (1:max_j) - lower + 1
    valid <- which(idx >= 1 & idx <= r)
    
    # if valid infections exist, add weighted basis function contributions
    if (length(valid) > 0) {
      X[i, ] <- colSums(X_tilde[idx[valid], , drop = FALSE] * pd[valid])
    }
  }
  
  S <- crossprod(diff(diag(k),diff=2)) # construct penalty matrix
  return(list(X_tilde=X_tilde,S=S,X=X,t_infection = lower:upper))
}

mats <- build_matrices(t,k = 80, pd)
X <- mats$X
X_tilde <- mats$X_tilde
S <- mats$S
head(X)

#============================================
# Step 2: Likelihood and gradient functions
#============================================
pnll <- function(gamma, y, X, S, lambda) {
  ## Computes the penalized negative log-likelihood and its gradient 
  ## for the Poisson deconvolution model used to estimate the daily infection curve f(t).
  ## The model assumes that observed daily deaths y_i follow a Poisson distribution
  ## A smoothness penalty on beta is included to control overfitting.
  ## The functions are designed for use with 'optim' for numerical optimization.
  beta <- exp(gamma)  # ensure positivity of beta
  mu <- as.vector(X %*% beta)  # Compute expected mean number of deaths under the model
  
  # Poisson log-likelihood 
  ll <- sum(dpois(y,mu,log = TRUE))
  
  # Add smoothness penalty term
  penalty <- lambda * sum(beta * (S %*% beta)) / 2
  # Return negative penalized log-likelihood
  -(ll - penalty)
}


pnll_grad <- function(gamma, y, X, S, lambda) {
  ## Compute the gradient of the penalized negative log-likelihood
  ## for the Poisson deconvolution model estimating daily infections.
  beta <- exp(gamma) # Transform to positive scale
  mu <- as.vector(X %*% beta) # Expected deaths
  
  # Gradient of log-likelihood w.r.t. gamma
  # dl/dgamma = F, where F = diag(y/mu - 1) * X * diag(beta)
  residual <- y / mu - 1  # residuals
  F <- (residual * X) * rep(beta, each = length(y))  # Element-wise multiplication
  dll_dgamma <- colSums(F)
  
  # Gradient of penalty w.r.t. gamma
  # dP/dgamma = diag(beta) * S * beta
  dP_dgamma <- beta * (S %*% beta)
  
  # Return negative gradient
  -(dll_dgamma - lambda * dP_dgamma)
}

# ---- Finite differencing test ----
test_gradient <- function(gamma, y, X, S, lambda) {
  ## Test the correctness of the analytic gradient function by finite differencing.
  # Analytic gradient from the implemented function
  grad_analytic <- pnll_grad(gamma, y, X, S, lambda)
  
  # Finite difference gradient
  eps <- 1e-7
  grad_fd <- numeric(length(gamma))
  for (i in 1:length(gamma)) {  # loop over each parameter in gamma
    gamma_plus <- gamma
    gamma_plus[i] <- gamma[i] + eps
    gamma_minus <- gamma
    gamma_minus[i] <- gamma[i] - eps
    # central difference formula
    grad_fd[i] <- (pnll(gamma_plus, y, X, S, lambda) - 
                     pnll(gamma_minus, y, X, S, lambda)) / (2 * eps)
  }
  
  # compare analytic vs numerical gradients
  max_diff <- max(abs(grad_analytic - grad_fd))
  cat("Max difference between analytic and finite difference gradient:", max_diff, "\n")
  
  invisible(max_diff < 1e-4)  # Should be very small
}

# Initialize gamma with small values (beta close to 1)
gamma_init <- rep(0, K)
# Test the gradient function
test_gradient(gamma_init, y, X, S, lambda = 5e-5)



#Q3
## ---- Optimize penalized likelihood ----
gamma_init <- rep(0, K)
lambda <- 5e-5
t_infection = mats$t_infection
# Optimize
fit_test <- optim(gamma_init, pnll, pnll_grad, 
                  lambda = lambda, X = X, y = y, S = S,
                  method = "BFGS", control = list(maxit = 1000))


## ---- Extract fitted parameters ----
gamma_hat <- fit_test$par
beta_hat <- exp(gamma_hat)

# Calculate fitted values
mu_hat <- as.vector(X %*% beta_hat) # fitted daily deaths
f_hat <- as.vector(X_tilde %*% beta_hat) # estimated infection curve


## ---- Plot the results ----
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))
plot(t, y, type = "p", pch = 16, col = "black",
     xlab = "Day of year 2020", ylab = "Count",
     main = "COVID-19 Deaths and Inferred Infections",
     ylim = c(0, max(y, f_hat) * 1.1))
lines(t, mu_hat, col = "red", lwd = 2)
lines(t_infection, f_hat, col = "blue", lwd = 2)
legend("topright", 
       legend = c("Observed deaths", "Fitted deaths", "Inferred infections"),
       col = c("black", "red", "blue"), 
       pch = c(16, NA, NA), lty = c(NA, 1, 1), lwd = 2)


## ---- Q4: Grid search for optimal lambda  ----

# Function to calculate BIC
calc_bic <- function(gamma, lambda, X, y, S) {
  beta <- exp(gamma)
  mu <- as.vector(X %*% beta)
  
  # Log-likelihood (without penalty)
  ll <- sum(dpois(y,mu,log = TRUE))
  
  # Calculate Hessian matrices for EDF
  # H_lambda = X^T * W * X + lambda * S, where W = diag(y/mu^2)
  W <- as.vector(y / mu^2)
  H_0 = crossprod(X, X * W)
  H_lambda <- H_0 + lambda * S
  
  #trace(A,B) = sum(A * t(B)) A and B are square matrices
  #H_0 <-  t(X) %*% (X * w)
  #Edf <- sum(diag(solve(H_lambda,H_0)))
  #Using Cholesky decomposition will make solving this equation more efficient,
  #but first, it's necessary to determine if the matrix satisfies the conditions for using this decomposition.
  cholH <- try(chol(H_lambda), silent = TRUE)
  if (inherits(cholH, "try-error")) {
    # Rollback plan(traditional way)
    EDF <- sum(diag(solve(H_lambda, H_0)))
  } else {
    # Avoid matrix inversion by using a two-step triangular decomposition method.
    Z <- backsolve(cholH, forwardsolve(t(cholH), H_0))
    EDF <- sum(diag(Z))
  }
  # BIC
  BIC <- -2 * ll + log(n) * EDF
  list(BIC = BIC, EDF = EDF, ll = ll)
}


# Grid search over log(lambda)
log_lambda <- seq(-13, -7, length = 50)
lambda_seq <- exp(log_lambda)
bic_values <- numeric(length(lambda_seq))
edf_values <- numeric(length(lambda_seq))

# Start from previous solution and update for each lambda
gamma_current <- gamma_init

for (i in 1:length(lambda_seq)) {
  lambda_i <- lambda_seq[i]
  
  # Optimize for this lambda
  fit_i <- optim(gamma_current, pnll, pnll_grad,
                 lambda = lambda_i, X = X, y = y, S = S,
                 method = "BFGS", control = list(maxit = 1000))
  
  gamma_current <- fit_i$par  # Use as starting point for next
  
  # Calculate BIC
  bic_info <- calc_bic(gamma_current, lambda_i, X, y, S)
  bic_values[i] <- bic_info$BIC
  edf_values[i] <- bic_info$EDF
  
}
# Find optimal lambda
idx_opt <- which.min(bic_values)
lambda_opt <- lambda_seq[idx_opt]
cat("\nOptimal lambda:", lambda_opt, "\n")
cat("Optimal BIC:", bic_values[idx_opt], "\n")
cat("Optimal EDF:", edf_values[idx_opt], "\n")

# Refit with optimal lambda
fit_opt <- optim(gamma_init, pnll, pnll_grad,
                 lambda = lambda_opt, X = X, y = y, S = S,
                 method = "BFGS", control = list(maxit = 1000))
#Back to function to find the optimal gamma and beta
gamma_opt <- fit_opt$par
beta_opt <- exp(gamma_opt)


# ==================== Q5: Non-parametric Bootstrap ====================

# Modify objective function to support weights
pnll_weighted <- function(gamma, y, X, S, lambda, weights) {
  ## Penalized negative log-likelihood function with weights
  beta <- exp(gamma)
  mu <- as.vector(X %*% beta)
  
  # Weighted Poisson log-likelihood
  ll <- sum(weights * dpois(y, mu, log = TRUE))
  
  # Smoothness penalty term
  penalty <- lambda * sum(beta * (S %*% beta)) / 2
  
  # Return negative penalized log-likelihood
  -(ll - penalty)
}

pnll_grad_weighted <- function(gamma, y, X, S, lambda, weights) {
  ## Gradient function with weights
  beta <- exp(gamma)
  mu <- as.vector(X %*% beta)
  
  # Weighted likelihood gradient
  residual <- y / mu - 1
  F_matrix <- (residual * X) * rep(beta, each = length(y))
  dll_dgamma <- colSums(weights * F_matrix)
  
  # Penalty gradient
  dP_dgamma <- beta * (S %*% beta)
  
  # Return negative gradient
  -(dll_dgamma - lambda * dP_dgamma)
}

# Use optimal lambda for bootstrap
lambda_fixed <- lambda_opt  # Use optimal lambda found in Q4
n_bootstrap <- 200  # Number of bootstrap replicates
n_obs <- length(y)  # Number of observations

# Store bootstrap results
bootstrap_f <- matrix(0, nrow = nrow(X_tilde), ncol = n_bootstrap)
bootstrap_beta <- matrix(0, nrow = K, ncol = n_bootstrap)

# Progress bar setup
cat("Starting bootstrap process...\n")
pb <- txtProgressBar(min = 0, max = n_bootstrap, style = 3)

# Use optimal solution as starting value for better convergence
gamma_start <- gamma_opt

# Bootstrap loop
for (b in 1:n_bootstrap) {
  # Generate bootstrap weights
  bootstrap_weights <- tabulate(sample(n_obs, replace = TRUE), n_obs)
  
  # Optimize using weighted functions
  fit_bootstrap <- optim(
    par = gamma_start,
    fn = pnll_weighted,
    gr = pnll_grad_weighted,
    y = y, X = X, S = S, lambda = lambda_fixed,
    weights = bootstrap_weights,
    method = "BFGS",
    control = list(maxit = 1000)
  )
  
  # Store results
  gamma_b <- fit_bootstrap$par
  beta_b <- exp(gamma_b)
  f_b <- as.vector(X_tilde %*% beta_b)
  
  bootstrap_beta[, b] <- beta_b
  bootstrap_f[, b] <- f_b
  
  # Update progress bar
  setTxtProgressBar(pb, b)
  
  # Update starting value every 50 iterations (avoid divergence)
  if (b %% 50 == 0) {
    gamma_start <- gamma_b
  }
}

close(pb)
cat("Bootstrap completed!\n")

# Calculate bootstrap confidence intervals
f_hat_opt <- as.vector(X_tilde %*% beta_opt)  # Optimal fitted infection curve

# Calculate bootstrap quantiles for each time point
f_bootstrap_ci <- apply(bootstrap_f, 1, function(x) {
  quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
})

f_lower <- f_bootstrap_ci[1, ]  # 95% CI lower bound
f_upper <- f_bootstrap_ci[2, ]  # 95% CI upper bound

# ==================== Q6: Final Plotting ====================

# Set graphics parameters
par(mfrow = c(2, 1), mar = c(4, 4, 2, 2))

# Plot 1: Death data and model fit
plot(t, y, type = "p", pch = 16, col = "black", cex = 0.7,
     main = "Daily Death Data and Model Fit",
     xlab = "Day of year 2020", ylab = "Daily Deaths",
     ylim = c(0, max(y) * 1.1))

# Add model fit line
mu_opt <- as.vector(X %*% beta_opt)
lines(t, mu_opt, col = "red", lwd = 2)

legend("topright", 
       legend = c("Observed deaths", "Fitted deaths"),
       col = c("black", "red"), 
       pch = c(16, NA), lwd = c(NA, 2),
       bty = "n")

# Plot 2: Estimated daily infection rate with confidence intervals
infection_dates <- mats$t_infection
plot(infection_dates, f_hat_opt, type = "l", col = "blue", lwd = 2,
     main = "Estimated Daily Infection Rate with 95% Confidence Intervals",
     xlab = "Day of year 2020", ylab = "Estimated Infections",
     ylim = c(0, max(f_upper) * 1.1))

# Add confidence intervals
polygon(c(infection_dates, rev(infection_dates)),
        c(f_lower, rev(f_upper)),
        col = rgb(0, 0, 1, 0.2), border = NA)

# Re-add central line for visibility
lines(infection_dates, f_hat_opt, col = "blue", lwd = 2)

# Add vertical line marking start of death data
death_start <- min(t)
abline(v = death_start, lty = 2, col = "red")
text(death_start, max(f_upper) * 0.9, 
     paste("Death Data Starts\n(Infection +", round(exp(edur)), "days)", sep = ""),
     pos = 2, cex = 0.8, col = "red")

legend("topright",
       legend = c("Inferred infections", "95% Confidence Interval"),
       col = c("blue", rgb(0, 0, 1, 0.3)),
       lwd = c(2, 8), lty = c(1, 1),
       bty = "n")

# Reset graphics parameters
par(mfrow = c(1, 1))

# Output bootstrap results summary
cat("\nBootstrap Results Summary:\n")
cat("Bootstrap replicates:", n_bootstrap, "\n")
cat("Lambda used:", lambda_fixed, "\n")
cat("Peak infection rate:", round(max(f_hat_opt), 2), "\n")
cat("Peak infection time:", infection_dates[which.max(f_hat_opt)], "\n")
cat("Infection curve length:", length(f_hat_opt), "days\n")

# Save important results
bootstrap_results <- list(
  lambda = lambda_fixed,
  f_hat = f_hat_opt,
  f_bootstrap = bootstrap_f,
  f_ci = list(lower = f_lower, upper = f_upper),
  infection_dates = infection_dates,
  death_dates = t,
  beta_opt = beta_opt,
  gamma_opt = gamma_opt
)

cat("\nBootstrap analysis completed! Results saved in 'bootstrap_results'.\n")


