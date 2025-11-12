#Group 6 :Xuyi Shi(S2796586), Chi Zhang(S2828433), Jiachen Guang(S2789777)
#Repo Link: https://github.com/daju1802-user/Extended-Statistical-Programming.git

## This code implements a smooth deconvolution model that reconstructs the daily infection trend 
## from reported COVID-19 deaths in hospitals.
## The model accounts for the infection-to-death delay distribution and 
## reconstructs the latent infection curve f(t) using a B-spline representation.
## By comparing model fits across a range of smoothing parameters and quantifying uncertainty through 
## non-parametric bootstrap resampling, we obtain robust estimates of the underlying infection dynamics along with confidence intervals.
## The analysis proceeds through four key computational stages: model specification and design matrix construction, 
## penalized likelihood optimization with gradient-based methods, smoothing parameter selection via information criteria, 
## and uncertainty quantification through bootstrap replication.
## These results give insight into how infection levels changed over time and demonstrate a practical, 
## data-driven approach for reconstructing infection dynamics from mortality data.

#Contribution:
#Xuyi Shi(%): 

#Chi Zhang(%): 

#Jiachen Guang(%): 


library(splines)
data <- read.table("engcov.txt", header = TRUE)
t <- data$julian  # Day of year
y <- data$nhs
n <- length(y)  # Number of observations


#============================================
# Step 1: Create model matrices X_tilde, X, and penalty matrix S
#============================================
d <- 1:80
edur <- 3.151  # Mean of log infection-to-death duration
sdur <- 0.469  # SD of log infection-to-death duration
pd <- dlnorm(d, edur, sdur)  # Log-normal distribution for infection-to-death
pd <- pd / sum(pd)  # Normalize to probability distribution
K <- 80  # Number of basis functions

build_matrices <- function(t,k,pd){
  ## Construct design matrices for the Poisson deconvolution model.
  ## This function prepares the matrices required to estimate the infection curve f(t),
  ## which is related to observed deaths through a convolution with the infection-to-death distribution.
  lower <- min(t) - 30  # Start infection curve 30 days before first death
  upper <- max(t)       # End infection curve at last death date
  
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
    max_j <- min(29 + i, length(pd))  # Limit lookback to available probability data
    # compute indices of infection times corresponding to death day
    idx <- t[i] - (1:max_j) - lower + 1
    valid <- which(idx >= 1 & idx <= r)
    
    # if valid infections exist, add weighted basis function contributions
    if (length(valid) > 0) {
      X[i, ] <- colSums(X_tilde[idx[valid], , drop = FALSE] * pd[valid])
    }
  }
  
  S <- crossprod(diff(diag(k),diff=2)) # construct penalty matrix (2nd differences)
  return(list(X_tilde=X_tilde,S=S,X=X,t_infection = lower:upper))
}

mats <- build_matrices(t,k = 80, pd)
X <- mats$X          # Design matrix for death model
X_tilde <- mats$X_tilde  # B-spline basis for infection curve
S <- mats$S          # Penalty matrix
head(X)

#============================================
# Step 2: Likelihood and gradient functions
#============================================
pnll <- function(gamma, y, X, S, lambda, w = rep(1, length(y))) {
  ## Computes the penalized negative log-likelihood and its gradient 
  ## for the Poisson deconvolution model used to estimate the daily infection curve f(t).
  ## The model assumes that observed daily deaths y_i follow a Poisson distribution
  ## A smoothness penalty on beta is included to control overfitting.
  ## The functions are designed for use with 'optim' for numerical optimization.
  beta <- exp(gamma)  # ensure positivity of beta via exponential transformation
  mu <- as.vector(X %*% beta)  # Compute expected mean number of deaths under the model
  
  # Poisson log-likelihood 
  ll <- sum(w * (y * log(mu) - mu - lgamma(y+1))) # w: bootstrap weights (default all 1s)
  
  # Add smoothness penalty term (penalizes roughness in beta coefficients)
  penalty <- lambda * sum(beta * (S %*% beta)) / 2
  # Return negative penalized log-likelihood
  -(ll - penalty)
}


pnll_grad <- function(gamma, y, X, S, lambda, w = rep(1, length(y))) {
  ## Compute the gradient of the penalized negative log-likelihood
  ## for the Poisson deconvolution model estimating daily infections.
  beta <- exp(gamma) # Transform to positive scale
  mu <- as.vector(X %*% beta) # Expected deaths
  
  # Gradient of log-likelihood w.r.t. gamma
  # dl/dgamma = F, where F = diag(y/mu - 1) * X * diag(beta)
  residual <- w * (y / mu - 1)  # Poisson residuals with weight
  F <- (residual * X) * rep(beta, each = length(y))  # Element-wise multiplication
  dll_dgamma <- colSums(F)
  
  # Gradient of penalty w.r.t. gamma
  # dP/dgamma = diag(beta) * S * beta
  dP_dgamma <- beta * (S %*% beta)
  
  # Return negative gradient (for minimization)
  -(dll_dgamma - lambda * dP_dgamma)
}

# ---- Finite differencing test ----
test_gradient <- function(gamma, y, X, S, lambda) {
  ## Test the correctness of the analytic gradient function by finite differencing.
  # Analytic gradient from the implemented function
  grad_analytic <- pnll_grad(gamma, y, X, S, lambda)
  
  # Finite difference gradient
  eps <- 1e-7  # Small perturbation for finite differences
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
lambda <- 5e-5  # Initial smoothing parameter for testing
t_infection = mats$t_infection
# Optimize using BFGS method
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
  ##Calculate the Bayesian Information Criterion (BIC) for the Poisson deconvolution model,
  ##used for selecting the smoothing parameter lambda. BIC balances model fit and complexity, 
  ##helping to choose the optimal degree of smoothing.
  beta <- exp(gamma)
  mu <- as.vector(X %*% beta)
  
  # Log-likelihood (without penalty)
  ll <- sum(y * log(mu) - mu)
  
  # Calculate Hessian matrices for EDF
  # H_lambda = X^T * W * X + lambda * S, where W = diag(y/mu^2)
  W <- as.vector(y / mu^2)  # Weights for Poisson GLM
  H_0 = crossprod(X, X * W)  # Hessian without penalty
  H_lambda <- H_0 + lambda * S  # Hessian with penalty
  
  # Calculate Effective Degrees of Freedom (EDF)
  # EDF = trace(H_lambda^{-1} * H_0)
  cholH <- try(chol(H_lambda), silent = TRUE)  # Try Cholesky decomposition
  if (inherits(cholH, "try-error")) {
    # Fallback: traditional matrix inversion
    EDF <- sum(diag(solve(H_lambda, H_0)))
  } else {
    # Efficient method using Cholesky decomposition
    Z <- backsolve(cholH, forwardsolve(t(cholH), H_0))
    EDF <- sum(diag(Z))
  }
  # BIC = -2 * log-likelihood + log(n) * EDF
  BIC <- -2 * ll + log(n) * EDF
  list(BIC = BIC, EDF = EDF, ll = ll)
}


# Grid search over log(lambda)
log_lambda <- seq(-13, -7, length = 50)  # Search range for log(lambda)
lambda_seq <- exp(log_lambda)  # Convert to actual lambda values
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
  
  gamma_current <- fit_i$par  # Use as starting point for next lambda
  
  # Calculate BIC
  bic_info <- calc_bic(gamma_current, lambda_i, X, y, S)
  bic_values[i] <- bic_info$BIC
  edf_values[i] <- bic_info$EDF
  
}
# Find optimal lambda (minimizing BIC)
idx_opt <- which.min(bic_values)
lambda_opt <- lambda_seq[idx_opt]
cat("\nOptimal lambda:", lambda_opt, "\n")
cat("Optimal BIC:", bic_values[idx_opt], "\n")
cat("Optimal EDF:", edf_values[idx_opt], "\n")

# Refit with optimal lambda to get final parameters
fit_opt <- optim(gamma_init, pnll, pnll_grad,
                 lambda = lambda_opt, X = X, y = y, S = S,
                 method = "BFGS", control = list(maxit = 1000))
# Extract optimal parameters
gamma_opt <- fit_opt$par
beta_opt <- exp(gamma_opt)


# ==================== Q5: Non-parametric Bootstrap ====================

# Use optimal lambda for bootstrap (fixed at optimal value from Q4)
lambda_fixed <- lambda_opt  
n_bootstrap <- 200  # Number of bootstrap replicates
n_obs <- length(y)  # Number of observations

# Store bootstrap results
bootstrap_f <- matrix(0, nrow = nrow(X_tilde), ncol = n_bootstrap)  # Infection curves
bootstrap_beta <- matrix(0, nrow = K, ncol = n_bootstrap)  # Beta parameters

# Use optimal solution as starting value for better convergence
gamma_start <- gamma_opt

# Bootstrap loop: generate 200 bootstrap samples
for (b in 1:n_bootstrap) {
  # Generate bootstrap weights (frequency counts for resampled observations)
  wb <- tabulate(sample(n_obs, replace = TRUE), n_obs)
  
  # Optimize using weighted functions for this bootstrap sample
  fit_bootstrap <- optim(
    par = gamma_start,
    fn = pnll,
    gr = pnll_grad,
    y = y, X = X, S = S, lambda = lambda_fixed,
    w = wb,
    method = "BFGS",
    control = list(maxit = 1000)
  )
  
  # Store results for this bootstrap sample
  gamma_b <- fit_bootstrap$par
  beta_b <- exp(gamma_b)
  f_b <- as.vector(X_tilde %*% beta_b)  # Infection curve for this bootstrap
  
  bootstrap_beta[, b] <- beta_b
  bootstrap_f[, b] <- f_b
  
}


# Calculate optimal fitted infection curve from original data
f_hat_opt <- as.vector(X_tilde %*% beta_opt)  

# Calculate bootstrap confidence intervals (2.5% and 97.5% quantiles)
f_bootstrap_ci <- apply(bootstrap_f, 1, function(x) {
  quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
})

f_lower <- f_bootstrap_ci[1, ]  # 95% CI lower bound (2.5% quantile)
f_upper <- f_bootstrap_ci[2, ]  # 95% CI upper bound (97.5% quantile)

# ==================== Q6: Combined Plot ====================

# Set graphics parameters for single plot
par(mar = c(5, 4, 4, 4))

# Create empty plot with appropriate ranges
plot_range_x <- range(c(t, infection_dates))
plot_range_y <- range(c(0, max(y) * 1.1, max(f_upper) * 1.1))

# Create main plot frame
plot(NA, type = "n", 
     xlim = plot_range_x, ylim = c(0, max(y) * 1.1),
     main = "Daily Deaths and Estimated Infections",
     xlab = "Day of year 2020", 
     ylab = "Daily Deaths")

# Plot death data and model fit (left y-axis)
points(t, y, pch = 16, col = "black", cex = 0.7)
lines(t, mu_opt, col = "red", lwd = 2)

# Add second y-axis for infections
par(new = TRUE)
plot(NA, type = "n", 
     xlim = plot_range_x, 
     ylim = c(0, max(f_upper) * 1.1),
     axes = FALSE, xlab = "", ylab = "")
axis(4, col = "blue", col.axis = "blue")
mtext("Estimated Infections", side = 4, line = 3, col = "blue")

# Plot infection rate with confidence intervals (right y-axis)
polygon(c(infection_dates, rev(infection_dates)),
        c(f_lower, rev(f_upper)),
        col = rgb(0, 0, 1, 0.2), border = NA)
lines(infection_dates, f_hat_opt, col = "blue", lwd = 2)

# Add vertical line marking start of death data
death_start <- min(t)
abline(v = death_start, lty = 2, col = "red")
text(death_start, max(f_upper) * 0.9, 
     paste("Death Data Starts\n(Infection +", round(exp(edur)), "days)", sep = ""),
     pos = 2, cex = 0.8, col = "red")

# Add legend
legend("topright",
       legend = c("Observed deaths", "Fitted deaths", 
                  "Inferred infections", "95% Confidence Interval"),
       col = c("black", "red", "blue", rgb(0, 0, 1, 0.3)),
       pch = c(16, NA, NA, NA),
       lwd = c(NA, 2, 2, 8), lty = c(NA, 1, 1, 1),
       bty = "n")

# Reset graphics parameters to default
par(mar = c(5, 4, 4, 2) + 0.1)

cat("Peak infection rate:", round(max(f_hat_opt), 2), "\n")
cat("Peak infection time:", infection_dates[which.max(f_hat_opt)], "\n")







