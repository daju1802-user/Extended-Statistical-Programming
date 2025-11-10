library(splines)
data <- read.table("engcov.txt", header = TRUE)
t <- data$julian  # Day of year
n <- length(y)  # Number of observations
y <- data$nhs
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


penalized_grad <- function(gamma, y, X, S, lambda) {
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
  # Analytic gradient
  grad_analytic <- pnll_grad(gamma, y, X, S, lambda)
  
  # Finite difference gradient
  eps <- 1e-7
  grad_fd <- numeric(length(gamma))
  for (i in 1:length(gamma)) {
    gamma_plus <- gamma
    gamma_plus[i] <- gamma[i] + eps
    gamma_minus <- gamma
    gamma_minus[i] <- gamma[i] - eps
    grad_fd[i] <- (pnll(gamma_plus, y, X, S, lambda) - 
                     pnll(gamma_minus, y, X, S, lambda)) / (2 * eps)
  }
  
  # Check agreement
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
fit <- optim(
  par = gamma0,
  fn = penalized_nll,
  gr = penalized_grad,
  y = y, X = mats$X, S = mats$S, lambda = lambda,
  method = "BFGS",
  control = list(maxit = 2000, trace = 1)
)

## ---- Extract fitted parameters ----
gamma_hat <- fit$par
beta_hat <- exp(gamma_hat)
mu_hat <- as.vector(mats$X %*% beta_hat)  # fitted daily deaths
f_hat <- as.vector(mats$tilde_X %*% beta_hat)  # estimated infection curve

## ---- Plot actual vs fitted daily deaths ----
plot(t_death, y, type = "p", pch = 16, col = "grey40",
     main = "Actual vs Fitted Daily Deaths",
     xlab = "Day (Julian date)", ylab = "Daily deaths")
lines(t_death, mu_hat, col = "red", lwd = 2)
legend("topright", legend = c("Observed deaths", "Fitted deaths"),
       col = c("grey40", "red"), lwd = c(NA, 2), pch = c(16, NA), bty = "n")

## ---- Plot estimated daily infection curve ----
plot(mats$t_f, f_hat, type = "l", col = "blue", lwd = 2,
     main = "Estimated Daily Infection Curve f(t)",
     xlab = "Day (Julian date)", ylab = "Estimated infections")

## add vertical line to show where death data starts
abline(v = min(t_death), lty = 2, col = "grey60")
text(min(t_death), max(f_hat)*0.9, "Deaths start", pos = 4, cex = 0.8)



## ---- Q4: Grid search for optimal lambda  ----

n <- length(y)
log_lambda_seq <- seq(-13, -7, length.out = 50)
lambda_seq <- exp(log_lambda_seq)
gamma_init <- rep(0, ncol(mats$X))

## Define helper function for one λ
# Function to calculate BIC
calc_bic <- function(gamma, lambda, X, y, S) {
  beta <- exp(gamma)
  mu <- as.vector(X %*% beta)
  
  # Log-likelihood (without penalty)
  ll <- sum(y * log(mu) - mu)
  
  # Calculate Hessian matrices for EDF
  # H_lambda = X^T * W * X + lambda * S, where W = diag(y/mu^2)
  W <- diag(as.vector(y / mu^2))
  H_0 <- crossprod(X, X * w)
  H_lambda <- H_0 + lambda * S
  
  # Effective degrees of freedom
  EDF <- sum(diag(solve(H_lambda) %*% H_0))
  
  cholH <- try(chol(H_lambda), silent = TRUE)
  if (inherits(cholH, "try-error")) {
    # Rollback plan
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






