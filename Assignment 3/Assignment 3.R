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
  lower <- min(t) - 30
  upper <- max(t)
  mid_squence <- seq(lower, upper, length.out = k-2)
  gap <- mid_squence[2] - mid_squence[1]
  knot <- c(
    seq(from = lower - 3*gap, to = lower - gap, length.out = 3),
    mid_squence,
    seq(from = upper + gap, to = upper + 3*gap, length.out = 3)
  )

  X_tilde <- splineDesign(knot,lower:upper,outer.ok = TRUE) #outer.ok = TRUE
  
  r <- nrow(X_tilde)

  X <- matrix(0, n, ncol(X_tilde))
  # Build X by convolution
  for (i in 1:n) {
    max_j <- min(29 + i, length(pd))
    idx <- t[i] - (1:max_j) - lower + 1
    valid <- which(idx >= 1 & idx <= r)
    if (length(valid) > 0) {
      X[i, ] <- colSums(X_tilde[idx[valid], , drop = FALSE] * pd[valid])
    }
  }

  S <- crossprod(diff(diag(k),diff=2))
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
  beta <- exp(gamma)
  mu <- as.vector(X %*% beta)
  mu[mu <= 1e-12] <- 1e-12
  
  grad_ll <- t(X) %*% (1 - y / mu)  # Gradient log-likelihood
  grad_pen <- lambda * (S %*% beta)  # Gradient penalty
  
  grad <- as.vector(beta * (grad_ll + grad_pen))  # Chain rule
  return(grad)
}

# ---- Finite differencing test ----
# Choose initial gamma
mats <- create_design_matrices(t_death, K = 80)
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



