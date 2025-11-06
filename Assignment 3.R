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
## 3. Constructs the design matrix X for μ_i:
##    Each row corresponds to one observed death day t_i.
##    Each column corresponds to a B-spline coefficient β_k.
##    Each μ_i is computed as a weighted sum of the basis functions
##    evaluated at the infection times t_i - j, weighted by π(j).
## 4. Computes the smoothing penalty matrix S:
##    S = t(D) %*% D, where D is the second-difference matrix.
##    Penalizes roughness in β, encouraging smooth f(t).