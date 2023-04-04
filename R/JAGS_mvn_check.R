# Header ------------------------------------------------------------------

# P-spline model in JAGS with robust specification of the roughness of the penalty.

# Mateus Maia & Andrew Parnell

# This file fits a spline regression model to data in JAGS, adding a roughness
# parameter into the precision of the splines to yielding a more robust
# hyperparameter specification

# Some boiler plate code to clear the workspace and load in required packages
rm(list = ls())
library(R2jags)
library(MASS) # Useful for mvrnorm function
library(splines) # Useful for creating the B-spline basis functions
library(boot) # For real example
source("R/other_functions.R")
source("R/sampler.R")
# Maths -------------------------------------------------------------------

# Notation:
# y: vector of all observations
# B: design matrix of spline basis functions
# beta; spline weights
# tau; residual precision
# tau_beta; spline random walk precision parameter
# delta; the roughness parameter for tau_b prior

# Likelihood
# y ~ N(B%*%beta, tau^-1)
# beta_j ~ N (beta_{j-1},tau_b^-1)
# beta_1 ~ N (0,tau_b_0^-1)

# Priors
# tau ~ gamma(a_tau, d_tau)
# tau_b ~ gamma(0.5*nu, 0.5*delta*nu)
# delta ~ gamma(a_delta,d_delta)

# Simulate data -----------------------------------------------------------


# Generating data
n_ <- 500
set.seed(42)
# Simulation 1
x <- matrix(seq(-pi,pi,length.out = n_))
x_new <- matrix(seq(-pi,pi,length.out = n_))
colnames(x) <- "x"
colnames(x_new) <- "x"
y <- sin(3*x) + rnorm(n = n_,sd = 0.1)
y <- y[x>0,,drop = FALSE]
x <- x[x>0,,drop = FALSE]

# Storing the original
x_train_original <- as.matrix(x)
x_train_scale <- as.matrix(x)

# Scaling x
x_min <- apply(as.matrix(x_train_original),2,min)
x_max <- apply(as.matrix(x_train_original),2,max)

# Normalising all the columns
# for(i in 1:ncol(x_train_scale)){
#   x_train_scale[,i] <- normalize_covariates_bart(y = x_train_scale[,i],a = x_min[i], b = x_max[i])
# }

# Scaling the y
min_y <- min(y)
max_y <- max(y)

absolut_min <- min(min(x_train_scale[,1]))
absolut_max <- max(max(x_train_scale[,1]))

# Setting the nIknots and the basis that will be used
nIknots <- 100
# Getting the internal knots
knots <- quantile(x_train_scale[,1],seq(0,1,length.out = nIknots+2))[-c(1,nIknots+2)]

# Creating the B spline
B_train <- as.matrix(splines::ns(x = x_train_scale[,1],knots = knots,
                                 intercept = FALSE,
                                 Boundary.knots = c(absolut_min,absolut_max)))

# Deciding if scale or not
# y_scale <- normalize_bart(y = y,a = min_y,b = max_y)
y_scale <- y

# Getting the naive sigma value
nsigma <- naive_sigma(x = x_train_scale,y = y_scale)

# Calculating tau hyperparam
df <- 3
sigquant <- 0.9
a_tau <- df/2

# Calculating lambda
qchi <- stats::qchisq(p = 1-sigquant,df = df,lower.tail = 1,ncp = 0)
lambda <- (nsigma*nsigma*qchi)/df
d_tau <- (lambda*df)/2

# Calculating penalisation matrix
n_diffs <- 1
D <- diff(diag(ncol(B_train)), diff = n_diffs )
P <- crossprod(D)

# Jags code ---------------------------------------------------------------

# Jags code to fit the model to the simulated data
model_code <- "
model {
  # Likelihood
  for (t in 1:N) {
    y[t] ~ dnorm(inprod(B[t,1:N_knots],beta[1:N_knots])+beta_intercept, tau)
  }

  # RW prior on beta
  beta_intercept ~ dnorm(0,tau_b_0)
  beta[1:N_knots] ~ dmnorm(rep(0,N_knots),tau_b*P)

  # Priors on beta values
  tau ~ dgamma(a_tau, d_tau)
  tau_b ~ dgamma(0.5 * nu, 0.5 * delta * nu)
  delta ~ dgamma(a_delta, d_delta)
}
"

# Set up the data
model_data <- list(
  N = nrow(B_train),
  y = c(y_scale),
  B = B_train,
  N_knots = ncol(B_train),
  a_tau = a_tau,
  d_tau = d_tau,
  a_delta = 0.0001, # Default values used in Jullion, A. and Lambert, P., 2007.
  d_delta = 0.0001, # Default values used in Jullion, A. and Lambert, P., 2007.
  tau_b_0 = 16/(diff(range(y_scale))^2),
  nu = 2,
  P = P
) # Default values used in Jullion, A. and Lambert, P., 2007.

# Choose the parameters to watch
model_parameters <- c("beta", "tau", "tau_b", "delta","beta_intercept")

# Run the model - can be slow
model_run <- jags(
  data = model_data,
  parameters.to.save = model_parameters,
  model.file = textConnection(model_code)
)
plot(model_run)

# Simulated results -------------------------------------------------------

# Results and output of the simulated example, to include convergence checking, output plots, interpretation etc
# print(model_run)

# Get the posterior betas and 50% CI
beta_post <- model_run$BUGSoutput$sims.list$beta
beta_quantile <- apply(beta_post, 2, quantile, prob = c(0.25, 0.5, 0.75))

# New prediction
# Plot the output with uncertainty bands
par(mfrow = c(1,1))
plot(x, y_scale)
lines(x, B_train %*% beta_quantile[2, ] + c(model_run$BUGSoutput$mean$beta_intercept), col = "blue") # Predicted line
lines(x, B_train %*% beta_quantile[1, ] + c(model_run$BUGSoutput$mean$beta_intercept), col = "blue", lty = 2) # Predicted low
lines(x, B_train %*% beta_quantile[3, ] + c(model_run$BUGSoutput$mean$beta_intercept), col = "blue", lty = 2) # Predicted high
legend("topleft", c(
  "True line",
  "Posterior lines (with 50% CI)",
  "Data"
),
lty = c(1, 1, -1),
pch = c(-1, -1, 1),
col = c("red", "blue", "black")
)



# Traceplots
par(mfrow=c(2,2))
plot(model_run$BUGSoutput$sims.list$beta_intercept,type = "l", main = expression(beta[0]))
plot(model_run$BUGSoutput$sims.list$tau_b,type = "l", main = expression(tau[b]))
plot(model_run$BUGSoutput$sims.list$delta,type = "l", main = expression(delta))
plot(model_run$BUGSoutput$sims.list$tau,type = "l", main = expression(tau))

