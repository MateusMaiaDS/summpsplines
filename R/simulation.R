rm(list=ls())
library(tidyverse)
source("R/other_functions.R")
source("R/sampler.R")
Rcpp::sourceCpp("src/sampler.cpp")

# Generating data
n_ <- 100
set.seed(42)

# Simulation 1
fried_sim <- mlbench::mlbench.friedman1(n = n_,sd = 0.01)
friedman_no_interaction <- function (n, sd = 1)
{
  x <- matrix(runif(4 * n), ncol = 4)
  y <- 10 * sin(pi * x[, 1] )
  y <- y + 20 * (x[, 2] - 0.5)^2 + 10 * x[, 3] + 5 * x[, 4]
  if (sd > 0) {
    y <- y + rnorm(n, sd = sd)
  }
  list(x = x, y = y)
}

sd_ <- 0.1
fried_sim <- friedman_no_interaction(n = n_,sd = sd_)
fried_sim_new_sample <- friedman_no_interaction(n = n_,sd = sd_)


x <- fried_sim$x[,,drop = FALSE]
x_new <- fried_sim_new_sample$x
y <- fried_sim$y


sp_mod <- rsp_sampler(x_train = x,y = y,
                      nIknots = 100,df = 3,
                      sigquant = 0.9,delta = 1,nu = 2,
                      a_delta = 0.0001,d_delta = 0.0001,
                      n_mcmc = 2500,n_burn = 500,
                      scale_y = TRUE)

# Formatting the sampler plot
par(mfrow=c(1,1))
plot(x,y,main = "P-Splines robust priors")
quantiles_y_hat <- apply(sp_mod$y_train_post,2,function(x){quantile(x,probs = c(0.025,0.5,0.975))})
lines(x,sin(3*x), col = "red")
lines(x,quantiles_y_hat[2,],col = "blue")
lines(x,quantiles_y_hat[1,],lty = "dashed", col = "blue")
lines(x,quantiles_y_hat[3,],lty = "dashed", col = "blue")

# Traceplots
par(mfrow=c(2,2))
plot(sp_mod$beta_0_post,type = "l", main = expression(beta[0]))
plot(sp_mod$tau_b_post,type = "l", main = expression(tau[b]))
plot(sp_mod$delta_post,type = "l", main = expression(delta))
plot(sp_mod$tau_post,type = "l", main = expression(tau))

