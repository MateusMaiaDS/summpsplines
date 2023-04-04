rm(list=ls())
library(tidyverse)

Rcpp::sourceCpp("src/sampler.cpp")
source("R/other_functions.R")
source("R/sampler.R")
n_ <- 500
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

# Transforming into data.frame
x <- as.data.frame(x)
x_test <- as.data.frame(x_new)


# Testing the mpsBART
bart_test <- rsp_sampler(x_train = x,y = unlist(c(y)),n_splines = 10,
                         n_mcmc = 2500,n_dif = 2,
                         nIknots = 30,delta = 1,df = 3,
                         sigquant = 0.9,nu = 2,
                         a_delta = 0.0001,d_delta = 0.0001,
                         n_burn = 1000,scale_y = TRUE)


# Running BART
bartmod <- dbarts::bart(x.train = x,y.train = unlist(c(y)),ntree = 200,x.test = x_test,keeptrees = TRUE)

# Convergence plots
par(mfrow = c(1,2))
plot(bart_test$tau_post,type = "l", main = expression(tau),ylab=  "")
plot(bartmod$sigma^-2, type = "l", main = paste0("BART: ",expression(tau)),ylab=  "")

par(mfrow = c(1,2))
plot(bart_test$y_train_post %>% colMeans(),y, main = 'mpsBART', xlab = "mpsBART pred", ylab = "y")
plot(bartmod$yhat.train.mean,y, main = "BART", xlab = "BART pred", ylab = "y")

# Comparing on the test set
pred_bart <- colMeans(predict(bartmod,fried_sim_new_sample$x))

# Storing the results
rmse(x = fried_sim$y,y = bartmod$yhat.train.mean)
rmse(x = fried_sim$y,y = colMeans(bart_test$y_train_post))


par(mfrow = c(1,1))
plot(bartmod$yhat.train.mean, colMeans(bart_test$y_train_post))

# Plotting all taus
par(mfrow=c(2,2))
plot(bart_test$tau_b_post[,1],type = "l", main = expression(tau[1]))
plot(bart_test$tau_b_post[,2],type = "l", main = expression(tau[2]))
plot(bart_test$tau_b_post[,3],type = "l", main = expression(tau[3]))
plot(bart_test$tau_b_post[,4],type = "l", main = expression(tau[4]))

