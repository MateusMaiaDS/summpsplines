# Creating the D (difference matrix)
D_gen <- function(p, n_dif){
  return(diff(diag(p),diff = n_dif))
}

# Wrapping the model with a R code
rsp_sampler <- function(x_train,
                   y,
                   nIknots,
                   df = 3,
                   sigquant = 0.9,
                   n_dif = 2,
                   delta,
                   nu,
                   a_delta,
                   d_delta,
                   n_mcmc,
                   n_burn,
                   scale_y = TRUE){

  # Scaling x
  x_min <- apply(as.matrix(x_train),2,min)
  x_max <- apply(as.matrix(x_train),2,max)

  # Storing the original
  x_train_original <- x_train
  x_train_scale <- x_train

  # Normalising all the columns
  for(i in 1:ncol(x_train)){
    x_train_scale[,i] <- normalize_covariates_bart(y = x_train_scale[,i],a = x_min[i], b = x_max[i])
  }

  # Scaling the y
  min_y <- min(y)
  max_y <- max(y)


  # Need to create an array for B_train
  B_train <- array(data = NA,
                   dim = c(nrow(x_train_scale),
                           nIknots+3,
                           ncol(x_train_scale)))

  absolut_min <- apply(x_train_scale,2,min)
  absolut_max <-apply(x_train_scale,2,max)


  # Iterating over all covariates
  for(i in 1:ncol(x_train_scale)){
    # Getting the internal knots
    knots <- quantile(x_train_scale[,i],seq(0,1,length.out = nIknots+2))[-c(1,nIknots+2)]

    # Creating the B spline
    B_train[,,i] <- as.matrix(splines::bs(x = x_train_scale[,i],knots = knots,
                                     intercept = FALSE,
                                     Boundary.knots = c(absolut_min[i],absolut_max[i])))
  }

  # Getting a penalized version of B_train
  if(n_dif!=0){
    D <- D_gen(p = ncol(B_train[,,1]),n_dif = n_dif)
  } else {
    D <- diag(nrow = ncol(B_train[,,1]))
  }
  # IN CASE WE WANT TO USE THE DIFFERENCE PENALISATION DIRECTLY OVER THE
  #BASIS FUNCTION
  # B_train <- B_train%*%crossprod(D,solve(tcrossprod(D)))

  # By adding P in the the beta prior
  P <- crossprod(D)

  # Scaling "y"
  if(scale_y){
    y_scale <- normalize_bart(y = y,a = min_y,b = max_y)
  } else {
    y_scale <- y
  }

  # Calculating \tau_{mu}
  if(scale_y){
    tau_b <- tau_b_0 <-  (4*1*(2^2))
    # tau_b <- (4*1*(2^2))*10

  } else {
    tau_b <- tau_b_0 <- (4*1*(2^2))/((min_y-max_y)^2)
    # tau_b <- (4*1*(2^2))*10
  }
  # Getting the naive sigma value
  nsigma <- naive_sigma(x = x_train_scale,y = y_scale)

  # Calculating tau hyperparam
  a_tau <- df/2

  # Calculating lambda
  qchi <- stats::qchisq(p = 1-sigquant,df = df,lower.tail = 1,ncp = 0)
  lambda <- (nsigma*nsigma*qchi)/df
  d_tau <- (lambda*df)/2

  # Call the bart function
  tau_init <- nsigma^(-2)
  tau_init <- 50

  sampler_list <- sp_sampler(B_train = B_train,
                             D_m = D,
                             y = as.matrix(y_scale),tau_b = tau_b,tau_b_intercept = tau_b_0,
                             tau = tau_init,a_tau = a_tau,d_tau = d_tau,
                             nu = nu,delta = delta,a_delta = a_delta,
                             d_delta = d_delta,n_mcmc = n_mcmc,n_burn = n_burn)

  # Tidying up the posterior elements
  if(scale_y){
    y_train_post <- unnormalize_bart(z = sampler_list[[3]],a = min_y,b = max_y)
    tau_b_post <- sampler_list[[4]]/((max_y-min_y)^2)
    tau_post <-  sampler_list[[6]]/((max_y-min_y)^2)
  } else {
    y_train_post <- sampler_list[[3]]
    tau_b_post <- sampler_list[[4]]
    tau_post <-  sampler_list[[6]]

  }

  beta_post <- sampler_list[[1]]
  beta_0_post <- sampler_list[[2]]
  delta_post <- sampler_list[[5]]

  return(list(beta_post = beta_post,
              beta_0_post = beta_0_post,
              y_train_post = y_train_post,
              tau_post = tau_post,
              tau_b_post = tau_b_post,
              delta_post = delta_post))
}





