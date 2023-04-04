#include "mpsplines_sampler.h"
#include <random>
#include <Rcpp.h>
using namespace std;

// Building the constructor
modelParam::modelParam(arma::cube B_train_,
                       arma::mat y_,
                       double tau_b_,
                       double tau_b_intercept_,
                       double tau_,
                       double a_tau_,
                       double d_tau_,
                       double nu_,
                       double delta_,
                       double a_delta_,
                       double d_delta_,
                       int n_mcmc_,
                       int n_burn_){

  // Create a matrix of ones
  arma::vec ones(y_.n_rows,arma::fill::ones);
  bt_ones = arma::mat(B_train_.n_cols,B_train_.n_slices);

  for(int i = 0;i<B_train_.n_slices;i++){
      bt_ones.col(i) = B_train_.slice(i).t()*ones;
  }

  // Filling the parameters
  y = y_;
  B_train = B_train_;
  d = B_train_.n_slices;
  tau_b = arma::vec(B_train_.n_slices,arma::fill::ones)*tau_b_;
  tau = tau_;
  a_tau = a_tau_;
  d_tau = d_tau_;
  nu = nu_;
  delta = arma::vec(B_train_.n_slices,arma::fill::ones)*delta_;
  a_delta = a_delta_;
  d_delta = d_delta_;
  n_mcmc = n_mcmc_;
  n_burn = n_burn_;
  p = B_train_.n_cols;


}

// // Create a function to generate matrix D (for penalisation)
arma::mat D(modelParam data){

        // Creating the matrix elements
        arma::mat D_m((data.p-2),data.p,arma::fill::zeros);

        for(int i=0;i<(data.p-2);i++){
                D_m(i,i) = 1;
                D_m(i,i+1) = -2;
                D_m(i,i+2) = 1;
        }

        return D_m;
}

// // Create a function to generate matrix D (for penalisation)
arma::mat D_diag(modelParam data){

  // Creating the matrix elements
  arma::mat D_m((data.p),data.p,arma::fill::zeros);

  for(int i=0;i<(data.p);i++){
    D_m(i,i) = 1;

  }

  return D_m;
}

// // Create a function to generate matrix D (for penalisation)
arma::mat D_first(modelParam data){

  // Creating the matrix elements
  arma::mat D_m((data.p-1),data.p,arma::fill::zeros);

  for(int i=0;i<(data.p-1);i++){
    D_m(i,i) = -1;
    D_m(i,i+1) = 1;
  }

  return D_m;
}

// Building the beta sampler
void beta_sampler(arma::mat& betas,
                  double& beta_0,
                  modelParam& data){


  // Iterating over all covariates
  for(int j = 0;j<data.d;j++){


      // Cov aux mat
      arma::mat cov_sum_aux(data.y.size(), 1,arma::fill::zeros);


      // Getting the sum element
      for(int k = 0; k < data.d; k++){
        if(k == j){
          continue;
        }
        cov_sum_aux = cov_sum_aux + data.B_train.slice(k)*betas.col(k);
      }

      // Getting the precision aux factor
      arma::mat inv_precision_aux = arma::inv(data.B_train.slice(j).t()*data.B_train.slice(j)+(data.tau_b(j)/data.tau)*data.P);
      arma::mat mean_aux = inv_precision_aux*(data.B_train.slice(j).t()*data.y-data.B_train.slice(j).t()*(beta_0+cov_sum_aux));
      arma::mat cov_aux = (1/data.tau)*inv_precision_aux;

      // cout << "Error sample BETA" << endl;
      arma::mat sample = arma::randn<arma::mat>(betas.n_rows);
      // cout << "Error variance" << endl;
      betas.col(j) = arma::chol(cov_aux,"lower")*sample + mean_aux;
  }

  return;
}

// Building the beta_0 sampler
void beta_0_sampler(arma::mat& betas,
                    double& beta_0,
                    modelParam& data){

  double mean_sum_aux = 0.0;

  // Getting the sum element
  for(int k = 0; k < data.d; k++){
      mean_sum_aux = mean_sum_aux + arma::as_scalar(betas.col(k).t()*(data.bt_ones.col(k)));
  }


  double s_gamma = data.y.size()+(data.tau_b_intercept/data.tau);
  double mean_aux = (1/s_gamma)*(arma::accu(data.y)-mean_sum_aux);
  double sd_aux = sqrt(1/(s_gamma*data.tau));

  beta_0 = arma::randn()*sd_aux + mean_aux;

}


// Building the \tau_b sampler
void tau_b_sampler(arma::mat& betas,
                   modelParam& data){

  for(int j; j< data.d; j++){
      // Calculating the shape and rate parameter
      double tau_b_shape = 0.5*data.p+0.5*data.nu;
      double rate_aux = arma::as_scalar(betas.col(j).t()*data.P*betas.col(j));
      double tau_b_rate = 0.5*rate_aux + 0.5*data.delta(j)*data.nu;

      data.tau_b(j) = R::rgamma(tau_b_shape,1/tau_b_rate);
  }

  return;
}

// Updating delta
void delta_sampler(modelParam& data){


  for(int j = 0; j<data.d; j++){
    // Calculating shape and rate parameter
    double delta_shape = 0.5*data.nu+data.a_delta;
    double delta_rate = 0.5*data.nu*data.tau_b(j) + data.d_delta;

    data.delta(j) = R::rgamma(delta_shape,1/delta_rate);
  }

  return;
}


// Updating the tau parameter
void tau_sampler(modelParam& data,
                 arma::vec& y_hat){

  double tau_res_sq_sum = dot((y_hat-data.y),(y_hat-data.y));

  data.tau = R::rgamma((0.5*data.y.n_rows+data.a_tau),1/(0.5*tau_res_sq_sum+data.d_tau));

  return;
}

// Generating the sample code for the sampler
//[[Rcpp::export]]
Rcpp::List sp_sampler(arma::cube B_train,
                     arma::mat y,
                     arma::mat D_m,
                     double tau_b,
                     double tau_b_intercept,
                     double tau,
                     double a_tau,
                     double d_tau,
                     double nu,
                     double delta,
                     double a_delta,
                     double d_delta,
                     int n_mcmc,
                     int n_burn){


    cout << "Error 0" << endl;

    // Initalising the data object
    modelParam data(    B_train,
                        y,
                        tau_b,
                        tau_b_intercept,
                        tau,
                        a_tau,
                        d_tau,
                        nu,
                        delta,
                        a_delta,
                        d_delta,
                        n_mcmc,
                        n_burn);

    // Generating the P matrix
    data.P = D_m.t()*D_m + arma::eye(B_train.n_cols,B_train.n_cols)*1e-8;

    // Initializing the vector of betas
    cout << "Error 1" << endl;
    arma::mat betas(data.p,data.d, arma::fill::ones);
    double beta_0 = 0;
    arma::vec y_hat(data.y.n_rows,arma::fill::zeros);

    // Storing the posteriors
    int n_post = n_mcmc-n_burn;
    arma::mat beta_post(n_post,data.p,arma::fill::ones);
    arma::vec beta_0_post(n_post,arma::fill::ones);
    arma::mat y_hat_post(n_post,data.y.n_rows,arma::fill::ones);

    arma::vec tau_post(n_post,arma::fill::ones);
    arma::mat tau_b_post(n_post,data.d,arma::fill::ones);
    arma::mat delta_post(n_post,data.d,arma::fill::ones);
    int post_iter = 0;


    // Initializing the sampling processes
    for(int i = 0; i < data.n_mcmc; i++){

      // cout << "Beta error" << endl;
      beta_sampler(betas,beta_0, data);
      // cout << "Beta_0 error" << endl;
      beta_0_sampler(betas, beta_0, data);
      // cout << "Tau_b error" << endl;
      tau_b_sampler(betas,data);
      cout << "Tau_b value corresponds to: ";
      for(int t = 0; t<data.d;t++){
        cout << data.tau_b(t) << " ";
      }
      cout << endl;
      // cout << "Delta error" << endl;
      delta_sampler(data);

      // Calculating the predictions
      // cout << "Y_hat error" << endl;
      arma::vec y_hat_aux(data.y.size(),arma::fill::zeros);
      for(int j = 0;j<data.d;j++){
        y_hat_aux = y_hat_aux + data.B_train.slice(j)*betas.col(j);
      }

      y_hat = y_hat_aux + beta_0;

      // cout << "Tau sampler error" << endl;
      tau_sampler(data,y_hat);


      // cout << " No error" << endl;
      // Iterating and storing the posterior samples
      if(i >= data.n_burn){
        // beta_post.row(post_iter) = betas.t();
        beta_0_post(post_iter) = beta_0;
        y_hat_post.row(post_iter) = y_hat.t();
        tau_b_post.row(post_iter) = data.tau_b.t();
        delta_post.row(post_iter) = data.delta.t();
        tau_post(post_iter) = data.tau;
        post_iter++;
      }


    }

    return Rcpp::List::create(beta_post,
                              beta_0_post,
                              y_hat_post,
                              tau_b_post,
                              delta_post,
                              tau_post);
}




