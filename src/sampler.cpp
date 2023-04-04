#include "mpsplines_sampler.h"
#include <random>
#include <Rcpp.h>
using namespace std;

// Building the constructor
modelParam::modelParam(arma::cube B_train_,
                       arma::mat y_,
                       int n_splines_,
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
  n_splines = n_splines_;
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

//[[Rcpp::export]]
arma::mat sum_exclude_col(arma::mat mat, int exclude_int){

  // Setting the sum matrix
  arma::mat m(mat.n_rows,1);

  if(exclude_int==0){
    m = sum(mat.cols(1,mat.n_cols-1),1);
  } else if(exclude_int == (mat.n_cols-1)){
    m = sum(mat.cols(0,mat.n_cols-2),1);
  } else {
    m = arma::sum(mat.cols(0,exclude_int-1),1) + arma::sum(mat.cols(exclude_int+1,mat.n_cols-1),1);
  }

  return m;
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
                  modelParam& data,
                  arma::vec& partial_residuals){


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
      arma::mat mean_aux = inv_precision_aux*(data.B_train.slice(j).t()*partial_residuals-data.B_train.slice(j).t()*(beta_0+cov_sum_aux));
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
                    modelParam& data,
                    arma::vec partial_residuals){

  double mean_sum_aux = 0.0;

  // Getting the sum element
  for(int k = 0; k < data.d; k++){
      mean_sum_aux = mean_sum_aux + arma::as_scalar(betas.col(k).t()*(data.bt_ones.col(k)));
  }

  // cout << "error 1" << endl;

  double s_gamma = data.y.size()+(data.tau_b_intercept/data.tau);
  // cout << "error 2" << endl;
  double mean_aux = (1/s_gamma)*(arma::accu(partial_residuals)-mean_sum_aux);
  // cout << "error 3" << endl;
  double sd_aux = sqrt(1/(s_gamma*data.tau));

  beta_0 = arma::randn()*sd_aux + mean_aux;
  // cout << "error 4" << endl;

}


// Building the \tau_b sampler
void tau_b_sampler(arma::cube& betas,
                   modelParam& data){

    arma::vec tau_b_shape_counter(data.d,arma::fill::zeros);
    arma::vec tau_b_rate_counter(data.d,arma::fill::zeros);

    for(int j = 0; j< data.d; j++){

            for(int k = 0; k < data.n_splines;k++){
                // Calculating the shape and rate parameter
                tau_b_shape_counter(j) = tau_b_shape_counter(j) + data.p;
                tau_b_rate_counter(j) =  tau_b_shape_counter(j) + arma::as_scalar(betas.slice(k).col(j).t()*data.P*betas.slice(k).col(j));
            }
            data.tau_b(j) = R::rgamma(0.5*data.nu+0.5*tau_b_shape_counter(j),1/(0.5*data.nu*data.delta(j)+0.5*tau_b_rate_counter(j)));
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

// // Generating the sample code for the sampler
// //[[Rcpp::export]]
// Rcpp::List sp_sampler(arma::cube B_train,
//                      arma::mat y,
//                      arma::mat D_m,
//                      double tau_b,
//                      double tau_b_intercept,
//                      double tau,
//                      double a_tau,
//                      double d_tau,
//                      double nu,
//                      double delta,
//                      double a_delta,
//                      double d_delta,
//                      int n_mcmc,
//                      int n_burn){
//
//
//     cout << "Error 0" << endl;
//
//     // Initalising the data object
//     modelParam data(    B_train,
//                         y,
//                         tau_b,
//                         tau_b_intercept,
//                         tau,
//                         a_tau,
//                         d_tau,
//                         nu,
//                         delta,
//                         a_delta,
//                         d_delta,
//                         n_mcmc,
//                         n_burn);
//
//     // Generating the P matrix
//     data.P = D_m.t()*D_m + arma::eye(B_train.n_cols,B_train.n_cols)*1e-8;
//
//     // Initializing the vector of betas
//     cout << "Error 1" << endl;
//     arma::mat betas(data.p,data.d, arma::fill::ones);
//     double beta_0 = 0;
//     arma::vec y_hat(data.y.n_rows,arma::fill::zeros);
//
//     // Storing the posteriors
//     int n_post = n_mcmc-n_burn;
//     arma::mat beta_post(n_post,data.p,arma::fill::ones);
//     arma::vec beta_0_post(n_post,arma::fill::ones);
//     arma::mat y_hat_post(n_post,data.y.n_rows,arma::fill::ones);
//
//     arma::vec tau_post(n_post,arma::fill::ones);
//     arma::mat tau_b_post(n_post,data.d,arma::fill::ones);
//     arma::mat delta_post(n_post,data.d,arma::fill::ones);
//     int post_iter = 0;
//
//
//     // Initializing the sampling processes
//     for(int i = 0; i < data.n_mcmc; i++){
//
//       // cout << "Beta error" << endl;
//       beta_sampler(betas,beta_0, data, data.y);
//       // cout << "Beta_0 error" << endl;
//       beta_0_sampler(betas, beta_0, data, data.y);
//       // cout << "Tau_b error" << endl;
//       tau_b_sampler(betas,data);
//       cout << "Tau_b value corresponds to: ";
//       for(int t = 0; t<data.d;t++){
//         cout << data.tau_b(t) << " ";
//       }
//       cout << endl;
//       // cout << "Delta error" << endl;
//       delta_sampler(data);
//
//       // Calculating the predictions
//       // cout << "Y_hat error" << endl;
//       arma::vec y_hat_aux(data.y.size(),arma::fill::zeros);
//       for(int j = 0;j<data.d;j++){
//         y_hat_aux = y_hat_aux + data.B_train.slice(j)*betas.col(j);
//       }
//
//       y_hat = y_hat_aux + beta_0;
//
//       // cout << "Tau sampler error" << endl;
//       tau_sampler(data,y_hat);
//
//
//       // cout << " No error" << endl;
//       // Iterating and storing the posterior samples
//       if(i >= data.n_burn){
//         // beta_post.row(post_iter) = betas.t();
//         beta_0_post(post_iter) = beta_0;
//         y_hat_post.row(post_iter) = y_hat.t();
//         tau_b_post.row(post_iter) = data.tau_b.t();
//         delta_post.row(post_iter) = data.delta.t();
//         tau_post(post_iter) = data.tau;
//         post_iter++;
//       }
//
//
//     }
//
//     return Rcpp::List::create(beta_post,
//                               beta_0_post,
//                               y_hat_post,
//                               tau_b_post,
//                               delta_post,
//                               tau_post);
// }


// Generating the sample code for the sampler
//[[Rcpp::export]]
Rcpp::List sum_sp_sampler(arma::cube B_train,
                      arma::mat y,
                      int n_splines,
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


  // Initalising the data object
  modelParam data(    B_train,
                      y,
                      n_splines,
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
  arma::cube betas(data.p,data.d,data.n_splines, arma::fill::ones);
  arma::vec beta_0(data.n_splines,arma::fill::zeros);
  arma::vec y_hat(data.y.n_rows,arma::fill::zeros);
  arma::vec y_hat_sum(data.y.n_rows,arma::fill::zeros);
  arma::vec partial_residuals(data.y.n_rows,arma::fill::zeros);
  arma::mat sum_splines(data.y.n_rows,data.n_splines,arma::fill::zeros);

  // Storing the posteriors
  int n_post = n_mcmc-n_burn;
  arma::mat beta_post(n_post,data.p,arma::fill::ones);
  arma::mat beta_0_post(n_post,data.n_splines,arma::fill::ones);
  arma::mat y_hat_post(n_post,data.y.n_rows,arma::fill::ones);

  arma::vec tau_post(n_post,arma::fill::ones);
  arma::mat tau_b_post(n_post,data.d,arma::fill::ones);
  arma::mat delta_post(n_post,data.d,arma::fill::ones);
  int post_iter = 0;


  // Initializing the sampling processes
  for(int i = 0; i < data.n_mcmc; i++){

    for(int j = 0; j < data.n_splines; j++){

          // Creating a vector for the single predictions
          arma::vec y_hat(data.y.n_rows,arma::fill::zeros);

          // Getting partial residuals
          partial_residuals = data.y - sum_exclude_col(sum_splines,j);

          // cout << "Beta error" << endl;
          beta_sampler(betas.slice(j),beta_0(j), data,partial_residuals);
          // cout << "Beta_0 error" << endl;
          beta_0_sampler(betas.slice(j), beta_0(j), data,partial_residuals);

          for(int d = 0 ; d<data.d;d++){
              y_hat = y_hat + data.B_train.slice(d)*(betas.slice(j)).col(d);
          }

          y_hat = y_hat + beta_0(j);

          sum_splines.col(j) = y_hat;
    }

    // Summing-up all predictions;
    y_hat_sum = sum(sum_splines,1);

    // cout << "Tau_b error" << endl;
    tau_b_sampler(betas,data);
    // // cout << "Tau_b value corresponds to: ";
    // for(int t = 0; t<data.d;t++){
    //   cout << data.tau_b(t) << " ";
    // }
    // cout << endl;
    // cout << "Delta error" << endl;
    delta_sampler(data);


    // cout << "Tau sampler error" << endl;
    tau_sampler(data,y_hat_sum);


    // cout << " No error" << endl;
    // Iterating and storing the posterior samples
    if(i >= data.n_burn){
      // beta_post.row(post_iter) = betas.t();
      beta_0_post.row(post_iter) = beta_0.t();
      y_hat_post.row(post_iter) = y_hat_sum.t();
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

