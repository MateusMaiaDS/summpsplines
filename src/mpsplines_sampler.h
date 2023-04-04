#include<RcppArmadillo.h>
#include<vector>
// Creating the struct
struct modelParam;

struct modelParam{

  arma::mat y;
  arma::cube B_train;

  // BART prior param specification
  int p;
  int d;
  int n_splines;
  arma::vec tau_b;
  double tau_b_intercept;
  double tau;
  double a_tau;
  double d_tau;
  double nu;
  arma::vec delta;
  double a_delta;
  double d_delta;
  arma::mat P;
  // MCMC spec.
  int n_mcmc;
  int n_burn;


  // Objects from the sampler calculator
  arma::mat bt_ones;

  // Defining the constructor for the model param
  modelParam(arma::cube B_train_,
             arma::mat y_,
             int n_splines,
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
             int n_burn_);

};
