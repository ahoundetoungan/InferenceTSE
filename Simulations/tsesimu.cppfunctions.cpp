// This file contains C++ functions that are used in the simulation study.

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;


double pi = acos(-1);

// [[Rcpp::export]]
// This function computes the cdf at x taking data as simulation from the distribution
// n in the length of data and k the length of x
arma::vec fcdf(const arma::vec& data, 
               arma::vec& x, 
               int& n, 
               int& k){
  arma::vec out(k, arma::fill::zeros);
  for(int i(0); i < n; ++ i){
    out.elem(arma::find(data(i) <= x)) += 1;
  }
  return out/n;
}

// [[Rcpp::export]]
// This function simulates data, given z
arma::mat filterdata(const int &N, 
                     const int &S,
                     const arma::vec &beta, 
                     const arma::vec &phi, 
                     const arma::mat &z) {
  arma::mat dt(N, S);
  arma::mat hi = arma::ones<arma::rowvec>(S)*beta(0)/(1 - beta(1) - beta(2)); //this is sigma^2 at i = 0 
  dt.row(0)    = phi(0)/(1 - phi(1)) + sqrt(hi)%z.row(0); // y at i = 0 because y at i - 1 = phi(0)/(1 - phi(1)), the mean
  for(int i(1); i < N; i++){
    hi         = beta(0) + beta(1)*pow(z.row(i-1), 2)%hi + beta(2)*hi;
    dt.row(i)  = phi(0) + phi(1)*dt.row(i-1) + sqrt(hi)%z.row(i);
  }
  return dt;
}

// [[Rcpp::export]]
// This function computes innovation z, given data
arma::vec filterz(const int &N, 
                  const arma::vec &beta, 
                  const arma::vec &phi, 
                  const arma::vec &dt) {
  arma::vec z(N);
  double hi(beta(0)/(1 - beta(1) - beta(2)));
  z(0)   = (dt(0) - phi(0)/(1 - phi(1)))/sqrt(hi);// because y at i - 1 = phi(0)/(1 - phi(1)), the mean
  for(int i(1); i < N; i++){
    hi   = beta(0) + beta(1)*pow(z(i-1), 2)*hi + beta(2)*hi;
    z(i) = (dt(i) - phi(0) - phi(1)*dt(i-1))/sqrt(hi);
  }
  return z;
}

// [[Rcpp::export]]
// This function computes z and sigma, given data
List filterzsigma(const int &N, 
                  const arma::vec &beta, 
                  const arma::vec &phi, 
                  const arma::vec &dt) {
  arma::vec z(N), sig(N);
  double hi(beta(0)/(1 - beta(1) - beta(2)));
  sig(0)   = sqrt(hi);
  z(0)     = (dt(0) - phi(0)/(1 - phi(1)))/sig(0);// because y at i - 1 = phi(0)/(1 - phi(1)), the mean
  for(int i(1); i < N; i++){
    hi     = beta(0) + beta(1)*pow(z(i-1), 2)*hi + beta(2)*hi;
    sig(i) = sqrt(hi);
    z(i)   = (dt(i) - phi(0) - phi(1)*dt(i-1))/sig(i);
  }
  return List::create(_["sigma"] = sig, _["z"] = z);
}

// This function computes unconstrained GARCH parameters as follows
// The parameter without constraint are
// phi1tilde = log(1 + phi1) - log(1 - phi1)
// omegatatilde = log(omega)
// alpha1tilde = log(alpha1) - log(1 - alpha1)
// beta1tilde = log(beta1) - log(1 - alpha1 - beta1)
// nutilde = log(nu - 2) if student
// The other parameters are without constraints
//[[Rcpp::export]]
arma::vec fbetat(const arma::vec& beta){
  double phi1t   = log(1 + beta(1)) - log(1 - beta(1));
  double omegat  = log(beta(2));
  double alpha1t = log(beta(3)) - log(1 - beta(3));
  double beta1t  = log(beta(4)) - log(1 - beta(3) - beta(4));
  
  arma::vec out;
  if(beta.n_elem == 6){
    double nut   = log(beta(5) - 2);
    out          = {beta(0), phi1t, omegat, alpha1t, beta1t, nut};
  } else{
    out          = {beta(0), phi1t, omegat, alpha1t, beta1t};
  }
  
  return out;
}

//[[Rcpp::export]]
// This function computes constrained GARCH parameters
arma::vec fbeta(const arma::vec& betat){
  double phi1    = (1 - exp(-betat(1)))/(1 + exp(-betat(1)));
  double omega   = exp(betat(2));
  double alpha1  = 1/(1 + exp(-betat(3)));
  double beta1   = (1 - alpha1)/(1 + exp(-betat(4)));
  
  
  arma::vec out;
  if(betat.n_elem == 6){
    double nu    = exp(betat(5)) + 2;
    out          = {betat(0), phi1, omega, alpha1, beta1, nu};
  } else{
    out          = {betat(0), phi1, omega, alpha1, beta1};
  }
  
  return out;
}


// This function computes The garch(1, 1) likelihood for each i with student distribution
//[[Rcpp::export]]
arma::vec fllhigarcht(const arma::vec& betat,
                      const int &N,
                      const arma::vec& dt){
  arma::vec beta = fbeta(betat);
  double nu      = beta(5);
  double scal    = sqrt(nu/(nu - 2));
  double lscal   = log(scal);
  arma::vec lh(N);
  
  // Compute sigma and z
  double h = beta(2)/(1 - beta(3) - beta(4));
  double z = (dt(0) - beta(0)/(1 - beta(1)))/sqrt(h);
  lh(0)    = R::dt(scal*z, nu, 1) + lscal - 0.5*log(h);
  for(int i(1); i < N; i++){
    h      = beta(2) + beta(3)*pow(z, 2)*h + beta(4)*h;
    z      = (dt(i) - beta(0) - beta(1)*dt(i-1))/sqrt(h);
    lh(i)  = R::dt(scal*z, nu, 1) + lscal - 0.5*log(h);
  }
  return lh;
}


// This function computes minus likelihood of garch(1, 1)
//[[Rcpp::export]]
double fllhgarcht(const arma::vec& betat,
                  const int &N,
                  const arma::vec& dt){
  return -sum(fllhigarcht(betat, N, dt));
}

// This function computes The garch(1, 1) likelihood for each i with normal distribution
//[[Rcpp::export]]
arma::vec fllhigarchn(const arma::vec& betat,
                      const int &N,
                      const arma::vec& dt){
  arma::vec beta = fbeta(betat);
  arma::vec lh(N);
  
  // Compute sigma and z
  double h = beta(2)/(1 - beta(3) - beta(4));
  double z = (dt(0) - beta(0)/(1 - beta(1)))/sqrt(h);
  lh(0)    = R::dnorm4(z, 0, 1, 1)- 0.5*log(h);
  for(int i(1); i < N; i++){
    h      = beta(2) + beta(3)*pow(z, 2)*h + beta(4)*h;
    z      = (dt(i) - beta(0) - beta(1)*dt(i-1))/sqrt(h);
    lh(i)  = R::dnorm4(z, 0, 1, 1)- 0.5*log(h);
  }
  return lh;
}


// This function computes minus likelihood of garch(1, 1)
//[[Rcpp::export]]
double fllhgarchn(const arma::vec& betat,
                  const int &N,
                  const arma::vec& dt){
  return -sum(fllhigarchn(betat, N, dt));
}

// This function compute the covariance matrix of unconstrained GARCH parameters using Delta method
//[[Rcpp::export]]
arma::mat fcovbetat(const arma::vec& beta, 
                    const arma::mat& covbeta,
                    const int& S){
  arma::mat jac = arma::eye<arma::mat>(6*S, 6*S);
  for(int s(0); s < S; ++ s){
    int tp = 6*s;
    jac(tp + 1, tp + 1) = 2/(1 - pow(beta(tp + 1), 2));
    jac(tp + 2, tp + 2) = 1/beta(tp + 2);
    jac(tp + 3, tp + 3) = 1/(beta(tp + 3)*(1 - beta(tp + 3)));
    jac(tp + 4, tp + 3) = 1/(1 - beta(tp + 3) - beta(tp + 4));
    jac(tp + 4, tp + 4) = (1 - beta(tp + 3))/(beta(tp + 4)*(1 - beta(tp + 3) - beta(tp + 4)));
    jac(tp + 5, tp + 5) = 1/(beta(tp + 5) - 2);
  }
  return jac*covbeta*jac.t();
}

// This function computes cdf of yi for each i with the t distribution
//[[Rcpp::export]]
arma::vec fGigarcht(const arma::vec& betat,
                    const int &N,
                    const arma::vec& dt){
  arma::vec beta = fbeta(betat);
  double nu      = beta(5);
  double scal    = sqrt(nu/(nu - 2));
  arma::vec out(N);
  
  // Compute sigma and z
  double h = beta(2)/(1 - beta(3) - beta(4));
  double z = (dt(0) - beta(0)/(1 - beta(1)))/sqrt(h);
  out(0)   = R::pt(scal*z, nu, 1, 0);
  for(int i(1); i < N; i++){
    h      = beta(2) + beta(3)*pow(z, 2)*h + beta(4)*h;
    z      = (dt(i) - beta(0) - beta(1)*dt(i-1))/sqrt(h);
    out(i) = R::pt(scal*z, nu, 1, 0);
  }
  return out;
}

// This function computes cdf of yi for each i with normal distribution
//[[Rcpp::export]]
arma::vec fGigarchn(const arma::vec& betat,
                    const int &N,
                    const arma::vec& dt){
  arma::vec beta = fbeta(betat);
  arma::vec out(N);
  
  // Compute sigma and z
  double h = beta(2)/(1 - beta(3) - beta(4));
  double z = (dt(0) - beta(0)/(1 - beta(1)))/sqrt(h);
  out(0)   = R::pnorm5(z, 0, 1, 1, 0);
  for(int i(1); i < N; i++){
    h      = beta(2) + beta(3)*pow(z, 2)*h + beta(4)*h;
    z      = (dt(i) - beta(0) - beta(1)*dt(i-1))/sqrt(h);
    out(i) = R::pnorm5(z, 0, 1, 1, 0);
  }
  return out;
}


//[[Rcpp::export]]
//log Clayton pdf for i for dimension n
arma::vec fllhiclayton(arma::mat& u, double& theta, int& dim){
  arma::vec tmp1 = arma::linspace(1, dim - 1, dim - 1);
  arma::vec tmp2 = sum(pow(u, -theta), 1) - dim + 1;
  return accu(log(tmp1*theta + 1)) - (theta + 1)*sum(log(u), 1) - (1/theta + dim)*log(tmp2);
}

//[[Rcpp::export]]
//minus log Clayton pdf for dimension n
double fllhclayton(arma::mat& u, double& ltheta, int& dim){
  double theta  = exp(ltheta);
  arma::vec tmp = sum(pow(u, -theta), 1) - dim + 1;
  if(tmp.min() < 0) return R_PosInf;
  return -sum(fllhiclayton(u, theta, dim));
}

//[[Rcpp::export]]
//First derivative of claytonpdf for dimension n
arma::vec d1lclaytonpdf(arma::mat& u, double& ltheta, int& dim){
  double theta  = exp(ltheta);
  arma::mat lu   = log(u);
  arma::vec tmp1 = arma::linspace(1, dim - 1, dim - 1);
  arma::mat uthe = pow(u, -theta);
  arma::vec tmp2 = sum(uthe, 1) - dim + 1;
  return ((dim + 1/theta)*sum(uthe%lu, 1)/tmp2 - sum(lu, 1) + 
          sum(tmp1/(tmp1*theta + 1)) + log(tmp2)/pow(theta, 2))*theta;
}  

//[[Rcpp::export]]
//Second derivative of claytonpdf for dimension n
arma::vec d2lclaytonpdf(arma::mat& u, double& ltheta, int& dim){
  double theta   = exp(ltheta);
  arma::mat lu   = log(u);
  arma::vec tmp1 = arma::linspace(1, dim - 1, dim - 1);
  arma::mat uthe = pow(u, -theta);
  arma::mat uthl = pow(u, -theta)%lu;
  arma::vec tmp2 = sum(uthe, 1) - dim + 1;
  arma::vec tmp3 = sum(uthl, 1);
  arma::vec tmp4 = sum(uthl%lu, 1);
  arma::vec tmp5 = sum(pow(uthl, 2), 1);
  return ((dim + 1/theta)*(pow(tmp3/tmp2, 2) - tmp4/tmp2) - 
          sum(pow(tmp1, 2)/pow(tmp1*theta + 1, 2)) - 2*tmp3/(pow(theta, 2)*tmp2) -
          2*log(tmp2)/pow(theta, 3))*pow(theta, 2) + d1lclaytonpdf(u, ltheta, dim)*theta;
}   


//[[Rcpp::export]]
//Second derivative of claytonpdf for dimension n with respect to theta and u
arma::mat dud1lclaytonpdf(arma::mat& u, double& ltheta, int& dim){
  double theta    = exp(ltheta);
  arma::mat lu    = log(u);
  arma::mat uthe  = pow(u, -theta);
  arma::mat utheo = pow(u, -theta - 1.0);
  arma::vec tmp2  = sum(uthe, 1) - dim + 1;
  arma::mat tmp3  = utheo.each_col()%((theta*dim + 1)*sum(uthe%lu, 1)/arma::square(tmp2) - 1/(theta*tmp2));
  arma::mat tmp4  = (utheo.each_col()%((dim + 1/theta)/tmp2))%(1 - theta*lu);
  return (-1/u + tmp3 + tmp4)*theta;
}   


//HAC estimation 
//Andrews, D. W. (1991). Heteroskedasticity and autocorrelation consistent covariance matrix estimation. 
//Econometrica: Journal of the Econometric Society, 817-858.
//Kernels 
double kTR(double x){
  if(abs(x) <= 1) return 1;
  return 0;
}

double kBT(double x){
  double abx = abs(x);
  if(abx <= 1) return 1 - abx;
  return 0;
}

double kPR(double x){
  double abx = abs(x);
  if(abx <= 0.5) return 1 -6*pow(x, 2) + 6*pow(abx, 3);
  if(abx <= 1) return 2*pow(1 - abx, 3);
  return 0;
}


double kTH(double x){
  if(abs(x) <= 0.5) return (1 + cos(pi*x))/2; 
  return 0;
}


//
double kQS(double x){
  double pix  = pi*x;
  double tmp1 = 12*pow(pix, 2);
  double tmp2 = 6*pix/5;
  return (25/tmp1) * (sin(tmp2)/tmp2 - cos(tmp2));
}

//[[Rcpp::export]]
arma::mat fHAC(const arma::mat& V, 
               const int& n, //nrow of V
               const int& r = 1, //number of parameters estimated
               const int& kernel = 5, //Quadratic spectral kernel
               const int& st = 0){ //default bandwith = ceil(0.75*pow(n, 1.0/3))
  double ST = st;
  if(st == 0){
    ST  = ceil(0.75*pow(n, 1.0/3));
  }
  arma::mat out = V.t()*V/n;
  arma::mat tp;
  
  switch(kernel)
  {
  case 1:
    for(int j(1); j < n; ++j){
      tp    = arma::trans(V.rows(0, n - 1 - j))*V.rows(j, n - 1)/n;
      out  += (kTR(j/ST)*(tp + tp.t()));
    }
    break;
  case 3:
    for(int j(1); j < n; ++j){
      tp    = arma::trans(V.rows(0, n - 1 - j))*V.rows(j, n - 1)/n;
      out  += (kPR(j/ST)*(tp + tp.t()));
    }
    break;
  case 4:
    for(int j(1); j < n; ++j){
      tp    = arma::trans(V.rows(0, n - 1 - j))*V.rows(j, n - 1)/n;
      out  += (kTH(j/ST)*(tp + tp.t()));
    }
    break;
  case 5:
    for(int j(1); j < n; ++j){
      tp    = arma::trans(V.rows(0, n - 1 - j))*V.rows(j, n - 1)/n;
      out  += (kQS(j/ST)*(tp + tp.t()));
    }
    break;
  default:
    for(int j(1); j < n; ++j){
      tp    = arma::trans(V.rows(0, n - 1 - j))*V.rows(j, n - 1)/n;
      out  += (kBT(j/ST)*(tp + tp.t()));
    }
    break;
  }
  return ((double)n/(n - r))*out;
}

//[[Rcpp::export]]
//This function compute weight for randomForest
arma::mat comp_weights(const arma::mat& lnew,
                       const arma::mat& ltrain,
                       const arma::mat& bagtree,
                       const int& ntrain,
                       const int& nnew,
                       const int& ntree){
  arma::mat out(ntrain, nnew, arma::fill::zeros);
  for(int t(0); t < ntree; ++ t){
    arma::umat tp1(ntrain, nnew);
    for(int s(0); s < nnew; ++ s){
      tp1.col(s) = (ltrain.col(t) == lnew(s, t));
    }
    arma::mat tp2(arma::conv_to<arma::mat>::from(tp1));
    out += arma::normalise(tp2.each_col()%bagtree.col(t), 1, 0);
  }
  return out.t()/ntree;
}

/***R
# Function to calculate weights for a given new observation x
fweights  <- function(model, train_x, new_x) {
  # Get leaf nodes for the new data point across all trees
  leaves_new   <- attr(predict(model, newdata = new_x, nodes = TRUE), "nodes")

  # Get leaf nodes for all training data
  leaves_train <- attr(predict(model, newdata = train_x, nodes = TRUE),"nodes")

  # Weights
  comp_weights(lnew = leaves_new, ltrain = leaves_train, bagtree = model$inbag,
               ntrain = nrow(train_x), nnew = nrow(new_x), ntree = model$ntree)
}
*/

//[[Rcpp::export]]
//This function compute many randomForest estimator by sampling the trees
arma::mat fsimRForest(const arma::mat& treepredict, 
                      const int& ntree,
                      const int& kappa,
                      const int& n){
  arma::mat out(n, kappa);
  for(int s(0); s < kappa; ++ s){
    out.col(s) = mean(treepredict.cols(randi<uvec>(ntree, distr_param(0, ntree - 1))), 1);
  }
  return out;
}
