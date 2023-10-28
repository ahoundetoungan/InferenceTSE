//' This file contains auxiliary functions for the application on peer effects.
//' Most functions come from the package PartialNetwork by Boucher and Houndetoungan.
//' We have adapted these functions to our framework.

// [[Rcpp::depends(RcppArmadillo, RcppEigen, RcppNumerical, RcppProgress)]]
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <RcppEigen.h>
#include <RcppNumerical.h>

typedef Eigen::Map<Eigen::MatrixXd> MapMatr;
typedef Eigen::Map<Eigen::VectorXd> MapVect;

using namespace Rcpp;
using namespace arma;
using namespace Numer;
using namespace std;
const double pi = acos(-1);

// A class to compute the objective function of the network formation model using censored Poisson
// In this implementation, we refer to the first stage estimator as "rho" instead 
// of "beta," as indicated in the paper.
class CENetPoi: public MFuncGrad
{
private:
  const arma::vec& d;
  const arma::mat& x;
  const arma::uvec& censure;
  const arma::vec& Nicum;
  const arma::vec& group;
  const int& n;
  const int& M;
  const int& Kx;
public:
  CENetPoi(const arma::vec& d_,
          const arma::mat& x_,
          const arma::uvec& censure_,
          const arma::vec& Nicum_,
          const arma::vec& group_,
          const int& n_,
          const int& M_,
          const int& Kx_) :
  d(d_),
  x(x_),
  censure(censure_),
  Nicum(Nicum_),
  group(group_),
  n(n_),
  M(M_),
  Kx(Kx_){}
  
  arma::vec Grad;
  
  double f_grad(Constvec& rho, Refvec grad)
  {
    Eigen::VectorXd rho0 = rho;  //make a copy
    arma::vec rhoa       = arma::vec(rho0.data(), rho0.size(), false, false); //converte into arma vec
    
    Rcpp::Rcout << "beta: \n";
    cout << arma::trans(rhoa) << endl;
    
    int nparms(M + Kx);
    double llh(0);
    arma::vec xb(x*rhoa.tail(Kx)), gd(nparms, arma::fill::zeros);

    for (int i(0); i < n; ++ i) {
      int n1          = Nicum(i);
      int n2          = Nicum(i + 1) - 1;
      int gi          = group(i);
      arma::vec exbm  = exp(-(rhoa(gi) + xb.subvec(n1, n2)));
      arma::vec tmp1  = 1/(1 + exbm);
      // arma::vec tmp1 = (rhoa(gi) + xb.subvec(n1, n2)); //This is for OLS
      arma::vec tmp2  =  exbm%pow(tmp1, 2);
      // arma::vec tmp2 =  arma::ones(n2 - n1 + 1); //This is for OLS
      double ld       = sum(tmp1);
      arma::vec dld(nparms, arma::fill::zeros);
      dld(gi)         = sum(tmp2);
      dld.tail(Kx)    = x.rows(n1, n2).t()*tmp2;
      
      if(censure(i) == 1){
        int rcens = n2 - n1 + 1;
        double l2 = R::ppois(rcens, ld, true, true);
        double l1 = R::ppois(d(i), ld, true, true);
        double tp = l2 + log(1 - exp(l1 - l2));
        llh      += tp;
        gd       +=  dld*(R::dpois(d(i), ld, false) - R::dpois(rcens, ld, false))*exp(-tp);
      } else {
        llh      += R::dpois(d(i), ld, true);
        gd       += d(i)*dld/ld -dld;
      }
    }
    if(!is_finite(sum(gd))){
      gd          = arma::ones(nparms)*1e8;
    }
    if(!is_finite(llh)){
      llh         = -1e8;
    }
    Rcpp::Rcout << "Gradient: " << "\n";
    cout<<gd.t()<<endl;
    
    grad = -Eigen::Map<Eigen::VectorXd>(gd.memptr(), nparms);
    Grad = gd;
    
    Rcpp::Rcout << "Objective: " << llh << "\n\n";
    Rcpp::Rcout << "***************\n";
    return -llh;
  }
};


// A class to compute the objective function of the network formation model using censored Normal
class CENetNor: public MFuncGrad
{
private:
  const arma::vec& d;
  const arma::mat& x;
  const arma::uvec& censure;
  const arma::vec& Nicum;
  const arma::vec& group;
  const int& n;
  const int& M;
  const int& Kx;
public:
  CENetNor(const arma::vec& d_,
           const arma::mat& x_,
           const arma::uvec& censure_,
           const arma::vec& Nicum_,
           const arma::vec& group_,
           const int& n_,
           const int& M_,
           const int& Kx_) :
  d(d_),
  x(x_),
  censure(censure_),
  Nicum(Nicum_),
  group(group_),
  n(n_),
  M(M_),
  Kx(Kx_){}
  
  arma::vec Grad;
  
  double f_grad(Constvec& rho, Refvec grad)
  {
    Eigen::VectorXd rho0 = rho;  //make a copy
    arma::vec rhoa       = arma::vec(rho0.data(), rho0.size(), false, false); //converte into arma vec
    
    Rcpp::Rcout << "beta: \n";
    cout << arma::trans(rhoa) << endl;
    
    int nparms(M + Kx + 1);
    double llh(0), se2(exp(rhoa(M + Kx))), se(sqrt(se2));
    arma::vec xb(x*rhoa.tail(Kx + 1).head(Kx)), gd(nparms, arma::fill::zeros);

    for (int i(0); i < n; ++ i) {
      int n1          = Nicum(i);
      int n2          = Nicum(i + 1) - 1;
      int gi          = group(i);
      arma::vec exbm  = exp(-(rhoa(gi) + xb.subvec(n1, n2)));
      arma::vec tmp1  = 1/(1 + exbm);
      // arma::vec tmp1 = (rhoa(gi) + xb.subvec(n1, n2)); //This is for OLS
      arma::vec tmp2  =  exbm%pow(tmp1, 2);
      // arma::vec tmp2 =  arma::ones(n2 - n1 + 1); //This is for OLS
      double ld       = sum(tmp1);
      arma::vec dld(M + Kx, arma::fill::zeros);
      dld(gi)         = sum(tmp2);
      dld.tail(Kx)    = x.rows(n1, n2).t()*tmp2;
      
      if(censure(i) == 1){
        int rcens = n2 - n1 + 1;
        double l2 = R::pnorm5(rcens, ld, se, true, true);
        double l1 = R::pnorm5(d(i), ld, se, true, true);
        double tp = l2 + log(1 - exp(l1 - l2));
        llh      += tp;
        gd.head(M + Kx) +=  dld*(R::dnorm4(d(i), ld, se, false) - R::dnorm4(rcens, ld, se, false))*exp(-tp);
        gd(M + Kx)      +=  ((d(i) - ld)*R::dnorm4(d(i), ld, se, false) - (rcens - ld)*R::dnorm4(rcens, ld, se, false))*exp(-tp)/2;
      } else {
        llh             += R::dnorm4(d(i), ld, se, true);
        gd.head(M + Kx) += dld*(d(i) - ld)/se2;
        gd(M + Kx)      += (0.5*pow(d(i) - ld, 2)/se2 - 0.5);
      }
    }
    if(!is_finite(sum(gd))){
      gd          = arma::ones(nparms)*1e10;
    }
    if(!is_finite(llh)){
      llh         = -1e10;
    }
    Rcpp::Rcout << "Gradient: " << "\n";
    cout<<gd.t()<<endl;
    
    grad = -Eigen::Map<Eigen::VectorXd>(gd.memptr(), nparms);
    Grad = gd;
    
    Rcpp::Rcout << "Objective: " << llh << "\n\n";
    Rcpp::Rcout << "***************\n";
    return -llh;
  }
};


// Estimation of the network formation model using censored Normal
//[[Rcpp::export]]
List ENetPoi(Eigen::VectorXd rho,
             const arma::vec& d,
             const arma::mat& x,
             const int& bound,
             const arma::vec& Nicum,
             const arma::vec& group,
             const int& n,
             const int& M,
             const int maxit = 300, 
             const double& eps_f = 1e-6, 
             const double& eps_g = 1e-5){
  int Kx(x.n_cols);
  arma::uvec censure = (d == bound);
  CENetPoi f(d, x, censure, Nicum, group, n, M, Kx);
  
  double fopt;
  int status = optim_lbfgs(f, rho, fopt, maxit, eps_f, eps_g);
  
  return Rcpp::List::create(
    Rcpp::Named("estimate") = rho,
    Rcpp::Named("value")    = fopt,
    Rcpp::Named("gradient") = f.Grad,
    Rcpp::Named("status")   = status);
}


// Estimation of the network formation model using censored Normal
//[[Rcpp::export]]
List ENetNor(Eigen::VectorXd rho,
             const arma::vec& d,
             const arma::mat& x,
             const int& bound,
             const arma::vec& Nicum,
             const arma::vec& group,
             const int& n,
             const int& M,
             const int maxit = 300, 
             const double& eps_f = 1e-6, 
             const double& eps_g = 1e-5){
  int Kx(x.n_cols);
  arma::uvec censure = (d == bound);
  CENetNor f(d, x, censure, Nicum, group, n, M, Kx);
  
  double fopt;
  int status = optim_lbfgs(f, rho, fopt, maxit, eps_f, eps_g);
  
  return Rcpp::List::create(
    Rcpp::Named("estimate") = rho,
    Rcpp::Named("value")    = fopt,
    Rcpp::Named("gradient") = f.Grad,
    Rcpp::Named("status")   = status);
}


// Starting point for the censored Normal
//[[Rcpp::export]]
arma::vec StartNetNor(arma::vec rho,
                      const arma::vec& d,
                      const arma::mat& x,
                      const int& bound,
                      const arma::vec& Nicum,
                      const arma::vec& group,
                      const int& n,
                      const int& M){
  int Kx(x.n_cols);
  arma::uvec censure(d == bound);
  arma::vec xb(x*rho.tail(Kx)), out(M + Kx + 1);
  out.head(M + Kx)    = rho;
  double se2(0);
  
  for (int i(0); i < n; ++ i) {
    if(censure(i) == 0){
      int n1          = Nicum(i);
      int n2          = Nicum(i + 1) - 1;
      int gi          = group(i);
      arma::vec exbm  = exp(-(rho(gi) + xb.subvec(n1, n2)));
      arma::vec tmp1  = 1/(1 + exbm);
      double ld       = sum(tmp1);
      se2            += pow(d(i) - ld, 2);
    }
  }
  out(M + Kx)         = log(se2/sum(censure == 0));
  return out;
}

// Polish the estimation of the network formation model using the NR algorithm
// Poisson model
//[[Rcpp::export]]
List NRNetPoi(const arma::vec& rho,
              const arma::vec& d,
              const arma::mat& x,
              const int& bound,
              const arma::vec& Nicum,
              const arma::vec& group,
              const int& n,
              const int& M,
              const int& maxit  = 100,
              const double& tol = 1e-5){
  int Kx(x.n_cols), nparms(M + Kx), iter(0);
  arma::uvec censure = (d == bound);
  double dist(R_PosInf), loglike;
  arma::vec rhk(rho);
  arma::mat Hess, grad;
  NumericVector rhR;
  
  while((iter < maxit) && (dist >= tol)) {
    ++ iter;
    double llh(0);
    arma::vec xb(x*rhk.tail(Kx));
    arma::mat gd(nparms, n, arma::fill::zeros), he(nparms, nparms, arma::fill::zeros);
    Progress prog(n, true);
    for (int i(0); i < n; ++ i) {
      int n1           = Nicum(i);
      int n2           = Nicum(i + 1) - 1;
      int gi           = group(i);
      arma::vec exbm   = exp(-(rhk(gi) + xb.subvec(n1, n2)));
      arma::vec tmp1   = 1/(1 + exbm);
      arma::vec tmp12  = pow(tmp1, 2);
      arma::vec tmp13  = pow(tmp1, 3);
      arma::vec tmp2   = exbm%tmp12;
      arma::vec tmp3   = -tmp2 + 2*pow(exbm, 2)%tmp13;
      
      // lambda
      double ld(sum(tmp1));
      
      // dlambda/drho
      arma::mat xi(n2 - n1 + 1, nparms, arma::fill::zeros);
      xi.tail_cols(Kx) = x.rows(n1, n2);
      xi.col(gi)      += 1;
      arma::vec dld    = xi.t()*tmp2;
      
      // ddlambda/ddrho
      arma::mat ddld   = xi.t()*(xi.each_col()%tmp3);
      
      if(censure(i) == 1){
        int rcens(n2 - n1 + 1);
        double l2(R::ppois(rcens, ld, true, true));
        double l1(R::ppois(d(i), ld, true, true));
        double tmp4(l2 + log(1 - exp(l1 - l2)));
        double tmp5(R::dpois(d(i), ld, false) - R::dpois(rcens, ld, false));
        double tmp6(R::dpois(d(i) - 1, ld, false) - R::dpois(rcens - 1, ld, false) - tmp5);
        double tmp7(tmp5*exp(-tmp4));
        llh      += tmp4;
        gd.col(i) =  dld*tmp7;
        he       += (ddld*tmp7 + (dld*dld.t())*(tmp6*exp(-tmp4) - pow(tmp7, 2)));
      } else {
        llh      += R::dpois(d(i), ld, true);
        double dl(d(i)/ld - 1), ddl(-d(i)/pow(ld, 2));
        gd.col(i) = dld*dl;
        he       += (ddld*dl + (dld*dld.t())*ddl);
      }
      prog.increment();
    }
    
    arma::vec rhl(rhk - solve(he, sum(gd, 1)));
    dist            = sum(abs(rhl - rhk));
    rhk             = rhl; 
    loglike         = llh;
    grad            = gd;
    Hess            = he;
    rhR             = wrap(rhk);
    rhR.attr("dim") = R_NilValue;
    Rcpp::print(rhR);
    Rprintf("iteration: %d *** loglik: %f *** distance: %f \n", iter, llh, dist);
  }
  
  arma::mat tmp8(grad*grad.t());
  arma::mat tmp9(arma::inv(Hess));
  arma::mat covm(tmp9*tmp8*tmp9.t());
  
  return List::create(Named("estimate")  = rhR,
                      Named("grad")      = grad,
                      Named("Hess")      = Hess,
                      Named("loglik")    = loglike,
                      Named("iteration") = iter,
                      Named("distance")  = dist,
                      Named("covmat")    = covm);
}

// Normal model
//[[Rcpp::export]]
List NRNetNor(const arma::vec& rho,
              const arma::vec& d,
              const arma::mat& x,
              const int& bound,
              const arma::vec& Nicum,
              const arma::vec& group,
              const int& n,
              const int& M,
              const int& maxit  = 100,
              const double& tol = 1e-5){
  int Kx(x.n_cols), nparms(M + Kx + 1), iter(0);
  arma::uvec censure = (d == bound);
  double dist(R_PosInf), loglike;
  arma::vec rhk(rho);
  arma::mat Hess, grad;
  NumericVector rhR;

  while((iter < maxit) && (dist >= tol)) {
    ++ iter;
    double llh(0), se2(exp(rhk(M + Kx))), se(sqrt(se2));
    arma::vec xb(x*rhk.tail(Kx + 1).head(Kx));
    arma::mat gd(nparms, n, arma::fill::zeros), he(nparms, nparms, arma::fill::zeros);
    Progress prog(n, true);
    for (int i(0); i < n; ++ i) {
      int n1           = Nicum(i);
      int n2           = Nicum(i + 1) - 1;
      int gi           = group(i);
      arma::vec exbm   = exp(-(rhk(gi) + xb.subvec(n1, n2)));
      arma::vec tmp1   = 1/(1 + exbm);
      arma::vec tmp12  = pow(tmp1, 2);
      arma::vec tmp13  = pow(tmp1, 3);
      arma::vec tmp2   = exbm%tmp12;
      arma::vec tmp3   = -tmp2 + 2*pow(exbm, 2)%tmp13;

      // lambda
      double ld(sum(tmp1));

      // dlambda/drho
      arma::mat xi(n2 - n1 + 1, M + Kx, arma::fill::zeros);
      xi.tail_cols(Kx) = x.rows(n1, n2);
      xi.col(gi)      += 1;
      arma::vec dld    = xi.t()*tmp2;

      // ddlambda/ddrho
      arma::mat ddld   = xi.t()*(xi.each_col()%tmp3);
      
      if(censure(i) == 1){
        // cout<<1<<endl;
        int rcens(n2 - n1 + 1);
        double l2 = R::pnorm5(rcens, ld, se, true, true);
        double l1 = R::pnorm5(d(i), ld, se, true, true);
        double tmp4(l2 + log(1 - exp(l1 - l2)));
        double tmp5b(R::dnorm4(d(i), ld, se, false) - R::dnorm4(rcens, ld, se, false));
        double tmp5s(((d(i) - ld)*R::dnorm4(d(i), ld, se, false) -
                     (rcens - ld)*R::dnorm4(rcens, ld, se, false))/2);
        double tmp6bb(2*tmp5s/se2);
        double tmp6ss((pow(d(i) - ld, 3)*R::dnorm4(d(i), ld, se, false) -
                      pow(rcens - ld, 3)*R::dnorm4(rcens, ld, se, false))/(4*se2));
        double tmp6bs((pow(d(i) - ld, 2)*R::dnorm4(d(i), ld, se, false) -
                      pow(rcens - ld, 2)*R::dnorm4(rcens, ld, se, false))/(2*se2));
        double tmp7b(tmp5b*exp(-tmp4));
        double tmp7s(tmp5s*exp(-tmp4));
        
        llh                   += tmp4;
        gd.submat(0, i, M + Kx - 1, i) = dld*tmp7b;
        gd(M + Kx, i)                  = tmp7s;
        he.submat(0, 0, M + Kx - 1, M + Kx -1) += (ddld*tmp7b + 
          (dld*dld.t())*(tmp6bb*exp(-tmp4) - pow(tmp7b, 2)));
        he(M + Kx, M + Kx)    += (-tmp7s/2 + tmp6ss*exp(-tmp4) - pow(tmp7s, 2));
        arma::vec tmp8         = dld*(-tmp7b/2 + tmp6bs*exp(-tmp4) - tmp7b*tmp7s);
        he.submat(0, M + Kx, M + Kx - 1, M + Kx) += tmp8;
        he.submat(M + Kx, 0, M + Kx, M + Kx - 1) += tmp8.t();
      } else {
        // cout<<2<<endl;
        llh      += R::dnorm4(d(i), ld, se, true);
        double dl1((d(i) - ld)/se2);
        gd.submat(0, i, M + Kx - 1, i) = dld*dl1;
        gd(M + Kx, i)                  = 0.5*pow(d(i) - ld, 2)/se2 - 0.5;
        he.submat(0, 0, M + Kx - 1, M + Kx - 1) += (ddld*dl1 - (dld*dld.t())/se2);
        he(M + Kx, M + Kx)    += (-pow(d(i) - ld, 2)/(2*se2));
        arma::vec tmp8         = dld*(-(d(i) - ld)/se2);
        he.submat(0, M + Kx, M + Kx - 1, M + Kx) += tmp8;
        he.submat(M + Kx, 0, M + Kx, M + Kx - 1) += tmp8.t();
      }
      prog.increment();
    }

    arma::vec rhl(rhk - solve(he, sum(gd, 1)));
    dist            = sum(abs(rhl - rhk));
    rhk             = rhl;
    loglike         = llh;
    grad            = gd;
    Hess            = he;
    rhR             = wrap(rhk);
    rhR.attr("dim") = R_NilValue;
    Rcpp::print(rhR);
    Rprintf("iteration: %d *** loglik: %f *** distance: %f \n", iter, llh, dist);
  }

  arma::mat tmp8(grad*grad.t());
  arma::mat tmp9(arma::inv(Hess));
  arma::mat covm(tmp9*tmp8*tmp9.t());

  return List::create(Named("estimate")  = rhR,
                      Named("grad")      = grad,
                      Named("Hess")      = Hess,
                      Named("loglik")    = loglike,
                      Named("iteration") = iter,
                      Named("distance")  = dist,
                      Named("covmat")    = covm);
}

// Estimation of the distribution of P
//[[Rcpp::export]]
arma::vec estimP(const arma::vec& rhofemale,
                 const arma::vec& rhomale,
                 const arma::mat& X,
                 const arma::uvec& groupj,
                 const arma::vec& femalej,
                 const arma::vec& G,
                 const arma::vec& Gobs,
                 const int& Kx,
                 const int& M){
  arma::vec xb = (rhofemale.elem(groupj) + X*rhofemale.subvec(M, M + Kx - 1))%femalej +
    (rhomale.elem(groupj) + X*rhomale.subvec(M, M + Kx - 1))%(1 - femalej);
  arma::vec P  = 1/(1 + exp(-xb));
  return Gobs%G + (1 - Gobs)%P;
}


// Estimation of the distribution of P when rho is simulated
//[[Rcpp::export]]
arma::vec estimPsim(const arma::vec& rho,
                    const arma::mat& cholcovrho,
                    const arma::mat& X,
                    const arma::uvec& groupj,
                    const arma::vec& femalej,
                    const arma::vec& G,
                    const arma::vec& Gobs,
                    const int& Kx,
                    const int& M,
                    const int& nparms){
  //simulate rhofemale and rhomale
  arma::vec tp1(arma::randn(2*nparms));
  tp1          = cholcovrho*tp1 + rho;
  arma::vec rhofemale(tp1.head(nparms)), rhomale(tp1.tail(nparms));
  arma::vec xb = (rhofemale.elem(groupj) + X*rhofemale.subvec(M, M + Kx - 1))%femalej +
    (rhomale.elem(groupj) + X*rhomale.subvec(M, M + Kx - 1))%(1 - femalej);
  arma::vec P  = 1/(1 + exp(-xb));
  return Gobs%G + (1 - Gobs)%P;
}


/***R
# This is the function summary.smmSAR adapted to our approach
simudist <- function(object, 
                     .fun, 
                     .args, 
                     sim    = 30,
                     dnetwork, 
                     data,
                     ...){
  stopifnot(inherits(object, "smmSAR"))
  details       <- object$details
  derM          <- details$av.grad.m
  aveMM         <- details$`av.m%*%t(m)`
  aveM          <- details$av.m
  iv.power      <- details$iv.power
  cond.var      <- !(is.null(derM)|is.null(aveMM)|is.null(aveM))
  W             <- details$W
  formula       <- object$formula
  fixed.effects <- object$fixed.effects
  contextual    <- object$contextual

  if(!missing(.args) & missing(.fun)){
    stop("`.args` is defined while `.fun` is missing.")
  }
  
  N             <- object$N
  M             <- object$n.group
  if(M < 2) stop("Inference is not possible with one group")
  Nsum          <- sum(N)
  
  seed          <- details$seed
  sysSeed       <- .GlobalEnv$.Random.seed
  on.exit({
    if (is.null(sysSeed)) {
      rm(".Random.seed", envir = .GlobalEnv)
    } else {
      .GlobalEnv$.Random.seed <- sysSeed 
    }
  })
  
  # data
  f.t.data     <- PartialNetwork:::formula.to.data.smm(formula = formula, data = data, fixed.effects = fixed.effects) 
  X1           <- f.t.data$X
  y            <- f.t.data$y
  GX2          <- f.t.data$GX
  Gy           <- f.t.data$Gy
  GXobs        <- !is.null(GX2)
  Gyobs        <- !is.null(Gy)
  col.x1       <- colnames(X1)
  col.gx2      <- colnames(GX2)
  intercept    <- ("(Intercept)" %in% col.x1)
  Kx1          <- length(col.x1) 
  if(!GXobs){
    col.gx2    <- paste0("G: ", col.x1[(1 + intercept):Kx1])
  }
  Kx2          <- length(col.gx2)  
  if((Kx1 - intercept) != Kx2) stop("The number of observed contextual variables does not suit")
  X2           <- X1[,(1 + intercept):Kx1, drop = FALSE]
  
  # controls
  nR           <- object$smm.ctr$R
  nS           <- 1
  nT           <- 1
  
  #sizes
  N            <- sapply(dnetwork, nrow)
  M            <- length(N)
  if(M < 2) stop("Inference is not possible with one group")
  Nsum         <- sum(N)
  Ncum         <- c(0, cumsum(N))
  Ilist        <- lapply(N, diag)
  Pm           <- iv.power - 1
  ninstr       <- Kx1 + iv.power*Kx2
  
  # Parameters
  alpha        <- object$estimates[1]
  beta         <- object$estimates[-1]
  
  # Functions fmvzeta and fmvzetaH
  fmvzeta      <- NULL
  fmvzetaH     <- NULL
  Afmvzeta     <- list(alpha = alpha, beta = beta, R = nR, distr = dnetwork, y = y, W = W, smoother = FALSE, 
                       hN = 0, Kx1 = Kx1, ninstr = ninstr, M = M, N = N, Pm = Pm, Ncum = Ncum)
  # Type = 3 // GX and Gy are unobserved
  if(!GXobs & !Gyobs){
    if(fixed.effects){
      fmvzeta   <- PartialNetwork:::fmvzeta3fe
      fmvzetaH  <- PartialNetwork:::fmvzetaH3fe
      Afmvzeta  <- c(Afmvzeta, list(S = nS, T = nT, Ilist = Ilist, X1 = X1, X2 = X2, Kx2 = Kx2))
    } else{
      fmvzeta   <- PartialNetwork:::fmvzeta3
      fmvzetaH  <- PartialNetwork:::fmvzetaH3
      Afmvzeta  <- c(Afmvzeta, list(S = nS, T = nT, Ilist = Ilist, X1 = X1, X2 = X2, Kx2 = Kx2))
    }
  }
  
  # Estimate H0 
  H0        <- NULL
  if(cond.var){
    tderM   <- t(derM)
    H0      <- -solve(tderM %*% W %*% derM, tderM)
  } else{
    on.exit({
      if (is.null(sysSeed)) {
        rm(".Random.seed", envir = .GlobalEnv)
      } else {
        .GlobalEnv$.Random.seed <- sysSeed 
      }
    })
    assign(".Random.seed", seed, envir = .GlobalEnv)#Note that I restore the system seed (it it exist) using on.exit, see the beginning of the function
    tmp     <- do.call(fmvzetaH, Afmvzeta)
    derM    <- tmp$derM
    tderM   <- t(derM)
    aveMM   <- tmp$sumMM/Nsum
    aveM    <- tmp$sumM/Nsum
    H0      <- -solve(tderM %*% W %*% derM, tderM)
  }
  
  # Simulate En and Vn
  tmp   <- lapply(1:sim, function(iteration) fOMEGA(iteration, .fun, .args, fmvzeta, Afmvzeta, M))
  psi   <-   H0 %*% W %*% sapply(1:sim, function(x)  t(chol(tmp[[x]]$VZ/Nsum)) %*% rnorm(ninstr) + tmp[[x]]$EZ/sqrt(Nsum))
  Sn    <-  apply(psi, 2, function(x) c(alpha, beta) - x/sqrt(Nsum))
  Sn
}


fOMEGA <- function(iteration, .fun, .args, fmvzeta, Afmvzeta, M) {
  cat("Iteration:", iteration, sep = "", "\n")
  Afmvzeta$distr  <- do.call(.fun, .args)
  tmp   <- do.call(fmvzeta, Afmvzeta)
  sumMM <- tmp$sumMM
  sumM  <- tmp$sumM
  list(VZ = sumMM - sumM %*% t(sumM)/M, EZ = sumM)
}

*/
