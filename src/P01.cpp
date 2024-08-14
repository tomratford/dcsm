// [[Rcpp::depends(RcppArmadillo, RcppEigen, RcppNumerical)]]

#include <RcppArmadillo.h>
using namespace arma;

#include <RcppNumerical.h>
using namespace Numer;

#include "P01.h"

// Need this else we get some very annoying error.
// This function evaluates only a single point `x`
// It must also be provided as it is used to calculate infinite integrals.
double P01int::operator()(const double &x) const {
  vec l(1, fill::value(left));
  vec r(1, fill::value(right));
  vec z(1, fill::value(arm));
  vec v(1, fill::value(x));
  vec res = fns->P01Integrand(v, l, r, z);
  return (res(0));
}
// This function evaluates every point `x` all at once, which makes it faster.
void P01int::eval(double *x, const int n) const {
  vec l(n, fill::value(left));
  vec r(n, fill::value(right));
  vec z(n, fill::value(arm));
  vec v(n);
  std::copy(x, x + n, v.begin());
  vec res = fns->P01Integrand(v, l, r, z);
  std::copy(res.begin(), res.end(), x);
}

//' Generic-ish function to compute the 0 \eqn{\to} 1 transition probability
//'
//' @description
//' Typically wrapped by a separate function so that only the first 3 arguments are provided.
//' @param l The lower bounds of the (untransformed) interval
//' @param r The upper bounds of the (untransformed) interval
//' @param z The arms (case/control etc...)
//' @param type A string of either \code{"RoystonParmar"} or \code{"JolyPenalised"}
//' @param theta01_,theta02_,theta12_,gammas01,knots01,gammas02,knots02,gammas12,knots12,boundaries Passed to the \link{IllnessDeath} constructors to create the relevant set of functions.
//' @details
//' If \code{type = "RoystonParmar"} is used then a transformed interval with the substitution \eqn{u = exp(v)} is used to avoid a discontinuity at 0.
//' This then becomes an infinite integral.
// [[Rcpp::export]]
vec P01(vec &l, vec &r, vec &z, String type, double theta01_, double theta02_,
        double theta12_, vec gammas01, vec knots01, vec gammas02, vec knots02,
        vec gammas12, vec knots12, vec boundaries) {

  vec lower(l.n_elem);
  vec upper(l.n_elem);
  IllnessDeath* fns;
  if (type == String("RoystonParmar")) {
    fns = new RoystonParmar(theta01_, theta02_, theta12_, gammas01, knots01,
                      gammas02, knots02, gammas12, knots12, boundaries);
    l.clamp(1e-9,datum::inf); // Do not allow zero => -inf
    lower = log(l);
    upper = log(r);
  } else if (type ==  "JolyPenalised") {
    fns = new JolyPenalised(theta01_, theta02_, theta12_, gammas01, knots01,
                      gammas02, knots02, gammas12, knots12, boundaries);
    lower = l;
    upper = r;
  } else {
    stop("type must be one of `RoystonParmar` or `JolyPenalised`");
  }
  
  vec res(l.n_elem);
  for (size_t i = 0; i < l.n_elem; i++) {
    double err_est;
    int err_code;

    P01int f(l(i), r(i), z(i), fns);
    res(i) = integrate(f, lower(i), upper(i), err_est, err_code);
    if (err_code > 0) {
      switch (err_code) {
      case 2:
        std::cout << "WARNING: Roundoff error, result = " << res(i) << " error estimate = " << err_est << "\n";
        res(i) = -1e10;
        break;
      default:
        stop("Integration code greater than zero: " + std::to_string(err_code));
      }
    }
  }
  
  delete fns;
  return (res);
}
