// [[Rcpp::depends(RcppArmadillo, RcppEigen, RcppNumerical)]]

#include <RcppArmadillo.h>
using namespace arma;

#include <RcppNumerical.h>
using namespace Numer;

#include "P01.h"

// Need this else we get some very annoying error
double P01int::operator()(const double &x) const { return (0); }

// [[Rcpp::export]]
vec P01(vec &l, vec &r, vec &z, String type, double theta01_, double theta02_,
        double theta12_, vec gammas01, vec knots01, vec gammas02, vec knots02,
        vec gammas12, vec knots12, vec boundaries) {

  IllnessDeath* fns;
  if (type == String("RoystonParmar")) {
    fns = new RoystonParmar(theta01_, theta02_, theta12_, gammas01, knots01,
                      gammas02, knots02, gammas12, knots12, boundaries);
  } else if (type ==  "JolyPenalised") {
    fns = new JolyPenalised(theta01_, theta02_, theta12_, gammas01, knots01,
                      gammas02, knots02, gammas12, knots12, boundaries);
  } else {
    stop("type must be one of `RoystonParmar` or `JolyPenalised`");
  }
  
  vec res(l.n_elem);
  vec logl = log(l);
  vec logr = log(r);
  for (size_t i = 0; i < l.n_elem; i++) {
    double err_est;
    int err_code;

    P01int f(l(i), r(i), z(i), fns);
    res(i) = integrate(f, logl(i), logr(i), err_est, err_code);
    if (err_code > 0) {
      exception e("Error code greater than zero", err_code);
      throw e;
      stop("This should have stopped by now");
    }
  }
  
  delete fns;
  return (res);
}
