// [[Rcpp::depends(RcppArmadillo, RcppEigen, RcppNumerical, splines2)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

#include <RcppNumerical.h>
using namespace Numer;

#include <splines2Armadillo.h>
using namespace splines2;

#include <functional>
#include "royston_parmar.h"

// [[Rcpp::plugins(cpp11)]]

NaturalCubicSpline::NaturalCubicSpline(vec &gammas_, vec &knots_,
                                       vec &boundaries_)
  : gammas(gammas_), knots(knots_), boundaries(boundaries_) {}

vec NaturalCubicSpline::S(vec &x) const {
  NaturalSpline spline(x, knots, boundaries);
  vec product = spline.basis() * gammas;
  return (product);
}

vec NaturalCubicSpline::dS(vec &x) const {
  NaturalSpline spline(x, knots, boundaries);
  vec product = spline.derivative() * gammas;
  return (product);
}

vec NaturalCubicSpline::intensity(vec &x, vec &z, double theta) {
  vec res = 1 / exp(x) * dS(x) * exp(S(x) + z * theta);
  return (res);
}

vec minus_Intensity(vec &l, vec &r, vec &z, double theta,
                    NaturalCubicSpline spline) {
  // Minus the integrated intensity (hence capitalised I)
  vec res = exp(spline.S(l)) - exp(spline.S(r));
  return (res * exp(z * theta));
}

RoystonParmarFns::RoystonParmarFns(double theta01_, double theta02_,
                                   double theta12_, vec gammas01, vec knots01,
                                   vec gammas02, vec knots02, vec gammas12,
                                   vec knots12, vec boundaries)
  : theta01(theta01_), theta02(theta02_), theta12(theta12_),
    spline01(gammas01, knots01, boundaries),
    spline02(gammas02, knots02, boundaries),
    spline12(gammas12, knots12, boundaries) {}

vec RoystonParmarFns::int02(vec &x, vec &z) {
  return (spline02.intensity(x, z, theta02));
}

vec RoystonParmarFns::int12(vec &x, vec &z) {
  return (spline12.intensity(x, z, theta12));
}

vec RoystonParmarFns::logP00(vec &l, vec &r, vec &z) const {
  return (minus_Intensity(l, r, z, theta01, spline01) +
          minus_Intensity(l, r, z, theta02, spline02));
}

vec RoystonParmarFns::P00(vec &r, vec &l, vec &z) const {
  return (exp(logP00(l, r, z)));
}

vec RoystonParmarFns::logP11(vec &l, vec &r, vec &z) const {
  return (minus_Intensity(l, r, z, theta12, spline12));
}

vec RoystonParmarFns::P11(vec &r, vec &l, vec &z) const {
  return (exp(logP11(l, r, z)));
}

P01Integrand::P01Integrand(double leftbound_, double rightbound_, vec &z_,
                           RoystonParmarFns fns_)
  : leftbound(leftbound_), rightbound(rightbound_), z(z_), fns(fns_) {}
void P01Integrand::eval(double *x, const int n) const {
  vec l(1);
  vec r(1);
  l(0) = leftbound;
  r(0) = rightbound;
  vec v(n);
  std::copy(x, x + n, v.begin());
  vec u = exp(v);
  vec res = fns.P00(l, u, z) * fns.spline01.dS(v) * exp(fns.spline01.S(v)) *
    fns.P11(u, r, z);
  x = res.memptr();
}

vec RoystonParmarFns::P01(vec &l, vec &r, vec &z) {
  vec res(l.n_elem);
  for (size_t i = 0; i < l.n_elem; i++) {
    double lower = l(i);
    double upper = r(i);
    double err_est;
    int err_code;

    P01Integrand f(lower, upper, z, *this);
    res(i) = integrate(f, lower, upper, err_est, err_code);
  }
  return (res);
}

RCPP_MODULE(mod_roystonparmar) {
  class_<RoystonParmarFns>("RoystonParmarFns")
  .constructor<double, double, double, vec, vec, vec, vec, vec, vec, vec>()
  .method("int02", &RoystonParmarFns::int02)
  .method("int12", &RoystonParmarFns::int12)
  .method("logP00", &RoystonParmarFns::logP00)
  .method("P00", &RoystonParmarFns::P00)
  .method("logP11", &RoystonParmarFns::logP11)
  .method("P11", &RoystonParmarFns::P11)
  .method("P01", &RoystonParmarFns::P01)
  ;
}
