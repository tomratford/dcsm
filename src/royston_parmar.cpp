// [[Rcpp::depends(RcppArmadillo, RcppEigen, RcppNumerical, splines2)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

#include <RcppNumerical.h>
using namespace Numer;

#include <splines2Armadillo.h>
using namespace splines2;

#include "royston_parmar.h"

// [[Rcpp::plugins(cpp11)]]

NaturalCubicSpline::NaturalCubicSpline(vec gammas_, vec knots_, vec boundaries_)
    : gammas(gammas_), knots(knots_), boundaries(boundaries_) {}

vec NaturalCubicSpline::S(vec &x) const {
  NaturalSpline spline(x, knots, boundaries);
  vec product = spline.basis() * gammas;
  return (product);
}

vec NaturalCubicSpline::dS(vec &x) const {
  NaturalSpline spline(x, knots, boundaries);
  // We do not want to include the intercept, force out using gammas
  vec gammas_(gammas);
  gammas_(0) = 0;
  vec product = spline.derivative() * gammas_;
  return (product);
}

vec NaturalCubicSpline::intensity(vec &t, vec &z, double theta) {
  vec x = log(t);
  vec res = (1 / exp(x)) % dS(x) % exp(S(x) + z * theta);
  return (res);
}

NaturalCubicSpline *makeSpline(vec gammas, vec knots, vec boundaries) {
  return new NaturalCubicSpline(gammas, knots, boundaries);
}

RCPP_MODULE(mySplines) {
  class_<NaturalCubicSpline>("NaturalCubicSpline")
      .factory<vec, vec, vec>(makeSpline)
      .method("S", &NaturalCubicSpline::S)
      .method("dS", &NaturalCubicSpline::dS)
      .method("intensity", &NaturalCubicSpline::intensity);
}

// minus the integrated intensity (hence capitalised I) between points l and r
vec minus_Intensity(vec &l, vec &r, vec &z, double theta,
                    NaturalCubicSpline spline) {
  vec res = exp(spline.S(l)) - exp(spline.S(r));
  return (res % exp(z * theta));
}

RoystonParmarFns::RoystonParmarFns(double theta01_, double theta02_,
                                   double theta12_, vec gammas01, vec knots01,
                                   vec gammas02, vec knots02, vec gammas12,
                                   vec knots12, vec boundaries)
    : theta01(theta01_), theta02(theta02_), theta12(theta12_),
      spline01(gammas01, knots01, boundaries),
      spline02(gammas02, knots02, boundaries),
      spline12(gammas12, knots12, boundaries) {}

vec RoystonParmarFns::int02(vec &t, vec &z) {
  return (spline02.intensity(t, z, theta02));
}

vec RoystonParmarFns::int12(vec &t, vec &z) {
  return (spline12.intensity(t, z, theta12));
}

vec RoystonParmarFns::logP00(vec &l, vec &r, vec &z) const {
  return (minus_Intensity(l, r, z, theta01, spline01) +
          minus_Intensity(l, r, z, theta02, spline02));
}

vec RoystonParmarFns::P00(vec &l, vec &r, vec &z) const {
  return (exp(logP00(l, r, z)));
}

vec RoystonParmarFns::logP11(vec &l, vec &r, vec &z) const {
  return (minus_Intensity(l, r, z, theta12, spline12));
}

vec RoystonParmarFns::P11(vec &l, vec &r, vec &z) const {
  return (exp(logP11(l, r, z)));
}

//' @name RoystonParmarFns$P01Integrand
//' @title Compute the integrand for the transition probability of progression.
//' @param v A vector of the log times
//' @param l The lower bound of the integral (not logged)
//' @param r The upper bound of the integral (not logged)
//' @param z The vector of patient treatment arms
vec RoystonParmarFns::P01Integrand(vec &v, vec &l, vec &r, vec &z) const {
  vec u = exp(v);
  return (P00(l, u, z) % spline01.dS(v) % exp(spline01.S(v)) % P11(u, r, z));
}

P01int::P01int(double left_, double right_, double arm_, RoystonParmarFns &fns_)
    : left(left_), right(right_), arm(arm_), fns(fns_) {}
double P01int::operator()(const double &x) const { return (0); }
void P01int::eval(double *x, const int n) const {
  vec l(n);
  l.fill(left);
  vec r(n);
  r.fill(right);
  vec z(n);
  z.fill(arm);
  vec v(n);
  std::copy(x, x + n, v.begin());
  vec res = fns.P01Integrand(v, l, r, z);
  std::copy(res.begin(), res.end(), x);
}

// [[Rcpp::export]]
vec P01(vec &l, vec &r, vec &z, double theta01_, double theta02_,
        double theta12_, vec gammas01, vec knots01, vec gammas02, vec knots02,
        vec gammas12, vec knots12, vec boundaries) {
  RoystonParmarFns fns(theta01_, theta02_, theta12_, gammas01, knots01,
                       gammas02, knots02, gammas12, knots12, boundaries);
  vec res(l.n_elem);
  vec logl = log(l);
  vec logr = log(r);
  for (size_t i = 0; i < l.n_elem; i++) {
    double err_est;
    int err_code;

    P01int f(l(i), r(i), z(i), fns);
    res(i) = integrate(f, logl(i), logr(i), err_est, err_code);
  }
  return (res);
}

RCPP_MODULE(RoystonParmar) {
  class_<RoystonParmarFns>("RoystonParmarFns")
      .constructor<double, double, double, vec, vec, vec, vec, vec, vec, vec>()
      .method("int02", &RoystonParmarFns::int02)
      .method("int12", &RoystonParmarFns::int12)
      .method("logP00", &RoystonParmarFns::logP00)
      .method("P00", &RoystonParmarFns::P00)
      .method("logP11", &RoystonParmarFns::logP11)
      .method("P11", &RoystonParmarFns::P11)
      .method("P01Integrand", &RoystonParmarFns::P01Integrand)
      .field("P01", &RoystonParmarFns::P01)
  ;
}
