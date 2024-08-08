// [[Rcpp::depends(RcppArmadillo, splines2)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

#include <splines2Armadillo.h>
using namespace splines2;

#include "splines.h"

// Natural cubic splines definitions

vec NaturalCubicSpline::S(const vec &x) const {
  NaturalSpline spline(x, knots, boundaries);
  vec product = spline.basis() * gammas;
  return (product);
}

vec NaturalCubicSpline::dS(const vec &x) const {
  NaturalSpline spline(x, knots, boundaries);
  // We do not want to include the intercept, force out using gammas
  vec gammas_(gammas);
  gammas_(0) = 0;
  vec product = spline.derivative() * gammas_;
  return (product);
}

vec NaturalCubicSpline::intensity(const vec &t, const vec &z, double theta) const {
  vec x = log(t);
  vec res = (1 / t) % dS(x) % exp(S(x) + z * theta);
  return (res);
}

vec NaturalCubicSpline::Intensity(const vec &t, const vec &z, double theta) const {
  vec x = log(t);
  return (exp(S(x) + z * theta));
}

// I-spline definitions

vec CubicISpline::S(const vec& t) const {
  ISpline spline(t, knots, 3, boundaries);
  vec product = spline.basis() * gammas;
  return (product);
}

vec CubicISpline::dS(const vec& t) const {
  ISpline spline(t, knots, 3, boundaries);
  vec product = spline.derivative() * gammas;
  return (product);
}

vec CubicISpline::intensity(const vec& t, const vec& z, double theta) const {
  vec res = dS(t) % exp(z * theta);
  return (res);
}

vec CubicISpline::Intensity(const vec& t, const vec& z, double theta) const {
  vec res = S(t) % exp(z * theta);
  return (res);
}


