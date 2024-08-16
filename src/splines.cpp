// [[Rcpp::depends(RcppArmadillo, splines2)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

#include <splines2Armadillo.h>
using namespace splines2;

#include "splines.h"

// Natural cubic splines definitions

// x in these function refers to log(t).

vec NaturalCubicSpline::S(const vec &x) const {
  NaturalSpline spline(x, knots, boundaries);
  vec product = spline.basis() * gammas;
  product.replace(datum::nan, -datum::inf); // Working under the assumption that NaN => x = -inf => t = 0
  return (product);
}

vec NaturalCubicSpline::dS(const vec &x) const {
  NaturalSpline spline(x, knots, boundaries);
  vec product = spline.derivative() * gammas;
  product.replace(datum::nan, -datum::inf); // Working under the assumption that NaN => x = -inf => t = 0
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

vec CubicISpline::dS2(const vec& t) const {
  ISpline spline(t, knots, 3, boundaries);
  vec product = spline.derivative(2) * gammas;
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

RCPP_MODULE(splines_mod) {
  class_<NaturalCubicSpline>("NaturalCubicSpline")
    .constructor<vec, vec, vec>()
    .method("S", &NaturalCubicSpline::S)
    .method("dS", &NaturalCubicSpline::dS)
    .method("int", &NaturalCubicSpline::intensity)
    .method("Int", &NaturalCubicSpline::Intensity)
    ;
  class_<CubicISpline>("CubicISpline")
    .constructor<vec, vec, vec>()
    .method("S", &CubicISpline::S)
    .method("dS", &CubicISpline::dS)
    .method("dS2", &CubicISpline::dS2)
    .method("int", &CubicISpline::intensity)
    .method("Int", &CubicISpline::Intensity)
    ;
}

