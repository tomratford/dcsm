#ifndef ROYSTON_PARMAR_H
#define ROYSTON_PARMAR_H

#include <RcppArmadillo.h>
using namespace arma;
#include <RcppNumerical.h>
using namespace Numer;

class NaturalCubicSpline {
private:
  const vec gammas;
  const vec knots;
  const vec boundaries;

public:
  NaturalCubicSpline(vec gammas_, vec knots_, vec boundaries_);
  vec S(vec &x) const;
  vec dS(vec &x) const;
  vec intensity(vec &x, vec &z, double theta);
  void debug();
};

vec minus_Intensity(vec &l, vec &r, vec &z, double theta,
                    NaturalCubicSpline spline);

class RoystonParmarFns {
private:
  const double theta01;
  const double theta02;
  const double theta12;

public:
  NaturalCubicSpline spline01;
  NaturalCubicSpline spline02;
  NaturalCubicSpline spline12;
  SEXP P01;

  RoystonParmarFns(double theta01_, double theta02_, double theta12_,
                   vec gammas01, vec knots01, vec gammas02, vec knots02,
                   vec gammas12, vec knots12, vec boundaries);

  vec int02(vec &x, vec &z);
  vec int12(vec &x, vec &z);
  vec logP00(vec &l, vec &r, vec &z) const;
  vec P00(vec &r, vec &l, vec &z) const;
  vec logP11(vec &l, vec &r, vec &z) const;
  vec P11(vec &r, vec &l, vec &z) const;
  vec P01Integrand(vec &v, vec &l, vec &r, vec &z) const;
};

class P01int : public Func {
private:
  double left;
  double right;
  double arm;
  RoystonParmarFns &fns;

public:
  P01int(double left_, double right_, double arm_, RoystonParmarFns &fns_);
  double operator()(const double &x) const;
  void eval(double *x, const int n) const;
};

vec P01(vec &l, vec &r, vec &z, RoystonParmarFns &fns_);

#endif
