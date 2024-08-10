#ifndef P01_H
#define P01_H

#include <RcppArmadillo.h>
using namespace arma;

#include <RcppNumerical.h>
using namespace Numer;

#include "illness_death.h"

class P01int : public Func {
 private:
  double left;
  double right;
  double arm;
  IllnessDeath* fns;

 public:
 P01int(double left_, double right_, double arm_, IllnessDeath *fns_)
   : left(left_), right(right_), arm(arm_), fns(fns_) {}
  double operator()(const double &x) const;
  void eval(double *x, const int n) const;
};

vec P01(const vec &l, const vec &r, const vec &z, String type, double theta01_, double theta02_,
        double theta12_, vec gammas01, vec knots01, vec gammas02, vec knots02,
        vec gammas12, vec knots12, vec boundaries);

#endif
