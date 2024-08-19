#ifndef ILLNESS_DEATH_H
#define ILLNESS_DEATH_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

#include "splines.h"

class IllnessDeath {
 public:
  virtual ~IllnessDeath() {}
  
  const double theta01, theta02, theta12;  
  std::unique_ptr<Spline> spline01, spline02, spline12;
  SEXP P01;

  IllnessDeath(double theta01_, double theta02_, double theta12_,
               std::unique_ptr<Spline> spline01_, std::unique_ptr<Spline> spline02_, std::unique_ptr<Spline> spline12_)
      : theta01(theta01_), theta02(theta02_), theta12(theta12_),
    spline01(std::move(spline01_)), spline02(std::move(spline02_)), spline12(std::move(spline12_)) {}

  vec int01(const vec &t, const vec &z) const {
    return (spline01->intensity(t, z, theta01));
  }
  vec int02(const vec &t, const vec &z) const {
    return (spline02->intensity(t, z, theta02));
  }
  vec int12(const vec &t, const vec &z) const {
    return (spline12->intensity(t, z, theta02));
  }
  vec logP00(const vec &l, const vec &r, const vec &z) const {
    return (spline01->minus_Intensity(l,r,z,theta01) +
            spline02->minus_Intensity(l,r,z,theta02));
  }
  vec P00(const vec &l, const vec &r, const vec &z) const {
    return (exp(logP00(l,r,z)));
  }
  vec logP11(const vec &l, const vec &r, const vec &z) const {
    return (spline12->minus_Intensity(l,r,z,theta12));
  }
  vec P11(const vec &l, const vec &r, const vec &z) const {
    return (exp(logP11(l,r,z)));
  }
  virtual vec P01Integrand(const vec &u, const vec &l, const vec &r, const vec &z) const = 0;
};

class RoystonParmar : public IllnessDeath {
public:
  RoystonParmar(double theta01_, double theta02_, double theta12_,
               const vec& gammas01, const vec& knots01,
               const vec& gammas02, const vec& knots02,
               const vec& gammas12, const vec& knots12, const vec& boundaries)
      : IllnessDeath(theta01_, theta02_, theta12_,
                     std::make_unique<NaturalCubicSpline>(gammas01, knots01, boundaries),
                     std::make_unique<NaturalCubicSpline>(gammas02, knots02, boundaries),
                     std::make_unique<NaturalCubicSpline>(gammas12, knots12, boundaries)) {}

  vec P01Integrand(const vec& v, const vec& l, const vec& r, const vec& z) const override;
};

class JolyPenalised : public IllnessDeath {
public:
  JolyPenalised(double theta01_, double theta02_, double theta12_,
               const vec& gammas01, const vec& knots01,
               const vec& gammas02, const vec& knots02, 
               const vec& gammas12, const vec& knots12, const vec& boundaries)
      : IllnessDeath(theta01_, theta02_, theta12_,
                     std::make_unique<CubicISpline>(gammas01, knots01, boundaries),
                     std::make_unique<CubicISpline>(gammas02, knots02, boundaries),
                     std::make_unique<CubicISpline>(gammas12, knots12, boundaries)) {}

  vec P01Integrand(const vec& u, const vec& l, const vec& r, const vec& z) const override;
};

IllnessDeath* MakeIllnessDeath(String type, double theta01_, double theta02_, double theta12_,
                               const vec& gammas01, const vec& knots01,
                               const vec& gammas02, const vec& knots02,
                               const vec& gammas12, const vec& knots12, const vec& boundaries);

#endif
