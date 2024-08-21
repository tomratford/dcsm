// [[Rcpp::depends(RcppArmadillo, splines2)]]

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

#include "illness_death.h"

vec RoystonParmar::P01Integrand(const vec& v, const vec& l, const vec& r, const vec& z) const {
  vec u = exp(v);
  vec res = P00(l, u, z) % spline01->dS(v) % exp(spline01->S(v)) % exp(z*theta01) % P11(u, r, z);
  return (res);
};

vec JolyPenalised::P01Integrand(const vec& u, const vec& l, const vec& r, const vec& z) const {
  return (P00(l, u, z) % spline01->intensity(u,z,theta01) % P11(u, r, z));
};

IllnessDeath* MakeIllnessDeath(String type, double theta01_, double theta02_, double theta12_,
                               const vec& gammas01, const vec& knots01,
                               const vec& gammas02, const vec& knots02,
                               const vec& gammas12, const vec& knots12, const vec& boundaries) {
  if (type == String("RoystonParmar")) {
    return new RoystonParmar(theta01_, theta02_, theta12_, gammas01, knots01,
                      gammas02, knots02, gammas12, knots12, boundaries);
  } else if (type ==  "JolyPenalised") {
    return new JolyPenalised(theta01_, theta02_, theta12_, gammas01, knots01,
                      gammas02, knots02, gammas12, knots12, boundaries);
  } else {
    stop("type must be one of `RoystonParmar` or `JolyPenalised`");
  }
}


RCPP_MODULE(illness_death_mod) {
  class_<IllnessDeath>("IllnessDeath")
    .factory<String, double, double, double, const vec&, const vec&, const vec&, const vec&, const vec&, const vec&, const vec&>(MakeIllnessDeath)
    .method("int01", &IllnessDeath::int01)
    .method("int02", &IllnessDeath::int02)
    .method("int12", &IllnessDeath::int12)
    .method("Int01", &IllnessDeath::Int01)
    .method("Int02", &IllnessDeath::Int02)
    .method("Int12", &IllnessDeath::Int12)
    .method("logP00", &IllnessDeath::logP00)
    .method("P00", &IllnessDeath::P00)
    .method("logP11", &IllnessDeath::logP11)
    .method("P11", &IllnessDeath::P11)
    .method("P01Integrand", &IllnessDeath::P01Integrand)
    .field("P01", &IllnessDeath::P01)
  ;
}
