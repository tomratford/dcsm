#ifndef SPLINES_H
#define SPLINES_H

#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

class Spline {
public:
  virtual ~Spline() {}
  virtual vec S(const vec& t) const = 0;
  virtual vec dS(const vec& t) const = 0;
  virtual vec intensity(const vec& t, const vec& z, double theta) const = 0;
  virtual vec Intensity(const vec& t, const vec& z, double theta) const = 0;

  vec minus_Intensity(const vec& l, const vec& r, const vec& z, double theta) const {
    vec res = Intensity(l, z, theta) - Intensity(r, z, theta);
    return res;
  }
};

class NaturalCubicSpline : public Spline {
private:
  const vec gammas;
  const vec knots;
  const vec boundaries;

public:
  NaturalCubicSpline(const vec& gammas_, const vec& knots_, const vec& boundaries_)
      : gammas(gammas_), knots(knots_), boundaries(boundaries_) {}

  vec S(const vec& x) const override;
  vec dS(const vec& x) const override;
  vec intensity(const vec& x, const vec& z, double theta) const override;
  vec Intensity(const vec& x, const vec& z, double theta) const override;
};

class CubicISpline : public Spline {
private:
  const vec gammas;
  const vec knots;
  const vec boundaries;

public:
  CubicISpline(const vec& gammas_, const vec& knots_, const vec& boundaries_)
      : gammas(gammas_), knots(knots_), boundaries(boundaries_) {}

  vec S(const vec& t) const override;
  vec dS(const vec& t) const override;
  vec intensity(const vec& t, const vec& z, double theta) const override;
  vec Intensity(const vec& t, const vec& z, double theta) const override;
};

#endif
