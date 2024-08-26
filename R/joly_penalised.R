#' Fit a dataset to a model with I-spline baseline transition intensities.
#' @description
#' Requires known interior knots and penalty parameters.
#'
#' @param data Data set to model from
#' @param k01 Number of interior knots in the 0 -> 1 transition
#' @param k02 Number of interior knots in the 0 -> 2 transition
#' @param k12 Number of interior knots in the 1 -> 2 transition
#' @param kappa01 Penalty parameter for the 0 -> 1 transition
#' @param kappa02 Penalty parameter for the 0 -> 1 transition
#' @param kappa12 Penalty parameter for the 0 -> 1 transition
#' @param compute_cross Compute the cross value (and relevant hessians)? EXPENSIVE!
#' @param debug Print every imputation?
#' @param initials optional list of initial values to replace call to \link{joly.initials}
#'
#' @return A \code{optim} object.
#' @export
joly.fit <- function(data,
                     k01,
                     k02,
                     k12,
                     kappa01 = 1,
                     kappa02 = 1,
                     kappa12 = 1,
                     compute_cross = T,
                     debug = F,
                     initials = NULL,
                     control = list(fnscale = -1, maxit = 500),
                     method = "BFGS",
                     ...) {
  if (is.null(initials))
    initials <- joly.initials(data, k01, k02, k12)
  p_init <- unlist(initials[c("theta01",
                              "theta02",
                              "theta12",
                              "gammas01",
                              "gammas02",
                              "gammas12")])
  opt_out <- optim(p_init, \(p) {
    pars <- make_pars2(p, initials)
    if (debug)
      print(pars[c("theta01",
                   "theta02",
                   "theta12",
                   "gammas01",
                   "gammas02",
                   "gammas12")])

    spline01 <- new(CubicISpline,
                    pars$gammas01,
                    pars$knots01,
                    pars$boundaries)
    spline02 <- new(CubicISpline,
                    pars$gammas02,
                    pars$knots02,
                    pars$boundaries)
    spline12 <- new(CubicISpline,
                    pars$gammas12,
                    pars$knots12,
                    pars$boundaries)

    # Get penalties
    penalty01 <- integrate(
      \(x) spline01$dS2(x) ^ 2,
      lower = pars$boundaries[1],
      upper = pars$boundaries[2]
    )$value
    penalty02 <- integrate(
      \(x) spline02$dS2(x) ^ 2,
      lower = pars$boundaries[1],
      upper = pars$boundaries[2]
    )$value
    penalty12 <- integrate(
      \(x) spline12$dS2(x) ^ 2,
      lower = pars$boundaries[1],
      upper = pars$boundaries[2]
    )$value
    penalty <- kappa01 * penalty01 + kappa02 * penalty02 + kappa12 * penalty12

    if (debug)
      print(penalty)
    res <- tryCatch(
      joly.ll(pars, data),
      error = \(e) {
        print(e)
        return(-1e10)
      }
    )
    if (debug)
      print(res)
    res - penalty
  }, control = control, method = method, ...)
  opt_out$loglik <- joly.ll(make_pars2(opt_out$par, initials), data)
  opt_out$AIC <- 2 * ((15 + k01 + k02 + k12) - opt_out$loglik)
  if (compute_cross) {
    opt_out$hessian2 <- optimHess(opt_out$par, \(p) {
      pars <- make_pars2(p, initials)
      res <- joly.ll(pars, data)
      res
    }, control = list(fnscale = -1))
    opt_out$hessian <- opt_out$hessian2 - optimHess(opt_out$par, \(p) {
      pars <- make_pars2(p, initials)
      spline01 <- new(CubicISpline,
                      pars$gammas01,
                      pars$knots01,
                      pars$boundaries)
      spline02 <- new(CubicISpline,
                      pars$gammas02,
                      pars$knots02,
                      pars$boundaries)
      spline12 <- new(CubicISpline,
                      pars$gammas12,
                      pars$knots12,
                      pars$boundaries)

      # Get penalties
      penalty01 <- integrate(
        \(x) spline01$dS2(x) ^ 2,
        lower = pars$boundaries[1],
        upper = pars$boundaries[2]
      )$value
      penalty02 <- integrate(
        \(x) spline02$dS2(x) ^ 2,
        lower = pars$boundaries[1],
        upper = pars$boundaries[2]
      )$value
      penalty12 <- integrate(
        \(x) spline12$dS2(x) ^ 2,
        lower = pars$boundaries[1],
        upper = pars$boundaries[2]
      )$value
      kappa01 * penalty01 + kappa02 * penalty02 + kappa12 * penalty12
    }, control = list(fnscale = -1))
    opt_out$cross_value <- opt_out$loglik - sum(diag(solve(opt_out$hessian) %*%
                                                       opt_out$hessian2))
  }
  attr(opt_out, "dist") <- "joly"
  attr(opt_out, "initials") <- initials
  class(opt_out) <- c("dcsm_mod")
  opt_out
}

#' Get a set of feasible initial values for a I-spline based model
#'
#' @param data Data set to model from, with columns `L`, `R`, `V`, `Delta1`, `Delta2`
#' @param k01 Number of interior knots in the 0 -> 1 transition
#' @param k02 Number of interior knots in the 0 -> 2 transition
#' @param k12 Number of interior knots in the 1 -> 2 transition
#'
#' @return A list of parameters, to be fed to \link{joly.ll}
#' @export
#'
#' @import dplyr
#' @importFrom survival survfit Surv coxph
joly.initials <- function(data, k01, k02, k12) {
  # set initial list
  initials = list(
    theta01 = 0,
    theta02 = 0,
    theta12 = 0,
    gammas01 = numeric(0),
    knots01 = numeric(0),
    gammas02 = numeric(0),
    knots02 = numeric(0),
    gammas12 = numeric(0),
    knots12 = numeric(0),
    #either day 1 or the smallest time in the data
    boundaries = c(0, max(data$V))
  )
  # work out knot values
  delta1 <-  pull(data, matches("d(elta)?0?1", perl=T))
  delta2 <-  pull(data, matches("d(elta)?0?2", perl=T))
  Rs <- pull(data, matches("^r(ight)?", perl=T))
  Vs <- pull(data, matches("^v$"))
  #important times
  itm <- coalesce(
    if_else(as.logical(delta1), Rs, NA),
    if_else(as.logical(delta2), Vs, NA)
  ) |> na.omit() |> sort()
  itmn <- length(itm)
  if (k01 > 0) {
    index <- round(seq(0,1,by=1/(k01+1))*itmn)
    initials$knots01 <- c(itm[index[-length(index)]]) #nb: index[length(index)] = boundary knot
  }
  if (k02 > 0) {
    index <- round(seq(0,1,by=1/(k02+1))*itmn)
    initials$knots02 <- c(itm[index[-length(index)]]) #nb: index[length(index)] = boundary knot
  }
  if (k12 > 0) {
    index <- round(seq(0,1,by=1/(k12+1))*itmn)
    initials$knots12 <- c(itm[index[-length(index)]]) #nb: index[length(index)] = boundary knot
  }

  # Set each gamma as max(cumhaz)/4
  cumhaz01 <- max(survfit(Surv(R, delta1) ~ 1, data)$cumhaz)
  cumhaz02 <- max(survfit(Surv(V, delta2) ~ 1, filter(data, delta1 == 0))$cumhaz)
  cumhaz12 <- max(survfit(Surv(V, delta2) ~ 1, filter(data, delta1 == 1))$cumhaz)

  initials$theta01 <- unname(coxph(Surv(R,delta1) ~ ATRTN, data)$coef)
  initials$theta02 <- unname(coxph(Surv(V,delta2) ~ ATRTN, filter(data, delta1 == 0))$coef)
  initials$theta12 <- unname(coxph(Surv(V,delta2) ~ ATRTN, filter(data, delta1 == 1))$coef)

  initials$gammas01 <- rep((cumhaz01/(k01+4)), k01+4)
  initials$gammas02 <- rep((cumhaz02/(k02+4)), k02+4)
  initials$gammas12 <- rep((cumhaz12/(k12+4)), k12+4)

  return(initials)
}

#' Compute the log likelihood for a model with illness death and I-spline cumulative intensities
#'
#' @param p parameters, as a named list based off \link{Joly}
#' @param data dataset as per instructions in \link{dc_loglik}
#' @param indiv return a vector of each individual log-likelihood instead of the sum?
#'
#' @return The log-likelihood for the data and parameters provided.
#' @export
joly.ll <- function(p, data, indiv=FALSE) {
  lik <- dc_loglik(data,
                   fns = do.call(
                     joly.fnBuilder,
                     p
                   ))
  lik[lik == -Inf] <- -1e10 # avoid non-finite diff
  if (indiv) {
    return(lik)
  } else {
    return(sum(lik))
  }
}

#' Functions for a illness-death model with I-spline cumulative intensities
#' @description
#' Based on research by Joly et al. (1998, 2002)
#'
#' @param theta01 A single numeric value representing treatment difference in the progression intensity
#' @param theta02 A single numeric value representing treatment difference in the pre-progression death intensity
#' @param theta12 A single numeric value representing treatment difference in the post-progression death- intensity
#' @param gammas01 Gamma values for the progression intensity spline
#' @param knots01 Interior knots for the progression intensity spline
#' @param gammas02 Gamma values for the pre-progression death intensity spline
#' @param knots02 Interior knots for the pre-progression death intensity spline
#' @param gammas12 Gamma values for the post-progression death intensity spline
#' @param knots12 Interior knots for the post-progression death intensity spline
#' @param boundaries Boundary knots for the data (typically day 1 and the administrative censoring time). Cannot be 0.
#'
#' @details
#' Each transition's cumulative intensity is modeled as
#' \deqn{
#'  \text{A}_{hj} = S_{hj}(t;\boldsymbol{\gamma})\exp(z\theta_{hj})
#' }
#' \code{gammasXX} must contain the \eqn{\boldsymbol{\gamma}_{hj}} values as parameters.
#' Using an example, if we have one knotpoint in the transition from 0 to 1: \code{gamma010}, \code{gamma011}, \code{gamma012} must be provided as parameters.
#'
#' @return A list of functions.
#'
#' @importFrom methods new
#'
#' @rdname Joly
joly.fnBuilder <- function(theta01, theta02, theta12, gammas01, knots01, gammas02, knots02, gammas12, knots12, boundaries) {

  fns <- new(IllnessDeath,
             "JolyPenalised",
             theta01, theta02, theta12,
             gammas01, knots01,
             gammas02, knots02,
             gammas12, knots12,
             boundaries)
  fns$P01 <- \(l,r,z) P01(l,r,z,"JolyPenalised",
                          theta01, theta02, theta12,
                          gammas01, knots01,
                          gammas02, knots02,
                          gammas12, knots12,
                          boundaries)
  # fns$P01 <- Vectorize(\(l, r, z) tryCatch(
  #   integrate(
  #     \(u) res <- fns$P01Integrand(u, rep(l, length(u)), rep(r, length(u)), rep(z, length(u))),
  #     lower = l,
  #     upper = r
  #   )$value,
  #   error = \(e) {
  #     print(e$message)
  #     curve(fns$P01Integrand(x, rep(l, length(x)), rep(r, length(x)), rep(z, length(x))),
  #           from = l,
  #           to = r)
  #     stop()
  #   }
  # ))
  fns
}
