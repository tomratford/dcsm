#' Get a set of (hopefully) feasible inital values for a Natural Cubic Spline based model
#'
#' @param data Data set to model from, with columns `L`, `R`, `V`, `Delta1`, `Delta2`
#' @param k01 Number of interior knots in the 0 -> 1 transition
#' @param k02 Number of interior knots in the 0 -> 2 transition
#' @param k12 Number of interior knots in the 1 -> 2 transition
#'
#' @return A list of parameters, to be fed to \link{royston_parmar.ll}
#' @export
#'
#' @importFrom dplyr pull matches coalesce if_else
#' @importFrom survival survfit Surv
royston_parmar.initials <- function(data, k01, k02, k12) {
  # set initial list
  initials = list(
    theta01 = numeric(0),
    theta02 = numeric(0),
    theta12 = numeric(0),
    gammas01 = numeric(0),
    knots01 = numeric(0),
    gammas02 = numeric(0),
    knots02 = numeric(0),
    gammas12 = numeric(0),
    knots12 = numeric(0),
    #either day 1 or the smallest time in the data, recall that royston & parmar defined on the log scale.
    boundaries = log(c(min(data$V, 1 / 365.25), max(data$V)))
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
    initials$knots01 <- log(c(itm[index[-length(index)]])) #nb: index[length(index)] = boundary knot
  }
  if (k02 > 0) {
    index <- round(seq(0,1,by=1/(k02+1))*itmn)
    initials$knots02 <- log(c(itm[index[-length(index)]])) #nb: index[length(index)] = boundary knot
  }
  if (k12 > 0) {
    index <- round(seq(0,1,by=1/(k12+1))*itmn)
    initials$knots12 <- log(c(itm[index[-length(index)]])) #nb: index[length(index)] = boundary knot
  }

  # try to get sensible initial values
  # transition 01 values
  fit01 <- survfit(Surv(L, R, delta1) ~ 1, data = data)
  cumhaz01 <- rep(fit01$cumhaz, fit01$n.event)
  coefs01 <- lm(
    log(cumhaz01) ~ 0 + nsp(
      log(PFSDY),
      intercept = T,
      knots = initials$knots01,
      Boundary.knots = initials$boundaries
    ) + ATRTN,
    filter(data,delta1 == 1)
  )$coef
  initials$gammas01 <- coefs01[-length(coefs01)]
  initials$theta01 <- coefs01[length(coefs01)]

  # transition 02 values
  data02 <- filter(data, delta1 == 0 & delta2 == 1)
  fit02 <- survfit(Surv(V, delta2) ~ 1, data = data02)
  cumhaz02 <- rep(fit02$cumhaz, fit02$n.event)
  coefs02 <- lm(
    log(cumhaz02) ~ 0 + nsp(
      log(PFSDY),
      intercept = T,
      knots = initials$knots02,
      Boundary.knots = initials$boundaries
    ) + ATRTN,
    data02
  )$coef
  initials$gammas02 <- coefs02[-length(coefs02)]
  initials$theta02 <- coefs02[length(coefs02)]

  # transition 12 values
  data12 <- filter(data, delta0 == 0 & delta2 == 1)
  fit12 <- survfit(Surv(V, delta2) ~ 1, data = data12)
  cumhaz12 <- rep(fit12$cumhaz, fit12$n.event)
  coefs12 <- lm(
    log(cumhaz12) ~ 0 + nsp(
      log(PFSDY),
      intercept = T,
      knots = initials$knots12,
      Boundary.knots = initials$boundaries
    ) + ATRTN,
    data12
  )$coef
  initials$gammas12 <- coefs12[-length(coefs12)]
  initials$theta12 <- coefs12[length(coefs12)]

  initials <- lapply(initials, unname)

  return(initials)
}

#' Compute the log likelihood for a model with illness death and natural cubic spline intensities
#'
#' @param p parameters, as a named list based off \link{RoystonParmar}
#' @param data dataset as per instructions in \link{dc_loglik}
#' @param indiv return a vector of each individual log-likelihood instead of the sum?
#'
#' @return The log-likelihood for the data and parameters provided.
#' @export
royston_parmar.ll <- function(p, data, indiv=FALSE) {
  lik <- dc_loglik(data,
            fns = do.call(
              royston_parmar.fnBuilder,
              p
            ))
  if (indiv) {
    return(lik)
  } else {
    return(sum(lik))
  }
}

#' Functions for a illness-death model with natural cubic spline intensities
#' @description
#' Based on research by Royston & Parmar (2002)
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
#' Each transition's baseline intensity is modeled as \eqn{\exp(S_{hj}(x;\boldsymbol{\gamma}_{hj}))}, where \eqn{x = \log(t)}.
#' We define \eqn{S_{hj}(x)} as a natural cubic spline such that it can be represented in terms of its basis functions \eqn{v_j(x)}
#' \deqn{
#'  S_{hj}(x) = \gamma_{hj0} + \gamma_{hj1}x + \sum^K_{i=1} \gamma_{hji} v_j(x)
#' }
#' with \eqn{K} the number of internal knot points.
#' \code{gammasXX} must contain the \eqn{\boldsymbol{\gamma}_{hj}} values as parameters.
#' Using an example, if we have one knotpoint in the transition from 0 to 1: \code{gamma010}, \code{gamma011}, \code{gamma012} must be provided as parameters.
#'
#' @return A list of functions.
#'
#' @importFrom methods new
#'
#' @rdname RoystonParmar
royston_parmar.fnBuilder <- function(theta01, theta02, theta12, gammas01, knots01, gammas02, knots02, gammas12, knots12, boundaries) {

  fns <- new(RoystonParmarFns,
             theta01, theta02, theta12,
             gammas01, knots01,
             gammas02, knots02,
             gammas12, knots12,
             boundaries)
  fns$P01 <- \(l,r,z) P01(l,r,z,
                          theta01, theta02, theta12,
                          gammas01, knots01,
                          gammas02, knots02,
                          gammas12, knots12,
                          boundaries)
  fns
}
