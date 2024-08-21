#' Fit a dataset to Royston & Parmar spline model.
#'
#' @param data Data set to model from
#' @param k01 Number of interior knots in the 0 -> 1 transition
#' @param k02 Number of interior knots in the 0 -> 2 transition
#' @param k12 Number of interior knots in the 1 -> 2 transition
#' @param debug Print every parameter tried?
#' @param initials Alternative list of initial values to try
#' @param control Passed to \link{optim}
#' @param method Passed to \link{optim}, defaults to \code{"BFGS"}
#' @param ... Passed to \link{optim}
#'
#' @return A \link{optim} object.
#' @export
royston_parmar.fit <- function(data,
                               k01,
                               k02,
                               k12,
                               debug = F,
                               initials = NULL,
                               control = list(fnscale = -1, maxit = 500),
                               method = "BFGS",
                               ...
) {
  if (is.null(initials))
    initials <- royston_parmar.initials(data, k01, k02, k12)
  p_init <- unlist(initials[c("theta01","theta02","theta12","gammas01","gammas02","gammas12")])
  if (debug) print(initials)
  opt_out <- optim(
    unlist(p_init),
    \(p) {
      pars <- make_pars2(p,initials)
      if (debug) print(pars[c("theta01", "theta02", "theta12", "gammas01","gammas02","gammas12")])
      res <- tryCatch(royston_parmar.ll(pars, data),
               error = \(e) {
                 print(e)
                 return(-1e10*length(data$V))
               })
      if (debug) print(res)
      if (!is.finite(res)) {
        return(-1e10*length(data$V))
      }
      res
    },
    control = control,
    method = method,
    ...
  )
  attr(opt_out, "dist") <- "royston_parmar"
  attr(opt_out, "initials") <- initials
  class(opt_out) <- c("dcsm_mod")
  opt_out
}

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
#' @import dplyr
#' @importFrom survival Surv survfit
#' @importFrom splines2 nsp
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
    boundaries = log(c(min(data$V, 1/365.25), max(data$V)))
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
  data01 <- filter(data,delta0 == 0)
  surv01 <- survfit(Surv(L,R,delta1,type="interval") ~ ATRTN,data01)
  cumhaz_raw <- -log(surv01$surv)
  cumhaz <- unlist(with(surv01, mapply(\(n,c) rep(c,n), n.event, cumhaz_raw)))
  cumhaz_is_finite <- is.finite(cumhaz)
  data01_events <- filter(data, delta1 == 1)
  coefs01 <- solve_qp_problem(cumhaz[cumhaz_is_finite],
                              log(data01_events$R)[cumhaz_is_finite],
                              data01_events$ATRTN[cumhaz_is_finite],
                              initials$knots01,
                              initials$boundaries)
  initials$gammas01 <- coefs01[-length(coefs01)]
  initials$theta01 <- coefs01[length(coefs01)]

  # transition 02 values
  data02 <- filter(data, delta1 == 1)
  surv02 <- survfit(Surv(V,delta2) ~ ATRTN,data02)
  cumhaz <- unlist(with(surv02, mapply(\(n,c) rep(c,n), n.event, cumhaz)))
  data02_events <- filter(data02, delta2 == 1)
  coefs02 <- solve_qp_problem(cumhaz,
                              log(data02_events$V),
                              data02_events$ATRTN,
                              initials$knots02,
                              initials$boundaries)
  initials$gammas02 <- coefs02[-length(coefs02)]
  initials$theta02 <- coefs02[length(coefs02)]

  # transition 12 values
  data12 <- filter(data, delta0 == 0)
  surv12 <- survfit(Surv(V,delta2) ~ ATRTN,data12)
  cumhaz <- unlist(with(surv12, mapply(\(n,c) rep(c,n), n.event, cumhaz)))
  data12_events <- filter(data12, delta2 == 1)
  coefs12 <- solve_qp_problem(cumhaz,
                              log(data12_events$V),
                              data12_events$ATRTN,
                              initials$knots12,
                              initials$boundaries)
  initials$gammas12 <- coefs12[-length(coefs12)]
  initials$theta12 <- coefs12[length(coefs12)]
  initials <- lapply(initials, unname)

  return(initials)
}

#' Find initial coefficients for Royston-Parmar
#' Taken nearly directly from flexsurv (c. Christopher Jackson)
#'
#' @param y cumhaz
#' @param x log vector of time (R/V)
#' @param z vector of coefficients
#' @param initial initial values so far (for knots/boundaries)
#'
#' @return vector of solutions
#' @importFrom quadprog solve.QP
solve_qp_problem <- function(y, x, z, knots, boundaries) {
  b <- splines2::nsp(
    x,
    knots = knots,
    Boundary.knots = boundaries,
    intercept = T
  )
  Xq <- cbind(b, z)

  kx <- x
  kr <- diff(range(kx))
  kx[1] <- x[1] - 0.01*kr
  db <- splines2::nsp(
    kx,
    knots = knots,
    Boundary.knots = boundaries,
    intercept = T,
    deriv=T
  )
  dXq <- cbind(db, rep(0, length(z)))

  Dmat <- t(Xq) %*% Xq
  if (!all(eigen(Dmat)$values > 0)) {
    stop("no sol")
  }
  solve.QP(
    Dmat = Dmat,
    dvec = t(t(y) %*% Xq),
    Amat = t(dXq),
    bvec = rep(1, length(y))
  )$solution
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

  fns <- new(IllnessDeath,
             "RoystonParmar",
             theta01, theta02, theta12,
             gammas01, knots01,
             gammas02, knots02,
             gammas12, knots12,
             boundaries)
  # fns$P01 <- \(l,r,z) P01(l,r,z,"RoystonParmar",
  #                         theta01, theta02, theta12,
  #                         gammas01, knots01,
  #                         gammas02, knots02,
  #                         gammas12, knots12,
  #                         boundaries)
  fns$P01 <- Vectorize(\(l, r, z) tryCatch(
    integrate(
      \(v) res <- fns$P01Integrand(v, rep(l, length(v)), rep(r, length(v)), rep(z, length(v))),
      lower = log(max(l, 1e-9)),
      upper = log(r)
    )$value,
    error = \(e) {
      # print(e)
      # curve(fns$P01Integrand(x, rep(l, length(x)), rep(r, length(x)), rep(z, length(x))),
      #       from = log(1e-9),
      #       to = log(r))
      #stop(e)
      return(-1e10) # I'm not sure if this is what it 'should' return, but it works sofar
    }
  ))
  fns
}
