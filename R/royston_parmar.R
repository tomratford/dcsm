#' Compute the log likelihood for a model with illness death and natural cubic spline intensities
#'
#' @param p parameters, as a named list based off \link{RoystonParmar}
#' @param data dataset as per instructions in \link{dc_loglik}
#'
#' @return The loglikelihood for the data and parameters provided.
#' @export
royston_parmar.ll <- function(p, data) {
  dc_loglik(data,
            fns = do.call(
              royston_parmar.fnBuilder,
              p
            )) |>
    sum()
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
