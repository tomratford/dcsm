#' Update an initial set of parameters with new value
#'
#' @param p Named numeric vector of new values
#' @param inits List of initial values, typically from \link{royston_parmar.initials} or \link{joly.initials}
#'
#' @return a list
#' @export
make_pars2 <- function(p,initials) {
  p2 <- unlist(initials)
  p2[names(p)] <- p
  pars <- make_pars(p2)
  lapply(pars, unname)
}

#' Bundle a vector of parameters into a list of parameters
#'
#' @param p A named vector of parameters, typically from a \link{joly.fit} or
#' \link{royston_parmar.fit} call.
#'
#' @return A list with entries \code{theta01}, \code{theta02}, \code{theta12},
#' \code{gammas01}, \code{knots01} \code{gammas02}, \code{knots02},
#' \code{gammas12}, \code{knots12}, \code{boundaries}.
#' @export
make_pars <- function(p) {
  pars <- list(
    theta01 = numeric(0),
    theta02 = numeric(0),
    theta12 = numeric(0),
    gammas01 = numeric(0),
    knots01 = numeric(0),
    gammas02 = numeric(0),
    knots02 = numeric(0),
    gammas12 = numeric(0),
    knots12 = numeric(0),
    boundaries = numeric(0)
  )
  pars$boundaries <- p[startsWith(names(p),"boundaries")]
  pars$theta01 <- p[["theta01"]]
  pars$theta02 <- p[["theta02"]]
  pars$theta12 <- p[["theta12"]]
  pars$gammas01 <- p[startsWith(names(p),"gammas01")]
  pars$gammas02 <- p[startsWith(names(p),"gammas02")]
  pars$gammas12 <- p[startsWith(names(p),"gammas12")]
  if (any(startsWith(names(p),"knots01"))) {
    pars$knots01 <- p[startsWith(names(p),"knots01")]
  }
  if (any(startsWith(names(p),"knots02"))) {
    pars$knots02 <- p[startsWith(names(p),"knots02")]
  }
  if (any(startsWith(names(p),"knots12"))) {
    pars$knots12 <- p[startsWith(names(p),"knots12")]
  }
  return(pars)
}
