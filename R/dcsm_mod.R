#' Print the treatment effect estimates for a dual-censored model
#'
#' @param x A \code{dcsm_mod} object
#' @param ... Not used
#'
#' @return Nothing, prints three estimates of the treatment effect for each transition
#' @export
print.dcsm_mod <- function(x, ...) {
  thetas <- x$par[c("theta01", "theta02", "theta12")]
  cat(paste(names(thetas), collapse = "\t"))
  cat("\n")
  cat(paste(format(thetas, digits = 3), collapse = "\t"))
}

#' Plot the survivor curve for a dual-censored model
#'
#' @param x A \code{dcsm_mod} object
#' @param treat_col The colour of the treatment effect line, set to \code{NULL} to only plot the baseline survivor function.
#' @param add Add this to a existing plot?
#' @param ... Passed to \link{curve}
#'
#' @return Nothing, plots to output
#' @export
plot.dcsm_mod <- function(x,
                          treat_col = "blue",
                          add = F,
                          CI = F,
                          ...) {
  if (attr(x, "dist") == "weibull") {
    pars <- as.list(x$par)
  } else {
    if (exists(attr(x, "initials"))) {
      pars <- make_pars2(x$par, attr(x, "initials"))
    } else {
      stop("No initials found, running this will brick your R session")
    }
  }
  fns <- switch(
    attr(x, "dist"),
    "weibull" = do.call(weib.fnBuilder, pars),
    "joly" = do.call(joly.fnBuilder, pars),
    "royston_parmar" = do.call(royston_parmar.fnBuilder, pars),
    stop("Invalid distribution in object")
  )
  if (CI) {
    se <- dcsm_se(x, fns)
  }
  curve(fns$P00(rep(0, length(x)), x, rep(0, length(x))), add = add, ...)
  if (CI) {
    curve(fns$P00(rep(0, length(x)), x, rep(0, length(x))) + qnorm(0.975) * se(x, 0), add = T, lty=2, ...)
    curve(fns$P00(rep(0, length(x)), x, rep(0, length(x))) - qnorm(0.975) * se(x, 0), add = T, lty=2, ...)
  }
  if (!is.null(treat_col)) {
    curve(fns$P00(rep(0, length(x)), x, rep(1, length(x))), add = T, col =
            treat_col, ...)
    if (CI) {
      curve(fns$P00(rep(0, length(x)), x, rep(1, length(x))) + qnorm(0.975) * se(x, 1), add = T, lty=2, ...)
      curve(fns$P00(rep(0, length(x)), x, rep(1, length(x))) - qnorm(0.975) * se(x, 1), add = T, lty=2, ...)
    }
  }
}

#' Return a function that calculate the standard error of the survivor curve
#'
#' @param x A \code{dscm_mod} object
#' @param fns List of functions used to generate SE
#'
#' @return A function to find the standard error at time \code{tm} with treatment arm \code{z}.
#' @importFrom splines2 nsp isp
dcsm_se <- function(x, fns) {
  if (!exists("hessian", weib_mod)) {
    stop("Hessian required, rerun with `hessian=T`")
  }
  Io <- solve(-x$hessian)
  switch(
    attr(x, "dist"),
    "weibull" = function(tm, z = 0) {
      delta <- matrix(0, ncol = 9, nrow = length(tm))

      delta[, 1] <- -tm ^ x$par[["gamma01"]] * exp(z * x$par[["theta01"]]) * fns$P00(0, tm, z)
      delta[, 2] <- -fns$int01(tm, z) * fns$P00(0, tm, z)
      delta[, 3] <- z * (-x$par[["lambda01"]] * (tm ^ x$par[["gamma01"]]) * exp(z * x$par[["theta01"]]) * fns$P00(0, tm, z))
      delta[, 4] <- -tm ^ x$par[["gamma02"]] * exp(z * x$par[["theta02"]]) * fns$P00(0, tm, z)
      delta[, 5] <- -fns$int02(tm, z) * fns$P00(0, tm, z)
      delta[, 6] <- z * (-x$par[["lambda02"]] * (tm ^ x$par[["gamma02"]]) * exp(z * x$par[["theta02"]]) * fns$P00(0, tm, z))

      diag(sqrt(delta %*% Io %*% t(delta)))
    },
    "royston_parmar" = function(tm, z = 0) {
      if (length(z) == 1) {
        z <- rep(z, length(tm))
      }
      zeros <- rep(0, length(tm))

      g01n <- length(startswith(names(x$par), "gammas01"))
      g02n <- length(startswith(names(x$par), "gammas02"))

      delta <- matrix(0, ncol = 3 + g01n + g02n + g12n, nrow = length(tm))

      P00s <- fns$P00(zeros, tm, z)
      A01 <- fns$Int01(tm, z)
      A02 <- fns$Int02(tm, z)

      delta[, 1] <- -z * A01 * P00s
      delta[, 2] <- -z * A02 * P00s

      delta[, 4:(4 + g01n)] <- nsp(
        tm,
        knots = pars$knots01,
        Boundary.knots = pars$boundaries,
        intercept = T,
        derivs = 1
      ) * A01 * P00s
      delta[, (4 + g01n):(4 + g01n + g02n)] <- nsp(
        tm,
        knots = pars$knots02,
        Boundary.knots = pars$boundaries,
        intercept = T,
        derivs = 1
      ) * A02 * P00s

      diag(sqrt(delta %*% Io %*% t(delta)))
    },
    "joly" = function(tm, z = 0) {
      if (length(z) == 1) {
        z <- rep(z, length(tm))
      }
      zeros <- rep(0, length(tm))

      g01n <- length(startswith(names(x$par), "gammas01"))
      g02n <- length(startswith(names(x$par), "gammas02"))

      delta <- matrix(0, ncol = 3 + g01n + g02n + g12n, nrow = length(tm))

      P00s <- fns$P00(zeros, tm, z)
      A01 <- fns$Int01(tm, z)
      A02 <- fns$Int02(tm, z)

      delta[, 1] <- -z * A01 * P00s
      delta[, 2] <- -z * A02 * P00s

      delta[, 4:(4 + g01n)] <- isp(
        tm,
        knots = pars$knots01,
        Boundary.knots = pars$boundaries,
        intercept = T,
        derivs = 1
      ) * exp(z * x$par[["theta01"]]) * P00s
      delta[, (4 + g01n):(4 + g01n + g02n)] <- isp(
        tm,
        knots = pars$knots02,
        Boundary.knots = pars$boundaries,
        intercept = T,
        derivs = 1
      ) * exp(z * x$par[["theta02"]]) * P00s

      diag(sqrt(delta %*% Io %*% t(delta)))
    },
    stop("Invalid distribution in object")
  )
}
