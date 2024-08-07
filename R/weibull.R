#' Fit a model based on starting parameters
#'
#' @param data Data set to fit to, must have entries as described in \link{dc_loglik}
#' @param initial Vector of initial parameters
#' @param control passed to \link{optim}, defaults to change the scale of the log-likelihood function
#' @param method passed to \link{optim}, defaults to Nelder-Mead.
#' @param ... additional parameters passed to \link{optim}
#'
#' @return A \link{optim} output
#' @export
weib.fit <- function(data,
                     initial = c(
                       "lambda01" = 1,
                       "gamma01" = 1,
                       "theta01" = 0,
                       "lambda02" = 1,
                       "gamma02" = 1,
                       "theta02" = 0,
                       "lambda12" = 1,
                       "gamma12" = 1,
                       "theta12" = 0
                     ),
                     control = list(fnscale = -1),
                     method = "Nelder-Mead",
                     ...) {
  optim(initial, weib.ll, data = data, control = control, ...)
}

#' Compute the log likelihood for a model with illness death and Weibull intensities
#'
#' @param p parameters, see \link{weibull}
#' @param data dataset as per instructions in \link{dc_loglik}
#'
#' @return The log-likelihood for the data and parameters provided.
#' @export
weib.ll <- function(p, data) {
  if (any(p[startsWith(names(p), "lambda")] < 0) || any(p[startsWith(names(p), "gamma")] < 0)) {
    return(-1e10)
  }
  dc_loglik(dat,
            fns = do.call(
              weib.fnBuilder,
              as.list(p)
            )) |>
    sum()
}

#' Functions for a illness-death model with Weibull intensities
#'
#' @param lambda see details, can be suffixed based on the transition
#' @param gamma see details, can be suffixed based on the transition
#' @param theta see details, can be suffixed based on the transition
#'
#' @details
#' Weibull is formulated such that its hazard is
#' \deqn{
#' H(t) = \lambda \gamma t^{\gamma-1}
#' }
#' this then models the baseline intensity function for each transition \eqn{h \to j}
#' \deqn{
#' \alpha_{hj}(t;z) = \alpha_{hj}(t;0) \exp(z\theta_{hj})
#' }
#'
#' @rdname weibull
#' @return A list of functions, to be used in `?weib.ll`
weib.fnBuilder <- function(lambda01,
                           gamma01,
                           theta01,
                           lambda02,
                           gamma02,
                           theta02,
                           lambda12,
                           gamma12,
                           theta12) {
  args <- as.list(environment()) # helpful for debugging
  rtn <- list(
    int12 = weib.int(lambda12, gamma12, theta12),
    int02 = weib.int(lambda02, gamma02, theta02),
    P00 = weib.P00(lambda01, gamma01, theta01, lambda02, gamma02, theta02),
    logP00 = weib.logP00(lambda01, gamma01, theta01, lambda02, gamma02, theta02),
    P01 = weib.P01(
      lambda01,
      gamma01,
      theta01,
      lambda02,
      gamma02,
      theta02,
      lambda12,
      gamma12,
      theta12
    ),
    P11 = weib.P11(lambda12, gamma12, theta12),
    logP11 = weib.logP11(lambda12, gamma12, theta12)
  )
  attr(rtn, "call") <- args
  rtn
}

#' @rdname weibull
weib.int <- function(lambda, gamma, theta) {
  function(t, z) {
    exp(z * theta) * lambda * gamma * (t ^ (gamma - 1))
  }
}

#' @rdname weibull
weib.cumint <- function(lambda, gamma, theta) {
  function(t, z)
    exp(z * theta) * lambda * (t ^ gamma)
}

#' @rdname weibull
weib.P00 <- function(lambda01,
                     gamma01,
                     theta01,
                     lambda02,
                     gamma02,
                     theta02) {
  inner <- function(u, z)
    (lambda01 * (u ^ gamma01)) * exp(z * theta01) + (lambda02 * (u ^
                                                                   gamma02)) * exp(z * theta02)

  function(l, r, z) {
    exp(inner(l, z) - inner(r, z))
  }
}

#' @rdname weibull
weib.logP00 <- function(lambda01,
                        gamma01,
                        theta01,
                        lambda02,
                        gamma02,
                        theta02) {
  inner <- function(u, z)
    (lambda01 * u ^ gamma01) * exp(z * theta01) + (lambda02 * u ^
                                                     gamma02) * exp(z * theta02)

  function(l, r, z) {
    inner(l, z) - inner(r, z)
  }
}

#' @rdname weibull
weib.P11 <- function(lambda12, gamma12, theta12) {
  inner <- function(u, z)
    (lambda12 * u ^ gamma12) * exp(z * theta12)

  function(l, r, z) {
    exp(inner(l, z) - inner(r, z))
  }
}

#' @rdname weibull
weib.logP11 <- function(lambda12, gamma12, theta12) {
  inner <- function(u, z)
    (lambda12 * u ^ gamma12) * exp(z * theta12)

  function(l, r, z) {
    inner(l, z) - inner(r, z)
  }
}

#' @rdname weibull
weib.P01 <- function(lambda01,
                     gamma01,
                     theta01,
                     lambda02,
                     gamma02,
                     theta02,
                     lambda12,
                     gamma12,
                     theta12) {
  # P00 <- weib.P00(lambda01, gamma01, theta01, lambda02, gamma02, theta02)
  # P11 <- weib.P11(lambda12, gamma12, theta12)
  # int01 <- weib.int(lambda01, gamma01, theta01)

  cumint12 <- weib.cumint(lambda12, gamma12, theta12)
  cumint01 <- weib.cumint(lambda01, gamma01, theta01)
  cumint02 <- weib.cumint(lambda02, gamma02, theta02)

  Vectorize(function(l, r, z) {
    factor <- gamma01 * lambda01 * exp(cumint01(l,z) + cumint02(l,z) - cumint12(r,z)) * exp(z*theta01)
    int <- integrate(\(v) exp(v)^(gamma01) * exp(cumint12(exp(v),z) - cumint01(exp(v),z) - cumint02(exp(v),z)), log(l), log(r))
    int$value * factor
  })
}
