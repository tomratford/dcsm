#' Compute the log likelihood for a set of dual censored observations
#'
#' @param data data frame containing required columns, see details
#' @param fns a list of functions to compute the likelihood
#'
#' @details
#' `data` must contain columns called `l` (optionally `left`), `r` (optionally `right`),
#' `v`, and either `dX`, `d0X`,`deltaX` or `delta0X` for all `x` \in (1,2,3).
#'
#' `fns` must contain entries `int12`, `int02`, `P00`, `logP00`, `P01`, `P11`,
#' `logP11`. These functions should be able to take a vector of inputs.
#'
#' @return vector of log likelihoods
#' @export
#'
#' @importFrom dplyr pull matches
dc_loglik <- function(data, z = "ATRTN", fns = list(
  int12 = \(x,z) 0,
  int02 = \(x,z) 0,
  P00 = \(l,r,z) 0,
  logP00 = \(l,r,z) 1,
  P01 = \(l,r,z) 0,
  P11 = \(l,r,z) 0,
  logP11 = \(l,r,z) 1
)) {
  delta0 <-  pull(data, matches("d(elta)?0?0", perl=T))
  delta1 <-  pull(data, matches("d(elta)?0?1", perl=T))
  delta2 <-  pull(data, matches("d(elta)?0?2", perl=T))
  L <- pull(data, matches("^l(eft)?", perl=T))
  R <- pull(data, matches("^r(ight)?", perl=T))
  V <- pull(data, matches("^v$"))
  Z <- pull(data, z)

  # Potential routes
  i <- delta0 == 0
  j <- delta1 == 0
  dth <- delta2 == 1
  bth <- i & j # these are the annoying ones, can't simplify log of exp
  i <- xor(i,bth) # known progressions
  j <- xor(j,bth) # known not to progress

  # For profiling purposes as this is commonly a bottleneck
  fnsP01 <- fns$P01

  # Factorising out P00(0,l), will be added back on at the end
  zeros <- rep(0,length(L)) # Avoids c++ errors
  logP00s <- fns$logP00(zeros, L, Z)

  # Known progression group
  if (any(i)) { # ensure we have at least one case
    lhs.P01s <- fnsP01(L[i], R[i], Z[i])
    lhs.logP11s <- fns$logP11(R[i], V[i], Z[i])
    tryCatch({
      lhs.loglik <- log(lhs.P01s) + lhs.logP11s
      nd_cnst <- which((i & dth)[i])
    },
    error = \(e) {
      print(e)
      browser()
    })
    if (any(nd_cnst)) { # i.e. did they die and need the intensity added
      tryCatch(
        lhs.loglik[nd_cnst] <- lhs.loglik[nd_cnst] + log(fns$int12(V[nd_cnst], Z[nd_cnst])),
        warning = \(w) {
          print(w)
          print(V[nd_cnst][which(is.nan(log(fns$int12(V[nd_cnst], Z[nd_cnst]))))])
        }
      )
    }
  } else {
    lhs.loglik <- 0
  }

  # known not to progress
  if (any(j)) { # Ensure we have at least one case
    rhs.loglik <- log(fns$P00(L[j], V[j], Z[j]))
    nd_cnst2 <- which((j & dth)[j])
    if (any(nd_cnst2)) {
      rhs.loglik[nd_cnst2] <- rhs.loglik[nd_cnst2] + log(fns$int02(V[nd_cnst2], Z[nd_cnst2]))
    }
  } else
  {
    rhs.loglik <- 0
  }

  # Dual censored observations
  both.lhs.P01s <- fnsP01(L[bth], R[bth], Z[bth])
  both.lhs.P11s <- fns$P11(R[bth], V[bth], Z[bth])
  both.lhs.lik <- both.lhs.P01s * both.lhs.P11s

  both.rhs.lik <- fns$P00(L[bth], V[bth], Z[bth])

  nd_cnst3 <- which((bth & dth)[bth])
  if (any(nd_cnst3)) {
    both.lhs.lik[nd_cnst3] <- both.lhs.lik[nd_cnst3] * fns$int12(V[nd_cnst3], Z[nd_cnst3])
    both.rhs.lik[nd_cnst3] <- both.rhs.lik[nd_cnst3] * fns$int02(V[nd_cnst3], Z[nd_cnst3])
  }
  both.lik <- both.lhs.lik + both.rhs.lik
  # in some cases we have negative probability
  # this typically happens in spline models so force to zero as I assume this
  # means shit optimisation with the gammas having dumb values.
  both.lik[both.lik < 0] <- 0
  both.loglik <- log(both.lik)

  # sum all together and return vector
  loglik <- numeric(length = nrow(data))
  loglik[i] <- lhs.loglik
  loglik[j] <- rhs.loglik
  loglik[bth] <- both.loglik
  loglik <- loglik + logP00s
  loglik
}
