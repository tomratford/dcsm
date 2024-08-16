#' Perform a simple Newton-Raphson procedure for generic f and df functions.
#'
#' @param f function taking two inputs, \code{x} and \code{z}.
#' @param df first derivative of function \code{f}.
#' @param z treatment arm
#' @param t1 Starting value for \code{f} and \code{df}
#'
#' @return single numeric of the true \eqn{t}
#' @export
newton_raphson <- function(f = \(x,z) (x-1)^2, df = \(x,z) 2*(x-1), z, t1 = 1) {
  x <- t1
  count <- 0
  maxiter <- 100
  while (abs(f(x,z)) > 1e-15) {
    x <- x - f(x,z)/df(x,z)
    count <- count + 1
    if (x < 0) {
      stop("negative t isn't possible")
    }
    if (count > maxiter) {
      stop("too many iterations")
    }
  }
  return(x)
}

sim_nonlinear <- function(N = 300,
                          # baseline intensity functions
                          a01 = \(x) x^(4/5),
                          A01 = \(x) (5/9)*x^(9/5),
                          a02 = \(x) 3*x/4,
                          A02 = \(x) (3*(x^2))/8,
                          a12 = \(x) (3*x/2)^(5/4),
                          A12 = \(x) ((2/3)^(3/4))*(x^(9/4)),
                          # treatment effects
                          theta01 = log(0.4),
                          theta02 = log(0.4),
                          theta12 = log(1),
                          K = 4,
                          pct_no_prg = 0.05,
                          silent=T) {
  # helper functions (these are)
  P00 <- \(x,z) exp(-A01(x)*exp(z*theta01) - A02(x)*exp(z*theta02))
  dP00 <- \(x,z) (-a01(x)*exp(z*theta01) - a02(x)*exp(z*theta02)) * P00(x,z)
  P11 <- \(s) (\(x,z) exp(-A12(x)*exp(z*theta12) - A12(s)*exp(z*theta12)))
  dP11 <- \(s) (\(x,z) (-a12(x)*exp(z*theta12)) * P11(s)(x,z))
  P01 <- \(x,z) Vectorize(\(x) integrate(\(v) P00(v,z)*a01(v)*exp(z*theta01)*P11(v)(x,z), 0, x)$value)(x)
  #dP01 <- \(x,z) P00(x,z) * a01(x) * exp(z*theta01)
  dP01 <- \(x,z) P00(x,z) * a01(x) * exp(z*theta01) - a12(2) * exp(z*theta12) * P01(x,z)
  P02 <- \(x,z) 1 - P01(x,z) - P00(x,z)
  dP02 <- \(x,z) -dP01(x,z) - dP00(x,z)

  # Assign treatment arm
  Zs <- rbinom(N, 1, 0.5)

  # calculate progression/pre-progression death

  # progression time

  u_01 <- runif(N)
  T_01s <- numeric(N)
  for (i in 1:N) {
    temp <- try(newton_raphson(\(x, z) - A01(x) * exp(z * theta01) - log(u_01[i]),
                               \(x, z) - a01(x) * exp(z * theta01),
                               Zs[i]), silent)
    if (class(temp) == "try-error") {
      T_01s[i] <- NA
    } else {
      T_01s[i] <- temp
    }
  }

  # pre-progression death time
  u_02 <- runif(N)
  T_02s <- numeric(N)
  for (i in 1:N) {
    temp <- try(newton_raphson(\(x, z) - A02(x) * exp(z * theta02) - log(u_02[i]),
                               \(x, z) - a02(x) * exp(z * theta02),
                               Zs[i]), silent)
    if (class(temp) == "try-error") {
      T_02s[i] <- NA
    } else {
      T_02s[i] <- temp
    }
  }

  # post progression death
  u_12 <- runif(N,max = 1 - u_01)
  T_12s <- numeric(N)
  for (i in 1:N) {
    if (is.na(T_01s[i]) || (!is.na(T_01s[i]) && !is.na(T_02s[i]) && T_02s[i] < T_01s[i])) {
      T_12s[i] <- NA
      next
    }
    temp <- try(newton_raphson(\(x, z) -A12(x)*exp(z*theta12) - log(u_12[i]),
                               \(x,z) -a12(x)*exp(z*theta12),
                               Zs[i]), silent)
    if (class(temp) == "try-error") {
      T_12s[i] <- NA
    } else {
      T_12s[i] <- T_01s[i] + temp
    }
  }

  D <- ifelse(T_01s < T_02s, T_12s, T_02s)

  C_a <- quantile(D[Zs == 1], 0.8, na.rm=T) # 20% of deaths censored
  names(C_a) <- NULL

  # Get visit times
  visits <- matrix(nrow = N, ncol = K + 2)
  visits[, 1] <- 0
  visits[, K + 2] <- C_a
  for (i in 1:K) {
    visits[, i + 1] <- i * (C_a / K) + rnorm(N, 0, C_a / (20 * K))
  }
  visits[visits < 0] <- 0
  visits <- t(apply(visits, 1, sort)) # ensure a sensible order of visits
  visits[visits > C_a] <- C_a # truncate any values greater censoring time

  # Inits
  delta0 <- vector("numeric", N)
  delta1 <- vector("numeric", N)
  delta2 <- vector("numeric", N)
  Ls <- apply(visits, 1, \(row) max(row[row < C_a]))
  Rs <- visits[, K + 2]
  Vs <- Rs

  # pre progression deaths
  pre_prog_death <- !is.na(T_01s) & !is.na(T_02s) & T_02s <= T_01s & T_02s < C_a
  delta2[pre_prog_death] <- 1
  Vs[pre_prog_death] <- T_02s[pre_prog_death]
  Rs[pre_prog_death] <- Vs[pre_prog_death]

  # post prog death
  post_prog_death <- !is.na(T_01s) & !is.na(T_02s) & !is.na(T_12s) & (T_01s < T_02s) &
    (T_12s < C_a) & (T_12s > T_01s)
  delta2[post_prog_death] <- 1
  Vs[post_prog_death] <- T_12s[post_prog_death]

  # get progression interval
  for (j in 1:(K + 1)) {
    in_this_int <- !is.na(T_01s) &
      (visits[, j] <= T_01s) &
      (T_01s < visits[, j + 1]) &
      !pre_prog_death & (visits[, j + 1] < C_a)
    delta1[in_this_int] <- 1
    Ls[in_this_int] <- visits[, j][in_this_int]
    Rs[in_this_int] <- visits[, j + 1][in_this_int]

    # Pre progression death lhs
    died_in_this_int <- !is.na(T_02s) & (pre_prog_death) & (visits[, j] < T_02s)
    Ls[died_in_this_int] <- visits[, j][died_in_this_int]

    # post progression death lhs
    died_in_this_int2 <- !is.na(T_12s) & (post_prog_death) &
      (in_this_int) & (visits[, j] < T_12s) & (T_12s < visits[, j + 1])
    delta1[died_in_this_int2] <- 0 # have to reset as we wouldn't know their progression status actually (would appear as if they just died pre-progression)
    Ls[died_in_this_int2] <- visits[, j][died_in_this_int2]
    Rs[died_in_this_int2] <- T_12s[died_in_this_int2]
  }
  # Hacky fix:
  # if R > V then set R to V provided it is less than L otherwise throw an error.
  # This is a weird edge case in the post progression case, don't really understand why it happens
  hack <- Rs > Vs
  if (any(Ls[hack] > Vs[hack])) {
    print(Ls[hack])
    print(Vs[hack])
    stop("This shouldn't happen. An L is greater than a V")
  }
  Rs[hack] <- Vs[hack]


  # potential known not progressed
  potentials <- which(delta0 == 0 & delta1 == 0)
  delta0[sample(potentials, ceiling(pct_no_prg * length(potentials)))] <- 1

  # get crude pfs time
  pfs_event <- as.numeric(delta1 | delta2)
  pfs_time <- rep(C_a, N)
  pfs_time[delta1 == 1] <- Rs[delta1 == 1]
  pfs_time[delta2 == 1] <- Vs[delta2 == 1]
  # true time for comparison
  true_time <- rep(C_a, N)
  true_time[delta1 == 1] <- T_01s[delta1 == 1]
  true_time[delta2 == 1] <- Vs[delta2 == 1]

  # Hide the true values that arent' valid (this will make checking easier?)
  # T_12s[(delta1 != 1) & (delta2 == 1)] <- NA
  # T_02s[(delta1 == 1) & (delta2 == 1)] <- NA
  # T_02s[delta2 != 1] <- NA
  # T_12s[delta2 != 1] <- NA

  dat <- data.frame(
    T01 = T_01s,
    T02 = T_02s,
    T12 = T_12s,
    PFS = pfs_event,
    PFSDY = pfs_time,
    TRUEDY = true_time,
    L = Ls,
    R = Rs,
    V = Vs,
    delta0 = delta0,
    delta1 = delta1,
    delta2 = delta2,
    ATRTN = Zs
  )
  attr(dat, "Ca") <- C_a
  attr(dat, "S") <- P00
  dat
}
