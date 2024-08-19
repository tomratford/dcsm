#' Simulate data based on exponential times to event
#'
#' @param N The number of observations to simulate
#' @param lambda01 The rate for time-to-progression
#' @param lambda02 The rate for time-to-death pre-progression
#' @param lambda12 The rate for time-to-death post-progression
#' @param theta01 The treatment effect for time-to-progression
#' @param theta02 The treatment effect for time-to-death pre-progression
#' @param theta12 The treatment effect for time-to-death post-progression
#' @param pct_no_prg The percent of observations assumed to have not progressed
#' @param K The number of visits
#'
#' @details
#' \eqn{T_{01} \sim \text{Exp}(\lambda_{01})} and
#' \eqn{T_{02} \sim \text{Exp}(\lambda_{02})} are drawn first. If \eqn{T_{01}}
#' is less than \eqn{T_02} then
#' \eqn{T_12 \sim T_01 + \text{Exp}(\lambda_{12} \equiv \text{Gamma}(2, \lambda_{01} \times \lambda_{02}))}.
#' These true values are then used to find \eqn{C_a}, the administrative
#' censoring time as the 80th percentile of \eqn{D = T_{02} \vee T_{12}}.
#'
#'
#' @return A data frame with 13 columns: \itemize{
#'  \item `T01`, `T02`, `T12` - The 'true' times of the event based on the `lambdaXX` parameter.
#'  \item `PFS` - Did the subject progress or die before censoring?
#'  \item `PFSDY` - The day of progression based on when it was first observed.
#'  \item `TRUEDY` - The true day of progression.
#'  \item `L`, `R` - Either the interval progression occurs in if `delta1 == 1`, or the last interval prior to censoring.
#'  \item `V` - Either the time of death if `delta2 == 1`, or the time of right-censoring.
#'  \item `delta0` - Did the patient definitely not progress?
#'  \item `delta1` - Did the patient definitely progress?
#'  \item `delta2` - Did the patient die?
#'  \item `ATRTN` - `0` => control arm, `1` => active arm.
#' }
#'
#' @export
#'
sim_exp <- function(N = 300,
                    # Rates to ensure certain censoring values
                    lambda01 = 4,
                    lambda02 = 1,
                    lambda12 = 1.5,
                    # treatment effects
                    theta01 = log(0.4),
                    theta02 = log(0.4),
                    theta12 = log(1),
                    # %s under certain conditions
                    pct_no_prg = 0.05,
                    # # of visits
                    K = 4) {
  # Assign treatment arm
  Zs <- rbinom(N, 1, 0.5)

  T_01s <- rexp(N, lambda01 * exp(Zs * theta01))
  T_02s <- rexp(N, lambda02 * exp(Zs * theta02))
  T_12s <- T_01s + rexp(N, lambda12 * exp(Zs * theta12)) #gamma

  D <- ifelse(T_01s < T_02s, T_12s, T_02s)

  C_a <- quantile(D[Zs == 1], 0.8) # 20% of deaths censored
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
  pre_prog_death <- T_02s <= T_01s & T_02s < C_a
  delta2[pre_prog_death] <- 1
  Vs[pre_prog_death] <- T_02s[pre_prog_death]
  Rs[pre_prog_death] <- Vs[pre_prog_death]

  # post prog death
  post_prog_death <- (T_01s < T_02s) &
    (T_12s < C_a) & (T_12s > T_01s)
  delta2[post_prog_death] <- 1
  Vs[post_prog_death] <- T_12s[post_prog_death]

  # get progression interval
  for (j in 1:(K + 1)) {
    in_this_int <- (visits[, j] <= T_01s) &
      (T_01s < visits[, j + 1]) &
      !pre_prog_death & (visits[, j + 1] < C_a)
    delta1[in_this_int] <- 1
    Ls[in_this_int] <- visits[, j][in_this_int]
    Rs[in_this_int] <- visits[, j + 1][in_this_int]

    # Pre progression death lhs
    died_in_this_int <- (pre_prog_death) & (visits[, j] < T_02s)
    Ls[died_in_this_int] <- visits[, j][died_in_this_int]

    # post progression death lhs
    died_in_this_int2 <- (post_prog_death) &
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

  # randomly allocate possible delta0s
  potentials <- sample(c(TRUE,FALSE), N, T, c(pct_no_prg,1 - pct_no_prg))
  known_not_prog <- potentials & delta1 == 0
  delta0[known_not_prog] <- 1

  # get crude pfs time
  pfs_event <- as.numeric(delta1 | delta2)
  pfs_time <- rep(C_a, N)
  pfs_time[delta1 == 1] <- Rs[delta1 == 1]
  pfs_time[delta2 == 1] <- Vs[delta2 == 1]
  # true time for comparison
  true_time <- rep(C_a, N)
  true_time[delta1 == 1] <- T_01s[delta1 == 1]
  true_time[delta2 == 1] <- Vs[delta2 == 1]

  # Hide the true values that arent' valid (this will make checking HARDER?)
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
  #attr(dat, "visits") <- visits
  dat
}
