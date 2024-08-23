dat <- sim_exp()

# Royston-Parmar model
# Forward fit from 0 knots
rknots <- c(0, 0, 0)
fit_royston <- \(ks) {
  v <- try(royston_parmar.fit(dat, ks[1], ks[2], ks[3]))
  if (class(v) == "try-error") {
    v <- list(value = NA)
  }
  v
}
royston_mod <- fit_royston(rknots)
history_royston_mods <- list(royston_mod)
history_rknots <- list(rknots)
history_pvals <- c()
print("Model fitting starts now")
while (TRUE) {
  possible_mods <- vector("list", 3)
  pvals <- vector("numeric", 3)
  possible_rknots <- vector("list", 3)

  possible_rknots[[1]] <- rknots + c(1, 0, 0)
  possible_rknots[[2]] <- rknots + c(0, 1, 0)
  possible_rknots[[3]] <- rknots + c(0, 0, 1)

  # restrict to a maximum of three knots
  possible_rknots <- lapply(possible_rknots,
                            \(rk) mapply(\(x,y) min(x,y), rk, c(3,3,3)))

  possible_mods <- parallel::mclapply(possible_rknots, fit_royston, mc.cores = 16)

  values <- sapply(possible_mods, \(x) x$value)
  Ts <- 2*(values - royston_mod$value)
  pvals <- pchisq(Ts,1,lower.tail=F)

  if (any(pvals <= 0.05)) {
    replace <- which(max(Ts) == Ts)
    royston_mod <- possible_mods[[replace]]
    rknots <- possible_rknots[[replace]]
    history_royston_mods[[length(history_royston_mods) + 1]] <- royston_mod
    history_rknots[[length(history_rknots) + 1]] <- rknots
    history_pvals <- c(history_pvals, pvals[[replace]])
  } else {
    break
  }
}
