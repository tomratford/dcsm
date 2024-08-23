dat <- sim_nonlinear()

# Penalised I-spline model
# Forward fitting with AIC
cat("Fitting Joly Model\n")
jknots <- c(0, 0, 0)
fit_joly <- \(ks) {
  v <- try(joly.fit(dat, ks[1], ks[2], ks[3], 100, 100, 100),T)
  if (class(v) == "try-error") {
    v <- list(AIC = NA)
  }
  v
}
joly_mod <- fit_joly(jknots)
history_joly_mods <- list(joly_mod)
history_jknots <- list(jknots)
history_aic <- c()
while (TRUE) {
  possible_mods <- vector("list", 3)
  AIC <- vector("numeric", 3)
  possible_jknots <- vector("list", 3)

  possible_jknots[[1]] <- jknots + c(1, 0, 0)
  possible_jknots[[2]] <- jknots + c(0, 1, 0)
  possible_jknots[[3]] <- jknots + c(0, 0, 1)

  possible_mods <- lapply(possible_jknots, fit_joly)

  AIC <- sapply(possible_mods, \(x) x$AIC) - joly_mod$AIC

  if (any(AIC >= 4)) {
    replace <- which(max(AIC) == AIC)
    cat("New model:",format(possible_jknots[[replace]]),"\n")
    joly_mod <- possible_mods[[replace]]
    jknots <- possible_jknots[[replace]]
    history_joly_mods[[length(history_joly_mods) + 1]] <- joly_mod
    history_jknots[[length(history_jknots) + 1]] <- jknots
    history_aic <- c(history_aic, AIC[[replace]])
  } else {
    break
  }
}
cat("Joly Model done\n")

save(
  list = c("weib_mod", "royston_mod", "joly_mod"),
  file = file.path(Sys.getenv("HOME"), "50203_mods.rda")
)
save(
  list = c(
    "history_royston_mods",
    "history_rknots",
    "history_pvals",
    "history_joly_mods",
    "history_jknots",
    "history_aic"
  ),
  file = file.path(Sys.getenv("HOME"), "50203_history.rda")
)
