# Run the model fitting procedure for each possible model,
# Save the final model objects to users home to manipulate later.

devtools::load_all()
dat <- AMGEN_20050203

# Weibull model
weib_mod <- weib.fit(dat)

# Royston-Parmar model
# Forward fit from 0 knots
cat("Fitting Royston Model\n")
rknots <- c(0, 0, 0)
fit_royston <- \(ks) {
  v <- try(royston_parmar.fit(dat, ks[1], ks[2], ks[3], control=list(fnscale=-1,maxit=500,trace=10)))
  if (class(v) == "try-error") {
    v <- list(value = NA)
  }
  cat("Fit model with knots", format(ks), "\n")
  v
}
royston_mod <- fit_royston(rknots)
history_royston_mods <- list(royston_mod)
history_rknots <- list(rknots)
history_pvals <- c()
while (TRUE) {
  cat("Trying new models\n")
  possible_mods <- vector("list", 3)
  pvals <- vector("numeric", 3)
  possible_rknots <- vector("list", 3)

  possible_rknots[[1]] <- rknots + c(1, 0, 0)
  possible_rknots[[2]] <- rknots + c(0, 1, 0)
  possible_rknots[[3]] <- rknots + c(0, 0, 1)
  cat(format(possible_rknots),"\n")
  # restrict to a maximum of three knots
  possible_rknots <- lapply(possible_rknots,
                            \(rk) mapply(\(x,y) min(x,y), rk, c(3,3,3)))

  possible_mods <- lapply(possible_rknots, fit_royston)
  cat("Fitted new models\n")

  values <- sapply(possible_mods, \(x) x$value)
  Ts <- 2*(values - royston_mod$value)
  pvals <- pchisq(Ts,1,lower.tail = F)

  if (any(pvals <= 0.05)) {
    replace <- which(max(Ts) == Ts)
    cat("Pvals:", format(pvals),"\n")
    cat("Replaced with new model:",format(possible_rknots[[replace]]),"\n")
    royston_mod <- possible_mods[[replace]]
    rknots <- possible_rknots[[replace]]
    history_royston_mods[[length(history_royston_mods) + 1]] <- royston_mod
    history_rknots[[length(history_rknots) + 1]] <- rknots
    history_pvals <- c(history_pvals, pvals[[replace]])
  } else {
    break
  }
}
cat("Royston Model done\n")

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
history_joly_mods <- list(royston_mod)
history_jknots <- list(jknots)
history_aic <- c()
while (TRUE) {
  possible_mods <- vector("list", 3)
  AIC <- vector("numeric", 3)
  possible_jknots <- vector("list", 3)

  possible_jknots[[1]] <- jknots + c(1, 0, 0)
  possible_jknots[[2]] <- jknots + c(0, 1, 0)
  possible_jknots[[3]] <- jknots + c(0, 0, 1)

  possible_mods <- parallel::mclapply(possible_rknots, fit_royston, mc.cores = 3)

  AIC <- sapply(possible_mods, \(x) x$AIC) - joly_mod$AIC

  if (any(AIC >= 4)) {
    replace <- which(max(AIC) == AIC)
    cat("New model:",format(possible_jknots[[replace]]),"\n")
    joly_mod <- possible_mods[[replace]]
    jknots <- possible_jknots[[replace]]
    history_joly_mods[[length(history_royston_mods) + 1]] <- joly_mod
    history_jknots[[length(history_rknots) + 1]] <- jknots
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
