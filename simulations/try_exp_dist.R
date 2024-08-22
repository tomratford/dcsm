dat <- sim_exp()

knots <- expand.grid(1:3,1:3,1:3)

fit_royston <- function(i) {
  v <- try(royston_parmar.fit(dat,knots[i,][[1]],knots[i,][[2]],knots[i,][[3]])$value)
  if (class(v) == "try-error") {
    v <- NA
  }
  v
}

parallel::mclapply(1:27, mc.cores = 16)
