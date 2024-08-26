devtools::load_all()

fit_model <- function(dat) {
  kappas <- c(100,1000,10000,100000,1000000)
  mods <- parallel::mclapply(kappas, \(k) {
    v <- try(joly.fit(dat, 5,5,5, k,k,k), silent = T)
    if (class(v) == 'try-error') {
      v <- list(cross_value = -Inf)
    }
    v
  }, mc.cores = 16)
  scores <- sapply(mods, \(m) m$cross_value)
  if (all(scores == -Inf)) {
    stop("All model fitting failed")
  }
  i <- which.max(scores)
  c(
    "final_mod" = mods[i],
    "scores" = scores
  )
}

joly_adhoc_50203 <- fit_model(AMGEN_20050203)
joly_adhoc_20408 <- fit_model(AMGEN_20020408)
