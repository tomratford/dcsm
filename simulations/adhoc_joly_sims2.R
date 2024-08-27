devtools::load_all()

lambda01 <- 4
lambda02 <- 1
lambda12 <- 0.5
theta01 <- log(0.4)
theta02 <- log(0.4)
theta12 <- log(1)

logS_t <- \(t,z) {
  -(lambda01*exp(z*theta01) + lambda02*exp(z*theta02))*t
}

run_exp_sim <- function(unused) {
  dat <- sim_exp(N = 300,
                 lambda01 = lambda01,
                 lambda02 = lambda02,
                 lambda12 = lambda12,
                 theta01 = theta01,
                 theta02 = theta02,
                 theta12 = theta12,
                 K = 8) # number of visits

  # penalised spline model
  joly_fit <- try(joly.fit(dat,5,5,5,1000,1000,1000,compute_cross = F),silent=T)
  if (class(joly_fit) == "try-error") {
    joly_fit <- list(
      par = c("theta01" = NA, "theta02" = NA, "theta12" = NA)
    )
    joly_error <- NA
  } else {
    joly_fns <- do.call(joly.fnBuilder, make_pars2(joly_fit$par, joly.initials(dat,5,5,5)))
    joly_error <- exp(logS_t(dat$PFSDY, dat$ATRTN)) - joly_fns$P00(rep(0,length(dat$PFSDY)), dat$PFSDY, dat$ATRTN)
  }

  #basic debugging
  cat("Completed Exp Sim #")
  cat(unused)
  cat("\n")

  return(
    c(
      "JolyMSE" = mean(joly_error ^ 2),
      "JolyTheta01" = joly_fit$par[["theta01"]],
      "JolyTheta02" = joly_fit$par[["theta02"]],
      "JolyTheta12" = joly_fit$par[["theta12"]]
    )
  )
}

save(list=c("joly_exp_dist_sims"), file=file.path(Sys.getenv("HOME"),"adhoc_joly_1000.rda"))

