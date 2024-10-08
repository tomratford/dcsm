# Simulate 100 exponential datasets, apply 3 model fitting procedures across
# all 100 and store the results in a data frame that is stored in the users home
# directory

devtools::load_all()

# true values
lambda01 <- 4
lambda02 <- 1
lambda12 <- 0.5
theta01 <- log(0.4)
theta02 <- log(0.4)
theta12 <- log(1)

# true survivor function
logS_t <- \(t,z) {
  -(lambda01*exp(z*theta01) + lambda02*exp(z*theta02))*t
}

set.seed(123)
run_simulation <- function(unused) {
  dat <- sim_exp(N = 300,
                 lambda01 = lambda01,
                 lambda02 = lambda02,
                 lambda12 = lambda12,
                 theta01 = theta01,
                 theta02 = theta02,
                 theta12 = theta12,
                 K = 8) # number of visits

  # cox model, for comparison
  cox_mod <- coxph(Surv(PFSDY, PFS) ~ ATRTN, data=dat)
  cox_error <- exp(logS_t(dat$PFSDY, dat$ATRTN)) - predict(cox_mod, type="survival")

  # weibull model
  weib_fit <- weib.fit(dat)
  weib_fns <- do.call(weib.fnBuilder,as.list(weib_fit$par))
  weib_error <-  exp(logS_t(dat$PFSDY, dat$ATRTN)) - weib_fns$P00(0, dat$PFSDY, dat$ATRTN)

  # penalised spline model
  joly_fit <- try(joly.fit(dat,0,0,0,100,100,100,compute_cross = F),silent=T)
  if (class(joly_fit) == "try-error") {
    joly_fit <- list(
      par = c("theta01" = NA, "theta02" = NA, "theta12" = NA)
    )
    joly_error <- NA
  } else {
    joly_fns <- do.call(joly.fnBuilder, make_pars2(joly_fit$par, joly.initials(dat,0,0,0)))
    joly_error <- exp(logS_t(dat$PFSDY, dat$ATRTN)) - joly_fns$P00(rep(0,length(dat$PFSDY)), dat$PFSDY, dat$ATRTN)
  }

  # royston parmar model
  royston_fit <- try(royston_parmar.fit(dat,0,0,2))
  if (class(royston_fit) == "try-error") {
    royston_fit <- list(
      par = c("theta01" = NA, "theta02" = NA, "theta12" = NA)
    )
    royston_error <- NA
  } else {
    royston_fns <- do.call(royston_parmar.fnBuilder, make_pars2(royston_fit$par, royston_parmar.initials(dat,0,0,2)))
    royston_error <- exp(logS_t(dat$PFSDY, dat$ATRTN)) - royston_fns$P00(rep(0,length(dat$PFSDY)), dat$PFSDY, dat$ATRTN)
  }

  #basic debugging
  cat("Completed Sim #")
  cat(unused)
  cat("\n")

  return(
    c(
      "CoxMSE" = mean(cox_error ^ 2),
      "CoxTheta" = unname(cox_mod$coef),
      "WeibMSE" = mean(weib_error ^ 2),
      "WeibTheta01" = weib_fit$par[["theta01"]],
      "WeibTheta02" = weib_fit$par[["theta02"]],
      "WeibTheta12" = weib_fit$par[["theta12"]],
      "JolyMSE" = mean(joly_error ^ 2),
      "JolyTheta01" = joly_fit$par[["theta01"]],
      "JolyTheta02" = joly_fit$par[["theta02"]],
      "JolyTheta12" = joly_fit$par[["theta12"]],
      "RoystonMSE" = mean(royston_error ^ 2),
      "RoystonTheta01" = royston_fit$par[["theta01"]],
      "RoystonTheta02" = royston_fit$par[["theta02"]],
      "RoystonTheta12" = royston_fit$par[["theta12"]]
    )
  )
}

exp_dist_sims <- parallel::mclapply(1:100, run_simulation, mc.cores=13)
save(list="exp_dist_sims", file=file.path(Sys.getenv("HOME"),"exp_dist_sims.rda"))
