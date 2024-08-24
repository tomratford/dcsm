# Simulate 100 datasets with non-linear log intensities, apply 3 model fitting
# procedures across all 100 and store the results in a data frame that is
# stored in the users home directory

devtools::load_all()

# true values
a01 <- \(x) x^(4/5)
A01 <- \(x) (5/9)*x^(9/5)
a02 <- \(x) 3*x/4
A02 <- \(x) (3*(x^2))/8
a12 <- \(x) (3*x/2)^(5/4)
A12 <- \(x) ((2/3)^(3/4))*(x^(9/4))
theta01 <- log(0.4)
theta02 <- log(0.4)
theta12 <- log(1)

# true survivor function
logS_t <- \(t,z) {
  -(A01(t)*exp(z*theta01) + A02(t)*exp(z*theta02))
}

set.seed(123)
run_simulation <- function(unused) {
  dat <- sim_nonlinear(N = 300,
                       a01 = a01,
                       A01 = A01,
                       a02 = a02,
                       A02 = A02,
                       a12 = a12,
                       A12 = A12,
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
  joly_fit <- try(joly.fit(dat,1,0,1,100,100,100,compute_cross = F),silent=T)
  if (class(joly_fit) == "try-error") {
    joly_fit <- list(
      par = c("theta01" = NA, "theta02" = NA, "theta12" = NA)
    )
    joly_error <- NA
  } else {
    joly_fns <- do.call(joly.fnBuilder, make_pars2(joly_fit$par, joly.initials(dat,1,0,1)))
    joly_error <- exp(logS_t(dat$PFSDY, dat$ATRTN)) - joly_fns$P00(rep(0,length(dat$PFSDY)), dat$PFSDY, dat$ATRTN)
  }

  # royston parmar model
  royston_fit <- try(royston_parmar.fit(dat,1,0,2), silent=T)
  if (class(royston_fit) == "try-error") {
    royston_fit <- list(
      par = c("theta01" = NA, "theta02" = NA, "theta12" = NA)
    )
    royston_error <- NA
  } else {
    royston_fns <- do.call(royston_parmar.fnBuilder, make_pars2(royston_fit$par, royston_parmar.initials(dat,1,0,2)))
    royston_error <- exp(logS_t(dat$PFSDY, dat$ATRTN)) - royston_fns$P00(rep(0,length(dat$PFSDY)), dat$PFSDY, dat$ATRTN)
  }

  #basic debugging
  cat("Completed Sim #")
  cat(unused)
  cat("\n")

  list("Cox" = cox_mod,
       "Weib" = weib_fit,
       "Royston" = royston_fit,
       "Joly" = joly_fit)
}

nonlinear_ex <- run_simulation()
save(list="nonlinear_ex", file=file.path(Sys.getenv("HOME"),"nonlinear_ex.rda"))
