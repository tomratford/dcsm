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
  joly_fit <- try(joly.fit(dat,5,5,5,10000,10000,10000,compute_cross = F),silent=T)
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

set.seed(123)
joly_exp_dist_sims <- parallel::mclapply(1:100, run_exp_sim, mc.cores=16)

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
run_nonlin_sim <- function(unused) {
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

  # penalised spline model
  joly_fit <- try(joly.fit(dat,5,5,5,100,100,100,compute_cross = F),silent=T)
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
  cat("Completed Nonlinear Sim #")
  cat(unused)
  cat("\n")

  c("JolyMSE" = mean(joly_error^2),
    "JolyTheta01" = joly_fit$par[["theta01"]],
    "JolyTheta02" = joly_fit$par[["theta02"]],
    "JolyTheta12" = joly_fit$par[["theta12"]])
}

set.seed(123)
joly_nonlin_sims <- parallel::mclapply(1:100, run_nonlin_sim, mc.cores=16)

save(list=c("joly_exp_dist_sims", "joly_nonlin_sims"), file=file.path(Sys.getenv("HOME"),"adhoc_joly.rda"))

