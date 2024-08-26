dat <- sim_nonlinear()

kappa <- c(100, 1000, 10000, 100000, 1000000)

parallel::mclapply(kappa, \(k) joly.fit(dat, 5,5,5, k,k,k), mc.cores = 5)
