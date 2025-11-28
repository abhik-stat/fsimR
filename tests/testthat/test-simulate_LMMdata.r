# library(testthat)

# context("Testing simulate_LMMdata()")


test_that("simulate_LMMdata basic output structure is valid", {
  set.seed(123)
  sim <- simulate_LMMdata(
    sim_settings = sim_spec(
      n_subj = 5,
      n_obs = 4,
      beta_coeff = c(1, 2)
    ),
    distr_settings = list(
      RE_distr = list(distr_name = "mvnorm", distr_params = list(sigma = diag(1))),
      error_distr = list(distr_name = "norm", distr_params = list(sd = 1))
    )
  )

  expect_s3_class(sim, "LMMdata")
  expect_true(all(c("y", "X", "Z", "RE") %in% names(sim)))
  expect_equal(length(sim$y), 20)
  expect_equal(nlevels(sim$ID), sim$sim_settings$n_subj)
})


test_that("simulate_LMMdata handles random slopes and intercepts", {
  set.seed(42)
  Sigma_RE <- matrix(c(1, 0.2, 0.2, 1), 2)
  sim <- simulate_LMMdata(
    sim_settings = sim_spec(
      n_subj = 10,
      n_obs = 6,
      beta_coeff = c(1, -0.5, 0.3),
      include.Zintercept = TRUE
    ),
    distr_settings = list(
      RE_distr = list(distr_name = "mvnorm", distr_params = list(sigma = Sigma_RE)),
      error_distr = list(distr_name = "norm", distr_params = list(sd = 1))
    )
  )

  expect_equal(ncol(sim$Z), ncol(sim$distr_settings$RE_distr$distr_params$sigma))
  expect_true(is.matrix(sim$y))
})


test_that("simulate_LMMdata works with correlated X and Z", {
  set.seed(77)
  Sigma_X <- matrix(c(1, 0.5, 0.5, 1), 2)
  Sigma_Z <- diag(2)
  Sigma_XZ <- matrix(c(0.3, 0.1, 0.2, 0.4), 2)
  Sigma <- rbind(
    cbind(Sigma_X, Sigma_XZ),
    cbind(t(Sigma_XZ), Sigma_Z)
  )

  sim <- simulate_LMMdata(
    sim_settings = sim_spec(
      n_subj = 10,
      n_obs = 5,
      beta_coeff = c(1, -0.5, 0.3),
      is.CorrelatedXZ = TRUE
    ),
    distr_settings = list(
      X_distr  = list(distr_name = "mvnorm", distr_params = list(sigma = Sigma_X)),
      Z_distr  = list(distr_name = "mvnorm", distr_params = list(sigma = Sigma_Z)),
      XZ_distr = list(distr_name = "mvnorm", distr_params = list(sigma = Sigma)),
      RE_distr = list(distr_name = "mvnorm", distr_params = list(sigma = diag(3))),
      error_distr = list(distr_name = "norm", distr_params = list(sd = 1))
    )
  )

  expect_true(is.matrix(sim$X))
  expect_true(is.matrix(sim$Z))
  expect_true(all(dim(sim$Z) > 0))

  # Check empirical covariance approximately matches Sigma
  emp_cov <- cov(cbind(sim$X[,-1], sim$Z[,-1]))
  expect_true(max(abs(emp_cov - Sigma)) < 0.5)
})


test_that("simulate_LMMdata supports non-Gaussian error distributions", {
  set.seed(99)
  sim <- simulate_LMMdata(
    sim_settings = sim_spec(
      n_subj = 5,
      n_obs = 5,
      beta_coeff = c(1, 1)
    ),
    distr_settings = list(
      RE_distr = list(distr_name = "mvnorm", distr_params = list(sigma = diag(1))),
      error_distr = list(distr_name = "t", distr_params = list(df = 5))
    )
  )

  expect_true(is.numeric(sim$y))
  expect_false(any(is.na(sim$y)))
  # Check for heavier tails than normal
  skip_if_not_installed("e1071")
  expect_true(abs(e1071::skewness(as.vector(sim$y))) < 2)
})


test_that("simulate_LMMdata supports multiple replications", {
  set.seed(2025)
  sim <- simulate_LMMdata(
    n_rep = 3,
    sim_settings = sim_spec(
      n_subj = 8,
      n_obs = 4,
      beta_coeff = c(1, 2)
    ),
    distr_settings = list(
      RE_distr = list(distr_name = "mvnorm", distr_params = list(sigma = diag(1))),
      error_distr = list(distr_name = "norm", distr_params = list(sd = 1))
    )
  )

  expect_type(sim$y, "list")
  expect_length(sim$y, 3)
  expect_equal(length(sim$RE), 3)
  # Check consistency of y dimensions
  for(rep_y in sim$y) {
    expect_equal(dim(rep_y), c(8*4, 1))
  }
})


test_that("simulate_LMMdata multiple replications with iterations have correct dimensions", {
  n_rep <- c(2, 3) # 2 independent reps × 3 iterations each
  sim <- simulate_LMMdata(
    n_rep = n_rep,
    sim_settings = sim_spec(
      n_subj = 5,
      n_obs = 4,
      beta_coeff = c(1, 2)
    ),
    distr_settings = list(
      RE_distr = list(distr_name = "mvnorm", distr_params = list(sigma = diag(2))),
      error_distr = list(distr_name = "norm", distr_params = list(sd = 1))
    )
  )

  expect_type(sim$y, "list")
  expect_length(sim$y, n_rep[1])

  for(rep_RE in sim$RE) {
    expect_true(is.list(rep_RE))
    expect_length(rep_RE, n_rep[2])
    for(iter_RE in rep_RE) {
      expect_equal(dim(iter_RE), c(5,2))
    }
  }

  for(rep_y in sim$y) {
    expect_true(is.matrix(rep_y))
    expect_equal(ncol(rep_y), n_rep[2])
    expect_equal(nrow(rep_y), 5*4)
  }

})



test_that("simulate_LMMdata works with user-supplied design matrices", {
  n_subj <- 5
  n_obs <- 6
  n <- n_subj * n_obs
  X_user <- matrix(rnorm(n * 3), ncol = 3)
  Z_user <- matrix(rnorm(n * 2), ncol = 2)

  sim <- simulate_LMMdata(sim_settings = sim_spec(
    n_subj = n_subj,
    n_obs = n_obs,
    X = X_user,
    Z = Z_user,
    beta_coeff = c(1, 0.5, -1)
  ))

  expect_equal(dim(sim$X), dim(X_user))
  expect_equal(dim(sim$Z), dim(Z_user))
})


test_that("simulate_LMMdata orthogonalizes X and Z when requested", {
  sim <- simulate_LMMdata(sim_settings = sim_spec(
    p = 5, q=2,
    include.Xintercept = FALSE,
    include.Zintercept = FALSE,
    orthogonalize.X = TRUE,
    orthogonalize.Z = TRUE,
    n_subj = 120,
    n_obs = 50
  ))
  cX=cor(sim$X)
  cZ=cor(sim$Z)
  expect_true(max(abs(cX[upper.tri(cX)])) < 1e-3)
  expect_true(max(abs(cZ[upper.tri(cZ)])) < 1e-3)
})


test_that("simulate_LMMdata handles copula-based X", {
  if (requireNamespace("copula", quietly = TRUE))
    skip("Required package 'copula' is not installed.")
  normal_cop <- copula::normalCopula(param = 0.6, dim = 2)
  distr_settings <- list(
    X_distr = list(
      distr_name = "copula",
      distr_params = distr_spec(
        copula = normal_cop,
        margins = list(list(dist="norm", params = list(mean = 0, sd = 1)),
                       list(dist="norm", params = list(mean = 0, sd = 2)))
      )
    )
  )

  sim <- simulate_LMMdata(
    sim_settings = sim_spec(n_subj = 10, n_obs = 3),
    distr_settings = distr_settings
  )

  expect_equal(ncol(sim$X), 3)
  # Check approximate correlation matches copula parameter
  emp_cor <- cor(sim$X[,-1])[1,2]
  expect_true(abs(emp_cor - 0.6) < 0.2)
})


test_that("simulate_LMMdata random effects match custom covariance", {
  Sigma_RE <- matrix(c(1, 0.3, 0.3,
                       0.3, 1, 0.2,
                       0.3, 0.2, 1), ncol = 3)
  RE_dist=list( distr_params = list(sigma = Sigma_RE))
  sim <- simulate_LMMdata(sim_settings = sim_spec(n_subj = 150, n_obs = 40),
                          distr_settings = distr_spec(RE_distr = RE_dist))

  emp_cov <- cov(sim$RE)
  expect_true(max(abs(emp_cov - Sigma_RE)) < 0.3)
})


test_that("simulate_LMMdata respects specified signal-to-noise ratio", {
  SNR_target <- 3
  RE_dist=list( distr_params = list(sigma = diag(2)))
  sim <- simulate_LMMdata(sim_settings = sim_spec(
    n_subj = 10,
    n_obs = 5,
    beta_coeff = c(1, 1),
    SNR = SNR_target),
    distr_settings = distr_spec(RE_distr = RE_dist)
    )

  signal_var <- var(as.vector(sim$X %*% sim$sim_settings$beta_coeff +
                                rowSums(sim$Z * sim$RE[sim$ID, ,drop = FALSE])))
  noise_var  <- var(as.vector(sim$y - (sim$X %*% sim$sim_settings$beta_coeff +
                                         rowSums(sim$Z * sim$RE[sim$ID, ,drop = FALSE]))))
  expect_true(abs(signal_var / noise_var - SNR_target) < 1)
})



test_that("simulate_LMMdata makes y approximately equals X*beta + Z*RE", {
  RE_dist=list( distr_params = list(sigma = diag(2)))
  sim <- simulate_LMMdata(sim_settings = sim_spec(
    n_subj = 80,
    n_obs = 40,
    beta_coeff = c(1, 2),
    SNR = 10),
    distr_settings = distr_spec(RE_distr = RE_dist)
  )

  y_hat <- sim$X %*% sim$sim_settings$beta_coeff +
    rowSums(sim$Z * sim$RE[sim$ID, ,drop = FALSE])
  expect_true(cor(as.vector(sim$y), as.vector(y_hat)) > 0.95)
})



test_that("SNR control fails for dependent errors", {
  n_subj <- 3
  n_obs <- rep(2, n_subj)
  n_total <- sum(n_obs)

  Sigma_dep <- matrix(0.5, nrow = n_total, ncol = n_total)
  diag(Sigma_dep) <- 1

  distr_settings <- list(
    error_distr = list(
      distr_name = "mvnorm",
      distr_params = list(sd = Sigma_dep)
    )
  )

  sim_settings <- list(n_subj = n_subj, n_obs = n_obs, SNR = 5)

  expect_error(
    simulate_LMMdata(n_rep = 1, sim_settings = sim_settings, distr_settings = distr_settings),
    "SNR control is allowed only for IID error distribution!"
  )
})




test_that("simulate_LMMdata works with dependent errors", {
  n_subj <- 5
  n_obs <- rep(3, n_subj)
  n_total <- sum(n_obs)

  Sigma_dep <- matrix(0.5, nrow = n_total, ncol = n_total)
  diag(Sigma_dep) <- 1

  distr_settings <- list(
    error_distr = list(
      distr_name = "mvnorm",
      distr_params = list(sd = Sigma_dep)
    )
  )

  sim_settings <- list(
    n_subj = n_subj,
    n_obs = n_obs,
    beta_coeff = 0,
    X = matrix(1, n_total, 1),
    Z = matrix(0, n_total, 1)
  )

  set.seed(123)
  sim <- simulate_LMMdata(n_rep = c(1, 1000), sim_settings = sim_settings, distr_settings = distr_settings)

  y <- sim$y
  emp_cov <- cov(t(y))
  max_diff <- max(abs(emp_cov - Sigma_dep))
  expect_true(max_diff < 0.2)

  expect_s3_class(sim, "LMMdata")
  expect_equal(nrow(sim$y), sum(n_obs))
})

test_that("simulate_LMMdata works with heterogeneous errors", {
  n_subj <- 4
  n_obs <- rep(2, n_subj)
  n_total <- sum(n_obs)

  hetero_sd <- diag(seq(0.5, 2, length.out = n_total))

  distr_settings <- list(
    error_distr = list(
      distr_name = "mvnorm",
      distr_params = list(sd = hetero_sd)
    )
  )

  sim_settings <- list(
    n_subj = n_subj,
    n_obs = n_obs,
    beta_coeff = 0,
    X = matrix(1, n_total, 1),
    Z = matrix(0, n_total, 1)
  )

  set.seed(456)
  sim <- simulate_LMMdata(n_rep = c(1, 1000), sim_settings = sim_settings, distr_settings = distr_settings)

  y <- sim$y
  emp_var <- diag(cov(t(y)))
  max_diff <- max(abs(emp_var - diag(hetero_sd)))
  expect_true(max_diff < 0.25)

  expect_s3_class(sim, "LMMdata")
  expect_equal(nrow(sim$y), sum(n_obs))
})


