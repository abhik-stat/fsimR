# library(testthat)

# context("Testing simulate_LMdata")


test_that("simulate.LMdata() basic simulation works", {
  set.seed(123)
  sim1 <- simulate_LMdata(sim_settings = sim_spec(
    n_subj = 20,
    beta_coeff = c(2, 1.5),
    sigma_e = 1
  ))

  expect_type(sim1, "list")
  expect_true(all(c("y", "X") %in% names(sim1)))
  expect_s3_class(sim1, "LMdata")

  expect_equal(nrow(sim1$X), 20)
  expect_equal(ncol(sim1$X), length(sim1$sim_settings$beta_coeff))
  expect_equal(length(sim1$y), 20)
})




test_that("simulate.LMdata() handles multiple replications correctly", {
  set.seed(42)
  sim2 <- simulate_LMdata(n_rep = 3, sim_settings = sim_spec(
    n_subj = 10,
    beta_coeff = c(1, -0.5)
  ))

  expect_length(sim2$X, 3)
  expect_length(sim2$y, 3)
  expect_true(all(sapply(sim2$X, nrow) == 10))

  expect_type(sim2$y, "list")
  expect_length(sim2$y, 3)
  # Check consistency of y dimensions
  for(rep_y in sim2$y) {
    expect_equal(dim(rep_y), c(10, 1))
  }
})



test_that("simulate.LMdata() adjusts sigma_e to match target SNR", {
  set.seed(99)
  sim <- simulate_LMdata(sim_settings = sim_spec(
    n_subj = 3000,
    beta_coeff = c(2, 1),
    SNR = 5
  ))

  snr_obs <- mean(unlist(sim$SNR))
  expect_true(abs(snr_obs - 5) < 0.5) # within tolerance

  y_hat <- sim$X %*% sim$sim_settings$beta_coeff
  signal_var <- var(as.vector(y_hat))
  noise_var  <- var(as.vector(sim$y - y_hat))
  expect_true(abs(signal_var / noise_var - 5) < 1)

  expect_true(cor(as.vector(sim$y), as.vector(y_hat)) > 0.9)


})



test_that("simulate_LMdata multiple replications with iterations have correct dimensions", {
  n_rep <- c(2, 3) # 2 independent reps × 3 iterations each
  sim <- simulate_LMdata(
    n_rep = n_rep,
    sim_settings = sim_spec(
      n_subj = 10,
      beta_coeff = c(1, 2)
    ),
    distr_settings = list(
      error_distr = list(distr_name = "norm", distr_params = list(sd = 1))
    )
  )

  expect_type(sim$y, "list")
  expect_length(sim$y, n_rep[1])


  for(rep_y in sim$y) {
    expect_true(is.matrix(rep_y))
    expect_equal(ncol(rep_y), n_rep[2])
    expect_equal(nrow(rep_y), 10)
  }

})




test_that("simulate.LMdata() throws error for bad inputs", {

  # --- 1. Negative sample size ---
  expect_error(
    simulate_LMdata(sim_settings = sim_spec(n_subj = -10)),
    regexp = "positive"
  )

  # --- 2. User-supplied X with wrong number of columns ---
  badX <- matrix(rnorm(10 * 5), ncol = 5)  # has 5 columns
  # but beta_coeff has length 3
  expect_error(
    simulate_LMdata(
      sim_settings = sim_spec(
        X = badX,
        beta_coeff = c(1, 2, 3)
      )
    ),
    regexp = "beta_coeff|X columns",
    ignore.case = TRUE
  )
})


test_that("simulate.LMdata() works with user-supplied X", {
  set.seed(55)
  X_custom <- cbind(1, matrix(rnorm(20 * 2), ncol = 2))
  sim5 <- simulate_LMdata(sim_settings = sim_spec(
    X = X_custom,
    beta_coeff = c(2, 1, -1),
    sigma_e = 0.5
  ))

  expect_equal(nrow(sim5$X), nrow(X_custom))
  expect_equal(unname(sim5$X), unname(as.matrix(X_custom)))
})

test_that("simulate.LMdata() supports orthogonalized predictors", {
  set.seed(10)
  sim4 <- simulate_LMdata(sim_settings = sim_spec(
    n_subj = 15,
    beta_coeff = c(2, -1, 0.5),
    orthogonalize.X = TRUE
  ))

  X <- sim4$X[,-1]
  cor_vals <- cor(X)
  off_diag <- cor_vals[upper.tri(cor_vals)]
  expect_true(mean(abs(off_diag)) < 0.1)
})



test_that("simulate.LMdata() supports non-normal errors", {
  set.seed(111)
  sim6 <- simulate_LMdata(sim_settings = sim_spec(
    n_subj = 20,
    beta_coeff = c(1, 2)
  ),
  distr_settings = list(error_distr = list(distr_name = "t", distr_params = list(df = 3)))
  )

  expect_equal(nrow(sim6$X), 20)
  expect_equal(length(sim6$y), 20)
  expect_true(is.numeric(sim6$y))
  expect_false(any(is.na(sim6$y)))
  skip_if_not_installed("e1071")
  expect_true(abs(e1071::skewness(as.vector(sim6$y))) < 2)
})




test_that("simulate.LMdata() supports custom generator", {
  set.seed(222)
  my_gen <- function(n, min = 0, max = 1) runif(n, min, max)

  sim7 <- simulate_LMdata(
    sim_settings = sim_spec(n_subj = 20, beta_coeff = c(1, 2)),
    distr_settings = list(
      X_distr = list(
        distr_name = "custom",
        distr_params = list(min = 0, max = 1),
        generator = my_gen
      ),
      error_distr = list(
        distr_name = "custom",
        distr_params = list(min = 0, max = 1),
        generator = my_gen
      )
    )
  )

  expect_equal(nrow(sim7$X), 20)
  expect_true(all(sim7$y >= 0 & sim7$y <= 5))
})


test_that("simulate.LMdata() supports copula-based X generation", {
  skip_if_not_installed("copula")
  library(copula)

  cop <- normalCopula(param = 0.5, dim = 3)
  margins <- list(
    list(dist = "norm", params = list(mean = 0, sd = 1)),
    list(dist = "norm", params = list(mean = 2, sd = 0.5)),
    list(dist = "unif", params = list(min = -1, max = 1))
  )

  sim8 <- simulate_LMdata(
    sim_settings = sim_spec(n_subj = 250, beta_coeff = c(1, 2, -1, 0.5)),
    distr_settings = list(
      X_distr = list(
        distr_name = "copula",
        distr_params = list(dim = 3, copula = cop, margins = margins)
      )
    )
  )

  expect_equal(nrow(sim8$X), 250)
  expect_true(all(is.finite(sim8$X)))

  expect_equal(ncol(sim8$X), 4)
  # Check approximate correlation matches copula parameter
  emp_cor <- cor(sim8$X[,-1])[1,2]
  expect_true(abs(emp_cor - 0.5) < 0.2)

})





test_that("SNR control fails for dependent errors", {
  n_subj <- 30

  Sigma_dep <- matrix(0.5, nrow = n_subj, ncol = n_subj)
  diag(Sigma_dep) <- 1

  distr_settings <- list(
    error_distr = list(
      distr_name = "mvnorm",
      distr_params = list(sd = Sigma_dep)
    )
  )

  sim_settings <- list(n_subj = n_subj, SNR = 5)

  expect_error(
    simulate_LMdata(n_rep = 1, sim_settings = sim_settings, distr_settings = distr_settings),
    "SNR control is allowed only for IID error distribution!"
  )
})




test_that("simulate_LMdata works with dependent errors", {
  n_subj <- 50

  Sigma_dep <- matrix(0.5, nrow = n_subj, ncol = n_subj)
  diag(Sigma_dep) <- 1

  distr_settings <- list(
    error_distr = list(
      distr_name = "mvnorm",
      distr_params = list(sd = Sigma_dep)
    )
  )

  sim_settings <- list(
    n_subj = n_subj,
    beta_coeff = 0,
    X = matrix(1, n_subj, 1))

  set.seed(123)
  sim <- simulate_LMdata(n_rep = c(1, 1000), sim_settings = sim_settings, distr_settings = distr_settings)

  y <- sim$y
  emp_cov <- cov(t(y))
  max_diff <- max(abs(emp_cov - Sigma_dep))
  expect_true(max_diff < 0.2)

  expect_s3_class(sim, "LMdata")
  expect_equal(nrow(sim$y), n_subj)
})

test_that("simulate_LMdata works with heterogeneous errors", {
  n_subj <- 40

  hetero_sd <- diag(seq(0.5, 2, length.out = n_subj))

  distr_settings <- list(
    error_distr = list(
      distr_name = "mvnorm",
      distr_params = list(sd = hetero_sd)
    )
  )

  sim_settings <- list(
    n_subj = n_subj,
    beta_coeff = 0,
    X = matrix(1, n_subj, 1)  )

  set.seed(456)
  sim <- simulate_LMdata(n_rep = c(1, 1000), sim_settings = sim_settings, distr_settings = distr_settings)

  y <- sim$y
  emp_var <- diag(cov(t(y)))
  max_diff <- max(abs(emp_var - diag(hetero_sd)))
  expect_true(max_diff < 0.25)

  expect_s3_class(sim, "LMdata")
  expect_equal(nrow(sim$y), n_subj)
})


