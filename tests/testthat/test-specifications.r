# library(testthat)

# context("Testing distr_spec, sim_spec, and match.fun.allR")


# ========== Tests for distr_spec ===================================
test_that("distr_spec default output", {
  out <- distr_spec()
  expect_type(out, "list")
  expect_equal(out$dim, 1)
  expect_equal(out$mean, 0)
  expect_equal(as.numeric(out$sd), 1)
  expect_equal(out$sd, as.matrix(1))
  expect_equal(out$alpha, 0)
})

test_that("distr_spec handles mean vector correctly", {
  out <- distr_spec(mean = c(1, 2, 3))
  expect_equal(out$mean, c(1, 2, 3))
  expect_equal(out$dim, 3)
})

test_that("distr_spec handles explicit dim with scalar mean", {
  out <- distr_spec(dim = 4, mean = 5)
  expect_equal(out$mean, rep(5, 4))
  expect_equal(out$dim, 4)
})

test_that("distr_spec resolves location aliases", {
  out1 <- distr_spec(mean = 5)
  out2 <- distr_spec(mu = 5)
  out3 <- distr_spec(location = 5)
  out4 <- distr_spec(xi = 5)
  expect_equal(out1$mean, 5)
  expect_equal(out2$mean, 5)
  expect_equal(out3$mean, 5)
  expect_equal(out4$mean, 5)
})


test_that("distr_spec scales scalar and matrix correctly", {
  out_scalar <- distr_spec(sd = 2, dim = 3)
  expect_equal(out_scalar$sd, diag(3) * 2)
  expect_equal(out_scalar$dim, 3)

  mat <- diag(3)
  out_matrix <- distr_spec(sigma = mat)
  expect_equal(out_matrix$sigma, mat)
})

test_that("distr_spec errors for incompatible mean and sd lengths", {
  expect_error(distr_spec(mean = c(1,2), sd = diag(3)))
})

test_that("distr_spec handles incompatible mu and location lengths as per preceedence", {
  out <- (distr_spec(mu = c(1,2), location = c(9,9,9)))
  expect_equal(out$dim, 2)
  expect_equal(out$location, out$mean)
  expect_equal(out$location, out$mu)
})

test_that("distr_spec handles mean as matrix or non-square sd matrix", {
  expect_error(distr_spec(mean = matrix(1:4, 2, 2)))
  expect_error(distr_spec(sd = matrix(1:6, 2, 3)))
})

test_that("distr_spec handles column vector mean matrix", {
  mat_mean <- matrix(1:4, 4, 1)
  out <- distr_spec(mean = mat_mean)
  expect_equal(out$mean, drop(mat_mean))
})

test_that("distr_spec does alpha vector based dimension inference", {
  out <- distr_spec(alpha = c(0.1, 0.2, 0.3))
  expect_equal(length(out$alpha), 3)
  expect_equal(out$dim, 3)
})

test_that("distr_spec adds extra arguments", {
  out <- distr_spec(extra1 = "hello", extra2 = 42)
  expect_equal(out$extra1, "hello")
  expect_equal(out$extra2, 42)
})

test_that("distr_spec throws error on incompatible dimensions", {
  expect_error(distr_spec(mean = 1:2, alpha = 1:3))
})


# ============ Tests for sim_spec =========================
test_that("sim_spec default output", {
  out <- sim_spec()
  expect_type(out, "list")
  expect_equal(out$n_subj, 10)
  expect_equal(out$n_obs, 5)
  expect_equal(out$p, 1)
  expect_equal(out$q, 1)
})

test_that("sim_spec infers p from beta_coeff", {
  out <- sim_spec(beta_coeff = c(1, 2, 3))
  expect_equal(out$p, 3)
  expect_equal(out$beta_coeff, c(1, 2, 3))
})

test_that("sim_spec stores user-supplied X and Z", {
  X <- matrix(1:6, 2, 3)
  Z <- matrix(1:4, 2, 2)
  out <- sim_spec(X = X, Z = Z)
  expect_equal(out$X, X)
  expect_equal(out$Z, Z)
})

test_that("sim_spec accepts additional arguments", {
  out <- sim_spec(custom_flag = TRUE)
  expect_true(out$custom_flag)
})

# ============= Tests for match.fun.allR ===========================
test_that("match.fun.allR returns function unchanged if input is function", {
  f <- function(x) x + 1
  expect_identical(match.fun.allR(f), f)
})

test_that("match.fun.allR finds base R function by name", {
  f <- match.fun.allR("rnorm")
  expect_true(is.function(f))
})

test_that("match.fun.allR throws error if function not found", {
  expect_error(match.fun.allR("nonexistent_fun"))
})

test_that("match.fun.allR can find explicit package function", {
  skip_if_not_installed("stats")
  f <- match.fun.allR("stats::rnorm")
  expect_true(is.function(f))
})

test_that("match.fun.allR with `list_all = TRUE` returns list or single function", {
  matches <- match.fun.allR("rnorm", list_all = TRUE)
  expect_true(is.list(matches) || is.function(matches))
  if (is.list(matches)) {
    expect_true(all(sapply(matches, is.function)))
  } else {
    expect_true(is.function(matches))
  }
})

test_that("match.fun.allR behaves identically to base match.fun for standard cases", {
  # Basic function name
  expect_identical(match.fun("mean"), match.fun.allR("mean"))
  expect_identical(match.fun("sum"), match.fun.allR("sum"))

  # Function passed directly
  f <- function(x) x + 1
  expect_identical(match.fun(f), match.fun.allR(f))

  # Function name as character string
  expect_identical(match.fun("sum"), match.fun.allR("sum"))

  # Function name as a symbol
  expect_identical(match.fun(quote(mean)), match.fun.allR(quote(mean)))
})


test_that("match.fun.allR and match.fun both error for non-existent functions", {
  err_base <- tryCatch(match.fun("non_existing_function"), error = function(e) e$message)
  err_custom <- tryCatch(match.fun.allR("non_existing_function"), error = function(e) e$message)
  expect_identical(class(err_base), class(err_custom))
})

test_that("match.fun.allR searches environments like match.fun", {
  env <- new.env()
  env$my_local_fun <- function(x) x * 10

  # Not in global env
  err_base <- tryCatch(match.fun("my_local_fun"), error = function(e) e$message)
  err_custom <- tryCatch(match.fun.allR("my_local_fun"), error = function(e) e$message)
  expect_identical(class(err_base), class(err_custom))

  # Assign to global env and retry
  assign("my_local_fun", env$my_local_fun, envir = .GlobalEnv)
  expect_identical(match.fun("my_local_fun"), match.fun.allR("my_local_fun"))

  # Clean up
  rm("my_local_fun", envir = .GlobalEnv)
})

test_that("match.fun.allR can find functions in installed packages", {
  skip_if_not_installed("dplyr")

  # 'filter' exists in dplyr
  err_base <- tryCatch(match.fun("filter"), error = function(e) e$message)
  err_custom <- tryCatch(match.fun.allR("filter"), error = function(e) e$message)

  # Either both error (if not loaded) or both succeed
  expect_identical(err_base, err_custom)
})


test_that("match.fun.allR errors for invalid input types", {
  err_base <- tryCatch(match.fun(list()), error = function(e) e$message)
  err_custom <- tryCatch(match.fun.allR(list()), error = function(e) e$message)
  expect_identical(err_base, err_custom)
})


