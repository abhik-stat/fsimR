# library(testthat)

# context("Testing add_contamination ")
#
# set.seed(123)

################################################################
# Updated helper function
################################################################

# Count contaminated elements/rows depending on type
count_diff <- function(original, contaminated) {

  # Vector case ----
  if (is.vector(original) && !is.list(original)) {
    return(sum(original != contaminated))
  }

  # Matrix case ----
  if (is.matrix(original)) {
    # count number of rows with at least one modified cell
    return(sum(apply(original != contaminated, 1, any)))
  }

  # List case (list of matrices or vectors) ----
  if (is.list(original)) {
    if (!is.list(contaminated))
      stop("Both original and contaminated must be lists")

    out <- mapply(
      count_diff,
      original,
      contaminated,
      SIMPLIFY = TRUE
    )
    return(out)
  }

  stop("Unsupported input type for count_diff()")
}

################################################################
# Original tests (unchanged)
################################################################

test_that("Additive contamination on numeric vector", {
  x <- rnorm(100)
  cont <- add_contamination(
    x,
    cont_pos = list(name="y", col=1),
    cont_settings = list(
      cont_prop = 0.1,
      cont_type = "additive",
      cont_distr = list(distr_name="norm", distr_params=list(mean=5, sd=1))
    )
  )
  expect_true(is.numeric(cont))
  expect_equal(length(cont), length(x))
  expect_gt(count_diff(x, cont), 0)  # Some rows must be modified
})

test_that("Multiplicative contamination on matrix", {
  mat <- matrix(rnorm(200), ncol=2)
  cont <- add_contamination(
    mat,
    cont_pos = list(name="y", col=1:2),
    cont_settings = list(
      cont_prop = 0.2,
      cont_type = "multiplicative",
      cont_value = 2
    )
  )
  expect_true(is.matrix(cont))
  expect_equal(dim(cont), dim(mat))
  expect_true(any(cont != mat))
  diff_rows <- which(mat != cont, arr.ind = TRUE)[,1]
  for (r in diff_rows) {
    for (c in 1:2) {
      if (mat[r, c] != cont[r, c]) {
        expect_equal(cont[r, c], mat[r, c]*2)
      }
    }
  }
})

test_that("Replace contamination with fixed values", {
  mat <- matrix(rnorm(100), ncol=2)
  cont <- add_contamination(
    mat,
    cont_pos = list(name="y"),
    cont_settings = list(
      cont_prop = 0.1,
      cont_type = "replace",
      cont_value = c(10, 20)
    )
  )
  expect_true(all(cont[cont != mat] %in% c(10,20)))
})

test_that("Outlier contamination generates extreme values", {
  mat <- matrix(rnorm(100, mean=0, sd=1), ncol=2)
  cont <- add_contamination(
    mat,
    cont_pos = list(name="y"),
    cont_settings = list(
      cont_prop = 0.1,
      cont_type = "outlier",
      outlier_factor = 5
    )
  )
  for (c in 1:2) {
    expect_true(any(cont[,c] > mean(mat[,c]) + 4*sd(mat[,c]) |
                      cont[,c] < mean(mat[,c]) - 4*sd(mat[,c])))
  }
})

test_that("Leverage contamination on list with X and y", {
  sim <- list(
    y = rnorm(100),
    X = cbind(1, rnorm(100), rnorm(100))
  )
  cont <- add_contamination(
    sim,
    cont_pos = list(name="X", col=2:3),
    cont_settings = list(
      cont_type = "leverage",
      cont_prop = 0.1,
      leverage_factor = 10,
      outlier_factor = 5
    )
  )
  expect_true(any(abs(cont$X[,2] - mean(sim$X[,2])) > 5*sd(sim$X[,2])))
  expect_true(any(abs(cont$X[,3] - mean(sim$X[,3])) > 5*sd(sim$X[,3])))
  row_idx <- which(cont$X[,2] != sim$X[,2] | cont$X[,3] != sim$X[,3])
  expect_true(all(sign(cont$y[row_idx] - mean(sim$y)) ==
                    -sign(rowSums(cont$X[row_idx,2:3] -
                                    colMeans(sim$X[,2:3])))))
})

test_that("cont_value overrides random generation", {
  x <- rnorm(50)
  cont <- add_contamination(
    x,
    cont_pos = list(name="y"),
    cont_settings = list(
      cont_prop = 0.2,
      cont_type = "additive",
      cont_value = 3
    )
  )
  expect_true(all((cont - x)[cont != x] == 3))
})

test_that("cont_prop=0 returns identical object", {
  x <- rnorm(50)
  cont <- add_contamination(
    x,
    cont_pos = list(name="y"),
    cont_settings = list(
      cont_prop = 0
    )
  )
  expect_equal(cont, x)
})

test_that("Multiple columns in y are all contaminated when col=NULL", {
  mat <- matrix(rnorm(100), ncol=4)
  cont <- add_contamination(
    mat,
    cont_pos = list(name="y", col=NULL),
    cont_settings = list(
      cont_prop = 0.1,
      cont_type = "replace",
      cont_value = 5
    )
  )
  expect_true(all(cont[cont != mat] == 5))
})

test_that("Error handling for invalid column index", {
  mat <- matrix(rnorm(10), ncol=2)
  expect_error(add_contamination(mat, cont_pos = list(name="y", col=5)))
})

test_that("Error handling for non-numeric data", {
  x <- letters[1:10]
  expect_error(add_contamination(x))
})

test_that("Casewise additive contamination modifies same rows", {
  mat <- matrix(rnorm(200), ncol=4)
  cont <- add_contamination(
    mat,
    cont_pos = list(name="y", col=1:4),
    cont_settings = list(
      cont_prop = 0.2,
      cont_type = "additive",
      cont_value = 5,
      cont_mode = "casewise"
    )
  )
  contaminated_rows <- lapply(1:4, function(c) which(mat[,c] != cont[,c]))
  for (k in 2:4) {
    expect_equal(contaminated_rows[[1]], contaminated_rows[[k]])
  }
})

test_that("Cellwise additive contamination modifies different rows", {
  mat <- matrix(rnorm(200), ncol=4)
  cont <- add_contamination(
    mat,
    cont_pos = list(name="y", col=1:4),
    cont_settings = list(
      cont_prop = 0.1,
      cont_type = "additive",
      cont_value = 3,
      cont_mode = "cellwise"
    )
  )
  sets <- lapply(1:4, function(c) which(mat[,c] != cont[,c]))
  expect_false(identical(sets[[1]], sets[[2]]))
})


test_that("Cellwise multiplicative contamination applies independently", {
  mat <- matrix(rnorm(150), ncol=3)
  cont <- add_contamination(
    mat,
    cont_pos = list(name="y", col=1:3),
    cont_settings = list(
      cont_prop = 0.1,
      cont_value = 4,
      cont_type = "multiplicative",
      cont_mode = "cellwise"
    )
  )
  for (c in 1:3) {
    changed <- which(mat[,c] != cont[,c])
    if (length(changed) > 0) {
      expect_true(all(cont[changed,c] == mat[changed,c] * 4))
    }
  }
})

test_that("Casewise leverage contamination shifts same rows in X", {
  sim <- list(
    y = rnorm(80),
    X = cbind(1, rnorm(80), rnorm(80))
  )
  cont <- add_contamination(
    sim,
    cont_pos = list(name="X", col=2:3),
    cont_settings = list(
      cont_prop = 0.25,
      cont_type = "leverage",
      leverage_factor = 12,
      outlier_factor = 4,
      cont_mode = "casewise"
    )
  )
  rows_2 <- which(sim$X[,2] != cont$X[,2])
  rows_3 <- which(sim$X[,3] != cont$X[,3])
  expect_equal(rows_2, rows_3)
})



test_that("Casewise replace contamination affects same rows", {
  nrows <- 34  # divisible by 3 columns
  ncols <- 3
  mat <- matrix(rnorm(nrows * ncols), nrow = nrows, ncol = ncols)

  cont <- add_contamination(
    mat,
    cont_pos = list(name="y", col=1:ncols),
    cont_settings = list(
      cont_prop = 0.3,
      cont_type = "replace",
      cont_value = 99,
      cont_mode = "casewise"
    )
  )

  contam_rows <- lapply(1:ncols, function(c) which(mat[,c] != cont[,c]))
  expect_equal(contam_rows[[1]], contam_rows[[2]])
  expect_equal(contam_rows[[1]], contam_rows[[3]])
})

test_that("Cellwise replace modifies different rows per column", {
  nrows <- 30  # divisible by 3 columns
  ncols <- 3
  mat <- matrix(rnorm(nrows * ncols), nrow = nrows, ncol = ncols)

  cont <- add_contamination(
    mat,
    cont_pos = list(name="y", col=1:ncols),
    cont_settings = list(
      cont_prop = 0.1,
      cont_type = "replace",
      cont_value = 99,
      cont_mode = "cellwise"
    )
  )

  sets <- lapply(1:ncols, function(c) which(mat[,c] != cont[,c]))
  expect_false(identical(sets[[1]], sets[[2]]))
})

test_that("Cellwise outlier contamination produces columnwise independent spikes", {
  nrows <- 66  # divisible by 3 columns
  ncols <- 3
  mat <- matrix(rnorm(nrows * ncols), nrow = nrows, ncol = ncols)

  cont <- add_contamination(
    mat,
    cont_pos = list(name="y", col=1:ncols),
    cont_settings = list(
      cont_prop = 0.05,
      cont_type = "outlier",
      outlier_factor = 8,
      cont_mode = "cellwise"
    )
  )

  for (c in 1:ncols) {
    expect_true(any(abs(cont[,c] - mean(mat[,c])) > 5 * sd(mat[,c])))
  }
})

