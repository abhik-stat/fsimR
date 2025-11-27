#' @title Add Contamination to Simulated Data
#'
#' @description
#' Introduces controlled contamination into simulated datasets.
#'
#' Contamination can be applied to:
#' - **Mixed Model datasets** (of class `LMMdata` or `GLMMdata`)
#' - **Regression datasets** (of class `LMdata` or `GLMdata`)
#' - **IID samples** (of class `IIDdata`)
#' - **Named lists** containing components such as:
#'   - `y`: response vector or matrix
#'   - `X`: fixed-effects design matrix
#'   - `Z`: random-effects design matrix (optional)
#'
#' Numeric vectors or matrices are treated as IID samples (columns = variables).
#'
#'
#' @param simData Numeric vector, matrix, named list, or simulation object
#'    (of class `LMMdata`, `GLMMdata`, `LMdata`, `GLMdata`, `IIDdata`).
#'
#'   - If a **list**, it must include the component specified in `cont_pos$name`.
#'   - If a **vector/matrix**, it is treated as the target component to contaminate.
#'
#' @param cont_pos Named list specifying the target component to contaminate:
#' \describe{
#'   \item{name}{Character string, one of `"y"`, `"X"`, or `"Z"` (default: `y`).
#'     Determines which component of the simulation object is modified.}
#'   \item{col}{Vector of column indices to contaminate.
#'     If `col = NULL` (default) or `name = "y"`, all columns are contaminated.
#'     For `y`, multiple columns usually represent replications in the simulation object,
#'     and all such columns are contaminated consistently.}
#' }
#'
#'
#' @param cont_settings Named list specifying contamination details:
#' \describe{
#'   \item{cont_mode}{Mode of contamination: `"casewise"` (default) or `"cellwise"`.
#'     - `"casewise"`: Same rows contaminated across all specified columns.
#'     - `"cellwise"`: Different rows may be contaminated independently per column.}
#'   \item{cont_prop}{Proportion of contamination (0–1).
#'      In `casewise` mode, this is the proportion of rows.
#'     In `cellwise` mode, this is the proportion of cells across all selected columns.
#'     (default: `0.1`)}
#'
#'   \item{cont_value}{Optional numeric scalar or vector specifying contamination values.
#'     If provided, this overrides `cont_distr` and no random generation occurs.
#'     *Length must be either 1 or equal to the number of contaminated columns*.}
#'
#'   \item{cont_distr}{Named list describing the distribution of contamination
#'     values (only when `cont_value` is not provided).
#'     Accepts:
#'     \describe{
#'       \item{distr_name}{Distribution name (e.g., `"norm"`, `"t"`, `"unif"`).}
#'       \item{distr_params}{List of distribution parameters (see [simulate_IIDdata]).}
#'       \item{generator}{Optional custom generator function.}
#'     }
#'     The default is the standard normal distribution.}
#'
#'   \item{cont_type}{Character string specifying the contamination mechanism:
#'     \describe{
#'       \item{"additive"}{Adds contamination values to selected rows/cells (default).}
#'
#'       \item{"multiplicative"}{Multiplies values in selected rows/cells
#'          by contamination values.}
#'
#'       \item{"replace"}{Replaces values in selected rows/cells with contamination values.}
#'
#'       \item{"outlier"}{Sets selected rows/cells to extreme values of the form:
#'     \deqn{\text{mean} \pm (\text{sd} \times \text{outlier_factor}).}
#'        This generates classical point outliers.}
#'
#'       \item{"leverage"}{Generates high-leverage points in `X` or `Z`.
#'         Corresponding rows in `y` (if present) are shifted in the opposite direction.
#'         See *Details* for full description.}
#'     }}
#'
#'   \item{leverage_factor}{Numeric multiplier controlling the magnitude
#'   of leverage contamination  (default: `5`).}
#'
#'   \item{outlier_factor}{Numeric multiplier controlling extreme values for `"outlier"`
#'   or `y` shifts in `"leverage"` (default: `5`).}
#'
#'   \item{outlier_dir}{Probability that an outlier is generated in the right
#'     (positive) tail versus the left (negative) tail.
#'     Default: `1` (all outliers on the right).}
#' }
#'
#'
#' @details
#' The function provides a general framework for injecting various forms of contamination
#' into simulated data, useful for robustness studies and sensitivity analyses.
#'
#' **Contamination mode:**
#' - *Casewise*: A total of `ceil(n × cont_prop)` rows are selected and
#' contaminated across all columns.
#' - *Cellwise*: A proportion of individual cells is contaminated independently
#' in each specified column. Rows are re-sampled independently per column.
#'
#'
#' *Leverage contamination:*
#' For `cont_type = "leverage"` applied to `X` or `Z`,
#' a subset of rows/cells is randomly selected for contamination.
#' Their values are set to the column mean plus a random shift,
#' drawn from a normal distribution with the column's standard deviation scaled by `leverage_factor`.
#' If a `y` component exists, the corresponding rows are shifted in the opposite direction,
#' scaled by `outlier_factor`, creating influential points within the regression framework.
#'
#' *Distribution-based contamination:*
#' When `cont_value = NULL`, random contamination values are generated using
#' `simulate_IIDdata()`, allowing flexible custom distributions or user-provided generators.
#'
#'
#'
#' @return
#' Modified simulated datasets of the same class/structure as the input,
#' with the requested contamination applied.
#' - If `simData` is a list, the contaminated component is replaced.
#' - If a vector or matrix is supplied, the returned object is the contaminated matrix.
#'
#'
#' @examples
#' # Additive contamination to an IID sample
#' set.seed(1)
#' x <- rnorm(100)
#' contaminated <- add_contamination(
#'   x,
#'   cont_pos = list(name="y", col=1),
#'   cont_settings = list(
#'     cont_prop = 0.1,
#'     cont_type = "additive",
#'     cont_distr = list(distr_name="norm", distr_params=list(mean=5, sd=1))
#'   )
#' )
#'
#'
#' # Outlier contamination in regression response
#' set.seed(2)
#' sim <- list(y = rnorm(200), X = cbind(1, rnorm(200)))
#'
#' sim_outliers <- add_contamination(
#'   sim,
#'   cont_pos = list(name="y"),
#'   cont_settings = list(
#'     cont_type="outlier",
#'     outlier_factor=8,
#'     outlier_dir=0.8
#'   )
#' )
#'
#'
#' # High-leverage contamination in design matrix X
#' set.seed(3)
#' sim <- list(
#'   y = rnorm(150),
#'   X = cbind(1, rnorm(150), rnorm(150))
#' )
#'
#' sim_lev <- add_contamination(
#'   sim,
#'   cont_pos = list(name="X", col=2:3),
#'   cont_settings = list(
#'     cont_type="leverage",
#'     leverage_factor=10,
#'     outlier_factor=6
#'   )
#' )
#'
#' # Cellwise contamination (different rows per column)
#' set.seed(4)
#' mat <- matrix(rnorm(50*3), ncol=3)
#' cont_cell <- add_contamination(
#'   mat,
#'   cont_pos = list(name="X", col=1:3),
#'   cont_settings = list(
#'     cont_prop=0.1,
#'     cont_type="replace",
#'     cont_value=99,
#'     cont_mode="cellwise"
#'   )
#' )
#'
#' @importFrom stats rnorm mean sd colMeans ncol nrow runif ceiling
#' @importFrom utils modifyList
#'
#' @export
add_contamination <- function(simData,
                              cont_pos = list(),
                              cont_settings = list()) {

  # Merge settings
  cont_pos <- modifyList(list(name = "y", col = NULL), cont_pos)

  default_settings <- list(
    cont_prop = 0.1,
    cont_value = NULL,
    cont_distr = list(distr_name="norm",
                      distr_params = distr_spec(),
                      generator = NULL),
    cont_type = "additive",
    cont_mode = "casewise",
    leverage_factor = 5,
    outlier_factor = 5,
    outlier_dir = 1
  )
  cont_settings <- modifyList(default_settings, cont_settings)

  # Extract contamination characteristics
  cont_prop <- cont_settings$cont_prop
  cont_value <- cont_settings$cont_value
  cont_type <- cont_settings$cont_type
  cont_mode       <- cont_settings$cont_mode
  leverage_factor <- cont_settings$leverage_factor
  outlier_factor <- cont_settings$outlier_factor
  cont_distr <- do.call(modifyList, list(default_settings$cont_distr, cont_settings$cont_distr))


  # Extract component to contaminate
  comp_name <- cont_pos$name

  if (is.list(simData)) {
    if (!comp_name %in% c("y","X","Z")) stop("Invalid component name")
    comp_data <- simData[[comp_name]]
  } else {
    comp_data <- simData
  }

  # Ensure list of matrices
  if (!is.list(comp_data)) {
    if (!is.matrix(comp_data)) {
      if (!is.numeric(comp_data)) stop("Target data must be numeric.")
      comp_data <- matrix(comp_data, ncol = 1)
    }
    comp_data1 <- list()
    comp_data1[[comp_name]] <- comp_data
    comp_data <- comp_data1
    # comp_data <- list(comp_data)
  }

  # Determine columns to contaminate
  cols <- cont_pos$col
  n_col0 <- ncol(comp_data[[1]])
  if (any(cols > n_col0)) stop("Some columns exceed data dimensions")
  col_idx <- if (comp_name == "y" || is.null(cols)) seq_len(n_col0) else cols
  n_col <- length(col_idx)

  # Number of contaminated rows/cells
  n <- nrow(comp_data[[1]])

  if (!(cont_mode %in% c("casewise","cellwise"))) {
    stop("cont_mode must be 'casewise' or 'cellwise'")
  }

  if (cont_mode == "casewise") {

    n_cont <- ceiling(n * cont_prop)                # rows to contaminate
    row_idx_casewise <- sample(n, n_cont, replace=FALSE)

  } else {  # CELLLWISE

    # proportion applies to CELLS in selected columns
    total_cells <- n * n_col
    target_cells <- ceiling(total_cells * cont_prop)
    n_cont <- min(ceiling(target_cells / n_col), n)
    row_idx_casewise <- NULL
  }

  # n_cont <- min(n, ceiling(n*cont_prop))
  if (n_cont == 0) return(simData)
  y_cont = list()

  # Apply contamination across list
  for (i in seq_along(comp_data)) {

    mat <- comp_data[[i]]

    # Outlier contamination ---
    if (cont_type == "outlier") {
      for (k in seq_along(col_idx)) {

        row_idx <- if (cont_mode == "casewise") row_idx_casewise else
                              sample(n, n_cont, replace = FALSE)

        col_values <- mat[, col_idx[k]]
        col_mean <- mean(col_values, na.rm = TRUE)
        col_sd <- sd(col_values, na.rm = TRUE)

        signs <- ifelse(runif(n_cont) < cont_settings$outlier_dir, 1, -1)

        mat[row_idx, col_idx[k]] <- col_mean + outlier_factor*col_sd*signs
      }
    }


    # Leverage-influential contamination ---
    if (cont_type == "leverage") {
      if (!comp_name %in% c("X","Z"))
        stop("Leverage-influential only applicable to X or Z")

      # col_idx <- if (is.null(cols)) seq_len(ncol(mat)) else cols  # columns to contaminate
      col_means <- colMeans(mat[, col_idx, drop = FALSE])
      col_sds <- apply(mat[, col_idx, drop = FALSE], 2, sd)

      # Contaminate each column separately
      for (k in seq_along(col_idx)) {

        row_idx <- if (cont_mode == "casewise") row_idx_casewise else
          sample(n, n_cont, replace = FALSE)

        mat_shift <- rnorm(n_cont) * leverage_factor * col_sds[k]
        mat[row_idx, col_idx[k]] <- col_means[k] + mat_shift
      }


      # Contaminate corresponding y in opposite direction
      if ("y" %in% names(simData)) {

        y_mat <- if(is.list(simData$y)) simData$y[[i]] else simData$y
        y_mean <- if (is.matrix(y_mat)) colMeans(y_mat) else mean(y_mat)
        y_sd <- if (is.matrix(y_mat)) apply(y_mat, 2, sd) else sd(y_mat)

        # Compute shift direction based on sum of contaminated columns
        shift_dir <- sign(rowSums(
                            mat[row_idx, col_idx, drop = FALSE] -
                            matrix(col_means, nrow = n_cont, ncol = n_col, byrow = TRUE)
                            ))

        if (is.matrix(y_mat)) {
          for (j in seq_len(ncol(y_mat))) {
            y_mat[row_idx, j] <- y_mean[j] - outlier_factor*y_sd[j]*shift_dir
          }
        } else {
          y_mat[row_idx] <- y_mean - outlier_factor*y_sd*shift_dir
        }
        y_cont[[i]] <- y_mat
      }
    }


    # Standard contamination: additive/multiplicative/replace ---
    if (cont_type %in% c("additive","multiplicative","replace")) {
      if (!is.null(cont_value)) {
        if (length(cont_value) == 1)
          cont_vals <- matrix(cont_value, nrow = n_cont, ncol = n_col)
        else if (length(cont_value) == n_col)
          cont_vals <- matrix(rep(cont_value, each = n_cont), nrow=n_cont)
        else stop("cont_value must be length 1 or equal to number of columns")
      } else {
        vals <- simulate_IIDdata(n_cont * n_col,
                                 distr_name = cont_distr$distr_name,
                                 distr_params = cont_distr$distr_params,
                                 generator=cont_distr$generator)
        cont_vals <- matrix(vals, nrow = n_cont, ncol = n_col)
      }

      for (k in seq_along(col_idx)) {

        row_idx <- if (cont_mode == "casewise") row_idx_casewise else
          sample(n, n_cont, replace = FALSE)

        if (cont_type == "additive")
          mat[row_idx, col_idx[k]] <- mat[row_idx, col_idx[k]] + cont_vals[, k]
        else if (cont_type == "multiplicative")
          mat[row_idx, col_idx[k]] <- mat[row_idx, col_idx[k]] * cont_vals[, k]
        else if (cont_type == "replace")
          mat[row_idx, col_idx[k]] <- cont_vals[, k]
      }
    }

    comp_data[[i]] <- mat
  }

  if (length(y_cont) > 0) simData$y <- if (length(y_cont) > 1) y_cont else y_mat

  if (length(comp_data) == 1) comp_data <- comp_data[[1]]

  if (!is.list(simData)) {simData <- comp_data} else {simData[[comp_name]] <- comp_data}

  return(simData)
}


