# =========== Common documentation for all summary methods ================
#' @title Summarize Simulated Data
#'
#' @description
#' Computes and prints a comprehensive summary for simulated datasets,
#' supporting multiple types:
#'
#' - **Mixed Models** (of class `LMMdata` or `GLMMdata`)
#' - **Regression Models** (of class `LMdata` or `GLMdata`)
#' - **IID samples** (of class `IIDdata`)
#'
#' The function also works for named lists containing:
#' - response vectors as `y`
#' - fixed-effects design matrix as `X`
#' - optional random-effects design matrix as `Z`
#' - optional random-effects realizations as `RE`
#'
#' Numeric vectors or matrices are treated as IID samples (columns = variables).
#'
#' For details of summary statistics, see Details.
#'
#' @param object Numeric vector, matrix, named list, or simulation object
#'    (of class `LMMdata`, `GLMMdata`, `LMdata`, `GLMdata`, `IIDdata`).
#'    If a list, it must include the components required for the selected `type`.
#'
#' @param replication Integer; which replication (dataset) to plot
#'    (LMM/GLMM, LM/GLM only; default `1`).
#' @param iteration Integer; which iteration within the replication to plot
#'    (LMM/GLMM, LM/GLM only; default `1`).
#'
#'
#' @details
#' The function prints a comprehensive summary for the simulated dataset:
#' - Simulation metadata such as number of subjects,
#'   observation counts, dimensions of `X`, `Z`, `RE` (for LMM/GLMM),
#'   fixed-effects coefficients, signal-to-noise ratio (`SNR`),
#'   and residual standard deviation (`sigma_e`).
#' - Optional distribution names if `distr_settings` is provided within object.
#' - Standard summary statistics (Mean, SD, Min, Q1, Median, Q3, Max)
#'   and total number of observations (`n`) for all variables in the object
#'   (`y`, `X`, `Z`, `RE` for mixed models; `y`, `X` for regression models;
#'   all columns for IID samples) for the specified `replication` and `iteration`.
#' - Correlation matrices for `X`, `Z`, and `RE` (if present),
#'   excluding intercept or any constant columns.
#'
#' @return
#' Invisibly returns a numeric matrix of computed summary statistics:
#' - **Rows:** variables present in the object (e.g., `y`, `X`, `Z`, `RE`).
#' - **Columns:** `Mean`, `SD`, `Min`, `Q1`, `Median`, `Q3`, `Max`, and total count `n`.
#'
#' @examples
#' \dontrun{
#' # Mixed Model
#' sim_settings <- sim_spec(p = 5)
#' lmm_data <- simulate_LMMdata(c(2,3), sim_settings)
#' summary(lmm_data, replication = 1, iteration = 1)
#'
#' # Regression Model
#' lm_data <- simulate_LMdata(c(2,3), sim_settings)
#' summary(lm_data, replication = 1, iteration = 1)
#'
#' # IID Samples (numeric matrix)
#' iid_data <- matrix(rnorm(100), ncol = 5)
#' summary(iid_data)
#'
#' # Named list with components like LMM
#' mylist <- list(
#'   y = lmm_data$y,
#'   X = lmm_data$X,
#'   Z = lmm_data$Z,
#'   RE = lmm_data$RE
#' )
#' summary(mylist, replication = 1, iteration = 1)
#' }
#'
#' @rdname summarySimData
NULL



# =========== Summary method for LMMdata/GLMMdata ================
#' @rdname summarySimData
#' @method summary LMMdata
#' @export
summary.LMMdata <- function(object, replication = 1, iteration = 1) {

  # Ensure required components exist
  required <- c("y", "X", "Z", "RE", "sigma_e", "SNR", "ID")
  missing_comp <- setdiff(required, names(object))
  if (length(missing_comp) > 0) {
    stop("Object is missing required components: ", paste(missing_comp, collapse = ", "))
  }

  # Ensure y, X, Z, RE, SNR are lists
  object[c("y", "X", "Z", "RE", "SNR", "sigma_e")] <- lapply(
                object[c("y", "X", "Z", "RE", "SNR", "sigma_e")],
                function(x) if (!is.list(x)) list(x) else x)

  # Number of simulations and iterations
  n_sim <- length(object$y)
  n_iter <- if (is.matrix(object$y[[1]])) ncol(object$y[[1]]) else length(object$y[[1]])

  # Bounds check
  if (replication > n_sim || replication < 1) stop("Replication index out of range")
  if (iteration > n_iter || iteration < 1) stop("Iteration index out of range")

  # Extract replication/iteration
  y_sel <- if (n_iter == 1) object$y[[replication]] else object$y[[replication]][, iteration, drop = TRUE]
  X_sel <- object$X[[replication]]
  Z_sel <- object$Z[[replication]]
  RE_sel <- if (n_iter == 1) object$RE[[replication]] else object$RE[[replication]][[iteration]]
  SNR_sel <- object$SNR[[replication]]
  sigma_e_sel <- object$sigma_e[[replication]]

  # beta coefficients
  beta_coeff <- if (!is.null(object$sim_settings$beta_coeff)) object$sim_settings$beta_coeff else rep(NA, ncol(X_sel))

  # Subjects and observations
  n_subj <- nlevels(object$ID)
  n_obs_vec <- as.integer(table(object$ID))
  n_total <- sum(n_obs_vec)

  # Ensure RE_sel is matrix and rename columns
  if (is.vector(RE_sel)) RE_sel <- matrix(RE_sel, ncol = 1)
  colnames(RE_sel) <- paste0("RE", seq_len(ncol(RE_sel)))

  # Convert all components to data frames and combine once
  all_df <- as.data.frame(cbind(
    y = y_sel,
    X_sel,
    Z_sel
  ))

  combined_summary <- rbind(
    Mean   = sapply(all_df, mean, na.rm = TRUE),
    SD     = sapply(all_df, sd, na.rm = TRUE),
    Min    = sapply(all_df, min, na.rm = TRUE),
    Q1     = sapply(all_df, function(x) quantile(x, 0.25, na.rm = TRUE)),
    Median = sapply(all_df, median, na.rm = TRUE),
    Q3     = sapply(all_df, function(x) quantile(x, 0.75, na.rm = TRUE)),
    Max    = sapply(all_df, max, na.rm = TRUE)
  )

  RE_sel <- as.data.frame(RE_sel)
  colnames(RE_sel) <- paste0("RE", seq_len(ncol(RE_sel)))

  RE_summary <- rbind(
    Mean   = sapply(RE_sel, mean, na.rm = TRUE),
    SD     = sapply(RE_sel, sd, na.rm = TRUE),
    Min    = sapply(RE_sel, min, na.rm = TRUE),
    Q1     = sapply(RE_sel, function(x) quantile(x, 0.25, na.rm = TRUE)),
    Median = sapply(RE_sel, median, na.rm = TRUE),
    Q3     = sapply(RE_sel, function(x) quantile(x, 0.75, na.rm = TRUE)),
    Max    = sapply(RE_sel, max, na.rm = TRUE)
  )


  combined_summary <- cbind(as.data.frame(combined_summary),
                            as.data.frame(RE_summary))

  # Transpose so variables are rows, stats are columns
  combined_summary <- t(combined_summary)

  # Add total n as extra column
  combined_summary <- cbind(combined_summary, n = n_total)

  # Print summary
  cat("==== LMM data Summary ====\n")
  cat("Total simulated datasets (replication):", n_sim, "\n")
  cat("Iterations per dataset:", n_iter, "\n")
  cat("Total observations per dataset:", n_total, "\n")
  cat("Number of subjects/clusters:", n_subj, "\n")
  if (length(unique(n_obs_vec)) > 1) {
    cat("Observations per subject (vector):", paste(n_obs_vec, collapse = ", "), "\n")
  } else {
    cat("Observations per subject (same for all):", n_obs_vec[1], "\n")
  }
  cat("Dimensions of X:", ncol(X_sel), "\n")
  cat("Dimensions of Z:", ncol(Z_sel), "\n")
  cat("Dimensions of RE:", paste(dim(RE_sel), collapse = " x "), "\n")
  cat("beta_coeff:", paste(round(beta_coeff, 3), collapse = ", "), "\n")
  cat("SNR (all iterations):", paste(round(SNR_sel, 3), collapse = ", "), "\n")
  cat("sigma_e (all iterations):", paste(round(sigma_e_sel, 3), collapse = ", "), "\n")

  if (!is.null(object$distr_settings)) {
    cat("\n---- Distribution Settings ----\n")
    for (nm in names(object$distr_settings)) {
      distr <- object$distr_settings[[nm]]
      if (!is.null(distr$distr_name)) cat(sprintf("%s: %s\n", nm, distr$distr_name))
    }
    cat("---------------------------------------\n")
  }

  cat(sprintf("\n---- Replication %d, Iteration %d ----\n", replication, iteration))
  cat("Combined summary (variables x statistics):\n")
  print(combined_summary)

  cat(sprintf("\nCorrelation matrix of X (excluding intercept):\n"))
  X_sel_nonzero <- X_sel[, sapply(as.data.frame(X_sel), sd) != 0, drop = FALSE]
  print(cor(X_sel_nonzero))

  cat(sprintf("\nCorrelation matrix of Z (excluding intercept):\n"))
  Z_sel_nonzero <- Z_sel[, sapply(as.data.frame(Z_sel), sd) != 0, drop = FALSE]
  print(cor(Z_sel_nonzero))

  cat(sprintf("\nCorrelation matrix of RE (excluding intercept):\n"))
  RE_sel_nonzero <- RE_sel[, sapply(as.data.frame(RE_sel), sd) != 0, drop = FALSE]
  print(cor(RE_sel_nonzero))

  invisible(combined_summary)
}


#' @rdname summarySimData
#' @method summary GLMMdata
#' @export
summary.GLMMdata <- function(object, replication = 1, iteration = 1) {
  summary.LMMdata(object, replication, iteration)
}



# =========== Summary method for LMdata/GLMdata ================
#' @rdname summarySimData
#' @method summary LMdata
#' @export
summary.LMdata <- function(object, replication = 1, iteration = 1) {

  # Ensure required components exist
  required <- c("y", "X", "sigma_e", "SNR")
  missing_comp <- setdiff(required, names(object))
  if (length(missing_comp) > 0) {
    stop("Object is missing required components: ", paste(missing_comp, collapse = ", "))
  }

  # Ensure y, X, Z, RE, SNR are lists
  object[c("y", "X", "SNR", "sigma_e")] <- lapply(object[c("y", "X", "SNR", "sigma_e")],
                                       function(x) if (!is.list(x)) list(x) else x)

  # Number of simulations and iterations
  n_sim <- length(object$y)
  n_iter <- if (is.matrix(object$y[[1]])) ncol(object$y[[1]]) else length(object$y[[1]])

  # Bounds check
  if (replication > n_sim || replication < 1) stop("Replication index out of range")
  if (iteration > n_iter || iteration < 1) stop("Iteration index out of range")

  # Extract replication/iteration
  y_sel <- if (n_iter == 1) object$y[[replication]] else object$y[[replication]][, iteration, drop = TRUE]
  X_sel <- object$X[[replication]]
  SNR_sel <- object$SNR[[replication]]
  sigma_e_sel <- object$sigma_e[[replication]]

  # beta coefficients
  beta_coeff <- if (!is.null(object$sim_settings$beta_coeff)) object$sim_settings$beta_coeff else rep(NA, ncol(X_sel))

  # Subjects and observations
  n_total <- length(y_sel)


  # Convert all components to data frames and combine once
  all_df <- as.data.frame(cbind(y = y_sel, X_sel))

  combined_summary <- rbind(
    Mean   = sapply(all_df, mean, na.rm = TRUE),
    SD     = sapply(all_df, sd, na.rm = TRUE),
    Min    = sapply(all_df, min, na.rm = TRUE),
    Q1     = sapply(all_df, function(x) quantile(x, 0.25, na.rm = TRUE)),
    Median = sapply(all_df, median, na.rm = TRUE),
    Q3     = sapply(all_df, function(x) quantile(x, 0.75, na.rm = TRUE)),
    Max    = sapply(all_df, max, na.rm = TRUE)
  )
  colnames(combined_summary)[1] <- "y"

  # Transpose so variables are rows, stats are columns
  combined_summary <- t(combined_summary)

  # Add total n as extra column
  combined_summary <- cbind(combined_summary, n = n_total)

  # Print summary
  cat("==== LMM data Summary ====\n")
  cat("Total simulated datasets (replication):", n_sim, "\n")
  cat("Iterations per dataset:", n_iter, "\n")
  cat("Total observations per dataset:", n_total, "\n")
  cat("Dimensions of X:", ncol(X_sel), "\n")
  cat("beta_coeff:", paste(round(beta_coeff, 3), collapse = ", "), "\n")
  cat("SNR (all iterations):", paste(round(SNR_sel, 3), collapse = ", "), "\n")
  cat("sigma_e (all iterations):", paste(round(sigma_e_sel, 3), collapse = ", "), "\n")

  if (!is.null(object$distr_settings)) {
    cat("\n---- Distribution Settings ----\n")
    for (nm in names(object$distr_settings)) {
      distr <- object$distr_settings[[nm]]
      if (!is.null(distr$distr_name)) cat(sprintf("%s: %s\n", nm, distr$distr_name))
    }
    cat("---------------------------------------\n")
  }

  cat(sprintf("\n---- Replication %d, Iteration %d ----\n", replication, iteration))
  cat("Combined summary (variables x statistics):\n")
  print(combined_summary)

  cat(sprintf("\nCorrelation matrix of X (excluding intercept):\n"))
  X_sel_nonzero <- X_sel[, sapply(as.data.frame(X_sel), sd) != 0, drop = FALSE]
  print(cor(X_sel_nonzero))

  invisible(combined_summary)
}


#' @rdname summarySimData
#' @method summary GLMdata
#' @export
summary.GLMdata <- function(object, replication = 1, iteration = 1) {
  summary.LMdata(object, replication, iteration)
}


# =========== Summary method for IIDdata ================
#' @rdname summarySimData
#' @method summary IIDdata
#' @export
summary.IIDdata <- function(data) {


  # format data suitably
  data <- as.matrix(data)
  n <- nrow(data)
  d <- ncol(data)
  if(is.null(colnames(data))) colnames(data) <- paste0("X", 1:d)

  # Convert all components to data frames and combine once
  all_df <- as.data.frame(data)

  combined_summary <- rbind(
    Mean   = sapply(all_df, mean, na.rm = TRUE),
    SD     = sapply(all_df, sd, na.rm = TRUE),
    Min    = sapply(all_df, min, na.rm = TRUE),
    Q1     = sapply(all_df, function(x) quantile(x, 0.25, na.rm = TRUE)),
    Median = sapply(all_df, median, na.rm = TRUE),
    Q3     = sapply(all_df, function(x) quantile(x, 0.75, na.rm = TRUE)),
    Max    = sapply(all_df, max, na.rm = TRUE)
  )

  # Transpose so variables are rows, stats are columns
  combined_summary <- t(combined_summary)

  # Add total n as extra column
  combined_summary <- cbind(combined_summary, n = n)

  cat("Combined summary (variables x statistics):\n")
  print(combined_summary)

  cat(sprintf("\nCorrelation matrix (excluding constant columns):\n"))
  X_sel_nonzero <- all_df[, sapply(all_df, sd) != 0, drop = FALSE]
  print(cor(X_sel_nonzero))

  invisible(combined_summary)
}



# =========== Summary method for Other classes ================

# Numeric wrapper -> IIDdata
#' @rdname summarySimData
#' @method summary numeric
#' @export
summary.numeric  <- function(object) {
  summary.IIDdata(object)
}


# Matrix wrapper -> IIDdata
#' @rdname summarySimData
#' @method summary matrix
#' @export
summary.matrix  <- function(object) {
  summary.IIDdata(object)
}


# List wrapper
#' @rdname summarySimData
#' @method summary list
#' @export
summary.list <- function(object, replication = 1, iteration = 1) {
  if (all(c("y","X","Z","RE") %in% names(object))) {
    summary.LMMdata(object, replication = 1, iteration = 1)
  } else if (all(c("y","X") %in% names(object))) {
    summary.LMdata(object, replication = 1, iteration = 1)
  } else {
    stop("Cannot determine summary method for this list. Ensure it contains valid components.")
  }
}
