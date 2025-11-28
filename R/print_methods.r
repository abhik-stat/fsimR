# =========== Common documentation for all print methods ================
#' @title Print Simulated Dataset
#'
#' @description
#' Prints a dataset from a simulation object, showing metadata and actual data.
#' Supports multiple types of simulated datasets:
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
#'
#' @param x Numeric vector, matrix, named list, or simulation object
#'    (of class `LMMdata`, `GLMMdata`, `LMdata`, `GLMdata`, `IIDdata`).
#'    If a list, it must include the components required for the selected `type`.
#'
#' @param replication Integer; which replication (dataset) to plot
#'    (LMM/GLMM, LM/GLM only; default `1`).
#' @param iteration Integer; which iteration within the replication to plot
#'    (LMM/GLMM, LM/GLM only; default `1`).
#' @param ... Any additional arguments.
#'
#' @details
#' The function prints:
#' - Metadata including number of subjects/observations, dimensions of design matrices,
#'    regression or fixed-effects coefficients, `SNR`, and `sigma_e`.
#' - Optional distribution names if `distr_settings` is provided within object.
#' - The actual dataset (`y`, `X`, `Z`, `RE`) for the specified replication and iteration.
#'
#' @examples
#' \dontrun{
#' # Mixed model
#' sim_data <- simulate_LMMdata(c(2,3), sim_spec(p=5))
#' print(sim_data, replication = 1, iteration = 2)
#'
#' # Regression model
#' sim_data_lm <- simulate_LMdata(c(2,3), sim_spec(p=5))
#' print(sim_data_lm)
#'
#' # IID samples
#' iid_data <- matrix(rnorm(50), ncol=5)
#' print(iid_data)
#'
#' # Named list
#' mylist <- list(y = rnorm(10), X = matrix(rnorm(20), ncol=2))
#' print(mylist)
#' }
#'
#' @name printSimData
NULL




# =========== Print method for LMMdata/GLMMdata ================
#' @rdname printSimData
#' @method print LMMdata
#' @export
print.LMMdata <- function(x, replication = 1, iteration = 1, ...) {

  # Ensure required components exist
  object <- x
  required <- c("y", "X", "Z", "RE", "sigma_e", "SNR", "ID")
  missing_comp <- setdiff(required, names(object))
  if (length(missing_comp) > 0) {
    stop("Object is missing required components: ", paste(missing_comp, collapse = ", "))
  }

  # Ensure y, X, Z, RE, SNR, sigma_e are lists
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

  # Print actual replication
  cat(sprintf("\n---- Replication %d, Iteration %d ----\n", replication, iteration))
  cat("Response vector (y):\n")
  print(y_sel)

  cat("\nFixed-effects design matrix (X):\n")
  print(X_sel)

  cat("\nRandom-effects design matrix (Z):\n")
  print(Z_sel)

  cat("\nRandom effects (RE):\n")
  print(RE_sel)

  invisible(NULL)
}



#' @rdname printSimData
#' @method print GLMMdata
#' @export
print.GLMMdata <- function(x, replication = 1, iteration = 1, ...) {
  print.LMMdata(x, replication, iteration, ...)
}



# =========== Print method for LMdata / GLMdata ================
#' @rdname printSimData
#' @method print LMdata
#' @export
print.LMdata <- function(x, replication = 1, iteration = 1, ...) {

  # Ensure required components exist
  object <- x
  required <- c("y", "X", "sigma_e", "SNR")
  missing_comp <- setdiff(required, names(object))
  if (length(missing_comp) > 0) {
    stop("Object is missing required components: ", paste(missing_comp, collapse = ", "))
  }

  # Ensure y, X, Z, RE, SNR are lists
  object[c("y", "X","SNR", "sigma_e")] <- lapply(object[c("y", "X", "SNR", "sigma_e")],
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

  # observations
  n_total <- length(y_sel)

  # Convert all components to data frames and combine once
  all_df <- as.data.frame(cbind(
    y = y_sel,
    X_sel))

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

  # Print actual replication
  cat(sprintf("\n---- Replication %d, Iteration %d ----\n", replication, iteration))
  cat("Response vector (y):\n")
  print(y_sel)

  cat("\nDesign matrix (X):\n")
  print(X_sel)


  invisible(NULL)
}



#' @rdname printSimData
#' @method print GLMMdata
#' @export
print.GLMMdata <- function(x, replication = 1, iteration = 1,...) {
  print.LMMdata(x, replication, iteration, ...)
}

