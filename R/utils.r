#' @title Fallback operator (Internal)
#'
#' @description
#' Returns `x` if not `NULL`, otherwise returns `y`.
#'
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) {
  if (!is.null(x)) x else y
}



#' ===================================================================================
#'
#' @title Resolve Alias Groups for Distribution Specification (Internal)
#'
#' @description
#' Internal helper used by `distr_spec()` to consolidate multiple alias
#' arguments for the same conceptual parameter (e.g., `mean` / `mu` /
#' `location` / `xi` for location parameter) into a single canonical representation.
#' The function selects the **first non-`NULL` alias** in the order
#' provided and interprets it according to `type`. It also updates the
#' implied dimension of the parameter when necessary.
#'
#'
#' @param values A named list of alias values. Entries may be `NULL`;
#'            the first non-`NULL` element determines the returned value.
#' @param dim Integer scalar specifying the current dimension,
#'            used for expanding scalars or validating shapes.
#' @param default Scalar default value to use when **no** alias is supplied.
#' @param type Character string, either `"vector"` or `"matrix"`,
#'        indicating the expected structure of the resolved parameter.
#'
#' @return A list with components:
#' \describe{
#'   \item{value}{The resolved vector or matrix.}
#'   \item{dim}{The updated dimension based on the provided alias value (may differ from input `dim`).}
#'   \item{has_alias}{Logical; `TRUE` if at least one alias was supplied, `FALSE` otherwise.}
#' }
#'
#'
#' @details
#' Depending on the `type` argument,
#' the user-supplied alias value is interpreted as:
#'
#'**Vector mode (`type = "vector"`)**
#' * Accepts scalars, atomic vectors, or matrices with one row or one column
#' * (which are then flattened into vectors).
#' * Multi-row and multi-column matrices are rejected.
#' * Scalar value is expanded to a vector of length `dim`.
#' * If a non-scalar vector is provided, `dim` is updated to its length.
#'
#'**Matrix mode (`type = "matrix"`)**
#' * Accepts scalars (expanded to `scalar \eqn{\times} diag(dim)`) or square matrices.
#' * If a square matrix is provided, `dim` is updated to its order.
#'
#' **Default behavior**
#' * If no alias is provided, the default is expanded to:
#'   * `rep(default, dim)`   for vectors
#'   * `default \eqn{\times} diag(dim)` for matrices
#'
#'
#' @keywords internal
#' @noRd
.resolve_alias_group <- function(values, dim, default = 0, type = c("vector", "matrix")) {

  type <- match.arg(type)

  # Extract non-null values in the order provided
  nonnull <- Filter(Negate(is.null), values)

  # -------------------------------------------------------------------
  # CASE 1: No alias given → return default value
  # -------------------------------------------------------------------
  if (length(nonnull) == 0) {
    if (type == "vector") {
      return(list(value = rep(default, dim),
                  dim   = dim,
                  has_alias = FALSE))
    } else {
      return(list(value = default * diag(dim),
                  dim   = dim,
                  has_alias = FALSE))
    }
  }

  # first alias wins
  v <- nonnull[[1]]

  # ------------------------------------------------------------
  # VECTOR CASE
  # ------------------------------------------------------------
  if (type == "vector") {

    # Matrix-like vector allowed (1 row or 1 col)
    if (is.matrix(v)) {
      if (!(nrow(v) == 1 || ncol(v) == 1)) {
        stop("Vector alias must be scalar, vector, or 1-row/1-col matrix.")
      }
      v <- as.vector(v)
    }

    # SCALAR → EXPAND
    if (length(v) == 1) {
      return(list(value = rep(v, dim),
                  dim   = dim,
                  has_alias = TRUE))
    }

    # NON-SCALAR vector updates dimension
    return(list(value = v,
                dim   = length(v),
                has_alias = TRUE))
  }

  # ------------------------------------------------------------
  # MATRIX CASE
  # ------------------------------------------------------------
  if (type == "matrix") {

    # SCALAR CASE → scalar*I(dim)
    if (!is.matrix(v)) {
      if (length(v) != 1) {
        stop("Matrix alias must be scalar or a square matrix.")
      }
      return(list(value = v * diag(dim),
                  dim   = dim,
                  has_alias = TRUE))
    }

    # MATRIX CASE
    if (nrow(v) != ncol(v))
      stop("Matrix alias must be square.")

    return(list(value = v,
                dim   = nrow(v),
                has_alias = TRUE))
  }

}



#' ===================================================================================
#' @title Strict Dimension Compatibility Checker (Internal)
#'
#' @description
#' Internal helper used by `distr_spec()` to enforce that any
#' user-supplied parameter has a dimension consistent with the final
#' resolved dimension of the distribution specification.
#'
#' @param name Character string naming the parameter (used in error messages).
#' @param user_dim Integer dimension of the user-supplied value, or `NULL`
#' if no explicit value was provided.
#' @param target_dim Integer scalar giving the final resolved dimension.
#'
#' @return Invisibly returns `TRUE` if dimensions are compatible; otherwise
#' throws an informative error.
#'

#' @details
#' If the parameter was not supplied (`user_dim = NULL`), the check is skipped.
#'
#' A mismatch between `user_dim` and `target_dim` produces an error of the form:
#' `"param_name has dimension k but resolved dimension is d."`
#'
#' @keywords internal
#' @noRd
.check_dim_strict <- function(name, user_dim, target_dim) {
  if (!is.null(user_dim) && user_dim != target_dim) {
    stop(sprintf("%s has dimension %d but resolved dimension is %d.",
                 name, user_dim, target_dim), call. = FALSE)
  }
  invisible(TRUE)
}

