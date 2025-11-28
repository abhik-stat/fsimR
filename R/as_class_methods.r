#' @title Coercion helpers for package-specific S3 classes
#'
#' @description
#' These functions provide lightweight coercion into the package’s custom classes,
#' such as `covMat`, `IIDdata`, `LMdata`, `LMMdata`, etc.
#' Each coercion function has the form: \code{as.<ClassName>(object)}
#' and simply prepends the corresponding class name to the object's
#' existing class vector.
#' This enables S3 dispatch while preserving the underlying data structure.
#'
#' @param object Any R object.
#'
#' @return The same object, with the target class prepended to its class vector
#' (if the object is compatible with the class).
#'
#' @name asClass
#'
#' @aliases
#' as.covMat
#' as.IIDdata
#' as.LMMdata
#' as.LMdata
#'
#' @export as.covMat
#' @export as.IIDdata
#' @export as.LMMdata
#' @export as.LMdata
NULL


##########################################################################
##  Individually defined coercion functions

#' @rdname asClass
as.covMat <- function(object) {
  if (is.matrix(object) || is.numeric(object))
    class(object) <- c("covMat", class(object))
  object
}

#' @rdname asClass
as.IIDdata <- function(object) {
  if (is.matrix(object) || is.data.frame(object) || is.numeric(object))
    class(object) <- c("IIDdata", class(object))
  object
}

#' @rdname asClass
as.LMMdata <- function(object) {
  if (is.list(object) || is.data.frame(object))
    class(object) <- c("LMMdata", class(object))
  object
}

#' @rdname asClass
as.LMdata <- function(object) {
  if (is.list(object) || is.data.frame(object))
    class(object) <- c("LMdata", class(object))
  object
}
