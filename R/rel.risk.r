#' Calculate the relative risk from a 2-by-2 table
#'
#' Computes the relative risk of a specified column of a two-by-two table.
#'
#' @param x A table or matrix containing frequencies.
#' @param col.num The column number upon which relative risk should be calculated.
#'
#' @return A numeric vector.
#' @author Paul W. Egeler, M.S.
#' @export
rel.risk <- function(x, col.num = 1) {
  stopifnot(identical(dim(x), c(2L, 2L)))
  y <- prop.table(x, 1)[, col.num]
  as.vector(y[1] / y[2])
}
