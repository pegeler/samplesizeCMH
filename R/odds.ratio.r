#' Create odds ratios estimates from a table
#'
#' @param x A two-dimensional matrix or table containing frequencies
#'
#' @return A numeric vector.
#' @author Paul W. Egeler, M.S.
#' @export
odds.ratio <- function (x) {

  stopifnot(identical(dim(x), c(2L,2L)))

  x[1,1] * x[2,2] / x[1,2] / x[2,1]

}
