#' Create odds ratios estimates from tables
#'
#' @param x A two-dimensional matrix or table containing frequencies
#'
#' @return A numeric vector.
#' @author Paul W. Egeler, M.S.
#' @export
odds.ratio <- function (x) {

  stopifnot(length(dim(x)) == 2L && all.equal(dim(x)[1:2],c(2,2)))

  x[1,1] * x[2,2] / x[1,2] / x[2,1]

}
