#' Create an odds ratio estimate from a 2-by-2 table of frequencies or proportions
#'
#' @param x A two-dimensional matrix or table containing frequencies or proportions.
#'
#' @return A numeric vector.
#' @examples
#' # Load in Titanic data from datasets package
#' data(Titanic, package = "datasets")
#'
#' # Get marginal table of survival by sex
#' marginal_table <- margin.table(Titanic, c(2,4))
#' marginal_table
#'
#' # Compute odds ratio of marginal table
#' odds.ratio(marginal_table)
#'
#' # Get partial tables of survival by sex, stratified by class
#' partial_tables <- margin.table(Titanic, c(2,4,1))
#' partial_tables
#'
#' # Compute odds ratio of each partial table
#' apply(partial_tables, 3, odds.ratio)
#'
#' @author Paul W. Egeler, M.S.
#' @export
odds.ratio <- function (x) {

  stopifnot(identical(dim(x), c(2L,2L)))

  as.vector(x[1,1] * x[2,2] / x[1,2] / x[2,1])

}
