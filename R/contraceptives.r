#' Oral contraceptive use and breast cancer rates
#'
#' This data summrizes the raw counts of the case-control study investigating
#' the link between breast cancer rates and oral contraceptive use, stratified
#' by age group. In toto, 10,890 subjects. See source for details.
#'
#' @source
#' Hennekens, C. H., F. E. Speizer, R. J. Lipnick, B. Rosner, C. Bain,
#' C. Belanger, M. J. Stampfer, W. Willett, and R. Peto. 1984. "A Case-Control
#' Study of Oral Contraceptive Use and Breast Cancer." \emph{Journal of the
#' National Cancer Institute} \strong{72} (1): 39–42. Table 1.
#' @docType data
#' @keywords datasets
#' @name contraceptives
#' @alias contraceptives_marginal
#' @usage data(contraceptives)
#' @format A 3-dimensional table.
#'
#' \enumerate{
#'   \item \code{OC Usage}: Subject exposure to oral contraceptives.
#'   \item \code{Disease Status}: Breast cancer present (case) or absent (control).
#'   \item \code{Age Group}: Age group of the subject.
#'
#' }
#'
#' \code{contraceptives_marginal} is the marginal table when age group is ignored.
NULL
