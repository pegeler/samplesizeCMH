#' samplesizeCMH: Power and Sample Size Calculation for the
#' Cochran-Mantel-Haenszel Test
#'
#' This package provides functions relating to power and sample size calculation
#' for the CMH test. There are also several helper functions for interconverting
#' probability, odds, relative risk, and odds ratio values.
#'
#' The \strong{Cochran-Mantel-Haenszel test} (CMH) is an inferential test for
#' the association between two binary variables, while controlling for a third
#' confounding nominal variable. Two variables of interest, \emph{X} and
#' \emph{Y}, are compared at each level of the confounder variable \emph{Z} and
#' the results are combined, creating a common odds ratio. Essentially, the CMH
#' test examines the \emph{weighted} association of \emph{X} and \emph{Y}. The
#' CMH test is a common technique in the field of biostatistics, where it is
#' often used for case-control studies.
#'
#' @section Sample Size Calculation:
#'
#' Given a target power which the researcher would like to achieve, a
#' calculation can be performed in order to estimate the appropriate number of
#' subjects for a study. The \code{\link{power.cmh.test}} function calculates
#' the required number of subjects per group to achieve a specified power for a
#' Cochran-Mantel-Haenszel test.
#'
#' @section Power Calculation:
#'
#' Researchers interested in estimating the probability of detecting a true
#' positive result from an inferential test must perform a power calculation
#' using a known sample size, effect size, significance level, \emph{et cetera}.
#' The \code{\link{power.cmh.test}} function can compute the power of a CMH test,
#' given parameters from the experiment.
#'
#' @docType package
#' @name samplesizeCMH
NULL

