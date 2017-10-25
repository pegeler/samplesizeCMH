#' samplesizeCMH: Sample Size Calcualtion for the Cochran-Mantel-Haenszel Test
#'
#' This package provides a function for calculating the required number of
#' subjects per group to achieve a specified power for a Cochran-Mantel-Haenszel
#' test using various methods.
#' There are also several helper functions for interconverting probability,
#' odds, relative risk, and odds ratio values.
#'
#' @section Sample Size Calculation:
#' This is where I talk about the main function that calculates sample size.
#'
#' @section Power Calculation:
#' I am thinking about adding a power calculator too.
#'
#' @section Helper Functions:
#' This is where I talk about the odds and proportions stuff
#'
#' @docType package
#' @name samplesizeCMH
NULL

# Warn the user that this package is still under construction
.onAttach <- function(libname, pkgname) {

  packageStartupMessage(
  "This package is under construction!\n",
  "------------------------------------\n",
  "It has not been tested for accuracy.\n",
  "Some features may be absent."
  )

}
