samplesizeCMH: Sample Size Calcualtion for the Cochran-Mantel-Haenszel Test
===============

by Paul. W. Egeler M.S. and Laura Kapitula, PhD

## Description
This R package provides a function for calculating the required number of
subjects per group to achieve a specified power for a Cochran-Mantel-Haenszel
test for stratified 2 x 2 tables using various methods.

There are also several helper functions for interconverting probability,
odds, relative risk, and odds ratio values.

Please see the package vingettes for more information on how this package is used.

## Installation

You must have [`devtools`](https://cran.r-project.org/web/packages/devtools/index.html)
installed in order to download and install this package. To do so, type the 
following into your R console:

    install.packages("devtools")

After you have `devtools`, type the following in the R console to download and 
install this package:

    devtools::install_github("pegeler/samplesizeCMH")
