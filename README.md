samplesizeCMH: Sample Size Calculation for the Cochran-Mantel-Haenszel Test
===============

by Paul W. Egeler M.S.

## Description
This package provides functions relating to power and sample size calculation
for the CMH test. There are also several helper functions for interconverting
probability, odds, relative risk, and odds ratio values.

Please see the package vingettes for more information on how this package is used.

### The Cochran Mantel Haenszel Test

The **Cochran-Mantel-Haenszel test** (CMH) is an inferential test for
the association between two binary variables, while controlling for a third
confounding nominal variable. Two variables of interest, *X* and
*Y*, are compared at each level of the confounder variable *Z* and
the results are combined, creating a common odds ratio. Essentially, the CMH
test examines the *weighted* association of *X* and *Y*. The
CMH test is a common technique in the field of biostatistics, where it is
often used for case-control studies.

### Sample Size Calculation

Given a target power which the researcher would like to achieve, a
calculation can be performed in order to estimate the appropriate number of
subjects for a study. The `power.cmh.test` function calculates
the required number of subjects per group to achieve a specified power for a
Cochran-Mantel-Haenszel test.

### Power Calculation

Researchers interested in estimating the probability of detecting a true
positive result from an inferential test must perform a power calculation
using a known sample size, effect size, significance level, *et cetera*.
The `power.cmh.test` function can compute the power of a CMH test,
given parameters from the experiment.

## Installation

Downloading and installing the package from GitHub is facilitated by
[`devtools`](https://CRAN.R-project.org/package=devtools).
To do so, type the following into your R console:

    install.packages("devtools")
    devtools::install_github("pegeler/samplesizeCMH")

