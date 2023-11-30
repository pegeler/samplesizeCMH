samplesizeCMH: Sample Size Calculation for the Cochran-Mantel-Haenszel Test
===============

[![CRAN\_version](http://www.r-pkg.org/badges/version/samplesizeCMH)](https://cran.r-project.org/package=samplesizeCMH)
[![Number\_of\_Downloads](https://cranlogs.r-pkg.org/badges/grand-total/samplesizeCMH)](https://cran.r-project.org/package=samplesizeCMH)

by Paul W. Egeler M.S.

## Description
This package provides functions relating to power and sample size calculation
for the CMH test. There are also several helper functions for interconverting
probability, odds, relative risk, and odds ratio values.

Please see the [package website](https://pegeler.github.io/samplesizeCMH/) for more information on how this package is used, including [documentation](https://pegeler.github.io/samplesizeCMH/reference/) and [vignettes](https://pegeler.github.io/samplesizeCMH/articles/).

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

Installation of the CRAN release can be done with `install.packages()`. From the R console:

```r
install.packages("samplesizeCMH")
```

Downloading and installing the latest version from GitHub is facilitated by
[`remotes`](https://CRAN.R-project.org/package=remotes).
To do so, type the following into your R console:

```r
if (!require("remotes")) install.packages("remotes")
remotes::install_github("pegeler/samplesizeCMH")
```
