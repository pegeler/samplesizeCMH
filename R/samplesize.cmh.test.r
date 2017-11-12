#' @title Sample size calculation for the Cochran Mantel Haenszel test
#'
#' @description
#' Compute the required number of subjects for the Cochran Mantel Haenszel test
#' for association in \emph{J} stratified 2 x 2 tables. The calculations can be
#' performed using various methods found in literature.
#'
#' @details
#' Terminology and symbolic conventions are borrowed from Woolson \emph{et al.}
#' (1986). The \code{p1} group is dubbed the \emph{Case} group and \code{p2}
#' group is called the \emph{Control} group.
#'
#' The \code{J} number of groups will be inferred by the length of vector
#' \code{t}.
#'
#' @section Methods:
#'
#' \describe{
#'
#'   \item{\code{cc.binomial}}{
#'     Continuity corrected weighted binomial sample size estimator described
#'     in the Nam (1992). This uses the weighted binomial sample size
#'     calculation described in Woolson \emph{et al.} (1986) but is enhanced
#'     for use with the continuity corrected Cochran's test.
#'   }
#'
#'   \item{\code{binomial}}{
#'     This sample size calculation is designed for case-control studies
#'     where one margin is fixed. The method is "based on the
#'     Cochran-Mantel-Haenszel statistic expressed as a weighted difference in
#'     binomial proportions."
#'     Described in the Woolson \emph{et al.} (1986).
#'   }
#'
#'   \item{\code{hypergeometric}}{
#'     Described in Munoz and Rosner (1984). This method uses the Cornfield
#'     approximation to determine asymptotic power/sample size. Since this uses
#'     the hypergeometric distribution for the Mantel-Haenszel test, both
#'     margins are assumed to be fixed. The \code{r} parameter represents the
#'     second margin allocation for each \emph{j}th stratum and is required.
#'
#'     Woolson \emph{et al.} (1986) state that this is the appropriate sample
#'     size estimator when both margins are fixed.
#'   }
#'
#'   \item{\code{least.squares}}{
#'     Described by Gail (1973). This has been superceded by subsequent methods
#'     described by Munoz and Rosner (1984), Woolson \emph{et al.} (1986), and
#'     later Nam (1992).
#'
#'     This method has several assumptions. Namely, the number subjects
#'     allocated to each case and control group within each \emph{j}th stratum
#'     must be equal. A warning is produced if this assumption is violated.
#'   }
#'
#'   \item{\code{unstratified}}{
#'     The use of this method ignores stratification of the \code{J} tables,
#'     essentially turning the data into a simple 2 x 2 table for
#'     calculation of sample size.
#'   }
#' }
#'
#' @section Arguments:
#' The arguments required differs based on the calculation method selected using
#' the \code{method} parameter. In all cases, \code{sig.level}, \code{power},
#' \code{alternative}, \code{t}, and \code{method} must be specified or the
#' default can be used.
#'
#' The remaining arguments may or may not be used based on the \code{method}
#' selected. See list below:
#'
#' \describe{
#'
#'   \item{\code{cc.binomial},\code{binomial}, and \code{unstratified}}{
#'
#'     The user must specify \code{s} in addition to \strong{one} of the
#'     following:
#'     \itemize{
#'       \item Both case and control proportion vectors, ex.,
#'         \itemize{\item \code{p1} and \code{p2} with \code{theta = NULL}.}
#'       \item One proportion vector and an effect size, ex.,
#'         \itemize{
#'           \item \code{p1} and \code{theta} with \code{p2 = NULL}, \strong{or}
#'           \item \code{p2} and \code{theta} with \code{p1 = NULL}.
#'       }
#'     }
#'   }
#'
#'   \item{\code{hypergeometric}}{
#'     The user must specify \code{r} and \code{theta}.
#'   }
#'
#'   \item{\code{least.squares}}{
#'     The user must specify \strong{one} of the following:
#'     \itemize{
#'       \item Both case and control proportion vectors, ex.,
#'         \itemize{\item \code{p1} and \code{p2} with \code{theta = NULL}.}
#'       \item One proportion vector and an effect size, ex.,
#'         \itemize{
#'           \item \code{p1} and \code{theta} with \code{p2 = NULL}, \strong{or}
#'           \item \code{p2} and \code{theta} with \code{p1 = NULL}.
#'       }
#'     }
#'   }
#' }
#'
#' @param p1 Vector of proportions of the \emph{J} case groups.
#' @param p2 Vector of proportions of the \emph{J} control groups.
#' @param theta Vector of odds ratios relating to the \emph{J} 2 x 2 tables.
#' @param sig.level Significance level (Type I error probability).
#' @param power Power of test (1 minus Type II error probability).
#' @param alternative Two- or one-sided test. Can be abbreviated.
#' @param r Proportion of exposed versus not exposed in \emph{J} stratum.
#'   (Hypergeometric method only)
#' @param s Proportion of case versus control in \emph{J} stratum.
#' @param t Proportion of total number of cases in \emph{J} stratum.
#' @param method Method for calculation of sample size. Can be abbreviated. (See
#'   \strong{Methods} section.)
#'
#' @return An object of class \code{"samplesize.cmh"}: a list of the arguments
#' (except theta, which is converted to the missing proportion if supplied)
#' as well as effective N, exact N, a vector of n's per each group,
#' a short description of the calculation method used, the original function
#' call, and a note.
#'
#' @seealso
#' \link[stats]{power.prop.test},
#' \link[stats]{mantelhaen.test},
#' \link{power.cmh.test},
#' \link[DescTools]{BreslowDayTest}
#'
#' @author Paul W. Egeler, M.S.
#'
#' @examples
#' # From "Sample size determination for case-control studies and the comparison
#' # of stratified and unstratified analyses", Nam (1992). See references.
#'
#' # Uncorrected sample size estimate first introduced
#' # by Woolson and others in 1986
#' sample_size_uncorrected <- samplesize.cmh.test(
#'   p2 = c(0.75,0.70,0.65,0.60),
#'   theta = 3,
#'   power = 0.9,
#'   t = c(0.10,0.40,0.35,0.15),
#'   alternative = "one",
#'   method = "bin"
#' )
#'
#' sample_size_uncorrected
#'
#' # We see that the N.exact is 171, the same as calculated by Nam
#' sample_size_uncorrected$N.exact
#'
#'
#' # Continuity corrected sample size estimate added by Nam
#' sample_size_corrected <- samplesize.cmh.test(
#'   p2 = c(0.75,0.70,0.65,0.60),
#'   theta = 3,
#'   power = 0.9,
#'   t = c(0.10,0.40,0.35,0.15),
#'   alternative = "one",
#'   method = "cc.bin"
#' )
#'
#' sample_size_corrected
#'
#' # We see that the N.exact is indeed equal to that which is reported in the paper
#' sample_size_corrected$N.exact
#'
#' @references
#' Gail, M. (1973). "The determination of sample sizes for trials involving
#' several 2 x 2 tables." \emph{Journal of Chronic Disease} \strong{26}: 669-673.
#'
#' Munoz, A. and B. Rosner. (1984). "Power and sample size for a collection of 2
#' x 2 tables." \emph{Biometrics} \strong{40}: 995-1004.
#'
#' Nam, J. (1992). "Sample size determination for case-control studies and the
#' comparison of stratified and unstratified analyses." \emph{Biometrics} \strong{48}: 389-395.
#'
#' Wittes, J. and S. Wallenstein. (1987). "The power of the Mantel-Haensel test."
#' \emph{Journal of the American Statistical Association} \strong{82}:
#' 1104-1109.
#'
#' Woolson, R. F., Bean, J. A., and P. B. Rojas. (1986).
#' "Sample size for case-control studies using Cochran's statistic."
#' \emph{Biometrics} \strong{42}: 927-932.
samplesize.cmh.test <- function(
  p1 = NULL,
  p2 = NULL,
  theta = NULL,
  sig.level = 0.05,
  power = 0.80,
  alternative = c("two.sided","one.sided"),
  r = NULL,
  s = 0.5,
  t,
  method = c(
    "cc.binomial",
    "binomial",
    "hypergeometric",
    "least.squares",
    "unstratified"
    )
  ) {

  # Determine the method of calculation to use
  method <- match.arg(method)
  methods <- c(
    "cc.binomial" = "Continuity corrected weighted difference between two binomial distributions",
    "binomial" = "Weighted difference between two binomial distributions (Case-Control)",
    "hypergeometric" = "Hypergeometric distribution with both margins fixed",
    "least.squares" = "Weighted least squares",
    "unstratified" = "Unstratified (ignoring confounding variable)"
  )

  # Determine upper, lower, or two-sided hypothesis
  alternative <- match.arg(alternative)
  tside <- switch(alternative, two.sided = 2, 1)

  # Process the expected proportions and/or effect size
  if (method == "hypergeometric" &&
      is.null(r) ||
      (is.null(theta) && (is.null(p1) && is.null(p2)))
  ) {
    stop("'r' and 'theta' must be specified for hypergeometric method")
    # NOTE: 'p1' and 'p2' can be used in place of 'theta' since they can be used
    #       to determine 'theta' using props2theta()
  } else if (sum(sapply(list(p1,p2,theta),is.null)) != 1L) {
    stop("exactly one of 'p1', 'p2', or 'theta' must be NULL")
  }

  # Infer J from length of t
  J <- length(t)

  # Check that 'r', 's', and 't' are reasonable
  if (method != "least.squares" &&
      !isTRUE(all.equal(s + (1 - s), rep(1,length(s))))) {
    stop("The variable 's' must be a decimal fraction")
  }
  if (method == "hypergeometric" &&
      !isTRUE(all.equal(r + (1 - r), rep(1,length(r))))
  ) {
    stop("The variable 'r' must be a decimal fraction")
  }
  if (!isTRUE(all.equal(sum(t),1))) {
    stop("The 't' levels must sum to 1")
  }

  # Determine missing p vector or theta
  if (is.null(theta)) {
    theta <- props2theta(p1,p2)
  } else if (method != "hypergeometric") {
    if (is.null(p1)) {
      p1 <- effect.size(p2,theta)
    } else {
      p2 <- effect.size(p1,theta)
    }
  }

  if (!isTRUE(all.equal(min(theta), max(theta)))) {
    warning("'theta' differs among groups. This may violate assumptions.")
  }

  # Calculate z values for alpha and beta
  z_a <- stats::qnorm(sig.level / tside)
  z_b <- stats::qnorm(1 - power)

  # Set up the different calcualtions based on method
  calculations <- list(
    "cc.binomial" = quote({
      N_uncorrected <- eval(calculations[["binomial"]])

      (1 + sqrt(1 + 2 / (Z * N_uncorrected) ))^2 * N_uncorrected / 4
    }),
    "binomial" = quote({
      pbar <- p1 * s + p2 * (1 - s)

      X <- sum(t * s * (1 - s) * pbar * (1 - pbar))

      Y <- sum(t * s * (1 - s) * ((1 - s) * p1 * (1 - p1) + s * p2 * (1 - p2)))

      Z <- sum(t * s * (1 - s) * (p1 - p2))

      (z_a * sqrt(X) + z_b * sqrt(Y))^2 / Z^2
    }),
    "hypergeometric" = quote({
      h <- (2 * r * s)^-1 * (((theta - 1)^-1 + r + s) -
           sqrt(((theta - 1)^-1 + r + s)^2 - 4 * r * s * (theta - 1)^-1 * theta))

      h_prime <- theta * (1 - r * h) * (1 - s * h) /
                 (1 + (theta - 1) * (r + s - 2 * r * s * h))

      (z_a * sqrt(sum(t * r * s * (1- r) * (1 - s))) +
          z_b * sqrt(sum(t * r * s * h_prime)))^2 /
      (sum(t * r * s * (h - 1)))^2

    }),
    "least.squares" =  quote({
      if (!isTRUE(all.equal(s,rep(0.5,J)))) {
        warning(
          "'least.squares' method assumes s = 0.5\n",
          "user input 's' will be ignored"
        )
      }

        sum_fg <-
          t *
          log(theta)^2 /
          (1/(p1 * (1 - p1))+ 1/(p2 * (1 - p2)))

        2(z_a + z_b)^2 / sum_fg

    }),
    "unstratified" =  quote({
      sum_ts <- sum(t * s)

      sum_t1ms <- sum(t * (1 - s))

      p2p <- sum(t * (1 - s) * p2 / sum_t1ms)

      p1p <- sum(t * s * p1 / sum_ts) # Note to confirm

      ppp <- p1p * sum_ts + p2p * sum_t1ms

      (z_a * sqrt(sum_ts * sum_t1ms * ppp * (1 - ppp)) +
       z_b * sqrt(sum_ts^2 * sum_t1ms^2 *
      (p1p * (1 - p1p) / sum_ts + p2p * (1 - p2p) / sum_t1ms)))^2 /
      ((p1p - p2p) * sum_ts * sum_t1ms)^2

    })
  )

  # Run the calculation
  N <- eval(calculations[[method]])

  # Return an object of class "samplesize.cmh"
  structure(
    list(
      n1 = ceiling(rep(N,J)*t*s),
      n2 = ceiling(rep(N,J)*t*(1 - s)),
      N.effective = sum(ceiling(rep(N,J)*t*s), ceiling(rep(N,J)*t*(1 - s))),
      N = N, p1 = p1, p2 = p2, sig.level = sig.level,
      power = power, alternative = alternative, r = r, s = s, t = t, J = J,
      note = paste(
        "N is the calculated *total* number of subjects;",
        "Effective N is the minimum number of subjects to satisfy whole number",
        "requirement in each cell",
        sep = "\n"
      ),
      method = method, method.desc = methods[method],
      call = match.call()
    ),
    class = "samplesize.cmh"
    )
}

# Print method so that "samplesize.cmh" will look nice
#' @export
print.samplesize.cmh <- function(x, ...) {
  with(x, cat(
    "Sample size calculation for the Cochran Mantel Haenszel test\n\n",

    "                 N = ",N,"\n",
    "       Effective N = ",N.effective,"\n",
    "Significance level = ",sig.level,"\n",
    "     Nominal Power = ",power,"\n",
    "             Sides = ",alternative,"\n\n",

    "Number of subjects per each group:\n",
    rep("_",10 + J * 7),"\n",
    "Group   |",sprintf(" %5i ",1:J),"\n",
    rep("=",10 + J * 7),"\n",
    "Case    |",sprintf(" %5i ",n1),"\n",
    "Control |",sprintf(" %5i ",n2),"\n\n",

    "METHOD: ",method.desc,"\n",
    "NOTE  : ",note,
    sep = ""
    )
  )

  invisible(x)

}
