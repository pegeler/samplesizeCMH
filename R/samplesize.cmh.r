#' @title Sample size calculation for the Cochran Mantel Haenszel test
#'
#' @description
#' Computes the sample size for the test for assocaition between
#' \emph{J 2 x 2} tables. The calculation can be performed using various
#' methods found in literature.
#'
#' @details
#' By convention, the \code{p1} group is dubbed the \emph{Case} group and \code{p2}
#' group is called the \emph{Control} group. This is carry-over language from the
#' Case-Control terminology used by Woolson \emph{et al.} (1986).
#'
#' You must specify one of the following:
#' \itemize{
#'   \item Both case and control proportion vectors, ex.,
#'     \itemize{\item \code{p1} and \code{p2} with \code{theta = NULL}.}
#'   \item One proportion vector and an effect size, ex.,
#'     \itemize{
#'       \item \code{p1} and \code{theta} with \code{p2 = NULL}, \strong{or}
#'       \item \code{p2} and \code{theta} with \code{p1 = NULL}.
#'     }
#' }
#' The \code{J} number of groups will be inferred by the maximum length of
#' \code{p1}, \code{p2}, or \code{theta}.
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
#'     \strong{Under construction}.
#'     Described in Munoz and Rosner (1984). Woolson \emph{et al.} (1986)
#'     state that this is the appropriate sample size estimator when both
#'     margins are fixed.
#'   }
#'
#'   \item{\code{least.squares}}{
#'     \strong{Under construction}.
#'     Described in Gail (1973). This may not be relevant, as subsequent
#'     Munoz and Rosner (1984), Woolson \emph{et al.} (1986), and later
#'     Nam (1992) described "improved" methods.
#'   }
#'
#'   \item{\code{unstratified}}{
#'     The use of this method ignores stratification of the \code{J} tables,
#'     essentially turning the data into a simple 2 x 2 table for
#'     calculation of sample size.
#'   }
#' }
#'
#' @param p1 Vector of proportions of the \emph{J} case groups.
#' @param p2 Vector of proportions of the \emph{J} control groups.
#' @param theta Vector of odds ratios relating to the \emph{J 2 x 2} tables.
#' @param sig.level Significance level (Type I error probability).
#' @param power Power of test (1 minus Type II error probability).
#' @param alternative Two- or one-sided test. Can be abbreviated.
#' @param s Proportion of case versus control in \emph{J} stratum.
#' @param t Proportion of total number of cases in \emph{J} stratum.
#' @param method Method for calculation of sample size.
#'  Can be abbreviated. (See \strong{Methods} section.)
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
#' \link{power.cmh},
#' \link[DescTools]{BreslowDayTest}
#'
#' @author Paul W. Egeler, M.S.; Laura Kapitula, PhD
#'
#' @examples
#' # From "Sample size determination for case-control studies and the comparison
#' # of stratified and unstratified analyses", Nam (1992). See references.
#'
#' # Uncorrected sample size estimate first introduced
#' # by Woolson and others in 1986
#' sample_size_uncorrected <- samplesize.cmh(
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
#' sample_size_corrected <- samplesize.cmh(
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
#'
#' @export
samplesize.cmh <- function(
  p1 = NULL,
  p2 = NULL,
  theta = NULL,
  sig.level = 0.05,
  power = 0.80,
  alternative = c("two.sided","one.sided"),
  s = 0.5,
  t = 1 / J,
  method = c(
    "cc.binomial",
    "binomial",
    "hypergeometric",
    "least.squares",
    "unstratified"
    )
  ) {

  # Process the expected proportions and/or effect size
  if (sum(sapply(list(p1,p2,theta),is.null)) != 1) {
    stop("exactly one of 'p1', 'p2', or 'theta' must be NULL")
  }

  # Infer J from vector lengths of first three args
  J <- max(sapply(list(p1,p2,theta),length))

  #Ensure that s and t are of correct lengths
  s <- rep(s, length.out = J)
  t <- rep(t, length.out = J)

  # Check that 's' and 't' are reasonable
  if (any(s >= 1)) stop("The variable 's' must be a decimal fraction")
  if (!isTRUE(all.equal(sum(t),1))) stop("The 't' levels must sum to 1")

  # If theta is used, determine the missing p vector
  if (is.null(p1)) {
    p1 <- (p2 * theta) / (1 - p2 + p2 * theta)
  }  else if (is.null(p2)) {
    p2 <- (p1 * theta) / (1 - p1 + p1 * theta)
  }

  # Determine upper, lower, or two-sided hypothesis
  alternative <- match.arg(alternative)
  tside <- switch(alternative, two.sided = 2, 1)

  # Calculate z values for alpha and beta
  z_a <- stats::qnorm(sig.level / tside)
  z_b <- stats::qnorm(1 - power)

  # Determine the method of calculation to use
  method <- match.arg(method)
  methods <- c(
    "cc.binomial" = "Continuity corrected weighted difference between two binomial distributions",
    "binomial" = "Weighted difference between two binomial distributions (Case-Control)",
    "hypergeometric" = "Hypergeometric distribution with both margins fixed",
    "least.squares" = "Weighted least squares",
    "unstratified" = "Unstratified (ignoring confounding variable)"
  )

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
      stop("Method under construction")
    }),
    "least.squares" =  quote({
      stop("Method under construction")
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
      N = sum(ceiling(rep(N,J)*t*s), ceiling(rep(N,J)*t*(1 - s))),
      N.exact = N, p1 = p1, p2 = p2, sig.level = sig.level,
      power = power, alternative = alternative, s = s, t = t, J = J,
      note = "N is *total* number of subjects",
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
