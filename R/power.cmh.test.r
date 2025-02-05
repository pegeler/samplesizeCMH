#' @title Power and sample size calculation for the Cochran-Mantel-Haenszel test
#'
#' @description
#' Compute the post-hoc power or required number of subjects for the
#' Cochran-Mantel-Haenszel test for association in \emph{J} stratified 2 x 2
#' tables.
#'
#' @details
#' This sample size calculation is based on the derivations described in the
#' Woolson \emph{et al.} (1986). It is designed for case-control studies where
#' one margin is fixed. The method is "based on the Cochran-Mantel-Haenszel
#' statistic expressed as a weighted difference in binomial proportions."
#'
#' Continuity corrected sample size is described in Nam's 1992 paper. This uses
#' the weighted binomial sample size calculation described in Woolson \emph{et
#' al.} (1986) but is enhanced for use with the continuity corrected Cochran's
#' test.
#'
#' Power calculations are based on the writings of Wittes and Wallenstein
#' (1987). They are functionally equivalent to the derivations of the sample
#' size calculation described by Woolson and others and Nam, but have slightly
#' added precision.
#'
#' Terminology and symbolic conventions are borrowed from Woolson \emph{et al.}
#' (1986). The \code{p1} group is dubbed the \emph{Case} group and \code{p2}
#' group is called the \emph{Control} group.
#'
#' @section Arguments:
#'
#' To calculate \strong{power}, the \code{power} parameter must be set to
#' \code{NULL}. To calculate \strong{sample size}, the \code{N} parameter must
#' be set to \code{NULL}.
#'
#' The \code{J} number of groups will be inferred by the maximum length of
#' \code{p1}, \code{p2}, or \code{theta}.
#'
#' Effect size must be specified using one of the following combinations of
#' arguments.
#' \itemize{
#'   \item Both case and control proportion vectors, ex.,
#'     \itemize{\item \code{p1} and \code{p2} with \code{theta = NULL}.}
#'   \item One proportion vector and an effect size, ex.,
#'     \itemize{
#'       \item \code{p1} and \code{theta} with \code{p2 = NULL}, \strong{or}
#'       \item \code{p2} and \code{theta} with \code{p1 = NULL}.
#'     }
#' }
#'
#' @param p1 Vector of proportions of the \emph{J} case groups.
#' @param p2 Vector of proportions of the \emph{J} control groups.
#' @param theta Vector of odds ratios relating to the \emph{J} 2 x 2 tables.
#' @param N Total number of subjects.
#' @param sig.level Significance level (Type I error probability).
#' @param power Power of test (1 minus Type II error probability).
#' @param alternative Two- or one-sided test. If one-sided, the direction of the
#'   association must be defined (less than 1 or greater than 1). Can be
#'   abbreviated.
#' @param s Proportion (weight) of case versus control in \emph{J} stratum.
#' @param t Proportion (weight) of total number of cases of \emph{J} stratum.
#' @param correct Logical indicating whether to apply continuity correction.
#'
#' @return
#' An object of class \code{"power.cmh"}: a list of the original arguments and
#' the calculated sample size or power. Also included are vectors of n's per
#' each group, an indicator or whether continuity correction was used, the
#' original function call, and \code{N.effective}.
#'
#' The vectors of n's per each group, \code{n1} and \code{n2}, are the
#' fractional n's required to achieve a final total N specified by the
#' calculation while satisfying the constraints of \code{s} and \code{t}.
#' However, the effective N, given the requirement of cell counts populated by
#' whole numbers is provided by \code{N.effective}. By default, the print method
#' is set to \code{n.frac = FALSE}, which will round each cell n up to the
#' nearest whole number.
#'
#' @seealso
#' \link[stats]{power.prop.test},
#' \link[stats]{mantelhaen.test},
#' \link[DescTools]{BreslowDayTest}
#'
#' @author Paul W. Egeler, M.S.
#'
#' @examples
#' # From "Sample size determination for case-control studies and the comparison
#' # of stratified and unstratified analyses", (Nam 1992). See references.
#'
#' # Uncorrected sample size estimate first introduced
#' # by Woolson and others in 1986
#' sample_size_uncorrected <- power.cmh.test(
#'   p2 = c(0.75, 0.70, 0.65, 0.60),
#'   theta = 3,
#'   power = 0.9,
#'   t = c(0.10, 0.40, 0.35, 0.15),
#'   alternative = "greater",
#'   correct = FALSE
#' )
#'
#' print(sample_size_uncorrected, detail = FALSE)
#'
#' # We see that the N is 171, the same as calculated by Nam
#' sample_size_uncorrected$N
#'
#'
#' # Continuity corrected sample size estimate added by Nam
#' sample_size_corrected <- power.cmh.test(
#'   p2 = c(0.75, 0.70, 0.65, 0.60),
#'   theta = 3,
#'   power = 0.9,
#'   t = c(0.10, 0.40, 0.35, 0.15),
#'   alternative = "greater",
#'   correct = TRUE
#' )
#'
#' print(sample_size_corrected, n.frac = TRUE)
#'
#' # We see that the N is indeed equal to that which is reported in the paper
#' sample_size_corrected$N
#'
#' @references
#' Gail, M. (1973). "The determination of sample sizes for trials involving
#' several 2 x 2 tables."
#' \emph{Journal of Chronic Disease} \strong{26}: 669-673.
#'
#' Munoz, A. and B. Rosner. (1984). "Power and sample size for a collection of 2
#' x 2 tables." \emph{Biometrics} \strong{40}: 995-1004.
#'
#' Nam, J. (1992). "Sample size determination for case-control studies and the
#' comparison of stratified and unstratified analyses."
#' \emph{Biometrics} \strong{48}: 389-395.
#'
#' Wittes, J. and S. Wallenstein. (1987). "The power of the Mantel-Haenszel
#' test." \emph{Journal of the American Statistical Association} \strong{82}:
#' 1104-1109.
#'
#' Woolson, R. F., Bean, J. A., and P. B. Rojas. (1986).
#' "Sample size for case-control studies using Cochran's statistic."
#' \emph{Biometrics} \strong{42}: 927-932.
#'
#' @export
power.cmh.test <- function(
  p1 = NULL,
  p2 = NULL,
  theta = NULL,
  N = NULL,
  sig.level = 0.05,
  power = 0.80,
  alternative = c("two.sided", "less", "greater"),
  s = 0.5,
  t = 1 / J,
  correct = TRUE
  ) {

  # Determine if the correct number of arguments are set to NULL
  if (count_nulls(p1, p2, theta) != 1L) {
    stop("exactly one of 'p1', 'p2', or 'theta' must be NULL")
  }

  if (count_nulls(N, power) != 1L) {
    stop("either 'N' or 'power' must be NULL")
  }

  # Determine upper, lower, or two-sided hypothesis
  alternative <- match.arg(alternative)
  tside <- switch(alternative, two.sided = 2, 1)
  lower.tail <- switch(alternative, less = TRUE, FALSE)

  # Infer J from vector lengths of first three args
  J <- max(lengths(list(p1, p2, theta)))

  # Ensure that 's' and 't' are of correct lengths
  s <- rep(as.vector(s), length.out = J)
  t <- rep(as.vector(t), length.out = J)

  # Check that 's' and 't' are reasonable
  if (any(s < 0) || any(s > 1)) {
    stop("The variable 's' must be a decimal fraction")
  }
  if (!isTRUE(all.equal(sum(t), 1))) {
    stop("The 't' levels must sum to 1")
  }

  # Determine missing 'p' or 'theta' vector
  if (is.null(theta)) {
    theta <- props2theta(p1, p2)
  } else if (is.null(p1)) {
    p1 <- effect.size(p2, theta)
  } else {
    p2 <- effect.size(p1, theta)
  }

  if (all(theta == 1)) {
    stop(
      "Effect size ('theta') is 1 across all strata.\n",
      "Calculation cannot be performed.\n"
    )
  }
  if (!all(theta >= 1) && !all(theta <= 1)) {
    message("NOTE: Effect size ('theta') differs in direction across strata.")
  }
  if (!isTRUE(all.equal(min(as.vector(theta)), max(as.vector(theta))))) {
    message("NOTE: Effect size ('theta') is not consistent between strata.")
  }

  # Calculate z value for alpha
  z_a <- stats::qnorm(sig.level / tside, lower.tail = lower.tail)

  # Determine if N or power is missing and then complete the calculation
  if (is.null(N)) {  # Get sample size
    z_b <- stats::qnorm(1 - power, lower.tail = lower.tail)

    pbar <- p1 * s + p2 * (1 - s)

    X <- sum(t * s * (1 - s) * pbar * (1 - pbar))

    Y <- sum(t * s * (1 - s) * ((1 - s) * p1 * (1 - p1) + s * p2 * (1 - p2)))

    Z <- sum(t * s * (1 - s) * (p1 - p2))

    # The calculation itself
    N <- (z_a * sqrt(X) + z_b * sqrt(Y))^2 / Z^2

    # Continuity correction
    if (correct) {
      if (!lower.tail)
        N <- (1 + sqrt(1 + 2 / (abs(Z) * N) ))^2 * N / 4
      else
        N <- (1 + sqrt(1 - 2 / (Z * N) ))^2 * N / 4
    }

    # Splitting N into groups
    n1 <- N*t*s
    n2 <- N*t*(1 - s)

    # Calculating effective N
    N.effective <- sum(ceiling(c(n1, n2)))

  } else {  # Get power

    std_dev <- sqrt(
      sum(t * s * (1 - s) * ((1 - s) * p1 * (1 - p1) + s * p2 * (1 - p2)))
    )

    pbar <- p1 * s + p2 * (1 - s)

    phi <- pbar * (1 - pbar) + (p1 - p2)^2 * s * (1 - s) / (N * t - 1)

    Eg <- sum(N * t * s * (1 - s) * (p1 - p2))
    numerator <- (
      ifelse(alternative == "two.sided", abs(Eg), Eg) +
      ifelse(!correct, 0, ifelse(alternative == "less", 0.5, -0.5))) / sqrt(N) -
      z_a * sqrt(sum(t * s * (1 - s) * phi))

    U <- numerator / std_dev

    power <- stats::pnorm(U, lower.tail = !lower.tail)

    # Splitting N into groups
    n1 <- round(N*t*s, 0)
    n2 <- round(N*t*(1 - s), 0)

    # N.effective will just be N
    N.effective <- N

  }

  # Return an object of class "power.cmh"
  structure(
    list(
      p1 = p1, p2 = p2, theta = theta, N = N, sig.level = sig.level,
      power = power, alternative = alternative, s = s, t = t, J = J,
      n1 = n1, n2 = n2, N.effective = N.effective, call = match.call()
    ),
    class = "power.cmh"
    )
}

#' Print \code{"power.cmh"} object
#'
#' The S3 print method for the \code{"power.cmh"} object
#'
#' @param x A \code{"power.cmh"} object.
#' @param detail Logical to toggle detailed or simple output.
#' @param n.frac Logical indicating whether sample n's should be rounded to the
#'   next whole number.
#' @param ... Ignored.
#' @export
print.power.cmh <- function(x, detail = TRUE, n.frac = FALSE, ...) {
  with(x, {
    cat(
      "Power and sample size calculation for the Cochran Mantel Haenszel test\n\n",

      "                 N = ", ifelse(n.frac, N, ceiling(N)), "\n",
      "       Effective N = ", N.effective, "\n",
      "Significance level = ", sig.level, "\n",
      "             Power = ", power, "\n",
      "       Alternative = ", alternative, "\n\n",
      sep = ""
    )

    if (detail) {

      if (n.frac) {
        cat(
          "Number of subjects per each group:\n",
          rep("_", 10 + J * 9), "\n",
          "Group   |", sprintf(" %7i ", 1:J), "\n",
          rep("=", 10 + J * 9), "\n",
          "Case    |", sprintf(" %7.2f ", n1), "\n",
          "Control |", sprintf(" %7.2f ", n2), "\n\n",
          sep = ""
        )
      } else {
        cat(
          "Number of subjects per each group:\n",
          rep("_", 10 + J * 7), "\n",
          "Group   |", sprintf(" %5i ", 1:J), "\n",
          rep("=", 10 + J * 7), "\n",
          "Case    |", sprintf(" %5i ", ceiling(n1)), "\n",
          "Control |", sprintf(" %5i ", ceiling(n2)), "\n\n",
          sep = ""
        )
      }
    }

    cat(
      "CALL: \n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = ""
    )

  })

  invisible(x)

}

count_nulls <- function(...) sum(vapply(list(...), is.null, logical(1)))
