#' @title Power and sample size calculation for the Cochran-Mantel-Haenszel test
#'
#' @description
#' Compute the post-hoc power or required number of subjects for the
#' Cochran-Mantel-Haenszel test for association in \emph{J} stratified 2 x 2
#' tables. The calculations can be performed using various methods found in
#' literature.
#'
#' @details
#' Terminology and symbolic conventions are borrowed from Woolson \emph{et al.}
#' (1986). The \code{p1} group is dubbed the \emph{Case} group and \code{p2}
#' group is called the \emph{Control} group.
#'
#' Power is solved iteratively.
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
#'   \item{\code{unstratified}}{
#'     The use of this method ignores stratification of the \code{J} tables,
#'     essentially turning the data into a simple 2 x 2 table for
#'     calculation of sample size.
#'   }
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
#' @param method Method for calculation of sample size. Can be abbreviated. (See
#'   \strong{Methods} section.)
#'
#' @return
#' An object of class \code{"power.cmh"}: a list of the original arguments and
#' the calculated sample size or power. Also included are vectors of n's per
#' each group, a short description of the calculation method used, the original
#' function call, and \code{N.effective}.
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
#'   p2 = c(0.75,0.70,0.65,0.60),
#'   theta = 3,
#'   power = 0.9,
#'   t = c(0.10,0.40,0.35,0.15),
#'   alternative = "greater",
#'   method = "bin"
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
#'   p2 = c(0.75,0.70,0.65,0.60),
#'   theta = 3,
#'   power = 0.9,
#'   t = c(0.10,0.40,0.35,0.15),
#'   alternative = "greater",
#'   method = "cc.bin"
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
#' Wittes, J. and S. Wallenstein. (1987). "The power of the Mantel-Haensel
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
  alternative = c("two.sided","less","greater"),
  s = 0.5,
  t = 1 / J,
  method = c(
    "cc.binomial",
    "binomial",
    "unstratified"
    )
  ) {

  # Process the expected proportions and/or effect size
  if (sum(sapply(list(p1,p2,theta),is.null)) != 1L) {
    stop("exactly one of 'p1', 'p2', or 'theta' must be NULL")
  }

  if (sum(sapply(list(N, power),is.null)) != 1L) {
    stop("either 'N'or 'power' must be NULL")
  }

  # Determine the method of calculation to use
  method <- match.arg(method)
  methods <- c(
    "cc.binomial" = paste(
      "Continuity corrected weighted difference",
      "between two binomial distributions"),
    "binomial" =
      "Weighted difference between two binomial distributions (Case-Control)",
    "unstratified" =
      "Unstratified (ignoring confounding variable)"
  )

  # Determine upper, lower, or two-sided hypothesis
  alternative <- match.arg(alternative)
  tside <- switch(alternative, two.sided = 2, 1)
  lower.tail <- switch(alternative, less = TRUE, FALSE)

  # Infer J from vector lengths of first three args
  J <- max(sapply(list(p1,p2,theta),length))

  #Ensure that s and t are of correct lengths
  s <- rep(as.vector(s), length.out = J)
  t <- rep(as.vector(t), length.out = J)

  # Check that 's' and 't' are reasonable
  if (!isTRUE(all.equal(s + (1 - s), rep(1,J)))) {
    stop("The variable 's' must be a decimal fraction")
  }
  if (!isTRUE(all.equal(sum(t),1))) {
    stop("The 't' levels must sum to 1")
  }

  # Determine missing p vector or theta
  if (is.null(theta)) {
    theta <- props2theta(p1,p2)
  } else if (is.null(p1)) {
    p1 <- effect.size(p2,theta)
  } else {
    p2 <- effect.size(p1,theta)
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

  # Calculate z values for alpha and beta
  z_a <- stats::qnorm(sig.level / tside, lower.tail = lower.tail)

  # Set up the different calcualtions based on method
  calculations <-
    list(
      sample.size = list(
        "cc.binomial" = quote({
          N_uncorrected <- eval(calculations[["sample.size"]][["binomial"]])

          if (!lower.tail)
            (1 + sqrt(1 + 2 / (abs(Z) * N_uncorrected) ))^2 * N_uncorrected / 4
          else
            (1 + sqrt(1 - 2 / (Z * N_uncorrected) ))^2 * N_uncorrected / 4
        }),
        "binomial" = quote({
          pbar <- p1 * s + p2 * (1 - s)

          X <- sum(t * s * (1 - s) * pbar * (1 - pbar))

          Y <- sum(t * s * (1 - s) * ((1 - s) * p1 * (1 - p1) + s * p2 * (1 - p2)))

          Z <- sum(t * s * (1 - s) * (p1 - p2))

          (z_a * sqrt(X) + z_b * sqrt(Y))^2 / Z^2
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
      ),
      power = list(
        "cc.binomial" = quote({
          std_dev <- sqrt(
            sum(t * s * (1 - s) * ((1 - s) * p1 * (1 - p1) + s * p2 * (1 - p2)))
          )

          pbar <- p1 * s + p2 * (1 - s)

          phi <- pbar * (1 - pbar) + (p1 - p2)^2 * s * (1 - s) / (N * t - 1)

          numerator <-
            if (alternative == "two.sided") {
              abs(sqrt(N) * sum(t * s * (1 - s) * (p1 - p2)) -
                     z_a * sqrt(sum(t * s * (1 - s) * phi))) - 0.5
            } else {
              # Do I need to change sign for lesser?
              (sqrt(N) * sum(t * s * (1 - s) * (p1 - p2)) -
                  z_a * sqrt(sum(t * s * (1 - s) * phi))) - 0.5
            }

          pnorm(numerator / std_dev, lower.tail = !lower.tail)
        }),
        "binomial" = quote({
          std_dev <- sqrt(
            sum(t * s * (1 - s) * ((1 - s) * p1 * (1 - p1) + s * p2 * (1 - p2)))
            )

          pbar <- p1 * s + p2 * (1 - s)

          phi <- pbar * (1 - pbar) + (p1 - p2)^2 * s * (1 - s) / (N * t - 1)

          numerator <-
            if (alternative == "two.sided") {
              abs(sqrt(N) * sum(t * s * (1 - s) * (p1 - p2)) -
                     z_a * sqrt(sum(t * s * (1 - s) * phi)))
            } else {
              sqrt(N) * sum(t * s * (1 - s) * (p1 - p2)) -
                  z_a * sqrt(sum(t * s * (1 - s) * phi))
            }

          pnorm(numerator / std_dev, lower.tail = !lower.tail)
        }),
        "unstratified" =  quote({
          stop("Under construction")
        })
      )
    )

  # Run the calculation
  if (is.null(N)) {
    z_b <- stats::qnorm(1 - power, lower.tail = lower.tail)
    N <- eval(calculations[["sample.size"]][[method]])
  } else {
    power <- eval(calculations[["power"]][[method]])
  }

  # Return an object of class "power.cmh"
  structure(
    list(
      p1 = p1, p2 = p2, theta = theta, N = N, sig.level = sig.level,
      power = power, alternative = alternative, s = s, t = t, J = J,
      n1 = N*t*s, n2 = N*t*(1 - s),
      N.effective = sum(ceiling(rep(N,J)*t*s), ceiling(rep(N,J)*t*(1 - s))),
      method = method, method.desc = methods[method],
      call = match.call()
    ),
    class = "power.cmh"
    )
}

# Print method so that "samplesize.cmh" will look nice
#' @export
print.power.cmh <- function(x, detail = TRUE, n.frac = FALSE, ...) {
  with(x, {
    cat(
      "Power and sample size calculation for the Cochran Mantel Haenszel test\n\n",

      "                 N = ",ceiling(N),"\n",
      "       Effective N = ",N.effective,"\n",
      "Significance level = ",sig.level,"\n",
      "             Power = ",power,"\n",
      "       Alternative = ",alternative,"\n\n",
      sep = ""
    )

    if (detail) {

      if (n.frac) {
        cat(
          "Number of subjects per each group:\n",
          rep("_",10 + J * 9),"\n",
          "Group   |",sprintf(" %7i ",1:J),"\n",
          rep("=",10 + J * 9),"\n",
          "Case    |",sprintf(" %7.2f ",n1),"\n",
          "Control |",sprintf(" %7.2f ",n2),"\n\n",
          sep = ""
        )
      } else {
        cat(
          "Number of subjects per each group:\n",
          rep("_",10 + J * 7),"\n",
          "Group   |",sprintf(" %5i ",1:J),"\n",
          rep("=",10 + J * 7),"\n",
          "Case    |",sprintf(" %5i ",ceiling(n1)),"\n",
          "Control |",sprintf(" %5i ",ceiling(n2)),"\n\n",
          sep = ""
        )
      }
    }

    cat(
      "CALL: \n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n",
      "METHOD: ", method.desc, "\n", sep = ""
    )

  })

  invisible(x)

}
