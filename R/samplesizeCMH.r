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
#'   \item{\code{binomial}}{
#'     This sample size calculation is designed for case-control studies
#'     where one margin is fixed. Described in the Woolson \emph{et al.} (1986).
#'   }
#'
#'   \item{\code{hypergeometric}}{
#'     Pass
#'   }
#'
#'   \item{\code{least.squares}}{
#'     Pass
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
#' a short description of the calculation method used, and a note.
#'
#' @seealso \link{power.prop.test}, \link{mantelhaen.test}
#'
#' @author Paul W. Egeler, M.S.
#'
#' @examples
#' ## Stata power cmh manual
#' ## Should make our own examples
#' #  Example 1
#' samplesizeCMH(
#' 	p1 = NULL,
#' 	p2 = c(.426,.444,.364),
#' 	theta = 2.5,
#' 	sig.level = 0.05,
#' 	power = .8,
#' 	alternative = "two"
#' )
#'
#' #  Example 2
#' samplesizeCMH(
#' 	p1 = NULL,
#' 	p2 = c(.426,.444,.364),
#' 	theta = 2.5,
#' 	sig.level = 0.05,
#' 	power = .8,
#' 	alternative = "two",
#' 	t = c(4,1,4) / sum(c(4,1,4))
#' )
#'
#' #  Example 3a
#' samplesizeCMH(
#' 	p1 = NULL,
#' 	p2 = c(.426,.444,.364),
#' 	theta = 2.5,
#' 	sig.level = 0.05,
#' 	power = .8,
#' 	alternative = "two",
#' 	t = c(4,1,4) / sum(c(4,1,4)),
#' 	s = 1 - c(.47,.57,.51)
#' )
#'
#' #  Example 3b
#' samplesizeCMH(
#' 	p1 = NULL,
#' 	p2 = c(.426,.444,.364),
#' 	theta = 2.5,
#' 	sig.level = 0.05,
#' 	power = .8,
#' 	alternative = "two",
#' 	t = c(4,1,4) / sum(c(4,1,4)),
#' 	s = 1 - c(.8,.7,.3)
#' )
#'
#' @source
#' Woolson, R. F., J. A. Bean, and P. B. Rojas. 1986.
#' Sample size for case-control studies using Cochran's statistic.
#' Biometrics 42: 927-932. PMID 3814733.
#'
#' Add Gail(1973) and Munoz and Rosner (1984)
#'
#' @export
samplesizeCMH <- function(
	p1 = NULL,
	p2 = NULL,
	theta = NULL,
	sig.level = 0.05,
	power = 0.80,
	alternative = c("two.sided","one.sided"),
	s = 0.5,
	t = 1 / J,
	method = c("binomial","hypergeometic","least.squares","unstratified")
	) {

	# Process the expected proportions and/or effect size
	if (sum(sapply(list(p1,p2,theta),is.null)) != 1) {
		stop("exactly one of 'p1', 'p2', or 'theta' must be NULL")
	}

	# Infer J from vector lengths of first three args
	J <- max(sapply(list(p1,p2,theta),length))

	# If theta is used, determine the missing p vector
	if (is.null(p1)) {
		p1 <- (p2 * theta) / (1 - p2 + p2 * theta)
	}	else if (is.null(p2)) {
		p2 <- (p1 * theta) / (1 - p1 + p1 * theta)
	}

	# Determine upper, lower, or two-sided hypothesis
	alternative <- match.arg(alternative)
	tside <- switch(alternative, two.sided = 2, 1)

	# Calculate z values for alpha and beta
	z_a <- qnorm(sig.level / tside)
	z_b <- qnorm(1 - power)

	# Determine the method of calculation to use
	method <- match.arg(method)
	methods <- c(
		"binomial" = "Weighted difference between two binomial distributions (Case-Control)",
		"hypergeometic" = "Hypergeometric distribution with both margins fixed",
		"least.squares" = "Weighted least squares",
		"unstratified" = "Unstratified (ignoring confounding variable)"
	)

	# Set up the different calcualtions based on method
	calculations <- list(
		"binomial" = quote({
			pbar <- p1 * s + p2 * (1 - s)

			X <- sum(t * s * (1 - s) * pbar * (1 - pbar))

			Y <- sum(t * s * (1 - s) * ((1 - s) * p1 * (1 - p1) + s * p2 * (1 - p2)))

			Z <- sum(t * s * (1 - s) * (p1 - p2))

			(z_a * sqrt(X) + z_b * sqrt(Y))^2 / Z^2
		}),
		"hypergeometric" = quote({
			999 # Need to add
		}),
		"least.squares" =  quote({
			888 # Need to add
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
			method = method, method.desc = methods[method]),
		class = "samplesize.cmh"
		)
}

# Print method so that "samplesize.cmh" will look nice
#' @export
print.samplesize.cmh <- function(x) {
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
}
