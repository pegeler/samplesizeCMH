#' Interconvert odds and proportion values
#'
#' These functions will create either odds for a given probability,
#' probability for a given odds, calculate the odds ratio between two
#' probabilities, or calculate effect size (raise a probability by theta)
#'
#' @param p,p1,p2 Proportion vector.
#' @param o Odds vector.
#' @param theta Odds ratio vector.
#' @param rr Relative risk vector (\code{p1 / p2}).
#' @examples
#' # Convert proportions of 0 through 1 to odds
#' props <- seq(0,1,0.1)
#' prop2odds(props)
#'
#' # Convert odds to proportions
#' odds2prop(1:3)
#'
#' # Raise a proportion by an effect size theta
#' effect.size(0.5, 2)
#'
#' # Find the odds ratio between two proportions
#' props2theta(0.75, 0.5)
#' @return A numeric vector.
#' @author Paul W. Egeler, M.S.
#' @name odds.and.proportions
NULL


#' @rdname odds.and.proportions
#' @export
prop2odds <- function(p) {p / (1 - p)}

#' @rdname odds.and.proportions
#' @export
odds2prop <- function(o) {o / (1 + o)}

#' @rdname odds.and.proportions
#' @export
effect.size <- function(p,theta) {(p * theta) / (1 - p + p * theta)}

#' @rdname odds.and.proportions
#' @export
props2theta <- function(p1,p2) {p1 * (1 - p2) / p2 / (1 - p1)}


# Should include? ---------------------------------------------------------

#' @rdname odds.and.proportions
#' @export
rr2theta <- function(rr,p1,p2) {rr * (1 - p2) / (1 - p1)}

#' @rdname odds.and.proportions
#' @export
theta2rr <- function(theta,p1,p2) {theta * (1 - p1) / (1 - p2)}
