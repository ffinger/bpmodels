##' Calculates the complementary log-probability
##'
##' Given x and norm, this calculates log(1-sum(exp(x)))
##' @param x log-probabilities
##' @return value
##' @author Sebastian Funk
##' @keywords internal
complementary_logprob <- function(x) {
    return(log1p(max(-sum(exp(x)), -1)))
}
