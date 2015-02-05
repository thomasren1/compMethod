#' tl method
#'
#' Calculate the order statistic for the a% non-parametric one-sided
#' tolerance limit with b% confidence
#'
#' @param mor is a vector of MOR data for all species
#' @param a is the percent of the population below the tolerance limit
#' @param b is the percent confidence of the tolerance limit
#' @return tolerance limit
#'
#' @export
#'
#'
tl <- function(mor, a = 5, b = 75, ...)
{

	N <- length(mor)
	b1 <- 0:N
	## Cumulative binomial probabilities for each value of b1 = 0,...,N
	b2 <- pbinom(b1, N, a / 100)
	## The index r of the order statistic so that at least b% of the
	## time the tolerance limit will lie below the value of interest.
	r <- max(b1[b2 < 1 - b / 100]) + 1
	return(sort(mor)[r])
}
