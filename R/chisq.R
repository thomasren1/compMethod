#' Chi Square test
#'
#' determine the number of pieces in each species group
#' with MOR values below the combined group's tolerance limit value and
#' conduct a Chi Square test at a significance level of 0.01.
#'
#' See Appendix X7, page 21-22 of ASTM D1990-07.
#'
#' @param mor is a vector of MOR data for all species
#' @param spfac is the species factor for mor
#'
#' @return a list giving the controlling subgroup (csg) and the tolerance
#' limit characteristic value (cval)
#'
#'
#' @note Follow with iterative Chi Square tests on species subsets if needed.
#' @note no continuity correction is applied when computing the test
#' statistic which results in a warning message
#'
#' @export


chisq <- function(mor, spfac)
{

	## Combined group's tolerance limit value
	tlval <- tl(mor)

	## Determine the number of pieces in each species group below the combined
	## group's tolerance limit value
	cta <- aggregate(mor, list(spfac), function(x) sum(x < tlval))

	## Contingency table
	sl <- levels(spfac)
	ct <- matrix(nrow = 2, ncol = length(sl))
	colnames(ct) <- cta[, 1]
	ct[1, ] <- cta[, 2]
	n <- aggregate(mor, list(spfac), length)[, 2]
	ct[2, ] <- n - ct[1, ]
	pchsq <- chisq.test(ct, correct = F)$p.value
	if(pchsq < 0.01)
		ichisq(ct, n, mor, spfac)
	else
		list(csg = sl, cval = tlval)
}
