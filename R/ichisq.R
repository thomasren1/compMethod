#' Iterative Chi Square tests
#'
#' determine the group of species that are
#' indistinguishable from the species with the highest percent of pieces
#' below the combined group's tolerance limit value.
#'
#' Apply at a significance level of 0.01.
#'
#' @param ct is the contingency table for all species
#' @param n is a vector of sample sizes for all species
#' @param mor is a vector of sample MOR data for all species
#' @param spfac is the species factor for mor
#'
#' @return a list giving the controlling subgroup (csg) and the tolerance
#' limit characteristic value (cval)
#'
#' @note no continuity correction is applied when computing the test
#' statistics which result in warning messages
#'
#' @export

ichisq <- function(ct, n, mor, spfac)
{
	## Iterative Chi Square tests - determine the group of species that are
	## indistinguishable from the species with the highest percent of pieces
	## below the combined group's tolerance limit value.  Apply at a
	## significance level of 0.01.
	##	 Input:
	##	    ct is the contingency table for all species
	##	    n is a vector of sample sizes for all species
	##	    mor is a vector of sample MOR data for all species
	##	    spfac is the species factor for mor
	##   Output:
	##      a list giving the controlling subgroup (csg) and the tolerance
	##	    limit characteristic value (cval)

	## NOTE : no continuity correction is applied when computing the test
	## statistics which result in warning messages

	## Sort the contingency table in decreasing order by percent of pieces
	## below the combined group tolerance limit
	## NOTE - tied species are not handled in any specific way
	sct <- ct[, order(ct[1, ] / n, decreasing = T)]
	l <- 2
	repeat
	{
		m <- sct[, 1:l]
		ifelse(chisq.test(m, correct = F)$p.value < 0.01, break, l <- l + 1)
	}
	csg <- colnames(m)[-l]
	smor <- mor[spfac%in%csg]
	cval <- tl(smor)
	list(csg = csg, cval = cval)
}
