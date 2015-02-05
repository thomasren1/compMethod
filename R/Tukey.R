#'
#' Tukey's multiple comparison test
#'
#' apply to the ranked data to
#' determine the group of species that are indistinguishable from the
#' species with the lowest median (the one with the lowest mean of the
#' ranked data) at a significance level of 0.01.
#'
#' See Appendix X6, page 20 of ASTM D1990-07.
#'
#'
#' @param moe is a vector of sample MOE data for all species
#' @param rmoe is the ranked sample MOE data for all species
#' @param spfac is the species factor for moe and rmoe
#' @param df1 is the treatment degrees of freedom for the non-parametric ANOVA test
#' @param df2 is the error degrees of freedom for the non-parametric ANOVA test
#' @param aovtab is the ANOVA table for the non-parametric ANOVA test
#'
#' @return a list giving the controlling subgroup (csg) and the mean
#' characteristic value (cval)
#'
#' @note  the mean characteristic value is not weighted by sample size
#' Mean differences between each species and the one with the lowest mean
#'
#' @export
#'
#'
#'



Tukey <- function(moe, rmoe, spfac, df1, df2, aovtab)
{

	gm <- aggregate(rmoe, list(spfac), mean)
	gmv <- gm[, 2]			# group mean values
	mmv <- min(gmv)			# minimum mean value
	gmd <- gmv - mmv			# group mean differences

	## Calculate w
	gn <- aggregate(rmoe, list(spfac), length)
	gnv <- gn[, 2]			# group sample sizes
	ind <- gmv == mmv
	nmmv <- gn[ind, 2]			# sample size of minimum mean value
	n <- 2 / ((1 / gnv) + (1 / nmmv))
	w <- qtukey(0.01, df1+1, df2, lower.tail = F) * sqrt(aovtab$"Mean Sq"[2]) * sqrt(1 / n)

	## Any two means more than w apart are significantly different and are not
	## in the controlling subgroup
	csg <- as.character(gm[!gmd > w, 1])
	cval <- mean(moe[spfac%in%csg])

	list(csg = csg, cval = cval)
}
