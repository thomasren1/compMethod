#' Non-parametric ANOVA
#'
#' Apply a parametric F test of the hypothesis of
#' equal means between species groups to the ranked MOE data at a
#' significance level of 0.01.
#'
#' See Appendix X5, page 20 of ASTM D1990-07.
#' Follow with a Tukey multiple comparison test if needed.
#'
#' @param moe is a vector of MOE data for all species
#' @param spfac is the species factor for moe
#' @return a list giving the controlling subgroup (csg) and the mean characteristic value (cval)
#' @note the mean characteristic value is not weighted by the species
#' sample sizes
#' @export
#'
#'
#'
#'



npaov <- function(moe, spfac)
{

	## Rank data, assign average ranks to ties
	rmoe <- rank(moe)
	fitlm <- lm(rmoe ~ spfac)
	aovtab <- anova(fitlm)
	fval <- aovtab$"F value"[1]
	df1 <- aovtab$"Df"[1]    		# Treatment degrees of freedom
	df2 <- aovtab$"Df"[2]    		# Error degrees of freedom
	if(fval > qf(0.01, df1, df2, lower.tail = F))
		Tukey(moe, rmoe, spfac, df1, df2, aovtab)
	else
		list(csg = levels(spfac), cval = mean(moe))
}
