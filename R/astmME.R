# Grouping method for Mean Properties as outlined in ASTM D1990-07 (file name
# "D1990 Attachment 1.pdf").  See Section 10.2.1, page 6
# General program that handles unequal sample sizes


npaov <- function(moe, spfac)
{
    ## Non-parametric ANOVA - apply a parametric F test of the hypothesis of
    ## equal means between species groups to the ranked MOE data at a
    ## significance level of 0.01.  See Appendix X5, page 20 of ASTM D1990-07.
    ## Follow with a Tukey multiple comparison test if needed.
    ##   Input:
    ##	    moe is a vector of MOE data for all species
    ##	    spfac is the species factor for moe 
    ##   Output:
    ##	    a list giving the controlling subgroup (csg) and the mean
    ##	    characteristic value (cval)

    ## NOTE : the mean characteristic value is not weighted by the species
    ## sample sizes

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


Tukey <- function(moe, rmoe, spfac, df1, df2, aovtab)
{
    ## Tukey's multiple comparison test - apply to the ranked data to 
    ## determine the group of species that are indistinguishable from the
    ## species with the lowest median (the one with the lowest mean of the
    ## ranked data) at a significance level of 0.01.  See Appendix X6, page 20
    ## of ASTM D1990-07.
    ##   Input:
    ##      moe is a vector of sample MOE data for all species
    ##      rmoe is the ranked sample MOE data for all species
    ##      spfac is the species factor for moe and rmoe 
    ##	    df1 is the treatment degrees of freedom for the non-parametric
    ##	    ANOVA test
    ##	    df2 is the error degrees of freedom for the non-parametric ANOVA
    ##	    test
    ##	    aovtab is the ANOVA table for the non-parametric ANOVA test
    ##   Output:
    ##      a list giving the controlling subgroup (csg) and the mean
    ##      characteristic value (cval)

    ## NOTE - the mean characteristic value is not weighted by sample size	    
    ## Mean differences between each species and the one with the lowest mean
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

