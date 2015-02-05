# Grouping method for Tolerance Limit Properties as outlined in ASTM D1990-07
# (file name "D1990 Attachment 1.pdf").  See Section 10.3.1, page 6-7
# General program that handles unequal sample sizes


tl <- function(mor, a = 5, b = 75, ...)
{
    ## Calculate the order statistic for the a% non-parametric one-sided
    ## tolerance limit with b% confidence
    ##   Input:
    ##      mor is a vector of MOR data for all species
    ##	    a is the percent of the population below the tolerance limit
    ##	    b is the percent confidence of the tolerance limit
    ##	 Output:
    ##	    tolerance limit

    N <- length(mor)
    b1 <- 0:N
    ## Cumulative binomial probabilities for each value of b1 = 0,...,N
    b2 <- pbinom(b1, N, a / 100)
    ## The index r of the order statistic so that at least b% of the 
    ## time the tolerance limit will lie below the value of interest.
    r <- max(b1[b2 < 1 - b / 100]) + 1
    return(sort(mor)[r])
}


chisq <- function(mor, spfac)
{
    ## Chi Square test - determine the number of pieces in each species group
    ## with MOR values below the combined group's tolerance limit value and
    ## conduct a Chi Square test at a significance level of 0.01.  See
    ## Appendix X7, page 21-22 of ASTM D1990-07.
    ## Follow with iterative Chi Square tests on species subsets if needed.
    ##	 Input:
    ##	    mor is a vector of MOR data for all species
    ##	    spfac is the species factor for mor
    ##   Output:
    ##	    a list giving the controlling subgroup (csg) and the tolerance
    ##	    limit characteristic value (cval)

    ## NOTE : no continuity correction is applied when computing the test
    ## statistic which results in a warning message

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

