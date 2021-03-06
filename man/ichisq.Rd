% Generated by roxygen2 (4.0.2): do not edit by hand
\name{ichisq}
\alias{ichisq}
\title{Iterative Chi Square tests}
\usage{
ichisq(ct, n, mor, spfac)
}
\arguments{
\item{ct}{is the contingency table for all species}

\item{n}{is a vector of sample sizes for all species}

\item{mor}{is a vector of sample MOR data for all species}

\item{spfac}{is the species factor for mor}
}
\value{
a list giving the controlling subgroup (csg) and the tolerance
limit characteristic value (cval)
}
\description{
determine the group of species that are
indistinguishable from the species with the highest percent of pieces
below the combined group's tolerance limit value.
}
\details{
Apply at a significance level of 0.01.
}
\note{
no continuity correction is applied when computing the test
statistics which result in warning messages
}

