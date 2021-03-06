\name{isomiR.summary}
\alias{isomiR.summary}
\title{The number of distinct isomiRs per locus.}
\usage{
isomiR.summary(x)
}
\arguments{x: A character vector of isomiRs and isomiRs should be same as our isomiRseeker tool generated isomiRs symbols.}
\description{
Counted the number of isoform/hairpin arms that produce
from a miRNA molecules and this distinct number of hairpin arms might be exhibits signifcant
variation compaire to archetype miRNA.
}
\examples{
#load the dataset
library(MDEHT)
data("ExampleData")
isomiRs = rownames(ExampleData)
utr.summary(x = isomiRs)
}
