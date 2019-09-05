\name{utr.summary}
\alias{utr.summary}
\title{View the genomic position relative to archetype miRNA}
\usage{
utr.summary(x)
}
\arguments{x: A character vector of isomiRs and isomiRs should be same as our isomiRseeker tool generated isomiRs symbols.}
\description{
Distribution of the number of isomiR endpoints at
each genomic position relative to the archetype’s coordinates.
}
\examples{
#load the dataset
library(MDEHT)
data("ExampleData")
isomiRs = rownames(ExampleData)
utr.summary(x = isomiRs)
}
