\name{utr.summary.plt}
\alias{utr.summary.plt}
\title{View the genomic position relative to archetype miRNA}
\usage{
utr.summary.plt(x)
}
\arguments{x: A character vector of isomiRs and isomiRs should be same as our isomiRseeker tool generated isomiRs symbols.}
\description{
Distribution of the number of isomiR endpoints at
each genomic position relative to the archetype’s coordinates.
}
\examples{
#load the dataset
data("ExampleData")
isomiRs = rownames(ExampleData)
utr.summary(x = isomiRs)
}
