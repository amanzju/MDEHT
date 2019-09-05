\name{MDEHT}
\alias{MDEHT}
\title{Multivariate Hotelling's T-square test in isomiRs RNA-Seq data.}
\usage{
MDEHT(X, Y, ...)
}
\arguments{
\item{X:}{ a matrix of one group (treat as tumor sample group) of RNA-Seq isomiRs counts (isomiRs are rows, samples are columns)}

\item{Y:}{ a matrix of another group (treat as normal sample group) of RNA-Seq isomiRs counts (isomiRs are rows, samples are columns) when Y =NULL then it will perform one sample Hotelling's test}

\item{inohf.isomiRs:}{ a logical parameter. if you want include only high frequency isoforms in a miRNA when it contains large number of isomiRs. Default FALSE}

\item{inohf.isomiRs.no:}{ only used when inohf.isomiRs is TRUE and default value 25}

\item{merge.olf.isomRs:}{a logical parameter. if you want merge only low frequency isoforms as a new isoform from a miRNA when it contains large number of isomiRs.Default FALSE}

\item{merge.olf.isomRs.no:}{only used when merge.olf.isomRs is TRUE and default value 25 }

\item{adjp.method:}{ p values adjusted method. Defaults to "BH".}

}

\description{
Executes multivariate Hotelling’s T-square test in isomiRs RNA-Seq data to detect differential expression miRNA. Where Hotelling’s T-square test is a generalization of Student's t- statistic and is widely used for testing the difference between two multivariate means.
}

\value{
 \item{tumor.mean:}{a vector of mean values from tumor group}
 \item{normal.mean:}{a vector of mean values from normal group}
 \item{T.value:}{a vector of Hotellig's T-statistic values}
 \item{p.value:}{a vector containing the raw p-values testing differential expression for each miRNA}
 \item{adj.pval:}{a vector containing the p-values after adjusting for multiple testing using the method of Benjamini-Hochberg}
 \item{FC:}{a vector containing the estimated log fold changes for each miRNA}
}
\examples{
################################
#required the following packages
install.packages("devtools")
library(devtools)
install_github("amanzju/MDEHT")
library(edgeR)
library(MDEHT)

#load an example dataset
data("ExampleData")

#use CPM normalization
NormalizedData = log2(cpm(ExampleData)+1)

#separate sample group
TumorData = NormalizedData[,grep("01A", colnames(NormalizedData))]
NormalData = NormalizedData[,grep("11A", colnames(NormalizedData))]

results = MDEHT(X = TumorData, Y = NormalData)
head(results)
###########################
}
