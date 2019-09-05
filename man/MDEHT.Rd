\name{MDEHT}
\alias{MDEHT}
\title{Multivariate Hotelling's T-square test in isomiRs RNA-Seq data.}
\usage{
MDEHT(X, Y, ...)
}
\arguments{
\item{X:}{ a matrix of isomiR read count data where rows represent isomiRs and columns represent samples (treat as tumor/case sample group)}

\item{Y:}{ a matrix of isomiR read count data where rows represent isomiRs and columns represent samples (treat as normal/control sample group). For Y = NULL, it will executes one sample Hotelling's test}

\item{inohf.isomiRs:}{ a logical parameter. include only high frequency isoforms in a miRNA when it contains large number of isomiRs. Default FALSE}

\item{inohf.isomiRs.no:}{ usable when inohf.isomiRs is TRUE and default value 25}

\item{merge.olf.isomRs:}{a logical parameter. merge only low-frequency isoforms into one group from a miRNA when it contains large number of isomiRs.Default FALSE}

\item{merge.olf.isomRs.no:}{usable when merge.olf.isomRs is TRUE and default value 25 }

\item{adjp.method:}{ p value adjusted method. Default to "BH".}

}

\description{
Executes multivariate Hotelling’s T-square test in isomiRs RNA-Seq data to detect differential expression miRNA. Hotelling’s T-square test is a generalization of Student's t- statistic and is widely used for testing the difference between two multivariate means. Also, Hotelling’s T-square test well controlled type I error rate and It can generate an overall conclusion even if multiple t-tests are inconsistent.
}

\value{
 \item{X.mean:}{a vector of mean values from tumor/case group}
 \item{Y.mean:}{a vector of mean values from normal/control group}
 \item{T.value:}{a vector of Hotellig's T-statistic values}
 \item{p.value:}{a vector the raw p-values from testing differential expression for each miRNA}
 \item{adj.pval:}{a vector of adjusted p-values after adjusting for multiple testing using the Benjamini-Hochberg  method}
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

results = MDEHT (X = TumorData, Y = NormalData)
head(results)

###########################
}
