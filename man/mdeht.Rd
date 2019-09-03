\name{MDEHT}
\alias{MDEHT}
\title{Multivariate Hotelling's T-square Test.}
\usage{
mdeht(X, Y, ...)
}
\description{
Executes multivariate Hotelling’s T-square test in isomiRs RNA-Seq data to detect differential expression miRNA. Where Hotelling’s T-square test is a generalization of Student's t- statistic and is widely used for testing the difference between two multivariate means.
}
\examples{
mdeht(X = TumorData, Y = NormalData)
}
