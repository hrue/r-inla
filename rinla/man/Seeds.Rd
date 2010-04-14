\name{Seeds}
\alias{Seeds}
\docType{data}
\title{ Factorial design}
\description{
  
  Proportion of seeds that germinated on each of 21 plates arranged according to a 2 by 2 factorial layout by seed and type of root extract
}
\usage{data(Seeds)}
\format{
  A data frame with 21 observations on the following 5 variables.
  \describe{
    \item{\code{r}}{number of germinated seeds per plate}
    \item{\code{n}}{number of total seeds per plate}
    \item{\code{x1}}{seed type}
    \item{\code{x2}}{root extracted}
    \item{\code{plate}}{indicator for the plate}
  }
}
\source{
 WinBUGS/OpenBUGS Manual Example, Vol. I 
}
\examples{
data(Seeds)

}
\keyword{datasets}
