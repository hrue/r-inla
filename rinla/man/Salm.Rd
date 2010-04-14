\name{Salm}
\alias{Salm}
\docType{data}
\title{ Extra-Poisson variation in dose-response study}
\description{
Breslow (1984) analyses some mutagenicity assay data (shown below) on salmonella in which three plates have been processed at each dose i of quinoline and the number of revertant colonies of TA98 Salmonella measured
}
\usage{data(Salm)}
\format{
  A data frame with 18 observations on the following 3 variables.
  \describe{
    \item{\code{y}}{number of salmonella bacteria}
    \item{\code{dose}}{dose of quinoline (mg per plate)}
    \item{\code{rand}}{indicator}
  }
}

\source{
WinBUGS/OpenBUGS manual Examples VOl.I
}

\examples{
data(Salm)

}
\keyword{datasets}
