\name{Scotland}
\alias{Scotland}
\docType{data}
\title{Conditional Autoregressive (CAR) model for disease mapping}
\description{
  
The rate of lip cancer in  56 counties in Scotland is recorder. The data
set includes the observed and expected cases (based on the population
and its age and sex distribution in the country), a covariate measuring
the percentage of the population engaged in agricolture, fishing or
forestry and the "position" of each county expressed as a list of
adjacent counties
  
  
}
\usage{data(Scotland)}
\format{
  A data frame with 56 observations on the following 4 variables.
  \describe{
    \item{\code{Counts}}{The number of lip cancer registered}
    \item{\code{E}}{The expected number of lip cancer }
    \item{\code{X}}{The percentage of the population engaged in agricolture, fishing or
forestry }
    \item{\code{Region}}{The county}
  }
}
\source{
OpenBUGS Example manual, GeoBUGS
}
\references{

  Clayton and Kaldor (1987) and Breslow and Clayton (1993)
}
\examples{
data(Scotland)

}
\keyword{datasets}
