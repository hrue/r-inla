\name{PRprec}
\alias{PRprec}
\docType{data}
\title{Parana' precipitation data}
\description{
  Daily accumulated precipitation on 616 stations on 
  each day of year 2012.
}
\usage{data(PRprec)}
\format{
  two objects:
  \code{PRprec}: a \code{data.frame} with the 616 rows and 
  368 columns. Each row represents each location. 
  The colums are
  \describe{
    \item{\code{Longitude}}{Longitude coordinate}
    \item{\code{Latitude}}{Latitude coordinate}
    \item{\code{Altitude}}{Altitude at the location}
    \item{\code{...}}{365 additinal columns, one each containing 
      the precipitation on each day from 2012-01-01 to 2012-03-31}
  }
  \code{PRborder}: a \code{matrix} with 2055 with a 
  Parana' state border coordinates 
  \describe{
    \item{\code{Longitude}}{Longitude coordinate}
    \item{\code{Latitude}}{Latitude coordinate}
  }
}

\source{
ANA - Agencia Nacional de Aguas 
\url{http://www.ana.gov.br}
}

\examples{
data(PRprec)

}
\keyword{datasets}
