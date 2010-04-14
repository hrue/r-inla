\name{Drivers}
\alias{Drivers}
\docType{data}
\title{ Time series with seasonal effect}
\description{
  Montly total of car drivers killed or several
  injuried in England from January 1969 to December 1984

  NB: The last 12 lines of the data set have the first column set to
  \code{NULL} since these data where not observed but we want to predict them.
}
\usage{data(Drivers)}
\format{
  A data frame with 204 observations on the following 4 variables.
  \describe{
    \item{\code{y}}{Number of deaths}
    \item{\code{belt}}{Indicator of weather the belt was compulsory to
      use (1) or not (0)}
    \item{\code{trend}}{time (in months)}
    \item{\code{seasonal}}{time (in months)}
  }
}


\references{
 Rue, H and Held, L. (2005)
  \emph{Gaussian Markov Random Fields - Theory and Applications}
  Chapman and Hall
}
\examples{
data(Drivers)

}
\keyword{datasets}
