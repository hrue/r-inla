\name{Epil}
\alias{Epil}
\docType{data}
\title{Repeated measures on Poisson counts}
\description{
  Seizure counts in a randomised trial of anti-convulsant therpay in
  epilepsy for 59 patients.
}
\usage{data(Epil)}
\format{
  A data frame with 236 observations on the following 7 variables.
  \describe{
    \item{\code{y}}{Number of seizures}
    \item{\code{Trt}}{indicator for the presence of treatment}
    \item{\code{Base}}{8-week baseline seizure counts}
    \item{\code{Age}}{Age of the patient}
    \item{\code{V4}}{indicator variable for the 4th visit.}
    \item{\code{rand}}{a numeric vector}
    \item{\code{Ind}}{indicator for the specific patient}
  }
}

\source{
WinBUGS/OpenBUGS Manual Examples Vol I
}

\examples{
data(Epil)

}
\keyword{datasets}
