\name{Kidney}
\alias{Kidney}
\docType{data}
\title{ Kidney infection data}
\description{
Times of infection from the time to insertion of the catheter for 38 kindey patients using portable dialysis equipment
}
\usage{data(Kidney)}
\format{
  A data frame with 76 observations on the following 9 variables.
  \describe{
    \item{\code{time}}{a numeric vector. Time to infection from the insertion of catheter}
    \item{\code{event}}{a numeric vector. 1: time of infection 0: time of censuring }
    \item{\code{age}}{a numeric vector. Age of the patient at the time of infection}
    \item{\code{sex}}{a numeric vector. Sex of the patient 0: male 1:female}
    \item{\code{disease}}{a numeric vector. Type of disease}
    \item{\code{dis1}}{a numeric vector. Dummy variable to codify the disease type.}
    \item{\code{dis2}}{a numeric vector. Dummy variable to codify the disease type.}
    \item{\code{dis3}}{a numeric vector. Dummy variable to codify the disease type.}
    \item{\code{ID}}{a numeric vector. Patient code.}
  }
}

\references{
  McGilchrist and C.W. Aisbett (1991),
   Regression with frailty in survival analysis,
  Biometrics,vol.47,pages 461--166.

  D.J. Spiegelhalter and A. Thomas and N.G. Best and W.R. Gilks (1995)
  BUGS: Bayesian Inference Using Gibbs sampling, Version 0.50.,
  MRC Biostatistics Unit, Cambridre, England.
}


\keyword{datasets}
