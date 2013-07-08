\name{Zambia}
\alias{Zambia}
\docType{data}
\title{ Semiparametric regression}
\description{
Undernutrition of children in each region of Zambia is measured through
a score computed on the basis of some anthropometric measures. The data
set contains also other infomation about each child.
}
\usage{data(Zambia)}
\format{
  A data frame with 4847 observations on the following 10 variables.
  \describe{
    \item{\code{hazstd}}{standardised Z score of stunting}
    \item{\code{bmi}}{body mass index of the mother}
    \item{\code{agc}}{age of the child in months}
    \item{\code{district}}{district where the child lives}
    \item{\code{rcw}}{mother employment status with categories "working"
    (1) and "not working" (-1)}
    \item{\code{edu1}}{mother's educations status with categories
      "complete primary but incomplete secondary " (\code{edu1=1)},
      "complete secondary or higher" (\code{edu2=1}) and "no education
      or incomplete primary" (\code{edu1=edu2=-1})}
    \item{\code{edu2}}{see above}
    \item{\code{tpr}}{locality of the domicile with categories "urban"
      (1) and "rural" (-1)}
    \item{\code{sex}}{gender of the child with categories "male" (1) and
    "female" (-1)}
    \item{\code{edu}}{DO NOT KNOW; check source}
}}

\source{
BayesX Manual \url{http://www.stat.uni-muenchen.de/~bayesx/bayesx.html}
}

\examples{
data(Zambia)

}
\keyword{datasets}
