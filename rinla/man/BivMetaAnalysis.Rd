\name{BivMetaAnalysis}
\alias{BivMetaAnalysis}
\docType{data}
\title{Bivariate Meta Analysis}
\description{
  Data are taken from a meta-analysis 
  to compare the utility of three types of diagnostic imaging -
  lymphangiography (LAG), computed tomography (CT) and magnetic
  resonance (MR) - to detect lymph node metastases in patients with
  cervical cancer. The dataset consists of a total of 46 studies:
  the first 17 for LAG, the following 19 for CT and the last
  10 for MR. 
}
\usage{data(BivMetaAnalysis)}
\format{
  A data frame with 92 observations on the following 9 variables.
  \describe{
    \item{\code{N}}{a numeric vector}
    \item{\code{Y}}{a numeric vector}
    \item{\code{diid}}{a numeric vector}
    \item{\code{lag.tp}}{a numeric vector}
    \item{\code{lag.tn}}{a numeric vector}
    \item{\code{ct.tp}}{a numeric vector}
    \item{\code{ct.tn}}{a numeric vector}
    \item{\code{mr.tp}}{a numeric vector}
    \item{\code{mr.tn}}{a numeric vector}
  }
}
%%\details{}
%%\source{}
\references{
   J. Scheidler and H. Hricak and K. K. Yu and L. Subak and M. R. 
Segal,"Radiological evaluation of lymph node metastases in patients with 
cervical cancer: a meta-analysis",JAMA 1997
 }
\examples{
data(BivMetaAnalysis)
}
\keyword{datasets}
