\name{inla.dBind}
\alias{inla.dBind}
\alias{inla.dBind-deprecated}
\title{
  Build a block-diagonal sparse matrix.
}
\description{
  Build a block-diagonal sparse matrix.
  Obsolete wrapper for \code{bdiag()}.
}
\usage{
inla.dBind(...)
}
\arguments{
  \item{\dots}{
    A list of square or rectangular matrices.
}
}
\value{
  A sparse matrix.
}
\author{
Finn Lindgren \email{finn.lindgren@gmail.com}
}
\examples{
inla.dBind(Matrix(1,2,1), Matrix(2,1,2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
