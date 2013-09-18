\name{inla.row.kron}
\alias{inla.row.kron}
\title{
  Row-wise Kronecker products
}
\description{
  Takes two Matrices and computes the row-wise Kronecker product.
  Optionally applies row-wise weights and/or applies an additional
  0/1 row-wise Kronecker matrix product, as needed by
  \code{\link{inla.spde.make.A}}.
}
\usage{
inla.row.kron(M1, M2, repl = NULL, n.repl = NULL, weights = NULL)
}
\arguments{
  \item{M1}{
    A matrix that can be transformed into a sparse Matrix.
}
  \item{M2}{
    A matrix that can be transformed into a sparse Matrix.
}
\item{repl}{
  An optional index vector.  For each entry, specifies which replicate
  the row belongs to, in the sense used in \code{\link{inla.spde.make.A}}.
}
\item{n.repl}{
  The maximum replicate index, in the sense used in
  \code{\link{inla.spde.make.A}}.
}
\item{weights}{
  Optional scaling weights to be applied row-wise to the resulting matrix.
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
  A \code{sparseMatrix} object.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Finn Lindgren \email{finn.lindgren@gmail.com}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{inla.spde.make.A}}
}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
