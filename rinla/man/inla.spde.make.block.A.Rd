\name{inla.spde.make.block.A}
\alias{inla.spde.make.block.A}
\title{
  Observation matrices for mesh models.
}
\description{
  Constructs observation/prediction weight matrices for numerical
  integration schemes for regional data problems.  Primarily intended
  for internal use by \code{\link{inla.spde.make.A}}.
}
\usage{
inla.spde.make.block.A(A, block, n.block = max(block),
                       weights = NULL,
                       rescale = c("none", "count", "weights", "sum"))
}
\arguments{
  \item{A}{
    A precomputed observation/prediction matrix for locations that are
    to be joined.
  }
  \item{block}{
    Indices specifying block groupings: Entries with the same
    \code{block} value are joined into a single row in the resulting
    matrix, and the \code{block} values are the row indices.
}
\item{n.block}{
  The number of blocks.
}
  \item{weights}{
    Optional scaling weights to be applied row-wise to the input
    \code{A} matrix.
  }
  \item{rescale}{
    Specifies what scaling method should be used when joining the
    rows of the \code{A} matrix as grouped by the \code{block}
    specification.
    \itemize{
    \item \code{'none'}: Straight sum, no rescaling.
    \item \code{'count'}: Divide by the number of entries in the block.
    \item \code{'weights'}: Divide by the sum of the \code{weight} values
    within each block.
    \item \code{'sum'}: Divide by the resulting row sums.
    }
  }
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
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
%%\keyword{ ~kwd1 }
