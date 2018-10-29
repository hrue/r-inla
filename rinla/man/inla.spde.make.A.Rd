\name{inla.spde.make.A}
\alias{inla.spde.make.A}
\title{
  Observation/prediction matrices for mesh models.
}
\description{
  Constructs observation/prediction weight matrices for models based on
  \code{\link{inla.mesh}} and \code{\link{inla.mesh.1d}} objects. 
}
\usage{
inla.spde.make.A(mesh = NULL, loc = NULL, index = NULL,
                 group = NULL, repl = 1L,
                 n.spde = NULL, n.group = NULL, n.repl = NULL,
                 group.mesh = NULL,
                 weights = NULL,
                 A.loc = NULL, A.group = NULL, group.index = NULL,
                 block = NULL, n.block = NULL,
                 block.rescale = c("none", "count", "weights", "sum"),
                 ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mesh}{
    An \code{\link{inla.mesh}} or \code{\link{inla.mesh.1d}} object
    specifying a function basis on a mesh domain.
    Alternatively, an
    \code{inla.spde} object that includes a mesh (e.g. from
    \code{\link{inla.spde2.matern}}).}
\item{loc}{
  Observation/prediction coordinates.  \code{mesh} and \code{loc}
  defines a matrix \code{A.loc} of mapping weights between basis
  function weights and field values.  If \code{loc} is \code{NULL},
  \code{A.loc} is defined as \code{Diagonal(n.spde, 1)}.
}
\item{index}{
  For each observation/prediction value, an index into \code{loc}.
  Default is \code{seq_len(nrow(A.loc))}.
}
\item{group}{
  For each observation/prediction value, an index into the group model.
}
  \item{repl}{
  For each observation/prediction value, the replicate index.
}
\item{n.spde}{
  The number of basis functions in the mesh model. (Note: may be
  different than the number of mesh vertices/nodes/knots.)
}
\item{n.group}{
  The size of the group model.
}
\item{n.repl}{
  The total number of replicates.
}
\item{group.mesh}{
  An optional \code{\link{inla.mesh.1d}} object for the group model.
}
\item{weights}{
  Optional scaling weights to be applied row-wise to the resulting matrix.
}
\item{A.loc}{
  Optional precomputed observation/prediction matrix. \code{A.loc} can
  be specified instead of \code{mesh}+\code{loc}, optionally with
  \code{index} supplied.
}
\item{A.group}{
  Optional precomputed observation/prediction matrix for the group
  model. \code{A.group} can be specified instead of \code{group} and/or
  \code{group.mesh}, optionally with \code{group.index} supplied.
}
  \item{group.index}{
  For each observation/prediction value, an index into the rows of
  \code{A.group}.
}
\item{block}{
  Optional indices specifying block groupings: Entries with the same
  \code{block} value are joined into a single row in the resulting
  matrix, and the \code{block} values are the row indices.  This is
  intended for construction of approximate integration schemes for
  regional data problems. See \code{\link{inla.spde.make.block.A}} for
  details.
}
\item{n.block}{
  The number of blocks.
}
  \item{block.rescale}{
  Specifies what scaling method should be used when joining entries as
  grouped by a \code{block} specification.  See
  \code{\link{inla.spde.make.block.A}} for details.
}
  \item{\dots}{
  Additional parameters.  Currently unused.
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
  \code{\link{inla.spde.make.index}}
}
\examples{
loc = matrix(runif(10000*2)*1000,10000,2)
mesh = inla.mesh.2d(loc=loc,
                    cutoff=50,
                    max.edge=c(50,500))
A = inla.spde.make.A(mesh, loc=loc)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
