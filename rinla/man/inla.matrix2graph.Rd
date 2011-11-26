\name{inla.matrix2graph}
\alias{inla.matrix2graph}
\alias{inla.graph2matrix}
\title{
  Convert between an adjecency matrix and the INLA graph description
  file format.
}
\description{
  Converts an adjecency matrix to the INLA graph description
  file format, and from a graph file to an adjecency matrix.
}
\usage{
inla.matrix2graph(Q, graph.file = "graph.txt", c.indexing = FALSE)
inla.graph2matrix(graph.file = "graph.txt")
}
\arguments{
  \item{Q}{
    Graph adjacency matrix
  }
  \item{graph.file}{
    Filename for the INLA graph file.
  }
  \item{c.indexing}{
    If \code{TRUE}, use zero-based indexing, if \code{FALSE} (default)
    use one-based indexing.
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  For \code{inla.matrix2graph}, the graph file is written to disk.

  For \code{inla.graph2matrix}, a sparse binary adjacency matrix is returned.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{inla.read.graph}},
  \code{\link{inla.write.graph}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
G = Matrix(c(1,1,0,0, 1,1,1,0, 0,1,1,1, 0,0,1,1),4,4)
inla.matrix2graph(G, "graph.txt")
graph = inla.read.graph("graph.txt")
summary(graph)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
