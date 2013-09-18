\name{inla.extract.el}
\alias{inla.extract.el}
\alias{inla.extract.el.data.frame}
\alias{inla.extract.el.list}
\alias{inla.extract.el.matrix}
\title{
  Extract elements by matching name from container objects.
}
\description{
  Extract elements by wildcard name matching from a \code{data.frame},
  \code{list}, or \code{matrix}.
}
\usage{
inla.extract.el(M, ...)

\method{inla.extract.el}{data.frame}(M, match, by.row = TRUE, ...) 
\method{inla.extract.el}{list}(M, match, ...) 
\method{inla.extract.el}{matrix}(M, match, by.row = TRUE, ...) 
}
\arguments{
  \item{M}{
    A container object.
  }
  \item{match}{
    A regex defining the matching criterion.
  }
  \item{by.row}{
    If \code{TRUE}, extract data by row, otherwise by column.
  }
  \item{\dots}{
    Additional arguments, not used.
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ something }
