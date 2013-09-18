\name{inla.generate.colors}
\alias{inla.generate.colors}
\title{
  Generate text RGB color specifications.
}
\description{
  Generates a tex RGB color specification matrix based on a color palette.
}
\usage{
inla.generate.colors(color,
                     color.axis = NULL,
                     color.n = 512,
                     color.palette = cm.colors,
                     color.truncate = FALSE,
                     alpha = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{color}{
    \code{character}, \code{matrix} or \code{vector}
}
\item{color.axis}{
  The min/max limit values for the color mapping.
}
\item{color.n}{
  The number of colors to use in the color palette.
}
\item{color.palette}{
  A color palette function.
}
\item{color.truncate}{
  If \code{TRUE}, truncate the colors at the color axis limits.
}
\item{alpha}{
  Transparency/opaqueness values.
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
%%\keyword{ ~kwd1 }
