\name{inla.spde.models}
\alias{inla.spde.models}
\alias{inla.spde1.models}
\alias{inla.spde2.models}
\title{
  List SPDE models supported by inla.spde objects
%%  ~~function to do ... ~~
}
\description{
  List SPDE models supported by inla.spde objects
}
\usage{
inla.spde.models(function.names=FALSE)
inla.spde1.models()
inla.spde2.models()
}
\arguments{
  \item{function.names}{
    If \code{FALSE}, return list model name lists.
    If \code{TRUE}, return list of model object constructor function
    names.
  }
}
\details{
  Returns a list of available SPDE model type name lists, one for each
  inla.spde model class (currently \code{\link{inla.spde1}} and
  \code{\link{inla.spde2}}).
}
\value{
  List of available SPDE model type name lists.
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

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
## Display help for each supported inla.spde2 model:
for (model in inla.spde2.models())
  print(help(paste("inla.spde2.", model, sep="")))

## Display help for each supported inla.spde* model:
models = inla.spde.models()
for (type in names(models))
  for (model in models[[type]])
    print(help(paste("inla.", type, ".", model, sep="")))

## Display help for each supported inla.spde* model (equivalent to above):
for (model in inla.spde.models(function.names=TRUE))
  print(help(model))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
