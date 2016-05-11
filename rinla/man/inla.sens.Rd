\name{inla.sens}
\alias{inla.sens}
\title{
Testing code for sensitivity
}
\description{
  TODO: Write a description
}
\usage{
inla.sens(inlaObj)
}
\arguments{
  \item{inlaObj}{
    The result from a run of \code{inla}.
  }
}
%\details{
  %%  ~~ If necessary, more details than the description above ~~
%}
\value{
  TODO: Describe what is returned
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Geir-Arne Fuglstad \email{geirarne.fuglstad@gmail.com}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{inla}}
}
\examples{
## Case 1: Simple linear regression on simulated data
  # Number of observations 
  nObs = 100

  # Measurement noise
  sdNoise = 0.1

  # Coefficients
  mu = 2
  beta = 1

  # Covariate
  x = runif(nObs)

  # Generate data
  y = mu + beta*x + rnorm(nObs)*sdNoise

  # Fit the model
  mod = inla(y ~ x,
             data = list(x = x, y = y))

  # Calculate sensitivites
  inla.sens(mod)

  # Close plot windows
  dev.off()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
