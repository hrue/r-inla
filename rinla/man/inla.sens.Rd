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
  TODO: This is an EXPERIMENTAL function!
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

  # Make some data unobserved
  nUnObs = 20
  y[(nObs-nUnObs+1):nObs] = NA
  
  # Fit the model
  mod = inla(y ~ x,
           data = list(x = x, y = y))
  
  # Calculate sensitivites
  inla.sens(mod)
  
## Case 2: Time series
  # Length of time-series
  nObs = 100
  
  # Measurement noise
  sdNoise = 0.1

  # Autoregressive process
  rho = 0.6
  sdProc = 0.1
  arP = matrix(0, nrow = nObs, ncol = 1)
  for(i in 2:nObs)
      arP[i] = rho*arP[i-1] + rnorm(1)*sdProc
  tIdx = 1:nObs

  # Coefficients
  mu = 2
  
  # Generate data
  y = mu + arP + rnorm(nObs)*sdNoise

  # Make some data unobserved
  nUnObs = 20
  y[(nObs-nUnObs+1):nObs] = NA
  idx = 1:nObs
  
  # Run INLA
  mod = inla(y ~ f(tIdx, model = "ar1"),
             data = list(y = y, tIdx = tIdx),
             control.inla = list(reordering = "metis"))
  
  # Calculate sensitivities
  inla.sens(mod)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
