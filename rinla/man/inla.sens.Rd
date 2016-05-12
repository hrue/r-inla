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

## Case 3: Epil dataset
  data(Epil)
  my.center = function(x) (x - mean(x))

  # make centered covariates
  Epil$CTrt    = my.center(Epil$Trt)
  Epil$ClBase4 = my.center(log(Epil$Base/4))
  Epil$CV4     = my.center(Epil$V4)
  Epil$ClAge   = my.center(log(Epil$Age))
  Epil$CBT     = my.center(Epil$Trt*Epil$ClBase4)

  # Define the model
  formula = y ~ ClBase4 + CTrt + CBT+ ClAge + CV4 +
            f(Ind, model="iid") + f(rand,model="iid")

  mod = inla(formula,family="poisson", data = Epil)
 
  # Calculate sensitivities
  inla.sens(mod)

## Case 4: Spatial data
  # Number of observations
  nObs = 100
  
  # Measurement noise
  sdNoise = 0.2

  # Spatial process
  sdProc = 1.0
  rho0 = 0.2

  # Coefficients
  beta0 = 1
  beta1 = 2
  
  # Generate spatial data + measurement noise
  loc = cbind(runif(nObs), runif(nObs))
  dd = as.matrix(dist(loc))
  Sig = sdProc^2*inla.matern.cov(nu = 1, kappa = sqrt(8)/rho0, dd, corr = TRUE)
  L = t(chol(Sig))
  u = L%*%rnorm(nObs) + rnorm(nObs)*sdNoise
  
  # Generate Covariate
  x = runif(nObs)-0.5

  # Combine to observations
  y = beta0 + beta1*x + u
  
  # Number of unobserved
  nUnObs = 2
  y[1:nUnObs] = NA
  
  # Mesh
  mesh = inla.mesh.2d(loc, max.edge = 0.05, cutoff = 0.05)
  
  # Make SPDE object
  spde  = inla.spde2.matern(mesh)
  spde2 = inla.spde2.matern(mesh, constr = TRUE)
  
  # Make A matrix
  A = inla.spde.make.A(mesh, loc)
  
  # Stack
  X = cbind(1, x)
  stk = inla.stack(data = list(y = y), A = list(A, 1),
                   effects = list(field = 1:spde$n.spde,
                                  X = X))
  
  # Run INLA
  mod1 = inla(y ~ -1 + X + f(field, model = spde),
              data = inla.stack.data(stk),
              control.predictor = list(A = inla.stack.A(stk)),
              control.family = list(prior = "pcprec",
                                    param = c(3, 0.05)))
  mod2 = inla(y ~ -1 + X + f(field, model = spde2),
              data = inla.stack.data(stk),
              control.predictor = list(A = inla.stack.A(stk)),
              control.family = list(prior = "pcprec",
                                    param = c(3, 0.05)))

  # Calculate sensitivities
  res1 = inla.sens(mod1)
  res2 = inla.sens(mod2)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }
