\name{inla.collect.results}
\alias{inla.collect.results}

\title{ Collect results from a inla run }
\description{\code{inla.collect.results} collect results  from a inla run.}

\usage{
inla.collect.results(results.dir,
                     control.results = inla.set.control.results.default(),
                     debug = FALSE, only.hyperparam=FALSE)
}

\arguments{
  \item{results.dir}{The directory where the results of the inla run are
  stored}
  \item{control.results}{a list of parameters controlling the
    output of the function. The list includes the following elements:
    
    return.marginals.random: a logical value indicating wheather
      the densities for the posterior marginals of the random elements
      should be retuned (default \code{FALSE}) 
    
      return.marginals.predictor: logical value indicating wheather posterior
      densities for linear predictor and fitted values should be returned
      (default \code{FALSE}).
  
  }
  \item{debug}{Binary variable, if TRUE some debugging information are printed  }

  \item{only.hyperparam}{Binary variable indicating wheather only the
    results for the hyperparameters should be collected}
}

\value{ The function returns an object of class \code{"inla"}, see the
help file for \code{inla} for details.}

\details{
  This function is mainly used inside \code{inla} and \code{inla.surv} to collect
  results after running the inla function. It can also be used to collect
  results into R after having runned a inla section outside R.
  
  In this last case some criteria have to be followed when writing the
  .ini file:

  1. All sections in the .ini file must have the same value in the field
  "compute". If "compute=0" then the argument \code{only.hyperparam} has
  to be set to \code{TRUE}

  2. The field "dir" in the sections of "type=ffield" has to start by
  "random.effect" and in the section of "type=linear" as to start by
  "fixed.effect"
      
  3. The inla program has to be run using the option "-b" which makes it
  print the output in binary format.
}


\author{Sara Martino \email{martino@math.ntnu.no}}

\keyword{IO}
