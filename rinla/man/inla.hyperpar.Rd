\name{inla.hyperpar}
\alias{inla.hyperpar}
\alias{hyperpar.inla}

\alias{hyperpar.inla.surv}
\alias{hyperpar.inla.inla}

\title{ Marginals for hyperparameters}
\description{
 Compute (more accurate) posterior marginals for the hyperparameters of
 the model
}
\usage{
`inla.hyperpar` =
    function(object,
             verbose = FALSE,
             dz = 0.75,
             diff.logdens = 7,
             quantiles = c(0.025, 0.5, 0.975),
             keep = FALSE
             )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{inla} or \code{surv.inla}. Usually a result of a call to
    \code{inla} or \code{surv.inla}}
  \item{verbose}{Boolean indicating wheather the inla program should run
    in a verbose mode. }
  \item{dz}{Step length used in the construction of the grid, see \code{inla}. }
  \item{diff.logdens}{Threshhold for accepting a configuration while
    building the grid, see \code{inla}.}
   \item{quantiles}{A vector of maximum 10 quantiles,
     to compute for each posterior marginal.}
   \item{keep}{A boolean variable indicating weather the working files
     should be kept.}
  }

\value{
  The object returned is of class \code{hyperpar.inla} and you may apply
  \code{plot} and \code{summary}.}
\references{
  See the references in \code{inla}}
\author{Havard Rue \email{hrue@math.ntnu.no}  }
\note{This function might take a long time if the number of hyperparameters in
  the model is large.}
\seealso{\code{\link{inla}} }
\examples{
##### Example epil #############################
### From WinBugs/OpenBugs Examples vol I #######

data(Epil)
attach(Epil)

#build centered covariates
M.Trt <- Trt-mean(Trt)
M.lBase4 <- log(Base/4)-mean(log(Base/4))
M.V4 <- V4-mean(V4)
M.lAge <- log(Age)-mean(log(Age))
M.BT <- Trt*M.lBase4-mean(Trt*M.lBase4)

formula <- y ~ M.lBase4 + M.Trt + M.BT+ M.lAge + M.V4 +
f(Ind,model="iid")+
f(rand,model="iid")

model=inla(formula,family="poisson")
hyper = hyperpar.inla(model)
summary(hyper)
plot(hyper)
}
\keyword{models, hyperparameters, inla}
