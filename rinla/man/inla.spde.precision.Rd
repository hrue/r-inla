\name{inla.spde.precision}
\alias{inla.spde.precision}
\alias{inla.spde1.precision}
\alias{inla.spde2.precision}
\alias{inla.spde.precision.inla.spde1}
\alias{inla.spde.precision.inla.spde2}
\title{Add title here
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
inla.spde2.precision(spde,
                     theta = NULL,
                     phi0 = inla.spde2.theta2phi0(spde, theta),
                     phi1 = inla.spde2.theta2phi1(spde, theta),
                     phi2 = inla.spde2.theta2phi2(spde, theta))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{spde}{
%%     ~~Describe \code{spde} here~~
}
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
  \item{phi0}{
%%     ~~Describe \code{phi0} here~~
}
  \item{phi1}{
%%     ~~Describe \code{phi1} here~~
}
  \item{phi2}{
%%     ~~Describe \code{phi2} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (spde, theta = NULL, phi0 = inla.spde2.theta2phi0(spde, 
    theta), phi1 = inla.spde2.theta2phi1(spde, theta), phi2 = inla.spde2.theta2phi2(spde, 
    theta)) 
{
    inla.require.inherits(spde, "inla.spde2", "'spde'")
    if (spde$f$model != "spde2") {
        stop("'inla.spde2.precision' only supports internal inla models 'spde2'")
    }
    D0 = Diagonal(spde$n.spde, phi0)
    D1 = Diagonal(spde$n.spde, phi1)
    D12 = Diagonal(spde$n.spde, phi1 * phi2)
    Q = (D0 \%*\% (D1 \%*\% spde$param.inla$M0 \%*\% D1 + D12 \%*\% spde$param.inla$M1 + 
        t(spde$param.inla$M1) \%*\% D12 + spde$param.inla$M2) \%*\% 
        D0)
    return(Q)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
