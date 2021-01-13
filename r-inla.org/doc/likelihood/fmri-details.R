### some details for family='fmri', by Ludger Benedikt Starke,
### Max Delbrück Center for Molecular Medicine, Berlin, Germany.

rm(list = ls())

## PDFs of the non-central chi via 'dchisq' ------------------------------------

# original parameters central
dchi <- function(z, nu) {
  out <- 2*z*dchisq(z^2, df = nu, ncp = 0, log = FALSE)
}


# non-central
dchi.nc <- function(z, nu, rho) {
  out <- 2*z*dchisq(z^2, df = nu, ncp = rho^2, log = FALSE)
}


# plus precision parameter
# z = sqrt(tau)*y
dchi.tau <- function(y, nu, lambda, tau) {
  out <- 2*tau*y*dchisq(tau*y^2, df = nu, ncp = tau*lambda^2, log = FALSE)
}


# plus exponential link
dchi.eta <- function(y, nu, eta, tau) {
  out <- 2*tau*y*dchisq(tau*y^2, df = nu, ncp = tau*exp(2*eta), log = FALSE)
}


# plus truncation
dchi.trunc <- function(y, nu, eta, tau, a) {
  pBelowA <- pchisq(tau*a^2, nu, ncp = tau*exp(2*eta))
  out <- 2*tau*y*dchisq(tau*y^2, df = nu, ncp = tau*exp(2*eta), log = FALSE)/(1 - pBelowA)
}



## tests: area = 1? ------------------------------------------------------------

# input parameters
nu <- 2.4                     # degree of freedom
rho <- 2.2                    # non-centrality parameter

tau = 0.5                     # precision parameter
lambda <- 1/sqrt(tau)*rho     # scaled non-centrality

eta = log(lambda)             # linear predictor
a = 2.9                       # truncation threshold

# original parameters non-central
integrand <- function(x) {dchi.nc(z = x, nu = nu, rho = rho)}
int <- integrate(integrand, lower = 10^(-10), upper = Inf)
area.dchi.ndc <- int$value

# plus precision parameter
integrand <- function(x) {dchi.tau(y = x, nu = nu, lambda = lambda, tau = tau)}
int <- integrate(integrand, lower = 10^(-10), upper = Inf)
area.dchi.tau <- int$value

# plus exponential link
integrand <- function(x) {dchi.eta(y = x, nu = nu, eta = eta, tau = tau)}
int <- integrate(integrand, lower = 10^(-10), upper = Inf)
area.dchi.eta <- int$value

# plus truncations
integrand <- function(x) {dchi.trunc(y = x, nu = nu, eta = eta, tau = tau, a = a)}
int <- integrate(integrand, a, upper = Inf)
area.dchi.trunc <- int$value

print(c(area.dchi.ndc, area.dchi.tau, area.dchi.eta, area.dchi.trunc))

 

## explicit PDFs of the non-central chi for comparison: ------------------------
# original parameters
ncChi.original <- function(z, nu, rho) {

  out <- rho*(exp(-(z^2 + rho^2)/2))*(z/rho)^(nu/2)*besselI(rho*z, nu/2 - 1)
  return(out)
}


# with precision parameter
ncChi.tau <- function(y, nu, lambda, tau) {

  out <- tau*lambda*exp(-tau*(y^2 + lambda^2)/2)*(y/lambda)^(nu/2)*besselI(tau*lambda*y, nu/2 - 1)
  return(out)
}


# exponential link
ncChi.eta <- function(y, nu, eta, tau) {

  out <- tau*exp(eta - tau*(y^2 + lambda^2)/2)*(y/exp(eta))^(nu/2)*besselI(tau*exp(eta)*y, nu/2 - 1)
  return(out)
}

# exponential link with truncation
ncChi.trunc <- function(y, nu, eta, tau, a) {

  pSmallerA <- pchisq((sqrt(tau)*a)^2, nu, ncp = (sqrt(tau)*exp(eta))^2)
  out <- tau*exp(eta - tau*(y^2 + lambda^2)/2)*(y/exp(eta))^(nu/2)*besselI(tau*exp(eta)*y, nu/2 - 1)/(1 - pSmallerA)
  return(out)
}



## comparison to dchisq versions -----------------------------------------------

test.vector <- seq(10^(-10), 4*lambda, length.out = 1000)

# original parameters
pdf.dchi.nc <- dchi.nc(test.vector, nu = nu, rho = rho)
pdf.ncChi.original <- ncChi.original(test.vector, nu = nu, rho = rho)

print(c('max difference dchi', max(abs(pdf.dchi.nc - pdf.ncChi.original))))

# with precision parameter
pdf.dchi.tau <- dchi.tau(test.vector, nu = nu, lambda = lambda, tau = tau)
pdf.ncChi.tau <- ncChi.tau(test.vector, nu = nu, lambda = lambda, tau = tau)

print(c('max difference, with scaling', max(abs(pdf.dchi.tau - pdf.ncChi.tau))))

# exponential link
pdf.dchi.eta <- dchi.eta(test.vector, nu = nu, eta = eta, tau = tau)
pdf.ncChi.eta <- ncChi.eta(test.vector, nu = nu, eta = eta, tau = tau)

print(c('max difference, with exponential link', max(abs(pdf.dchi.eta - pdf.ncChi.eta))))

# exponential link with truncation
pdf.dchi.trunc <- dchi.trunc(test.vector, nu = nu, eta = eta, tau = tau, a = a)
pdf.ncChi.trunc <- ncChi.trunc(test.vector, nu = nu, eta = eta, tau = tau, a = a)

print(c('max difference, with truncation', max(abs(pdf.dchi.trunc - pdf.ncChi.trunc))))
