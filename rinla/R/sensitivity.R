## Export: inla.sens

##!\name{inla.sens}
##!
##!\title{Calculate sensitivity measurements}
##!
##!\description{TODO}
##!\usage{
##!   inla.sens(inlaRes)
##!}
##!
##!\arguments{
##!  \item{inlaRes}{Object returned by \code{inla} function.}
##!}
##!\value{
##!  \code{inla.sens}  plots robustness and returns object with different robustnesses
##!}
##!\author{Geir-Arne Fuglstad \email{geirarne.fuglstad@gmail.com}}
##!\examples{
##! TODO
##!}

inla.sens = function(inlaObj, lambda = 1){
    # Ensure that $misc$configs information is available
    if(!inlaObj$.args$control.compute$config){
        # Turn on storage of x|theta distributions
        inlaObj$.args$control.compute$config = TRUE

        # Use the previously calculated mode
        inlaObj$.args$control.mode$result = NULL
        inlaObj$.args$control.mode$restart = FALSE
        inlaObj$.args$control.mode$theta = inlaObj$mode$theta
        inlaObj$.args$control.mode$x = inlaObj$mode$x

        # Run INLA with all the settings
        inlaObj = do.call("inla", args = inlaObj$.args)
    }

    # Number of different theta values
    nTheta = inlaObj$misc$configs$nconfig

    # Number of latent components
    nLatent = length(inlaObj$misc$configs[[1]]$improved.mean)

    # Store parameters of distributions in matrices
    muMarg = matrix(0, nrow = nLatent, ncol = nTheta)
    sdMarg = muMarg
    skMarg = muMarg
    for(idxC in 1:nTheta){
        muMarg[, idxC] = inlaObj$misc$configs$config[[idxC]]$improved.mean
        sdMarg[, idxC] = sqrt(diag(inlaObj$misc$configs$config[[idxC]]$Qinv))
        skMarg[, idxC] = inlaObj$misc$configs$config[[idxC]]$skewness
    }

    ## Pre-compute robustification (Could be tabulated)
        # Resolution
        nSamples = 2e4

        # Correct lambda for number of parameters
        nPar = p + p*(p+1)/2
        lambda = lambda/sqrt(nPar)

        # Sample distances
        R = rexp(n = nSamples, rate = lambda)

        # Approximate the draws from the ellipsiod with draws from n-sphere
        z12 = matrix(rnorm(n = 2*nSamples), ncol = 2)
        zRest2 = rchisq(n = nSamples, df = nPar-2)
        sphereRadius = sqrt(rowSums(z12^2) + zRest2)
        sphereRadius = cbind(sphereRadius, sphereRadius)
        z12 = z12/sphereRadius

        # Scale with distances
        z12 = z12*R

        # Correct for different semi-axis
        z12[, 1] = z12[, 1]/1
        z12[, 2] = z12[, 2]/2

        # Extract mean and standard deviations
        muRobust = z12[, 1]
        sdRobust = exp(z12[, 1])

        # Choose interval to compute robustification of Gaussian on
        len = sqrt(var(muRobust) + mean(sdRobust^2))

        # Choose resolution for grid
        nGrid = 1e4

        # Calculate distribution
        

    return(inlaObj)
}
