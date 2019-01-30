## NoExport: inla.sens
## NoExport: inla.sens.distance
## NoExport: inla.sens.distance.skew
## NoExport: inla.sens.skewMap
## NoExport: inla.sens.post

##!\name{inla.sens}
##!
##!\title{Calculate sensitivity measurements}
##!
##!\description{TODO}
##!\usage{
##!inla.sens(inlaObj, lambda = 0.3, nThreads = NULL, seed = NULL,
##!          nGrid = 1e4, nSamples = 2e4, nIntGrid = 1e4, useSkew = FALSE,
##!          calcPriorSens = FALSE, makePlots = TRUE)
##!}
##!
##!\arguments{
##!  \item{inlaRes}{Object returned by \code{inla} function.}
##!  \item{lambda}{TODO}
##!  \item{nThreads}{TODO}
##!  \item{seed}{TODO}
##!  \item{nGrid}{TODO}
##!  \item{nSamples}{TODO}
##!  \item{nIntGrid}{TODO}
##!  \item{useSkew}{TODO}
##!  \item{calcPriorSens}{TODO}
##!  \item{makePlots}{TODO}
##!}
##!\value{
##!  \code{inla.sens}  plots robustness and returns object with different robustnesses
##!}
##!\author{Geir-Arne Fuglstad \email{geirarne.fuglstad@gmail.com}}
##!\examples{
##! TODO
##!}

inla.sens = function(inlaObj, lambda = 0.3, nThreads = NULL, seed = NULL,
                     nGrid = 1e4, nSamples = 2e4, nIntGrid = 1e4, useSkew = FALSE,
                     calcPriorSens = FALSE, makePlots = TRUE)
{
    ## Ensure reproducability
    if(!is.null(seed))
        set.seed(seed)

    ## Ensure that $misc$configs information is available
    if(!inlaObj$.args$control.compute$config){
        ## Turn on storage of x|theta distributions
        inlaObj$.args$control.compute$config = TRUE

        ## Use the previously calculated mode
        inlaObj$.args$control.mode$result = NULL
        inlaObj$.args$control.mode$restart = FALSE
        inlaObj$.args$control.mode$theta = inlaObj$mode$theta
        inlaObj$.args$control.mode$x = inlaObj$mode$x

        ## Run INLA with all the settings
        inlaObj = do.call("inla", args = inlaObj$.args)
    }

    ## Re-run with all data set to NA to get prior-level
    if(calcPriorSens){
        inlaObjRef = inlaObj
        inlaObjRef$.args$data$y = NA*inlaObjRef$.args$data$y
        inlaObjRef = do.call("inla", args = inlaObjRef$.args)
        priorSens = inla.sens(inlaObjRef, lambda = lambda, nThreads = nThreads,
                              seed = seed, nGrid = nGrid, nSamples = nSamples,
                              nIntGrid = nIntGrid, useSkew = useSkew, calcPriorSens = FALSE)
    }

    ## Number of different theta values
    nTheta = inlaObj$misc$configs$nconfig

    ## Number of latent components
    nLatent = length(inlaObj$misc$configs$config[[1]]$improved.mean)

    ## Store parameters of distributions in matrices
    muMarg = matrix(0, nrow = nTheta, ncol = nLatent)
    sdMarg = muMarg
    skMarg = muMarg
    prob = vector(mode = "numeric", length = nTheta)
    for(idxC in 1:nTheta){
        muMarg[idxC, ] = inlaObj$misc$configs$config[[idxC]]$improved.mean
        sdMarg[idxC, ] = sqrt(diag(inlaObj$misc$configs$config[[idxC]]$Qinv))
        skMarg[idxC, ] = inlaObj$misc$configs$config[[idxC]]$skewness
        prob[idxC] = exp(inlaObj$misc$configs$config[[idxC]]$log.posterior)
    }
    prob = prob/sum(prob)

    ## Pre-compute robustification (Could be tabulated)
    ## Using the exponential distribution too easily leads to a distribution for the standard deviations
    ## that has no mean

    ## Correct lambda for number of parameters
    ## nPar = nLatent + nLatent*(nLatent+1)/2
    ## lambda = lambda/sqrt(nPar)

    ## ## Make sure mean and variance of standard deviations will exist


    ## Sample distances
    ## R = rexp(n = nSamples, rate = lambda)

    ## With the Gaussian distribution we can ensure the correct behaviour
    ## Scale variance according to number of samples
    nPar = nLatent + nLatent*(nLatent+1)/2
    sigDist = lambda*sqrt(nPar-3)

    ## Simulate distances
    R = abs(rnorm(nSamples, sd = sigDist))

    ## Approximate the draws from the ellipsiod with draws from n-sphere
    z12 = matrix(rnorm(n = 2*nSamples), ncol = 2)
    zRest2 = rchisq(n = nSamples, df = nPar-2)
    sphereRadius = sqrt(rowSums(z12^2) + zRest2)
    sphereRadius = cbind(sphereRadius, sphereRadius)
    z12 = z12/sphereRadius

    ## Scale with distances
    z12 = z12*R

    ## Correct for different semi-axis
    z12[, 1] = z12[, 1]/1
    z12[, 2] = z12[, 2]/2

    ## Extract mean and standard deviations
    muRobust = z12[, 1]
    sdRobust = exp(z12[, 2])

    ## Choose interval to compute robustification of Gaussian on
    len = sqrt(var(muRobust) + mean(sdRobust^2))

    ## Calculate distribution
    xR = seq(-20*len, 20*len, length.out = nGrid)
    yR = vector(mode = "numeric", length = nGrid)
    for(idxS in 1:nSamples){
        yR = yR + dnorm(xR, mean = muRobust[idxS], sd = sdRobust[idxS])
    }
    yR = log(yR)-log(nSamples)
    robMarg = list(x = xR, y = yR)

    ## Calculate distance between original marginal posterior and
    ## posterior with uncertainty added in each conditional density
    ## x_i | \theta, y
    inla.require("doParallel")
    if(is.null(nThreads)){
        doParallel::registerDoParallel()
    } else{
        doParallel::registerDoParallel(nThreads)
    }
    inla.require("foreach")
    nWorkers = foreach::getDoParWorkers()
    breaks = floor(seq(1, nLatent+1, length.out = nWorkers+1))
    ds = foreach::foreach(idxW = 1:nWorkers, .combine = 'c') %dopar%{
        tmpRes = vector(mode = "numeric", length = breaks[idxW+1]-breaks[idxW])
        for(idx in breaks[idxW]:(breaks[idxW+1]-1)){
            if(useSkew){
                tmpRes[idx-breaks[idxW]+1] = inla.sens.distance.skew(
                    muMarg[, idx], sdMarg[, idx], skMarg[, idx], prob, robMarg, nIntGrid)
            } else{
                tmpRes[idx-breaks[idxW]+1] = inla.sens.distance(
                    muMarg[, idx], sdMarg[, idx], skMarg[, idx], prob, robMarg, nIntGrid)
            }
        }

        tmpRes
    }

    ## Calculate max distance
    dMax = inla.sens.distance(0, 1, 0, 1, robMarg, nIntGrid)

    ## Standardize against max distance
    stdDist = (dMax-ds)/dMax

    ## Store results in groups
    res = list()

    ## Make one plot for each group of variables
    groups = inlaObj$misc$configs$contents
    nGroups = length(groups$tag)
    xLab = c()
    val  = c()
    cex.names = 1.2
    cex.axis  = 1.2
    lwd       = 1.2
    for(idxP in 1:nGroups){
        sIdx = groups$start[idxP]
        eIdx = sIdx + groups$length[idxP]-1
        if(groups$length[idxP] == 1){
            xLab = c(xLab, groups$tag[idxP])
            val  = c(val,  stdDist[sIdx])
        } else{
            if(makePlots){
                inla.dev.new()
                barplot(stdDist[sIdx:eIdx], 
                        space = 2,
                        ylim = c(0, 1), 
                        main = groups$tag[idxP], 
                        names.arg = 1:(eIdx-sIdx+1), 
                        xlab = "Index", 
                        ylab = "Uncertainty",
                        cex.names = cex.names,
                        cex.axis  = cex.axis,
                        lwd       = lwd)

                ## Add prior level to plot if calculated
                if(calcPriorSens){
                    barplot(priorSens[[length(res)+1]]$val, 
                            add = TRUE, 
                            col = rgb(1, 0, 0, alpha = .20),
                            space = 2,
                            cex.axis = cex.axis)
                }
            }

            ## Add groups into result object
            res = c(res, list(list(tag = groups$tag[idxP],
                                   val = stdDist[sIdx:eIdx])))
        }
    }
    if(length(val) >= 1){
        if(makePlots){
            inla.dev.new()
            barplot(val, 
                    space = 2, 
                    ylim = c(0, 1), 
                    main = "Fixed effects", 
                    names.arg = xLab, 
                    ylab = "Uncertainty",
                    cex.names = cex.names,
                    cex.axis  = cex.axis,
                    lwd       = lwd)

            ## Add prior level to plot if calculated
            if(calcPriorSens){
                barplot(priorSens[[length(res)+1]]$val, 
                        add = TRUE, 
                        col = rgb(1, 0, 0, alpha = .20),
                        space = 2,
                        cex.axis = cex.axis)
            }
        }

        ## Add fixed effects to result object
        res = c(res, list(list(tag = "Fixed effects",
                               val = val,
                               names = xLab)))
    }

    return(res)
}

inla.sens.distance = function(muMarg, sdMarg, skMarg, prob, robMarg, nGrid, extraLen = 20)
{
    ## Estimate required integration grid
    sdMax = max(sdMarg)
    xs = seq(-1, 1, length.out = nGrid)*sdMax*extraLen + mean(muMarg)

    ## Iterate through \theta values
    yO = vector(mode = "numeric", length = nGrid)
    yR = yO
    for(idxT in 1:length(muMarg)){
        ## Use precomputed table of standard robust distribution
        xx = (xs - muMarg[idxT])/sdMarg[idxT]
        yy = exp(spline(x = robMarg$x, y = robMarg$y, xout = xx)$y)/sdMarg[idxT]

        ## Remove the extrapolated values
        yy[(xx < min(robMarg$x)) | (xx > max(robMarg$x))] = 0

        ## Add original and robust to their respective mixtures
        yO = yO + prob[idxT]*dnorm(xs, mean = muMarg[idxT], sd = sdMarg[idxT])
        yR = yR + prob[idxT]*yy
    }

    ##  Calculate KLD between original and added uncertainty
    intG = yR*log(yR/yO)
    intG[yR == 0] = 0
    KLD = sum((intG[-nGrid] + intG[-1])*(xs[2]-xs[1])/2)

    ## Convert to distance and return value
    return(sqrt(2*KLD))
}

inla.sens.distance.skew = function(muMarg, sdMarg, skMarg, prob, robMarg, nGrid, extraLen = 20)
{
    inla.require("sn")
    ## Estimate required integration grid
    sdMax = max(sdMarg)
    xs = seq(-1, 1, length.out = nGrid)*sdMax*extraLen + mean(muMarg)

    ## Iterate through \theta values
    yO = vector(mode = "numeric", length = nGrid)
    yR = yO
    for(idxT in 1:length(muMarg)){
        ## Extract parameters of skew normal marginal
        mu  = muMarg[idxT]
        sig = sdMarg[idxT]
        gamma = skMarg[idxT]
        if(is.nan(gamma))
            gamma = 0

        ## Convert to standard parametrization of skew normal
        p = inla.sens.skewMap(c(mu, sig, gamma))

        ## Map to standard distribution
        sIdx = 1
        mIdx = floor(nGrid/2)
        eIdx = nGrid
        tmpXX1 = sn::psn(xs[sIdx:(mIdx-1)], xi = p[1], omega = p[2], alpha = p[3])
        tmpXX2 = sn::psn(-xs[mIdx:eIdx], xi = -p[1], omega = p[2], alpha = -p[3])
        xx = qnorm(tmpXX1)
        xx2 = -qnorm(tmpXX2)
        xx = c(xx, xx2)

        ## Use precomputed table of standard robust distribution
        yy = exp(spline(x = robMarg$x, y = robMarg$y, xout = xx)$y)

        ## Correct for transformation
        yy = yy*exp(sn::dsn(xs, xi = p[1], omega = p[2], alpha = p[3], log = TRUE)-dnorm(xx, log = TRUE))

        ## Remove the extrapolated values
        yy[(xx < min(robMarg$x)) | (xx > max(robMarg$x))] = 0

        ## Add original and robust to their respective mixtures
        yO = yO + prob[idxT]*sn::dsn(xs, xi = p[1], omega = p[2], alpha = p[3])
        yR = yR + prob[idxT]*yy
    }

    ##  Calculate KLD between original and added uncertainty
    intG = yR*log(yR/yO)
    intG[yR == 0] = 0
    KLD = sum((intG[-nGrid] + intG[-1])*(xs[2]-xs[1])/2)

    ## Convert to distance and return value
    return(sqrt(2*KLD))
}

inla.sens.skewMap = function(x){
    ## Extract desired parameters of skew normal
    mu = x[1]
    s = x[2]
    g1 = x[3]

    ## Transform to usual parametrization of skew normal
    d = sign(g1)*sqrt(abs(g1)^(2/3)/(2/pi*(((4-pi)/2)^(2/3)+abs(g1)^(2/3))))
    alpha = d/sqrt(1-d^2)
    w = s*(1-2*d^2/pi)^(-0.5)
    ksi = mu - w*d*sqrt(2/pi)

    return(c(ksi, w, alpha))
}

inla.sens.post = function(res.sens, rIdx, X)
{
    ## Number of components
    nC = length(res.sens)
    nP = length(res.sens[[1]]$val)

    ## Iterate through each posterior and plot
    par(ask = TRUE)
    for(idxP in 1:nP){
        ## Get values and names
        val = NULL
        nam = NULL
        for(idxC in 1:nC){
            if(res.sens[[idxC]]$tag == "Fixed effects"){
                tmpIdx = (X[idxP, ] != 0)
                val = c(val, res.sens[[idxC]]$val[tmpIdx])
                nam = c(nam, res.sens[[idxC]]$names[tmpIdx])
            } else{
                if(!is.na(rIdx[idxP, idxC])){
                    val = c(val, res.sens[[idxC]]$val[rIdx[idxP, idxC]])
                    nam = c(nam, res.sens[[idxC]]$tag)
                }
            }
        }

        ## Standardize with linear predictor
        val = val/val[1]
        inla.dev.new()
        cex.names = 1.5
        cex.axis = 1.5
        lwd = 1.5
        barplot(val, 
                space = 2,
                ylim = c(1/exp(3), exp(3)), 
                main = paste("Linear predicter", idxP), 
                names.arg = nam, 
                xlab = "Index", 
                ylab = "Resiliency",
                cex.names = cex.names,
                cex.axis  = cex.axis,
                lwd       = lwd,
                log = c("y"))
    }
    par(ask = FALSE)
}
