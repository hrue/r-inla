`inla.create.data.file` =
    function(y.orig = NULL,
             MPredictor = NULL,
             mf=NULL,
             scale=NULL,
             E=NULL,
             Ntrials=NULL,
             event=NULL,
             family=NULL,
             data.dir=NULL,
             file=NULL,
             debug=FALSE)
{
    if (is.null(y.orig)) {
        y.orig = c(mf[,1])
    } else {
        if (!inherits(y.orig,"inla.surv"))
            y.orig = as.data.frame(y.orig)
    }
    
    ##FIXME
    ##n.data=length(mf[,1])
    n.data = dim(y.orig)[1]
    ind=seq(0,n.data-1)

    if (debug)
        cat("inla.create.data.file: n.data = ", n.data, "\n")

    ## RESPONSE (gaussian and T can have the weights for the precision: y.i\sim N(\eta,prec\w.i))
    if (inla.one.of(family, c("gaussian", "normal", "t", "laplace", "sn", "skewnormal", "gev"))) {
        if(is.null(scale)) {
            response=cbind(ind,rep(1,n.data), y.orig)
            null.dat = is.na(response[,3])
            response = response[!null.dat,]
        } else {
            response = cbind(ind, y.orig)
            null.dat = is.na(response[,2])
            response = response[!null.dat,]
            scale = scale[!null.dat]

            if (debug)
                cat(paste("length(scale), nrow(response)",length(scale),nrow(response),"\n",sep=" "  ))
            
            if(length(scale)!=nrow(response)) {
                file.remove(file)
                file.remove(data.dir)
                stop("Length of scale has to be the same as the length of data")
            }
            else
                response=cbind(response[,1],scale,response[,2])
        }   
    } else if (inla.one.of(family, c("poisson",
                                     "zeroinflatedpoisson0",
                                     "zeroinflatedpoisson1", 
                                     "zeroinflatedpoisson2", 
                                     "nbinomial",
                                     "zeroinflatednbinomial0",
                                     "zeroinflatednbinomial1",
                                     "zeroinflatednbinomial2"))) {
        ## Poisson/nbinomial family has the E field
        if (is.null(E)) {
            response = cbind(ind,rep(1,n.data),y.orig)
            null.dat = is.na(response[,3])
            response = response[!null.dat,]
        } else {
            response = cbind(ind, y.orig)
            null.dat = is.na(response[,2])
            response = response[!null.dat,]
            E = E[!null.dat]
            if(length(E)!=nrow(response)) {
                file.remove(file)
                file.remove(data.dir)
                stop("Length of E has to be the same as the length of data")
            }
            else
                response=cbind(response[,1], E, response[,2])
        }       
    } else if (inla.one.of(family, c("poissonext"))) {
        ## PoissonExt family has the E field, which is a n x 3 matrix
        if (is.null(E)) {
            ## this is the default value '[1,0,0]'
            response = cbind(ind, rep(1,n.data), rep(0, n.data), rep(0, n.data), y.orig)
            null.dat = is.na(response[,3])
            response = response[!null.dat,]
        } else {
            response = cbind(ind, y.orig)
            null.dat = is.na(response[,2])
            response = response[!null.dat,]
            if (!is.matrix(E) || dim(E)[2] != 3)
                 stop(paste("Family 'poissonext': E has to be a 'n x 3' matrix;", dim(E)))
            E = E[!null.dat, ]
            if (dim(E)[1] != nrow(response)) {
                file.remove(file)
                file.remove(data.dir)
                stop(paste("Family 'poissonext': dim(E)[1] has to be the same as the length of data;", dim(E)[1], nrow(response)))
            }
            else
                response=cbind(response[,1], E, response[,2])
        }       
    } else if(inla.one.of(family,
                          c("binomial",
                            "zeroinflatedbinomial0",
                            "zeroinflatedbinomial1",
                            "zeroinflatedbinomial2",
                            "zeroinflatedbetabinomial2"))) {
        response = cbind(ind, y.orig)
        null.dat = is.na(response[,2])
        response = response[!null.dat,]
        if (is.null(Ntrials)) {
            Ntrials = rep(1,n.data)
        } else {
            Ntrials =  Ntrials[!null.dat]
        }
        if(is.null(Ntrials)) {
            file.remove(file)
            file.remove(data.dir)
            stop("Number of binomial trials has to be provided")
        } else {
            if(sum(Ntrials >= response[,2]) < nrow(response)) {
                file.remove(file)
                file.remove(data.dir)
                stop("Nnumber of trials have to be larger than y")
            }
            else
                response=cbind(response[,1],Ntrials,response[,2])
        }
    } else if (inla.one.of(family, c("exponential", "weibull", "weibullcure"))) {
        if (!inla.lmodel.properties(family)$survival)
            stop("This should not happen.")

        if (is.null(y.orig$truncation))
            stop("Responce does not contain variable `truncation'.")
        if (is.null(y.orig$lower))
            stop("Responce does not contain variable `lower'.")
        if (is.null(y.orig$upper))
            stop("Responce does not contain variable `upper'.")
        if (is.null(y.orig$event))
            stop("Responce does not contain variable `event'.")
        if (is.null(y.orig$time))
            stop("Responce does not contain variable `time'.")

        idx = !is.na(y.orig$time)
        response = cbind(ind[idx], y.orig$event[idx], y.orig$truncation[idx], y.orig$lower[idx], y.orig$upper[idx], y.orig$time[idx])

        if (any(is.na(response)))
            stop("NA in truncation/event/lower/upper/time is not allowed")
    } else if (inla.one.of(family, c("stochvol", "stochvolt", "stochvolnig"))) {
        response = cbind(ind, y.orig)
        null.dat = is.na(response[,2])
        response = response[!null.dat,]
    } else {
        file.remove(file)
        file.remove(data.dir)
        stop(paste("response family:", family, ", not recognised"))
    }
    ##create the DATA file
    file.data = inla.tempfile(tmpdir=data.dir)
    file.create(file.data)
    write(t(response),ncolumns=ncol(response),file=file.data,append=FALSE)

    file.data = gsub(data.dir, "$inladatadir", file.data, fixed=TRUE)

    return(file.data)
}
