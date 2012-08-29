`inla.rgeneric.define` = function(model.def = NULL, ...)
{
    model = list(
            definition = model.def,
            fifo = list(c2R = tempfile(), R2c = tempfile()), 
            args = list(...)
            )
    class(model) = "inla-rgeneric"
    return (model)
}

`inla.rgeneric.loop` = function(model, debug = FALSE)
{
    debug.cat = function(...) {
        if (debug)
            cat("R", ...)
    }
     
    stopifnot(inherits(model, "inla-rgeneric"))
    while (TRUE) {
        if (file.exists(model$fifo$R2c)) {
            R2c = fifo(model$fifo$R2c, "wb", blocking=TRUE)
            break
        }
    }
    debug.cat("Open file R2c", model$fifo$R2c, "\n")

    while (TRUE) {
        if (file.exists(model$fifo$c2R)) {
            c2R = fifo(model$fifo$c2R, "rb", blocking=TRUE)
            break
        }
    }
    debug.cat("Open file c2R", model$fifo$c2R, "\n")
    
    Id = NA
    n = -1L
    ntheta = -1L
    
    ## Enter the loop
    while (TRUE) {
        Id = readBin(c2R, what = integer())
        cmd = readBin(c2R, what = character())
        debug.cat("Id[", Id, "]. Got command:", cmd, "\n")

        if (cmd %in% "Q") {
            debug.cat("ntheta", ntheta, "\n")
            if (ntheta > 0L) {
                theta = readBin(c2R, what = numeric(), n = ntheta)
           } else {
                theta = numeric(0)
            }
            debug.cat("Got theta", theta, "\n")
            Q = do.call(model$definition, args = list(cmd = cmd, theta = theta, args = model$args))
            Q = inla.as.dgTMatrix(Q)
            stopifnot(dim(Q)[1L] == n && dim(Q)[2L] == n)
            len = length(Q@i)
            writeBin(as.integer(len), R2c)
            writeBin(as.integer(Q@i), R2c)
            writeBin(as.integer(Q@j), R2c)
            writeBin(Q@x, R2c)
        } else if (cmd %in% "graph") {
            G = do.call(model$definition, args = list(cmd = cmd, theta = NULL, args = model$args))
            G = inla.as.dgTMatrix(G)
            n = dim(G)[1L]
            len = length(G@i)
            debug.cat("Put len =", len, "\n")
            writeBin(as.integer(len), R2c)
            writeBin(as.integer(G@i), R2c)
            writeBin(as.integer(G@j), R2c)
        } else if (cmd %in% "initial") {
            init = do.call(model$definition, args = list(cmd = cmd, theta = NULL, args = model$args))
            ntheta = length(init)
            debug.cat("Put ntheta =", ntheta, "\n")
            debug.cat("Put initial", init, "\n")
            writeBin(as.integer(length(init)), R2c)
            if (length(init) > 0L) {
                writeBin(init, R2c)
            }
        } else if (cmd %in% "log.norm.const") {
            if (ntheta > 0L) {
                theta = readBin(c2R, what = numeric(), n = ntheta)
            } else {
                theta = numeric(0)
            }
            debug.cat("Got theta",  theta, "\n")
            lnc = do.call(model$definition, args = list(cmd = cmd, theta = theta, args = model$args))
            debug.cat("Put log.norm.const", lnc, "\n")
            writeBin(lnc, R2c)
        } else if (cmd %in% "log.prior") {
            if (ntheta > 0L) {
                theta = readBin(c2R, what = numeric(), n = ntheta)
            } else {
                theta = numeric(0)
            }
            debug.cat("Got theta",  theta, "\n")
            lp = do.call(model$definition, args = list(cmd = cmd, theta = theta, args = model$args))
            debug.cat("Put log.prior", lp, "\n")
            writeBin(lp, R2c)
        } else if (cmd %in% c("quit", "exit")) {
            debug.cat("Got a cleanup-and-exit-loop-message\n")
            close(R2c)
            unlink(R2c)
            close(c2R)
            unlink(c2R)
            break
        } else {
            stop(paste("Unknown command", cmd))
        }
    }

    return (invisible())
}

`inla.rgeneric.ar1.model` = function(
        cmd = c("Q", "graph", "initial", "log.norm.const", "log.prior"),
        theta = NULL,
        args = NULL)
{
    ## this is an example of the 'rgeneric' model. here we implement
    ## the AR-1 model as described in inla.doc("ar1"), where 'rho' is
    ## the lag-1 correlation and 'prec' is the *marginal* (not
    ## conditional) precision.

    interpret.theta = function(n, ntheta, theta)
    {
        ## internal helper-function to map the parameters from the internal-scale to the user-scale
        return (list(prec = exp(theta[1L]),
                     rho = 2*exp(theta[2L])/(1+exp(theta[2L])) - 1.0))
    }

    Q = function(n, ntheta, theta)
    {
        ## returns the precision matrix for given parameters
        param = interpret.theta(n, ntheta, theta)
        Q = param$prec/(1-param$rho^2) * toeplitz(c(1+param$rho^2, -param$rho, rep(0, n-2L)))
        ## first and last diagonal term is 1*marg.prec
        Q[1, 1] = Q[n, n] = param$prec/(1-param$rho^2)
        ## make it a sparse-matrix
        Q = inla.as.dgTMatrix(Q)
        return (Q)
    }

    graph = function(n, ntheta, theta)
    {
        ## return the graph of the model. the values of Q is only
        ## interpreted as zero or non-zero.
        Q = toeplitz(c(1, 1, rep(0, n-2L)))
        ## make it a sparse-matrix
        Q = inla.as.dgTMatrix(Q)
        return (Q)
    }

    log.norm.const = function(n, ntheta, theta)
    {
        ## return the log(normalising constant) for the model.
        param = interpret.theta(n, ntheta, theta)
        prec.innovation  = param$prec / (1.0 - param$rho^2)
        val = n * (- 0.5 * log(2*pi) + 0.5 * log(prec.innovation)) + 0.5 * log(1.0 - param$rho^2)
        return (val)
    }

    log.prior = function(n, ntheta, theta)
    {
        ## return the log-prior for the hyperparameters. the
        ## '+theta[1L]' is the log(Jacobian) for having a gamma prior
        ## on the precision and convert it into the prior for the
        ## log(precision).
        param = interpret.theta(n, ntheta, theta)
        val = (dgamma(param$prec, shape = 1, rate = 1, log=TRUE) + theta[1L] + 
               dnorm(theta[2L], mean = 0, sd = 1, log=TRUE))
        return (val)
    }

    initial = function(n, ntheta, theta)
    {
        ## return the initial values
        return (rep(1, ntheta))
    }

    cmd = match.arg(cmd)
    val = do.call(cmd, args = list(n = args$n, ntheta = args$ntheta, theta = theta))
    return (val)
}
