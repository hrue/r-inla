`inla.R.generic.define` = function(model.def = NULL, ...)
{
    model = list(
            definition = model.def,
            fifo = list(c2R = tempfile(), R2c = tempfile()), 
            args = list(...)
            )
    class(model) = "inla-R-generic"
    return (model)
}

`inla.R.generic.loop` = function(model, debug = TRUE)
{
    debug.cat = function(...) {
        if (debug)
            cat("R", ...)
    }
    
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
            Q = model$definition(cmd, input = theta, args = model$args)
            Q = inla.as.dgTMatrix(Q)
            stopifnot(dim(Q)[1L] == n && dim(Q)[2L] == n)
            len = length(Q@i)
            writeBin(as.integer(len), R2c)
            writeBin(as.integer(Q@i), R2c)
            writeBin(as.integer(Q@j), R2c)
            writeBin(Q@x, R2c)
        } else if (cmd %in% "graph") {
            G = model$definition(cmd, input = NULL, args = model$args)
            G = inla.as.dgTMatrix(G)
            n = dim(G)[1L]
            len = length(G@i)
            debug.cat("Put len =", len, "\n")
            writeBin(as.integer(len), R2c)
            writeBin(as.integer(G@i), R2c)
            writeBin(as.integer(G@j), R2c)
        } else if (cmd %in% "initial") {
            initial = model$definition(cmd, input = NULL, args = model$args)
            ntheta = length(initial)
            debug.cat("Put ntheta =", ntheta, "\n")
            debug.cat("Put initial", initial, "\n")
            writeBin(as.integer(length(initial)), R2c)
            if (length(initial) > 0L) {
                writeBin(initial, R2c)
            }
        } else if (cmd %in% "extra") {
            if (ntheta > 0L) {
                theta = readBin(c2R, what = numeric(), n = ntheta)
            } else {
                theta = numeric(0)
            }
            debug.cat("Got theta",  theta, "\n")
            extra = model$definition("extra", input = theta, args = model$args)
            debug.cat("Put extra", c(extra$model, extra$prior), "\n")
            writeBin(extra$model, R2c)
            writeBin(extra$prior, R2c)
        } else if (cmd %in% c("quit", "exit", "break")) {
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

`inla.R.generic.ar1.model` = function(
        what = list("Q", "graph", "initial", "extra"),
        input = NULL,
        args)
{
    ## AR1 model with 

    interpret.theta = function(theta)
    {
        return (list(prec = exp(theta[1L]),
                     rho = 2*exp(theta[2L])/(1+exp(theta[2L])) - 1.0))
    }

    Qfunc = function(theta, n)
    {
        param = interpret.theta(theta)
        Q = toeplitz(c(1+param$rho^2, -param$rho, rep(0, n-3L), -param$rho))
        Q[1, n] = Q[n, 1] = 0
        Q[1, 1] = Q[n, n] = 1
        prec.noise = param$prec/(1-param$rho^2)
        Q = prec.noise * Q
        Q = inla.as.dgTMatrix(Q)
        ##print("AR1  Q-matrix"); print(Q)
        return (Q)
    }

    graph = function(n)
    {
        Q = toeplitz(c(1, 1, rep(0, n-3L), 1))
        Q[1, n] = Q[n, 1] = 0
        Q = inla.as.dgTMatrix(Q)
        return (Q)
    }

    extra = function(theta, n)
    {
        param = interpret.theta(theta)
        prec.innovation  = param$prec / (1.0 - param$rho^2)

        val = list(
                model = n * (- 0.5 * log(2*pi) + 0.5 * log(prec.innovation)) + 0.5 * log(1.0 - param$rho^2), 
                prior = (dgamma(param$prec, shape = 1, rate = 1, log=TRUE) + theta[1L] + 
                         dnorm(theta[2L], mean = 0, sd = 1, log=TRUE)
                         )
                )
        return (val)
    }

    initial = function()
    {
        return (c(-1, 1))
    }

    what = match.arg(what)

    n = args$n
    ntheta = args$ntheta
    
    if (what %in% "Q") {
        return (Qfunc(input, n))
    } else if (what %in% "graph") {
        return (graph(n))
    } else if (what %in% "initial") {
        return (initial())
    } else if (what %in% "extra") {
        aa = extra(input, n)
        return (extra(input, n))
    } else {
        stop("This shuld not happen.")
    }
}
