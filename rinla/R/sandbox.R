## Nothing to export

inla.extract.Q <- function(what, result) {
    stopifnot(!missing(what))
    stopifnot(!missing(result))

    conf <- result$misc$configs
    cont <- conf$contents
    i <- which(cont$tag == what)
    if (length(i) != 1L) {
        stop(paste0("Cannot find '", what, "'"))
    }
    p.len <- cont$length[1]
    i.start <- cont$start[i]
    i.len <- cont$length[i]

    k <- NA
    for (j in seq_along(conf$config)) {
        if (conf$config[[j]]$log.posterior.orig == 0.0) {
            k <- j
            break
        }
    }
    stopifnot(!is.na(k))
    Q <- conf$config[[k]]$Q
    diag(Q) <- conf$config[[k]]$Qprior.diag
    d <- diag(Q)
    log.prec <- inla.models()$predictor$predictor$hyper$theta$initial
    browser()
    for (i in i.start:(i.start + i.len - 1L)) {
        Q[i, i] <- Q[i, i] - sum(Q[1:p.len, i]^2 / (d[1:p.len] - exp(log.prec)))
    }
    Q <- Q[i.start:(i.start + i.len - 1), i.start:(i.start + i.len - 1)]
    Q <- Q + t(Q)
    diag(Q) <- diag(Q) / 2.0

    return(Q)
}
