inla.get.hyperparameters.default = function(model) {
    return (list(theta1 = list(initial = 1,  fixed = TRUE, prior = "normal",  param = c(0, 1),
                         name = "Precision",  short.name = "prec"), 
                 theta2 = list(initial = 2,  fixed = FALSE, prior = "NORMAL",  param = c(10, 11),
                         name = "Lag 1 correlation",  short.name = "rho")))
}


inla.set.hyperparameters = function(model = NULL,  hyper = NULL, 
        initial = NULL, fixed = NULL,  prior = NULL,  param = NULL)
{
    ## name of the list can be in any CASE, and a longer name is ok,
    ## like 'list(precision=list(...))' where the 'name="prec"'.

    hyper.new = inla.get.hyperparameters.default(model)
    nhyper = length(hyper.new)
    if (nhyper == 0) {
        return (hyper.new)
    }
    
    if (!is.null(hyper)) {
        for (nm in names(hyper)) {
            valid.keywords = c(
                    "theta",
                    paste("theta",  1:nhyper,  sep=""),
                    as.character(sapply(hyper.new, function(x) x$short.name)), 
                    as.character(sapply(hyper.new, function(x) x$name)))

            if (is.null(nm)) {
                stop(paste("Missing name/keyword in `hyper'; must be one of: ", inla.paste(valid.keywords),  ".",  sep=""))
            }
            if (!any(inla.strncasecmp(nm, valid.keywords))) {
                stop(paste("Unknown keyword in `hyper': `",  nm,
                           "'. Must be one of: ",  inla.paste(valid.keywords),  ".",  sep=""))
            }
        }
    }

    ## need `prior' to be before `param'!
    keywords = c("initial",  "fixed", "prior", "param")

    for (key in keywords) {
        cmd = "<KEY>.off = 0L"
        inla.eval(gsub("<KEY>", key, cmd))
    }
    
    ## do each hyperparameter one by time. fill into `hyper.new'
    for(ih in 1:nhyper) {
        
        ## find idx in the new one
        idx.new = which(names(hyper.new) == paste("theta", ih, sep=""))
        if (ih == 1 && length(idx.new) == 0) {
            idx.new = which(names(hyper.new) == "theta")
        }
        stopifnot(length(idx.new) > 0)
        name = hyper.new[[idx.new]]$name         ## full name
        short.name = hyper.new[[idx.new]]$short.name ## short name
        stopifnot(!is.null(short.name))

        if (!is.null(hyper)) {
            ## find idx in the given one
            idx = which(names(hyper) == paste("theta", ih, sep=""))
            if (length(idx) == 0 && nhyper == 1) {
                idx = which(names(hyper) == "theta")
            }
            if (length(idx) == 0) {
                idx = which(inla.strncasecmp(names(hyper), name))
            }
            if (length(idx) == 0) {
                idx = which(inla.strncasecmp(names(hyper), short.name))
            }

            if (length(idx) == 0) {
                h = NULL
                idx = NULL
            } else {
                h = hyper[[idx]]
                if (!is.list(h)) {
                    stop(paste("Argument `hyper$theta", ih, "' or `",  short.name,
                               "' or `", name, "', is not a list: ", h,  sep=""))
                }
            }
        } else {
            h = NULL
            idx = NULL
        }

        for (key in keywords) {

            ## start of cmd is here
            ## check if there is an argument <KEY> as well
            if (!is.null(key)) {
                ## known length
                if (\"<KEY>\" == \"param\") {

                    len = inla.prior.properties(hyper.new[[idx.new]]$prior)$nparameters
                    if (len < length(hyper.new[[idx.new]]$<KEY>) && len > 0) {
                        hyper.new[[idx.new]]$<KEY> = hyper.new[[idx.new]]$<KEY>[1:len]
                    } else if (len < length(hyper.new[[idx.new]]$<KEY>) && len == 0) {
                        hyper.new[[idx.new]]$<KEY> = numeric(0)
                    } else if (len > length(hyper.new[[idx.new]]$<KEY>)) {
                        hyper.new[[idx.new]]$<KEY> = c(hyper.new[[idx.new]]$<KEY>, rep(NA, len - length(hyper.new[[idx.new]]$<KEY>)))
                    }
                } else {
                    len = length(hyper.new[[idx.new]]$<KEY>)
                }
                ## given length
                llen = length(<KEY>) - <KEY>.off
            
                if (llen < len) {
                    <KEY> = c(<KEY>,  rep(NA,  len - llen))
                }

                if (len > 0) {
                    ## set those NA's to the default ones

                    ii = <KEY>[<KEY>.off + 1:len]
                    idxx = which(!is.na(ii))
                    hyper.new[[idx.new]]$<KEY>[idxx] = ii[idxx]
                }
                <KEY>.off = <KEY>.off + len
            }

            if (!is.null(h) && !is.null(h$<KEY>)) {
                ## known length
                if (\"<KEY>\" == \"param\") {
                    len = inla.prior.properties(hyper.new[[idx.new]]$prior)$nparameters
                    if (len < length(hyper.new[[idx.new]]$<KEY>) && len > 0) {
                        hyper.new[[idx.new]]$<KEY> = hyper.new[[idx.new]]$<KEY>[1:len]
                    } else if (len < length(hyper.new[[idx.new]]$<KEY>) && len == 0) {
                        hyper.new[[idx.new]]$<KEY> = NULL
                    } else if (len > length(hyper.new[[idx.new]]$<KEY>)) {
                        hyper.new[[idx.new]]$<KEY> = c(hyper.new[[idx.new]]$<KEY>, rep(NA, len - length(hyper.new[[idx.new]]$<KEY>)))
                    }
                } else {
                    len = length(hyper.new[[idx.new]]$<KEY>)
                }
                ## given length
                llen = length(h$<KEY>)
            
                if (llen > len) {
                    stop(paste(\"model\",  model,  \", hyperparam\", ih,  \", length(hyper$<KEY>) =\",
                               llen,  \">\",  len,  sep = \" \"))
                } else if (llen < len) {
                    h$<KEY> = c(h$<KEY>,  rep(NA,  len - llen))
                }

                if (len > 0) {
                    ## set those NA's to the default ones
                    idxx = which(!is.na(h$<KEY>))
                    hyper.new[[idx.new]]$<KEY>[idxx] = h$<KEY>[idxx]
                }
            }
            ## end of cmd is here
            "

            ## evaluate it after doing the replacement
            inla.eval(gsub("<KEY>", key, cmd))

            ans = inla.prior.properties(hyper.new[[idx.new]]$prior, stop.on.error = TRUE)

            if (key == "param") {
                if (length(hyper.new[[idx.new]]$param) != ans$nparameters) {
                    stop(paste("Wrong length of prior-parameters, prior `", hyper.new[[idx.new]]$prior,  "' needs ",
                               ans$nparameters,  " parameters, you have ",  length(hyper.new[[idx.new]]$param),  ".",  sep=""))
                }
            }
        }
    }

    for (key in keywords) {
        cmd = "
        if (!is.null(<KEY>) && (length(<KEY>) > <KEY>.off)) {
            stop(paste(\"Length of argument `<KEY>':\", length(<KEY>),
                       \", does not match the total length of `<KEY>' in `hyper':\", <KEY>.off, sep=\"\"))
        }
        ## cmd ends here
        "

        ## evaluate it after doing the replacement
        inla.eval(gsub("<KEY>", key, cmd))
    }

    return (hyper.new)
}
