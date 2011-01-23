inla.get.hyper.default = function(model) {
    return (list(theta1 = list(initial = 1:6,  fixed = TRUE, prior = "normal",  param = c(0, 1)), 
                 theta2 = list(initial = 10:11,  fixed = FALSE, prior = "NORMAL",  param = c(10, 11))))
}


inla.set.hyper = function(model = NULL,  hyper = NULL, 
        initial = NULL, fixed = NULL,  prior = NULL,  param = NULL) {

    hyper.new = inla.get.hyper.default(model)
    nhyper = length(hyper.new)
    if (nhyper == 0) {
        return (hyper.new)
    }
    
    if (!is.null(hyper)) {
        for (nm in names(hyper)) {
            if (is.null(nm))
                stop("Missing name in `hyper'; must be `theta' or `theta1', `theta2', ...")
            if (!is.element(nm, c("theta",  paste("theta",  1:100,  sep=""))))
                stop(paste("Unknown name in `hyper': `",  nm,
                           "'. Must be `theta' or `theta1', `theta2', ...",  sep=""))
        }
    }

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
        hn = hyper.new[[idx.new]]

        if (!is.null(hyper)) {
            ## find idx in the given one
            idx = which(names(hyper) == paste("theta", ih, sep=""))
            if (nhyper == 1 && length(idx) == 0) {
                idx = which(names(hyper) == "theta")
            }
            if (length(idx) == 0) {
                ## cannot find name theta`ih'
                ## stop(paste("Cannot find name: `theta", ih,  "' in the `hyper' argument.", sep=""))
                h = NULL
                idx = NULL
            } else {
                h = hyper[[idx]]
                if (!is.list(h)) {
                    stop(paste("Argument `hyper$theta", ih, " is not a list: ", h,  sep=""))
                }
            }
        } else {
            h = NULL
            idx = NULL
        }

        for (key in keywords) {

            ## just write this using <KEY> and then evaluate it using <KEY> = key
            cmd = "
            ## start of cmd is here
            ## check if there is an argument <KEY> as well
            if (!is.null(<KEY>)) {
                ## known length
                len = length(hn$<KEY>)
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
                len = length(hn$<KEY>)
                ## given length
                llen = length(h$<KEY>)
            
                if (llen > len) {
                    stop(paste(\"model\",  model,  \"hyperparam\", ih,  \"length(hyper$<KEY>) =\",
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
        }
    }

    for (key in keywords) {
        cmd = "
        if (!is.null(<KEY>) && (length(<KEY>) > <KEY>.off)) {
            stop(paste(\"length of argument `<KEY>':\", length(<KEY>),
                       \", does not match the total length of `<KEY>' in `hyper':\", <KEY>.off, sep=\"\"))
        }
        ## cmd ends here
        "

        ## evaluate it after doing the replacement
        inla.eval(gsub("<KEY>", key, cmd))
    }

    return (hyper.new)
}
