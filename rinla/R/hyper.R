## this is a bit messy. (model, section) picks out the right model in
## the correct section. this model has a default setting for
## 'hyper'. simply replace these entries with given ones, first in
## 'hyper', then in 'initial', 'fixed', 'prior' and 'param'. somewhat
## tricky is that the length in 'param' depends on the argument(s) of
## 'prior' as various priors can have different number of
## parameters. also, if NA or NULL is given, then we use the default
## ones.

inla.set.hyper = function(
        model = NULL,
        section = names(inla.models()), 
        hyper = NULL, 
        initial = NULL,
        fixed = NULL,
        prior = NULL,
        param = NULL)
{
    ## name of the list can be in any CASE, and a longer name is ok,
    ## like 'list(precision=list(...))' where the 'name="prec"'.

    section = match.arg(section)
    hyper.new = inla.model.properties(model, section)$hyper
    nhyper = length(hyper.new)

    if (nhyper == 0) {
        if (!is.null(initial) && length(initial) > 0) {
            stop(inla.paste(c("Model:", model, "[", section, "], has none hyperparameters, but 'initial' is: ",
                              initial), sep=""))
        }
        if (!is.null(fixed) && length(fixed) > 0) {
            stop(inla.paste(c("Model:", model, "[", section, "], has none hyperparameters, but 'fixed' is:",
                              fixed), sep=""))
        }
        if (!is.null(param) && length(param) > 0) {
            stop(inla.paste(c("Model:", model, "[", section, "], has none hyperparameters, but 'param' is:",
                              param), sep=""))
        }
        if (!is.null(prior) && length(prior) > 0) {
            stop(inla.paste(c("Model:", model, "[", section, "], has none hyperparameters, but 'prior' is:",
                              param), sep=""))
        }
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
                stop(paste("Missing name/keyword in `hyper'; must be one of: ", inla.paste(valid.keywords),  ".",
                           sep=""))
            }
            if (!any(inla.strcasecmp(nm, valid.keywords))) {
                stop(paste("Unknown keyword in `hyper': `",  nm,
                           "'. Must be one of: ",  inla.paste(valid.keywords),  ".",  sep=""))
            }
        }
    }

    ## need `prior' to be before `param'!
    keywords = c("initial",  "fixed", "prior", "param")
    off = list(0, 0, 0, 0)
    names(off) = keywords

    ## do each hyperparameter one by time. fill into `hyper.new'
    for(ih in 1:nhyper) {
        
        ## find idx in the new one
        idx.new = which(names(hyper.new) == paste("theta", ih, sep=""))
        if (ih == 1 && length(idx.new) == 0) {
            idx.new = which(names(hyper.new) == "theta")
        }
        stopifnot(length(idx.new) > 0)
        name = hyper.new[[idx.new]]$name             ## full name
        short.name = hyper.new[[idx.new]]$short.name ## short name
        stopifnot(!is.null(short.name))

        if (!is.null(hyper)) {
            ## find idx in the given one
            idx = which(names(hyper) == paste("theta", ih, sep=""))
            if (length(idx) == 0 && nhyper == 1) {
                idx = which(names(hyper) == "theta")
            }
            if (length(idx) == 0) {
                idx = which(inla.strcasecmp(names(hyper), name))
            }
            if (length(idx) == 0) {
                idx = which(inla.strcasecmp(names(hyper), short.name))
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

            key.val = inla.eval(key)

            ## start of cmd is here
            if (!is.null(key.val)) {
                ## known length
                if (key == "param") {
                    len = inla.model.properties(hyper.new[[idx.new]]$prior, "prior")$nparameters
                    if (len < length(hyper.new[[idx.new]][[key]]) && len > 0) {
                        hyper.new[[idx.new]][[key]] = hyper.new[[idx.new]][[key]][1:len]
                    } else if (len < length(hyper.new[[idx.new]][[key]]) && len == 0) {
                        hyper.new[[idx.new]][[key]] = numeric(0)
                    } else if (len > length(hyper.new[[idx.new]][[key]])) {
                        hyper.new[[idx.new]][[key]] = c(hyper.new[[idx.new]][[key]],
                                                    rep(NA, len - length(hyper.new[[idx.new]][[key]])))
                    }
                } else {
                    len = length(hyper.new[[idx.new]][[key]])
                }
                ## given length
                llen = length(key.val) - off[[key]]
            
                if (llen < len) {
                    key.val = c(key.val,  rep(NA,  len - llen))
                }

                if (len > 0) {
                    ## set those NA's to the default ones
                    ii = key.val[off[[key]] + 1:len]
                    idxx = which(!(is.na(ii) | is.null(ii)))
                    if (length(idxx) > 0) {
                        hyper.new[[idx.new]][[key]][idxx] = ii[idxx]
                    }
                }
                off[[key]] = off[[key]] + len
            }

            if (!is.null(h) && !is.null(h[[key]]) && !is.na(h[[key]])) {
                ## known length
                if (key == "param") {
                    len = inla.model.properties(hyper.new[[idx.new]]$prior, "prior")$nparameters
                    if (len < length(hyper.new[[idx.new]][[key]]) && len > 0) {
                        hyper.new[[idx.new]][[key]] = hyper.new[[idx.new]][[key]][1:len]
                    } else if (len < length(hyper.new[[idx.new]][[key]]) && len == 0) {
                        hyper.new[[idx.new]][[key]] = NULL
                    } else if (len > length(hyper.new[[idx.new]][[key]])) {
                        hyper.new[[idx.new]][[key]] = c(hyper.new[[idx.new]][[key]],
                                                    rep(NA, len - length(hyper.new[[idx.new]][[key]])))
                    }
                } else {
                    len = length(hyper.new[[idx.new]][[key]])
                }
                ## given length
                llen = length(h[[key]])
            
                if (llen > len) {
                    stop(paste("model",  model, "[", section, "], hyperparam", ih,  ", length(hyper[[key]]) =",
                               llen,  ">",  len,  sep = " "))
                } else if (llen < len) {
                    h[[key]] = c(h[[key]],  rep(NA,  len - llen))
                }
                
                if (len > 0) {
                    ## set those NA's to the default ones
                    idxx = which(!is.na(h[[key]]) & !is.null(h[[key]]))
                    if (length(idxx) > 0) {
                        hyper.new[[idx.new]][[key]][idxx] = h[[key]][idxx]
                    }
                }
            }

            ans = inla.model.properties(hyper.new[[idx.new]]$prior, "prior", stop.on.error = TRUE)
        }
        
        if (key == "param") {
            if (length(hyper.new[[idx.new]]$param) != ans$nparameters) {
                stop(paste("Wrong length of prior-parameters, prior `", hyper.new[[idx.new]]$prior,  "' needs ",
                           ans$nparameters,  " parameters, you have ",
                           length(hyper.new[[idx.new]]$param),  ".",  sep=""))
            }
        }
    }
    

    for (key in keywords) {
        key.val = inla.eval(key)
        if (!is.null(key.val) && (length(key.val) > off[[key]])) {
            stop(paste("Length of argument `", key, "':", length(key.val),
                       ", does not match the total length of `", key, "' in `hyper':", off[[key]], sep=""))
        }
    }

    return (hyper.new)
}
