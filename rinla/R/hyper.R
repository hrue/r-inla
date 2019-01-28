## Nothing to Export.

## this is a bit messy. (model, section) picks out the right model in
## the correct section. this model has a default setting for
## 'hyper'. simply replace these entries with given ones, first in
## 'hyper', then in 'initial', 'fixed', 'prior' and 'param'. somewhat
## tricky is that the length in 'param' depends on the argument(s) of
## 'prior' as various priors can have different number of
## parameters. also, if NA or NULL is given, then we use the default
## ones.

## I think its a better strategy to take initial, fixed, prior, param
## and make a hyper-similar structure and then merge these two.  the
## current code is messy but does not have to be.

inla.set.hyper = function(
        model = NULL,
        section = NULL, 
        hyper = NULL,
        initial = NULL,
        fixed = NULL,
        prior = NULL,
        param = NULL,
        debug = FALSE,
        hyper.default = NULL)
{
    ## name of the list can be in any CASE, and a longer name is ok,
    ## like 'list(precision=list(...))' where the 'name="prec"'.
    skip.final.check = FALSE

    if (is.null(section)) {
        stop("No section given; please fix...")
    }
    mm = inla.models()
    section = match.arg(section, names(mm))    

    if (is.null(hyper.default)) {
        ## default values are given in the inla.models()
        hyper.new = inla.model.properties(model, section)$hyper
    } else {
        ## default values are given by the user. process these
        hyper.new = inla.set.hyper(model, section, hyper = hyper.default,  debug = debug)
    }
    nhyper = length(hyper.new)

    if (debug) {
        cat(paste("* Get default hyper from model",  model,  "in section",  section), "\n")
    }

    if (nhyper == 0L) {
        if (!is.null(hyper) && !identical(hyper,  list())) {
            stop(inla.paste(c("Model", model, "[", section,
                              "], has none hyperparameters, but 'hyper' is ",
                              inla.paste(list(hyper))), sep = " "))
        }
        if (!is.null(initial) && length(initial) > 0L) {
            stop(inla.paste(c("Model", model, "[", section, "], has none hyperparameters, but 'initial' is ",
                              initial), sep=" "))
        }
        if (!is.null(fixed) && length(fixed) > 0L) {
            stop(inla.paste(c("Model", model, "[", section, "], has none hyperparameters, but 'fixed' is",
                              fixed), sep=" "))
        }
        if (!is.null(param) && length(param) > 0L) {
            stop(inla.paste(c("Model", model, "[", section, "], has none hyperparameters, but 'param' is",
                              param), sep=" "))
        }
        if (!is.null(prior) && length(prior) > 0L) {
            stop(inla.paste(c("Model", model, "[", section, "], has none hyperparameters, but 'prior' is",
                              param), sep=" "))
        }
        return (hyper.new)
    }

    if (!is.null(hyper)) {
        for (nm in names(hyper)) {
            if (debug) {
                cat(paste("* Check contents of given hyper$", nm,  sep=""), "\n")
            }
            valid.keywords = c(
                    "theta",
                    paste("theta",  1L:nhyper,  sep=""),
                    as.character(sapply(hyper.new, function(x) x$short.name)),
                    as.character(sapply(hyper.new, function(x) x$name)))

            if (is.null(nm)) {
                stop(paste("Missing name/keyword in `hyper'; must be one of ", inla.paste(valid.keywords),  ".",
                           sep=" "))
            }
            if (!any(inla.strcasecmp(nm, valid.keywords))) {
                stop(paste("Unknown keyword in `hyper' `",  nm,
                           "'. Must be one of ",  inla.paste(valid.keywords),  ".",  sep=" "))
            }
        }
    }

    ## need `prior' to be before `param'!
    keywords = c("initial",  "fixed", "prior", "param")
    off = list(0L, 0L, 0L, 0L)
    names(off) = keywords

    ## do each hyperparameter one by time. fill into `hyper.new'
    for(ih in 1L:nhyper) {

        if (debug) {
            cat(paste("** Check hyperparameter",  ih,  "with name", names(hyper.new)[ih]), "\n")
        }

        ## find idx in the new one
        idx.new = which(names(hyper.new) == paste("theta", ih, sep=""))
        if (ih == 1L && length(idx.new) == 0L) {
            idx.new = which(names(hyper.new) == "theta")
        }
        stopifnot(length(idx.new) > 0L)
        name = hyper.new[[idx.new]]$name             ## full name
        short.name = hyper.new[[idx.new]]$short.name ## short name
        stopifnot(!is.null(short.name))

        if (debug) {
            cat(paste("** idx.new =",  idx.new), "\n")
        }

        if (!is.null(hyper)) {
            ## find idx in the given one
            idx = which(names(hyper) == paste("theta", ih, sep=""))
            if (length(idx) == 0L && nhyper == 1L) {
                idx = which(names(hyper) == "theta")
            }
            if (length(idx) == 0L) {
                idx = which(inla.strcasecmp(names(hyper), name))
            }
            if (length(idx) == 0L) {
                idx = which(inla.strcasecmp(names(hyper), short.name))
            }

            if (length(idx) == 0L) {
                h = NULL
                idx = NULL
            } else {
                h = hyper[[idx]]
                if (!is.list(h)) {
                    stop(paste("Argument `hyper$theta", ih, "' or `",  short.name,
                               "' or `", name, "', is not a list: ", h,  sep=""))
                }
            }

            if (debug) {
                cat(paste("** Corresponding idx in given hyper =",  idx), "\n")
                if (is.null(h)) {
                    cat(paste("*** h is NULL"), "\n")
                } else {
                    cat(paste("*** h is !NULL"), "\n")
                }
            }

        } else {
            h = NULL
            idx = NULL
        }

        for (key in keywords) {

            key.val = inla.eval(key)
            if (debug) {
                cat(paste("*** Process key =",  key), "\n")
                cat(inla.paste("*** key =", inla.ifelse(is.null(key.val) || (is.numeric(key.val) && length(key.val)==0),
                                                        "NULL",  as.character(key.val))), "\n")
            }

            ## special case. if attribute is set to 'read only' then refuse to change value
            read.only = attr(hyper.new[[idx.new]][[key]], "inla.read.only")
            if (!is.null(read.only) && read.only == TRUE) {
                if (!is.null(h[[key]]) && h[[key]] != hyper.new[[idx.new]][[key]])
                    stop(paste("Setting hyperparameter `", name,
                               "'. Key `", key, "' = '", hyper.new[[idx.new]][[key]],
                               "' is 'read-only', cannot change to `", h[[key]], "'.", sep=""))
                ii = idx.new + off[[key]]
                if (!is.null(key.val) && !is.na(key.val[ii])) {
                    if (key.val[ii] != hyper.new[[idx.new]][[key]]) {
                        stop(paste("Setting hyperparameter `", name,
                                   "'. Key `", key,"' = '", hyper.new[[idx.new]][[key]],
                                   "' is 'read-only', cannot change to `", key.val[ii], "'.", sep=""))
                    }
                }
            }

            ## start of cmd is here
            if (!is.null(key.val)) {
                ## known length
                if (key == "param") {
                    len = inla.model.properties(hyper.new[[idx.new]]$prior, "prior")$nparameters
                    if (len < 0L) {
                        ## this is a special case, where the number of
                        ## parameters in the prior vary. first example
                        ## is the spde2/spde3 model. in this case, we just
                        ## eat all the remaining parameters.
                        skip.final.check = TRUE
                    } else {
                        ## this is the normal case
                        if (len <= length(hyper.new[[idx.new]][[key]])) {
                            if (len > 0L) {
                                hyper.new[[idx.new]][[key]] = hyper.new[[idx.new]][[key]][1L:len]
                            } else if (len == 0L) {
                                hyper.new[[idx.new]][[key]] = numeric(0L)
                            } else {
                                stop("SHOULD NOT HAPPEN")
                            }
                        } else if (len > length(hyper.new[[idx.new]][[key]])) {
                            hyper.new[[idx.new]][[key]] = c(hyper.new[[idx.new]][[key]],
                                                        rep(NA, len - length(hyper.new[[idx.new]][[key]])))
                        } else {
                            stop("SHOULD NOT HAPPEN")
                        }
                    }
                } else {
                    len = length(hyper.new[[idx.new]][[key]])
                }
                ## given length
                llen = length(key.val) - off[[key]]

                if (llen < len) {
                    key.val = c(key.val,  rep(NA,  len - llen))
                }

                if (llen >= 0L) {
                    if (len < 0L) {
                        if (debug) {
                            cat(inla.paste(c("*** Replace hyper.new|", idx.new, "|", key, "|", idxx, " with ", key.val)), "\n")
                        }
                        hyper.new[[idx.new]][[key]] = key.val
                    } else {
                        if (len > 0L) {
                            ## set those NA's to the default ones
                            ii = key.val[off[[key]] + 1L:len]
                            idxx = which(!(is.na(ii) | is.null(ii)))
                            if (length(idxx) > 0L) {
                                if (debug) {
                                    cat(inla.paste(c("*** Replace hyper.new|", idx.new, "|", key, "|", idxx, " with ", ii[idxx])), "\n")
                                }
                                hyper.new[[idx.new]][[key]][idxx] = ii[idxx]
                            }
                        }
                    }
                    off[[key]] = off[[key]] + len
                }
            }

            test.val = (!is.null(h) && !is.null(h[[key]]) && !(is.na(h[[key]])))

            if (debug) {
                cat(paste("*** test.val =", test.val), "\n")
            }

            if (is.na(test.val) || test.val) {
                ## known length
                if (key == "param") {
                    len = inla.model.properties(hyper.new[[idx.new]]$prior, "prior")$nparameters
                    if (debug) {
                        cat(paste("*** nparam =", len), "\n")
                    }
                    if (len < 0L) {
                        ## see explanation above
                        skip.final.check = TRUE
                        hyper.new[[idx.new]][[key]] = h[[key]]
                    } else {
                        ## normal case
                        if (len < length(hyper.new[[idx.new]][[key]]) && len > 0L) {
                            hyper.new[[idx.new]][[key]] = hyper.new[[idx.new]][[key]][1L:len]
                        } else if (len < length(hyper.new[[idx.new]][[key]]) && len == 0L) {
                            hyper.new[[idx.new]][[key]] = NULL
                        } else if (len > length(hyper.new[[idx.new]][[key]])) {
                            hyper.new[[idx.new]][[key]] = c(hyper.new[[idx.new]][[key]],
                                                        rep(NA, len - length(hyper.new[[idx.new]][[key]])))
                        }
                    }
                } else {
                    len = length(hyper.new[[idx.new]][[key]])
                }
                ## given length
                llen = length(h[[key]])

                if (len >= 0L) {
                    if (llen > len) {
                        stop(paste("model",  model, "[", section, "], hyperparam", ih,  ", length(hyper[[key]]) =",
                                   llen,  ">",  len,  sep = " "))
                    } else if (llen < len) {
                        h[[key]] = c(h[[key]],  rep(NA,  len - llen))
                    }
                }
                if (len != 0L) {
                    ## set those NA's to the default ones
                    idxx = which(!is.na(h[[key]]) & !is.null(h[[key]]))
                    if (length(idxx) > 0L) {
                        if (debug) {
                            cat(inla.paste(c("*** Replace hyper.new|", idx.new, "|", key, "|", idxx, " with h|", key, "|",  idxx,
                                             "=",  h[[key]][idxx])), "\n")
                        }
                        hyper.new[[idx.new]][[key]][idxx] = h[[key]][idxx]
                    }
                }
            }

            if (key == "param") {
                ans = inla.model.properties(hyper.new[[idx.new]]$prior, "prior", stop.on.error = TRUE)
                if (ans$nparameters >= 0L) {
                    if (length(hyper.new[[idx.new]]$param) != ans$nparameters) {
                        stop(paste("Wrong length of prior-parameters, prior `", hyper.new[[idx.new]]$prior,  "' needs ",
                                   ans$nparameters,  " parameters, you have ",
                                   length(hyper.new[[idx.new]]$param),  ".",  sep=""))
                    }
                }
            }
        }
    }

    for (key in keywords) {
        key.val = inla.eval(key)
        if (key != "param" || (key == "param" && !skip.final.check)) {
            if (!is.null(key.val) && (length(key.val) > off[[key]])) {
                stop(paste("Length of argument `", key, "' is ", length(key.val),
                           ", does not match the total length of `", key, "' in `hyper' which is ", off[[key]], sep=""))
            }
        }
    }

    return (hyper.new)
}
 
