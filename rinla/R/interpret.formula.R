## nothing to export

`inla.interpret.formula` = function (
    gf,
    debug=FALSE,
    data.same.len = NULL,
    data=NULL,
    data.model = NULL,
    parent.frame = NULL)
{
    ## use a copy if data.model is !NULL, as we have to assign an
    ## entry to it. otherwise, just use the
    if (!is.null(data.model)) {
        ## make a copy of the environment
        p.env = new.env(hash = TRUE, parent = environment(gf))
        ## then add the entries in data.model to that environment
        for(nm in names(data.model)) {
            idx = which(names(data.model) == nm)
            assign(nm, data.model[[idx]], envir = p.env)
        }
    } else {
        ## otherwise, just use (for read-only) this environment
        if (missing(parent.frame)) {
            p.env = environment(gf)
        } else {
            p.env = parent.frame
        }
    }

    ## argument data is here only used for infereing '.'
    tf = terms.formula(gf, specials = c("f"), data=NULL)
    terms = attr(tf, "term.labels")
    nt = length(terms)

    if (attr(tf, "response") > 0) {
        ## fixf formula with ONLY fixed effects.  randf formula with
        ## ONLY random effect.  weightf formula where are the names of
        ## the (possible) weigths for the covariates
        response = as.character(attr(tf, "variables")[2])
        fixf = randf = weightf= paste(response, "~", sep = "")
    } else {
        stop("\n\tA response variable has to be present")
    }

    if (length(attr(tf, "offset")) > 0) {
        i = grep("^offset[(]", attr(tf, "variables"))
        offset = as.character(attr(tf, "variables")[i])
        if (debug)
            cat("found offset\n")
    } else {
        offset = NULL
        if (debug)
            cat("no offset\n")
    }

    rt = attr(tf, "specials")$f
    vtab = attr(tf, "factors")
    if (length(rt) > 0) {
        for (i in 1:length(rt)) {
            ind = (1:nt)[as.logical(vtab[rt[i], ])]
            rt[i] = ind
        }
    }

    k =  ks = kp = 1
    len.rt = length(rt)
    random.spec = list()
    if (nt>0) {
        for (i in 1:nt) {
            if (k <= len.rt && ((ks <= len.rt && rt[ks] == i))) {
                st = eval(parse(text = gsub("^f\\(","INLA::f(", terms[i])), envir = data, enclos = p.env)
                random.spec[[k]] = st
                if (ks <= len.rt && rt[ks] == i) {
                    ks = ks + 1
                } else {
                    kt = kt + 1
                }
                k = k + 1
            } else {
                if (kp > 1) {
                    fixf = paste(fixf, " + ", terms[i], sep = "")
                } else {
                    fixf = paste(fixf, terms[i], sep = "")
                }
                kp = kp + 1
            }
        }
    }

    ##if -1 the intercept is not included
    if (attr(tf, "intercept") == 0) {
        intercept=FALSE
    } else {
        intercept=TRUE
    }

    if (length(random.spec) > 0) {
        ## check that we are not using the same covariate more than
        ## one time.
        all.terms = c()
        for (i in 1:length(random.spec)) {
            all.terms = c(all.terms, random.spec[[i]]$term)
        }
        for (i in 1:length(random.spec)) {
            if (sum(random.spec[[i]]$term == all.terms) > 1) {
                ## duplicate names!
                stop(inla.paste(c("\n",
                           "The covariate `", random.spec[[i]]$term, "' are used in the following f()-terms: ",
                                  inla.paste(which(random.spec[[i]]$term == all.terms), sep=", "),
                           "\n",
                           "Only one is allowed since the naming convention of the f()-terms depends on unique names of the covariates",
                           "\n\n",
                           "You can resolve this as follows:",
                           "\n",
                           "Not allowed: formula = y ~ f(x, <etc>) + f(x, <etc>) + <etc>",
                           "\n",
                           "Fix: data$xx = data$x; formula = y ~ f(x, <etc>) + f(xx, <etc>) + <etc>",
                           "\n\n"), sep=""))
            }
        }
    }

    ##number of covariate which have weights attached!
    n.weights = 0
    if (length(random.spec) > 0)
        for (i in 1:length(random.spec)) {
            ff1 = paste(random.spec[[i]]$term, collapse = "+")

            if(!is.null(random.spec[[i]]$weights)) {
                ww1 = paste(random.spec[[i]]$weights, collapse = "+")
                n.weights = n.weights+1
            }
            if(i==1) {
                randf = paste(randf, ff1)
                if(!is.null(random.spec[[i]]$weights))
                    weightf = paste(weightf, ww1)
            } else {
                randf = paste(randf, "+", ff1)
                if(!is.null(random.spec[[i]]$weights))
                    weightf = paste(weightf,"+", ww1)
            }
        }

    if(len.rt>0) {
        randf = as.formula(randf, p.env)
    } else {
        randf = NULL
    }

    if (intercept) {
        fixf = paste(fixf, " +1", sep="")
        nt = nt + 1
    } else {
        fixf = paste(fixf, " -1", sep="")
        nt = nt + 1
    }

    if ((nt-len.rt)>0) {
        fixf = as.formula(fixf, p.env)
    } else {
        fixf = NULL
    }

    if(n.weights> 0) {
        weightf = as.formula(weightf, p.env)
    } else {
        weightf = NULL
    }

    ret = list(randf=randf,
        random.spec = random.spec, n.random = len.rt,
        fixf=fixf, n.fix = (nt-len.rt),
        weightf=weightf, n.weights=n.weights,
        offset=offset, response = response)

    return (ret)
}
