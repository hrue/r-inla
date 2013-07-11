##! \name{inla.collect.results}
##! \alias{inla.collect.results}
##! \alias{collect.results}
##! \title{Collect results from a inla-call}
##! \description{\code{inla.collect.results} collect results  from a inla-call}
##! \usage{
##! inla.collect.result(
##!              results.dir,
##!              control.results = inla.set.control.results.default(),
##!              debug=FALSE,
##!              only.hyperparam=FALSE,
##!              file.log = NULL)
##!}
##! \arguments{

`inla.collect.results` =
    function(
        ##! \item{results.dir}{The directory where the results of the inla run are stored}
        results.dir, 

        ##! \item{control.results}{a list of parameters controlling the
        ##! output of the function; see \code{?control.results}}
        control.results = inla.set.control.results.default(),

        ##! \item{debug}{Logical. If \code{TRUE} some debugging information are printed}
        debug=FALSE,

        ##! \item{only.hyperparam}{Binary variable indicating wheather only the
        ##! results for the hyperparameters should be collected}
        only.hyperparam=FALSE, 

        ##! \item{file.log}{Character. The filename, if any, of the logfile for
        ##! the internal calculations}
        file.log = NULL)
{
    ##! }
    ##! \value{ The function returns an object of class \code{"inla"}, see the
    ##! help file for \code{inla} for details.}
    ##! 
    ##! \details{This function is mainly used inside \code{inla} and
    ##!  \code{inla.surv} to collect results after running the inla
    ##!  function. It can also be used to collect results into R after having
    ##!  runned a inla section outside R.  }

    if (is.na(file.info(results.dir)$isdir) || 
        !file.info(results.dir)$isdir) {
        stop(paste("This is not a directory: ", results.dir, "\n"))
    }

    filename = paste(results.dir,  "/.ok",  sep="")
    res.ok = file.exists(filename)
    if (!res.ok) {
        inla.inlaprogram.has.crashed()
    }
 
    if (!only.hyperparam) {
        res.fixed = inla.collect.fixed(results.dir, debug)
        res.lincomb = inla.collect.lincomb(results.dir, debug, derived=FALSE)
        res.lincomb.derived = inla.collect.lincomb(results.dir, debug, derived = TRUE)
        res.dic = inla.collect.dic(results.dir, debug)
        res.cpo.pit = inla.collect.cpo(results.dir, debug)
        res.po = inla.collect.po(results.dir, debug)
        res.random = inla.collect.random(results.dir, control.results$return.marginals.random, debug)
        res.predictor = inla.collect.predictor(results.dir, control.results$return.marginals.predictor, debug)
        res.configurations = inla.collect.configurations(results.dir, debug)
        res.spde2.blc = inla.collect.spde2.blc(results.dir, control.results$return.marginals.random, debug)

        file=paste(results.dir,.Platform$file.sep,"neffp",.Platform$file.sep,"neffp.dat", sep="")
        neffp = matrix(inla.read.binary.file(file), 3, 1)
        rownames(neffp) = inla.trim(c("Expectected  number of parameters",
                        "Stdev of the number of parameters",
                        "Number of equivalent replicates"))
    } else {
        res.fixed=NULL
        res.lincomb = NULL
        res.lincomb.derived = NULL
        res.dic=NULL
        res.cpo.pit =NULL
        res.po = NULL
        res.random=NULL
        res.predictor =NULL
        res.configurations = NULL
        res.spde2.blc = NULL
        neffp =NULL
    }
    res.mlik = inla.collect.mlik(results.dir, debug)
    res.q = inla.collect.q(results.dir, debug)
    res.graph = inla.collect.graph(results.dir, debug)
    res.offset = inla.collect.offset.linear.predictor(results.dir, debug)

    ##get the hyperparameters
    theta.mode = inla.read.binary.file(paste(results.dir,.Platform$file.sep,".theta_mode", sep=""))[-1]
    x.mode = inla.read.binary.file(paste(results.dir,.Platform$file.sep,".x_mode", sep=""))[-1]
    hgid = readLines(paste(results.dir,.Platform$file.sep,".hgid", sep=""))
    
    if (length(theta.mode)>0) {
        res.hyper = inla.collect.hyperpar(results.dir, debug)

        ##get the joint (if printed)
        alldir = dir(results.dir)
        if (length(grep("joint.dat", alldir))==1) {
            if (debug) {
                print("inla.collect.joint hyperpar")
            }
            fnm = paste(results.dir,"/joint.dat", sep="")
            if (file.info(fnm)$size > 0) {
                joint.hyper = read.table(fnm)
            } else {
                joint.hyper = NULL
            }
        } else {
            joint.hyper = NULL
        }
    } else {
        res.hyper = NULL
        joint.hyper = NULL
    }

    logfile = inla.collect.logfile(file.log, debug)
    misc = inla.collect.misc(results.dir, debug)
    theta.tags = NULL
    mode.status = NA
    if (!is.null(misc)) {
        ## put also theta.mode in here
        misc$theta.mode = theta.mode
        ## we need theta.tags for later usage
        if (!is.null(misc$theta.tags)) {
            theta.tags = misc$theta.tags
        }
        mode.status = misc$mode.status

        if (!is.null(misc$lincomb.derived.correlation.matrix)) {
            stopifnot(!is.null(res.lincomb.derived))
            id = res.lincomb.derived$summary.lincomb.derived$ID
            tag = rownames(res.lincomb.derived$summary.lincomb.derived)
            R = misc$lincomb.derived.correlation.matrix
            rownames(R) = colnames(R) = tag[id]
            misc$lincomb.derived.correlation.matrix = R
        }
        if (!is.null(misc$lincomb.derived.covariance.matrix)) {
            ## same code as above
            stopifnot(!is.null(res.lincomb.derived))
            id = res.lincomb.derived$summary.lincomb.derived$ID
            tag = rownames(res.lincomb.derived$summary.lincomb.derived)
            R = misc$lincomb.derived.covariance.matrix
            rownames(R) = colnames(R) = tag[id]
            misc$lincomb.derived.covariance.matrix = R
        }
    }

    ## add the names of the theta's here, as they are available.
    if (!is.null(misc) && !is.null(joint.hyper)) {
        colnames(joint.hyper) = c(misc$theta.tags, "Log posterior density")
    }

    res = c(res.fixed, res.lincomb, res.lincomb.derived, res.mlik,
            list(cpo=res.cpo.pit), list(po = res.po), 
            res.random, res.predictor, res.hyper,
            res.configurations, res.offset, res.spde2.blc, logfile, 
            list(misc = misc,
                 dic=res.dic, mode = list(theta = theta.mode, x = x.mode,
                                      theta.tags = theta.tags, mode.status = mode.status,
                                      log.posterior.mode = misc$log.posterior.mode),
                 neffp=neffp,
                 joint.hyper=joint.hyper, nhyper=length(theta.mode),
                 version = list(inla.call = hgid, inla.call.builtin = hgid, R.INLA=inla.version("hgid"))), 
            list(Q=res.q),
            res.graph, ok = res.ok)
    class(res) = "inla"

    if (inla.getOption("internal.experimental.mode")) {
        if (debug)
            print("...Fix marginals")
        
        ## set the inla.marginal class to all the marginals, and add tag
        ## used for plotting.  all these have two levels:
        idxs = grep("marginals[.](fixed|linear[.]predictor|lincomb[.]derived|lincomb|hyperpar|fitted[.]values)", names(res))
        if (length(idxs) > 0) {
            for(idx in idxs) {
                if (!is.null(res[[idx]])) {
                    name.1 = names(res)[idx]
                    attr(res[[idx]], "inla.tag") = name.1
                    class(res[[idx]]) = "inla.marginals"

                    if (length(res[[idx]])>0) {
                        for(i in 1:length(res[[idx]])) {
                            name.2 = names(res[[idx]])[i]
                            if (!is.null(res[[idx]][[i]])) {
                                attr(res[[idx]][[i]], "inla.tag") = paste(name.1, name.2)
                                class(res[[idx]][[i]]) = "inla.marginal"
                            }
                        }
                    }
                }
            }
        }

        if (debug)
            print("...Fix marginals 1")

        ## all these have three levels:
        idxs = grep("marginals[.]random", names(res))
        if (length(idxs) > 0) {
            for(idx in idxs) {
                if (!is.null(res[[idx]])) {
                    name.1 = names(res)[idx]
                    name.2 = names(res[[idx]])

                    if (length(res[[idx]])>0) {
                        for(i in 1:length(res[[idx]])) {
                            name.3 = name.2[i]
                            name.4 = names(res[[idx]][[i]])

                            attr(res[[idx]][[i]], "inla.tag") = paste(name.1, name.3)
                            class(res[[idx]][[i]]) = "inla.marginals"

                            if (length(res[[idx]][[i]]) > 0) {
                                for(j in 1:length(res[[idx]][[i]])) {
                                    name.5 = name.4[j]
                                    if (!is.null(res[[idx]][[i]][[j]])) {
                                        attr(res[[idx]][[i]][[j]], "inla.tag") = paste(name.1, name.3, name.5)
                                        class(res[[idx]][[i]][[j]]) = "inla.marginal"
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        if (debug)
            print("...Fix marginals done.")
    }

    return(res)
}
