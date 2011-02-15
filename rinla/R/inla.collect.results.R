
`inla.collect.results` =
    function(results.dir, control.results = inla.set.control.results.default(),
             debug=FALSE, only.hyperparam=FALSE)
{
    if (is.na(file.info(results.dir)$isdir) || 
        !file.info(results.dir)$isdir) {
        stop(paste("This is not a directory: ", results.dir, "\n"))
    }

    if (!only.hyperparam) {
        res.fixed = inla.collect.fixed(results.dir, debug)
        res.lincomb = inla.collect.lincomb(results.dir, debug)
        res.dic = inla.collect.dic(results.dir, debug)
        res.cpo.pit = inla.collect.cpo(results.dir, debug)
        res.random = inla.collect.random(results.dir, control.results$return.marginals.random, debug)
        res.predictor = inla.collect.predictor(results.dir, control.results$return.marginals.predictor, debug)
        res.configurations = inla.collect.configurations(results.dir, debug)

        file=paste(results.dir,.Platform$file.sep,"neffp",.Platform$file.sep,"neffp.dat", sep="")
        neffp = matrix(inla.read.binary.file(file), 3, 1)
        rownames(neffp) = inla.trim(c("Expectected  number of parameters",
                        "Stdev of the number of parameters",
                        "Number of equivalent replicates"))

    } else {
        res.fixed=NULL
        res.lincomb = NULL
        res.dic=NULL
        res.cpo.pit =NULL
        res.random=NULL
        res.predictor =NULL
        res.configurations = NULL
        neffp =NULL
    }
    res.mlik = inla.collect.mlik(results.dir, debug)
    res.q = inla.collect.q(results.dir, control.results$image.dim, debug)
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

    misc = inla.collect.misc(results.dir, debug)
    theta.tags = NULL
    if (!is.null(misc)) {
        ## put also theta.mode in here
        misc$theta.mode = theta.mode
        ## we need theta.tags for later usage
        if (!is.null(misc$theta.tags)) {
            theta.tags = misc$theta.tags
        }
    }

    res = c(res.fixed, res.lincomb, res.mlik, res.cpo.pit, res.random, res.predictor, res.hyper,
            res.configurations, res.offset,
            list(misc = misc,
                 dic=res.dic, mode = list(theta = theta.mode, x=x.mode,
                                      theta.tags = theta.tags),
                 neffp=neffp,
                 joint.hyper=joint.hyper, nhyper=length(theta.mode),
                 version = list(inla = hgid, Rinla=inla.version(hgid=TRUE))),
            res.q)
    class(res) = "inla"

    return(res)
}
