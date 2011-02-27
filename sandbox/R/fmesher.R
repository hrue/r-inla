mesh.segm = function(loc=NULL,idx=NULL,grp=NULL,is.bnd=TRUE)
{
    if ((missing(loc) || is.null(loc)) &&
        (missing(idx) || is.null(idx)))
        stop("At most one of 'loc' and 'idx' may be missing or null.")
    if (!missing(loc) && !is.null(loc)) {
        if (!is.matrix(loc))
            stop("'loc' must be a a matrix")
        if (missing(idx) || is.null(idx))
            idx = (inla.ifelse(is.bnd,
                               c(1:nrow(loc),1),
                               c(1:nrow(loc))))
    } else {
        loc = NULL
    }

    if (!missing(idx) && !is.null(idx)) {
        if (!is.vector(idx) && !is.matrix(idx))
            stop("'idx' must be a vector or a matrix")
        if (is.vector(idx))
            idx = as.matrix(idx,nrow=length(idx),ncol=1)
        if (ncol(idx)==1) {
            if (nrow(idx)<2)
                warning("Segment specification must have at least 2 idxices.")
            idx = matrix(c(idx[-nrow(idx)],idx[-1]), nrow(idx)-1, 2)
        }
        storage.mode(idx) <- "integer"
        if (!is.null(loc) && (max(idx, na.rm=TRUE)>nrow(loc)))
            warning("Segment indices (max=", max(idx, na.rm=TRUE),
                    ") exceed specified location list length (",
                    nrow(loc), ").")
    }

    if (!missing(grp) && !is.null(grp)) {
        if (!is.vector(grp) && !is.matrix(grp))
            stop("'grp' must be a vector or a matrix")
        grp = matrix(grp, min(length(grp), nrow(idx)), 1)
        if (nrow(grp)<nrow(idx))
            grp = (matrix(c(as.vector(grp),
                            rep(grp[nrow(grp)], nrow(idx)-length(grp))),
                            nrow(idx), 1))
        storage.mode(grp) <- "integer"
    } else
        grp = NULL

    ## Filter away NAs in loc and idx
    if (!is.null(loc)) {
        idx[is.na(idx)] = 0L ## Avoid R annoyances with logical+NA indexing
        while (sum(is.na(loc))>0) {
            i = min(which(rowSums(is.na(loc))>0))
            loc = loc[-i,,drop=FALSE]
            idx[idx==i] = 0L
            idx[idx>i] = idx[idx>i]-1L
        }
        idx[idx==0L] = NA
    }
    while (sum(is.na(idx))>0) {
        i = min(which(rowSums(is.na(idx))>0))
        idx = idx[-i,,drop=FALSE]
        if (!is.null(grp))
            grp = grp[-i,,drop=FALSE]
    }

    ret = list(loc=loc, idx=idx, grp=grp, is.bnd=is.bnd)
    class(ret) <- "fmesher.segm"
    return(ret)
}


lines.fmesher.segm = function (segm, loc=NULL, ...)
{
    if (!is.null(segm$loc))
        loc = segm$loc
    stopifnot(!is.null(loc), ncol(loc)>=2)

    for (grp in unique(segm$grp)) {
        idx = which(segm$grp==grp)
        lines(loc[t(cbind(segm$idx[idx,, drop=FALSE], NA)), 1],
              loc[t(cbind(segm$idx[idx,, drop=FALSE], NA)), 2],
              type="l",
              col=c("black", "blue", "red", "green")[1+(grp%%4)],
              ...)
    }
}

plot.fmesher.mesh = function (mesh, add=FALSE, lwd=1, ...)
{
    if (!inherits(mesh, "fmesher.mesh"))
        stop("'mesh' bust be an 'fmesher.mesh' object.")

    idx = cbind(mesh$graph$tv[,c(1:3,1), drop=FALSE], NA)
    x = mesh$loc[t(idx), 1]
    y = mesh$loc[t(idx), 2]

    if (!add) {
        plot.new()
        plot.window(xlim=range(mesh$loc[,1]), ylim=range(mesh$loc[,2]), "")
    }
    lines(x, y, type="l")
    if (!is.null(mesh$segm$bnd))
        lines(mesh$segm$bnd, mesh$loc, lwd=lwd+1, ...)
    if (!is.null(mesh$segm$int))
        lines(mesh$segm$int, mesh$loc, lwd=lwd+1, ...)
    if (!add)
        title("Constrained refined Delaunay triangulation",
              deparse(substitute(mesh)))
    return(invisible())
}


############################
##    parse.grid.input = function(grid)
##    {
##        default = (list(dims=c(10,10),
##                        xlim=c(0,1),
##                        ylim=c(0,1),
##                        units=""))
##        for (name in names(grid))
##            default[name] = grid[name]
##        return(default)
##    }

mesh.grid = function(dims=c(2,2), xlim=c(0,1), ylim=c(0,1), x=NULL, y=NULL)
{
    if (is.null(x)) {
        x = seq(xlim[1], xlim[2], length.out=dims[1])
    } else {
        x = as.vector(x)
        dims[1] = length(x)
        xlim = range(x)
    }
    if (is.null(y)) {
        y = seq(ylim[1], ylim[2], length.out=dims[2])
    } else {
        y = as.vector(y)
        dims[2] = length(y)
        ylim = range(x)
    }
    loc = (cbind(rep(x, times = dims[2]),
                 rep(y, each = dims[1])))
    ## Construct grid boundary
    segm.idx = (c(1:(dims[1]-1),
                  dims[1]*(1:(dims[2]-1)),
                  dims[1]*dims[2]-(0:(dims[1]-2)),
                  dims[1]*((dims[2]-1):1)+1))
    segm = mesh.segm(loc[segm.idx,, drop=FALSE], is.bnd=TRUE)

    grid = list(loc=loc, segm=segm)
    class(grid) = "fmesher.grid"
    return(grid)
}




inla.mesh.parse.segm.input = function(boundary=NULL, interior=NULL, n=0)
{
###########################################
    homogenise.segm.input = function(x, is.bnd)
    {
        if (is.matrix(x) || is.vector(x)) { ## Coordinates or indices
            x = (inla.ifelse(is.matrix(x),x,
                             as.matrix(x,nrow=length(x),ncol=1)))
            if (is.integer(x)) { ## Indices
                ret = mesh.segm(NULL,x,NULL,is.bnd)
            } else if (is.numeric(x)) { ## Coordinates
                ret = mesh.segm(x,NULL,NULL,is.bnd)
            } else {
                stop("Segment info matrix must be numeric or integer.")
            }
        } else if (inherits(x, "fmesher.segm")) {
            ## Override x$is.bnd
            ret = mesh.segm(x$loc,x$idx,x$grp,is.bnd)
        } else if (!is.null(x)) {
            stop("Segment info must be matrix or 'fmesher.segm' object.")
        } else {
            ret = NULL
        }
        return(ret)
    }
##################################################
    homogenise.segm.grp = function(input) {
        grp.idx = 0L
        for (k in 1:length(input)) if (!is.null(input[[k]])) {
            if (!inherits(input[[k]], "fmesher.segm")) {
                stop("Segment info is not a segment list. ",
                     "This should not happen.")}
            if (is.null(input[[k]]$grp)) {
                grp.idx = grp.idx+1L
                input[[k]] = (mesh.segm(input[[k]][[1]],
                                        input[[k]][[2]],
                                        grp.idx,
                                        input[[k]][[4]]))
            } else {
                grp.idx = max(grp.idx, input[[k]]$grp, na.rm=TRUE)
            }
        }
        return(input)
    }
##################################################
    parse.segm.input = function(input, loc.n) {
        loc = NULL
        bnd = list(idx = matrix(,0,2), grp = matrix(,0,1))
        int = list(idx = matrix(,0,2), grp = matrix(,0,1))
        storage.mode(bnd$idx) <- "integer"
        storage.mode(bnd$grp) <- "integer"
        storage.mode(int$idx) <- "integer"
        storage.mode(int$grp) <- "integer"
        for (k in 1:length(input)) if (!is.null(input[[k]])) {
            if (!inherits(input[[k]],"fmesher.segm")) {
                stop("Segment info is not a segment list. ",
                     "This should not happen.")
            }
            prev.loc.n = loc.n
            if (!is.null(input[[k]]$loc)) {
                extra.loc.n = nrow(input[[k]]$loc)
                idx.offset = loc.n
                loc.n = loc.n + extra.loc.n
                loc = (inla.ifelse(is.null(loc),
                                   input[[k]]$loc,
                                   rbind(loc, input[[k]]$loc)))
            } else {
                idx.offset = 0L
            }
            if (input[[k]]$is.bnd) {
                bnd$idx = rbind(bnd$idx, input[[k]]$idx+idx.offset)
                bnd$grp = rbind(bnd$grp, input[[k]]$grp)
            } else {
                int$idx = rbind(int$idx, input[[k]]$idx+idx.offset)
                int$grp = rbind(int$grp, input[[k]]$grp)
            }
        }
        if (nrow(bnd$idx)==0)
            bnd = NULL
        if (nrow(int$idx)==0)
            int = NULL
        return(list(loc = loc,
                    bnd = (inla.ifelse(is.null(bnd),
                                       NULL,
                                       mesh.segm(bnd$loc,
                                                 bnd$idx,
                                                 bnd$grp,
                                                 TRUE))),
                    int = (inla.ifelse(is.null(int),
                                       NULL,
                                       mesh.segm(int$loc,
                                                 int$idx,
                                                 int$grp,
                                                 FALSE)))))
    }
###########################

    segm = (c(lapply(inla.ifelse(inherits(boundary, "list"),
                                 boundary, list(boundary)),
                     function(x){homogenise.segm.input(x, TRUE)}),
              lapply(inla.ifelse(inherits(interior, "list"),
                                 interior, list(interior)),
                     function(x){homogenise.segm.input(x, FALSE)})))
    segm = homogenise.segm.grp(segm)
    return(parse.segm.input(segm, n))
}



`inla.mesh` = function(loc=NULL, tv=NULL,
                       boundary=NULL, interior=NULL,
                       extend=(missing(tv) || is.null(tv)),
                       refine=FALSE,
                       grid=NULL, manifold=NULL,
                       cutoff.distance = 0,
                       plot.delay = NULL,
                       dir = tempdir(), prefix = NULL)
{

################################
    filter.locations = function(loc, cutoff.distance)
    {
        ## Map locations to nodes, avoiding near-duplicates.
        loc.n = nrow(loc)
        loc.dim = ncol(loc)

        loc.is.na = (rowSums(is.na(loc))>0)
        if (sum(loc.is.na)>0)
            stop("NAs in locations not yet supported.")

        node.coord = matrix(nrow=loc.n, ncol=loc.dim)
        map.loc.to.node = rep(0L, nrow=loc.n)
        excluded = c()
        loc.i = 1L
        node.i.max = 1L
        node.coord[node.i.max,] = loc[loc.i,]
        map.loc.to.node[[loc.i]] = node.i.max
        for (loc.i in 2:loc.n) {
            loc.to.node.dist =
                sqrt(rowSums((as.matrix(rep(1, node.i.max)) %*%
                              loc[loc.i,] -
                              node.coord[1:node.i.max,])^2))
            if (min(loc.to.node.dist) > cutoff.distance) {
                node.i.max = node.i.max+1L
                node.coord[node.i.max,] = loc[loc.i,]
                map.loc.to.node[[loc.i]] = node.i.max
            } else {
                excluded = c(excluded, loc.i)
            }
        }
        ## Remove excess nodes.
        node.coord = node.coord[1:node.i.max,]

        ## Identify nearest nodes for excluded locations.
        for (loc.i in excluded) {
            loc.to.node.dist =
                sqrt(rowSums((as.matrix(rep(1, node.i.max)) %*%
                              loc[loc.i,] -
                              node.coord)^2))
            node.i = which.min(loc.to.node.dist);
            map.loc.to.node[loc.i] = node.i
        }

        return(list(loc = node.coord, node.idx = map.loc.to.node))
    }
############################

    if (!missing(manifold))
        warning("Option 'manifold' not implemented.")
    if (is.logical(extend) && extend) extend = list()
    if (is.logical(refine) && refine) refine = list()

    if (missing(grid) || is.null(grid)) {
        grid = list(loc=NULL, segm=NULL)
        grid.n = 0
    } else {
        if (!inherits(grid, "fmesher.grid"))
            stop("'grid' must be an 'fmesher.grid' object.")
        if (!is.null(tv)) {
            warning("Both 'grid' and 'tv' specified.  Ignoring 'tv'.")
            tv = NULL
        }
        if (!inherits(extend, "list")) {
            boundary = (c(inla.ifelse(inherits(boundary, "list"),
                                      boundary, list(boundary)),
                          list(grid$segm)))
        }
        grid.n = max(0,nrow(grid$loc))
    }
    loc.n = max(0,nrow(loc))

    segm = inla.mesh.parse.segm.input(boundary, interior, loc.n+grid.n)
    segm.n = max(0,nrow(segm$loc))

    tmp = filter.locations(rbind(loc, grid$loc, segm$loc), cutoff.distance)
    loc0 = tmp$loc
    idx.all = tmp$node.idx
    idx = (list(loc = inla.ifelse(loc.n>0, idx.all[1:loc.n], NULL),
                grid = inla.ifelse(grid.n>0, idx.all[loc.n+(1:grid.n)], NULL),
                segm = (inla.ifelse(segm.n>0,
                                    idx.all[loc.n+grid.n+(1:segm.n)],
                                    NULL))))

    ## Remap indices
    if (!is.null(tv)) {
        tv = matrix(idx.all[as.vector(tv)], nrow=nrow(tv), ncol=ncol(tv))
    }
    if (!is.null(segm$bnd)) {
        segm$bnd$idx = (matrix(idx.all[as.vector(segm$bnd$idx)],
                               nrow=nrow(segm$bnd$idx),
                               ncol=ncol(segm$bnd$idx)))
    }
    if (!is.null(segm$int)) {
        segm$int$idx = (matrix(idx.all[as.vector(segm$int$idx)],
                               nrow=nrow(segm$int$idx),
                               ncol=ncol(segm$int$idx)))
    }

    ## Where to put the files?
    prefix = inla.fmesher.make.prefix(dir, prefix)

    ## The default name for the input locations in fmesher is s0:
    loc.file = inla.write.fmesher.file(loc0,
      filename = paste(prefix, "input.s", sep=""))
    all.args = "--input=input.s"

    if (!is.null(tv)) {
        inla.write.fmesher.file(tv-1L, filename = paste(prefix, "input.tv", sep=""))
        all.args = paste(all.args, ",input.tv", sep="")
    }
    if (!is.null(segm$bnd)) {
        inla.write.fmesher.file(segm$bnd$idx-1L,
                                filename = paste(prefix, "input.segm.bnd.idx", sep=""))
        inla.write.fmesher.file(segm$bnd$grp,
                                filename = paste(prefix, "input.segm.bnd.grp", sep=""))
        all.args = paste(all.args," --boundary=input.segm.bnd.idx")
        all.args = paste(all.args," --boundarygrp=input.segm.bnd.grp")
    }
    if (!is.null(segm$int)) {
        inla.write.fmesher.file(segm$int$idx-1L,
                                filename = paste(prefix, "input.segm.int.idx", sep=""))
        inla.write.fmesher.file(segm$int$grp,
                                filename = paste(prefix, "input.segm.int.grp", sep=""))
        all.args = paste(all.args," --interior=input.segm.int.idx")
        all.args = paste(all.args," --interiorgrp=input.segm.int.grp")
    }


    if (inherits(extend,"list")) {
        cet = c(0,0)
        cet[1] = inla.ifelse(is.null(extend$n), 8, extend$n)
        cet[2] = inla.ifelse(is.null(extend$offset), -0.15, extend$offset)
        all.args = (paste(all.args," --cet=",
                          cet[1],",", cet[2], sep=""))
    }
    if (inherits(refine,"list")) {
        rcdt = c(0,0,0)
        rcdt[1] = inla.ifelse(is.null(refine$min.angle), 21, refine$min.angle)
        rcdt[2] = inla.ifelse(is.null(refine$max.edge), -1.0, refine$max.edge)
        rcdt[3] = inla.ifelse(is.null(refine$max.edge), -0.5, refine$max.edge)
        rcdt[2] = (inla.ifelse(is.null(refine$max.edge.extra),
                               rcdt[2], refine$max.edge.extra))
        rcdt[3] = (inla.ifelse(is.null(refine$max.edge.data),
                               rcdt[3], refine$max.edge.data))
        all.args = (paste(all.args," --rcdt=",
                          rcdt[1],",", rcdt[2],",", rcdt[3], sep=""))
    }

    if (!is.null(plot.delay)) {
        all.args = paste(all.args," --x11=", plot.delay, sep="")
    }
    all.args = paste(all.args, inla.getOption("fmesher.arg"))

    ## Call fmesher:
    time.used = system.time({
        echoc = (inla.fmesher.call(all.args=all.args,
                                   prefix=prefix))})

    ## Read the mesh:
    loc = inla.read.fmesher.file(paste(prefix, "s", sep=""))
    graph = (list(tv = 1L+inla.read.fmesher.file(paste(prefix, "tv", sep="")),
                  vt = 1L+inla.read.fmesher.file(paste(prefix, "vt", sep="")),
                  tt = 1L+inla.read.fmesher.file(paste(prefix, "tt", sep="")),
                  tti = 1L+inla.read.fmesher.file(paste(prefix, "tti", sep="")),
                  vv = inla.read.fmesher.file(paste(prefix, "vv", sep=""))))
    graph$tv[graph$tv==0L] = NA
    graph$vt[graph$vt==0L] = NA
    graph$tt[graph$tt==0L] = NA
    graph$tti[graph$tti==0L] = NA

    ## Read constraint segment information:
    segm.bnd = (mesh.segm(NULL,
                          1L+inla.read.fmesher.file(paste(prefix,
                                                          "segm.bnd",
                                                          sep="")),
                          inla.read.fmesher.file(paste(prefix,
                                                       "segm.bnd.grp",
                                                       sep="")),
                          TRUE))
    if (!is.null(segm$int)) {
        segm.int = (mesh.segm(NULL,
                              1L+inla.read.fmesher.file(paste(prefix,
                                                              "segm.int",
                                                              sep="")),
                              inla.read.fmesher.file(paste(prefix,
                                                           "segm.int.grp",
                                                           sep="")),
                              FALSE))
    } else {
        segm.int = NULL
    }

    mesh = (list(meta = (list(call=match.call(),
                              fmesher.args = all.args,
                              time = time.used,
                              prefix = prefix)),
                 loc = loc,
                 graph = graph,
                 segm = list(bnd=segm.bnd, int=segm.int),
                 idx = idx))

    class(mesh) <- "fmesher.mesh"

    return(mesh)
}

`summary.fmesher.mesh` = function(x, verbose=FALSE, ...)
{
    ## provides a summary for a mesh object
    if (!inherits(x, "fmesher.mesh"))
        stop("'x' must inherit from class \"fmesher.mesh\"")

    ret = list(verbose=verbose)
    if (verbose) {
        ret = c(ret, list(call=x$meta$call))
        ret = c(ret, list(fmesher.args=x$meta$fmesher.args))
    }
    ret = c(ret, list(nV=nrow(mesh$loc)))
    ret = c(ret, list(nT=nrow(mesh$graph$tv)))
    ret = c(ret, list(xlim=range(mesh$loc[,1])))
    ret = c(ret, list(ylim=range(mesh$loc[,2])))
    ret = c(ret, list(zlim=range(mesh$loc[,3])))

    my.segm = function(x) {
        if (is.null(x))
            return(list(n=0, grps=NULL))
        n = max(0, nrow(x$idx))
        if (max(0, length(unique(x$grp)))>0) {
            grps = unique(x$grp)
        } else {
            grps = NULL
        }
        return(list(n=n, grps=grps))
    }
    if(!is.null(x$segm)) {
        ret = c(ret, list(segm.bnd=my.segm(x$segm$bnd)))
        ret = c(ret, list(segm.int=my.segm(x$segm$int)))
    } else {
        ret = c(ret, list(segm.bnd=my.segm(NULL)))
        ret = c(ret, list(segm.int=my.segm(NULL)))
    }

    class(ret) <- "summary.fmesher.mesh"
    return (ret)
}

`print.summary.fmesher.mesh` = function(x, ...)
{
    if (!inherits(x, "summary.fmesher.mesh"))
        stop("'x' must inherit from class \"summary.fmesher.mesh\"")

    if (x$verbose) {
        cat("\nCall:\n")
        print(x$call)

        cat("\nfmesher arguments:\n", x$fmesher.args, "\n", sep = "")
    }

    cat("\nVertices:\t", as.character(x$nV), "\n", sep="")
    cat("Triangles:\t", as.character(x$nT), "\n", sep="")

    my.print.segm = function(x) {
        cat(as.character(x$n))
        if (!is.null(x$grps)) {
            n = length(x$grps)
            cat(" (", n, " group", inla.ifelse(n==1, "", "s"), sep="")
            if (n <= 10) {
                cat(":", x$grps, sep=" ")
            } else {
                cat(":", x$grps[1:10], "...", sep=" ")
            }
            cat(")")
        }
        cat("\n", sep="")
        return(invisible())
    }
    cat("Boundary segm.:\t")
    my.print.segm(x$segm.bnd)
    cat("Interior segm.:\t")
    my.print.segm(x$segm.int)
    cat("xlim:\t", x$xlim[1], " ", x$xlim[2], "\n", sep="")
    cat("ylim:\t", x$ylim[1], " ", x$ylim[2], "\n", sep="")
    cat("zlim:\t", x$zlim[1], " ", x$zlim[2], "\n", sep="")
    cat("\n")

}
