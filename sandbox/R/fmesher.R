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
    class(ret) <- c("list", "inla.fmesher.mesh.segm")
    return(ret)
}


plot.inla.fmesher.mesh.segm = function (segm, loc=NULL)
{
    if (!is.null(segm$loc))
        loc = segm$loc
    stopifnot(!is.null(loc), ncol(loc)>=2)

    for (grp in unique(segm$grp)) {
        idx = which(segm$grp==grp)
        lines(loc[t(cbind(segm$idx[idx,,drop=FALSE],NA)),1],
              loc[t(cbind(segm$idx[idx,,drop=FALSE],NA)),2],
              type="l",
              col=c("black","blue","red","green")[1+(grp%%4)],lwd=5)
    }
}


`inla.mesh` = function(s,tv=NULL,
                       boundary=NULL, interior=NULL,
                       extend=FALSE, refine=FALSE,
                       grid=NULL, manifold=NULL,
                       cutoff.distance = 0,
                       plot.delay = NULL,
                       dir = tempdir(), prefix = NULL)
{

##################################################
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
                stop("Unexpected input.  Segment info matrix must be numeric or integer.")
            }
        } else if (is.list(x) || inherits(x, "inla.fmesher.mesh.segm")) {
            if (length(x)==1) {
                ret = mesh.segm(x[[1]],NULL,NULL,is.bnd)
            } else if (length(x)==2) {
                ret = mesh.segm(x[[1]],x[[2]],NULL,is.bnd)
            } else if (length(x)==3) {
                ret = mesh.segm(x[[1]],x[[2]],x[[3]],is.bnd)
            } else if (length(x)==4) {
                if (x[[4]] != is.bnd)
                    warning("Boundary/interior mismatch in constraint input.")
                ret = mesh.segm(x[[1]],x[[2]],x[[3]],is.bnd)
            } else {
                stop("Unexpected input.  Segment info must have length 1, 2, 3, or 4.")
            }
        } else if (!is.null(x)) {
            stop("Unexpected input.  Segment info must be matrix or list.")
        } else {
            ret = NULL
        }
        return(ret)
    }
##################################################
    homogenise.segm.grp = function(input) {
        grp.idx = 0L
        for (k in 1:length(input)) if (!is.null(input[[k]])) {
            if (!inherits(input[[k]], "inla.fmesher.mesh.segm")) {
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
            if (!inherits(input[[k]],"inla.fmesher.mesh.segm")) {
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
        return(list(loc = loc, bnd = bnd, int = int))
    }
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

    if (!missing(grid))
        warning("Option 'grid' not implemented.")
    if (!missing(manifold))
        warning("Option 'manifold' not implemented.")

    segm = (c(lapply(inla.ifelse(inherits(boundary, "list") &&
                                 !inherits(boundary,
                                           "inla.fmesher.mesh.segm"),
                                 boundary, list(boundary)),
                     function(x){homogenise.segm.input(x, TRUE)}),
              lapply(inla.ifelse(inherits(interior, "list") &&
                                 !inherits(interior,
                                           "inla.fmesher.mesh.segm"),
                                 interior, list(interior)),
                     function(x){homogenise.segm.input(x, FALSE)})))
    segm = homogenise.segm.grp(segm)
    segm = parse.segm.input(segm, nrow(s))

    tmp = filter.locations(rbind(s, segm$loc), cutoff.distance)
    loc0 = tmp$loc
    idx = tmp$node.idx

    ## Remap indices
    if (!is.null(tv)) {
        tv = matrix(idx[as.vector(tv)], nrow=nrow(tv), ncol=ncol(tv))
    }
    if (!is.null(segm$bnd)) {
        segm$bnd$idx = (matrix(idx[as.vector(segm$bnd$idx)],
                               nrow=nrow(segm$bnd$idx),
                               ncol=ncol(segm$bnd$idx)))
    }
    if (!is.null(segm$int)) {
        segm$int$idx = (matrix(idx[as.vector(segm$int$idx)],
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


##    if (!inherits(cet)) {
##        all.args = paste(all.args," --cet=8,-0.15", sep="")
##        all.args = paste(all.args," --cet=", cet[1],",", cet[2], sep="")
##    }
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
                               rcdt[2], refine$max.edge.data))
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
    mesh = (list(s = inla.read.fmesher.file(paste(prefix, "s", sep="")),
                 tv = 1L+inla.read.fmesher.file(paste(prefix, "tv", sep=""))))

    ## Read constraint segment information:
    segm.bnd = (mesh.segm(NULL,
                          1L+inla.read.fmesher.file(paste(prefix, "segm.bnd", sep="")),
                          inla.read.fmesher.file(paste(prefix, "segm.bnd.grp", sep="")),
                          TRUE))
    if (!is.null(segm$int)) {
        segm.int = (mesh.segm(NULL,
                              1L+inla.read.fmesher.file(paste(prefix, "segm.int", sep="")),
                              inla.read.fmesher.file(paste(prefix, "segm.int.grp", sep="")),
                              FALSE))
    } else {
        segm.int = NULL
    }

    mesh = (list(meta = (list(call=match.call(),
                              fmesher.args = all.args,
                              time = time.used,
                              prefix = prefix)),
                 s = mesh$s, tv = mesh$tv,
                 segm = list(bnd=segm.bnd, int=segm.int),
                 node.idx=idx))

    class(mesh) <- "fmesher"

    return(mesh)
}

`summary.fmesher` = function(x, ...)
{
    ## provides a summary for a mesh object
    if (!inherits(x, "fmesher"))
        stop("'x' must inherit from class \"fmesher\"")

    ret = list()
    ret = c(ret, list(call=x$meta$call))
    ret = c(ret, list(fmesher.args=x$meta$fmesher.args))
    ret = c(ret, list(mesh=x))

    class(ret) <- "summary.fmesher"
    return (ret)
}

`print.summary.fmesher` = function(x, ...)
{
    if (!inherits(x, "summary.fmesher"))
        stop("'x' must inherit from class \"summary.fmesher\"")

    cat("\nCall:\n", deparse(x$call), "\n", sep = "")

    cat("\nfmesher arguments:\n", x$fmesher.args, "\n", sep = "")

    if(!is.null(x$mesh)) {
        cat("\nVertices:\t", as.character(nrow(x$mesh$s)), "\n", sep="")
        cat("Triangles:\t", as.character(nrow(x$mesh$tv)), "\n", sep="")
        cat("\n")
    } else {
        cat("\nThe mesh is empty\n\n")
    }

    if(!is.null(x$mesh$segm)) {
        cat("Boundary (groups):\t",
            as.character(max(0,nrow(x$mesh$segm$bnd$idx))),
            " (", as.character(length(unique(x$mesh$segm$bnd$grp))), ")",
            "\n", sep="")
        cat("Interior (groups):\t",
            as.character(max(0,nrow(x$mesh$segm$int$idx))),
            " (", as.character(length(unique(x$mesh$segm$int$grp))), ")",
            "\n", sep="")
        cat("\n")
    } else {
        cat("No segment information\n\n")
    }

}
