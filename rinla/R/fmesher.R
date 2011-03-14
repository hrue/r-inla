inla.mesh.segm = function(...)
{
    warning("'inla.mesh.segm' is deprecated.  Use 'inla.mesh.segment' instead.")
    return(inla.mesh.segment(...))
}


inla.mesh.segment =
    function(loc=NULL,
             idx=NULL,
             grp=NULL,
             is.bnd=TRUE)
{
    if ((missing(loc) || is.null(loc)) &&
        (missing(idx) || is.null(idx)))
        stop("At most one of 'loc' and 'idx' may be missing or null.")
    if (!missing(loc) && !is.null(loc)) {
        if (!is.matrix(loc))
            stop("'loc' must be a a matrix")
        if (!is.double(loc))
            storage.mode(loc) = "double"
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
            idx = as.matrix(idx, nrow=length(idx), ncol=1)
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

    if (!is.null(loc)) {
        ## Identify unused locations and remap indices accordingly.
        idx.new = rep(0L, nrow(loc))
        idx.new[as.vector(idx)] = 1L
        loc = loc[idx.new==1L,, drop=FALSE]
        idx.new[idx.new==1L] = 1L:sum(idx.new)
        idx = (matrix(idx.new[as.vector(idx)],
                      nrow=nrow(idx),
                      ncol=ncol(idx)))
    }

    ret = list(loc=loc, idx=idx, grp=grp, is.bnd=is.bnd)
    class(ret) <- "inla.mesh.segment"
    return(ret)
}


lines.inla.mesh.segment = function(segm, loc=NULL, ...)
{
    if (!is.null(segm$loc))
        loc = segm$loc
    stopifnot(!is.null(loc), ncol(loc)>=2)

    grps = inla.ifelse(is.null(segm$grp), rep(0L,nrow(segm$idx)), segm$grp)
    for (grp in unique(grps)) {
        idx = which(grps==grp)
        lines(loc[t(cbind(segm$idx[idx,, drop=FALSE], NA)), 1],
              loc[t(cbind(segm$idx[idx,, drop=FALSE], NA)), 2],
              type="l",
              col=c("black", "blue", "red", "green")[1+(grp%%4)],
              ...)
    }
}

plot.inla.mesh =
    function(mesh,
             t.sub=1:nrow(mesh$graph$tv),
             add=FALSE,
             lwd=1,
             col="black",
             ...)
{
    inla.require.inherits(mesh, "inla.mesh", "'mesh'")

    idx = cbind(mesh$graph$tv[t.sub,c(1:3,1), drop=FALSE], NA)
    x = mesh$loc[t(idx), 1]
    y = mesh$loc[t(idx), 2]

    if (!add) {
        plot.new()
        plot.window(xlim=range(mesh$loc[,1]), ylim=range(mesh$loc[,2]), "")
    }
    lines(x, y, type="l", col=col, lwd=lwd)
    if (!is.null(mesh$segm$bnd))
        lines(mesh$segm$bnd, mesh$loc, lwd=lwd+1, ...)
    if (!is.null(mesh$segm$int))
        lines(mesh$segm$int, mesh$loc, lwd=lwd+1, ...)
    if (!add)
        title("Constrained refined Delaunay triangulation",
              deparse(substitute(mesh)))
    return(invisible())
}



inla.mesh.lattice =
    function(x=seq(0, 1, length.out=2),
             y=seq(0, 1, length.out=2),
             z=NULL,
             dims = (inla.ifelse(is.matrix(x),
                                 dim(x),
                                 c(length(x), length(y)))),
             units = NULL)
{
    units = match.arg(units, c("default", "longlat","longsinlat"))

    if (is.matrix(x)) {
        if (!identical(dims, dim(x)) ||
            !identical(dims, dim(y)) ||
            (is.matrix(z) && !identical(dims, dim(z))))
            stop("The size of matrices 'x', 'y', and 'z' must match 'dims'.")
        loc = cbind(as.vector(x), as.vector(y), as.vector(z))
        x = NULL
        y = NULL
    } else {
        loc = (cbind(rep(x, times = dims[2]),
                     rep(y, each = dims[1])))
    }
    if (!is.double(loc))
        storage.mode(loc) = "double"

    if (identical(units, "longlat")) {
        ## Transform onto a sphere
        loc = (cbind(cos(loc[,1]*pi/180)*cos(loc[,2]*pi/180),
                     sin(loc[,1]*pi/180)*cos(loc[,2]*pi/180),
                     sin(loc[,2]*pi/180)))
    } else if (identical(units, "longsinlat")) {
        ## Transform onto a sphere
        coslat = sapply(loc[,2], function(x) sqrt(max(0, 1-x^2)))
        loc = (cbind(cos(loc[,1]*pi/180)*coslat,
                     sin(loc[,1]*pi/180)*coslat,
                     loc[,2]
                     ))
    }

    ## Construct lattice boundary
    segm.idx = (c(1:(dims[1]-1),
                  dims[1]*(1:(dims[2]-1)),
                  dims[1]*dims[2]-(0:(dims[1]-2)),
                  dims[1]*((dims[2]-1):1)+1))
    segm.grp = (c(rep(1L, dims[1]-1),
                  rep(2L, dims[2]-1),
                  rep(3L, dims[1]-1),
                  rep(4L, dims[2]-1)))

    segm = (inla.mesh.segment(loc=loc[segm.idx,, drop=FALSE],
                              grp=segm.grp,
                              is.bnd=TRUE))

    lattice = list(dims=dims, x=x, y=y, loc=loc, segm=segm)
    class(lattice) = "inla.mesh.lattice"
    return(lattice)
}


extract.groups = function(...)
{
    UseMethod("extract.groups")
}

extract.groups.inla.mesh.segment =
    function(segm,
             groups,
             groups.new=groups,
             ...)
{
    inla.require.inherits(segm, "inla.mesh.segment", "'segm'")

    if (length(groups.new)==1L) {
        groups.new = rep(groups.new, length(groups))
    }
    if (length(groups.new)!=length(groups)) {
        stop("Length of 'groups.new' (", length(groups.new),
             ") does not match length of 'groups' (",length(groups),")")
    }

    idx = c()
    segm.grp = c()
    for (k in 1:length(groups)) {
        extract.idx = which(segm$grp==groups[k])
        idx = c(idx, extract.idx)
        segm.grp = c(segm.grp, rep(groups.new[k], length(extract.idx)))
    }
    segm.idx = segm$idx[idx,, drop=FALSE]

    return(inla.mesh.segment(loc=segm$loc,
                             idx=segm.idx,
                             grp=segm.grp,
                             segm$is.bnd))
}



inla.mesh.parse.segm.input =
    function(boundary=NULL,
             interior=NULL,
             segm.offset=0L,
             loc.offset=0L)
{
###########################################
    homogenise.segm.input = function(x, is.bnd)
    {
        if (is.matrix(x) || is.vector(x)) { ## Coordinates or indices
            x = (inla.ifelse(is.matrix(x),x,
                             as.matrix(x,nrow=length(x),ncol=1)))
            if (is.integer(x)) { ## Indices
                ret = inla.mesh.segment(NULL, x, NULL, is.bnd)
            } else if (is.numeric(x)) { ## Coordinates
                ret = inla.mesh.segment(x, NULL, NULL, is.bnd)
            } else {
                stop("Segment info matrix must be numeric or integer.")
            }
        } else if (inherits(x, "inla.mesh.segment")) {
            ## Override x$is.bnd:
            ret = inla.mesh.segment(x$loc, x$idx, x$grp, is.bnd)
        } else if (!is.null(x)) {
            inla.require.inherits(NULL,
                                  c("matrix", "inla.mesh.segment"),
                                  "Segment info")
        } else {
            ret = NULL
        }
        return(ret)
    }
##################################################
    homogenise.segm.grp = function(input) {
        grp.idx = 0L
        for (k in 1:length(input)) if (!is.null(input[[k]])) {
            inla.require.inherits(input[[k]],
                                  "inla.mesh.segment",
                                  "Segment info list members ")
            if (is.null(input[[k]]$grp)) {
                grp.idx = grp.idx+1L
                input[[k]] = (inla.mesh.segment(input[[k]][[1]],
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
    parse.segm.input = function(input, segm.offset=0L, loc.offset=0L)
    {
        loc = NULL
        bnd = list(loc=NULL, idx = matrix(,0,2), grp = matrix(,0,1))
        int = list(loc=NULL, idx = matrix(,0,2), grp = matrix(,0,1))
        storage.mode(bnd$idx) <- "integer"
        storage.mode(bnd$grp) <- "integer"
        storage.mode(int$idx) <- "integer"
        storage.mode(int$grp) <- "integer"
        for (k in 1:length(input)) if (!is.null(input[[k]])) {
            inla.require.inherits(input[[k]],
                                  "inla.mesh.segment",
                                  "Segment info list members ")
            if (!is.null(input[[k]]$loc)) {
                extra.loc.n = nrow(input[[k]]$loc)
                idx.offset = segm.offset
                segm.offset = segm.offset + extra.loc.n
                loc = (inla.ifelse(is.null(loc),
                                   input[[k]]$loc,
                                   rbind(loc, input[[k]]$loc)))
            } else {
                idx.offset = loc.offset
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
                                       inla.mesh.segment(bnd$loc,
                                                         bnd$idx,
                                                         bnd$grp,
                                                         TRUE))),
                    int = (inla.ifelse(is.null(int),
                                       NULL,
                                       inla.mesh.segment(int$loc,
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
    return(parse.segm.input(segm, segm.offset, loc.offset))
}






################################
##
## Old code.  Filtering is now done in fmesher itself.
## Retained for now so that we can check if the results are the same.
##
inla.mesh.filter.locations = function(loc, cutoff)
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
                              loc[loc.i,, drop=FALSE] -
                              node.coord[1:node.i.max,, drop=FALSE])^2))
            if (min(loc.to.node.dist) > cutoff) {
                node.i.max = node.i.max+1L
                node.coord[node.i.max,] = loc[loc.i,, drop=FALSE]
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
                              loc[loc.i,, drop=FALSE] -
                              node.coord)^2))
            node.i = which.min(loc.to.node.dist);
            map.loc.to.node[[loc.i]] = node.i
        }

        return(list(loc = node.coord, node.idx = map.loc.to.node))
    }
############################




inla.mesh = function(...)
{
    UseMethod("inla.mesh")
}


inla.mesh.default =
    function(loc=NULL, tv=NULL,
             boundary=NULL, interior=NULL,
             extend = (missing(tv) || is.null(tv)),
             refine=FALSE,
             lattice=NULL,
             globe=NULL,
             cutoff = 0,
             plot.delay = NULL,
             data.dir,
             keep = (!missing(data.dir) && !is.null(data.dir)))
{

    time.total = system.time({ ## Entire function timing start

    time.pre = system.time({ ## Pre-processing timing start

    if (!is.null(loc)) {
        if (!is.matrix(loc)) {
            loc = as.matrix(loc)
        }
        if (!is.double(loc)) {
            storage.mode(loc) = "double"
        }
    }

    if (is.logical(extend) && extend) extend = list()
    if (is.logical(refine) && refine) refine = list()

    if (missing(lattice) || is.null(lattice)) {
        lattice = list(loc=NULL, segm=NULL)
        lattice.n = 0L
    } else {
        inla.require.inherits(lattice, "inla.mesh.lattice", "'lattice'")
        if (!is.null(tv)) {
            warning("Both 'lattice' and 'tv' specified.  Ignoring 'tv'.")
            tv = NULL
        }
        if (!inherits(extend, "list")) {
            boundary = (c(inla.ifelse(inherits(boundary, "list"),
                                      boundary, list(boundary)),
                          list(lattice$segm)))
        }
        lattice.n = max(0L,nrow(lattice$loc))
    }
    loc.n = max(0L,nrow(loc))

    segm = (inla.mesh.parse.segm.input(boundary,
                                       interior,
                                       loc.n,
                                       0L))
    segm.n = max(0,nrow(segm$loc))
    ## Run parse again now that we know where the indices should point:
    segm = (inla.mesh.parse.segm.input(boundary,
                                       interior,
                                       0L,
                                       segm.n+lattice.n))

    if (TRUE) {
        loc0 = rbind(segm$loc, lattice$loc, loc)
        if ((!is.null(loc0)) && (nrow(loc0)>0))
            idx0 = 1:nrow(loc0)
        else
            idx0 = c()
    } else {
        ##
        ## Old code.  Filtering is now done in fmesher itself.
        ## Retained for now so that we can check if the results are the same.
        ##

        tmp = (inla.mesh.filter.locations(rbind(segm$loc, lattice$loc, loc),
                                          cutoff))
        loc0 = tmp$loc
        idx0 = tmp$node.idx

        ## Remap indices
        if (!is.null(segm$bnd)) {
            segm$bnd$idx = (matrix(idx0[(as.vector(segm$bnd$idx)+segm.n-1L) %%
                                        (segm.n+lattice.n+loc.n)+1L],
                                   nrow=nrow(segm$bnd$idx),
                                   ncol=ncol(segm$bnd$idx)))
        }
        if (!is.null(segm$int)) {
            segm$int$idx = (matrix(idx0[segm.n+as.vector(segm$int$idx)],
                                   nrow=nrow(segm$int$idx),
                                   ncol=ncol(segm$int$idx)))
        }
        if (!is.null(tv)) {
            tv = (matrix(idx0[segm.n+lattice.n+as.vector(tv)],
                         nrow=nrow(tv),
                         ncol=ncol(tv)))
        }
    }


    ## Where to put the files?
    if (keep) {
        if (missing(data.dir)) {
            data.dir = inla.fmesher.make.dir("inla.mesh.data")
        }
        prefix = paste(data.dir, "/mesh.", sep="")
        keep.dir = TRUE
    } else {
        if (missing(data.dir)) {
            prefix = inla.tempfile(pattern="fmesher", tmpdir=tempdir())
            prefix = paste(prefix, ".", sep="")
            keep.dir = TRUE
        } else {
            data.dir = inla.fmesher.make.dir(data.dir)
            prefix = paste(data.dir, "/mesh.", sep="")
            ## We should not try to delete the session tempdir...
            keep.dir = identical(tempdir(), dirname(prefix))
        }
    }
    prefix = inla.fmesher.make.prefix(NULL, prefix)

    if (is.null(loc0)) {
        all.args = ""
    } else {
        loc.file = fmesher.write(loc0, prefix, "input.s")
        all.args = "--input=input.s"

        if (!is.null(tv)) {
            fmesher.write(tv-1L, prefix, "input.tv")
            all.args = paste(all.args, ",input.tv", sep="")
        }
    }
    if (!missing(cutoff)) {
        all.args = paste(all.args, " --cutoff=", cutoff, sep="")
    }
    if (!is.null(segm$bnd)) {
        fmesher.write(segm$bnd$idx-1L, prefix, "input.segm.bnd.idx")
        fmesher.write(segm$bnd$grp, prefix, "input.segm.bnd.grp")
        all.args = paste(all.args," --boundary=input.segm.bnd.idx")
        all.args = paste(all.args," --boundarygrp=input.segm.bnd.grp")
    }
    if (!is.null(segm$int)) {
        fmesher.write(segm$int$idx-1L, prefix, "input.segm.int.idx")
        fmesher.write(segm$int$grp, prefix, "input.segm.int.grp")
        all.args = paste(all.args," --interior=input.segm.int.idx")
        all.args = paste(all.args," --interiorgrp=input.segm.int.grp")
    }
    if (!missing(globe) && !is.null(globe)) {
        all.args = paste(all.args, " --globe=", globe, sep="")
    }


    if (inherits(extend,"list")) {
        cet = c(0,0)
        cet[1] = inla.ifelse(is.null(extend$n), 16, extend$n)
        cet[2] = inla.ifelse(is.null(extend$offset), -0.1, extend$offset)
        all.args = (paste(all.args," --cet=",
                          cet[1],",", cet[2], sep=""))
    }
    if (inherits(refine,"list")) {
        rcdt = c(0,0,0)
        if (!missing(globe) && !is.null(globe)) {
            max.edge.default=10
        } else {
            max.edge.default = (sqrt(diff(range(loc0[,1]))^2+
                                     diff(range(loc0[,2]))^2+
                                     inla.ifelse(ncol(loc0)<3,
                                                 0,
                                                 diff(range(loc0[,3]))^2)
                                     ))
        }
        if ((inherits(extend,"list")) && (!is.null(extend$offset))) {
            max.edge.default = (max.edge.default +
                                max(0,2*extend$offset))
            max.edge.default = (max.edge.default *
                                (1+max(0,-2*extend$offset)))
        }
        max.edge.default = max.edge.default*10 ## "*10": Better to be safe.
        rcdt[1] = inla.ifelse(is.null(refine$min.angle), 21, refine$min.angle)
        rcdt[2] = (inla.ifelse(is.null(refine$max.edge),
                               max.edge.default,
                               refine$max.edge))
        rcdt[3] = (inla.ifelse(is.null(refine$max.edge),
                               max.edge.default,
                               refine$max.edge))
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

    }) ## Pre-processing timing end

    ## Call fmesher:
    time.fmesher = system.time({
        echoc = (inla.fmesher.call(all.args=all.args,
                                   prefix=prefix))
    })

    time.post = system.time({ ## Post-processing timing start

    ## Read the mesh:
    manifold = 1L+fmesher.read(prefix, "manifold")
    manifold = list("M", "R2", "S2")[[manifold]]

    loc = fmesher.read(prefix, "s")
    graph = (list(tv = 1L+fmesher.read(prefix, "tv"),
                  vt = 1L+fmesher.read(prefix, "vt"),
                  tt = 1L+fmesher.read(prefix, "tt"),
                  tti = 1L+fmesher.read(prefix, "tti"),
                  vv = fmesher.read(prefix, "vv")))
    graph$tv[graph$tv==0L] = NA
    graph$vt[graph$vt==0L] = NA
    graph$tt[graph$tt==0L] = NA
    graph$tti[graph$tti==0L] = NA

    ## Read the vertex input/output mapping:
    idx.all = 1L+fmesher.read(prefix, "idx")
    idx = (list(loc = (inla.ifelse(loc.n>0,
                                   idx.all[idx0[segm.n+lattice.n+(1:loc.n)]],
                                   NULL)),
                lattice = (inla.ifelse(lattice.n>0,
                                    idx.all[idx0[segm.n+(1:lattice.n)]],
                                    NULL)),
                segm = (inla.ifelse(segm.n>0,
                                    idx.all[idx0[(1:segm.n)]],
                                    NULL))))

    ## Read constraint segment information:
    segm.bnd = (inla.mesh.segment(NULL,
                                  1L+fmesher.read(prefix, "segm.bnd.idx"),
                                  fmesher.read(prefix, "segm.bnd.grp"),
                                  TRUE))
    segm.int = (inla.mesh.segment(NULL,
                                  1L+fmesher.read(prefix, "segm.int.idx"),
                                  fmesher.read(prefix, "segm.int.grp"),
                                  FALSE))

    if (!keep)
        unlink(paste(prefix, "*", sep=""), recursive=FALSE)
    if (!keep.dir) {
        unlink(dirname(prefix), recursive=TRUE)
    }

    }) ## Post-processing timing end

    time.object = system.time({ ## Object construction timing start
    mesh = (list(meta = (list(call=match.call(),
                              fmesher.args = all.args,
                              time = (rbind(pre = time.pre,
                                            fmesher = time.fmesher,
                                            post = time.post)),
                              prefix = prefix)),
                 manifold = manifold,
                 n = nrow(loc),
                 loc = loc,
                 graph = graph,
                 segm = list(bnd=segm.bnd, int=segm.int),
                 idx = idx))
    class(mesh) <- "inla.mesh"

    }) ## Object construction timing end

    }) ## Entire function timing start

    mesh$meta$time = (rbind(mesh$meta$time,
                            object=time.object,
                            total=time.total))
    return(mesh)
}


inla.mesh.inla.mesh = function(mesh, ...)
{
    inla.require.inherits(mesh, "inla.mesh", "'mesh'")

    not.known = function (mesh, queryname)
    {
        stop(paste("Query '", queryname,
                   "' unknown.", sep=""))
    }
    not.implemented = function (mesh, queryname)
    {
        stop(paste("Query '", queryname,
                   "' not implemented for inla.mesh.", sep=""))
##        stop(paste("Query '", queryname,
##                   "' not implemented for inla.mesh for mesh type '",
##                   model$type, "'.", sep=""))
    }

    result = list()
    queries = inla.parse.queries(...)
    if (length(queries)==0L)
        return(result)

    for (query.idx in 1:length(queries)) {
        query = names(queries)[query.idx]
        param = queries[[query.idx]]
        answer = NULL
        query = (match.arg(query, c("tt.neighbours",
                                    "vt.neighbours"
                                    )))
        if (identical(query, "tt.neighbours")) {
            if (is.null(param))
                param = list(c(1))
            if (length(param)<2)
                param = c(as.list(param), list(c(1)))
            nT = nrow(mesh$graph$tv)
            i = rep(1:nT,3)
            j = as.vector(mesh$graph$tt)
            i = i[!is.na(j)]
            j = j[!is.na(j)]
            tt = sparseMatrix(i=i, j=j, x=rep(1,length(i)), dims=c(nT,nT))
            answer = (sparseMatrix(i=param[[1]],
                                   j=rep(1,length(param[[1]])),
                                   x=1,
                                   dims=c(nT,1)))
            tt0 = (answer == 0.5)*1
            tt1 = answer
            order = 0
            while (order<min(param[[2]])) {
                order = order+1
                tt0 = tt1
                tt1 = answer
                answer = ((((tt %*% answer) > 0) - tt1 - tt0) > 0)
            }
            while (order<max(param[[2]])) {
                order = order+1
                answer = ((((tt %*% answer) > 0) - tt1 - tt0) > 0)
            }

            ##                not.implemented(mesh,query)
        } else if (identical(query, "vt.neighbours")) {
            if (is.null(param))
                param = list(1)
            if (length(param)<2)
                param = c(as.list(param), list(c(1)))
            nV = nrow(mesh$loc)
            nT = nrow(mesh$graph$tv)
            i = rep(1:nT,3)
            j = as.vector(mesh$graph$tv)
#            i = i[!is.na(j)]
#            j = j[!is.na(j)]
            tv = sparseMatrix(i=i, j=j, x=rep(1,length(i)), dims=c(nT,nV))
            vv = (sparseMatrix(i=param[[1]],
                               j=rep(1,length(param[[1]])),
                               x=1,
                               dims=c(nV,1)))
            vt = (tv %*% vv ) > 0
            vv0 = (vv == 0.5)*1
            vv1 = vv
            vt0 = (tv %*% vv0 ) > 0
            vt1 = vt
            order = 0
            while (order<min(param[[2]])) {
                order = order+1
                vv0 = vv1
                vv1 = vv
                vv = ((((mesh$graph$vv %*% vv) > 0) - vv1 - vv0) > 0)
                vt0 = vt1
                vt1 = vt
                vt = (((tv %*% vv) > 0) - vt1 - vt0 ) > 0
            }
            while (order<max(param[[2]])) {
                order = order+1
                vv = ((((mesh$graph$vv %*% vv) > 0) - vv1 - vv0) > 0)
                vt = (((((tv %*% vv) > 0) - vt1 - vt0 ) > 0 ) + vt) > 0
            }
            answer = vt

            ##                not.implemented(mesh,query)
        } else if (!identical(query, "")) {
            not.known(mesh,query)
        }
        ## Expand the result list:
        result[query.idx] = list(NULL)
        names(result)[query.idx] = query
        ## Set the answer:
        if (!is.null(answer))
            result[[query.idx]] = answer
    }

    return(result)
}

summary.inla.mesh = function(x, verbose=FALSE, ...)
{
    ## provides a summary for a mesh object
    inla.require.inherits(x, "inla.mesh", "'x'")

    ret = list(verbose=verbose)
    if (verbose) {
        ret = (c(ret, list(call=x$meta$call,
                           fmesher.args=x$meta$fmesher.args,
                           prefix=x$meta$prefix,
                           time = x$meta$time)))
    }
    ret = (c(ret, list(manifold=x$manifold,
                       nV=x$n,
                       nT=nrow(x$graph$tv),
                       xlim=range(x$loc[,1]),
                       ylim=range(x$loc[,2]),
                       zlim=range(x$loc[,3]))))

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
        ret = (c(ret, list(segm.bnd=my.segm(x$segm$bnd),
                           segm.int=my.segm(x$segm$int))))
    } else {
        ret = (c(ret, list(segm.bnd=my.segm(NULL),
                           segm.int=my.segm(NULL))))
    }

    class(ret) <- "summary.inla.mesh"
    return (ret)
}

print.summary.inla.mesh = function(x, ...)
{
    my.print.proc_time = function (x, ...)
    {
        if (!is.matrix(x)) {
            y = matrix(x,1,5)
        } else {
            y = x
        }
        for (k in 1:nrow(x)) {
            if (!is.na(y[k,4L]))
                y[k,1L] <- y[k,1L] + y[k,4L]
            if (!is.na(y[k,5L]))
                y[k,2L] <- y[k,2L] + y[k,5L]
        }
        y <- y[,1L:3L]
        colnames(y) <- c(gettext("user"), gettext("system"), gettext("elapsed"))
        print(y, ...)
        invisible(x)
    }


    inla.require.inherits(x, "summary.inla.mesh", "'x'")

    if (x$verbose) {
        cat("\nCall:\n")
        print(x$call)

        cat("\nfmesher:\t", x$fmesher.args, "\n", sep = "")
        cat("prefix:\t\t", x$prefix, "\n", sep = "")

        cat("\nTimings:\n")
        my.print.proc_time(x$time)
    }

    cat("\nManifold:\t", x$manifold, "\n", sep="")
    cat("Vertices:\t", as.character(x$nV), "\n", sep="")
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

    invisible(x)
}



inla.mesh.project = function(...)
{
    UseMethod("inla.mesh.project")
}

inla.mesh.project.inla.mesh = function(mesh, loc)
{
    inla.require.inherits(mesh, "inla.mesh", "'mesh'")

    smorg = (inla.fmesher.smorg(mesh$loc,
                                mesh$graph$tv,
                                points2mesh=loc))
    ti = smorg$p2m.t
    b = smorg$p2m.b

    ok = (ti[,1] > 0L)
    ti[ti[,1] == 0L,1] = NA

    ii = which(ok)
    A = (sparseMatrix(dims=c(nrow(loc),mesh$n),
                      i = rep(ii, 3),
                      j = as.vector(mesh$graph$tv[ti[ii,1],]),
                      x = as.vector(b[ii,]) ))

    return (list(t=ti, bary=b, A=A, ok=ok))
}


inla.mesh.project.inla.mesh.projector =
    function(projector, field)
{
    inla.require.inherits(projector, "inla.mesh.projector", "'projector'")

    if (is.data.frame(field)) {
        field = as.matrix(field)
    }

    if (is.null(dim(field))) {
        data = as.vector(projector$proj$A %*% as.vector(field))
        data[!projector$proj$ok] = NA
        return(matrix(data,
                      projector$lattice$dims[1],
                      projector$lattice$dims[2]))
    } else {
        data = projector$proj$A %*% field
        data[!projector$proj$ok,,drop=FALSE] = NA
        return(data)
    }
}


inla.mesh.projector =
    function(mesh,
             xlim=range(mesh$loc[,1]),
             ylim=range(mesh$loc[,2]),
             dims=c(100,100),
             projection=NULL)
{
    inla.require.inherits(mesh, "inla.mesh", "'mesh'")

    if (identical(mesh$manifold, "R2")) {
        units = "default"
        x = seq(xlim[1], xlim[2], length.out=dims[1])
        y = seq(ylim[1], ylim[2], length.out=dims[2])
    } else if (identical(mesh$manifold, "S2")) {
        projection = match.arg(projection, c("longlat", "longsinlat"))
        units = projection
        if (missing(xlim) || is.null(xlim)) {
            xlim = c(-180,180)
        }
        if (missing(ylim) || is.null(ylim)) {
            ylim = c(-90,90)
        }
        x = seq(xlim[1], xlim[2], length.out=dims[1])
        if (identical(projection, "longlat")) {
            y = seq(ylim[1], ylim[2], length.out=dims[2])
        } else {
            y = (seq(sin(ylim[1]*pi/180),
                     sin(ylim[2]*pi/2),
                     length.out=dims[2]))
        }
    }

    lattice = (inla.mesh.lattice(x=x, y=y, units = units))

    proj = inla.mesh.project(mesh, lattice$loc)
    projector = list(x=x, y=y, lattice=lattice, proj=proj)
    class(projector) = "inla.mesh.projector"

    return (projector)
}


inla.mesh.basis =
    function(mesh,
             type="b.spline",
             n=3,
             degree=2,
             knot.placement="uniform.area",
             rot.inv=TRUE)
{
    inla.require.inherits(mesh, "inla.mesh", "'mesh'")

    type = match.arg(type, c("b.spline", "sph.harm"))
    knot.placement = (match.arg(knot.placement,
                                c("uniform.area",
                                  "uniform.latitude")))

    if (identical(type, "b.spline")) {
        if (identical(mesh$manifold, "R2")) {
            long = ((mesh$loc[,1]-min(mesh$loc[,1]))/
                    diff(range(mesh$loc[,1])))*90*pi/180
            sinlat =  (((mesh$loc[,2]-min(mesh$loc[,2]))/
                       diff(range(mesh$loc[,2])))*2-1)
            coslat = sapply(sinlat, function(x) sqrt(max(0.0,1.0-x^2)))
            loc = (matrix(c(cos(long)*coslat,
                            sin(long)*coslat,
                            sinlat), mesh$n, 3))
            knots = 0
        } else if (identical(mesh$manifold, "S2")) {
            loc = mesh$loc
            knots = identical(knot.placement, "uniform.latitude")
        } else {
            stop("Only know how to make B-splines on R2 and S2.")
        }
        degree = max(0L, min(n-1L, degree))
        basis = (inla.fmesher.smorg(loc,
                                    mesh$graph$tv,
                                    bspline = c(n, degree, knots),
                                    fem=-1)$bspline)
        if (!rot.inv) {
            warning("Currently only 'rot.inv=TRUE' is supported for B-splines.")
        }
    } else if (identical(type, "sph.harm")) {
        if (!identical(mesh$manifold, "S2")) {
            stop("Only know how to make spherical harmonics on S2.")
        }
        if (rot.inv) {
            basis = (inla.fmesher.smorg(mesh$loc,
                                        mesh$graph$tv,
                                        sph0=n)$sph0)
        } else {
            basis = (inla.fmesher.smorg(mesh$loc,
                                        mesh$graph$tv,
                                        sph=n)$sph)
        }
    }

    return(basis)
}





old.mesh.class = function(...)
{
    UseMethod("old.mesh.class")
}

old.mesh.class.inla.mesh = function(mesh, ...)
{
    fmesh=list(mesh=mesh)
    fmesh$mesh$s = mesh$loc
    fmesh$mesh$tv = mesh$graph$tv
    fmesh$locations.idx = mesh$idx$loc
    class(fmesh)="inla.fmesher.mesh"
    return(fmesh)
}




inla.spde = function(...)
{
    UseMethod("inla.spde")
}

inla.spde.inla.mesh =
    function(mesh,
             model=c("matern", "imatern", "matern.osc"),
             param=NULL,
             ...)
{
    inla.require.inherits(mesh, "inla.mesh", "'mesh'")

    model = match.arg(model)
    if (is.null(param))
        param = list()

    spde.prefix = inla.fmesher.make.prefix(NULL, NULL)

    spde = (list(model = model,
                 mesh = mesh,
                 internal = list(),
                 f = (list(model="spde",
                           spde.prefix=spde.prefix,
                           n=nrow(mesh$loc)
                           ))
                 ))
    class(spde) = "inla.spde"

    if (identical(model, "matern") ||
        identical(model, "imatern") ||
        identical(model, "matern.osc")
        ) {

        if (is.null(param))
            param = list()
        if (is.null(param$alpha))
            param$alpha = 2
        if (is.null(param$basis.T))
            param$basis.T = matrix(1, mesh$n, 1)
        else if (!is.matrix(param$basis.T)) {
            len = length(as.vector(param$basis.T))
            if (len == 1L)
                param$basis.K = matrix(as.vector(param$basis.T), mesh$n, 1)
            else
                param$basis.T = as.matrix(param$basis.T)
        }
        if (nrow(param$basis.T) != mesh$n)
            stop(paste("'basis.T' has ", nrow(basis.T),
                       " rows; expected ", mesh$n, ".", sep=""))
        if (identical(model, "matern") ||
            identical(model, "matern.osc")
            ) {
            if (is.null(param$basis.K)) {
                param$basis.K = matrix(1, mesh$n, 1)
            } else if (!is.matrix(param$basis.K)) {
                len = length(as.vector(param$basis.K))
                if (len == 1L)
                    param$basis.K = matrix(as.vector(param$basis.K), mesh$n, 1)
                else
                    param$basis.K = as.matrix(param$basis.K)
            }
        } else {
            param$basis.K = matrix(0, mesh$n, 1)
        }
        if (nrow(param$basis.K) != mesh$n)
            stop(paste("'basis.K' has ", nrow(basis.K),
                       " rows; expected ", mesh$n, ".", sep=""))
        spde$internal = (c(spde$internal,
                           list(alpha = param$alpha,
                                basis.T = param$basis.T,
                                basis.K = param$basis.K)
                           ))

        mesh.range = (max(c(diff(range(mesh$loc[,1])),
                            diff(range(mesh$loc[,2])),
                            diff(range(mesh$loc[,3]))
                            )))

        spde$internal = (c(spde$internal,
                           inla.fmesher.smorg(mesh$loc,
                                              mesh$graph$tv,
                                              fem=2,
                                              output=list("c0", "g1", "g2"))))

        if (param$alpha==2) {
            kappa0 = sqrt(8)/(mesh.range*0.2)
            tau0 = 1/sqrt(4*pi*kappa0^2)/1.0
        } else if (param$alpha==1) {
            spde$internal$g2 = spde$internal$g1
            spde$internal$g1 = spde$internal$g1*0.0

            kappa0 = sqrt(sqrt(8)/(mesh.range*0.2))
            tau0 = 1/sqrt(4*pi)/1.0
        }
        ## inla checks PREFIX valididy by looking for "s":
        fmesher.write(spde$mesh$loc, spde.prefix, "s")
        ## Write the precision building blocks:
        fmesher.write(spde$internal$c0, spde.prefix, "c0")
        fmesher.write(spde$internal$g1, spde.prefix, "g1")
        fmesher.write(spde$internal$g2, spde.prefix, "g2")
        fmesher.write(spde$internal$basis.T, spde.prefix, "basisT")
        fmesher.write(spde$internal$basis.K, spde.prefix, "basisK")

        if (identical(model, "matern")) {
            spde$f$hyper = (list(theta1=(list(initial=log(tau0),
                                              param=c(log(tau0), 0.01))),
                                 theta2=(list(initial=log(kappa0^2),
                                              param=c(log(kappa0^2), 0.1))),
                                 theta4=(list(initial=log(tau0),
                                              param=c(log(tau0), 0.01)))
                                 ))
        } else if (identical(model, "imatern")) {
            spde$f$hyper = (list(theta1=(list(initial=log(tau0),
                                              param=c(log(tau0), 0.01))),
                                 theta2=(list(initial=-20,
                                              fixed=TRUE)),
                                 theta4=(list(initial=log(tau0),
                                              param=c(log(tau0), 0.01)))
                                 ))
        } else if (identical(model, "matern.osc")) {
            spde$f$hyper = (list(theta1=(list(initial=log(tau0),
                                              param=c(log(tau0), 0.01))),
                                 theta2=(list(initial=log(kappa0^2),
                                              param=c(log(kappa0^2), 0.1))),
                                 theta3=(list(fixed=FALSE,
                                              initial=0,
                                              param=c(0, 0.01))),
                                 theta4=(list(initial=log(tau0),
                                              param=c(log(tau0), 0.01)))
                                 ))
        }
    } else {
        stop(paste("Model '", model, "' unknown or not implemented.", sep=""))
    }

    return(invisible(spde))
}




inla.parse.queries = function(...)
{
    queries = list(...)
    if (length(queries)==0)
        return(queries)

    ## Make sure that we have a list of names, empty or not:
    if (is.null(names(queries))) {
        q.names = rep("", length(queries))
    } else {
        q.names = names(queries)
    }

    ## All nameless entries must be strings with query names.  Replace
    ## empty names with those names, and set those entries to NULL.
    for (query.idx in 1:length(queries)) {
        if (q.names[[query.idx]]=="") {
            if (is.character(queries[[query.idx]])) {
                names(queries)[query.idx] = queries[[query.idx]]
                queries[query.idx] = list(NULL)
            } else {
                queries[query.idx] = list(NULL)
                warning(paste("Unnamed query ignored.  Check query #",
                              query.idx, sep=""))
            }
        }
    }

    return(queries)
}



inla.spde.inla.spde = function(spde, ...)
{
    inla.require.inherits(spde, "inla.spde", "'spde'")

    not.known = function (spde, queryname)
    {
        stop(paste("Query '", queryname,
                   "' unknown.", sep=""))
    }
    not.implemented = function (spde, queryname)
    {
        stop(paste("Query '", queryname,
                   "' not implemented for inla.spde model '",
                   spde$model, "'.", sep=""))
    }
    param.to.fcn =
        function(basis, theta, values, n, name)
        {
            if (!is.null(values)) {
                if ((is.vector(values) &&
                     (length(values)==n)) ||
                    (is.matrix(values) &&
                     (nrow(values)==n))) {
                    fcn = as.vector(values)
                } else {
                    fcn = as.vector(values)
                }
                if (length(fcn) == 1L)
                    fcn = rep(fcn, n)
                else if (length(fcn) != n)
                    stop(paste("Length of '", name, "' is ", length(fcn),
                               ", should be ", n, ".", sep=""))
            } else {
                fcn = exp(basis %*% theta)
            }
            return(fcn)
        }

    result = list()
    queries = inla.parse.queries(...)
    if (length(queries)==0L)
        return(result)

    for (query.idx in 1:length(queries)) {
        query = names(queries)[query.idx]
        param = queries[[query.idx]]
        answer = NULL
        query = (match.arg(query, c("precision",
                                    "sample"
                                    )))
        if (identical(query, "precision")) {
            if (identical(spde$model, "matern")) {
                tau = (param.to.fcn(spde$internal$basis.T,
                                    param$theta.T, param$tau,
                                    spde$mesh$n, "tau"))
                kappa2 = (param.to.fcn(spde$internal$basis.K,
                                       param$theta.K, param$kappa2,
                                       spde$mesh$n, "kappa2"))
                dT = Diagonal(spde$mesh$n, tau)
                dK2 = Diagonal(spde$mesh$n, kappa2)
                tmp = dK2 %*% spde$internal$g1
                answer = (dT %*% (dK2 %*% spde$internal$c0 %*% dK2+
                                  tmp + t(tmp) +
                                  spde$internal$g2) %*% dT)
            } else if (identical(spde$model, "matern.osc")) {
                tau = (param.to.fcn(spde$internal$basis.T,
                                    param$theta.T, param$tau,
                                    spde$mesh$n, "tau"))
                kappa2 = (param.to.fcn(spde$internal$basis.K,
                                       param$theta.K, param$kappa2,
                                       spde$mesh$n, "kappa2"))
                osc = param$osc
                dT = Diagonal(spde$mesh$n, tau)
                dK2 = Diagonal(spde$mesh$n, kappa2)
                tmp = dK2 %*% spde$internal$g1
                answer = (dT %*% (dK2 %*% spde$internal$c0 %*% dK2+
                                  cos(pi*osc)*(tmp + t(tmp)) +
                                  spde$internal$g2) %*% dT)
            } else if (identical(spde$model, "imatern")) {
                tau = (exp(param.to.fcn(spde$internal$basis.T,
                                        param$theta.T, param$tau,
                                        spde$mesh$n, "tau")))
                dT = Diagonal(spde$mesh$n, tau)
                answer = (dT %*% ((1e-10)*spde$internal$c0+
                                  spde$internal$g2) %*% dT)
            } else {
                not.implemented(spde,query)
            }
        } else if (identical(query, "sample")) {
            Q = inla.spde(spde, precision=param)$precision
            finn = (inla.finn(Q, seed=(inla.ifelse(is.null(param$seed),
                                                   0L,
                                                   param$seed))))
            answer = finn$sample
        } else if (!identical(query, "")) {
            not.known(spde,query)
        }
        ## Expand the result list:
        result[query.idx] = list(NULL)
        names(result)[query.idx] = query
        ## Set the answer:
        if (!is.null(answer))
            result[[query.idx]] = answer
    }

    return(result)
}


inla.spde.inla = function(inla, name, spde, ...)
{
    inla.require.inherits(inla, "inla", "'inla'")
    warning("'inla.spde.inla' is not implemented yet.")

    return(spde)
}







`inla.fmesher.smorg` =
    function(loc, tv,
             fem=NULL,
             aniso=NULL,
             gradients=FALSE,
             sph0=NULL,
             sph=NULL,
             bspline=NULL,
             points2mesh=NULL,
             output=NULL,
             keep=FALSE)
{
    prefix = inla.fmesher.make.prefix(NULL, NULL)

    n = nrow(loc)
    s.dim = ncol(loc)

    if (s.dim==1)
        stop("1-D models not implemented yet.")
    stopifnot(s.dim>=2, s.dim<=3)

    if (missing(output) || is.null(output)) {
        output.given = FALSE
        output = NULL
    } else {
        output.given = TRUE
    }
    output.fem = list("c0", "g1", "g2")
    output.aniso = list("g1aniso", "g2aniso")
    output.gradients = list("dx", "dy", "dz")
    output.sph0 = list("sph0")
    output.sph = list("sph")
    output.bspline = list("bspline")
    output.p2m = list("p2m.t", "p2m.b")

    fmesher.write(loc, prefix, "s")
    fmesher.write(tv-1L, prefix, "tv")
    all.args = "--smorg --input=s,tv"

    ## additional arguments
    if (!is.null(fem)) {
        all.args = paste(all.args," --fem=", max(2,fem), sep="")
        if (!output.given) output = c(output, output.fem)
    }
    if (!is.null(aniso)) {
        fmesher.write(aniso[[1]], prefix, "aniso.gamma")
        fmesher.write(aniso[[2]], prefix, "aniso.vec")

        all.args = paste(all.args," --aniso=aniso.gamma,aniso.vec", sep="")
        if (!output.given) output = c(output, output.aniso)
    }
    if (!is.null(gradients) && gradients) {
        all.args = paste(all.args," --grad", sep="")
        if (!output.given) output = c(output, output.gradients)
    }
    if (!is.null(sph0)) {
        all.args = paste(all.args," --sph0=", sph0, sep="")
        if (!output.given) output = c(output, output.sph0)
    }
    if (!is.null(sph)) {
        all.args = paste(all.args," --sph=", sph, sep="")
        if (!output.given) output = c(output, output.sph)
    }
    if (!is.null(bspline)) {
        all.args = (paste(all.args, " --bspline=",
                          bspline[1], ",", bspline[2], ",", bspline[3],
                          sep=""))
        if (!output.given) output = c(output, output.bspline)
    }
    if (!is.null(points2mesh)) {
        fmesher.write(points2mesh, prefix, "p2m")

        all.args = paste(all.args," --points2mesh=p2m", sep="")
        if (!output.given) output = c(output, output.p2m)
    }
    all.args = paste(all.args, inla.getOption("fmesher.arg"))

    echoc = inla.fmesher.call(all.args=all.args, prefix=prefix)

    result = list()
    for (name in output) {
        if (identical(name, "p2m.t"))
            if (!file.exists(paste(prefix, name, sep="")))
                result[[name]] = fmesher.read(prefix, "points2mesh.t")+1L
            else
                result[[name]] = fmesher.read(prefix, name)+1L
        else if (identical(name, "p2m.b"))
            if (!file.exists(paste(prefix, name, sep="")))
                result[[name]] = fmesher.read(prefix, "points2mesh.b")
            else
                result[[name]] = fmesher.read(prefix, name)
        else
            result[[name]] = fmesher.read(prefix, name)
    }

    if (!keep)
        unlink(paste(prefix, "*", sep=""), recursive=FALSE)

    return (result)
}





##inla.spde(mesh, model=list("matern"), ...)
##inla.spde(mesh, model=list("heat", Qw=..., t=...), ...)
##inla.spde(mesh, model=list("imatern"), ...)


