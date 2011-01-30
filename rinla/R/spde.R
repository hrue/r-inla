## create a spde model using fmesher.

`inla.fmesher.call` = function(
        fmesher.call = inla.getOption("fmesher.call"),
        all.args, prefix)
{
    if (inla.os("linux") || inla.os("mac")) {
        echoc = system(paste(shQuote(fmesher.call), all.args, prefix))
    }
    else if(inla.os("windows")) {
        if (!remote) {
            echoc = try(system(paste(shQuote(fmesher.call),
                    all.args,
                    prefix)), silent=TRUE)
            echoc = 0
        } else {
            echoc = try(inla.cygwin.run.command(
                    paste(inla.cygwin.map.filename(fmesher.call),
                          all.args,
                          inla.cygwin.map.filename(prefix))), silent=TRUE)
            echoc = 0
        }
    }
    else
        stop("\n\tNot supported architecture.")
}

`inla.fmesher.make.prefix` = function(dir = tempdir(), prefix = NULL)
{
    if (is.null(prefix)) {
        if (is.null(dir))
            prefix = paste(tempfile(), ".", sep="")
        else {
            inla.dir.create(dir)
            prefix = paste(tempfile(tmpdir=dir), ".", sep="")
        }
    } else {
        inla.dir.create(inla.dirname(prefix))
    }
    return (prefix)
}

`inla.create.spde` = function(mesh=NULL, prefix = NULL,
        fem=NULL, sph0=NULL, sph=NULL, bspline=NULL, points2mesh=NULL)
{
    ## create files using fmesher

    if (is.null(mesh) && is.null(prefix))
        stop("At least one of mesh and prefix must be specified.")

    prefix = inla.fmesher.make.prefix(NULL, prefix)

    if (is.null(mesh)) {
        ## Need to know the size of the graph.
        s = inla.read.fmesher.file(paste(prefix, "s", sep=""))
        stopifnot(is.matrix(s))
        n = nrow(s)
        s.dim = ncol(s)
    } else {
        stopifnot(is.matrix(mesh$s))
        n = nrow(mesh$s)
        s.dim = ncol(mesh$s)
    }

    if (s.dim==1)
        stop("1-D models not implemented yet.")
    stopifnot(s.dim>=2, s.dim<=3)

    if (!is.null(mesh)) {
        stopifnot(is.matrix(mesh$tv))
        inla.write.fmesher.file(mesh$s, filename = paste(prefix, "s", sep=""))
        inla.write.fmesher.file(inla.affirm.integer(mesh$tv-1L),
                                filename = paste(prefix, "tv", sep=""))
    }
    all.args = "--input=s,tv"

    ## additional arguments
    if (!is.null(fem)) {
        all.args = paste(all.args," --fem=", fem, sep="")
    }
    if (!is.null(sph0)) {
        all.args = paste(all.args," --sph0=", sph0, sep="")
    }
    if (!is.null(sph)) {
        all.args = paste(all.args," --sph=", sph, sep="")
    }
    if (!is.null(bspline)) {
        all.args = paste(all.args, " --bspline=",
                bspline[1], ",", bspline[2], ",", bspline[3],
                sep="")
    }
    if (!is.null(points2mesh)) {
        inla.write.fmesher.file(points2mesh,
                                filename = paste(prefix,
                                                 "points2mesh",
                                                 sep=""))

        all.args = paste(all.args," --points2mesh=points2mesh", sep="")
    }
    all.args = paste(all.args, inla.getOption("fmesher.arg"))

    echoc = inla.fmesher.call(all.args=all.args, prefix=prefix)

    if (is.null(mesh)) {
      mesh = list(s=s)
    }

    return (list(prefix = prefix, locations = mesh$s, n = n))
}

`inla.fmesher.mesh` = function(
        locations,
        dir = tempdir(),
        prefix = NULL,
        rcdt = c(21,-1.0,-0.5),
        cet = c(8,-0.1),
        boundary = NULL,
        interior = NULL,
        x11.delay = NULL,
        min.input.distance = 0)
{

    ## Internal helper function:
    handle.segment.input = function (input,is.boundary,loc.n,group.index)
    {
        k = 1
        previous.was.locations = FALSE
        locations = NULL
        indices = NULL
        groups = NULL
        while (k <= length(input)) {
            if ((class(input[[k]]) != "matrix") &&
                (class(input[[k]]) != "integer")) {
                stop("Input information must be a list of matrices or node indices.")
            }
            
            if ((class(input[[k]]) == "matrix") &&
                (class(input[[k]][1,1]) == "numeric")) {
                ## Input coordinates
                previous.locations.length = loc.n
                loc.n = loc.n + nrow(input[[k]])
                locations = rbind(locations,input[[k]])
                ## Look-ahead: if the next matrix is not an index vector,
                ## insert one.
                if (!((k < length(input)) &&
                      ((class(input[[k+1]]) == "integer") ||
                       ((class(input[[k+1]]) == "matrix") &&
                        (class(input[[k+1]][1,1]) == "integer"))))) {
                    if (k < length(input))
                        if (is.boundary)
                            input = (c(input[1:k],
                                       list(c(1:nrow(input[[k]]),1L)),
                                       input[(k+1):length(input)]))
                        else
                            input = (c(input[1:k],
                                       list(c(1:nrow(input[[k]]))),
                                       input[(k+1):length(input)]))
                    else
                        if (is.boundary)
                            input = (c(input[1:k],
                                       list(c(1:nrow(input[[k]]),1L))))
                        else
                            input = (c(input[1:k],
                                       list(c(1:nrow(input[[k]])))))
                }
                previous.was.locations = TRUE
                k = k+1
            }
            if ((class(input[[k]]) == "integer") ||
                ((class(input[[k]]) == "matrix") &&
                 (class(input[[k]][1,1]) == "integer"))) {
                ## Input node indices
                if (class(input[[k]]) != "matrix") {
                    input[[k]] = matrix(input[[k]],length(input[[k]]),1)
                }
                ## If the previous element was locations, shift indices.
                if (previous.was.locations) {
                    input[[k]] = (input[[k]]+
                                  previous.locations.length);
                }

                indices.n = nrow(input[[k]])
                indices = rbind(indices,input[[k]])

                ## Look-ahead: if the next matrix is not a group vector,
                ## insert one.
                if (!((k < length(input)) &&
                      ((class(input[[k+1]]) == "integer") ||
                       ((class(input[[k+1]]) == "matrix") &&
                        (class(input[[k+1]][1,1]) == "integer"))))) {
                    group.index = group.index+1L
                    if (k < length(input))
                        input = (c(input[1:k],
                                   list(group.index),
                                   input[(k+1):length(input)]))
                    else
                        input = (c(input[1:k],
                                   list(group.index)))
                }
                previous.was.locations = FALSE
                k = k+1

                ## Groups:
                if (class(input[[k]]) != "matrix") {
                    input[[k]] = matrix(input[[k]],length(input[[k]]),1)
                }
                if ((indices.n>1) && (nrow(input[[k]])!=(indices.n-1))) {
                    input[[k]] = matrix(input[[k]],indices.n-1,1)
                }
                    
                groups = rbind(groups,input[[k]])

            } else {
                stop("Unexpected input information type, element ",
                     as.character(k),".")
            }

            ## Advance to next input element.
            k = k+1
        }

        return (list(locations=locations,indices=indices,
                     groups=groups,group.index=group.index))
    }
    ## End of helper.

    stopifnot(is.matrix(locations))
    loc.n = nrow(locations)
    loc.dim = ncol(locations)
    stopifnot(loc.dim>=2, loc.dim<=3)

    prefix = inla.fmesher.make.prefix(dir, prefix)

    segm.bnd.input = (handle.segment.input(boundary,TRUE,nrow(locations),0))
    locations = rbind(locations,segm.bnd.input$locations)
    segm.int.input = (handle.segment.input(interior,FALSE,nrow(locations),
                                           segm.bnd.input$group.index))
    locations = rbind(locations,segm.int.input$locations)

    loc.n = nrow(locations)
    loc.dim = ncol(locations)

    ## Map locations to nodes, avoiding near-duplicates.
    node.coord = matrix(nrow=loc.n, ncol=loc.dim)
    map.loc.to.node = rep(0, nrow=loc.n)
    excluded = c()
    loc.i = 1
    node.i.max = 1
    node.coord[node.i.max,] = locations[loc.i,]
    map.loc.to.node[[loc.i]] = node.i.max
    for (loc.i in 2:loc.n) {
        loc.to.node.dist = 
            sqrt(rowSums((as.matrix(rep(1, node.i.max)) %*% locations[loc.i,] -
                          node.coord[1:node.i.max,])^2))
        if (min(loc.to.node.dist) > min.input.distance) {
            node.i.max = node.i.max+1
            node.coord[node.i.max,] = locations[loc.i,]
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
            sqrt(rowSums((as.matrix(rep(1, node.i.max)) %*% locations[loc.i,] -
                          node.coord)^2))
        node.i = which.min(loc.to.node.dist);
        map.loc.to.node[loc.i] = node.i
    }

    s0 = node.coord
    idx = map.loc.to.node

    ## The default name for the input locations in fmesher is s0:
    loc.file = inla.write.fmesher.file(s0,
            filename = paste(prefix, "s0", sep=""))

    ## additional arguments
    all.args = ""
    if (!is.null(boundary)) {
        bnd.idx = (matrix(idx[segm.bnd.input$indices],
                          nrow(segm.bnd.input$indices),
                          ncol(segm.bnd.input$indices)))
        bnd.idx[is.na(bnd.idx)] = 0L
        storage.mode(bnd.idx) = "integer"
        bnd.grp = segm.bnd.input$groups
        storage.mode(bnd.grp) = "integer"
        
        loc.file = inla.write.fmesher.file(bnd.idx-1L,
            filename = paste(prefix, "segm.bnd0", sep=""))
        loc.file = inla.write.fmesher.file(bnd.grp,
            filename = paste(prefix, "segm.bndgrp0", sep=""))
        all.args = paste(all.args," --boundary=segm.bnd0")
        all.args = paste(all.args," --boundarygrp=segm.bndgrp0")
    }
    if (!is.null(interior)) {
        int.idx = (matrix(idx[segm.int.input$indices],
                          nrow(segm.int.input$indices),
                          ncol(segm.int.input$indices)))
        int.idx[is.na(int.idx)] = 0L
        storage.mode(int.idx) = "integer"
        int.grp = segm.int.input$groups
        storage.mode(int.grp) = "integer"
        
        loc.file = inla.write.fmesher.file(int.idx-1L,
            filename = paste(prefix, "segm.int0", sep=""))
        loc.file = inla.write.fmesher.file(int.grp,
            filename = paste(prefix, "segm.intgrp0", sep=""))
        all.args = paste(all.args," --interior=segm.int0")
        all.args = paste(all.args," --interiorgrp=segm.intgrp0")
    }
    if (!is.null(cet)) {
        all.args = paste(all.args," --cet=", cet[1],",", cet[2], sep="")
    }
    if (!is.null(rcdt)) {
        all.args = paste(all.args," --rcdt=", rcdt[1],",", rcdt[2],",", rcdt[3], sep="")
    }
    if (!is.null(x11.delay)) {
        all.args = paste(all.args," --x11=", x11.delay, sep="")
    }
    all.args = paste(all.args, inla.getOption("fmesher.arg"))

    echoc = inla.fmesher.call(all.args=all.args, prefix=prefix)
    
    mesh = (list(tv = 1L+inla.read.fmesher.file(paste(prefix, "tv", sep="")),
                 s = inla.read.fmesher.file(paste(prefix, "s", sep=""))))
    
    bnd.info = (list(idx = 1L+inla.read.fmesher.file(paste(prefix, "segm.bnd", sep="")),
                     grp = inla.read.fmesher.file(paste(prefix, "segm.bnd.grp", sep=""))
                     ))
    if (!is.null(interior)) {
        int.info  = (list(idx = 1L+inla.read.fmesher.file(paste(prefix, "segm.int", sep="")),
                          grp = inla.read.fmesher.file(paste(prefix, "segm.int.grp", sep=""))
                          ))
    } else {
        int.info = NULL
    }
    mesh$segm = list(bnd=bnd.info,int=int.info)
    
    m = list(prefix = prefix, mesh = mesh, locations.idx = idx,
            call = match.call(expand.dots=TRUE))
    class(m) = "inla.fmesher.mesh"
    return (m)
}

`plot.inla.fmesher.mesh` = function(m, color = "green", size = 2, lwd=2, ...)
{
    ## a simple function that plots the mesh from inla.fmesher.mesh()

    if (length(color) == 1)
        color = rep(color, dim(m$mesh$s)[1])

    require(rgl)
    open3d()
    view3d(0, 0, fov=0)
    rgl.points(m$mesh$s[m$locations.idx, ], size=2*size, lwd=lwd, color = "blue", ...)
    plot.inla.trimesh(m$mesh$tv, m$mesh$s, color = color, size=size, lwd=lwd, ...)

    return (invisible())
}
