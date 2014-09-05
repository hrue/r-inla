## Internal: inla.internal.sp2segment.join
##
## S3methods; also export some methods explicitly
## Export: inla.sp2segment
## Export: inla.sp2segment.Polygon inla.sp2segment.Polygons
## Export: inla.sp2segment.SpatialPolygons
## Export: inla.sp2segment.SpatialPolygonsDataFrame
## Export: inla.sp2segment!Polygon inla.sp2segment!Polygons
## Export: inla.sp2segment!SpatialPolygons
## Export: inla.sp2segment!SpatialPolygonsDataFrame
## Export: inla.sp2segment.Line inla.sp2segment.Lines
## Export: inla.sp2segment.SpatialLines
## Export: inla.sp2segment.SpatialLinesDataFrame
## Export: inla.sp2segment!Line inla.sp2segment!Lines
## Export: inla.sp2segment!SpatialLines
## Export: inla.sp2segment!SpatialLinesDataFrame
## Export: as.inla.mesh.segment
## Export: as.inla.mesh.segment!Polygon as.inla.mesh.segment!Polygons
## Export: as.inla.mesh.segment!SpatialPolygons
## Export: as.inla.mesh.segment!SpatialPolygonsDataFrame
## Export: as.inla.mesh.segment!Line as.inla.mesh.segment!Lines
## Export: as.inla.mesh.segment!SpatialLines
## Export: as.inla.mesh.segment!SpatialLinesDataFrame


## Input: list of segments, all closed polygons.
inla.internal.sp2segment.join <- function(inp, grp=NULL, closed=TRUE) {
    if (length(inp) > 0) {
        out.loc = matrix(0,0,ncol(inp[[1]]$loc))
    } else {
        out.loc = matrix(0,0,2)
    }
    out.idx = matrix(0L,0,2)
    if (is.null(grp)) {
        out.grp = NULL
    } else {
        out.grp = integer(0)
    }
    for (k in seq_along(inp)) {
        inp.loc = inp[[k]]$loc
        inp.idx = inp[[k]]$idx
        inp.grp = inp[[k]]$grp
        offset = nrow(out.loc)
        n = nrow(as.matrix(inp.idx))
        if (closed) {
            if (!is.null(grp) && is.null(inp.grp)) {
                inp.grp = rep(grp[k], n)
            }
            if (ncol(as.matrix(inp.idx))==1) {
                inp.idx = cbind(inp.idx, inp.idx[c(2:n,1)])
            }
        } else {
            if (!is.null(grp) && is.null(inp.grp)) {
                inp.grp = rep(grp[k], n-1)
            }
            if (ncol(as.matrix(inp.idx))==1) {
                inp.idx = cbind(inp.idx[-n], inp.idx[-1])
            }
        }
        out.loc = rbind(out.loc, inp.loc)
        out.idx = rbind(out.idx, inp.idx+offset)
        if (!is.null(grp)) {
            out.grp = c(out.grp, inp.grp)
        }
    }
    out = inla.mesh.segment(loc=out.loc,idx=out.idx,grp=out.grp,is.bnd=FALSE)
}


as.inla.mesh.segment <-
    function(sp, ...)
{
    UseMethod("as.inla.mesh.segment")
}

as.inla.mesh.segment.Line <-
    function(sp, reverse=FALSE, ...)
{
    loc = sp@coords
    n = dim(loc)[1L]
    if (reverse) {
        idx <- seq(n, 1L, length=n)
    } else {
        idx <- seq_len(n)
    }
    return(inla.mesh.segment(loc = loc, idx = idx, is.bnd = FALSE))
}

as.inla.mesh.segment.Lines <-
    function (sp, join = TRUE, ...)
{
    segm <- as.list(lapply(sp@Lines,
                           function(x) as.inla.mesh.segment(x, ...)))
    if (join)
        segm = inla.internal.sp2segment.join(segm, grp = NULL, closed=FALSE)
    return(segm)
}

as.inla.mesh.segment.SpatialLines <-
    function (sp, join = TRUE, grp = NULL, ...)
{
    segm = list()
    for (k in 1:length(sp@lines))
        segm[[k]] = as.inla.mesh.segment(sp@lines[[k]], join = TRUE, ...)
    if (join) {
        if (missing(grp)) {
            grp = 1:length(segm)
        }
        segm = inla.internal.sp2segment.join(segm, grp = grp, closed=FALSE)
    }
    return(segm)
}

as.inla.mesh.segment.SpatialLinesDataFrame <-
    function (sp, ...)
{
    as.inla.mesh.segment.SpatialLines(sp, ...)
}

as.inla.mesh.segment.SpatialPolygons <-
    function(sp, join=TRUE, grp=NULL, ...)
{
    segm = list()
    for (k in 1:length(sp@polygons))
        segm[[k]] = as.inla.mesh.segment(sp@polygons[[k]], join=TRUE)
    if (join) {
        if (missing(grp)) {
            grp = 1:length(segm)
        }
        segm = inla.internal.sp2segment.join(segm, grp=grp)
    }
    return(segm)
}

as.inla.mesh.segment.SpatialPolygonsDataFrame <-
    function(sp, ...)
{
    as.inla.mesh.segment.SpatialPolygons(sp, ...)
}

as.inla.mesh.segment.Polygons <-
    function(sp, join=TRUE, ...)
{
    segm = as.list(lapply(sp@Polygons, function (x) as.inla.mesh.segment(x)))
    if (join)
        segm = inla.internal.sp2segment.join(segm, grp=NULL)
    return(segm)
}

as.inla.mesh.segment.Polygon <-
    function(sp, ...)
{
    loc = sp@coords[-dim(sp@coords)[1L],,drop=FALSE]
    n = dim(loc)[1L]
    if (sp@hole)
        if (sp@ringDir==1)
            idx = c(1L:n,1L)
        else
            idx = c(1L,seq(n,1L,length.out=n))
    else
        if (sp@ringDir==1)
            idx = c(1L,seq(n,1L,length.out=n))
        else
            idx = c(1L:n,1L)
    return(inla.mesh.segment(loc=loc, idx=idx, is.bnd=TRUE))
}








inla.sp2segment.Line <-
    function (sp, reverse=FALSE, ...)
{
    as.inla.mesh.segment.Line(sp, reverse=reverse, ...)
}

inla.sp2segment.Lines <-
    function (sp, join = TRUE, ...)
{
    as.inla.mesh.segment.Lines(sp, join=join, ...)
}

inla.sp2segment.SpatialLines <-
    function (sp, join = TRUE, grp = NULL, ...)
{
    as.inla.mesh.segment.SpatialLines(sp, join=join, grp=grp, ...)
}

inla.sp2segment.SpatialLinesDataFrame <-
    function (sp, ...)
{
    as.inla.mesh.segment.SpatialLines(sp, ...)
}

inla.sp2segment.SpatialPolygons <-
    function(sp, join=TRUE, grp=NULL, ...)
{
    as.inla.mesh.segment.SpatialPolygons(sp, join=join, grp=grp, ...)
}

inla.sp2segment.SpatialPolygonsDataFrame <-
    function(sp, ...)
{
    as.inla.mesh.segment.SpatialPolygons(sp, ...)
}

inla.sp2segment.Polygons <-
    function(sp, join=TRUE, ...)
{
    as.inla.mesh.segment.Polygons(sp, join=join, ...)
}

inla.sp2segment.Polygon <-
    function(sp, ...)
{
    as.inla.mesh.segment.Polygons(s, ...)
}
