## Internal: internal.update.crs
## Internal: inla.internal.sp2segment.join
## Internal: inla.crs.transform.oblique
## Internal: inla.rotmat3213
## Internal: inla.rotmat3123
## Internal: inla.crs.graticule
## Internal: inla.crs.tissot
## Internal: internal.clip
## Internal: inla.crs.bounds
## Internal: inla.crs.bounds.check
##
## S3methods; also export some methods explicitly
## Export: inla.sp2segment
## Export: as.inla.mesh.segment
## Export: as.inla.mesh.segment!Polygon
## Export: as.inla.mesh.segment!Polygons
## Export: as.inla.mesh.segment!SpatialPolygons
## Export: as.inla.mesh.segment!SpatialPolygonsDataFrame
## Export: as.inla.mesh.segment!Line
## Export: as.inla.mesh.segment!Lines
## Export: as.inla.mesh.segment!SpatialLines
## Export: as.inla.mesh.segment!SpatialLinesDataFrame
## Export: as.inla.mesh.segment!SpatialPoints
## Export: as.inla.mesh.segment!SpatialPointsDataFrame
## Export: inla.CRS
## Export: inla.CRSargs
## Export: inla.as.list.CRSargs
## Export: inla.as.CRSargs.list
## Export: inla.as.list.CRS
## Export: inla.as.CRS.list
## Export: inla.identical.CRS
## Export: inla.spTransform
## Export: inla.spTransform!SpatialPoints
## Export: inla.spTransform!inla.mesh
## Export: inla.spTransform!inla.mesh.segment
## Export: inla.spTransform!inla.mesh.lattice
## Export: inla.spTransform!default
## Export: plot!CRS plot!inla.CRS


internal.clip <- function(bounds, coords, eps=0.05) {
  ## Clip 2D coordinate matrix of polylines and generate a list of Line objects
  ## bounds is from inla.crs.bounds
  ## This implementation only removes "long" line segments.
  thelines <- list()
  ## Rudimentary cutting:
  toolong <-
    which(c(TRUE,
            (diff(coords[,1]) / diff(bounds$xlim))^2
            + (diff(coords[,2]) / diff(bounds$ylim))^2
            > eps^2,
            TRUE))
  start <- toolong[-length(toolong)]
  ending <- toolong[-1] - 1
  for (i in seq_along(start)) {
    if (start[i] < ending[i]) {
      thelines <-
        c(thelines,
          list(sp::Line(coords[start[i]:ending[i], 1:2])))
    }
  }
  thelines
}


inla.crs.graticule <- function(x, by=c(15, 15, 45), add=FALSE, do.plot=TRUE,
                               eps=0.05, ...)
{
  ## Graticule
  if (is.null(by)) {
    return(invisible(list()))
  }
  if (length(by) < 2) {
    by <- by[c(1,1,1)]
  } else if (length(by) < 3) {
    by <- by[c(1,2,1)]
  }
  n <- c(floor(180/by[1]), ceiling(90/by[2])-1, floor(180/by[3]))
  bounds <- inla.crs.bounds(x)
  if (by[1] > 0) {
    special.poles <- (by[1] != by[3]) && (by[2] > 0)
    lon <- ((1-n[1]):n[1]) * by[1]
    if (special.poles) {
      lat <- seq( -n[2]*by[2], n[2]*by[2], length=91)
    } else {
      lat <- seq( -90+1e-6, 90-1e-6, length=91)
    }
    meridians <- as.matrix(expand.grid(lat, lon)[,2:1])
    proj.mer.coords <- inla.spTransform(meridians, inla.CRS("longlat"), x)
    proj.mer.coords1 <- matrix(proj.mer.coords[,1], length(lat),
                               length(lon))
    proj.mer.coords2 <- matrix(proj.mer.coords[,2], length(lat),
                               length(lon))

    mer.coords <-
      unlist(lapply(seq_along(lon),
                    function(k) {
        internal.clip(bounds, cbind(proj.mer.coords1[,k,drop=FALSE],
                                    proj.mer.coords2[,k,drop=FALSE]),
                      eps=eps)
      }),
      recursive=FALSE)

    if (special.poles) {
      if (by[3] > 0) {
        lon <- ((1-n[3]):n[3]) * by[3]
        lat <- seq( -90+1e-6, -n[2]*by[2], length=ceiling((90-n[2]*by[2])/2)+1)
        meridians <- as.matrix(expand.grid(lat, lon)[,2:1])
        proj.mer.coords <- inla.spTransform(meridians, inla.CRS("longlat"), x)
        proj.mer.coords1 <- matrix(proj.mer.coords[,1], length(lat),
                                   length(lon))
        proj.mer.coords2 <- matrix(proj.mer.coords[,2], length(lat),
                                   length(lon))
        mer.coords <-
          c(mer.coords,
            unlist(lapply(seq_along(lon),
                          function(k) {
              internal.clip(bounds, cbind(proj.mer.coords1[,k,drop=FALSE],
                                          proj.mer.coords2[,k,drop=FALSE]),
                            eps=eps)
            }),
            recursive=FALSE))

        lat <- seq( n[2]*by[2], 90-1e-6, length=ceiling((90-n[2]*by[2])/2)+1)
        meridians <- as.matrix(expand.grid(lat, lon)[,2:1])
        proj.mer.coords <- inla.spTransform(meridians, inla.CRS("longlat"), x)
        proj.mer.coords1 <- matrix(proj.mer.coords[,1], length(lat),
                                   length(lon))
        proj.mer.coords2 <- matrix(proj.mer.coords[,2], length(lat),
                                   length(lon))
        mer.coords <-
          c(mer.coords,
            unlist(lapply(seq_along(lon),
                          function(k) {
              internal.clip(bounds, cbind(proj.mer.coords1[,k,drop=FALSE],
                                          proj.mer.coords2[,k,drop=FALSE]),
                            eps=eps)
            }),
            recursive=FALSE))
      }
    }

    proj.mer <-
      sp::SpatialLines(list(sp::Lines(mer.coords, ID="meridians")),
                       proj4string=CRS(inla.CRSargs(x)))
    if (do.plot) {
      args <- list(x=proj.mer, ...)
      args <- args[intersect(names(args),
                             union(names(formals(plot.default)),
                                   union(names(formals(sp:::plot.Spatial)),
                                         names(formals(sp:::plotSpatialLines)))
                                   ))]
      if (add) {
        do.call(plot, c(list(add=TRUE), args))
      } else {
        do.call(plot, args)
        add <- TRUE
      }
    }
  } else {
    proj.mer <- NULL
  }
  if (by[2] > 0) {
    lon <- seq(-180+1e-6, 180-1e-6, length=181)
    lat <- ((-n[2]):n[2]) * by[2]
    parallels <- as.matrix(expand.grid(lon, lat))
    proj.par.coords <- inla.spTransform(parallels, inla.CRS("longlat"), x)
    proj.par.coords1 <- matrix(proj.par.coords[,1], length(lon),
                               length(lat))
    proj.par.coords2 <- matrix(proj.par.coords[,2], length(lon),
                               length(lat))
    proj.par <-
      sp::SpatialLines(list(sp::Lines(
        unlist(lapply(seq_along(lat),
                      function(k) {
          internal.clip(bounds, cbind(proj.par.coords1[,k],
                                      proj.par.coords2[,k]), eps=eps)
        }),
        recursive=FALSE),
        ID="parallels")),
        proj4string=CRS(inla.CRSargs(x)))
    if (do.plot) {
      args <- list(x=proj.par, ...)
      args <- args[intersect(names(args),
                             union(names(formals(plot.default)),
                                   union(names(formals(sp:::plot.Spatial)),
                                         names(formals(sp:::plotSpatialLines)))
                                   ))]
      if (add) {
        do.call(plot, c(list(add=TRUE), args))
      } else {
        do.call(plot, args)
      }
    }
  } else {
    proj.par <- NULL
  }
  invisible(list(meridians=proj.mer, parallels=proj.par))
}

inla.crs.tissot <- function(x, by=c(30, 30, 30), add=FALSE, do.plot=TRUE,
                            eps=0.05, diff.eps=1e-2, ...)
{
  if (is.null(by)) {
    return(invisible(list()))
  }
  if (length(by) < 2) {
    by <- c(by[1], by[1], 30)
  } else if (length(by) < 3) {
    by <- c(by[1:2], 30)
  }
  bounds <- inla.crs.bounds(x)
  n <- c(floor(180/by[1]), ceiling(90/by[2])-1)
  lon <- ((1-n[1]):n[1]) * by[1]
  lat <- ((-n[2]):n[2]) * by[2]
  loc0.lon <- loc0.lat <- loc0 <-
    cbind(as.matrix(expand.grid(lat, lon)[,2:1]), 0)
  loc0.lon[,1] <- loc0.lon[,1] + diff.eps / cos(loc0.lat[,2]*pi/180)
  loc0.lat[,2] <- loc0.lat[,2] + diff.eps
  crs.longlat <- inla.CRS("longlat")

  loc1 <- inla.spTransform(loc0, crs.longlat, x)
  loc1.lon <- inla.spTransform(loc0.lon, crs.longlat, x)
  loc1.lat <- inla.spTransform(loc0.lat, crs.longlat, x)
  ok <- (rowSums(is.na(loc1)) +
         rowSums(is.na(loc1.lon)) +
         rowSums(is.na(loc1.lat)) == 0)
  loc1 <- loc1[ok,,drop=FALSE]
  loc1.lon <- loc1.lon[ok,,drop=FALSE]
  loc1.lat <- loc1.lat[ok,,drop=FALSE]

  diff.lon <- (loc1.lon - loc1)/eps
  diff.lat <- (loc1.lat - loc1)/eps

  scale <- by[3]
  theta <- seq(0, 2*pi, length=181)
  ct <- cos(theta) * scale
  st <- sin(theta) * scale

  collection <-
      sp::SpatialLines(list(sp::Lines(
        unlist(lapply(seq_len(nrow(loc1)),
                      function(k) {
          loc1.ellipse <-
            cbind(loc1[k,1] + diff.lon[k,1] * ct + diff.lat[k,1] * st,
                  loc1[k,2] + diff.lon[k,2] * ct + diff.lat[k,2] * st)
          internal.clip(bounds, loc1.ellipse, eps=eps)
        }),
        recursive=FALSE),
        ID="parallels")),
        proj4string=CRS(inla.CRSargs(x)))
  if (do.plot) {
    args <- list(x=collection, ...)
    args <- args[intersect(names(args),
                           union(names(formals(plot.default)),
                                 union(names(formals(sp:::plot.Spatial)),
                                       names(formals(sp:::plotSpatialLines)))
                                 ))]
    if (add) {
      do.call(plot, c(list(add=TRUE), args))
    } else {
      do.call(plot, args)
    }
  }

  invisible(list(tissot=collection))
}

plot.inla.CRS <- function(x, xlim=NULL, ylim=NULL,
                          outline=TRUE,
                          graticule=c(15, 15, 45),
                          tissot=c(30,30,30),
                          asp=1,
                          add=FALSE,
                          eps=0.05,
                          ...)
{
  bounds <- inla.crs.bounds(x)
  if (is.null(xlim)) xlim <- bounds$xlim
  if (is.null(ylim)) ylim <- bounds$ylim
  if (!add) {
    args <- list(x=NA, type="n", xlim=xlim, ylim=ylim, asp=asp, ...)
    args <- args[intersect(names(args), names(formals(plot.default)))]
    do.call(plot, args)
  }
  ## Outline
  if (outline) {
    args <- list(x=bounds$polygon, ...)
    args <- args[intersect(names(args), union(names(formals(lines.default)),
                                              names(formals(plot.xy))))]
    do.call(lines, args)
  }
  ## Graticule
  inla.crs.graticule(x, by=graticule, add=TRUE, do.plot=TRUE, eps=eps, ...)
  ## Tissot
  inla.crs.tissot(x, by=tissot, add=TRUE, do.plot=TRUE, eps=eps, ...)
  invisible(NULL)
}

plot.CRS <- function(x, xlim=NULL, ylim=NULL,
                     outline=TRUE,
                     graticule=c(15, 15, 45),
                     tissot=c(30,30,30),
                     asp=1,
                     add=FALSE,
                     eps=0.05,
                     ...) {
  invisible(plot.inla.CRS(x, xlim=xlim, ylim=ylim,
                          outline=outline, graticule=graticule, tissot=tissot,
                          asp=asp, add=add, ...))
}



inla.CRS <- function(projargs = NA_character_, doCheckCRSArgs = TRUE,
                     args=NULL, oblique=NULL, ...) {
  halfroot <- "+a=0.7071067811865476 +b=0.7071067811865476"
  predef <- list(
    hammer = paste("+proj=hammer +ellps=sphere +units=m", halfroot),
    lambert = "+proj=cea +ellps=sphere +lat_ts=0 +units=m +a=1 +b=1",
    longlat = "+proj=longlat +ellps=sphere +a=1 +b=1",
    mollweide = paste("+proj=moll +ellps=sphere +units=m", halfroot),
    sphere = "+proj=geocent +ellps=sphere +a=1 +b=1 +units=m")
  if (is.character(projargs)) {
    if (projargs %in% names(predef)) {
      projargs <- predef[[projargs]]
    }
    x <- CRS(projargs, doCheckCRSArgs=doCheckCRSArgs)
  } else if (inherits(projargs, "CRS")) {
    x <- projargs
  } else {
    stop(paste("Unsupported projargs input class",
               paste(class(projargs), collapse=",")))
  }
  if (!is.null(args)) {
    if (typeof(args) != "list") {
      stop("'args' must be NULL or a list of name=value pairs.")
    }
    xargs <- inla.as.list.CRS(x)
    for (name in names(args)) {
      xargs[[name]] <- args[[name]]
    }
    x <- CRS(inla.as.CRSargs.list(xargs), doCheckCRSArgs=doCheckCRSArgs)
  }

  if (!is.null(oblique)) {
    stopifnot(is.vector(oblique))
    if (length(oblique) < 4) {
      oblique <- c(oblique, rep(0, 4-length(oblique)))
    }
    x <- list(crs=x, oblique=oblique)
    class(x) <- "inla.CRS"
  }
  x
}

inla.as.list.CRS <- function(x, ...) {
  inla.as.list.CRSargs(inla.CRSargs(x))
}

inla.as.CRS.list <- function(x, ...) {
  inla.CRS(args=x)
}

inla.CRSargs <- function(x, ...) {
  if (inherits(x, "inla.CRS")) {
    x <- x[["crs"]]
  }
  if (is.null(x)) {
    as.character(NA)
  } else {
    stopifnot(inla.require("rgdal"))
    rgdal::CRSargs(x)
  }
}


## CRS proj4 string for name=value pair list
inla.as.CRSargs.list <- function(x, ...) {
  paste(lapply(names(x),
               function(xx) {
    if (is.na(x[[xx]])) {
      paste("+", xx, sep="")
    } else {
      paste("+", xx, "=", x[[xx]], sep="")
    }
  }),
  collapse=" ")
}

## List of name=value pairs from CRS proj4 string
inla.as.list.CRSargs <- function(x, ...) {
  if (is.na(x)) {
    return(list())
  }
  if (!is.character(x)) {
    stop("proj4string must be of class character")
  }
  do.call(c, lapply(strsplit(x=strsplit(x=paste(" ", x, sep=""),
                                        split=" \\+")[[1]][-1],
                             split="="),
                    function(x) {
               xx <- list(x[2])
               names(xx) <- x[1]
               xx
             }))
}


inla.rotmat3213 <- function(rot)
{
  cs <- cos(rot[1])
  sn <- sin(rot[1])
  R <- matrix(c(cs, -sn, 0,
                sn, cs, 0,
                0, 0, 1), 3, 3)
  cs <- cos(rot[2])
  sn <- sin(rot[2])
  R <- R %*% matrix(c(cs, 0, sn,
                      0, 1, 0,
                      -sn, 0, cs), 3, 3)
  cs <- cos(rot[3])
  sn <- sin(rot[3])
  R <- R %*% matrix(c(1, 0, 0,
                      0, cs, -sn,
                      0, sn, cs), 3, 3)
  cs <- cos(rot[4])
  sn <- sin(rot[4])
  R <- R %*% matrix(c(cs, -sn, 0,
                      sn, cs, 0,
                      0, 0, 1), 3, 3)
  R
}

inla.rotmat3123 <- function(rot)
{
  cs <- cos(rot[4])
  sn <- sin(rot[4])
  R <- matrix(c(cs, -sn, 0,
                sn, cs, 0,
                0, 0, 1), 3, 3)
  cs <- cos(rot[3])
  sn <- sin(rot[3])
  R <- R %*% matrix(c(1, 0, 0,
                      0, cs, -sn,
                      0, sn, cs), 3, 3)
  cs <- cos(rot[2])
  sn <- sin(rot[2])
  R <- R %*% matrix(c(cs, 0, sn,
                      0, 1, 0,
                      -sn, 0, cs), 3, 3)
  cs <- cos(rot[1])
  sn <- sin(rot[1])
  R <- R %*% matrix(c(cs, -sn, 0,
                      sn, cs, 0,
                      0, 0, 1), 3, 3)
  R
}


inla.crs.transform.oblique <- function(x, oblique, to.oblique=TRUE) {
  if (to.oblique) {
    ## Transform to oblique orientation
    ## 1) Rotate -oblique[1] around (0,0,1)
    ## 2) Rotate +oblique[2] around (0,1,0)
    ## 3) Rotate -oblique[3] around (1,0,0)
    ## 3) Rotate -oblique[4] around (0,0,1)
    x %*% inla.rotmat3213(c(-1,1,-1,-1) * oblique * pi/180)
  } else {
    ## Transform back from oblique orientation
    ## 1) Rotate +oblique[4] around (0,0,1)
    ## 2) Rotate +oblique[3] around (1,0,0)
    ## 3) Rotate -oblique[2] around (0,1,0)
    ## 4) Rotate +oblique[1] around (0,0,1)
    x %*% inla.rotmat3123(c(1,-1,1,1) *oblique * pi/180)
  }
}



## +proj=longlat in (-180,180)x(-90,90)
## +proj=moll in (-2,2)x(-1,1) scaled by +a and +b, and +units
## +proj=lambert in (-pi,pi)x(-1,1) scaled by +a and +b, and +units
inla.crs.bounds <- function(crs, warn.unknown = FALSE) {
  args <- inla.as.list.CRS(crs)
  if (args[["proj"]] == "longlat") {
    bounds <- list(type="rectangle", xlim=c(-180,180), ylim=c(-90,90))
  } else if (args[["proj"]] == "cea") {
    axis <- c(pi, 1)
    if (!is.null(args[["a"]])) {
      axis[1] <- axis[1] * as.numeric(args$a)
    }
    if (!is.null(args[["b"]])) {
      axis[2] <- axis[2] * as.numeric(args$a)^0.5 * as.numeric(args$b)^0.5
    }
    ## TODO: Handle "lat_ts" and "units"
    bounds <- list(type="rectangle",
                   xlim=c(-1,1)*axis[1], ylim=c(-1,1)*axis[2])
  } else if (args[["proj"]] %in% c("moll", "hammer")) {
    axis <- c(2, 1)
    center <- c(0,0)
    if (!is.null(args[["a"]])) {
      axis[1] <- axis[1] * as.numeric(args$a) / sqrt(1/2)
      axis[2] <- axis[2] * as.numeric(args$a) / sqrt(1/2)
    }
    ## TODO: Handle "units"
    bounds <- list(type="ellipse", axis=axis, center=center,
                   xlim=center[1]+c(-1,1)*axis[1],
                   ylim=center[2]+c(-1,1)*axis[2])
  } else if (args[["proj"]] == "tmerc") {
    bounds <- list(type="rectangle", xlim=c(-Inf, Inf), ylim=c(-Inf, Inf))
  } else if (args[["proj"]] == "geocent") {
    bounds <- list(type="rectangle", xlim=c(-Inf, Inf), ylim=c(-Inf, Inf))
  } else {
    if (warn.unknown) {
      warning("Could not determine shape of transformation bounds. Using infinite rectangle.")
    }
    bounds <- list(type="rectangle", xlim=c(-Inf, Inf), ylim=c(-Inf, Inf))
  }
  if (bounds$type == "rectangle") {
    bounds$polygon <- cbind(bounds$xlim[c(1,2,2,1,1)],
                            bounds$ylim[c(1,1,2,2,1)])
  } else if (bounds$type == "ellipse") {
    theta <- seq(0, 2*pi, length=1000)
    bounds$polygon <- cbind(bounds$center[1] + bounds$axis[1]*cos(theta),
                            bounds$center[2] + bounds$axis[2]*sin(theta))
  } else {
    stop("Unknown transformation type. This should not happen.")
  }
  bounds
}
## TRUE/FALSE for points inside/outside projection domain.
inla.crs.bounds.check <- function(x, bounds) {
  inla.require.inherits(x, "matrix")
  if (all(is.finite(bounds$xlim)) && all(is.finite(bounds$ylim))) {
    (sp::point.in.polygon(x[,1], x[,2],
                          bounds$polygon[,1], bounds$polygon[,2])
      > 0)
  } else {
    ok <- rep(TRUE, nrow(x))
  }
}



internal.update.crs <- function(crs, newcrs, mismatch.allowed) {
  if (is.null(crs)) {
    newcrs
  } else {
    if (!mismatch.allowed && !identical(crs, newcrs)) {
      show(crs)
      show(newcrs)
      stop("CRS information mismatch.")
    }
    crs
  }
}


inla.identical.CRS <- function(crs0, crs1, crsonly=FALSE) {
  if (!crsonly) {
    identical(crs0, crs1)
  } else {
    if (inherits(crs0, "inla.CRS")) {
      crs0 <- crs0$crs
    }
    if (inherits(crs1, "inla.CRS")) {
      crs1 <- crs1$crs
    }
    identical(crs0, crs1)
  }
}

## Low level transformation of raw coordinates.
inla.spTransform.default <- function(x, crs0, crs1, passthrough=FALSE, ...) {
  ok0 <- (!is.null(crs0) &&
          ((inherits(crs0, "CRS") && !is.na(inla.CRSargs(crs0))) ||
           (inherits(crs0, "inla.CRS"))))
  ok1 <- (!is.null(crs1) &&
          ((inherits(crs1, "CRS") && !is.na(inla.CRSargs(crs1))) ||
           (inherits(crs1, "inla.CRS"))))
  if (ok0 && ok1) {
    if (ncol(x) == 2) {
      x <- cbind(x, 0)
    }
    onsphere <- inla.identical.CRS(crs0, inla.CRS("sphere"), crsonly=TRUE)
    isgeocentric <- identical(inla.as.list.CRS(crs0)[["proj"]], "geocent")
    if (isgeocentric) {
      ok <- TRUE
    } else {
      bounds <- inla.crs.bounds(crs0)
      if (identical(inla.as.list.CRS(crs0)[["proj"]], "longlat")) {
        ## Wrap longitudes to [-180,180]
        needswrap <- (x[,1] < -180) | (x[,1] > 180)
        if (any(needswrap)) {
          x[needswrap,1] <- ((x[needswrap,1] + 180) %% 360) - 180
        }
      }
      ok <- inla.crs.bounds.check(x, bounds)
      if (!all(ok)) {
        xx <- x
      }
    }
    if (inherits(crs0, "inla.CRS")) {
      if (!onsphere) {
        x <- spTransform(SpatialPoints(x[ok,,drop=FALSE], proj4string=crs0$crs),
                         inla.CRS("sphere"))
      }
      if (!is.null(crs0$oblique)) {
        x <- SpatialPoints(inla.crs.transform.oblique(coordinates(x),
                                                      crs0$oblique,
                                                      to.oblique=FALSE),
                           proj4string=inla.CRS("sphere"))
      }
      onshpere <- TRUE
    } else {
      x <- SpatialPoints(x[ok,,drop=FALSE], proj4string=crs0)
    }
    if (inherits(crs1, "inla.CRS")) {
      if (!onsphere) {
        x <- spTransform(x, inla.CRS("sphere"))
      }
      if (!is.null(crs1$oblique)) {
        x <- SpatialPoints(inla.crs.transform.oblique(coordinates(x),
                                                      crs1$oblique,
                                                      to.oblique=TRUE),
                           proj4string=inla.CRS("sphere"))
      }
      x <- spTransform(x, crs1$crs)
    } else {
      x <- spTransform(x, crs1)
    }
    if (!all(ok)) {
      xx[ok,] <- coordinates(x)
      xx[!ok,] <- NA
      x <- xx
    }
  } else if (!passthrough) {
    if (!ok0) {
      stop("'crs0' is an invalid coordinate reference object.")
    }
    if (!ok1) {
      stop("'crs1' is an invalid coordinate reference object.")
    }
  }
  if (is.matrix(x)) {
    invisible(x)
  } else {
    invisible(coordinates(x))
  }
}

inla.spTransform.SpatialPoints <- function(x, CRSobj, passthrough=FALSE, ...) {
  ok0 <- !is.na(proj4string(x))
  ok1 <- (!missing(CRSobj) && !is.null(CRSobj) &&
          (inherits(CRSobj, "CRS") && !is.na(inla.CRSargs(CRSobj))))
  if (ok0 && ok1) {
    invisible(SpatialPoints(inla.spTransform(coordinates(x),
                                             CRS(proj4string(x)),
                                             CRSobj),
                            proj4string=CRSobj))
  } else if (ok1) { ## Know: !ok0 && ok1
    if (!passthrough) {
      stop("Invalid origin CRS for SpatialPoints")
    }
    invisible(SpatialPoints(coordinates(x), proj4string=CRSobj))
  } else { ## Know: (ok0 || !ok0) && !ok1
    if (!passthrough) {
      stop("Invalid target CRS for SpatialPoints")
    }
    invisible(SpatialPoints(coordinates(x), proj4string=inla.CRS()))
  }
}

inla.spTransform.inla.mesh.lattice <- function(x, CRSobj, passthrough=FALSE, ...) {
  x$segm <- inla.spTransform(x$segm, CRSobj, passthrough=passthrough)
  x$loc <- inla.spTransform(x$loc, x$crs, CRSobj, passthrough=passthrough)
  x$crs <- CRSobj
  invisible(x)
}

inla.spTransform.inla.mesh.segment <- function(x, CRSobj, passthrough=FALSE, ...) {
  x$loc <- inla.spTransform(x$loc, x$crs, CRSobj, passthrough=passthrough)
  x$crs <- CRSobj
  invisible(x)
}

inla.spTransform.inla.mesh <- function(x, CRSobj, passthrough=FALSE, ...) {
  x$loc <- inla.spTransform(x$loc, x$crs, CRSobj, passthrough=passthrough)
  args <- inla.as.list.CRS(CRSobj)
  if (identical(args[["proj"]], "geocent")) {
    x$manifold <- "S2"
  } else {
    x$manifold <- "R2"
  }
  x$crs <- CRSobj
  invisible(x)
}

inla.spTransform <- function(x, ...) {
  UseMethod("inla.spTransform")
}


## Input: list of segments, all closed polygons.
inla.internal.sp2segment.join <- function(inp, grp=NULL, closed=TRUE) {
    crs <- NULL
    if (length(inp) > 0) {
      out.loc = matrix(0,0,ncol(inp[[1]]$loc))
      for (k in seq_along(inp)) {
        crs <- internal.update.crs(crs, inp[[k]]$crs, mismatch.allowed=FALSE)
      }
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
    inla.mesh.segment(loc=out.loc, idx=out.idx, grp=out.grp, is.bnd=FALSE,
                      crs=crs)
}


as.inla.mesh.segment <-
    function(sp, ...)
{
    UseMethod("as.inla.mesh.segment")
}

inla.sp2segment <-
    function(sp, ...)
{
    UseMethod("as.inla.mesh.segment")
}



as.inla.mesh.segment.SpatialPoints <-
    function (sp, reverse=FALSE, grp = NULL, is.bnd=TRUE, ...)
  {
    crs <- CRS(proj4string(sp))
    loc <- coordinates(sp)

    n = dim(loc)[1L]
    if (reverse) {
        idx <- seq(n, 1L, length=n)
    } else {
        idx <- seq_len(n)
    }
    inla.mesh.segment(loc = loc, idx = idx, grp = grp, is.bnd = is.bnd,
                      crs=crs)
  }

as.inla.mesh.segment.SpatialPointsDataFrame <-
    function (sp, ...)
{
    as.inla.mesh.segment.SpatialLines(sp, ...)
}



as.inla.mesh.segment.Line <-
    function(sp, reverse=FALSE, crs=NULL, ...)
{
    loc = sp@coords
    n = dim(loc)[1L]
    if (reverse) {
        idx <- seq(n, 1L, length=n)
    } else {
        idx <- seq_len(n)
    }
    inla.mesh.segment(loc = loc, idx = idx, is.bnd = FALSE, crs=crs)
}

as.inla.mesh.segment.Lines <-
    function (sp, join = TRUE, crs=NULL, ...)
{
    segm <- as.list(lapply(sp@Lines,
                           function(x) as.inla.mesh.segment(x, crs=crs, ...)))
    if (join)
        segm = inla.internal.sp2segment.join(segm, grp = NULL, closed=FALSE)
    segm
}

as.inla.mesh.segment.SpatialLines <-
    function (sp, join = TRUE, grp = NULL, ...)
  {
    crs <- CRS(proj4string(sp))
    segm = list()
    for (k in 1:length(sp@lines))
      segm[[k]] = as.inla.mesh.segment(sp@lines[[k]], join = TRUE,
                                       crs = crs, ...)
    if (join) {
        if (missing(grp)) {
            grp = 1:length(segm)
        }
        segm = inla.internal.sp2segment.join(segm, grp = grp, closed=FALSE)
    }
    segm
}

as.inla.mesh.segment.SpatialLinesDataFrame <-
    function (sp, ...)
{
    as.inla.mesh.segment.SpatialLines(sp, ...)
}

as.inla.mesh.segment.SpatialPolygons <-
    function(sp, join=TRUE, grp=NULL, ...)
{
    crs <- CRS(proj4string(sp))
    segm = list()
    for (k in 1:length(sp@polygons))
        segm[[k]] = as.inla.mesh.segment(sp@polygons[[k]], join=TRUE, crs=crs)
    if (join) {
        if (missing(grp)) {
            grp = 1:length(segm)
        }
        segm = inla.internal.sp2segment.join(segm, grp=grp)
    }
    segm
}

as.inla.mesh.segment.SpatialPolygonsDataFrame <-
    function(sp, ...)
{
    as.inla.mesh.segment.SpatialPolygons(sp, ...)
}

as.inla.mesh.segment.Polygons <-
    function(sp, join=TRUE, crs=NULL, ...)
{
    segm = as.list(lapply(sp@Polygons,
                          function (x) as.inla.mesh.segment(x, crs=crs)))
    if (join)
        segm = inla.internal.sp2segment.join(segm, grp=NULL)
    segm
}

as.inla.mesh.segment.Polygon <-
    function(sp, crs=NULL, ...)
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
    inla.mesh.segment(loc=loc, idx=idx, is.bnd=TRUE, crs=crs)
}
