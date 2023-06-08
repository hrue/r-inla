## Internal: internal.update.crs
## Internal: inla.internal.sp2segment.join
## Internal: inla.crs.graticule
## Internal: inla.crs.tissot
## Internal: internal.clip
## Internal: inla.crs.bounds
## Internal: inla.crs.bounds.check




#' PROJ6 detection
#'
#' Detect whether PROJ6 is available for INLA. Deprecated and always returns `TRUE`..
#'
#' @aliases inla.has_PROJ6 inla.not_for_PROJ6 inla.not_for_PROJ4
#' inla.fallback_PROJ6 inla.requires_PROJ6
#' @details `inla.has_PROJ6` is called to check if PROJ6&GDAL3 are available.
#' @return For `inla.has_PROJ6`, always returns `TRUE`. Previously: logical; `TRUE` if PROJ6 is available,
#' `FALSE` otherwise
#' @examples
#'
#' inla.has_PROJ6()
#' @export inla.has_PROJ6
#' @rdname inla.has_PROJ6
inla.has_PROJ6 <- function() {
    lifecycle::deprecate_warn(
        "23.01.01",
        "inla.has_PROJ6()",
        details = "PROJ4 is no longer supported, and inla.has_PROJ6() always returns TRUE"
    )
    TRUE
}

#' @describeIn inla.has_PROJ6 Called to warn about using old PROJ4
#' features even though PROJ6 is available
#' @param fun The name of the calling function

inla.not_for_PROJ6 <- function(fun) {
    lifecycle::deprecate_soft(when = "2023.04.19",
                              what = "inla.not_for_PROJ6()",
                              details = "Should no longer be used.")
    if (inla.has_PROJ6()) {
        stack <- sys.calls()
        stack <- lapply(as.list(stack), function(x) as.character(deparse(x)))
        warning(paste0(
            "'",
            fun,
            "' should not be used with PROJ6\n",
            "Call stack for developer debugging:\n",
            paste0(lapply(stack, function(x) substr(x, 1, 70)),
                collapse = "\n"
            )
        ))
    }
}

#' @describeIn inla.has_PROJ6 Called to give an error when
#' calling methods that are only available for PROJ6

inla.not_for_PROJ4 <- function(fun) {
    lifecycle::deprecate_soft(when = "2023.04.19",
                              what = "inla.not_for_PROJ4()",
                              details = "Should no longer be used.")
    if (!inla.has_PROJ6()) {
        stop(paste0(
            "'",
            fun,
            "' is not supported for PROJ4"
        ))
    }
}

#' @describeIn inla.has_PROJ6 Called to warn about falling back
#' to using old PROJ4 methods when a RPOJ6 method hasn't been implemented

inla.fallback_PROJ6 <- function(fun) {
    lifecycle::deprecate_soft(when = "2023.04.19",
                              what = "inla.fallback_PROJ6()",
                              details = "Should no longer be used.")
    if (inla.has_PROJ6()) {
        stop(paste0(
            "'",
            fun,
            "' method for PROJ6 not implemented. PROJ4 is no longer supported."
        ))
    }
}

#' @describeIn inla.has_PROJ6 Called to give an error when PROJ6
#' is required but not available

inla.requires_PROJ6 <- function(fun) {
    lifecycle::deprecate_soft(when = "2023.04.19",
                              what = "inla.requires_PROJ6()",
                              details = "Should no longer be used.")
    if (!inla.has_PROJ6()) {
        stop(paste0(
            "'",
            fun,
            "' requires PROJ6"
        ))
    }
}






#' @title Extract CRS information
#'
#' @description
#' Wrapper for CRS(projargs) (PROJ4) and CRS(wkt) for `sp::Spatial`
#' objects.
#'
#' This function is a convenience method to workaround PROJ4/PROJ6 differences,
#' and the lack of a crs extraction method for Spatial objects.
#'
#' @param x A `sp::Spatial` object
#' @return A `CRS` object, or NULL if no valid CRS identified
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @examples
#'
#' \dontrun{
#' if (interactive()) {
#'     s <- sp::SpatialPoints(matrix(1:6, 3, 2), proj4string = inla.CRS("sphere"))
#'     inla.sp_get_crs(s)
#' }
#' }
#'
#' @export inla.sp_get_crs
inla.sp_get_crs <- function(x) {
    if (inla.getOption("fmesher.evolution") >= 2L) {
        lifecycle::deprecate_soft(when = "2023.06.06",
                                  what = "inla.sp_get_crs()",
                                  with = "inlabru::fm_CRS()")
    }
    if (is.null(x)) {
        return(NULL)
    }
    fm_CRS(x)
}





internal.clip <- function(bounds, coords, eps = 0.05) {
    ## Clip 2D coordinate matrix of polylines and generate a list of Line objects
    ## bounds is from inla.crs.bounds
    ## This implementation only removes "long" line segments.
    thelines <- list()
    ## Rudimentary cutting:
    toolong <-
        which(c(
            TRUE,
            (diff(coords[, 1]) / diff(bounds$xlim))^2
                + (diff(coords[, 2]) / diff(bounds$ylim))^2
            > eps^2,
            TRUE
        ))
    start <- toolong[-length(toolong)]
    ending <- toolong[-1] - 1
    for (i in seq_along(start)) {
        if (start[i] < ending[i]) {
            thelines <-
                c(
                    thelines,
                    list(sp::Line(coords[start[i]:ending[i], 1:2]))
                )
        }
    }
    thelines
}




#' Handling CRS/WKT
#'
#' Get and set CRS object or WKT string properties.
#'
#'
#' @aliases inla.wkt_is_geocent inla.crs_is_geocent
#' inla.wkt_get_ellipsoid_radius inla.crs_get_ellipsoid_radius
#' inla.wkt_set_ellipsoid_radius inla.crs_set_ellipsoid_radius
#' inla.wkt_unit_params crs_wkt inla.wkt_get_lengthunit inla.wkt_set_lengthunit
#' inla.crs_get_wkt inla.crs_get_lengthunit inla.crs_set_lengthunit
#' @param wkt A WKT2 character string
#' @param crs A `sp::CRS` or `inla.CRS` object
#' @param unit character, name of a unit. Supported names are "metre",
#' "kilometre", and the aliases "meter", "m", International metre",
#' "kilometer", and "km", as defined by `inla.wkt_unit_params` or the
#' `params` argument. (For legacy PROJ4 use, only "m" and "km" are
#' supported)
#' @param params Length unit definitions, in the list format produced by
#' `inla.wkt_unit_params()`, Default: NULL, which invokes
#' `inla.wkt_unit_params()`
#' @return For `inla.wkt_unit_params`, a list of named unit definitions
#'
#' For `inla.wkt_get_lengthunit`, a list of length units used in the wkt
#' string, excluding the ellipsoid radius unit.
#'
#' For `inla.wkt_set_lengthunit`, a WKT2 string with altered length units.
#' Note that the length unit for the ellipsoid radius is unchanged.
#'
#' For `inla.crs_get_wkt`, WKT2 string.
#'
#' For `inla.crs_get_lengthunit`, a list of length units used in the wkt
#' string, excluding the ellipsoid radius unit. (For legacy PROJ4 code, the raw
#' units from the proj4string are returned, if present.)
#'
#' For `inla.crs_set_lengthunit`, a `sp::CRS` object with altered
#' length units. Note that the length unit for the ellipsoid radius is
#' unchanged.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.sp_get_crs()]
#' @examples
#'
#' \dontrun{
#'   c1 <- inla.CRS("globe")
#'   inla.crs_get_lengthunit(c1)
#'   c2 <- inla.crs_set_lengthunit(c1, "km")
#'   inla.crs_get_lengthunit(c2)
#' }
#'
#' @rdname crs_wkt
#' @export inla.wkt_is_geocent
inla.wkt_is_geocent <- function(wkt) {
    if (inla.getOption("fmesher.evolution") >= 2L) {
        lifecycle::deprecate_soft(when = "2023.06.06",
                                  what = "inla.wkt_is_geocent()",
                                  with = "inlabru::fm_wkt_is_geocent()")
    }
    inlabru::fm_wkt_is_geocent(wkt)
}

#' @export
#' @rdname crs_wkt

inla.crs_is_geocent <- function(crs) {
    if (inla.getOption("fmesher.evolution") >= 2L) {
        lifecycle::deprecate_soft(when = "2023.06.06",
                                  what = "inla.crs_is_geocent()",
                                  with = "inlabru::fm_crs_is_geocent()")
    }
    inlabru::fm_crs_is_geocent(crs)
}


#' @rdname crs_wkt
#' @export

inla.wkt_get_ellipsoid_radius <- function(wkt) {
    if (inla.getOption("fmesher.evolution") >= 2L) {
        lifecycle::deprecate_soft(when = "2023.06.06",
                                  what = "inla.wkt_get_ellipsoid_radius()",
                                  with = "inlabru::fm_ellipsoid_radius()")
    }
    inlabru::fm_ellipsoid_radius(wkt)
}

#' @rdname crs_wkt
#' @export

inla.crs_get_ellipsoid_radius <- function(crs) {
    if (inla.getOption("fmesher.evolution") >= 2L) {
        lifecycle::deprecate_soft(when = "2023.06.06",
                                  what = "inla.crs_get_ellipsoid_radius()",
                                  with = "inlabru::fm_ellipsoid_radius()")
    }
    inlabru::fm_ellipsoid_radius(crs)
}


#' @rdname crs_wkt
#' @param radius numeric
#' @export

inla.wkt_set_ellipsoid_radius <- function(wkt, radius) {
    if (inla.getOption("fmesher.evolution") >= 2L) {
        lifecycle::deprecate_soft(when = "2023.06.06",
                                  what = "inla.wkt_set_ellipsoid_radius()",
                                  with = "`inlabru::fm_ellipsoid_radius<-`()")
    }
    inlabru::fm_wkt_get_ellisoid_radius(wkt) <- radius
    wkt
}

#' @rdname crs_wkt
#' @export

inla.crs_set_ellipsoid_radius <- function(crs, radius) {
    if (inla.getOption("fmesher.evolution") >= 2L) {
        lifecycle::deprecate_soft(when = "2023.06.06",
                                  what = "inla.crs_set_ellipsoid_radius()",
                                  with = "`inlabru::fm_ellipsoid_radius<-`()")
    }
    inlabru::fm_wkt_get_ellisoid_radius(crs) <- radius
    crs
}





inla.crs.graticule <- function(x, by = c(15, 15, 45), add = FALSE, do.plot = TRUE,
                               eps = 0.05, ...) {
    inla.fallback_PROJ6("inla.crs.graticule")

    ## Graticule
    if (is.null(by)) {
        return(invisible(list()))
    }
    if (length(by) < 2) {
        by <- by[c(1, 1, 1)]
    } else if (length(by) < 3) {
        by <- by[c(1, 2, 1)]
    }
    n <- c(floor(180 / by[1]), ceiling(90 / by[2]) - 1, floor(180 / by[3]))
    bounds <- inla.crs.bounds(x)
    if (by[1] > 0) {
        special.poles <- (by[1] != by[3]) && (by[2] > 0)
        lon <- ((1 - n[1]):n[1]) * by[1]
        if (special.poles) {
            lat <- seq(-n[2] * by[2], n[2] * by[2], length = 91)
        } else {
            lat <- seq(-90 + 1e-6, 90 - 1e-6, length = 91)
        }
        meridians <- as.matrix(expand.grid(lat, lon)[, 2:1])
        proj.mer.coords <- fm_transform(meridians, crs0 = inla.CRS("longlat"), crs = x)
        proj.mer.coords1 <- matrix(
            proj.mer.coords[, 1], length(lat),
            length(lon)
        )
        proj.mer.coords2 <- matrix(
            proj.mer.coords[, 2], length(lat),
            length(lon)
        )

        mer.coords <-
            unlist(lapply(
                seq_along(lon),
                function(k) {
                    internal.clip(bounds, cbind(
                        proj.mer.coords1[, k, drop = FALSE],
                        proj.mer.coords2[, k, drop = FALSE]
                    ),
                    eps = eps
                    )
                }
            ),
            recursive = FALSE
            )

        if (special.poles) {
            if (by[3] > 0) {
                lon <- ((1 - n[3]):n[3]) * by[3]
                lat <- seq(-90 + 1e-6, -n[2] * by[2], length = ceiling((90 - n[2] * by[2]) / 2) + 1)
                meridians <- as.matrix(expand.grid(lat, lon)[, 2:1])
                proj.mer.coords <- fm_transform(meridians, crs0 = inla.CRS("longlat"), crs = x)
                proj.mer.coords1 <- matrix(
                    proj.mer.coords[, 1], length(lat),
                    length(lon)
                )
                proj.mer.coords2 <- matrix(
                    proj.mer.coords[, 2], length(lat),
                    length(lon)
                )
                mer.coords <-
                    c(
                        mer.coords,
                        unlist(lapply(
                            seq_along(lon),
                            function(k) {
                                internal.clip(bounds, cbind(
                                    proj.mer.coords1[, k, drop = FALSE],
                                    proj.mer.coords2[, k, drop = FALSE]
                                ),
                                eps = eps
                                )
                            }
                        ),
                        recursive = FALSE
                        )
                    )

                lat <- seq(n[2] * by[2], 90 - 1e-6, length = ceiling((90 - n[2] * by[2]) / 2) + 1)
                meridians <- as.matrix(expand.grid(lat, lon)[, 2:1])
                proj.mer.coords <- fm_transform(meridians, crs0 = inla.CRS("longlat"), crs = x)
                proj.mer.coords1 <- matrix(
                    proj.mer.coords[, 1], length(lat),
                    length(lon)
                )
                proj.mer.coords2 <- matrix(
                    proj.mer.coords[, 2], length(lat),
                    length(lon)
                )
                mer.coords <-
                    c(
                        mer.coords,
                        unlist(lapply(
                            seq_along(lon),
                            function(k) {
                                internal.clip(bounds, cbind(
                                    proj.mer.coords1[, k, drop = FALSE],
                                    proj.mer.coords2[, k, drop = FALSE]
                                ),
                                eps = eps
                                )
                            }
                        ),
                        recursive = FALSE
                        )
                    )
            }
        }

        proj.mer <-
            sp::SpatialLines(list(sp::Lines(mer.coords, ID = "meridians")),
                proj4string = CRS(inla.CRSargs(x))
            )
        if (do.plot) {
            args <- list(x = proj.mer, ...)
            args <- args[intersect(
                names(args),
                union(
                    names(formals(plot.default)),
                    union(
                        names(formals(sp:::plot.Spatial)),
                        names(formals(sp:::plotSpatialLines))
                    )
                )
            )]
            if (add) {
                do.call(plot, c(list(add = TRUE), args))
            } else {
                do.call(plot, args)
                add <- TRUE
            }
        }
    } else {
        proj.mer <- NULL
    }
    if (by[2] > 0) {
        lon <- seq(-180 + 1e-6, 180 - 1e-6, length = 181)
        lat <- ((-n[2]):n[2]) * by[2]
        parallels <- as.matrix(expand.grid(lon, lat))
        proj.par.coords <- fm_transform(parallels, crs0 = inla.CRS("longlat"), crs = x)
        proj.par.coords1 <- matrix(
            proj.par.coords[, 1], length(lon),
            length(lat)
        )
        proj.par.coords2 <- matrix(
            proj.par.coords[, 2], length(lon),
            length(lat)
        )
        proj.par <-
            sp::SpatialLines(list(sp::Lines(
                unlist(lapply(
                    seq_along(lat),
                    function(k) {
                        internal.clip(bounds, cbind(
                            proj.par.coords1[, k],
                            proj.par.coords2[, k]
                        ), eps = eps)
                    }
                ),
                recursive = FALSE
                ),
                ID = "parallels"
            )),
            proj4string = CRS(inla.CRSargs(x))
            )
        if (do.plot) {
            args <- list(x = proj.par, ...)
            args <- args[intersect(
                names(args),
                union(
                    names(formals(plot.default)),
                    union(
                        names(formals(sp:::plot.Spatial)),
                        names(formals(sp:::plotSpatialLines))
                    )
                )
            )]
            if (add) {
                do.call(plot, c(list(add = TRUE), args))
            } else {
                do.call(plot, args)
            }
        }
    } else {
        proj.par <- NULL
    }
    invisible(list(meridians = proj.mer, parallels = proj.par))
}

inla.crs.tissot <- function(x, by = c(30, 30, 30), add = FALSE, do.plot = TRUE,
                            eps = 0.05, diff.eps = 1e-2, ...) {
    inla.fallback_PROJ6("inla.crs.graticule")

    if (is.null(by)) {
        return(invisible(list()))
    }
    if (length(by) < 2) {
        by <- c(by[1], by[1], 30)
    } else if (length(by) < 3) {
        by <- c(by[1:2], 30)
    }
    bounds <- inla.crs.bounds(x)
    n <- c(floor(180 / by[1]), ceiling(90 / by[2]) - 1)
    lon <- ((1 - n[1]):n[1]) * by[1]
    lat <- ((-n[2]):n[2]) * by[2]
    loc0.lon <- loc0.lat <- loc0 <-
        cbind(as.matrix(expand.grid(lat, lon)[, 2:1]), 0)
    loc0.lon[, 1] <- loc0.lon[, 1] + diff.eps / cos(loc0.lat[, 2] * pi / 180)
    loc0.lat[, 2] <- loc0.lat[, 2] + diff.eps
    crs.longlat <- inla.CRS("longlat")

    loc1 <- fm_transform(loc0, crs0 = crs.longlat, crs = x)
    loc1.lon <- fm_transform(loc0.lon, crs0 = crs.longlat, crs = x)
    loc1.lat <- fm_transform(loc0.lat, crs0 = crs.longlat, crs = x)
    ok <- (rowSums(is.na(loc1)) +
        rowSums(is.na(loc1.lon)) +
        rowSums(is.na(loc1.lat)) == 0)
    loc1 <- loc1[ok, , drop = FALSE]
    loc1.lon <- loc1.lon[ok, , drop = FALSE]
    loc1.lat <- loc1.lat[ok, , drop = FALSE]

    diff.lon <- (loc1.lon - loc1) / eps
    diff.lat <- (loc1.lat - loc1) / eps

    scale <- by[3]
    theta <- seq(0, 2 * pi, length = 181)
    ct <- cos(theta) * scale
    st <- sin(theta) * scale

    collection <-
        sp::SpatialLines(list(sp::Lines(
            unlist(lapply(
                seq_len(nrow(loc1)),
                function(k) {
                    loc1.ellipse <-
                        cbind(
                            loc1[k, 1] + diff.lon[k, 1] * ct + diff.lat[k, 1] * st,
                            loc1[k, 2] + diff.lon[k, 2] * ct + diff.lat[k, 2] * st
                        )
                    internal.clip(bounds, loc1.ellipse, eps = eps)
                }
            ),
            recursive = FALSE
            ),
            ID = "parallels"
        )),
        proj4string = CRS(inla.CRSargs(x))
        )
    if (do.plot) {
        args <- list(x = collection, ...)
        args <- args[intersect(
            names(args),
            union(
                names(formals(plot.default)),
                union(
                    names(formals(sp:::plot.Spatial)),
                    names(formals(sp:::plotSpatialLines))
                )
            )
        )]
        if (add) {
            do.call(plot, c(list(add = TRUE), args))
        } else {
            do.call(plot, args)
        }
    }

    invisible(list(tissot = collection))
}



#' Plot CRS and inla.CRS objects
#'
#' Plot the outline of a CRS or inla.CRS projection, with optional graticules
#' (transformed parallels and meridians) and Tissot indicatrices.
#'
#'
#' @aliases plot.inla.CRS plot.CRS
#' @param x A `CRS` or [inla.CRS()] object.
#' @param xlim Optional x-axis limits.
#' @param ylim Optional y-axis limits.
#' @param outline Logical, if `TRUE`, draw the outline of the projection.
#' @param graticule Vector of length at most 3, to plot meridians with spacing
#' `graticule[1]` degrees and parallels with spacing `graticule[2]`
#' degrees. `graticule[3]` optionally specifies the spacing above and
#' below the first and last parallel.  When `graticule[1]==0` no meridians
#' are drawn, and when `graticule[2]==0` no parallels are drawn. Use
#' `graticule=NULL` to skip drawing a graticule.
#' @param tissot Vector of length at most 3, to plot Tissot's indicatrices with
#' spacing `tissot[1]` degrees and parallels with spacing `tissot[2]`
#' degrees. `tissot[3]` specifices a scaling factor.  Use
#' `tissot=NULL` to skip drawing a Tissot's indicatrices.
#' @param asp The aspect ratio for the plot, default 1.
#' @param add If `TRUE`, add the projecton plot to an existing plot.
#' @param eps Clipping tolerance for rudimentary boundary clipping
#' @param \dots Additional arguments passed on to the internal calls to
#' `plot` and `lines`.
#' @author Finn Lindgren <finn.lindgren@@gmail.com>
#' @seealso [inla.CRS()]
#' @examples
#'
#' if (require("sf") && require("sp")) {
#'     oblique <- c(0, 45, 45, 0)
#'     for (projtype in c("longlat", "lambert", "mollweide", "hammer")) {
#'         plot(inla.CRS(projtype), main = projtype)
#'         plot(inla.CRS(projtype, oblique = oblique), main = paste("oblique", projtype))
#'     }
#' }
#' @method plot inla.CRS
#' @export
plot.inla.CRS <- function(x, xlim = NULL, ylim = NULL,
                          outline = TRUE,
                          graticule = c(15, 15, 45),
                          tissot = c(30, 30, 30),
                          asp = 1,
                          add = FALSE,
                          eps = 0.05,
                          ...) {
    inla.fallback_PROJ6("plot.inla.CRS")

    bounds <- inla.crs.bounds(x)
    if (is.null(xlim)) xlim <- bounds$xlim
    if (is.null(ylim)) ylim <- bounds$ylim
    if (!add) {
        args <- list(x = NA, type = "n", xlim = xlim, ylim = ylim, asp = asp, ...)
        args <- args[intersect(names(args), names(formals(plot.default)))]
        do.call(plot, args)
    }
    ## Outline
    if (outline) {
        args <- list(x = bounds$polygon, ...)
        args <- args[intersect(names(args), union(
            names(formals(lines.default)),
            names(formals(plot.xy))
        ))]
        do.call(lines, args)
    }
    ## Graticule
    inla.crs.graticule(x, by = graticule, add = TRUE, do.plot = TRUE, eps = eps, ...)
    ## Tissot
    inla.crs.tissot(x, by = tissot, add = TRUE, do.plot = TRUE, eps = eps, ...)
    invisible(NULL)
}

#' @export
#' @method plot CRS
#' @rdname plot.inla.CRS
plot.CRS <- function(x, xlim = NULL, ylim = NULL,
                     outline = TRUE,
                     graticule = c(15, 15, 45),
                     tissot = c(30, 30, 30),
                     asp = 1,
                     add = FALSE,
                     eps = 0.05,
                     ...) {
    invisible(plot.inla.CRS(x,
        xlim = xlim, ylim = ylim,
        outline = outline, graticule = graticule, tissot = tissot,
        asp = asp, add = add, ...
    ))
}




#' Create a coordinate reference system object
#'
#' Creates either a CRS object or an inla.CRS object, describing a coordinate
#' reference system. Deprecated in favour of `inlabru::fm_CRS()`.
#'
#' @param \dots Arguments passed on to `inlabru::fm_CRS(...)`.
#' @return Either an `sp::CRS` object or an `inla.CRS` object,
#' depending on if the coordinate reference system described by the parameters
#' can be expressed with a pure `sp::CRS` object or not.
#'
#' An S3 `inla.CRS` object is a list, usually (but not necessarily)
#' containing at least one element: \item{crs }{The basic `sp::CRS`
#' object}
#'
#' `inla.wkt_predef` returns a WKT2 string defining a projection
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [sp::CRS()], [crs_wkt()],
#' [inla.sp_get_crs()] [plot.CRS()],
#' [inla.identical.CRS()]
#' @examples
#'
#'
#' if (require("sf") && require("sp")) {
#'     crs1 <- inla.CRS("longlat_globe")
#'     crs2 <- inla.CRS("lambert_globe")
#'     crs3 <- inla.CRS("mollweide_norm")
#'     crs4 <- inla.CRS("hammer_globe")
#'     crs5 <- inla.CRS("sphere")
#'     crs6 <- inla.CRS("globe")
#' }
#' \dontrun{
#' names(inla.wkt_predef())
#' }
#'
#' @export inla.CRS
inla.CRS <- function(...) {
    if (inla.getOption("fmesher.evolution") >= 2L) {
        lifecycle::deprecate_soft(when = "2023.06.06",
                              what = "inla.CRS()",
                              with = "inlabru::fm_CRS()")
    }
    inlabru::fm_CRS(...)
}

#' @return `inla.wkt_predef` returns a WKT2 string defining a projection
#' @examples
#' \dontrun{
#' names(inla.wkt_predef())
#' }
#' @export
#' @rdname inla.CRS

inla.wkt_predef <- function() {
    list(
        hammer_norm = 'PROJCRS["unknown",BASEGEOGCRS["unknown",DATUM["unknown",ELLIPSOID["unknown",0.707106781186548,0,LENGTHUNIT["metre",1,ID["EPSG",9001]]]],PRIMEM["Reference meridian",0,ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]]],CONVERSION["unknown",METHOD["PROJ hammer"]],CS[Cartesian,2],AXIS["(E)",east,ORDER[1],LENGTHUNIT["metre",1,ID["EPSG",9001]]],AXIS["(N)",north,ORDER[2],LENGTHUNIT["metre",1,ID["EPSG",9001]]]]',
        lambert_norm = 'PROJCRS["unknown",BASEGEOGCRS["unknown",DATUM["unknown",ELLIPSOID["unknown",1,0,LENGTHUNIT["metre",1,ID["EPSG",9001]]]],PRIMEM["Reference meridian",0,ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]]],CONVERSION["unknown",METHOD["Lambert Cylindrical Equal Area (Spherical)",ID["EPSG",9834]],PARAMETER["Latitude of 1st standard parallel",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8823]],PARAMETER["Longitude of natural origin",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8802]],PARAMETER["False easting",0,LENGTHUNIT["metre",1],ID["EPSG",8806]],PARAMETER["False northing",0,LENGTHUNIT["metre",1],ID["EPSG",8807]]],CS[Cartesian,2],AXIS["(E)",east,ORDER[1],LENGTHUNIT["metre",1,ID["EPSG",9001]]],AXIS["(N)",north,ORDER[2],LENGTHUNIT["metre",1,ID["EPSG",9001]]]]',
        longlat_norm = 'GEOGCRS["unknown",DATUM["unknown",ELLIPSOID["unknown",1,0,LENGTHUNIT["metre",1,ID["EPSG",9001]]]],PRIMEM["Reference meridian",0,ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]],CS[ellipsoidal,2],AXIS["longitude",east,ORDER[1],ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]],AXIS["latitude",north,ORDER[2],ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]]]',
        mollweide_norm = 'PROJCRS["unknown",BASEGEOGCRS["unknown",DATUM["unknown",ELLIPSOID["unknown",0.707106781186548,0,LENGTHUNIT["metre",1,ID["EPSG",9001]]]],PRIMEM["Reference meridian",0,ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]]],CONVERSION["unknown",METHOD["Mollweide"],PARAMETER["Longitude of natural origin",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8802]],PARAMETER["False easting",0,LENGTHUNIT["metre",1],ID["EPSG",8806]],PARAMETER["False northing",0,LENGTHUNIT["metre",1],ID["EPSG",8807]]],CS[Cartesian,2],AXIS["(E)",east,ORDER[1],LENGTHUNIT["metre",1,ID["EPSG",9001]]],AXIS["(N)",north,ORDER[2],LENGTHUNIT["metre",1,ID["EPSG",9001]]]]',
        hammer_globe = 'PROJCRS["unknown",BASEGEOGCRS["unknown",DATUM["Unknown based on Normal Sphere (r=6370997) ellipsoid",ELLIPSOID["Normal Sphere (r=6370997)",6370997,0,LENGTHUNIT["metre",1,ID["EPSG",9001]]]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8901]]],CONVERSION["unknown",METHOD["PROJ hammer"]],CS[Cartesian,2],AXIS["(E)",east,ORDER[1],LENGTHUNIT["kilometre",1000,ID["EPSG",9036]]],AXIS["(N)",north,ORDER[2],LENGTHUNIT["kilometre",1000,ID["EPSG",9036]]]]',
        lambert_globe = 'PROJCRS["unknown",BASEGEOGCRS["unknown",DATUM["Unknown based on Normal Sphere (r=6370997) ellipsoid",ELLIPSOID["Normal Sphere (r=6370997)",6370997,0,LENGTHUNIT["metre",1,ID["EPSG",9001]]]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8901]]],CONVERSION["unknown",METHOD["Lambert Cylindrical Equal Area (Spherical)",ID["EPSG",9834]],PARAMETER["Latitude of 1st standard parallel",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8823]],PARAMETER["Longitude of natural origin",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8802]],PARAMETER["False easting",0,LENGTHUNIT["kilometre",1000],ID["EPSG",8806]],PARAMETER["False northing",0,LENGTHUNIT["kilometre",1000],ID["EPSG",8807]]],CS[Cartesian,2],AXIS["(E)",east,ORDER[1],LENGTHUNIT["kilometre",1000,ID["EPSG",9036]]],AXIS["(N)",north,ORDER[2],LENGTHUNIT["kilometre",1000,ID["EPSG",9036]]]]',
        longlat_globe = 'GEOGCRS["unknown",DATUM["Unknown based on Normal Sphere (r=6370997) ellipsoid",ELLIPSOID["Normal Sphere (r=6370997)",6370997,0,LENGTHUNIT["metre",1,ID["EPSG",9001]]]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8901]],CS[ellipsoidal,2],AXIS["longitude",east,ORDER[1],ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]],AXIS["latitude",north,ORDER[2],ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]]]',
        mollweide_globe = 'PROJCRS["unknown",BASEGEOGCRS["unknown",DATUM["Unknown based on Normal Sphere (r=6370997) ellipsoid",ELLIPSOID["Normal Sphere (r=6370997)",6370997,0,LENGTHUNIT["metre",1,ID["EPSG",9001]]]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8901]]],CONVERSION["unknown",METHOD["Mollweide"],PARAMETER["Longitude of natural origin",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8802]],PARAMETER["False easting",0,LENGTHUNIT["kilometre",1000],ID["EPSG",8806]],PARAMETER["False northing",0,LENGTHUNIT["kilometre",1000],ID["EPSG",8807]]],CS[Cartesian,2],AXIS["(E)",east,ORDER[1],LENGTHUNIT["kilometre",1000,ID["EPSG",9036]]],AXIS["(N)",north,ORDER[2],LENGTHUNIT["kilometre",1000,ID["EPSG",9036]]]]',
        sphere = 'GEODCRS["unknown",DATUM["unknown",ELLIPSOID["unknown",1,0,LENGTHUNIT["metre",1,ID["EPSG",9001]]]],PRIMEM["Reference meridian",0,ANGLEUNIT["degree",0.0174532925199433,ID["EPSG",9122]]],CS[Cartesian,3],AXIS["(X)",geocentricX,ORDER[1],LENGTHUNIT["metre",1,ID["EPSG",9001]]],AXIS["(Y)",geocentricY,ORDER[2],LENGTHUNIT["metre",1,ID["EPSG",9001]]],AXIS["(Z)",geocentricZ,ORDER[3],LENGTHUNIT["metre",1,ID["EPSG",9001]]]]',
        globe = 'GEODCRS["unknown",DATUM["Unknown based on Normal Sphere (r=6370997) ellipsoid",ELLIPSOID["Normal Sphere (r=6370997)",6370997,0,LENGTHUNIT["metre",1,ID["EPSG",9001]]]],PRIMEM["Greenwich",0,ANGLEUNIT["degree",0.0174532925199433],ID["EPSG",8901]],CS[Cartesian,3],AXIS["(X)",geocentricX,ORDER[1],LENGTHUNIT["kilometre",1000,ID["EPSG",9036]]],AXIS["(Y)",geocentricY,ORDER[2],LENGTHUNIT["kilometre",1000,ID["EPSG",9036]]],AXIS["(Z)",geocentricZ,ORDER[3],LENGTHUNIT["kilometre",1000,ID["EPSG",9036]]]]'
    )
}





#' Internal WKT handling
#'
#' Conversion between WKT and a tree representation
#'
#'
#' @aliases inla.as.wkt_tree.wkt inla.as.wkt.wkt_tree inla.wkt_tree_get_item
#' inla.wkt_tree_set_item
#' @param x A WKT2 string, or a `wkt_tree` list structure
#' @param \dots Unused
#' @param item character vector with item labels identifying a parameter item
#' entry.
#' @param duplicate For items that have more than one match, `duplicate`
#' indicates the index number of the desired version. Default: 1
#' @param item_tree An item tree identifying a parameter item entry
#' @export inla.as.wkt_tree.wkt
#' @rdname wkt_tree
inla.as.wkt_tree.wkt <- function(x, ...) {
    if (inla.getOption("fmesher.evolution") >= 2L) {
        lifecycle::deprecate_soft(when = "2023.06.06",
                                  what = "inla.as.wkt_tree.wkt()",
                                  with = "inlabru::fm_wkt_as_wkt_tree()")
    }
    inlabru::fm_wkt_as_wkt_tree(x, ...)
}

#' @rdname wkt_tree
#' @param pretty logical
#' @export

inla.as.wkt.wkt_tree <- function(x, pretty = FALSE, ...) {
    if (inla.getOption("fmesher.evolution") >= 2L) {
        lifecycle::deprecate_soft(when = "2023.06.06",
                                  what = "inla.as.wkt.wkt_tree()",
                                  with = "inlabru::fm_wkt_tree_as_wkt()")
    }
    inlabru::fm_wkt_tree_as_wkt(x, pretty = pretty, ...)
}

#' @param item character vector with item labels identifying a parameter item
#' entry.
#' @param duplicate For items that have more than one match, `duplicate`
#' indicates the index number of the desired version. Default: 1
#' @rdname wkt_tree
#' @export

inla.wkt_tree_get_item <- function(x, item, duplicate = 1) {
    if (inla.getOption("fmesher.evolution") >= 2L) {
        lifecycle::deprecate_soft(when = "2023.06.06",
                                  what = "inla.wkt_tree_get_item()",
                                  with = "inlabru::fm_wkt_tree_get_item()")
    }
    inlabru::fm_wkt_tree_get_item(x, item = item, duplicate = duplicate)
}

#' @param item_tree An item tree identifying a parameter item entry
#' @rdname wkt_tree
#' @export

inla.wkt_tree_set_item <- function(x, item_tree, duplicate = 1) {
    if (inla.getOption("fmesher.evolution") >= 2L) {
        lifecycle::deprecate_soft(when = "2023.06.06",
                                  what = "inla.wkt_tree_set_item()",
                                  with = "inlabru::fm_wkt_tree_set_item()")
    }
    inlabru::fm_wkt_tree_set_item(x, item_tree = item_tree, duplicate = duplicate)
}



#' Show expanded CRS arguments
#'
#' Wrapper for `sp::CRS` and `inla.CRS` objects to extract the
#' coordinate reference system argument string.
#' Should no longer be used with PROJ6;
#' see [inla.crs_get_wkt()]
#'
#' @aliases inla.CRSargs inla.as.list.CRS inla.as.list.CRSargs inla.as.CRS.list
#' inla.as.CRSargs.list
#' @param x An `sp::CRS` or `inla.CRS` object (for
#' `inla.CRSargs` and `inla.as.list.CRS`), a character string (for
#' `inla.as.list.CRSargs`), or a list (for `inla.as.CRS.list` and
#' `inla.as.CRSargs.list`).
#' @param \dots Additional arguments passed on to other methods.
#' @return For `inla.CRSargs` and `inla.as.CRSargs.list`, a character
#' string with PROJ.4 arguments.
#'
#' For `inla.as.list.CRS` and `inla.as.list.CRSargs`, a list of
#' name/value pairs.
#'
#' For `inla.as.CRS.list`, a `CRS` or `inla.CRS` object.
#' @author Finn Lindgren <finn.lindgren@@gmail.com>
#' @seealso [inla.CRS()]
#' @examples
#'
#' if (require("sf") && require("sp")) {
#'     crs0 <- inla.CRS("longlat")
#'     p4s <- inla.CRSargs(crs0)
#'     lst <- inla.as.list.CRSargs(p4s)
#'     crs1 <- inla.as.CRS.list(lst)
#'     lst$a <- 2
#'     crs2 <- inla.CRS(p4s, args = lst)
#'     print(inla.CRSargs(crs0))
#'     print(inla.CRSargs(crs1))
#'     print(inla.CRSargs(crs2))
#' }
#' @export
#' @rdname CRSargs
inla.CRSargs <- function(x, ...) {
    if (inla.getOption("fmesher.evolution") >= 2L) {
        lifecycle::deprecate_soft(when = "2023.06.06",
                                  what = "inla.CRSargs()",
                                  with = "inlabru::fm_proj4string()")
    }

    if (inherits(x, "inla.CRS")) {
        x <- x[["crs"]]
    }
    if (is.null(x)) {
        as.character(NA)
    } else {
        inla.require("sf", stop.on.error = TRUE)
        sf::st_crs(x)$proj4string
    }
}


#' @details * `inla.as.CRSargs.list`: CRS proj4 string for name=value pair list
#'
#' @export
#' @rdname CRSargs
inla.as.CRSargs.list <- function(x, ...) {
    paste(lapply(
        names(x),
        function(xx) {
            if (is.na(x[[xx]])) {
                paste("+", xx, sep = "")
            } else {
                paste("+", xx, "=", x[[xx]], sep = "")
            }
        }
    ),
    collapse = " "
    )
}

#' @details * `inla.as.list.CRSargs`: List of name=value pairs from CRS proj4 string
#'
#' @export
#' @rdname CRSargs
inla.as.list.CRSargs <- function(x, ...) {
    if (is.na(x)) {
        return(list())
    }
    if (!is.character(x)) {
        stop("proj4string must be of class character")
    }
    do.call(c, lapply(
        strsplit(
            x = strsplit(
                x = paste(" ", x, sep = ""),
                split = " \\+"
            )[[1]][-1],
            split = "="
        ),
        function(x) {
            xx <- list(x[2])
            names(xx) <- x[1]
            xx
        }
    ))
}


#' @export
#' @rdname CRSargs
inla.as.list.CRS <- function(x, ...) {
    inla.as.list.CRSargs(inla.CRSargs(x))
}

#' @export
#' @rdname CRSargs
inla.as.CRS.list <- function(x, ...) {
    inla.CRS(args = x)
}





#' @param crs A `sp::CRS` or `inla.CRS` object
#' @param wkt A WKT2 character string
#' @param unit character, name of a unit. Supported names are
#' "metre", "kilometre", and the aliases "meter", "m", International metre",
#' "kilometer", and "km", as defined by `inla.wkt_unit_params` or the
#' `params` argument. (For legacy PROJ4 use, only "m" and "km" are
#' supported)
#' @param params Length unit definitions, in the list format produced by
#' `inla.wkt_unit_params()`, Default: NULL, which invokes
#' `inla.wkt_unit_params()`
#' @return For `inla.wkt_unit_params`, a
#' list of named unit definitions
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @examples
#' \dontrun{
#'   c1 <- inla.CRS("globe")
#'   inla.crs_get_lengthunit(c1)
#'   c2 <- inla.crs_set_lengthunit(c1, "km")
#'   inla.crs_get_lengthunit(c2)
#' }
#' @export
#' @seealso [inla.sp_get_crs()]
#' @aliases crs_wkt
#' @rdname crs_wkt

inla.wkt_unit_params <- function() {
    params <- list(
        "metre" =
            list(
                '"metre"',
                "1",
                list(
                    label = "ID",
                    params = list('"EPSG"', "9001")
                )
            ),
        "kilometre" =
            list(
                '"kilometre"',
                "1000",
                list(
                    label = "ID",
                    params = list('"EPSG"', "9036")
                )
            )
    )
    params[["meter"]] <- params[["metre"]]
    params[["m"]] <- params[["metre"]]
    params[["International metre"]] <- params[["metre"]]
    params[["kilometer"]] <- params[["kilometre"]]
    params[["km"]] <- params[["kilometre"]]
    params
}

#' @export
#' @rdname crs_wkt
#' @return For `inla.wkt_get_lengthunit`, a
#' list of length units used in the wkt string, excluding the ellipsoid radius
#' unit.
#' @importFrom inlabru fm_length_unit

inla.wkt_get_lengthunit <- function(wkt) {
    if (inla.getOption("fmesher.evolution") >= 2L) {
        lifecycle::deprecate_soft(when = "2023.06.06",
                                  what = "inla.wkt_get_lengthunit()",
                                  with = "inlabru::fm_length_unit()")
    }
    inlabru::fm_length_unit(wkt)
}

#' @export
#' @rdname crs_wkt
#' @return For `inla.wkt_set_lengthunit`, a
#' WKT2 string with altered length units.
#' Note that the length unit for the ellipsoid radius is unchanged.
#' @importFrom inlabru `fm_length_unit<-`

inla.wkt_set_lengthunit <- function(wkt, unit, params = NULL) {
    if (inla.getOption("fmesher.evolution") >= 2L) {
        lifecycle::deprecate_soft(when = "2023.06.06",
                                  what = "inla.wkt_set_lengthunit()",
                                  with = "`inlabru::fm_length_unit<-`()")
    }

    if (!is.null(params)) {
        warning("Non-NULL 'params' for inla.wkt_set_lengthunit ignored.")
    }
    inlabru::fm_length_unit(wkt) <- unit
    wkt
}

#' @return For `inla.crs_get_wkt`, WKT2 string.
#' @export
#' @rdname crs_wkt

inla.crs_get_wkt <- function(crs) {
    if (inla.getOption("fmesher.evolution") >= 2L) {
        lifecycle::deprecate_soft(when = "2023.06.06",
                                  what = "inla.crs_get_wkt()",
                                  with = "inlabru::fm_wkt()")
    }
    
    inlabru::fm_wkt(crs)
}

#' @return For `inla.crs_get_lengthunit`, a
#' list of length units used in the wkt string, excluding the ellipsoid radius
#' unit. (For legacy PROJ4 code, the raw units from the proj4string are
#' returned, if present.)
#' @export
#' @rdname crs_wkt

inla.crs_get_lengthunit <- function(crs) {
    if (inla.getOption("fmesher.evolution") >= 2L) {
        lifecycle::deprecate_soft(when = "2023.06.06",
                                  what = "inla.crs_get_lengthunit()",
                                  with = "inlabru::fm_length_unit()")
    }
    inlabru::fm_length_unit(crs)
}

#' @return For `inla.crs_set_lengthunit`, a `sp::CRS` object with
#' altered length units.
#' Note that the length unit for the ellipsoid radius is unchanged.
#' @export
#' @rdname crs_wkt

inla.crs_set_lengthunit <- function(crs, unit, params = NULL) {
    if (inla.getOption("fmesher.evolution") >= 2L) {
        lifecycle::deprecate_soft(when = "2023.06.06",
                                  what = "inla.crs_set_lengthunit()",
                                  with = "`inlabru::fm_length_unit<-`()")
    }
    if (!is.null(params)) {
        warning("Non-NULL 'params' for inla.crs_set_lengthunit ignored.")
    }
    inlabru::fm_length_unit(crs) <- unit
    crs
}





inla.wkt_tree_projection_type <- function(wt) {
    axis1 <- inla.wkt_tree_get_item(wt, "AXIS", 1)
    axis2 <- inla.wkt_tree_get_item(wt, "AXIS", 2)
    if (identical(axis1[["params"]][[1]], '"longitude"') &&
        identical(axis2[["params"]][[1]], '"latitude"')) {
        return("longlat")
    }
    conversion <- inla.wkt_tree_get_item(wt, "CONVERSION")
    if (!is.null(conversion)) {
        method <- inla.wkt_tree_get_item(conversion, "METHOD")
        if (identical(method[["params"]][[1]], '"Lambert Cylindrical Equal Area (Sherical)"')) {
            return("lambert")
        }
        if (identical(method[["params"]][[1]], '"Mollweide"')) {
            return("mollweide")
        }
        if (identical(method[["params"]][[1]], '"PROJ hammer"')) {
            return("hammer")
        }
        if (identical(method[["params"]][[1]], '"tmerc"')) {
            return("tmerc")
        }
    }
    NULL
}

inla.wkt_projection_type <- function(wkt) {
    wt <- inlabru::fm_wkt_as_wkt_tree(wkt)
    inla.wkt_tree_projection_type(wt)
}

inla.crs_projection_type <- function(crs) {
    wkt <- inlabru::fm_wkt(crs)
    inla.wkt_projection_type(wkt)
}

## +proj=longlat in (-180,180)x(-90,90)
## +proj=moll in (-2,2)x(-1,1) scaled by +a and +b, and +units
## +proj=lambert in (-pi,pi)x(-1,1) scaled by +a and +b, and +units
inla.crs.bounds <- function(crs, warn.unknown = FALSE) {
    wkt <- inla.crs_get_wkt(crs)
    wt <- inla.as.wkt_tree.wkt(wkt)
    type <- inla.wkt_tree_projection_type(wt)
    
    if (is.null(type)) {
        if (inla.wkt_is_geocent(wkt)) {
            bounds <-
                list(
                    type = "rectangle",
                    xlim = c(-Inf, Inf),
                    ylim = c(-Inf, Inf)
                )
        } else {
            if (warn.unknown) {
                warning(
                    "Could not determine shape of transformation bounds. Using infinite rectangle."
                )
            }
            bounds <-
                list(
                    type = "rectangle",
                    xlim = c(-Inf, Inf),
                    ylim = c(-Inf, Inf)
                )
        }
    } else if (type == "longlat") {
        bounds <-
            list(
                type = "rectangle",
                xlim = c(-180, 180),
                ylim = c(-90, 90)
            )
    } else if (type == "lambert") {
        axis <- c(pi, 1)
        radius <- inla.wkt_get_ellipsoid_radius(wkt)
        axis[1] <- axis[1] * radius
        # TODO: handle eccentricity
        axis[2] <- axis[2] * sqrt(radius) * sqrt(radius)
        # TODO: Handle units"
        bounds <- list(
            type = "rectangle",
            xlim = c(-1, 1) * axis[1],
            ylim = c(-1, 1) * axis[2]
        )
    } else if (type %in% c("mollweide", "hammer")) {
        axis <- c(2, 1)
        center <- c(0, 0)
        radius <- inla.wkt_get_ellipsoid_radius(wkt)
        axis[1] <- axis[1] * radius / sqrt(1 / 2)
        axis[2] <- axis[2] * radius / sqrt(1 / 2)
        # TODO: Handle "units"
        bounds <- list(
            type = "ellipse",
            axis = axis,
            center = center,
            xlim = center[1] + c(-1, 1) * axis[1],
            ylim = center[2] + c(-1, 1) * axis[2]
        )
    } else if (type == "tmerc") {
        bounds <-
            list(
                type = "rectangle",
                xlim = c(-Inf, Inf),
                ylim = c(-Inf, Inf)
            )
    } else {
        stop("'inla.crs.bounds' internal error: transformation detected but not handled.")
    }
    
    if (bounds$type == "rectangle") {
        bounds$polygon <- cbind(
            bounds$xlim[c(1, 2, 2, 1, 1)],
            bounds$ylim[c(1, 1, 2, 2, 1)]
        )
    } else if (bounds$type == "ellipse") {
        theta <- seq(0, 2 * pi, length = 1000)
        bounds$polygon <- cbind(
            bounds$center[1] + bounds$axis[1] * cos(theta),
            bounds$center[2] + bounds$axis[2] * sin(theta)
        )
    } else {
        stop("Unknown transformation type. This should not happen.")
    }
    bounds
}

## TRUE/FALSE for points inside/outside projection domain.
inla.crs.bounds.check <- function(x, bounds) {
    inla.require.inherits(x, "matrix")
    if (all(is.finite(bounds$xlim)) && all(is.finite(bounds$ylim))) {
        (sp::point.in.polygon(
            x[, 1], x[, 2],
            bounds$polygon[, 1], bounds$polygon[, 2]
        )
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




#' Test CRS and inla.CRS for equality
#'
#' Wrapper for identical, optionally testing only the CRS part of two objects
#' Deprecated in favour of `inlabru::fm_identical_CRS()`
#'
#' @export inla.identical.CRS
inla.identical.CRS <- function(...) {
    inlabru::fm_identical_CRS(...)
}






#' Wrapper method for `inlabru::fm_transform`
#'
#' Handles transformation of various inla objects accorting to coordinate
#' reference systems of `sf::crs`, `sp::CRS` or `inla.CRS` class.
#'
#'
#' @param x The object that should be transformed from it's current CRS to a
#' new CRS
#' @param crs0 The source crs object, only needed for raw coordinate inputs
#' @param CRSobj The target crs object
#' @param passthrough default FALSE. Setting to TRUE allows objects with no CRS
#' information to be passed through without transformation.
#' @param \dots Potential additional arguments
#' @return The object is returned with its coordinates transformed
#' @author Finn Lindgren <finn.lindgren@@gmail.com>
#' @seealso [inla.CRS()]
#' @examples
#'
#' if (require("sf") && require("sp") && require("inlabru")) {
#'     latt <- inla.mesh.lattice(-10:10, 40:60)
#'     mesh1 <- inla.mesh.create(
#'         lattice = latt, extend = FALSE, refine = FALSE,
#'         crs = fm_CRS("longlat_norm")
#'     )
#'     mesh2 <- fm_transform(mesh1, fm_crs("lambert_globe"))
#'     summary(mesh1)
#'     summary(mesh2)
#' }
#' @export inla.spTransform
inla.spTransform <- function(x, CRSobj, ...) {
    inlabru::fm_transform(x, crs = CRSobj, ...)
}




inla.crs_detect_manifold <- function(crs) {
    inlabru::fm_crs_detect_manifold(crs)
}





#' @title Convert `sp` objects to `inla.mesh.segment` objects.
#'
#' @description Wrapper for `inlabru::fm_as_inla_mesh_segment`
#'
#' @param sp An `sp` polygon object of class `Polygon`,
#' `Polygons`, `SpatialPolygons`, or `SpatialPolygonsDataFrame`.
#' @param join If `TRUE`, join multiple polygons into a single segment
#' (possibly non-simply connected).
#' @param grp Group ID specification for each polygon, as used by
#' [inla.mesh.segment()], one ID per polygon.
#' @param reverse Logical, indicating if the line sequence should be traversed
#' backwards.
#' @param crs An optional `CRS` or `inla.CRS` object
#' @param \dots Additional arguments passed on to other methods.
#' @param is.bnd logical
#' @return A [inla.mesh.segment()] object, or a list of
#' [inla.mesh.segment()] objects.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.mesh.segment()]
#' @export as.inla.mesh.segment
#' @importFrom inlabru fm_as_inla_mesh_segment
as.inla.mesh.segment <-
    function(sp, ...) {
        inlabru::fm_as_inla_mesh_segment(sp, ...)
    }

#' @export
#' @rdname as.inla.mesh.segment
inla.sp2segment <-
    function(sp, ...) {
        inlabru::fm_as_inla_mesh_segment(sp, ...)
    }
