% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spmesh.R
\name{crs_wkt}
\alias{crs_wkt}
\alias{inla.wkt_is_geocent}
\alias{inla.crs_is_geocent}
\alias{inla.wkt_get_ellipsoid_radius}
\alias{inla.crs_get_ellipsoid_radius}
\alias{inla.wkt_set_ellipsoid_radius}
\alias{inla.crs_set_ellipsoid_radius}
\alias{inla.wkt_unit_params}
\alias{inla.wkt_get_lengthunit}
\alias{inla.wkt_set_lengthunit}
\alias{inla.crs_get_wkt}
\alias{inla.crs_get_lengthunit}
\alias{inla.crs_set_lengthunit}
\title{Handling CRS/WKT}
\usage{
inla.wkt_is_geocent(wkt)

inla.crs_is_geocent(crs)

inla.wkt_get_ellipsoid_radius(wkt)

inla.crs_get_ellipsoid_radius(crs)

inla.wkt_set_ellipsoid_radius(wkt, radius)

inla.crs_set_ellipsoid_radius(crs, radius)

inla.wkt_unit_params()

inla.wkt_get_lengthunit(wkt)

inla.wkt_set_lengthunit(wkt, unit, params = NULL)

inla.crs_get_wkt(crs)

inla.crs_get_lengthunit(crs)

inla.crs_set_lengthunit(crs, unit, params = NULL)
}
\arguments{
\item{wkt}{A WKT2 character string}

\item{crs}{A \code{sp::CRS} or \code{inla.CRS} object}

\item{radius}{numeric}

\item{unit}{character, name of a unit. Supported names are
"metre", "kilometre", and the aliases "meter", "m", International metre",
"kilometer", and "km", as defined by \code{inla.wkt_unit_params} or the
\code{params} argument. (For legacy PROJ4 use, only "m" and "km" are
supported)}

\item{params}{Length unit definitions, in the list format produced by
\code{inla.wkt_unit_params()}, Default: NULL, which invokes
\code{inla.wkt_unit_params()}}
}
\value{
For \code{inla.wkt_unit_params}, a list of named unit definitions

For \code{inla.wkt_get_lengthunit}, a list of length units used in the wkt
string, excluding the ellipsoid radius unit.

For \code{inla.wkt_set_lengthunit}, a WKT2 string with altered length units.
Note that the length unit for the ellipsoid radius is unchanged.

For \code{inla.crs_get_wkt}, WKT2 string.

For \code{inla.crs_get_lengthunit}, a list of length units used in the wkt
string, excluding the ellipsoid radius unit. (For legacy PROJ4 code, the raw
units from the proj4string are returned, if present.)

For \code{inla.crs_set_lengthunit}, a \code{sp::CRS} object with altered
length units. Note that the length unit for the ellipsoid radius is
unchanged.

For \code{inla.wkt_unit_params}, a
list of named unit definitions

For \code{inla.wkt_get_lengthunit}, a
list of length units used in the wkt string, excluding the ellipsoid radius
unit.

For \code{inla.wkt_set_lengthunit}, a
WKT2 string with altered length units.
Note that the length unit for the ellipsoid radius is unchanged.

For \code{inla.crs_get_wkt}, WKT2 string.

For \code{inla.crs_get_lengthunit}, a
list of length units used in the wkt string, excluding the ellipsoid radius
unit. (For legacy PROJ4 code, the raw units from the proj4string are
returned, if present.)

For \code{inla.crs_set_lengthunit}, a \code{sp::CRS} object with
altered length units.
Note that the length unit for the ellipsoid radius is unchanged.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} in favour of \code{\link[fmesher:fm_crs_wkt]{fmesher::fm_wkt()}} and related
methods.

Get and set CRS object or WKT string properties.
}
\section{Functions}{
\itemize{
\item \code{inla.wkt_is_geocent()}: \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} in favour of
\code{\link[fmesher:fm_crs_wkt]{fmesher::fm_wkt_is_geocent()}}

\item \code{inla.crs_is_geocent()}: \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} in favour of
\code{\link[fmesher:fm_crs_wkt]{fmesher::fm_crs_is_geocent()}}

\item \code{inla.wkt_get_ellipsoid_radius()}: \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} in favour of
\code{\link[fmesher:fm_crs_wkt]{fmesher::fm_ellipsoid_radius()}}

\item \code{inla.crs_get_ellipsoid_radius()}: \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} in favour of
\code{\link[fmesher:fm_crs_wkt]{fmesher::fm_ellipsoid_radius()}}

\item \code{inla.wkt_set_ellipsoid_radius()}: \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} in favour of
\code{\link[fmesher:fm_crs_wkt]{fmesher::fm_wkt_set_ellipsoid_radius()}}

\item \code{inla.crs_set_ellipsoid_radius()}: \ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#deprecated}{\figure{lifecycle-deprecated.svg}{options: alt='[Deprecated]'}}}{\strong{[Deprecated]}} in favour of
\code{\link[fmesher:fm_crs_wkt]{fmesher::fm_ellipsoid_radius<-()}}

}}
\examples{
\dontrun{
c1 <- fmesher::fm_CRS("globe")
inla.crs_get_lengthunit(c1)
c2 <- inla.crs_set_lengthunit(c1, "km")
inla.crs_get_lengthunit(c2)
}

\dontrun{
c1 <- inla.CRS("globe")
inla.crs_get_lengthunit(c1)
c2 <- inla.crs_set_lengthunit(c1, "km")
inla.crs_get_lengthunit(c2)
}
}
\seealso{
\code{\link[=inla.sp_get_crs]{inla.sp_get_crs()}}

\code{\link[=inla.sp_get_crs]{inla.sp_get_crs()}}
}
\author{
Finn Lindgren \email{finn.lindgren@gmail.com}
}
