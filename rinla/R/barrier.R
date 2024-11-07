#' Functions for defining the Barrier models
#'
#' Functions for defining Barrier models as an `inla rgeneric` model
#'
#' This model is described in the ArXiv preprint arXiv:1608.03787.  For
#' examples, see <https://haakonbakkagit.github.io/btopic128.html>
#'
#' @aliases inla.barrier barrier inla.barrier.pcmatern barrier.pcmatern
#' inla.barrier.polygon barrier.polygon inla.barrier.q barrier.q
#' inla.barrier.fem barrier.fem
#' @param mesh The mesh to build the model on, from inla.mesh.2d
#' @param barrier.triangles The numerical ids of the triangles that make up the
#' barrier area
#' @param prior.range 2 parameters `(range0,Prange)` for the prior spatial
#' range.  If `Prange` is `NA`, then `range0` is used as a fixed
#' range value (not tested).
#' @param prior.sigma 2 parameters `(sig0,Psig)` for the prior marginal
#' standard deviation sigma.  If `Psig` is `NA`, then `sig0` is
#' used as a fixed sigma value (not tested).
#' @param range.fraction The length of the spatial range inside the barrier
#' area, as a fraction of the range parameter.
#' @param Omega Advanced option for creating a set of permeable barriers (not
#' documented)
#' @param enable.INLAspacetime Use the implentation in the package `INLAspacetime`
#' instead if available (default TRUE) if its available.
#' You may need set this option to `FALSE` if you want to
#' extract properties of the model for other use,  like for example `inla.rgeneric.q`.
#' @return
#' * `inla.barrier.pcmatern` gives the (rgeneric) model object
#' for fitting the model in INLA
#' * `inla.barrier.polygon` gives the polygon
#' around the barrier (mainly for plotting)
#' * `inla.barrier.q` is an
#' internal method producing the Q matrix from a result of inla.barrier.fem,
#' * `inla.barrier.fem` is an internal method producing the Finite Element
#' matrices.
#' @author Haakon Bakka \email{bakka@@r-inla.org}
#' @seealso inla.spde2.pcmatern
#'
#' @export
#' @rdname inla.barrier
#' @details * `inla.barrier.pcmatern`
#' This function creates the model component used in inla(...)
`inla.barrier.pcmatern` <- function(mesh, barrier.triangles, prior.range, prior.sigma,
                                    range.fraction = 0.2, enable.INLAspacetime = TRUE) {
    ## Give default values if absolutely needed
    if (missing(prior.range)) {
        warning("Arbitrary prior values chosen automatically. This may suffice for a first attempt, 
            but should be changed in any serious analysis.")
        prior.range <- c(diff(range(mesh$loc[, 1])) / 5, 0.5)
    }

    if (missing(prior.sigma)) {
        prior.sigma <- c(1, 0.5)
    }

    ## Input verification
    stopifnot(inherits(mesh, "inla.mesh"))
    stopifnot(range.fraction > 0.000001)

    if (enable.INLAspacetime && requireNamespace("INLAspacetime")) {
	warning("Using implementation from the `INLAspacetime` package")
	return(INLAspacetime::barrierModel.define(
                                  mesh = mesh,
                                  barrier.triangles = barrier.triangles,
                                  prior.range = prior.range, 
                                  prior.sigma = prior.sigma,
                                  range.fraction = range.fraction))
    } else {
	warning(paste(
            "Please install the `INLAspacetime` package\n",
            "which contains an implementation that runs faster!"))
    }

    ## ## ## FUNCTIONS FOR RGENERIC MODEL SETUP ## ## ##

    ## This function is the model component definition in the rgeneric inla framework
    ## See inla.doc('rgeneric')
    barrier.rgeneric.model <- function(cmd = c(
                                           "graph", "Q", "mu", "initial",
                                           "log.norm.const", "log.prior", "quit"
                                       ),
                                       theta = NULL) {

        envir = parent.env(environment())

        ## Dynamic input
        ## theta = log(sigma), log(range1), log(range2), ...

        ## Static input
        prior.sigma <- obj$prior.sigma
        prior.range <- obj$prior.range

        ## The finite element matrices
        fem <- obj$fem

        ## The function for creating the precision matrix
        inla.barrier.q <- obj$inla.barrier.q

        ## The fraction used for the barrier model
        range.fraction <- obj$range.fraction

        ## Initial values on the internal scale
        ## Can be overrridden in inla(...)
        ## Warning: Not implemented for multiple ranges
        initial <- function(theta) {
            initial.theta <- c()
            if (!is.na(prior.sigma[2])) {
                initial.theta <- c(0, initial.theta)
            }
            if (!is.na(prior.range[2])) {
                initial.theta <- c(initial.theta, 0)
            }

            return(initial.theta)
        }

        ## This is the log-prior for the internal parametrisation theta
        ## Log of probability for log transform of exponential distribution
        log.prior <- function(theta) {
            ## Lambdas:
            ## The prior parameters are the lambdas in the exponential
            ## priors for standard deviation and inverse-range
            ## - the first is log(prob)/exceed, the second log(prob)*exceed
            ## - the second is exponential for inverse range, therefore multiplication!

            val <- 0
            ## If the sigma is not fixed
            if (!is.na(prior.sigma[2])) {
                lambda0 <- -log(prior.sigma[2]) / prior.sigma[1]
                ## Prior for standard deviation
                val <- val + log(lambda0) - lambda0 * exp(theta[1]) + theta[1]
                ## Which of the thetas are ranges ...
                ## ... not the first
                theta.ran <- theta[-1]
            } else {
                ## ... all of them
                theta.ran <- theta
            }

            ## If the range(s) is/are not fixed
            if (!is.na(prior.range[2])) {
                lambda1 <- -log(prior.range[2]) * prior.range[1]

                ## Prior for range(s)
                ## All are forced to have the same prior (if there are more than 1)
                for (logrange in theta.ran) {
                    val <- val + log(lambda1) - lambda1 * exp(-logrange) + -logrange
                }
            }
            return(val)
        }

        Q <- function(theta) {
            ## Make a theta that contains all variables
            theta.full <- theta
            ## If there is no theta param for sigma
            ## add the fixed sigma value
            if (is.na(prior.sigma[2])) {
                theta.full <- c(log(prior.sigma[1]), theta.full)
            }
            ## If there is no theta param for sigma
            ## add the fixed range value
            if (is.na(prior.range[2])) {
                theta.full <- c(theta.full, log(prior.range[1]))
                ## How it would be for multiple ranges (not barrier model)
                                        # all.ranges = rep(prior.range[1], fem$hdim)
                                        # theta.full = c(theta.full, log(all.ranges))
            }

            ## Not implemented for multiple ranges
            ## Hence, we should always have 2 theta param at this stage
            stopifnot(length(theta.full) == 2)

            ## Construct the precision matrix
            ## Not from the raw theta, but using both parameters
            Q <- inla.barrier.q(fem = fem, ranges = exp(theta.full[2]) * c(1, range.fraction),
                                sigma = exp(theta.full[1]), envir = envir)
            return(Q)
        }

        ## Here we create the graph by just calling the Q function a few times
        ## on some arbitrary inputs (several: to be sure we do not get accidental 0s)
        graph <- function(theta) {
            requireNamespace(methods)
            ntheta <- 2 # only for barrier model
            theta.full <- (1:ntheta) / 3.217233456
            G1 <- inla.barrier.q(fem = fem, ranges = exp(theta.full[2]) * c(1, range.fraction), sigma = exp(theta.full[1]))
            G1[G1 != 0] <- 1
            theta.full <- (1:ntheta)^2 / 12.1543534
            G2 <- inla.barrier.q(fem = fem, ranges = exp(theta.full[2]) * c(1, range.fraction), sigma = exp(theta.full[1]))
            G2[G2 != 0] <- 1

            return(G1 + G2)
        }

        ## These can always be like this for all rgeneric models
        mu <- function(theta) numeric(0)
        log.norm.const <- function(theta) numeric(0)
        quit <- function(theta) invisible()

        val <- do.call(match.arg(cmd), args = list(theta))
        return(val)
    }

    ## Create a valid Omega from the barrier.triangles
    barrier.triangles <- unique(barrier.triangles)

    ## Create the barrier model
    obj <- list()
    obj$prior.sigma <- prior.sigma
    obj$prior.range <- prior.range
    obj$range.fraction <- range.fraction
    obj$inla.barrier.q <- inla.barrier.q

    obj$fem <- inla.barrier.fem(mesh, barrier.triangles = barrier.triangles)
    barrier.model <- inla.rgeneric.define(model = barrier.rgeneric.model, optimize = TRUE, obj = obj)

    if (!is.na(prior.sigma[2]) && !is.na(prior.range[2])) {
        ## All ok
    } else {
        warning("Not properly tested, let us know if you have problems.")
    }

    return(barrier.model)
}

#' @export
#' @rdname inla.barrier
#' @details * `inla.barrier.polygon` This function constructs SpatialPolygons for the different subdomains (areas)
`inla.barrier.polygon` <- function(mesh, barrier.triangles, Omega = NULL) {
    ## Requires an inla mesh to work
    stopifnot(inherits(mesh, "inla.mesh"))
    ## Requires sf for combining polygons
    inla.require("sf", stop.on.error = TRUE)

    if (missing(barrier.triangles)) {
        ## Use Omega
    } else {
        ## Create a valid Omega
        barrier.triangles <- unique(barrier.triangles)
        ## Number of triangles in total
        t <- length(mesh$graph$tv[, 1])
        remaining <- setdiff(1:t, barrier.triangles)
        if (!is.null(Omega)) warning("Omega is replaced by barrier.triangles")
        Omega <- list(remaining, barrier.triangles)
    }

    Omega.SP.list <- list()
    for (j in 1:length(Omega)) {
        poly.list <- list()
        for (tri in Omega[[j]]) {
            px <- mesh$graph$tv[tri, ]
            temp <- mesh$loc[px, ] # is a 3 by 3 matrix of node locations
            poly.list <-
                c(poly.list, sp::Polygon(rbind(temp[c(3,2,1), 1:2], temp[3, 1:2]), hole = FALSE))
        }
        mesh.polys <- sp::SpatialPolygons(list(sp::Polygons(poly.list, ID = "noid")))
        if (interactive() &&
            compareVersion(getNamespaceVersion("sf"), "1.0-10") < 0) {
            warning(
                paste0(
                    "The sf::st_as_sfc.SpatialPolygons() function is broken for some inputs in versions before 1.0-10.\n",
                    "Please check your output carefully, and upgrade to sf version 1.0-10 or later."),
                immediate. = TRUE
            )
        }
        inla.require("sf", stop.on.error = TRUE)
        mesh.polys <- sf::st_as_sfc(mesh.polys)
        Omega.SP.list[[j]] <- sf::as_Spatial(sf::st_union(mesh.polys))
    }

    if (missing(barrier.triangles)) {
        return(Omega.SP.list)
    } else {
        ## Only return the polygon describing the barrier
        return(Omega.SP.list[[2]])
    }
}

### PRECISION MATRIX FUNCTIONS ###
## I.e. Solve the differential equation (the SPDE)

#' @details * `inla.barrier.q`: This function computes a specific precision matrix
#' @param fem represents the Barrier model or the Different Terrains (DT) model,
#' by containing all the needed matrices to solve the SPDE
#' @param ranges,sigma the hyperparameters that determine Q
#' @param envir the environment used for caching (with optimize=TRUE), if any
#' @export
#' @rdname inla.barrier
`inla.barrier.q` <- function(fem, ranges, sigma = 1, envir = NULL) {
    if (is.null(ranges)) stop("ranges cannot be NULL")
    if (any(is.na(ranges))) stop("No range can be NA")

    xi <- length(ranges)
    if (xi != length(fem$D)) {
        print("inla.barrier.q has encountered an error. Will stop.")
        stop("Ranges do no correspond to fem")
    }
    if (any(ranges < 0.001)) {
        warning("This hyper parameter value may fail. A very small maximum edge length needed in the mesh.")
    }

    Cdiag <- ranges[1]^2 * fem$C[[1]] # already raised to power 2
    if (xi > 1) {
        for (k in 2:xi) {
            Cdiag <- Cdiag + ranges[k]^2 * fem$C[[k]]
        }
    }
    N <- length(Cdiag)
    Cinv <- sparseMatrix(i = 1:N, j = 1:N, x = 1.0 / Cdiag, dims = c(N, N), repr = "T")

    A <- fem$I
    for (k in 1:xi) {
        A <- A + (ranges[k]^2 / 8.0) * fem$D[[k]]
    }

    ## Build the precision matrix
    ## The last multiplication factor is because the C matrix
    ## is different by a factor of 3 compared to the corresponding matrix
    ## in the stationary spde approach
    ## (factor of 3 commented as C from inla.barrier.fem() is now as in the stationary case)
    Q <- inla.as.sparse(t(A) %*% Cinv %*% A * (1 / sigma^2) / pi * 2) ## * 3)

    ## if it exists, then this is the environment used for caching and optimize=TRUE.
    ## otherwise, use optimize=FALSE format
    if (is.environment(envir)) {
        if (!exists("cache.done", envir = envir)) {
            Qx.idx <- which(Q@i <= Q@j)
            assign("Qx.idx", Qx.idx, envir = envir)
            assign("cache.done", TRUE, envir = envir)
        } else {
            Qx.idx <- get("Qx.idx", envir = envir)
        }
        return(Q@x[Qx.idx])
    } else {
        return (Q)
    }
}

#' @details * `inla.barrier.fem` This function computes the Finite Element
#' matrices that are needed to compute the precision matrix Q later
#' @export
#' @rdname inla.barrier
`inla.barrier.fem` <- function(mesh, barrier.triangles, Omega = NULL) {
    stopifnot(inherits(mesh, "inla.mesh"))

    if (missing(barrier.triangles) && is.null(Omega)) stop("Input barrier triangles")

    if (missing(barrier.triangles)) {
        ## Use Omega
    } else {
        ## Create a valid Omega
        barrier.triangles <- unique(barrier.triangles)
        t <- length(mesh$graph$tv[, 1])
        ## - all triangles
        remaining <- setdiff(1:t, barrier.triangles)
        if (!is.null(Omega)) warning("Omega is replaced by barrier.triangles")
        Omega <- list(remaining, barrier.triangles)
    }

    ## This function computes the Finite Element matrix of the white noise
    ## on one subdomain (area)
    ## Returns a diagonal matrix, represented as a vector
    dt.fem.white <- function(mesh, subdomain) {

        ## Pre-allocation
        Ck <- rep(0, mesh$n)

        for (t in subdomain) {
            ## Node indexes for triangle t:
            px <- mesh$graph$tv[t, ]
            temp <- mesh$loc[px, ] # is a 3 by 3 matrix of node locations
            p1 <- t(t(temp[1, c(1, 2)]))
            p2 <- t(t(temp[2, c(1, 2)]))
            p3 <- t(t(temp[3, c(1, 2)]))
            Ts <- cbind(p2 - p1, p3 - p1) # is the transformation to reference triangle
            area <- abs(det(Ts)) * 0.5
            for (i in 1:3) {
                Ck[px[i]] <- Ck[px[i]] + area
            }
        }
        return(Ck/3) ## divide by 3 to match with C in inla.mesh.fem()
    }

    ## This function computes the Finite Element matrix for the '1' in
    ## the SPDE (the identity operator) <phi, phi>
    ## Result does not depend on the subdomains
    dt.fem.identity <- function(mesh) {
        ## Preallocation
        len <- length(mesh$graph$tv[, 1])
        index.i <- rep(0, len * 6)
        index.j <- rep(0, len * 6)
        Aij <- rep(0, len * 6) # matrix values
        counter <- 1

        for (t in 1:len) {
            ## Node indexes for triangle t:
            px <- mesh$graph$tv[t, ]
            temp <- mesh$loc[px, ] # is a 3 by 3 matrix of node locations
            p1 <- t(t(temp[1, c(1, 2)]))
            p2 <- t(t(temp[2, c(1, 2)]))
            p3 <- t(t(temp[3, c(1, 2)]))
            Ts <- cbind(p2 - p1, p3 - p1) # is the transformation to reference triangle
            twiceArea <- abs(det(Ts))

            for (i in 1:3) {
                index.i[counter] <- px[i]
                index.j[counter] <- px[i] # same
                Aij[counter] <- (twiceArea) * 1 / 12
                counter <- counter + 1
            }

            for (i in 1:2) {
                for (j in (i + 1):3) {
                    index.i[counter] <- px[i]
                    index.j[counter] <- px[j]
                    Aij[counter] <- (twiceArea) * 1 / 24
                    counter <- counter + 1
                    ## symmetry:
                    index.i[counter] <- px[j]
                    index.j[counter] <- px[i]
                    Aij[counter] <- (twiceArea) * 1 / 24
                    counter <- counter + 1
                }
            }
        }

        I <- sparseMatrix(i = index.i, j = index.j, x = Aij,
                          dims = c(mesh$n, mesh$n), repr = "T")
        return(I)
    }

    ## This function computes the Finite Element matrix of the laplace operator
    ## on one subdomain (area)
    ## The resulting matrix is very sparse
    dt.fem.laplace <- function(mesh, subdomain) {
        ## The nabla phi's
        Nphix <- rbind(c(-1, -1), c(1, 0), c(0, 1))
        len <- length(subdomain)
        index.i <- rep(0, len * 9)
        index.j <- rep(0, len * 9)
        Aij <- rep(0, len * 9) # matrix values
        counter <- 1

        for (tri in subdomain) {
            px <- mesh$graph$tv[tri, ]
            temp <- mesh$loc[px, ] # is a 3 by 3 matrix of node locations
            p1 <- t(t(temp[1, c(1, 2)]))
            p2 <- t(t(temp[2, c(1, 2)]))
            p3 <- t(t(temp[3, c(1, 2)]))
            Ts <- cbind(p2 - p1, p3 - p1) # is the transformation to reference triangle
            TTTinv <- solve(t(Ts) %*% Ts)
            area <- abs(det(Ts)) * 0.5

            for (k in 1:3) {
                for (m in 1:3) {
                    tmp <- (3 * m + k - 4) * length(subdomain)
                    index.i[(tmp + counter)] <- px[k]
                    index.j[(tmp + counter)] <- px[m]
                    Aij[(tmp + counter)] <- area * Nphix[k, c(1, 2)] %*% TTTinv %*% as.matrix(Nphix[m, c(1, 2)])
                }
            }
            counter <- counter + 1
        }

        Dk <- sparseMatrix(i = index.i, j = index.j, x = Aij, dims = c(mesh$n, mesh$n), repr = "T")
        return(Dk)
    }

    ## Call all the functions to get the resulting list of fem matrices

    xi <- length(Omega)
    
    if (requireNamespace("INLAspacetime")) {
        warning("Using implementation from the `INLAspacetime` package")
        fem <- INLAspacetime::mesh2fem.barrier(
                                  mesh = mesh, 
                                  barrier.triangles = Omega[[2L]])
    } else {
        warning(paste(
            "Please install the `INLAspacetime` package\n",
            "which contains an implementation that runs faster!"))
        fem <- list()
        fem$I <- dt.fem.identity(mesh)
        fem$D <- list()
        fem$C <- list()
        for (k in 1:xi) {
            fem$D[[k]] <- dt.fem.laplace(mesh, Omega[[k]])
        }
        for (k in 1:xi) {
            fem$C[[k]] <- dt.fem.white(mesh, Omega[[k]])
        }
        fem$hdim <- xi
    }

    return(fem)
}
