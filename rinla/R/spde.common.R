
#' Build a block-diagonal sparse matrix.
#'
#' Build a block-diagonal sparse matrix.  Obsolete wrapper for `bdiag()`.
#'
#'
#' @aliases inla.dBind inla.dBind-deprecated
#' @param \dots A list of square or rectangular matrices.
#' @return A sparse matrix.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @keywords internal
#' @examples
#' \dontrun{
#'   inla.dBind(Matrix(1, 2, 1), Matrix(2, 1, 2))
#' }
#' bdiag(Matrix(1, 2, 1), Matrix(2, 1, 2))
#' @export inla.dBind
inla.dBind <- function(...) {
    .Deprecated("Matrix::bdiag")
    return(bdiag(...))
}



#' Extract elements by matching name from container objects.
#'
#' Extract elements by wildcard name matching from a `data.frame`,
#' `list`, or `matrix`.
#'
#'
#' @aliases inla.extract.el inla.extract.el.data.frame inla.extract.el.list
#' inla.extract.el.matrix
#' @param M A container object.
#' @param match A regex defining the matching criterion.
#' @param by.row If `TRUE`, extract data by row, otherwise by column.
#' @param \dots Additional arguments, not used.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @export inla.extract.el
inla.extract.el <- function(M, ...) {
    if (is.null(M)) {
          return(NULL)
      }
    UseMethod("inla.extract.el", M)
}

inla.regex.match <- function(x, match) {
    return(strsplit(x, match)[[1]][1] == "")
}

#' @export
#' @rdname inla.extract.el
inla.extract.el.matrix <- function(M, match, by.row = TRUE, ...) {
    if (by.row) {
        return(M[sapply(rownames(M), inla.regex.match, match = match), , drop = FALSE])
    } else {
        return(M[, sapply(colnames(M), inla.regex.match, match = match), drop = FALSE])
    }
}

#' @export
#' @rdname inla.extract.el
inla.extract.el.data.frame <- function(M, match, by.row = TRUE, ...) {
    if (by.row) {
        return(M[sapply(rownames(M), inla.regex.match, match = match), , drop = FALSE])
    } else {
        return(M[, sapply(colnames(M), inla.regex.match, match = match), drop = FALSE])
    }
}

#' @export
#' @rdname inla.extract.el
inla.extract.el.list <- function(M, match, ...) {
    return(M[sapply(names(M), inla.regex.match, match = match)])
}



inla.spde.homogenise_B_matrix <- function(B, n.spde, n.theta) {
    if (!is.numeric(B)) {
          stop("B matrix must be numeric.")
      }
    if (is.matrix(B)) {
        if ((nrow(B) != 1) && (nrow(B) != n.spde)) {
            stop(inla.paste(list("B matrix has",
                as.character(nrow(B)),
                "rows but should have 1 or",
                as.character(n.spde),
                sep = " "
            )))
        }
        if ((ncol(B) != 1) && (ncol(B) != 1 + n.theta)) {
            stop(inla.paste(list("B matrix has",
                as.character(ncol(B)),
                "columns but should have 1 or",
                as.character(1 + n.theta),
                sep = " "
            )))
        }
        if (ncol(B) == 1) {
            return(cbind(as.vector(B), matrix(0.0, n.spde, n.theta)))
        } else if (ncol(B) == 1 + n.theta) {
            if (nrow(B) == 1) {
                return(matrix(as.vector(B), n.spde, 1 + n.theta, byrow = TRUE))
            } else if (nrow(B) == n.spde) {
                return(B)
            }
        }
    } else { ## !is.matrix(B)
        if ((length(B) == 1) || (length(B) == n.spde)) {
            return(cbind(B, matrix(0.0, n.spde, n.theta)))
        } else if (length(B) == 1 + n.theta) {
            return(matrix(B, n.spde, 1 + n.theta, byrow = TRUE))
        } else {
            stop(inla.paste(list(
                "Length of B vector is",
                as.character(length(B)),
                "but should be 1,",
                as.character(1 + n.theta), "or",
                as.character(n.spde)
            ),
            sep = " "
            ))
        }
    }
    stop(inla.paste(list("Unrecognised structure for B matrix"),
        sep = " "
    ))
}





#' Numerical evaluation of Matern and related covariance functions.
#'
#' Calculates covariance and correlation functions for Matern models and
#' related oscillating SPDE models, on \eqn{R^d}{R^d} and on the sphere,
#' \eqn{S^2}{S^2}.
#'
#' On \eqn{R^d}{R^d}, the models are *defined* by the spectral density
#' given by
#' \deqn{S(w) = \frac{1}{(2\pi)^d (\kappa^4 + 2 \kappa^2 \cos(\pi
#' \theta) |w|^2 + |w|^4)^{(\nu + d/2)/2}}
#' }{S(w) = 1 / ( (2\pi)^d * (kappa^4 + 2 kappa^2 * cos(pi * theta) * |w|^2 +
#' |w|^4)^((nu + d/2)/2) )}
#'
#' On \eqn{S^2}{S^2}, the models are *defined* by the spectral
#' coefficients
#' \deqn{S(k) = \frac{2k+1}{4\pi (\kappa^4 + 2 \kappa^2 \cos(\pi
#' \theta) k(k+1) + k^2(k+1)^2)^{(\nu +
#' 1)/2}}
#' }{S(k) = (2k+1) / (4 pi (kappa^4 + 2 kappa^2 cos(pi theta) k(k+1) +
#' k^2(k+1)^2)^((\nu + 1)/2) )}
#'
#' @aliases inla.matern.cov inla.matern.cov.s2
#' @param nu The Matern smoothness parameter.
#' @param kappa The spatial scale parameter.
#' @param x Distance values.
#' @param d Space dimension; the domain is \eqn{R^d}{R^d}.
#' @param corr If `TRUE`, calculate correlations, otherwise calculate
#' covariances.  Only used for pure Matern models (i.e. with
#' \eqn{\theta=0}{theta=0}).
#' @param norm.corr If `TRUE`, normalise by the estimated variance, giving
#' approximate correlations.
#' @param theta Oscillation strength parameter.
#' @param epsilon Tolerance for detecting points close to distance zero.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @export inla.matern.cov
inla.matern.cov <- function(nu, kappa, x, d = 1, corr = FALSE, norm.corr = FALSE, theta, epsilon = 1e-8) {
    if (missing(theta)) { ## Ordinary Matern
        y <- kappa * abs(x)
        if (corr) {
            ok <- (y >= epsilon)
            if (nu <= 0) {
                covariance <- y * 0
                covariance[!ok] <- 1 - y / epsilon
            } else {
                covariance <- y * 0
                covariance[ok] <-
                    2^(1 - nu) / gamma(nu) * (y[ok])^nu * besselK(y[ok], nu)
                if (any(!ok)) {
                    scale <- 2^(1 - nu) / gamma(nu)
                    ## corr = scale y^nu K_nu(y)
                    ## = 1 - b y^(2 nu) + o(y^(2 nu), 0 < \nu < 1
                    ## = 1 - b y^2 + o(y^2), 1 <= \nu
                    ## (1-corr(eps)/vari)/eps^(2nu) = b
                    if (nu < 1) {
                        exponent <- 2 * nu
                    } else {
                        exponent <- 2
                    }
                    corr.eps <-
                        scale * epsilon^nu * besselK(epsilon, nu)
                    b <- (1 - corr.eps) / epsilon^exponent
                    covariance[!ok] <- (1 - b * y[!ok]^exponent)
                }
            }
            return(covariance)
        } else {
            ok <- (y >= epsilon)
            covariance <- y * 0
            covariance[ok] <-
                2^(1 - nu) / gamma(nu + d / 2) / (4 * pi)^(d / 2) / kappa^(2 * nu) *
                    (y[ok])^nu * besselK(y[ok], nu)
            if (any(!ok)) {
                if (nu > 0) { ## Regular Matern case
                    vari <- gamma(nu) / gamma(nu + d / 2) / (4 * pi)^(d / 2) / kappa^(2 * nu)
                    scale <- 2^(1 - nu) / gamma(nu)
                    ## corr = scale y^nu K_nu(y)
                    ## = 1 - b y^(2 nu) + o(y^(2 nu), 0 < \nu < 1
                    ## = 1 - b y^2 + o(y^2), 1 <= \nu
                    ## (1-corr(eps)/vari)/eps^(2nu) = b
                    if (nu < 1) {
                        exponent <- 2 * nu
                    } else {
                        exponent <- 2
                    }
                    corr.eps <-
                        scale * epsilon^nu * besselK(epsilon, nu)
                    b <- (1 - corr.eps) / epsilon^exponent
                    covariance[!ok] <- vari * (1 - b * y[!ok]^exponent)
                } else if (nu == 0) { ## Limiting Matern case
                    g <- 0.577215664901484 ## Euler's constant
                    covariance[!ok] <-
                        2 / gamma(d / 2) / (4 * pi)^(d / 2) *
                            (-log(y[!ok] / 2) - g)
                } else { ## (nu<0)
                    ## TODO: check this...
                    covariance[!ok] <-
                        ((2^(1 - nu) / gamma(nu + d / 2) / (4 * pi)^(d / 2) / kappa^(2 * nu) *
                            gamma(nu) * 2^(nu - 1)) * (1 - (y[!ok] / epsilon)) +
                            (2^(1 - nu) / gamma(nu + d / 2) / (4 * pi)^(d / 2) / kappa^(2 * nu) *
                                epsilon^nu * besselK(epsilon, nu)) * (y[!ok] / epsilon))
                }
            }
            return(covariance)
        }
    } else { ## Oscillating covariances
        y <- abs(x)
        if (d > 2L) {
            warning("Dimension > 2 not implemented for oscillating models.")
        }
        freq.max <- 1000 / max(y)
        freq.n <- 10000
        w <- seq(0, freq.max, length.out = freq.n)
        dw <- w[2] - w[1]
        spec <- 1 / (2 * pi)^d / (kappa^4 + 2 * kappa^2 * cos(pi * theta) * w^2 + w^4)^((nu + d / 2) / 2)
        if (d == 1L) {
            covariance <- y * 0 + spec[1] * dw
        } else {
            covariance <- y * 0
        }
        for (k in 2:freq.n) {
            if (d == 1L) {
                covariance <- covariance + 2 * cos(y * w[k]) * spec[k] * dw
            } else {
                covariance <- covariance + w[k] * besselJ(y * w[k], 0) * spec[k] * dw
            }
        }

        if (norm.corr) {
            noise.variance <- 1 / covariance[1]
        } else {
            noise.variance <- 1
        }

        return(covariance * noise.variance)
    }
}

#' @export
#' @param freq.max The maximum allowed harmonic order. Current default 40, to
#'   be changed to a dynamic choice based on error bounds.
#' @rdname inla.matern.cov
inla.matern.cov.s2 <- function(nu, kappa, x, norm.corr = FALSE, theta = 0,
                               freq.max = NULL) {
    inla.require("gsl", stop.on.error = TRUE)
    y <- cos(abs(x))
    
    # TODO: use the error bounds to pick a freq.max; ok since gsl Legendre
    # polynomial evaluations are stable to high order
    if (is.null(freq.max)) {
        freq.max <- 40L
    }
    freq.n <- freq.max + 1L
    w <- 0L:freq.max
    spec <- 1 / (kappa^4 +
                     2 * kappa^2 * cos(pi * theta) * w * (w + 1) +
                     w^2 * (w + 1)^2)^((nu + 1) / 2) / (4 * pi)
    covariance <- y * 0
    # TODO: blockwise calculations with legendre_Pl_array
    for (k in w) {
        covariance <- (covariance + (2 * k + 1) * spec[k + 1] *
                           gsl::legendre_Pl(l = k, x = y))
    }
    
    if (norm.corr) {
        noise.variance <- 1 / covariance[1]
    } else {
        noise.variance <- 1
    }
    
    return(covariance * noise.variance)
}

# @export
# @details `inla.matern.cov.s2.unstable` is an old implementation using
# `orthopolynom::legendre.polynomials` and `orthopolynom::polynomial.values`
# @rdname inla.matern.cov
# inla.matern.cov.s2.unstable <- function(nu, kappa, x, norm.corr = FALSE, theta = 0,
#                                freq.max = NULL) {
#     # Example to show how much better the new version is:
#     # library(tidyverse)
#     # library(ggplot2)
#     # ggplot(
#     #     data.frame(x=rep(seq(0,0.001,length.out=1000),times=2),
#     #                Type = rep(c("Unstable","Stable"), times=2)) %>%
#     #         mutate(Cov=case_when(
#     #             Type=="Unstable" ~ inla.matern.cov.s2.unstable(1,1,x,freq.max=40),
#     #             TRUE ~ inla.matern.cov.s2(1,1,x,freq.max=40)))
#     # ) +
#     #   geom_line(aes(x,Cov,col=Type))
#     
#     inla.require("orthopolynom", stop.on.error = TRUE)
#     y <- cos(abs(x))
#     
#     # TODO: use the error bounds to pick a freq.max; ok since gsl Legendre
#     # polynomial evaluations are stable to high order
#     if (is.null(freq.max)) {
#         freq.max <- 40L
#     }
#     freq.n <- freq.max + 1L
#     w <- 0L:freq.max
#     spec <- 1 / (kappa^4 +
#                      2 * kappa^2 * cos(pi * theta) * w * (w + 1) +
#                      w^2 * (w + 1)^2)^((nu + 1) / 2) / (4 * pi)
#     covariance <- y * 0
#     leg <- orthopolynom::legendre.polynomials(freq.max)
#     for (k in w) {
#         covariance <- (covariance + (2 * k + 1) * spec[k + 1] *
#                            orthopolynom::polynomial.values(leg[k + 1], y)[[1]])
#     }
#     
#     if (norm.corr) {
#         noise.variance <- 1 / covariance[1]
#     } else {
#         noise.variance <- 1
#     }
#     
#     return(covariance * noise.variance)
# }





#' List SPDE models supported by inla.spde objects
#'
#' List SPDE models supported by inla.spde objects
#'
#' Returns a list of available SPDE model type name lists, one for each
#' inla.spde model class (currently [inla.spde1()] and
#' [inla.spde2()]).
#'
#' @aliases inla.spde.models inla.spde1.models inla.spde2.models
#' @param function.names If `FALSE`, return list model name lists.  If
#' `TRUE`, return list of model object constructor function names.
#' @return List of available SPDE model type name lists.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @examples
#' \dontrun{
#' ## Display help for each supported inla.spde2 model:
#' for (model in inla.spde2.models()) {
#'     print(help(paste("inla.spde2.", model, sep = "")))
#' }
#'
#' ## Display help for each supported inla.spde* model:
#' models <- inla.spde.models()
#' for (type in names(models)) {
#'    for (model in models[[type]]) {
#'        print(help(paste("inla.", type, ".", model, sep = "")))
#'    }
#' }
#'
#' ## Display help for each supported inla.spde* model (equivalent to above):
#' for (model in inla.spde.models(function.names = TRUE)) {
#'     print(help(model))
#' }
#' }
#' @export inla.spde.models
inla.spde.models <- function(function.names = FALSE) {
    types <- c("spde1", "spde2")
    models <- list()
    for (t in types) {
        models[[t]] <-
            do.call(
                what = paste("inla.", t, ".models", sep = ""),
                args = list()
            )
        if (function.names) {
            models[[t]] <- paste("inla.", t, ".", models[[t]], sep = "")
        }
    }

    if (function.names) {
        models <- as.vector(do.call(c, models))
    }

    return(models)
}




#' Sample from SPDE models
#'
#' Old methods fo sampling from a SPDE model.  For new code, use
#' [inla.spde.precision()] and [inla.qsample()] instead.
#'
#'
#' @aliases inla.spde.sample inla.spde.sample.default
#' inla.spde.sample.inla.spde
#' @param precision A precision matrix.
#' @param seed The seed for the pseudo-random generator.
#' @param spde An `inla.spde` object.
#' @param \dots Parameters passed on to other methods.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.spde.precision()], [inla.qsample()]
#' @export inla.spde.sample
inla.spde.sample <- function(...) {
    warning("inla.spde.sample is deprecated.  Please use inla.qsample() in combination with inla.spde.precision() instead.")
    UseMethod("inla.spde.sample")
}

#' @export
#' @rdname inla.spde.sample
inla.spde.sample.default <-
    function(precision, seed = NULL, ...) {
        return(inla.finn(precision,
            seed = (inla.ifelse(
                is.null(seed),
                0L,
                seed
            ))
        )$sample)
    }

#' @export
#' @rdname inla.spde.sample
inla.spde.sample.inla.spde <-
    function(spde, seed = NULL, ...) {
        precision <- inla.spde.precision(spde, ...)
        return(inla.spde.sample(precision, seed = seed))
    }





#' @title Precision matrices for SPDE models
#'
#' @description Calculates the precision matrix for given parameter
#'  values based on an `inla.spde` model object.
#'
#' @aliases inla.spde.precision inla.spde1.precision inla.spde2.precision
#' inla.spde.precision.inla.spde1 inla.spde.precision.inla.spde2
#' @param spde An `inla.spde` object.
#' @param theta The parameter vector.
#' @param phi0 Internal parameter for a generic model.  Expert option only.
#' @param phi1 Internal parameter for a generic model.  Expert option only.
#' @param phi2 Internal parameter for a generic model.  Expert option only.
#' @param \dots Additional parameters passed on to other methods.
#' @return A sparse precision matrix.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.spde.models()], [inla.spde2.generic()],
#' [inla.spde2.theta2phi0()], [inla.spde2.theta2phi1()],
#' [inla.spde2.theta2phi2()]
#' @export
#' @rdname inla.spde.precision
inla.spde.precision <- function(...) {
    UseMethod("inla.spde.precision")
}



#' SPDE result extraction from INLA estimation results
#'
#' Exctract field and parameter values and distributions for an
#' `inla.spde` SPDE effect from an inla result object.
#'
#'
#' @aliases inla.spde.result inla.spde.result.inla.spde1
#' inla.spde.result.inla.spde2 inla.spde1.result inla.spde2.result
#' @param inla An `inla` object obtained from a call to [inla()]
#' @param name A character string with the name of the SPDE effect in the inla
#' formula.
#' @param spde The `inla.spde` object used for the effect in the inla
#' formula. (Note: this could have been stored in the inla output, but isn't.)
#' Usually the result of a call to [inla.spde2.matern()].
#' @param do.transform If `TRUE`, also calculate marginals transformed to
#' user-scale.  Setting to `FALSE` is useful for large non-stationary
#' models, as transforming many marginal densities is time-consuming.
#' @param \dots Further arguments passed to and from other methods.
#' @return For `inla.spde2` models, a list, where the nominal range and
#' variance are defined as the values that would have been obtained with a
#' stationary model and no boundary effects: \item{marginals.kappa }{Marginal
#' densities for kappa} \item{marginals.log.kappa }{Marginal densities for
#' log(kappa)} \item{marginals.log.range.nominal }{Marginal densities for
#' log(range)} \item{marginals.log.tau }{Marginal densities for log(tau)}
#' \item{marginals.log.variance.nominal }{Marginal densities for log(variance)}
#' \item{marginals.range.nominal }{Marginal densities for range}
#' \item{marginals.tau }{Marginal densities for tau} \item{marginals.theta
#' }{Marginal densities for the theta parameters} \item{marginals.values
#' }{Marginal densities for the field values} \item{marginals.variance.nominal
#' }{Marginal densities for variance} \item{summary.hyperpar }{The SPDE related
#' part of the inla hyperpar output summary} \item{summary.log.kappa }{Summary
#' statistics for log(kappa)} \item{summary.log.range.nominal }{Summary
#' statistics for log(range)} \item{summary.log.tau }{Summary statistics for
#' log(tau)} \item{summary.log.variance.nominal }{Summary statistics for
#' log(kappa)} \item{summary.theta }{Summary statistics for the theta
#' parameters} \item{summary.values }{Summary statistics for the field values}
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.spde.models()], [inla.spde2.matern()]
#' @examples
#'
#' loc <- matrix(runif(100 * 2), 100, 2)
#' mesh <- fmesher::fm_mesh_2d_inla(loc.domain = loc, max.edge = c(0.1, 0.5))
#' spde <- inla.spde2.matern(mesh)
#' index <- inla.spde.make.index("spatial", mesh$n, n.repl = 2)
#' spatial.A <- inla.spde.make.A(mesh, loc,
#'     index = rep(1:nrow(loc), 2),
#'     repl = rep(1:2, each = nrow(loc))
#' )
#' ## Toy example with no spatial correlation (range=zero)
#' y <- 10 + rnorm(100 * 2)
#' stack <- inla.stack(
#'     data = list(y = y),
#'     A = list(spatial.A),
#'     effects = list(c(index, list(intercept = 1))),
#'     tag = "tag"
#' )
#' data <- inla.stack.data(stack, spde = spde)
#' formula <- y ~ -1 + intercept + f(spatial,
#'     model = spde,
#'     replicate = spatial.repl
#' )
#' result <- inla(formula,
#'     family = "gaussian", data = data,
#'     control.predictor = list(A = inla.stack.A(stack))
#' )
#' spde.result <- inla.spde.result(result, "spatial", spde)
#' plot(spde.result$marginals.range.nominal[[1]], type = "l")
#' @export inla.spde.result
inla.spde.result <- function(...) {
    inla.require.inherits(list(...)[[1]], "inla", "First parameter")
    inla.require.inherits(list(...)[[2]], "character", "Second parameter")
    UseMethod("inla.spde.result", list(...)[[3]])
}









#' SPDE model index vector generation
#'
#' Generates a list of named index vectors for an SPDE model.
#'
#'
#' @param name A character string with the base name of the effect.
#' @param n.spde The size of the model, typically from `spde$n.spde`.
#' @param n.group The size of the `group` model.
#' @param n.repl The number of model replicates.
#' @param ...  Additional parameters.  Currently unused.
#' @return A list of named index vectors. \item{name }{Indices into the vector
#' of latent variables} \item{name.group }{'group' indices} \item{name.repl
#' }{Indices for replicates}
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.spde.make.A()], [inla.spde2.result()]
#' @examples
#'
#' loc <- matrix(runif(100 * 2), 100, 2)
#' mesh <- fmesher::fm_mesh_2d_inla(loc.domain = loc, max.edge = c(0.1, 0.5))
#' spde <- inla.spde2.matern(mesh)
#' index <- inla.spde.make.index("spatial", spde$n.spde, n.repl = 2)
#' spatial.A <- inla.spde.make.A(mesh, loc,
#'     index = rep(1:nrow(loc), 2),
#'     repl = rep(1:2, each = nrow(loc))
#' )
#' y <- 10 + rnorm(100 * 2)
#' stack <- inla.stack(
#'     data = list(y = y),
#'     A = list(spatial.A),
#'     effects = list(c(index, list(intercept = 1))),
#'     tag = "tag"
#' )
#' data <- inla.stack.data(stack, spde = spde)
#' formula <- y ~ -1 + intercept + f(spatial,
#'     model = spde,
#'     replicate = spatial.repl
#' )
#' result <- inla(formula,
#'     family = "gaussian", data = data,
#'     control.predictor = list(A = inla.stack.A(stack))
#' )
#' spde.result <- inla.spde2.result(result, "spatial", spde)
#' @export inla.spde.make.index
inla.spde.make.index <- function(name, n.spde, n.group = 1, n.repl = 1, ...) {
    if ("n.mesh" %in% names(list(...))) {
        warning("'n.mesh' is deprecated, please use 'n.spde' instead.")
        if (missing(n.spde) || is.null(n.spde)) {
              n.spde <- list(...)$n.mesh
          }
    }

    name.group <- paste(name, ".group", sep = "")
    name.repl <- paste(name, ".repl", sep = "")
    out <- list()
    out[[name]] <- rep(rep(1:n.spde, times = n.group), times = n.repl)
    out[[name.group]] <- rep(rep(1:n.group, each = n.spde), times = n.repl)
    out[[name.repl]] <- rep(1:n.repl, each = n.spde * n.group)
    return(out)
}







#' Row-wise Kronecker products
#' 
#' `r lifecycle::badge("deprecated")` in favour of `fmesher::fm_row_kron()`, which
#' is typically an order of magnitude faster than the old `inla.row.kron()` implementation.
#'  `inla.row.kron()` now calls `fm_row_kron()` internally instead.
#'
#' Takes two Matrices and computes the row-wise Kronecker product.  Optionally
#' applies row-wise weights and/or applies an additional 0/1 row-wise Kronecker
#' matrix product, as needed by [inla.spde.make.A()].
#'
#' @param M1 A matrix that can be transformed into a sparse Matrix.
#' @param M2 A matrix that can be transformed into a sparse Matrix.
#' @param repl An optional index vector.  For each entry, specifies which
#' replicate the row belongs to, in the sense used in
#' [inla.spde.make.A()].
#' @param n.repl The maximum replicate index, in the sense used in
#' [inla.spde.make.A()].
#' @param weights Optional scaling weights to be applied row-wise to the
#' resulting matrix.
#' @return A `sparseMatrix` object.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.spde.make.A()]
#' @export inla.row.kron
#' @keywords internal
inla.row.kron <- function(M1, M2, repl = NULL, n.repl = NULL, weights = NULL) {
    lifecycle::deprecate_soft(
      when = "24.05.28",
      what = "inla.row.kron()",
      with = "fmesher::fm_row_kron()",
      details = "fmesher::fm_row_kron() is typically an order of magnitude faster than the old inla.row.kron() implementation.")
    fmesher::fm_row_kron(M1,
                         M2,
                         repl = repl,
                         n.repl = n.repl,
                         weights = weights)
}






#' Observation matrices for mesh models.
#'
#' Constructs observation/prediction weight matrices for numerical integration
#' schemes for regional data problems.  Primarily intended for internal use by
#' [inla.spde.make.A()].
#'
#'
#' @param A A precomputed observation/prediction matrix for locations that are
#' to be joined.
#' @param block Indices specifying block groupings: Entries with the same
#' `block` value are joined into a single row in the resulting matrix, and
#' the `block` values are the row indices.
#' @param n.block The number of blocks.
#' @param weights Optional scaling weights to be applied row-wise to the input
#' `A` matrix.
#' @param rescale Specifies what scaling method should be used when joining the
#' rows of the `A` matrix as grouped by the `block` specification.
#' \itemize{ \item `'none'`: Straight sum, no rescaling.  \item
#' `'count'`: Divide by the number of entries in the block.  \item
#' `'weights'`: Divide by the sum of the `weight` values within each
#' block.  \item `'sum'`: Divide by the resulting row sums.  }
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.spde.make.A()]
#' @export inla.spde.make.block.A
inla.spde.make.block.A <-
    function(A,
             block,
             n.block = max(block),
             weights = NULL,
             rescale = c("none", "count", "weights", "sum")) {
        ## Add A-matrix rows belonging to the same "block" with optional
        ## weights and optional rescaling, by dividing row i by a_i, for |a_i|>0:
        ## B(i) = {j; block(i)==block(j) and sum_k |A_jk| > 0 }
        ## "count":   a_i = #{j in B(i)} }
        ## "weights": a_i = \sum_{j in B(i)} weights_j
        ## "sum":     a_i = \sum_{j in B(i)} \sum_k A_jk weights_j
        ##
        ## This function makes use of the feature of sparseMatrix to sum all
        ## values for multiple instances of the same (i,j) pair.

        A <- inla.as.dgTMatrix(A)
        N <- nrow(A)
        ## length(block) should be == N or 1
        if (length(block) == 1L) {
            block <- rep(block, N)
        }
        if (is.null(weights)) {
            weights <- rep(1, N)
        }

        rescale <- match.arg(rescale)

        if (!(rescale == "none")) {
            if (rescale == "count") {
                ## Count the non-zero rows within each block
                sums <- (sparseMatrix(
                    i = block,
                    j = rep(1L, N),
                    x = (rowSums(abs(A)) > 0) * 1.0,
                    dims = c(n.block, 1L)
                ))[block]
            } else if (rescale == "weights") {
                ## Sum the weights within each block
                sums <- (sparseMatrix(
                    i = block,
                    j = rep(1L, N),
                    x = (rowSums(abs(A)) > 0) * weights,
                    dims = c(n.block, 1L)
                ))[block]
            } else { ## (rescale == "sum"){
                ## Sum the weighted values within each block
                sums <- (sparseMatrix(
                    i = block,
                    j = rep(1L, N),
                    x = rowSums(A) * weights,
                    dims = c(n.block, 1L)
                ))[block]
            }
            ## Normalise:
            ok <- (abs(sums) > 0)
            weights[ok] <- weights[ok] / sums[ok]
        }

        return(inla.as.dgTMatrix(sparseMatrix(
            i = block[1L + A@i],
            j = 1L + A@j,
            x = A@x * weights[1L + A@i],
            dims = c(n.block, ncol(A))
        )))
    }





#' Observation/prediction matrices for mesh models.
#'
#' Constructs observation/prediction weight matrices for models based on
#' [inla.mesh()] and [inla.mesh.1d()] objects.
#'
#'
#' @param mesh An [inla.mesh()] or [inla.mesh.1d()] object
#' specifying a function basis on a mesh domain.  Alternatively, an
#' `inla.spde` object that includes a mesh (e.g. from
#' [inla.spde2.matern()]).
#' @param loc Observation/prediction coordinates.  `mesh` and `loc`
#' defines a matrix `A.loc` of mapping weights between basis function
#' weights and field values.  If `loc` is `NULL`, `A.loc` is
#' defined as `Diagonal(n.spde, 1)`.
#' @param index For each observation/prediction value, an index into
#' `loc`.  Default is `seq_len(nrow(A.loc))`.
#' @param group For each observation/prediction value, an index into the group
#' model.
#' @param repl For each observation/prediction value, the replicate index.
#' @param n.spde The number of basis functions in the mesh model. (Note: may be
#' different than the number of mesh vertices/nodes/knots.)
#' @param n.group The size of the group model.
#' @param n.repl The total number of replicates.
#' @param group.mesh An optional [inla.mesh.1d()] object for the
#' group model.
#' @param weights Optional scaling weights to be applied row-wise to the
#' resulting matrix.
#' @param A.loc Optional precomputed observation/prediction matrix.
#' `A.loc` can be specified instead of `mesh`+`loc`, optionally
#' with `index` supplied.
#' @param A.group Optional precomputed observation/prediction matrix for the
#' group model. `A.group` can be specified instead of `group` and/or
#' `group.mesh`, optionally with `group.index` supplied.
#' @param group.index For each observation/prediction value, an index into the
#' rows of `A.group`.
#' @param block Optional indices specifying block groupings: Entries with the
#' same `block` value are joined into a single row in the resulting
#' matrix, and the `block` values are the row indices.  This is intended
#' for construction of approximate integration schemes for regional data
#' problems. See [inla.spde.make.block.A()] for details.
#' @param n.block The number of blocks.
#' @param block.rescale Specifies what scaling method should be used when
#' joining entries as grouped by a `block` specification.  See
#' [inla.spde.make.block.A()] for details.
#' @param \dots Additional parameters.  Currently unused.
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [inla.spde.make.index()]
#' @examples
#'
#' loc <- matrix(runif(10000 * 2) * 1000, 10000, 2)
#' mesh <- inla.mesh.2d(
#'     loc = loc,
#'     cutoff = 50,
#'     max.edge = c(50, 500)
#' )
#' A <- inla.spde.make.A(mesh, loc = loc)
#' @export inla.spde.make.A
inla.spde.make.A <-
    function(mesh = NULL,
             loc = NULL,
             index = NULL,
             group = NULL,
             repl = 1L,
             n.spde = NULL,
             n.group = NULL,
             n.repl = NULL,
             group.mesh = NULL,
             weights = NULL,
             A.loc = NULL,
             A.group = NULL,
             group.index = NULL,
             block = NULL,
             n.block = NULL,
             block.rescale = c("none", "count", "weights", "sum"),
             ...)
    ## Deprecated/obsolete parameters: n.mesh, group.method
    {
        ## A.loc can be specified instead of mesh+loc, optionally with
        ## index supplied.
        ## A.group can be specified instead of group and/or group.mesh,
        ## optionally with group.index supplied.

        if ("n.mesh" %in% names(list(...))) {
            warning("'n.mesh' is deprecated, use 'n.spde' instead.")
            n.spde <- list(...)$n.mesh
        }
        if ("group.method" %in% names(list(...))) {
            group.method <-
                match.arg(list(...)$group.method, c("nearest", "S0", "S1"))
            warning(paste("'group.method=", group.method,
                "' is deprecated.  Specify 'degree=",
                switch(group.method, nearest = "0", S0 = "0", S1 = "1"),
                "' in inla.mesh.1d() instead.",
                sep = ""
            ))
        }

        if (is.null(mesh)) {
            if (is.null(A.loc) && is.null(n.spde)) {
                  stop("At least one of 'mesh', 'n.spde', and 'A.loc' must be specified.")
              }
            if (!is.null(A.loc)) {
                n.spde <- ncol(A.loc)
            }
        } else if (inherits(mesh, "inla.spde")) {
            spde <- mesh
            mesh <- spde$mesh
            inla.require.inherits(spde$mesh, c("inla.mesh", "inla.mesh.1d"), "'spde$mesh'")
            n.spde <- spde$n.spde
        } else {
            inla.require.inherits(mesh, c("inla.mesh", "inla.mesh.1d"), "'mesh'")
            if (inherits(mesh, "inla.mesh.1d")) {
                n.spde <- mesh$m
            } else {
                n.spde <- mesh$n
            }
        }
        if (!is.null(group.mesh)) {
            inla.require.inherits(group.mesh, "inla.mesh.1d", "'mesh'")
        }

        n.group <-
            ifelse(!is.null(n.group),
                n.group,
                ifelse(!is.null(A.group),
                    nrow(A.group),
                    ifelse(!is.null(group.mesh),
                        group.mesh$m,
                        max(
                            1,
                            ifelse(is.null(group),
                                1,
                                ifelse(length(group) == 0,
                                    1,
                                    max(group)
                                )
                            )
                        )
                    )
                )
            )
        n.repl <-
            ifelse(!is.null(n.repl),
                n.repl,
                max(1, ifelse(is.null(repl),
                    1,
                    ifelse(length(repl) == 0,
                        1,
                        max(repl)
                    )
                ))
            )


        ## Handle loc and index input semantics:
        if (is.null(loc)) {
            if (is.null(A.loc)) {
                A.loc <- Diagonal(n.spde, 1)
            }
        } else {
            if (is.null(mesh)) {
                  stop("'loc' specified but 'mesh' is NULL.")
              }

            A.loc <- fmesher::fm_evaluator(mesh, loc = loc)$proj$A
        }
        if (is.null(index)) {
            index <- seq_len(nrow(A.loc))
        }
        ## Now 'index' points into the rows of 'A.loc'

        if (is.null(n.block)) {
            n.block <- ifelse(is.null(block), length(index), max(block))
        }
        block.rescale <- match.arg(block.rescale)

        ## Handle group semantics:
        ## TODO: FIXME!!! group, group.index, group.mesh, A.group, etc
        if (!is.null(A.group)) {
            if (!is.null(group) || !is.null(group.mesh)) {
                warning("'A.group' has been specified; ignoring non-NULL 'group' or 'group.mesh'.")
            }
        } else if (!is.null(group.mesh)) {
            if (is.null(group)) {
                group <- rep(group.mesh$mid[1], length(index))
            }
        } else if (is.null(group)) {
            group <- rep(1L, length(index))
        } else if (length(group) == 1) {
            group <- rep(group, length(index))
        }
        if (is.null(group.index)) {
            group.index <- seq_len(length(group))
        }
        ## Now 'group.index' points into the rows of 'A.group' or 'group'
        if (length(group.index) != length(index)) {
            stop(paste("length(group.index) != length(index): ",
                length(group.index), " != ", length(index),
                sep = ""
            ))
        }

        if (!is.null(group.mesh) && is.null(A.group)) {
            A.group <- inla.mesh.1d.A(group.mesh, loc = group)
        }
        ## Now 'group.index' points into the rows of 'A.group' or 'group'

        ## Handle repl semantics:
        if (is.null(repl)) {
            repl <- rep(1, length(index))
        } else if (length(repl) == 1) {
            repl <- rep(repl, length(index))
        } else if (length(repl) != length(index)) {
            stop(paste("length(repl) != length(index): ",
                length(repl), " != ", length(index),
                sep = ""
            ))
        }

        if (length(index) > 0L) {
            A.loc <- inla.as.dgTMatrix(A.loc[index, , drop = FALSE])

            if (length(A.loc@i) > 0L) {
                if (is.null(weights)) {
                    weights <- rep(1, length(index))
                } else if (length(weights) == 1L) {
                    weights <- rep(weights[1], length(index))
                }
                if (!is.null(block)) {
                    ## Leave the rescaling until the block phase,
                    ## so that the proper rescaling can be determined.
                    block.weights <- weights
                    weights <- rep(1, length(index))
                }

                if (!is.null(A.group)) {
                    A.group <- inla.as.dgTMatrix(A.group[group.index, , drop = FALSE])
                    A <- fmesher::fm_row_kron(
                        A.group,
                        A.loc,
                        repl = repl,
                        n.repl = n.repl,
                        weights = weights
                    )
                    ## More general version:
                    ## A = fm_row_kron(A.repl,
                    ##   fm_row_kron(A.group, A.loc),
                    ##   weights=weights))
                } else {
                    i <- 1L + A.loc@i
                    group.i <- group[group.index[i]]
                    repl.i <- repl[i]
                    weights.i <- weights[i]
                    A <- (sparseMatrix(
                        i = i,
                        j = (1L + A.loc@j +
                            n.spde * (group.i - 1L) +
                            n.spde * n.group * (repl.i - 1L)),
                        x = weights.i * A.loc@x,
                        dims = (c(
                            length(index),
                            n.spde * n.group * n.repl
                        ))
                    ))
                }
                if (!is.null(block)) {
                    A <- (inla.spde.make.block.A(
                        A = A,
                        block = block,
                        n.block = n.block,
                        weights = block.weights,
                        rescale = block.rescale
                    ))
                }
            } else {
                A <- (sparseMatrix(
                    i = integer(0),
                    j = integer(0),
                    x = numeric(0),
                    dims = c(n.block, n.spde * n.group * n.repl)
                ))
            }
        } else {
            A <- (sparseMatrix(
                i = integer(0),
                j = integer(0),
                x = numeric(0),
                dims = c(0L, n.spde * n.group * n.repl)
            ))
        }
        return(A)
    }

#' @noRd
rbind.inla.data.stack.info <- function(...) {
  l <- list(...)
  names(l) <- NULL
  names.tmp <- do.call(c, lapply(l, function(x) x$names))
  ncol.tmp <- do.call(c, lapply(l, function(x) x$ncol))
  
  ncol <- c()
  names <- list()
  for (k in 1:length(names.tmp)) {
    name <- names(names.tmp)[k]
    if (!is.null(names[[name]])) {
      if (!identical(
        names[[name]],
        names.tmp[[k]]
      )) {
        stop("Name mismatch.")
      }
    }
    names[[name]] <- names.tmp[[k]]
    
    if (!is.null(as.list(ncol)[[name]])) {
      if (ncol[name] != ncol.tmp[[k]]) {
        stop("ncol mismatch.")
      }
    }
    ncol[name] <- ncol.tmp[[k]]
  }
  
  data <- dplyr::bind_rows(lapply(l, function(x) x[["data"]]))
  
  offset <- 0
  index <- list()
  for (k in 1:length(l)) {
    for (j in 1:length(l[[k]]$index)) {
      if (is.null(index[[names(l[[k]]$index)[j]]])) {
        index[[names(l[[k]]$index)[j]]] <- l[[k]]$index[[j]] + offset
      } else {
        index[[names(l[[k]]$index)[j]]] <-
          c(
            index[[names(l[[k]]$index)[j]]],
            l[[k]]$index[[j]] + offset
          )
      }
    }
    offset <- offset + l[[k]]$nrow
  }
  
  info <-
    list(
      data = data,
      nrow = nrow(data),
      ncol = ncol,
      names = names,
      index = index
    )
  class(info) <- "inla.data.stack.info"
  
  return(info)
}

#' @noRd
rbind.inla.stack.responses <- function(l) {
  null.l <- vapply(l, is.null, logical(1))
  l <- l[!null.l]
  if (length(l) == 0) {
    return(list(NULL))
  }
  
  l <- lapply(l, function(x) {
    if (inherits(x, "inla.mdata")) {
      # Make sure mdata is a data.frame
      attribs <- attributes(x)
      x <- as.data.frame(x)
      attr(x, "inla.ncols") <- attribs$inla.ncols
      attr(x, "names.ori") <- attribs$names.ori
      class(x) <- c("inla.mdata", "data.frame")
      x
    } else {
      x
    }
  })
  
  classes <- lapply(l, class)
  if (length(unique(classes)) > 1) {
    stop("Cannot rbind responses with different classes.")
  }
  
  attribs <- lapply(l, attributes)
  
  if (all(vapply(l, is.data.frame, logical(1)))) {
    response <- dplyr::bind_rows(l)
    class(response) <- classes[[1]]
    if (inherits(response, "inla.mdata")) {
      attr(response, "inla.ncols") <- attribs[[1]]$inla.ncols
      attr(response, "names.ori") <- attribs[[1]]$names.ori
    }
  } else {
    response <- do.call(c, l)
  }
  
  return(list(response))
}

#' @noRd
expand.inla.stack.responses <- function(l) {
  nrows <- vapply(l, NROW, integer(1))
  
  responses <- lapply(
    seq_along(l),
    function(k) {
      x <- l[[k]]
      if (is.data.frame(x)) {
        y <- dplyr::bind_rows(
          as.data.frame(
            matrix(
              NA,
              nrow = sum(nrows[seq_len(k - 1)]),
              ncol = 0
            )
          ),
          x,
          as.data.frame(
            matrix(
              NA,
              nrow = sum(nrows[-seq_len(k)]),
              ncol = 0
            )
          )
        )
        if (inherits(x, "inla.mdata")) {
          attr(y, "inla.ncols") <- attr(x, "inla.ncols")
          attr(y, "names.ori") <- attr(x, "names.ori")
          class(y) <- c("inla.mdata", "data.frame")
        }
      } else if (is.matrix(x)) {
        NA_ <- x[1, 1]
        is.na(NA_) <- TRUE
        y <- rbind(
          matrix(
            NA_,
            nrow = sum(nrows[seq_len(k - 1)]),
            ncol = ncol(x)
          ),
          x,
          matrix(
            NA_,
            nrow = sum(nrows[-seq_len(k)]),
            ncol = ncol(x)
          )
        )
      } else if (is.vector(x)) {
        NA_ <- x[1]
        is.na(NA_) <- TRUE
        y <- c(
          rep(NA_, sum(nrows[seq_len(k - 1)])),
          x,
          rep(NA_, sum(nrows[-seq_len(k)])))
      } else {
        stop("Don't know how to expand responses of class '",
             paste0(class(x), collapse = ", "), "'.")
      }
      y
    })
  

  return(responses)
}

#' @describeIn inla.stack Remove unused entries from an existing stack
#'
#' @export
inla.stack.remove.unused <- function(stack) {
    inla.require.inherits(stack, "inla.data.stack", "'stack'")

    if (stack$effects$nrow < 2) {
        return(stack)
    }

    ## Remove components with no effect:
    remove <- rep(FALSE, stack$effects$nrow)
    remove.unused.indices <-
        which(colSums(abs(stack$A[, , drop = FALSE])) == 0)
    remove[remove.unused.indices] <- TRUE

    index.new <- rep(as.integer(NA), stack$effects$nrow)

    ncol.A <- sum(!remove)
    if (ncol.A > 0) {
          index.new[!remove] <- seq_len(ncol.A)
      }
    index.new[remove] <- index.new[index.new[remove]]

    for (k in 1:length(stack$effects$index)) {
        stack$effects$index[[k]] <- index.new[stack$effects$index[[k]]]
    }

    A <- inla.as.dgTMatrix(stack$A)
    j.new <- index.new[A@j + 1L]
    ## Check for any zero-elements in remove.unused-columns:
    ok <- !is.na(j.new)
    stack$A <-
        sparseMatrix(
            i = A@i[ok] + 1L,
            j = j.new[ok],
            x = A@x[ok],
            dims = c(nrow(A), ncol.A)
        )

    stack$effects$data <- stack$effects$data[!remove, , drop = FALSE]
    stack$effects$nrow <- ncol.A

    return(stack)
}

#' @describeIn inla.stack Compress an existing stack by removing duplicates
#'
#' @export
inla.stack.compress <- function(stack, remove.unused = TRUE) {
    inla.require.inherits(stack, "inla.data.stack", "'stack'")

    if (stack$effects$nrow < 2) {
        return(stack)
    }

    ii <- do.call(order, as.list(stack$effects$data))
    jj.dupl <-
        which(1L ==
            diff(c(
                duplicated(stack$effects$data[ii, , drop = FALSE]),
                FALSE
            )))
    kk.dupl <-
        which(-1L ==
            diff(c(
                duplicated(stack$effects$data[ii, , drop = FALSE]),
                FALSE
            )))
    ## ii[jj.dupl] are the rows that have duplicates.
    ## ii[(jj.dupl[k]+1):kk.dupl[k]] are the duplicate rows for each k

    remove <- rep(FALSE, stack$effects$nrow)
    index.new <- rep(as.integer(NA), stack$effects$nrow)

    if (length(jj.dupl) > 0) {
        for (k in 1:length(jj.dupl)) {
            i <- ii[jj.dupl[k]]
            j <- ii[(jj.dupl[k] + 1):kk.dupl[k]]

            remove[j] <- TRUE
            index.new[j] <- i
        }
    }

    ncol.A <- sum(!remove)
    if (ncol.A > 0) {
          index.new[!remove] <- seq_len(ncol.A)
      }
    index.new[remove] <- index.new[index.new[remove]]

    for (k in 1:length(stack$effects$index)) {
        stack$effects$index[[k]] <- index.new[stack$effects$index[[k]]]
    }

    A <- inla.as.dgTMatrix(stack$A)
    j.new <- index.new[A@j + 1L]
    ## Check for any zero-elements in remove.unused-columns:
    ok <- !is.na(j.new)
    stack$A <-
        sparseMatrix(
            i = A@i[ok] + 1L,
            j = j.new[ok],
            x = A@x[ok],
            dims = c(nrow(A), ncol.A)
        )

    stack$effects$data <- stack$effects$data[!remove, , drop = FALSE]
    stack$effects$nrow <- ncol.A

    if (remove.unused) {
        return(inla.stack.remove.unused(stack))
    } else {
        return(stack)
    }
}





#' Data stacking for advanced INLA models
#'
#' Functions for combining data, effects and observation matrices into
#' `inla.stack` objects, and extracting information from such objects.
#'
#' For models with a single effects collection, the outer list container for
#' `A` and `effects` may be omitted.
#'
#' Component size definitions:
#' * \eqn{n_l}{n_l} effect blocks
#' * \eqn{n_k}{n_k} effects
#' * \eqn{n_i}{n_i} data values
#' * \eqn{n_{j,l}}{n_jl} effect size for block \eqn{l}{l}
#' * \eqn{n_j}{n_j} \eqn{= \sum_{l=1}^{n_l} n_{j,l}}{sum_l n_jl} total
#' effect size
#'
#' Input: \describe{
#' \item{`data`}{\eqn{(y^1, \ldots, y^p)}{(y1,\dots,y2)} \eqn{p}{p}
#' vectors, each of length \eqn{n_i}{n_i}}
#' \item{`A`}{\eqn{(A^1, \ldots, A^{n_l})}{(A1,\dots,A2)} matrices of size
#' \eqn{n_i \times n_{j,l}}{n_i by n_jl}}
#' \item{`effects`}{\eqn{\left((x^{1,1},\ldots,x^{n_k,1}), \ldots,
#' (x^{1,n_l},\ldots,x^{n_k,n_l})\right)}{((x_[1,1],\dots,x_[n_k,1]),\dots(x_[1,n_l],\dots,x_[n_k,n_l]))}
#' collections of effect
#' vectors of length \eqn{n_{j,l}}{n_jl} }
#' }
#'
#' \deqn{\mbox{predictor}(y^1, \ldots, y^p) \sim
#' \sum_{l=1}^{n_l} A^l \sum_{k=1}^{n_k} g(k, x^{k,l})
#' = \tilde{A} \sum_{k=1}^{n_k} g(k, \tilde{x}^k)
#' }{ predictor(y^1, \ldots, y^p) ~ sum_{l=1}^{n_l} A^l sum_{k=1}^{n_k}
#' g(k, x^{k,l}) = tilde{A} sum_{k=1}^{n_k} g(k, tilde{x}^k) }
#' where
#' \deqn{\tilde{A} = \mbox{cbind}\left( A^1, \ldots, A^{n_l} \right)
#' }{ tilde{A} = cbind( A^1, ..., A^{n_l} ) }
#' and
#' \deqn{\tilde{x}^k = \mbox{rbind}\left( x^{k,1}, \ldots, x^{k,n_l} \right)
#' }{ tilde{x}^k = rbind( x^{k,1}, ..., x^{k,n_l} ) }
#' and for each block \eqn{l}{l}, any missing
#' \eqn{x^{k,l}} is replaced by an `NA` vector.
#'
#' @aliases inla.stack inla.stack.remove.unused inla.stack.compress
#' inla.stack.sum inla.stack.join inla.stack.index inla.stack.LHS
#' inla.stack.RHS inla.stack.data inla.stack.A
#' @param stack A `inla.data.stack` object, created by a call to
#' `inla.stack`, `inla.stack.sum`, or `inla.stack.join`.
#' @param remove.unused If `TRUE`, compress the model by removing rows of
#' effects corresponding to all-zero columns in the `A` matrix (and
#' removing those columns).
#' @param multi.family logical or character; in `inla.data.join`, if `TRUE`,
#' the `response` part of the stack is joined as a `list`. If `character`,
#' denotes the name of a `data` element that should be joined as a multi-column
#' matrix. Default is `FALSE`, which joins both the `data` and `response`
#' with regular row binding with `dplyr::bind_rows`.
#' @param ... For `inla.stack.join`, two or more data stacks of class
#' `inla.data.stack`, created by a call to `inla.stack`,
#' `inla.stack.sum`, or `inla.stack.join`. For
#' `inla.stack.data`, a list of variables to be joined with the data list.
#' @param compress If `TRUE`, compress the model by removing duplicated
#' rows of effects, replacing the corresponding A-matrix columns with a single
#' column containing the sum.
#' @param data A list or codedata.frame of named data vectors. Scalars are
#' expanded to match the number of rows in the A matrices, or any non-scalar
#' data vectors. An error is given if the input is inconsistent.
#' @param A A list of observation matrices. Scalars are expanded to diagonal
#' matrices matching the effect vector lengths. An error is given if the input
#' is inconsistent or ambiguous.
#' @param effects A collection of effects/predictors.  Each list element
#' corresponds to an observation matrix, and must either be a single vector, a
#' list of vectors, or a `data.frame`. Single-element effect vectors are
#' expanded to vectors matching the number of columns in the corresponding A
#' matrix.  An error is given if the input is inconsistent or ombiguous.
#' @param tag A string specifying a tag for later identification.
#' @return A data stack of class `inla.data.stack`.
#' Elements: \itemize{
#' \item`data` \eqn{=(y^1, \ldots, y^p, \tilde{x}^1, \ldots,
#' \tilde{x}^{n_k})}{=(y^1, \dots, y^p, tilde{x}^1, \dots,
#' tilde{x}^{n_k})}
#' \item`A` \eqn{=\tilde{A}}{=tilde{A}}
#' \item`data.names` List
#' of data names, length \eqn{p}
#' \item`effect.names` List of effect names,
#' length \eqn{n_k}
#' \item`n.data` Data length, \eqn{n_i}
#' \item`index`
#' List indexed by `tag`s, each element indexing into \eqn{i=1, \ldots,
#' n_i} }
#' @section Functions: \itemize{
#' \item `inla.stack.remove.unused`: Remove
#' unused entries from an existing stack
#'
#' \item `inla.stack.compress`: Compress an existing stack by removing
#' duplicates
#'
#' \item `inla.stack`: Shorthand for inla.stack.join and inla.stack.sum
#'
#' \item `inla.stack.sum`: Create data stack as a sum of predictors
#'
#' \item `inla.stack.join`: Join two or more data stacks
#'
#' \item `inla.stack.index`: Extract tagged indices
#'
#' \item `inla.stack.LHS`: Extract data associated with the "left hand
#' side" of the model (e.g. the data itself, `Ntrials`, `link`,
#' `E`)
#'
#' \item `inla.stack.RHS`: Extract data associated with the "right hand
#' side" of the model (all the covariates/predictors)
#'
#' \item `inla.stack.data`: Extract data for an inla call, and optionally
#' join with other variables
#'
#' \item `inla.stack.A`: Extract the "A matrix" for control.predictor }
#' @seealso [inla.spde.make.A()], [inla.spde.make.index()]
#' @keywords fmesher
#' @examples
#'
#' n <- 200
#' loc <- matrix(runif(n * 2), n, 2)
#' mesh <- inla.mesh.2d(
#'     loc.domain = loc,
#'     max.edge = c(0.05, 0.2)
#' )
#' proj.obs <- inla.mesh.projector(mesh, loc = loc)
#' proj.pred <- inla.mesh.projector(mesh, loc = mesh$loc)
#' spde <- inla.spde2.pcmatern(mesh,
#'     prior.range = c(0.01, 0.01),
#'     prior.sigma = c(10, 0.01)
#' )
#'
#' covar <- rnorm(n)
#' field <- inla.qsample(n = 1, Q = inla.spde.precision(spde, theta = log(c(0.5, 1))))[, 1]
#' y <- 2 * covar + inla.mesh.project(proj.obs, field)
#'
#' A.obs <- inla.spde.make.A(mesh, loc = loc)
#' A.pred <- inla.spde.make.A(mesh, loc = proj.pred$loc)
#' stack.obs <-
#'     inla.stack(
#'         data = list(y = y),
#'         A = list(A.obs, 1),
#'         effects = list(c(
#'             list(Intercept = 1),
#'             inla.spde.make.index("spatial", spde$n.spde)
#'         ),
#'         covar = covar
#'         ),
#'         tag = "obs"
#'     )
#' stack.pred <-
#'     inla.stack(
#'         data = list(y = NA),
#'         A = list(A.pred),
#'         effects = list(c(
#'             list(Intercept = 1),
#'             inla.spde.make.index("spatial", mesh$n)
#'         )),
#'         tag = "pred"
#'     )
#' stack <- inla.stack(stack.obs, stack.pred)
#'
#' formula <- y ~ -1 + Intercept + covar + f(spatial, model = spde)
#' result1 <- inla(formula,
#'     data = inla.stack.data(stack.obs, spde = spde),
#'     family = "gaussian",
#'     control.predictor = list(
#'         A = inla.stack.A(stack.obs),
#'         compute = TRUE
#'     )
#' )
#'
#' plot(y, result1$summary.fitted.values[inla.stack.index(stack.obs, "obs")$data, "mean"],
#'     main = "Observations vs posterior predicted values at the data locations"
#' )
#'
#' result2 <- inla(formula,
#'     data = inla.stack.data(stack, spde = spde),
#'     family = "gaussian",
#'     control.predictor = list(
#'         A = inla.stack.A(stack),
#'         compute = TRUE
#'     )
#' )
#'
#' field.pred <- inla.mesh.project(
#'     proj.pred,
#'     result2$summary.fitted.values[inla.stack.index(stack, "pred")$data, "mean"]
#' )
#' field.pred.sd <- inla.mesh.project(
#'     proj.pred,
#'     result2$summary.fitted.values[inla.stack.index(stack, "pred")$data, "sd"]
#' )
#'
#' plot(field, field.pred, main = "True vs predicted field")
#' abline(0, 1)
#' image(inla.mesh.project(mesh,
#'     field = field,
#'     dims = c(200, 200)
#' ),
#' main = "True field"
#' )
#' image(inla.mesh.project(mesh,
#'     field = field.pred,
#'     dims = c(200, 200)
#' ),
#' main = "Posterior field mean"
#' )
#' image(inla.mesh.project(mesh,
#'     field = field.pred.sd,
#'     dims = c(200, 200)
#' ),
#' main = "Prediction standard deviation"
#' )
#' plot(field, (field.pred - field) / 1,
#'     main = "True field vs standardised prediction residuals"
#' )
#' @export inla.stack
inla.stack <- function(..., compress = TRUE, remove.unused = TRUE, multi.family = FALSE) {
    if (all(sapply(list(...), function(x) inherits(x, "inla.data.stack")))) {
        return(do.call(
            inla.stack.join,
            c(list(...),
                compress = compress,
                remove.unused = remove.unused,
                multi.family = multi.family
            )
        ))
    } else {
        return(do.call(
            inla.stack.sum,
            c(list(...),
                compress = compress,
                remove.unused = remove.unused
            )
        ))
    }
}


#' @describeIn inla.stack Create data stack as a sum of predictors
#' @param response A list of response vectors, matrices, data.frame, or other special
#' response objects, such as `inla.mdata`. In ordinary user-side code,
#' the list has length 1.
#'
#' @export
inla.stack.sum <- function(data, A, effects, response = NULL,
                           tag = "",
                           compress = TRUE,
                           remove.unused = TRUE) {
    input.nrow <- function(x) {
        return(inla.ifelse(
          is.matrix(x) || is(x, "Matrix"),
          nrow(x),
          inla.ifelse(
            inherits(x, "inla.mdata"),
            stop("inla.mdata objects must be given as 'response' input"),
            inla.ifelse(
              is.data.frame(x),
              rep(nrow(x), ncol(x)),
              length(x)
            )
          )
        ))
    }
    input.ncol <- function(x) {
      return(inla.ifelse(
        is.matrix(x) || is(x, "Matrix"),
        ncol(x),
        inla.ifelse(
          inherits(x, "inla.mdata"),
          rep(1L, length(x)),
          inla.ifelse(
            is.data.frame(x),
            rep(1L, ncol(x)),
            1L
          )
        )
      ))
    }

    input.list.nrow <- function(l) {
        if (is.data.frame(l)) {
              return(input.nrow(l))
          }
        return(do.call(c, lapply(l, input.nrow)))
    }
    input.list.ncol <- function(l) {
        if (is.data.frame(l)) {
              return(input.ncol(l))
          }
        return(do.call(c, lapply(l, input.ncol)))
    }
    input.list.names <- function(l) {
        if (is.data.frame(l)) {
              return(colnames(l))
          }
        is.df <- vapply(l, is.data.frame, logical(1))
        name <- vector("list", length(l))
        if (!is.null(names(l))) {
              name[!is.df] <-
                  lapply(
                      names(l)[!is.df],
                      function(x) list(x)
                  )
          } else {
              name[!is.df] <- ""
          }
        name[is.df] <-
            lapply(
                l[is.df],
                function(x) as.list(colnames(x))
            )

        return(do.call(c, name))
    }


    parse.input.list <- function(l, n.A, error.tag, tag = "", n.A.strict = FALSE) {
        ncol <- input.list.ncol(l)
        nrow <- input.list.nrow(l)
        names <- input.list.names(l)
        if ((n.A > 1) && any(nrow == 1)) {
            for (k in which(nrow == 1)) {
                if (ncol[k] == 1) {
                    l[[k]] <- rep(l[[k]], n.A)
                    nrow[k] <- n.A
                } else {
                    stop(paste(error.tag,
                        "Automatic expansion only available for scalars.",
                        sep = ""
                    ))
                }
            }
        }

        if (length(unique(c(names, ""))) < length(c(names, ""))) {
            stop(paste(error.tag,
                "All variables must have unique names\n",
                "Names: ('",
                paste(names, collapse = "', '", sep = ""),
                "')",
                sep = ""
            ))
        }

        for (k in seq_len(length(names))) {
            if (ncol[k] == 1) {
                names(names)[k] <- names[[k]][[1]]
                names[[k]] <- c(names[[k]][[1]])
            } else {
                names(names)[k] <- names[[k]][[1]]
                names[[k]] <- paste(names[[k]][[1]], ".", seq_len(ncol[k]), sep = "")
            }
        }

        names(nrow) <- names(names)
        names(ncol) <- names(names)

        ## data = as.data.frame(do.call(cbind, l))
        data <- as.data.frame(l)
        if (!is.null(names)) {
          names(data) <- do.call(c, names)
        }
        nrow <- nrow(data)
        if ((n.A > 1 || n.A.strict) && (nrow != n.A)) {
            stop(paste(error.tag,
                "Mismatching row sizes: ",
                paste(nrow, collapse = ",", sep = ""),
                ", n.A=", n.A,
                sep = ""
            ))
        }

        index <- list(seq_len(nrow))
        if (!is.null(tag)) {
            names(index) <- tag
        }

        info <- list(data = data, nrow = nrow, ncol = ncol, names = names, index = index)
        class(info) <- "inla.data.stack.info"

        return(info)
    }

    if (is.null(tag)) {
          stop("'tag' must not be 'NULL'")
      }

    ## Check if only a single block was specified.
    if (!is.list(A)) {
        A <- list(A)
        effects <- list(effects)
    }
    if (length(A) != length(effects)) {
          stop(paste("length(A)=", length(A),
              " should be equal to length(effects)=", length(effects),
              sep = ""
          ))
      }

    n.effects <- length(effects)

    eff <- list()
    for (k in 1:n.effects) {
        if (is.data.frame(effects[[k]])) {
            eff[[k]] <-
                parse.input.list(
                    list(effects[[k]]),
                    input.ncol(A[[k]]),
                    paste("Effect block ", k, ":\n", sep = ""),
                    tag
                )
        } else {
            if (!is.list(effects[[k]])) {
                tmp <-
                    inla.ifelse(
                        is.null(names(effects)[k]),
                        "",
                        names(effects)[k]
                    )
                effects[[k]] <- list(effects[[k]])
                names(effects[[k]]) <- tmp
            }
            eff[[k]] <-
                parse.input.list(
                    effects[[k]],
                    input.ncol(A[[k]]),
                    paste("Effect block ", k, ":\n", sep = ""),
                    tag
                )
        }
    }

    for (k in 1:n.effects) {
        if (is.vector(A[[k]])) {
            A[[k]] <- Matrix(A[[k]], input.nrow(A[[k]]), 1)
        }
        if ((input.ncol(A[[k]]) == 1) && (eff[[k]]$nrow > 1)) {
            if (input.nrow(A[[k]]) != 1) {
                  stop(paste("ncol(A) does not match nrow(effect) for block ",
                      k, ": ",
                      input.ncol(A[[k]]), " != ", eff[[k]]$nrow,
                      sep = ""
                  ))
              }
            A[[k]] <- Diagonal(eff[[k]]$nrow, A[[k]][1, 1])
        } else if (input.ncol(A[[k]]) != eff[[k]]$nrow) {
            stop(paste("ncol(A) does not match nrow(effect) for block ",
                k, ": ",
                input.ncol(A[[k]]), " != ", eff[[k]]$nrow,
                sep = ""
            ))
        }
    }
    if (length(unique(input.list.nrow(A))) > 1) {
        stop(paste("Row count mismatch for A: ",
            paste(input.list.nrow(A), collapse = ",", sep = ""),
            sep = ""
        ))
    }
    A.nrow <- nrow(A[[1]])
    A.ncol <- input.list.ncol(A)

    data <-
        parse.input.list(inla.ifelse(
            is.data.frame(data),
            as.list(data),
            data
        ),
        A.nrow,
        paste("Data block:\n", sep = ""),
        tag,
        n.A.strict = TRUE
        )

    effects <- do.call(rbind.inla.data.stack.info, eff)

    A.matrix <- do.call(cbind, A)
    A.nrow <- nrow(A.matrix)
    A.ncol <- ncol(A.matrix)

    if (length(unique(c(names(data$names), names(effects$names)))) <
        length(c(names(data$names), names(effects$names)))) {
        stop(paste("Names for data and effects must not coincide.\n",
            "Data names:   ",
            paste(names(data$names), collapse = ", ", sep = ""),
            "\n",
            "Effect names: ",
            paste(names(effects$names), collapse = ", ", sep = ""),
            sep = ""
        ))
    }

    stack <- list(
      A = A.matrix,
      data = data,
      effects = effects,
      response = response
    )
    class(stack) <- "inla.data.stack"

    if (compress) {
        return(inla.stack.compress(stack, remove.unused = remove.unused))
    } else if (remove.unused) {
        return(inla.stack.remove.unused(stack))
    } else {
        return(stack)
    }
}

# Expand observation vectors/matrices in stacks into to a multicolumn matrix for multiple likelihoods
#
# @aliases inla.stack.mexpand
# @name inla.stack.mexpand
# @export
# @param ... List of stacks that contain vector observations
#            (existing multilikelihood observation matrices are also permitted)
# @param old.names A vector of strings with the names of the observation vector/matrix for each stack.
#        If a single string, this is assumed for all the stacks. (default "BRU.response")
# @param new.name The name to be used for the expanded observation matrix,
#        possibly the same as an old name. (default "BRU.response")
# @return a list of modified stacks with multicolumn observations
# @author Fabian E. Bachl \email{f.e.bachl@@bath.ac.uk} and Finn Lindgren \email{finn.lindgren@@gmail.com}
#

inla.stack.mexpand <- function(...,
                               old.names = "response",
                               new.name = "response") {
  stacks <- list(...)
  if (length(old.names) == 1) {
    old.names <- rep(old.names, length(stacks))
  }
  y.cols <- unlist(lapply(seq_along(stacks),
                          function(x, stacks, old.names) {
                            LHS <- INLA::inla.stack.LHS(stacks[[x]])[[old.names[x]]]
                            ifelse(is.vector(LHS), 1, ncol(LHS))
                          },
                          stacks = stacks, old.names = old.names
  ))
  y.offset <- c(0, cumsum(y.cols))
  y.cols.total <- sum(y.cols)
  for (j in seq_along(stacks)) {
    LHS <- INLA::inla.stack.LHS(stacks[[j]])
    RHS <- INLA::inla.stack.RHS(stacks[[j]])
    A <- INLA::inla.stack.A(stacks[[j]])
    response <- stacks[[j]][["response"]]
    # Access the raw tag indexing information
    tags <- list(
      data = stacks[[j]]$data$index,
      effects = stacks[[j]]$effects$index
    )
    
    # Expand the observation vector/matrix into a multilikelihood observation matrix:
    y.rows <- ifelse(is.vector(LHS[[old.names[j]]]),
                     length(LHS[[old.names[j]]]),
                     nrow(LHS[[old.names[j]]])
    )
    LHS[[new.name]] <-
      cbind(
        matrix(NA, nrow = y.rows, ncol = y.offset[j]),
        LHS[[old.names[j]]],
        matrix(NA, nrow = y.rows, ncol = y.cols.total - y.offset[j + 1])
      )
    
    # Create the modified stack, with model compression disabled to prevent modifications:
    stacks[[j]] <-
      INLA::inla.stack.sum(data = LHS, A = A, effects = RHS, compress = FALSE, remove.unused = FALSE,
                           response = response)
    # Since the row indexing is unchanged, copy the tag index information:
    stacks[[j]]$data$index <- tags$data
    stacks[[j]]$effects$index <- tags$effects
  }
  stacks
}

#' @describeIn inla.stack Join two or more data stacks
#'
#' @export
inla.stack.join <- function(..., compress = TRUE, remove.unused = TRUE, multi.family = FALSE) {
    if (is.character(multi.family)) {
        # NA expand response variable in the data part as matrix
        S.input <- inla.stack.mexpand(...,
                                      old.names = multi.family,
                                      new.name = multi.family)
        # Nothing more to do for multi-likelihoods
        # For non-multi, need to join the matrices (see below)
        multi.family <- FALSE
    } else {
      S.input <- list(...)
    }
  
    # Join responses as a list
    responses <- do.call(c,
                        lapply(S.input,
                               function(x) {
                                 if (is.null(x[["response"]])) {
                                   list(NULL)
                                 } else {
                                   x[["response"]]
                                 }
                               }))
    if (isFALSE(multi.family)) {
        responses <- rbind.inla.stack.responses(responses)
    }
    if (isTRUE(multi.family)) {
      responses <- expand.inla.stack.responses(responses)
    }

    data <- do.call(
        rbind.inla.data.stack.info,
        lapply(S.input, function(x) x$data)
    )
    effects <- do.call(
        rbind.inla.data.stack.info,
        lapply(S.input, function(x) x$effects)
    )
    ## The .bdiag form of bdiag takes a list as input.
    A <- .bdiag(lapply(S.input, function(x) x$A))

    S.output <- list(A = A, data = data, effects = effects, responses = responses)
    class(S.output) <- "inla.data.stack"

    if (length(unique(c(names(data$names), names(effects$names)))) <
        length(c(names(data$names), names(effects$names)))) {
        stop(paste("Names for data and effects must not coincide.\n",
            "Data names:   ",
            paste(names(data$names), collapse = ", ", sep = ""),
            "\n",
            "Effect names: ",
            paste(names(effects$names), collapse = ", ", sep = ""),
            sep = ""
        ))
    }

    if (compress) {
        return(inla.stack.compress(S.output, remove.unused = remove.unused))
    } else if (remove.unused) {
        return(inla.stack.remove.unused(S.output))
    } else {
        return(S.output)
    }
}





#' @describeIn inla.stack Extract tagged indices
#'
#' @export
inla.stack.index <- function(stack, tag) {
    inla.require.inherits(stack, "inla.data.stack", "'stack'")

    if (is.null(tag)) {
        return(list(
            data = as.vector(do.call(c, stack$data$index)),
            effects = (stack$data$nrow +
                as.vector(do.call(c, stack$effects$index)))
        ))
    } else {
        return(list(
            data = as.vector(do.call(c, stack$data$index[tag])),
            effects = (stack$data$nrow +
                as.vector(do.call(c, stack$effects$index[tag])))
        ))
    }
}

## Extract an object list from an inla.stack internal structure
inla.stack.do.extract <- function(dat) {
    inla.require.inherits(dat, "inla.data.stack.info", "'dat'")

    handle.entry <- function(x) {
        if (dat$ncol[[x]] > 1) {
            return(matrix(
                do.call(
                    c,
                    dat$data[dat$names[[x]]]
                ),
                dat$nrow,
                dat$ncol[[x]]
            ))
        } else if (is.factor(dat$data[[dat$names[[x]]]])) {
            return(dat$data[[dat$names[[x]]]])
        }
        return(as.vector(dat$data[[dat$names[[x]]]]))
    }

    out <- lapply(names(dat$names), handle.entry)
    names(out) <- names(dat$names)

    return(out)
}


#' @describeIn inla.stack Extract data associated with the "left hand side" of the model
#' (e.g. the data itself, `Ntrials`, `link`, `E`)
#'
#' @export
inla.stack.LHS <- function(stack) {
    inla.require.inherits(stack, "inla.data.stack", "'stack'")

    return(inla.stack.do.extract(stack$data))
}

#' @describeIn inla.stack Extract data associated with the "right hand side" of the model
#' (all the covariates/predictors)
#'
#' @export
inla.stack.RHS <- function(stack) {
    inla.require.inherits(stack, "inla.data.stack", "'stack'")

    return(inla.stack.do.extract(stack$effects))
}

#' @describeIn inla.stack Extract data for an inla call, and optionally join with other variables
#'
#' @export
inla.stack.data <- function(stack, ...) {
    inla.require.inherits(stack, "inla.data.stack", "'stack'")

    return(c(
        inla.stack.do.extract(stack$data),
        inla.stack.do.extract(stack$effects),
        list(...)
    ))
}

#' @describeIn inla.stack Extract the "A matrix" for control.predictor
#'
#' @export
inla.stack.A <- function(stack) {
  inla.require.inherits(stack, "inla.data.stack", "'stack'")
  return(stack$A)
}

#' @describeIn inla.stack Extract the response variable or list of response objects
#'
#' @export
inla.stack.response <- function(stack) {
  inla.require.inherits(stack, "inla.data.stack", "'stack'")
  return(stack$response)
}

#' @describeIn inla.stack Print information about an `inla.data.stack`
#'
#' @method print inla.data.stack
#' @export
print.inla.data.stack <- function(x, ...) {
  inla.require.inherits(x, "inla.data.stack", "'stack'")
  LHS <- inla.stack.LHS(x)
  RHS <- inla.stack.RHS(x)
  A <- inla.stack.A(x)
  response <- inla.stack.response(x)
  LHS_n <- if (is.data.frame(LHS)) {
    NROW(LHS)
  } else if (length(LHS) > 0) {
    unique(vapply(LHS, NROW, 1L))
  } else {
    0
  }
  RHS_n <- if (is.data.frame(RHS)) {
    NROW(RHS)
  } else if (length(RHS) > 0) {
    unique(vapply(RHS, NROW, 1L))
  } else {
    0
  }
  cat("Data stack with\n  ",
      "data:    ", paste0(names(LHS), collapse = ", "),
      ", size: ", paste0(LHS_n, collapse = ", "), "\n  ",
      "effects: ", paste0(names(RHS), collapse = ", "),
      ", size: ", paste0(RHS_n, collapse = ", "), "\n  ",
      "A:       ", nrow(A), " times ", ncol(A), "\n",
      sep = "")
  if (isTRUE(x[["multi_family"]])) {
    cat("  response: ", length(response), "response objects\n",
        sep = "")
  }
  return(invisible(x))
}



