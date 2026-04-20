## Handle some imports specially to avoid warnings
## from clashes with Matrix
#' @rawNamespace import(stats, except = c(cov2cor, toeplitz, update))
#' @rawNamespace import(graphics, except = c(image))
#' @import grDevices
#' @rawNamespace import(utils, except= c(tail, head))
#' @import methods
#' @import splines
#' @import Matrix
NULL
