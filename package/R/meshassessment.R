## meshassessment.R
##
## Copyright (C) 2016, 2018, Finn Lindgren
##
## Permission is hereby granted, free of charge, to any person obtaining a copy of this software
## and associated documentation files (the “Software”), to deal in the Software without
## restriction, including without limitation the rights to use, copy, modify, merge, publish,
## distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the
## Software is furnished to do so, subject to the following conditions:
## 
## The above copyright notice and this permission notice shall be included in all copies or
## substantial portions of the Software.
## 
## THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
## BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
## NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
## DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
## FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.



#' Interactive mesh building and diagnostics
#'
#' Assess the finite element approximation errors in a mesh for interactive R
#' sessions.  More detailed assessment tools are in [meshbuilder()].
#'
#' @param mesh An `inla.mesh`
#' @param spatial.range numeric; the spatial range parameter to use for the assessment
#' @param alpha numeric; A valid `inla.spde2.pcmatern` `alpha` parameter
#' @param dims 2-numeric; the grid size
#' @author Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @seealso [fmesher::fm_mesh_2d()], [fmesher::fm_rcdt_2d()], meshbuilder
#' @examples
#'
#' library(fmesher)
#' bnd <- fm_segm(cbind(
#'     c(0, 10, 10, 0, 0),
#'     c(0, 0, 10, 10, 0)
#' ), is.bnd = TRUE)
#' mesh <- fm_mesh_2d_inla(boundary = bnd, max.edge = 1)
#' out <- inla.mesh.assessment(mesh, spatial.range = 3, alpha = 2)
#' 
#' @rdname meshassessment
#' @export inla.mesh.assessment
inla.mesh.assessment <- function(mesh, spatial.range, alpha = 2,
                                 dims = c(500, 500)) {
    mesh.edgelengths <- function(mesh, proj) {
        i <- c()
        val <- c()
        num <- c()
        tri_idx <- rbind(c(1, 2), c(2, 3), c(3, 1))
        for (k in 1:3) {
            ti1 <- tri_idx[k, 1]
            ti2 <- tri_idx[k, 2]
            len <- Matrix::rowSums((mesh$loc[mesh$graph$tv[, ti2], ] -
                mesh$loc[mesh$graph$tv[, ti1], ])^2)^0.5
            i <- c(i, mesh$graph$tv[, ti1], mesh$graph$tv[, ti2])
            val <- c(val, rep(as.vector(len), times = 2))
            num <- c(num, rep(1, 2 * length(as.vector(len))))
        }
        avg_len <-
            as.vector(sparseMatrix(i = i, j = rep(1, length(i)), x = val)) /
                as.vector(sparseMatrix(i = i, j = rep(1, length(i)), x = num))

        b <- fmesher::fm_basis(proj, full = TRUE)
        proj_len <- as.vector(b$A %*% avg_len)
        proj_len[!b$ok] <- NA
        proj_len
    }
    mesh.proj <- function(mesh, dims) {
        fmesher::fm_evaluator(mesh, dims = dims)
    }
    mesh.spde <- function(mesh, alpha) {
        INLA::inla.spde2.pcmatern(mesh,
            alpha = alpha,
            prior.range = c(1, 0.5),
            prior.sigma = c(1, 0.5)
        )
    }
    mesh.Q <- function(spde, spatial.range) {
        INLA::inla.spde.precision(spde, theta = log(c(spatial.range, 1)))
    }
    mesh.S <- function(Q) {
        INLA::inla.qinv(Q, reordering = INLA::inla.reorderings())
    }
    mesh.sd <- function(proj, S) {
        b <- fmesher::fm_basis(proj, full = TRUE)
        v <- Matrix::rowSums(b$A * (b$A %*% S))
        v[!b$ok] <- NA
        matrix(v^0.5, length(proj$x), length(proj$y))
    }
    mesh.sd.deviation.approx <- function(proj, S, sd0) {
        b <- fmesher::fm_basis(proj, full = TRUE)
        val <- b$A %*% (
            as.vector(t(b$A[b$ok, , drop = FALSE]) %*%
                as.vector(sd0)[b$ok]) /
                colSums(b$A[b$ok, , drop = FALSE]))
        val[!b$ok] <- NA
        matrix(
            1 + (as.vector(sd0) - val),
            length(proj$x), length(proj$y)
        )
    }
    mesh.sd.bound <- function(proj, S) {
        fmesher::fm_evaluate(proj, field = diag(mesh.S())^0.5)
    }


    if (mesh$manifold != "R2") {
        warning("inla.mesh.assessment has only been tested on flat 2D manifolds.")
    }

    spde <- mesh.spde(mesh = mesh, alpha = alpha)
    Q <- mesh.Q(spde = spde, spatial.range = spatial.range)
    S <- mesh.S(Q = Q)
    proj <- mesh.proj(mesh = mesh, dims = dims)
    sd0 <- mesh.sd(proj = proj, S = S)
    sd.deviation <- mesh.sd.deviation.approx(proj = proj, S = S, sd0 = sd0)
    edgelengths <- mesh.edgelengths(mesh, proj)

    out <- data.frame(
        x = proj$lattice$loc[, 1],
        y = proj$lattice$loc[, 2],
        sd = as.vector(sd0),
        sd.dev = as.vector(sd.deviation),
        edge.len = edgelengths
    )
    sp::coordinates(out) <- ~ x + y
    out
}
