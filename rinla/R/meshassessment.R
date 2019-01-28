## meshassessment.R
##
##   Copyright (C) 2016, 2018, Finn Lindgren
##
##   This program is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Only one function is exported:
## Export: inla.mesh.assessment

##!\name{inla.mesh.assessment}
##!\alias{inla.mesh.assessment}
##!\title{Interactive mesh building and diagnostics}
##!\description{Assess the finite element
##!  approximation errors in a mesh for interactive R sessions.
##!  More detailed assessment tools are in \code{\link{meshbuilder}}.
##!}
##!\usage{
##!inla.mesh.assessment(mesh)
##!}
##!\author{Finn Lindgren \email{finn.lindgren@gmail.com}}
##!\seealso{inla.mesh.2d, inla.mesh.create, meshbuilder}
##!\examples{
##!bnd <- inla.mesh.segment(cbind(c(0, 10, 10, 0, 0),
##!                               c(0, 0, 10, 10, 0)), bnd = TRUE)
##!mesh <- inla.mesh.2d(boundary = bnd, max.edge = 1)
##!out <- inla.mesh.assessment(mesh, spatial.range = 3, alpha = 2)
##!}



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
                              mesh$loc[mesh$graph$tv[, ti1], ]) ^ 2) ^ 0.5
      i <- c(i, mesh$graph$tv[, ti1], mesh$graph$tv[, ti2])
      val <- c(val, rep(as.vector(len), times = 2))
      num <- c(num, rep(1, 2 * length(as.vector(len))))
    }
    avg_len <-
      as.vector(sparseMatrix(i = i, j = rep(1, length(i)), x = val)) /
        as.vector(sparseMatrix(i = i, j = rep(1, length(i)), x = num))

    proj_len <- as.vector(proj$proj$A %*% avg_len)
    proj_len[!proj$proj$ok] <- NA
    proj_len
  }
  mesh.proj <- function(mesh, dims) {
    INLA::inla.mesh.projector(mesh, dims = dims)
  }
  mesh.spde <- function(mesh, alpha) {
    INLA::inla.spde2.pcmatern(mesh, alpha = alpha,
                              prior.range = c(1, 0.5),
                              prior.sigma = c(1, 0.5))
  }
  mesh.Q <- function(spde, spatial.range) {
    INLA::inla.spde.precision(spde, theta = log(c(spatial.range, 1)))
  }
  mesh.S <- function(Q) {
    INLA::inla.qinv(Q, reordering = INLA::inla.reorderings())
  }
  mesh.sd <- function(proj, S) {
    v <- Matrix::rowSums(proj$proj$A * (proj$proj$A %*% S))
    v[!proj$proj$ok] <- NA
    matrix(v ^ 0.5, length(proj$x), length(proj$y))
  }
  mesh.sd.deviation.approx <- function(proj, S, sd0) {
    val <- proj$proj$A %*% (
             as.vector(t(proj$proj$A[proj$proj$ok, , drop = FALSE]) %*%
                       as.vector(sd0)[proj$proj$ok]) /
             colSums(proj$proj$A[proj$proj$ok, , drop = FALSE]))
    val[!proj$proj$ok] <- NA
    matrix(
          1 + (as.vector(sd0) - val),
          length(proj$x), length(proj$y)
          )
  }
  mesh.sd.bound <- function(proj, S) {
    INLA::inla.mesh.project(proj, field = diag(mesh.S()) ^ 0.5)
  }

  
  if (mesh$manifold != "R2") {
    warning("inla.mesh.assessment has only been tested on flat 2D manifolds.")
  }
  
  spde <- mesh.spde(mesh = mesh, alpha=alpha)
  Q <- mesh.Q(spde = spde, spatial.range = spatial.range)
  S <- mesh.S(Q = Q)
  proj <- mesh.proj(mesh = mesh, dims = dims)
  sd0 <- mesh.sd(proj = proj, S = S)
  sd.deviation <- mesh.sd.deviation.approx(proj = proj, S = S, sd0 = sd0)
  edgelengths <- mesh.edgelengths(mesh, proj)

  out <- data.frame(x = proj$lattice$loc[,1],
                    y = proj$lattice$loc[,2],
                    sd = as.vector(sd0),
                    sd.dev = as.vector(sd.deviation),
                    edge.len = edgelengths)
  coordinates(out) <- ~ x+y
  out
}
