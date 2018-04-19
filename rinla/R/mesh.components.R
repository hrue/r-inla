## Export: inla.mesh.components


##!\name{inla.mesh.components}
##!\alias{inla.mesh.components}
##!\title{Compute connected mesh subsets}
##!\description{Compute subsets of vertices and triangles in an inla.mesh object
##!  that are connected by edges.
##!}
##!\usage{
##!inla.mesh.components(mesh)
##!}
##!\value{A list with elements \code{vertex} and \code{triangle}, vectors of
##!  integer labels for which connected component they belong, and \code{info},
##!  a \code{data.frame} with columns
##!  \item{component}{Connected component integer label.} 
##!  \item{nV}{The number of vertices in the component.}
##!  \item{nT}{The number of triangles in the component.}
##!  \item{area}{The surface area associated with the component.
##!  Component lables are not comparable across different meshes, but
##!  some ordering stability is guaranteed by initiating each component
##!  from the lowest numbered triangle whenever a new component is initiated.}
##!}
##!\author{Finn Lindgren \email{finn.lindgren@gmail.com}}
##!\seealso{inla.mesh.2d, inla.mesh.create}
##!\examples{
##!# Construct two simple meshes:
##!loc <- matrix(c(0,1,0,1), 2, 2)
##!mesh1 <- inla.mesh.2d(loc = loc, max.edge=0.1)
##!bnd <- inla.nonconvex.hull(loc, 0.3)
##!mesh2 <- inla.mesh.2d(boundary = bnd, max.edge=0.1)
##!
##!# Compute connectivity information:
##!conn1 <- inla.mesh.components(mesh1)
##!conn2 <- inla.mesh.components(mesh2)
##!# One component, simply connected mesh
##!conn1$info
##!# Two disconnected components
##!conn2$info
##!
##!# Extract the subset mesh for the largest component:
##!# (Note: some information is lost, such as fixed segments,
##!#  and boundary edge labels.)
##!maxi <- conn2$info$component[which.max(conn2$info$area)]
##!mesh3 <- inla.mesh.create(loc = mesh2$loc,
##!                          tv = mesh2$graph$tv[conn2$triangle == maxi,,drop=FALSE])
##!}
inla.mesh.components <- function(mesh) {
  vertex <- integer(mesh$n)
  Nt <- nrow(mesh$graph$tv)
  triangle <- integer(Nt)
  ok <- !is.na(mesh$graph$tt)
  tt <- Matrix::sparseMatrix(i = c(which(ok[,1]),
                                   which(ok[,2]),
                                   which(ok[,3])),
                             j = c(mesh$graph$tt[ok[,1], 1],
                                   mesh$graph$tt[ok[,2], 2],
                                   mesh$graph$tt[ok[,3], 3]),
                             x = rep(1, sum(ok)),
                             dims = c(Nt, Nt))
  component <- 0
  while (any(triangle == 0)) {
    component <- component + 1
    tri <- integer(Nt)
    tri[min(which(triangle == 0))] <- 1
    tri_prev <- integer(Nt)
    while (any(tri_prev != tri)) {
      tri_prev <- tri
      tri <- tri | as.vector(tt %*% tri > 0)
    }
    triangle[tri > 0] <- component
    vtx <- sort(unique(as.vector(mesh$graph$tv[tri > 0, ])))
    if (any(vertex[vtx] > 0)) {
      warning(paste0("Corner-only connected triangles detected.\n",
                     "  Vertices = ", paste0(vtx[vertex[vtx] > 0], collapse = ", "), "\n",
                     "  Components = ", paste0(vertex[vtx[vertex[vtx] > 0]], collapse = ", "), "\n",
                     "  New component = ", component, "\n",
                     "  Vertex component information will be inconsistent.", "\n",
                     "  Triangle component information will ignore corner-only connections."))
    }
    vertex[vtx] <- component
  }

  if (component == 0) {
    info <- data.frame(component = integer(0), nV = integer(0), nT = integer(0), area = numeric(0))
  } else {
    fem <- inla.mesh.fem(mesh, order = 1)
    info <- do.call(rbind,
      lapply(sort(unique(triangle)),
        function(x) data.frame(component = x,
                               nV = sum(vertex == x),
                               nT = sum(triangle == x),
                               area = sum(fem$ta[triangle == x]))))
  }

  list(
    vertex = vertex,
    triangle = triangle,
    info = info
  )
}
