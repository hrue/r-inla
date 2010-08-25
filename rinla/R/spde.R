## create a spde model using fmesher.

`inla.fmesher.call` = function(
  fmesher.call = inla.getOption("fmesher.call"),
  all.args,prefix)
  {
    if (inla.os("linux") || inla.os("mac")) {
        echoc = system(paste(shQuote(fmesher.call), all.args, prefix))
    }
    else if(inla.os("windows")) {
        if (!remote) {
            echoc = try(system(paste(shQuote(fmesher.call),
                                     all.args,
                                     prefix)), silent=TRUE)
            echoc = 0
        } else {
            echoc = try(inla.cygwin.run.command(
                    paste(inla.cygwin.map.filename(fmesher.call),
                          all.args,
                          inla.cygwin.map.filename(prefix))), silent=TRUE)
            echoc = 0
        }
    }
    else
        stop("\n\tNot supported architecture.")
  }

`inla.fmesher.make.prefix` = function(dir = tempdir(), prefix = NULL)
{
  if (is.null(prefix)) {
    if (is.null(dir))
      prefix = paste(tempfile(), ".", sep="")
    else {
      ## create the directory if not already there.
      if (!inla.is.dir(dir)) {
        if (!dir.create(dir, recursive=TRUE)) {
          stop(paste("Failed to create directory:", dir))
        }
      }
      prefix = paste(tempfile(tmpdir=dir), ".", sep="")
    }
  } else {
    ## Make sure the prefix directory exists.
    if (!inla.is.dir(inla.dirname(prefix))) {
      if (!dir.create(inla.dirname(prefix), recursive=TRUE)) {
        stop(paste("Failed to create directory:", inla.dirname(prefix)))
      }
    }
  }
  return (prefix)
}

`inla.create.spde` = function(mesh=NULL, prefix = NULL,
  fem=NULL, sph0=NULL, sph=NULL)
  {
    ## create files using fmesher

    if (is.null(mesh) && is.null(prefix))
      stop("At least one of mesh and prefix must be specified.")

    if (is.null(mesh)) {
      ## Need to know the size of the graph.
      s = inla.read.fmesher.file(paste(prefix, "s", sep=""))
      stopifnot(is.matrix(s))
      n = nrow(s)
      s.dim = ncol(s)
    } else {
      stopifnot(is.matrix(mesh$s))
      n = nrow(mesh$s)
      s.dim = ncol(mesh$s)
    }

    if (s.dim==1)
      stop("1-D models not implemented yet.")
    stopifnot(s.dim>=2, s.dim<=3)

    prefix = inla.fmesher.make.prefix(NULL,prefix)

    if (!is.null(mesh)) {
      stopifnot(is.matrix(mesh$tv))
      inla.write.fmesher.file(mesh$s, filename = paste(prefix, "s", sep=""))
      inla.write.fmesher.file(inla.affirm.integer(mesh$tv-1L),
                              filename = paste(prefix, "tv", sep=""))
    }
    all.args = "--input=s,tv"

    ## additional argumets
    if (!is.null(fem)) {
      all.args = paste(all.args," --fem=",fem,sep="")
    }
    if (!is.null(sph0)) {
      all.args = paste(all.args," --sph0=",sph0,sep="")
    }
    if (!is.null(sph)) {
      all.args = paste(all.args," --sph=",sph,sep="")
    }
    all.args = paste(all.args, inla.getOption("fmesher.arg"))

    echoc = inla.fmesher.call(all.args=all.args,prefix=prefix)
    
    return (list(prefix = prefix, locations = mesh$s, n = n))
}





`inla.fmesher.mesh` = function(
  locations,
  boundary = NULL,
  dir = tempdir(),
  prefix = NULL,
  rcdt = c(21,-1.0,-0.5),
  x11.delay = NULL,
  min.input.distance = 0)
{

    ## create files using fmesher

    stopifnot(is.matrix(locations))
    loc.n = dim(locations)[1]
    loc.dim = dim(locations)[2]
    stopifnot(loc.dim>=2, loc.dim<=3)

    if (!is.null(boundary))
      stop("Boundary handling is not implemented in the R interface yet.")

    prefix = inla.fmesher.make.prefix(dir,prefix)

    ## Map locations to nodes, avoiding near-duplicates.
    node.coord = matrix(nrow=loc.n,ncol=loc.dim)
    map.loc.to.node = rep(0,nrow=loc.n)
    excluded = c()
    loc.i = 1
    node.i.max = 1
    node.coord[node.i.max,] = locations[loc.i,]
    map.loc.to.node[[loc.i]] = node.i.max
    for (loc.i in 2:loc.n) {
      loc.to.node.dist = 
        sqrt(rowSums((as.matrix(rep(1,node.i.max)) %*% locations[loc.i,] -
                      node.coord[1:node.i.max,])^2))
      if (min(loc.to.node.dist) > min.input.distance) {
        node.i.max = node.i.max+1
        node.coord[node.i.max,] = locations[loc.i,]
        map.loc.to.node[[loc.i]] = node.i.max
      } else {
        excluded = c(excluded,loc.i)
      }
    } 
    ## Remove excess nodes.
    node.coord = node.coord[1:node.i.max,]
    ## Identify nearest nodes for excluded locations.
    for (loc.i in excluded) {
      loc.to.node.dist = 
        sqrt(rowSums((as.matrix(rep(1,node.i.max)) %*% locations[loc.i,] -
                      node.coord)^2))
      node.i = which.min(loc.to.node.dist);
      map.loc.to.node[loc.i] = node.i
    }

    s0 = node.coord
    idx = map.loc.to.node

    ## The default name for the input locations in fmesher is s0:
    loc.file = inla.write.fmesher.file(s0,
      filename = paste(prefix, "s0", sep=""))

    ## additional arguments
    all.args = ""
    if (!is.null(rcdt)) {
      all.args = paste(all.args," --rcdt=",rcdt[1],",",rcdt[2],",",rcdt[3],sep="")
    }
    if (!is.null(x11.delay)) {
      all.args = paste(all.args," --x11=",x11.delay,sep="")
    }
    all.args = paste(all.args, inla.getOption("fmesher.arg"))

    echoc = inla.fmesher.call(all.args=all.args,prefix=prefix)
    
    mesh = list(
      tv = 1L+inla.read.fmesher.file(paste(prefix, "tv", sep="")),
      s = inla.read.fmesher.file(paste(prefix, "s", sep="")))
    
    return (list(prefix = prefix, mesh = mesh, locations.idx = idx))
}
