## This demo is similar to demo(Bym), but demonstrate various ways to
## spesify and deal with the 'graph'

data(Germany)
g = system.file("demodata/germany.graph", package="INLA")
source(system.file("demodata/Bym-map.R", package="INLA"))
summary(Germany)

## just make a duplicated column
Germany = cbind(Germany,region.struct=Germany$region)

## graph given as a filename

## we can view the file as a binary matrix, where G[i, j] = 1 if i~j
## and i!=j,  and G[i, i] = 1

## view it as a matrix
inla.dev.new()
inla.spy(g)

formula1 = Y ~ f(region.struct,model="besag",graph=g) +
               f(region,model="iid")
result1  =  inla(formula1,family="poisson",data=Germany,E=E)

inla.dev.new()
Bym.map(result1$summary.random$region.struct$mean)
title("Using graph from file")

## graph given as a inla.graph-object. we can create this object, by
## either reading a graph definition from file, or to specify the
## neigbours in a (sparse) symmetric matrix

go = inla.read.graph(g)
## short summary of its properties
summary(go) 
## view it as a matrix
inla.spy(go)

## view it as a graph. this require the Rgraphviz package from the
## bioconductor project.
if (require(Rgraphviz)) {
    inla.dev.new()
    plot(go)
}

# standard BYM model. the graph is spesified as a filename
formula1 = Y ~ f(region.struct,model="besag",graph=go) +
               f(region,model="iid")
result1  =  inla(formula1,family="poisson",data=Germany,E=E)

inla.dev.new()
Bym.map(result1$summary.random$region.struct$mean)
title("Using graph from an inla.graph-object")

## graph given as a binary matrix (with zero or non-zero term). we can
## also specify the graph as a binary (sparse) symmetric matrix.

## here,  we use the graph we already have
gm = inla.graph2matrix(g)
## or
gm = inla.graph2matrix(go)

## and we can view the result
inla.dev.new()
inla.spy(gm)

## 'gm' can also be a an ordinary dense matrix, but does not have to
if (FALSE)
    gm = as.matrix(gm)

# standard BYM model. the graph is spesified as a filename
formula1 = Y ~ f(region.struct,model="besag",graph=gm) +
               f(region,model="iid")
result1  =  inla(formula1,family="poisson",data=Germany,E=E)

inla.dev.new()
Bym.map(result1$summary.random$region.struct$mean)
title("Using graph from an neighbour-matrix")

## the matrix 'gm' is easier to modify manually, like adding
## neigbhours. recall that this matrix has to be symmetric.

## these two nodes are selected by 'random'. 
stopifnot(gm[1, 10] == 0)
gm[1, 10] = gm[10, 1] = 1

## to rebuild the inla.graph-object, we just re-read it from the
## matrix.
go = inla.read.graph(gm)

## write a new file
fnm = tempfile()
inla.write.graph(gm, filename = fnm)

## and we can use the modified matrix, 'gm', or the inla.graph-object,
## 'go', or the filename, 'fnm', as the as the graph-argument.
formula1 = Y ~ f(region.struct,model="besag",graph=gm) +
               f(region,model="iid")
result1  =  inla(formula1,family="poisson",data=Germany,E=E)

inla.dev.new()
Bym.map(result1$summary.random$region.struct$mean)
title("Using modified graph")

## lets make a random graph, using a dense matrix, just to illustrate...
n = dim(gm)[1L]
gr = matrix(0, n, n)
for(i in 1L:n) {
    ## just sample at random some neigbhours, can include 'i'
    ## itself. so we remove those afterwards.
    js = sample(1L:n, rpois(1, lambda = 10))
    wi = which(js == i)
    if (length(wi) > 0L) 
        js = js[-wi]
    gr[i, js] = 1
    gr[js, i] = 1
}

inla.dev.new()
inla.spy(gr)

formula1 = Y ~ f(region.struct,model="besag",graph=gr) +
               f(region,model="iid")
result1  =  inla(formula1,family="poisson",data=Germany,E=E)

## clean up
unlink(fnm)
