## various graph convert utilities that I got from someone. I don't
## know how useful they are really...

`geobugs2inla` = function(adj, num, graph.file="graph.dat")
{
    return (inla.graph.convert.2(adj, num, graph.file))
}

`inla.graph.convert.2` = function(adj, num, graph.file="graph.dat")
{
    ## A function for converting GeoBUGS adjacency data into the INLA
    ## graph format. Kindly provided by Aki Havunlinna tkk.fi; thanks.

    fd = file(graph.file,  "w")
    len <- length(num)
    cat(len, '\n', file=fd)
    k = 1L
    for(i in 1L:len) {
        if (num[i] > 0) {
            cat(i, num[i], adj[k:(k+num[i]-1L)], "\n", file = fd)
        } else {
            cat(i, num[i], "\n", file=fd)
        }
        k = k + num[i]
    }
    close(fd)
    return (graph.file)
}
