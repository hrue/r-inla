## Export: inla.read.graph inla.write.graph
## Export: summary!inla.graph
## Export: plot!inla.graph
## Export: print!inla.graph.summary

##!\name{read.graph}
##!\alias{read.graph}
##!\alias{write.graph}
##!\alias{inla.read.graph}
##!\alias{inla.write.graph}
##!\alias{inla.graph}
##!\alias{summary.inla.graph}
##!\alias{plot.inla.graph}
##!\alias{print.inla.graph.summary}
##!\title{Read and write a graph-object}
##!\description{Construct a graph-object from a file or a matrix; write graph-object to file}
##!\usage{
##!inla.read.graph(..., size.only = FALSE)
##!inla.write.graph(graph, filename = "graph.dat", mode = c("binary", "ascii"), ...)
##!
##!\method{summary}{inla.graph}(object, ...)
##!\method{plot}{inla.graph}(x, y, ...)
##!\method{print}{inla.graph.summary}(x, ...)
##!}
##!\arguments{
##!    \item{filename}{The filename of the graph.}
##!    \item{graph}{An \code{inla.graph}-object, a (sparse) symmetric matrix, a filename containing the graph,
##!                 a list or collection of characters and/or numbers defining the graph, 
##!                 or a neighbours list with class \code{nb} (see \code{spdep::card} and
##!                 \code{spdep::poly2nb} for for details of \code{nb} and an example a function
##!                 returning an \code{nb} object}
##!    \item{mode}{The mode of the file; ascii-file or a (gzip-compressed) binary. Default value depends on 
##!                the inla.option \code{internal.binary.mode} which is default \code{TRUE}; see \code{inla.setOption}.}
##!    \item{object}{An \code{inla.graph} -object} 
##!    \item{x}{An \code{inla.graph} -object} 
##!    \item{y}{Not used}
##!    \item{size.only}{Only read the size of the graph}
##!    \item{...}{Additional arguments. In \code{inla.read.graph},
##!               then it is the graph definition (object, matrix, character, filename),  plus extra arguments.
##!               In \code{inla.write.graph} it is extra arguments to \code{inla.read.graph}.}
##!}
##!\value{
##!    The output of \code{inla.read.graph}, is an \code{inla.graph} object, with elements
##!    \item{n}{is the size of the graph}
##!    \item{nnbs}{is a vector with the number of neigbours}
##!    \item{nbs}{is a list-list with the neigbours}
##!    \item{cc}{list with connected component information
##!        \itemize{
##!            \item{\code{id}}{is a vector with the connected component id for each node (starting from 1)}
##!            \item{\code{n}}{is the number of connected components}
##!            \item{\code{nodes}}{is a list-list of nodes belonging to each connected component}
##!        }
##!    }
##!    Methods implemented for \code{inla.graph} are \code{summary} and \code{plot}.
##!    The method \code{plot} require the libraries \code{Rgraphviz} and \code{graph} from the Bioconductor-project,
##!    see \url{https://www.bioconductor.org}.
##!}
##!\author{Havard Rue \email{hrue@r-inla.org}}
##!\seealso{
##!    \code{\link{inla.spy}}
##!}
##!\examples{
##!## a graph from a file
##!cat("3 1 1 2 2 1 1 3 0\n", file="g.dat")
##!g = inla.read.graph("g.dat")
##!## writing an inla.graph-object to file
##!g.file = inla.write.graph(g, mode="binary")
##!## re-reading it from that file
##!gg = inla.read.graph(g.file)
##!summary(g)
##!##
##!Not run:
##!plot(g)
##!inla.spy(g)
##!## when defining the graph directly in the call, 
##!## we can use a mix of character and numbers
##!g = inla.read.graph(c(3, 1, "1 2 2 1 1 3", 0))
##!inla.spy(c(3, 1, "1 2 2 1 1 3 0"))
##!inla.spy(c(3, 1, "1 2 2 1 1 3 0"),  reordering=3:1)
##!inla.write.graph(c(3, 1, "1 2 2 1 1 3 0"))
##!
##!## building a graph from adjacency matrix
##!adjacent = matrix(0, nrow = 4, ncol = 4)
##!adjacent[1,4] = adjacent[4,1] = 1
##!adjacent[2,4] = adjacent[4,2] = 1
##!adjacent[2,3] = adjacent[3,2] = 1
##!adjacent[3,4] = adjacent[4,3] = 1
##!g = inla.read.graph(adjacent)
##!plot(g)
##!summary(g)
##!End(Not run)
##!}

`inla.graph.binary.file.magic` = function()
{
    ## the value of the first integer (read binary) in a binary
    ## filename. this value must be the same as
    ## 'GMRFLib_BINARY_GRAPH_FILE_MAGIC' in GMRFLib/graph.h

    return (-1L)
}

`inla.add.graph.cc` = function(...)
{
    ## add the cc information to a graph

    args = list(...)
    if (length(args) == 0L) {
        return (NULL)
    }

    ## need this test to avoid infinite recursion, as inla.read.graph also call inla.add.graph.cc!
    if (class(args[[1L]]) == "inla.graph") {
        graph = args[[1L]]
    } else {
        graph = inla.read.graph(...)
    }

    cc = list(id = NA, n = NA, nodes = NA)
    n = graph$n

    if (TRUE) {
        ## need these {...} as we want to define 's' globally inside
        ## these, to prevent copying the object all the time duing the
        ## recursions.
    
        s = integer(n)
        s[] = 0L
        k = 1L
        
        do.visit = function(x) x ## just to avoid warning for missing
                                 ## function 'do.visit' during compile
        do.visit = inla.cmpfun(function(idxs) {
            if (any(s[idxs] == 0L)) {
                which.idxs = idxs[which(s[idxs] == 0L)]
                s[which.idxs] <<- k
                ## its ok to refer to 'graph' here:
                visit.nodes = unique(unlist(lapply(which.idxs, function(x) graph$nbs[[x]])))
                
                ## check which of the visit.nodes that needs to be
                ## visited. although this is done already in the
                ## beginning of this routine, but we do that also here
                ## to reduce the depth of the recursive call
                visit.nodes = visit.nodes[which(s[visit.nodes] == 0L)]
                if (length(visit.nodes) > 0L) {
                    do.visit(visit.nodes)
                }
            }
            return (invisible())
        })

        ## need to allow for larger recursion depth,  temporary
        ex.save = getOption("expressions")
        options(expressions=500000L)
        for (i in 1L:n) {
            if (s[i] == 0L) {
                do.visit(i)
                k = k + 1L
            }
        }
        options(paste("expressions=", ex.save,  sep=""))
        
        cc$id = s
        cc$n = max(s)
        cc$nodes = lapply(1L:cc$n, function(cc.id, cs) sort(which(cc.id == cs)), cs = s)
    }

    graph$cc = cc
    return(graph)
}

`inla.graph.size` = function(...)
{
    return(inla.read.graph(..., size.only = TRUE))
}
    

`inla.read.graph` = function(..., size.only = FALSE)
{
    ## graph is either a filename, a graph-object, a (sparse) matrix,
    ## or a list of integers or strings defining the graph.

    `inla.read.graph.ascii.internal` = function(filename, offset = 0L, size.only = FALSE)
    {
        ## offset it needed if the graph is zero-based, then offset is
        ## set to 1.
        stopifnot(file.exists(filename))

        s = readLines(filename)
        if (length(s) == 0L) {
            return (NULL)
        }
        
        ## remove comment lines
        s = sapply(s, function(x) return (gsub("#.*$", "", x))) #
        ## convert "1 2 3" into 1 2 3
        s = as.integer(unlist(sapply(s, function(x) strsplit(x, "[ \t]+"))))
        ## remove possibe NA's that might appear due to spaces at the end of the file
        s = s[!is.na(s)]
    
        n = s[1L]
        if (size.only) {
            return(n)
        }

        g = list(n = n, nnbs = numeric(n), nbs = rep(list(numeric()), n))

        k = 2L
        for(i in 1L:n) {

            if (s[k] + offset == 0L) {
                ## this is a zero-based graph
                return (inla.read.graph.ascii.internal(filename, offset=1L))
            }

            stopifnot(s[k] + offset >= 1L && s[k] + offset <= n)
            idx = s[k] + offset
            k = k+1L

            g$nnbs[idx] = s[k]
            k = k+1L

            if (g$nnbs[idx] > 0L) {
                g$nbs[[idx]] = s[k:(k + g$nnbs[idx] -1L)] + offset
                k = k + g$nnbs[idx]
            }
        }

        stopifnot(k -1L == length(s))
        class(g) = "inla.graph"
        if (length(g$nbs) < g$n) {
            g$nbs = c(g$nbs, rep(list(numeric()), g$n - length(g$nbs)))
        }
        g = inla.add.graph.cc(g)

        return (g)
    }

    `inla.read.graph.binary.internal` = function(filename, offset=0L, size.only = FALSE)
    {
        ## offset it needed if the graph is zero-based, then offset is
        ## set to 1.

        ## read the binary filename, which is the output from inla().
        stopifnot(file.exists(filename))

        ## read the first int,  and check that its the key.
        fp = gzfile(filename, "rb")
        s = as.integer(readBin(fp, integer(), n = 1L))
        close(fp)
        if (length(s) == 0L || s[1L] != inla.graph.binary.file.magic()) {
            ## then its not a binary filename
            return (NULL)
        }
    
        ## since we're using gzfiles (as GMRFLib do that), we don't know
        ## how many elements this file contains from looking at the
        ## size. so we got to try to read to many simply...
        n.try = 2^12
        while(TRUE) {
            fp = gzfile(filename, "rb")
            s = as.integer(readBin(fp, integer(), n = n.try))
            close(fp)
            if (length(s) < n.try) {
                break
            } else {
                n.try = n.try * 4L
            }
        }

        ## remove the key
        s = s[-1L]

        ## then the rest is the graph
        n = s[1L]
        if (size.only) {
            return (n)
        }
        g = list(n = n, nnbs = numeric(n), nbs = rep(list(numeric()), n))

        ## graphs are always 1-based by definition
        k = 2L
        for(i in 1L:n) {

            if (s[k] + offset == 0L) {
                ## this is a zero-based graph
                return (inla.read.graph.binary.internal(filename, offset=1L))
            }

            stopifnot(s[k] + offset >= 1L && s[k] + offset <= n)
            idx = s[k] + offset
            k = k+1L
        
            g$nnbs[idx] = s[k]
            k = k+1L

            if (g$nnbs[idx] > 0L) {
                g$nbs[[idx]] = s[k:(k + g$nnbs[idx] -1L)] + offset
                k = k + g$nnbs[idx]
            }
        }
        stopifnot(k -1L == length(s))
        class(g) = "inla.graph"
        if (length(g$nbs) < g$n) {
            g$nbs = c(g$nbs, rep(list(numeric()), g$n - length(g$nbs)))
        }
        g = inla.add.graph.cc(g)

        return (g)
    }

    `inla.matrix2graph.internal` = function(Q, size.only = FALSE)
    {
        if (missing(Q)) {
            return (NULL)
        }
    
        n = dim(Q)
        if (size.only) {
            return(n[1L])
        }
        Q = inla.as.dgTMatrix(Q)

        if (n[1] != n[2]) {
            stop(paste("Matrix must be a square matrix, dim(Q) =", dim(Q)))
        }

        n = dim(Q)[1]
        g = list(n = n, nnbs = numeric(n), nbs = rep(list(numeric()), n), graph.file = NA)

        if (TRUE) {
            diag(Q) = 1
            Q = inla.as.sparse(Q) ## to avoid possible duplicates
            ord = order(Q@i)
            Q@i = Q@i[ord]
            Q@j = Q@j[ord]
            Q@x = Q@x[ord]
            hash.len = table(Q@i)
            hash.idx = c(1L, 1L+cumsum(hash.len))
            stopifnot(length(hash.len) == ncol(Q))

            for(i in 1L:n) {
                if (hash.len[i] > 1L) {
                    idx = hash.idx[i]:(hash.idx[i] + hash.len[i] - 1L)
                    j = Q@j[idx] + 1L
                    x = Q@x[idx]
                    j = j[ (x != 0.0) & (j != i) ]
                } else {
                    j = NULL
                }
                g$nbs[[i]] = j
                g$nnbs[i] = length(j)
            }
        } else {
            if (TRUE) {
                ## new improved version, using apply. DO NOT PASS 'Q' as argument,  slower...
                g$nbs = lapply(1L:n,
                        inla.cmpfun(function(i) {
                            ## inline: row = inla.sparse.get(Q, row = i)
                            idx = which(Q@i == i-1L)
                            row = list(i = i, j = Q@j[idx] + 1L, values = Q@x[idx])
                            if (length(row$j) > 0) {
                                row$j = row$j[ (row$values != 0.0) & (row$j != i) ]
                            }
                            return (row$j)
                        }))
                g$nnbs = sapply(g$nbs, length)
            } else {
                ## keep old version...
                for(i in 1L:n) {
                    row = inla.sparse.get(Q, row = i)
                    nb = length(row$j)
                    if (nb > 0) {
                        ## setting elements of a sparse-matrix to 0 does not
                        ## necessarily remove that entry.
                        row$j = row$j[ (row$values != 0.0) & (row$j != i) ]
                        nb = length(row$j)
                    }
                    g$nnbs[i] = nb
                    if (g$nnbs[i] > 0L) {
                        g$nbs[[i]] = row$j
                    }
                }
            }
        }        
        class(g) = "inla.graph"
        if (length(g$nbs) < g$n) {
            g$nbs = c(g$nbs, rep(list(numeric()), g$n - length(g$nbs)))
        }
        g = inla.add.graph.cc(g)
        return (g)
    }

    ##
    ## code starts here, really...
    ##

    args = list(...)
    graph = args[[1L]]

    if (is.character(graph) || length(args) > 1L ||
        (is.numeric(graph) && !(is.matrix(graph) || is(graph, "Matrix")))) {
        graph = paste(as.character(graph))

        ## if the file exists, its a file
        if (length(graph) == 1L && file.exists(graph)) {
            ## try binary first, if it fail, try ascii
            g = inla.read.graph.binary.internal(..., size.only = size.only)
            if (is.null(g)) {
                g = inla.read.graph.ascii.internal(..., size.only = size.only)
            }
            return (g)
        } else {
            ## otherwise, its the definition itself
            tfile = tempfile()
            cat(unlist(args), sep = "\n", file=tfile,  append=FALSE)
            ## recursive call
            g = inla.read.graph(tfile, size.only = size.only)
            unlink(tfile)
            return (g)
        }
    } else if (inherits(graph, "inla.graph")) {
        ## no need to do anything. 
        if (size.only) {
            return (graph$n)
        } else {
            return (graph)
        }
    } else if (inherits(graph, "nb")) {
        ## a neigbour-graph from spdep with class="nb".
        ## this can replace spdep::nb2INLA.
        ## call spdep::nb2listw and use spdep coercion.
        inla.require("spdep")
        Q = spdep::nb2listw(graph, style="B", zero.policy=TRUE)
        Q = inla.as.sparse(as(Q, "symmetricMatrix"))
        return (inla.matrix2graph.internal(Q, size.only = size.only))
    } else {
        return (inla.matrix2graph.internal(..., size.only = size.only))
    }

    stopifnot(FALSE)
    return (NULL)
}

`inla.write.graph` = function(graph, filename = "graph.dat", mode = c("binary", "ascii"), ...)
{
    `inla.write.graph.ascii.internal` = function(graph, filename = "graph.dat")
    {
        ## write a graph read from inla.read.graph, or in that format, to
        ## file.
        fd = file(filename , "w")
        cat(graph$n, "\n", file = fd)
        for(i in 1:graph$n) {
            cat(i, graph$nnbs[i], graph$nbs[[i]], "\n", file = fd)
        }
        close(fd)
        return (filename)
    }

    `inla.write.graph.binary.internal` = function(graph, filename = "graph.dat")
    {
        ## write a graph to file,  1-based binary format.
        fd = file(filename , "wb")
        writeBin(as.integer(inla.graph.binary.file.magic()), fd)
        writeBin(as.integer(graph$n), fd)
        if (graph$n > 0L) {
            for(i in 1:graph$n) {
                writeBin(as.integer(i), fd)
                writeBin(as.integer(graph$nnbs[i]), fd)
                if (graph$nnbs[i] > 0L) {
                    writeBin(as.integer(graph$nbs[[i]]), fd)
                }
            }
        }
        close(fd)
        return (filename)
    }

    ##
    ## code starts here
    ##
    if (missing(mode)) {
        ## if nothing is giving, use the global option depending on
        ## the internal.binary.mode. this option is default TRUE, but
        ## can be set to FALSE to ease debugging.
        mode = inla.ifelse(inla.getOption("internal.binary.mode"), "binary", "ascii")
    }
    mode = match.arg(mode)
    g = inla.read.graph(graph, ...)

    if (mode == "binary") {
        return (invisible(inla.write.graph.binary.internal(g, filename)))
    } else if (mode == "ascii") {
        return (invisible(inla.write.graph.ascii.internal(g, filename)))
    } else {
        stopifnot(FALSE)
    }
}

`plot.inla.graph` = function(x, y, ...)
{
    ## these are default options to plot for class inla.graph
    filter = filter.args = c("neato", "dot", "fdp", "twopi")
    attrs = NULL
    scale = 0.5
    node.names = NULL
    ## we evaluate them here,  as they are set in '...'
    inla.eval.dots(...)

    ## I add here some tools to view and summarize a such graphs...
    inla.require("Rgraphviz") || stop("Need library 'Rgraphviz' from Bioconductor: see https://www.bioconductor.org")
    inla.require("graph") || stop("Need library 'graph' from Bioconductor: see https://www.bioconductor.org")

    filter = match.arg(filter, filter.args)
    if (is.null(attrs)) {
        attrs = Rgraphviz::getDefaultAttrs(layoutType = filter)
    }
    if (!is.null(node.names)) {
        stopifnot(length(node.names) == x$n)
    } else {
        node.names = as.character(1:x$n)
    }
    g <- new("graphNEL", nodes = node.names, edgemode = "undirected")
    for (i in 1L:x$n) {
        if (x$nnbs[i] > 0L) {
            j = x$nbs[[i]]
            j = j[j > i]
            if (length(j) > 0L) {
                g = graph::addEdge(node.names[i], node.names[j], g)
            }
        }
    }
    attrs$node$height = as.numeric(attrs$node$height) * scale
    attrs$node$width = as.numeric(attrs$node$width) * scale
    plot(g, filter, attrs = attrs, ...)
}

`summary.inla.graph` = function(object, ...)
{
    ret = list()
    ret = c(ret, list(n = object$n))
    if (!is.null(object$cc)) {
        ret = c(ret, list(ncc = object$cc$n))
    } else {
        ret = c(ret, list(ncc = NA))
    }
    ret = c(ret, list(nnbs = table(object$nnbs)))

    class(ret) = "inla.graph.summary"
    return(ret)
}

`print.inla.graph.summary` = function(x, ...)
{
    cat(paste("\tn = ",  x$n, "\n"))
    cat(paste("\tncc = ",  x$ncc, "\n"))
    w = max(nchar(names(x$nnbs)))
    cat(inla.paste(c("\tnnbs = (names) ",  format(names(x$nnbs), width = w, justify = "right"), "\n")))
    cat(inla.paste(c("\t       (count) ",  format(x$nnbs, width = w, justify = "right"), "\n")))
    return(invisible())
}
