##!\name{read.graph}
##!\alias{read.graph}
##!\alias{write.graph}
##!\alias{inla.read.graph}
##!\alias{inla.write.graph}
##!\alias{inla.graph}
##!\title{Read and write a graph-object}
##!\description{Reads a graph-object to a file and write graph-object to file}
##!\usage{
##!g = inla.read.graph(graph)
##!inla.write.graph(graph, filename, mode = c("binary", "ascii""))
##!summary(graph)
##!plot(graph)
##!\arguments{
##!    \item{filename}{The filename of the graph.}
##!    \item{graph}{An \code{inla.graph}-object, a (sparse) symmetric matrix or a filename containing the graph.}
##!    \item{mode}{The mode of the file; ascii-file or a (gzip-compressed) binary. Default value depends on 
##!                the inla.option \code{internal.binary.mode} which is default \code{TRUE}; see \code{inla.setOption}.}
##!}
##!\value{
##!    The output of \code{inla.read.graph}, is an \code{inla.graph} object, with elements
##!    \item{n}{is the size of the graph}
##!    \item{nnbs}{is a vector with the number of neigbours}
##!    \item{nbs}{is a list-list with the neigbours}
##!    \item{filename}{is the name of the file if this \code{inla.graph}-object was read from a file}
##!    \item{cc}{list with connected component information (this entry can be auto-generated; see below)
##!        \itemize{
##!            \item{\code{id}}{is a vector with the connected component id for each node (starting from 1)}
##!            \item{\code{n}}{is the number of connected components}
##!            \item{\code{nodes}}{is a list-list of nodes belonging to each connected component}
##!        }
##!    }
##!    The connected component information,  can be generated from the rest of the graph-structure,
##!    using \code{graph = inla.add.graph.cc(graph)}.
##!    Methods implemented for \code{inla.graph} are \code{summary} and \code{plot}.
##!    The method \code{plot} require the libraries \code{Rgraphviz} and \code{graph} from the Bioconductor-project,
##!    see \code{http://www.bioconductor.org}.
##!
##!    The \code{graph} is either a filename containing the \code{inla.graph}-object (previously written using
##!    \code{inla.write.graph()},  or an \code{inla.graph}-object,  or a (sparse) symmetric matrix whos non-zero pattern
##!    determines the graph.
##!}
##!\author{Havard Rue \email{hrue@math.ntnu.no}}
##!\seealso{
##!    \code{\link{inla.matrix2graph}}
##!    \code{\link{inla.graph2matrix}} 
##!    \code{\link{inla.spy}}
##!}
##!\examples{
##!cat("3 1 1 2 2 1 1 3 0\n", file="g.dat")
##!g = inla.read.graph("g.dat")
##!g.file = inla.write.graph(g, mode="binary")
##!gg = inla.read.graph(g.file)
##!summary(g)
##!summary(gg)
##!plot(gg)
##!inla.spy(gg)
##!ggg = inla.graph2matrix(gg)
##!ggg.file = inla.write.graph(ggg)
##!inla.spy(ggg.file)
##!}

`inla.graph.binary.file.magic` = function()
{
    ## the value of the first integer (read binary) in a binary
    ## filename. this value must be the same as the 'define
    ## GMRFLib_BINARY_GRAPH_FILE_MAGIC' in GMRFLib/graph.h

    return (-1L)
}

`inla.add.graph.cc` = function(graph)
{
    ## add the cc information to a graph
    
    cc = list(id = NA, n = NA, nodes = NA)
    n = graph$n

    if (TRUE) {
        ## need these {...} as we want to define 's' globally inside
        ## these, to prevent copying the object all the time duing the
        ## recursions.
    
        s = integer(n)
        s[] = 0L
        k = 1L

        do.visit = function(idx) {
            if (s[idx] == 0L) {
                ## this is ok,  as it picks up 's' from within the if(TRUE) {...} 
                s[idx] <<- k
                if (graph$nnbs[idx] > 0L) {
                    for(iidx in graph$nbs[[idx]]) {
                        do.visit(iidx)
                    }
                }
            } 
            return (invisible())
        }

        for (i in 1:n) {
            if (s[i] == 0L) {
                do.visit(i)
                k = k + 1L
            }
        }
    
        cc$id = s
        cc$n = max(s)
        cc$nodes = lapply(1L:cc$n, function(cc.id, cs) which(cc.id == cs), cs = s)
    }

    graph$cc = cc
    return(graph)
}

`inla.read.graph` = function(graph)
{
    ## graph is either a filename, a graph-object, or a (sparse) matrix.

    `inla.read.graph.ascii` = function(filename)
    {
        stopifnot(file.exists(filename))

        s = readLines(filename)
        if (length(s) == 0L) {
            return (NULL)
        }
        
        ## remove comment lines
        s = sapply(s, function(x) return (gsub("#.*$", "", x))) #
        ## convert "1 2 3" into 1 2 3
        s = as.integer(unlist(sapply(s, function(x) strsplit(x, " "))))
        ## remove possibe NA's that might appear due to spaces at the end of the file
        s = s[!is.na(s)]
    
        n = s[1L]
        g = list(n = n, nnbs = numeric(n), nbs = rep(list(numeric()), n), filename = filename)

        k = 2L
        for(i in 1L:n) {

            if (s[k] == 0L) {
                stop("The nodes of the graph must be numbered from 1...n,  not 0...n-1. Please fix.")
            }

            stopifnot(s[k] >= 1L && s[k] <= n)
            idx = s[k]
            k = k+1L

            g$nnbs[idx] = s[k]
            k = k+1L

            if (g$nnbs[idx] > 0L) {
                g$nbs[[idx]] = s[k:(k + g$nnbs[idx] -1L)]
                k = k + g$nnbs[idx]
            }
        }

        stopifnot(k -1L == length(s))
        class(g) = "inla.graph"
        g = inla.add.graph.cc(g)

        return (g)
    }

    `inla.read.graph.binary` = function(filename)
    {
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
        g = list(n = n, nnbs = numeric(n), nbs = rep(list(numeric()), n), filename = filename)

        ## graphs are always 1-based by definition
        k = 2L
        for(i in 1L:n) {

            if (s[k] == 0L) {
                stop("The nodes of the graph must be numbered from 1...n,  not 0...n-1. Please fix.")
            }

            stopifnot(s[k] >= 1L && s[k] <= n)
            idx = s[k]
            k = k+1L
        
            g$nnbs[idx] = s[k]
            k = k+1L

            if (g$nnbs[idx] > 0L) {
                g$nbs[[idx]] = s[k:(k + g$nnbs[idx] -1L)] 
                k = k + g$nnbs[idx]
            }
        }
        stopifnot(k -1L == length(s))
        class(g) = "inla.graph"
        g = inla.add.graph.cc(g)

        return (g)
    }

    `inla.matrix2graph` = function(Q)
    {
        if (missing(Q)) {
            return (NULL)
        }
    
        if (!is(Q, "dgTMatrix")) {
            Q = inla.as.dgTMatrix(Q)
        }
        n = dim(Q)

        if (n[1] != n[2]) {
            stop(paste("Matrix must be a square matrix, dim(Q) =", dim(Q)))
        }

        n = dim(Q)[1]
        g = list(n = n, nnbs = numeric(n), nbs = rep(list(numeric()), n), graph.file = NA)

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
        class(g) = "inla.graph"
        g = inla.add.graph.cc(g)
    
        return (g)
    }

    ##
    ## code starts here, really...
    ##
    if (is.character(graph)) {
        stopifnot(file.exists(graph))
         ## try binary first, if it fail, try ascii
        g = inla.read.graph.binary(graph)
        if (is.null(g)) {
            g = inla.read.graph.ascii(graph)
        }
        return (g)
    } else if (class(graph) == "inla.graph") {
        ## no need to do anything. 
        return (graph)
    } else {
        return (inla.matrix2graph(graph))
    }

    stopifnot(FALSE)
    return (NULL)
}

`inla.write.graph` = function(graph,  filename = "graph.dat", mode = c("binary", "ascii"))
{
    `inla.write.graph.ascii` = function(graph, filename = "graph.dat")
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

    `inla.write.graph.binary` = function(graph, filename = "graph.dat")
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

    if (is.character(graph)) {
        g = inla.read.graph(graph)
    } else if (class(graph) == "inla.graph") {
        g = graph
    } else {
        g = inla.matrix2graph(graph)
    }

    if (mode == "binary") {
        return (invisible(inla.write.graph.binary(g, filename)))
    } else if (mode == "ascii") {
        return (invisible(inla.write.graph.ascii(g, filename)))
    } else {
        stopifnot(FALSE)
    }
}

`plot.inla.graph` = function(
        graph,
        filter = c("neato", "fdp"),
        attrs = NULL,
        scale = 0.5,
        node.names = NULL,
        ...)
{
    ## I add here some tools to view and summarize a such graphs...
    require(Rgraphviz) || stop("Need library 'Rgraphviz' from Bioconductor: see http://www.bioconductor.org")
    require(graph) || stop("Need library 'graph' from Bioconductor: see http://www.bioconductor.org")

    filter = match.arg(filter)
    if (is.null(attrs)) {
        attrs = getDefaultAttrs(layoutType = filter)
    }
    if (!is.null(node.names)) {
        stopifnot(length(node.names) == graph$n)
    } else {
        node.names = as.character(1:graph$n)
    }
    g <- new("graphNEL", nodes = node.names, edgemode = "undirected")
    for (i in 1L:graph$n) {
        if (graph$nnbs[i] > 0L) {
            j = graph$nbs[[i]]
            j = j[j > i]
            if (length(j) > 0L) {
                g = addEdge(node.names[i], node.names[j], g)
            }
        }
    }
    attrs$node$height = as.numeric(attrs$node$height) * scale
    attrs$node$width = as.numeric(attrs$node$width) * scale
    plot(g, filter, attrs = attrs, ...)
}

`summary.inla.graph` = function(graph, ...)
{
    ret = list()
    ret = c(ret, list(n = graph$n))
    if (!is.null(graph$cc)) {
        ret = c(ret, list(ncc = graph$cc$n))
    } else {
        ret = c(ret, list(ncc = NA))
    }
    ret = c(ret, list(filename = graph$filename))
    ret = c(ret, list(nnbs = table(graph$nnbs)))

    class(ret) = "inla.graph.summary"
    return(ret)
}

`print.inla.graph.summary` = function(go, ...)
{
    cat(paste("\tn = ",  go$n, "\n"))
    cat(paste("\tncc = ",  go$ncc, "\n"))
    cat(paste("\tfilename = ", shQuote(go$filename), "\n"))
    w = max(nchar(names(go$nnbs)))
    cat(inla.paste(c("\tnnbs = (names) ",  format(names(go$nnbs), width = w, justify = "right"), "\n")))
    cat(inla.paste(c("\t       (count) ",  format(go$nnbs, width = w, justify = "right"), "\n")))
    return(invisible())
}
