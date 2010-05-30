## obsolete functions goes here.

`inla.obsolete` = function(old, new) {
    stop(inla.paste(c("INLA-function `", as.character(old)[1], "' is obsolete. Use function `",
                      as.character(new), "' instead."), sep=""))
}

`inla.Cmatrix2matrix` = function(C, symmetric = TRUE)
{
    inla.obsolete(match.call(), "inla.sparse2matrix")
    return (invisible())
}

`inla.matrix2Cmatrix` = function(Q, symmetric = TRUE)
{
    inla.obsolete(match.call(), "inla.matrix2sparse")
    return (invisible())
}
    
`inla.Cmatrix2file` = function(Cmatrix, filename = NULL, c.indexing = FALSE)
{
    inla.obsolete(match.call(), "inla.sparse2file")
    return (invisible())
}

