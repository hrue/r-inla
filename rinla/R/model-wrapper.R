## Nothing to export

## I am not sure if this is a good idea after all...

`inla.model.wrapper` = function(
        covariates = NULL,
        model = c("joint"),
        data = NULL,
        Q = NULL,
        precision = 1.0E8,
        hyper = NULL,
        initial = NULL,
        prior = NULL,
        param = NULL,
        fixed = NULL,
        ...)
{
    model = match.arg(model, several.ok = FALSE)

    if (inla.one.of(model,  "joint")) {

        res = list()
        res$model = "generic0"

        if (is.character(covariates)) {
            if (is.null(data)) {
                covariates = eval.parent(parse(text=covariates))
            } else {
                covariates = eval(parse(text=covariates), data)
            }
        }
        covariates = as.matrix(covariates)
        n = dim(covariates)[1]
        nb = dim(covariates)[2]
        N = n+nb

        if (is.null(Q) || any(dim(Q) != nb)) {
            stop("Matrix Q is either null or has wrong dimension")
        }

        QQ = sparseMatrix(i=1:n, j=1:n, x=rep(precision, n), dims = c(N, N))
        for(i in 1:n) {
            cov.row = -precision * covariates[i, 1:nb, drop=FALSE]
            QQ[(n+1):N, i] = cov.row
            QQ[i, (n+1):N] = t(cov.row)
        }
        QQ[(n+1):N, (n+1):N] = precision * t(covariates) %*% covariates + Q

        hyper = inla.set.hyper(model, "wrapper", hyper, 
                initial, fixed, prior, param)         

        ret = list(f = list(Cmatrix = QQ, hyper.default = hyper,
                           model = "generic", n = N,  values = 1:N, rankdef = n))
        class(ret) = "inla.wrapper.model"
        
        return (ret)
    } else {
        stop("This should not happen.")
    }
}
