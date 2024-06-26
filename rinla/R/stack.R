#' @noRd
rbind.inla.data.stack.info <- function(...) {
    l <- list(...)
    names(l) <- NULL
    names.tmp <- do.call(c, lapply(l, function(x) x$names))
    ncol.tmp <- do.call(c, lapply(l, function(x) x$ncol))
    
    ncol <- c()
    names <- list()
    for (k in 1:length(names.tmp)) {
        name <- names(names.tmp)[k]
        if (!is.null(names[[name]])) {
            if (!identical(
                names[[name]],
                names.tmp[[k]]
            )) {
                stop("Name mismatch.")
            }
        }
        names[[name]] <- names.tmp[[k]]
        
        if (!is.null(as.list(ncol)[[name]])) {
            if (ncol[name] != ncol.tmp[[k]]) {
                stop("ncol mismatch.")
            }
        }
        ncol[name] <- ncol.tmp[[k]]
    }
    
    data <- dplyr::bind_rows(lapply(l, function(x) x[["data"]]))
    
    offset <- 0
    index <- list()
    for (k in 1:length(l)) {
        for (j in 1:length(l[[k]]$index)) {
            if (is.null(index[[names(l[[k]]$index)[j]]])) {
                index[[names(l[[k]]$index)[j]]] <- l[[k]]$index[[j]] + offset
            } else {
                index[[names(l[[k]]$index)[j]]] <-
                    c(
                        index[[names(l[[k]]$index)[j]]],
                        l[[k]]$index[[j]] + offset
                    )
            }
        }
        offset <- offset + l[[k]]$nrow
    }
    
    info <-
        list(
            data = data,
            nrow = nrow(data),
            ncol = ncol,
            names = names,
            index = index
        )
    class(info) <- "inla.data.stack.info"
    
    return(info)
}

ensure.canonical.response.objects <- function(responses) {
    null.l <- vapply(responses, is.null, logical(1))
    responses <- responses[!null.l]
    if (length(responses) == 0) {
        return(list(NULL))
    }
    responses <- lapply(responses, function(x) {
        if (inherits(x, "inla.mdata")) {
            # Make sure mdata is a data.frame
            attribs <- attributes(x)
            x <- as.data.frame(x)
            attr(x, "inla.ncols") <- attribs$inla.ncols
            attr(x, "names.ori") <- attribs$names.ori
            class(x) <- c("inla.mdata", "data.frame")
            x
        } else if (inherits(x, "inla.surv")) {
            # Make sure inla.surv is a data.frame
            attribs <- attributes(x)
            y <- unclass(x)
            y <- x[!vapply(y, is.null, logical(1))]
            x <- as.data.frame(y)
            attr(x, "names.ori") <- attribs$names.ori
            class(x) <- c("inla.surv", "data.frame")
            x
        } else {
            x
        }
    })
    
    responses
}

#' @noRd
rbind.inla.stack.responses <- function(l) {
    l <- ensure.canonical.response.objects(l)    
    
    classes <- lapply(l, class)
    if (length(unique(classes)) > 1) {
        stop("Cannot rbind responses with different classes.")
    }
    
    attribs <- lapply(l, attributes)
    
    if (all(vapply(l, is.data.frame, logical(1)))) {
        response <- dplyr::bind_rows(l)
        class(response) <- classes[[1]]
        if (inherits(response, "inla.mdata")) {
            attr(response, "inla.ncols") <- attribs[[1]]$inla.ncols
            attr(response, "names.ori") <- attribs[[1]]$names.ori
        } else if (inherits(response, "inla.surv")) {
            attr(response, "names.ori") <- attribs[[1]]$names.ori
        }
    } else {
        response <- do.call(c, l)
    }
    
    return(list(response))
}

#' @noRd
expand.inla.stack.responses <- function(l) {
    nrows <- vapply(l, NROW, integer(1))
    
    responses <- lapply(
        seq_along(l),
        function(k) {
            x <- l[[k]]
            if (is.data.frame(x)) {
                y <- dplyr::bind_rows(
                    as.data.frame(
                        matrix(
                            NA,
                            nrow = sum(nrows[seq_len(k - 1)]),
                            ncol = 0
                        )
                    ),
                    x,
                    as.data.frame(
                        matrix(
                            NA,
                            nrow = sum(nrows[-seq_len(k)]),
                            ncol = 0
                        )
                    )
                )
                if (inherits(x, "inla.mdata")) {
                    attr(y, "inla.ncols") <- attr(x, "inla.ncols")
                    attr(y, "names.ori") <- attr(x, "names.ori")
                    class(y) <- c("inla.mdata", "data.frame")
                } else if (inherits(x, "inla.surv")) {
                    attr(y, "names.ori") <- attr(x, "names.ori")
                    class(y) <- c("inla.surv", "data.frame")
                }
            } else if (is.matrix(x)) {
                NA_ <- x[1, 1]
                is.na(NA_) <- TRUE
                y <- rbind(
                    matrix(
                        NA_,
                        nrow = sum(nrows[seq_len(k - 1)]),
                        ncol = ncol(x)
                    ),
                    x,
                    matrix(
                        NA_,
                        nrow = sum(nrows[-seq_len(k)]),
                        ncol = ncol(x)
                    )
                )
            } else if (is.vector(x)) {
                NA_ <- x[1]
                is.na(NA_) <- TRUE
                y <- c(
                    rep(NA_, sum(nrows[seq_len(k - 1)])),
                    x,
                    rep(NA_, sum(nrows[-seq_len(k)])))
            } else {
                stop("Don't know how to expand responses of class '",
                     paste0(class(x), collapse = ", "), "'.")
            }
            y
        })
    
    
    return(responses)
}

#' @describeIn inla.stack Remove unused entries from an existing stack
#'
#' @export
inla.stack.remove.unused <- function(stack) {
    inla.require.inherits(stack, "inla.data.stack", "'stack'")
    
    if (stack$effects$nrow < 2) {
        return(stack)
    }
    
    ## Remove components with no effect:
    remove <- rep(FALSE, stack$effects$nrow)
    remove.unused.indices <-
        which(colSums(abs(stack$A[, , drop = FALSE])) == 0)
    remove[remove.unused.indices] <- TRUE
    
    index.new <- rep(as.integer(NA), stack$effects$nrow)
    
    ncol.A <- sum(!remove)
    if (ncol.A > 0) {
        index.new[!remove] <- seq_len(ncol.A)
    }
    index.new[remove] <- index.new[index.new[remove]]
    
    for (k in 1:length(stack$effects$index)) {
        stack$effects$index[[k]] <- index.new[stack$effects$index[[k]]]
    }
    
    A <- inla.as.dgTMatrix(stack$A)
    j.new <- index.new[A@j + 1L]
    ## Check for any zero-elements in remove.unused-columns:
    ok <- !is.na(j.new)
    stack$A <-
        sparseMatrix(
            i = A@i[ok] + 1L,
            j = j.new[ok],
            x = A@x[ok],
            dims = c(nrow(A), ncol.A)
        )
    
    stack$effects$data <- stack$effects$data[!remove, , drop = FALSE]
    stack$effects$nrow <- ncol.A
    
    return(stack)
}

#' @describeIn inla.stack Compress an existing stack by removing duplicates
#'
#' @export
inla.stack.compress <- function(stack, remove.unused = TRUE) {
    inla.require.inherits(stack, "inla.data.stack", "'stack'")
    
    if (stack$effects$nrow < 2) {
        return(stack)
    }
    
    ii <- do.call(order, as.list(stack$effects$data))
    jj.dupl <-
        which(1L ==
                  diff(c(
                      duplicated(stack$effects$data[ii, , drop = FALSE]),
                      FALSE
                  )))
    kk.dupl <-
        which(-1L ==
                  diff(c(
                      duplicated(stack$effects$data[ii, , drop = FALSE]),
                      FALSE
                  )))
    ## ii[jj.dupl] are the rows that have duplicates.
    ## ii[(jj.dupl[k]+1):kk.dupl[k]] are the duplicate rows for each k
    
    remove <- rep(FALSE, stack$effects$nrow)
    index.new <- rep(as.integer(NA), stack$effects$nrow)
    
    if (length(jj.dupl) > 0) {
        for (k in 1:length(jj.dupl)) {
            i <- ii[jj.dupl[k]]
            j <- ii[(jj.dupl[k] + 1):kk.dupl[k]]
            
            remove[j] <- TRUE
            index.new[j] <- i
        }
    }
    
    ncol.A <- sum(!remove)
    if (ncol.A > 0) {
        index.new[!remove] <- seq_len(ncol.A)
    }
    index.new[remove] <- index.new[index.new[remove]]
    
    for (k in 1:length(stack$effects$index)) {
        stack$effects$index[[k]] <- index.new[stack$effects$index[[k]]]
    }
    
    A <- inla.as.dgTMatrix(stack$A)
    j.new <- index.new[A@j + 1L]
    ## Check for any zero-elements in remove.unused-columns:
    ok <- !is.na(j.new)
    stack$A <-
        sparseMatrix(
            i = A@i[ok] + 1L,
            j = j.new[ok],
            x = A@x[ok],
            dims = c(nrow(A), ncol.A)
        )
    
    stack$effects$data <- stack$effects$data[!remove, , drop = FALSE]
    stack$effects$nrow <- ncol.A
    
    if (remove.unused) {
        return(inla.stack.remove.unused(stack))
    } else {
        return(stack)
    }
}





#' Data stacking for advanced INLA models
#'
#' Functions for combining data, effects and observation matrices into
#' `inla.stack` objects, and extracting information from such objects.
#'
#' For models with a single effects collection, the outer list container for
#' `A` and `effects` may be omitted.
#'
#' Component size definitions:
#' * \eqn{n_l}{n_l} effect blocks
#' * \eqn{n_k}{n_k} effects
#' * \eqn{n_i}{n_i} data values
#' * \eqn{n_{j,l}}{n_jl} effect size for block \eqn{l}{l}
#' * \eqn{n_j}{n_j} \eqn{= \sum_{l=1}^{n_l} n_{j,l}}{sum_l n_jl} total
#' effect size
#'
#' Input: \describe{
#' \item{`data`}{\eqn{(y^1, \ldots, y^p)}{(y1,\dots,y2)} \eqn{p}{p}
#' vectors, each of length \eqn{n_i}{n_i}}
#' \item{`A`}{\eqn{(A^1, \ldots, A^{n_l})}{(A1,\dots,A2)} matrices of size
#' \eqn{n_i \times n_{j,l}}{n_i by n_jl}}
#' \item{`effects`}{\eqn{\left((x^{1,1},\ldots,x^{n_k,1}), \ldots,
#' (x^{1,n_l},\ldots,x^{n_k,n_l})\right)}{((x_[1,1],\dots,x_[n_k,1]),\dots(x_[1,n_l],\dots,x_[n_k,n_l]))}
#' collections of effect
#' vectors of length \eqn{n_{j,l}}{n_jl} }
#' }
#'
#' \deqn{\mbox{predictor}(y^1, \ldots, y^p) \sim
#' \sum_{l=1}^{n_l} A^l \sum_{k=1}^{n_k} g(k, x^{k,l})
#' = \tilde{A} \sum_{k=1}^{n_k} g(k, \tilde{x}^k)
#' }{ predictor(y^1, \ldots, y^p) ~ sum_{l=1}^{n_l} A^l sum_{k=1}^{n_k}
#' g(k, x^{k,l}) = tilde{A} sum_{k=1}^{n_k} g(k, tilde{x}^k) }
#' where
#' \deqn{\tilde{A} = \mbox{cbind}\left( A^1, \ldots, A^{n_l} \right)
#' }{ tilde{A} = cbind( A^1, ..., A^{n_l} ) }
#' and
#' \deqn{\tilde{x}^k = \mbox{rbind}\left( x^{k,1}, \ldots, x^{k,n_l} \right)
#' }{ tilde{x}^k = rbind( x^{k,1}, ..., x^{k,n_l} ) }
#' and for each block \eqn{l}{l}, any missing
#' \eqn{x^{k,l}} is replaced by an `NA` vector.
#'
#' @aliases inla.stack inla.stack.remove.unused inla.stack.compress
#' inla.stack.sum inla.stack.join inla.stack.index inla.stack.LHS
#' inla.stack.RHS inla.stack.data inla.stack.A
#' @param stack A `inla.data.stack` object, created by a call to
#' `inla.stack`, `inla.stack.sum`, or `inla.stack.join`.
#' @param remove.unused If `TRUE`, compress the model by removing rows of
#' effects corresponding to all-zero columns in the `A` matrix (and
#' removing those columns).
#' @param multi.family logical or character. For `inla.data.join`, if `TRUE`,
#' the `response` part of the stack is joined as a `list`. If `character`,
#' denotes the name of a `data` element that should be joined as a multi-column
#' matrix. Default is `FALSE`, which joins both the `data` and `responses`
#' elements with regular row binding with `dplyr::bind_rows`.
#' @param ... For `inla.stack.join`, two or more data stacks of class
#' `inla.data.stack`, created by a call to `inla.stack`,
#' `inla.stack.sum`, or `inla.stack.join`. For
#' `inla.stack.data`, a list of variables to be joined with the data list.
#' @param compress If `TRUE`, compress the model by removing duplicated
#' rows of effects, replacing the corresponding A-matrix columns with a single
#' column containing the sum.
#' @param data A list or codedata.frame of named data vectors. Scalars are
#' expanded to match the number of rows in the A matrices, or any non-scalar
#' data vectors. An error is given if the input is inconsistent.
#' @param A A list of observation matrices. Scalars are expanded to diagonal
#' matrices matching the effect vector lengths. An error is given if the input
#' is inconsistent or ambiguous.
#' @param effects A collection of effects/predictors.  Each list element
#' corresponds to an observation matrix, and must either be a single vector, a
#' list of vectors, or a `data.frame`. Single-element effect vectors are
#' expanded to vectors matching the number of columns in the corresponding A
#' matrix.  An error is given if the input is inconsistent or ombiguous.
#' @param tag A string specifying a tag for later identification.
#' @return A data stack of class `inla.data.stack`.
#' Elements: \itemize{
#' \item`data` \eqn{=(y^1, \ldots, y^p, \tilde{x}^1, \ldots,
#' \tilde{x}^{n_k})}{=(y^1, \dots, y^p, tilde{x}^1, \dots,
#' tilde{x}^{n_k})}
#' \item`A` \eqn{=\tilde{A}}{=tilde{A}}
#' \item`data.names` List
#' of data names, length \eqn{p}
#' \item`effect.names` List of effect names,
#' length \eqn{n_k}
#' \item`n.data` Data length, \eqn{n_i}
#' \item`index`
#' List indexed by `tag`s, each element indexing into \eqn{i=1, \ldots,
#' n_i} }
#' @section Functions: \itemize{
#' \item `inla.stack.remove.unused`: Remove
#' unused entries from an existing stack
#'
#' \item `inla.stack.compress`: Compress an existing stack by removing
#' duplicates
#'
#' \item `inla.stack`: Shorthand for inla.stack.join and inla.stack.sum
#'
#' \item `inla.stack.sum`: Create data stack as a sum of predictors
#'
#' \item `inla.stack.join`: Join two or more data stacks
#'
#' \item `inla.stack.index`: Extract tagged indices
#'
#' \item `inla.stack.LHS`: Extract data associated with the "left hand
#' side" of the model (e.g. the data itself, `Ntrials`, `link`,
#' `E`)
#'
#' \item `inla.stack.RHS`: Extract data associated with the "right hand
#' side" of the model (all the covariates/predictors)
#'
#' \item `inla.stack.data`: Extract data for an inla call, and optionally
#' join with other variables
#'
#' \item `inla.stack.A`: Extract the "A matrix" for control.predictor }
#' @seealso [inla.spde.make.A()], [inla.spde.make.index()]
#' @keywords fmesher
#' @examples
#'
#' n <- 200
#' loc <- matrix(runif(n * 2), n, 2)
#' mesh <- inla.mesh.2d(
#'     loc.domain = loc,
#'     max.edge = c(0.05, 0.2)
#' )
#' proj.obs <- inla.mesh.projector(mesh, loc = loc)
#' proj.pred <- inla.mesh.projector(mesh, loc = mesh$loc)
#' spde <- inla.spde2.pcmatern(mesh,
#'     prior.range = c(0.01, 0.01),
#'     prior.sigma = c(10, 0.01)
#' )
#'
#' covar <- rnorm(n)
#' field <- inla.qsample(n = 1, Q = inla.spde.precision(spde, theta = log(c(0.5, 1))))[, 1]
#' y <- 2 * covar + inla.mesh.project(proj.obs, field)
#'
#' A.obs <- inla.spde.make.A(mesh, loc = loc)
#' A.pred <- inla.spde.make.A(mesh, loc = proj.pred$loc)
#' stack.obs <-
#'     inla.stack(
#'         data = list(y = y),
#'         A = list(A.obs, 1),
#'         effects = list(c(
#'             list(Intercept = 1),
#'             inla.spde.make.index("spatial", spde$n.spde)
#'         ),
#'         covar = covar
#'         ),
#'         tag = "obs"
#'     )
#' stack.pred <-
#'     inla.stack(
#'         data = list(y = NA),
#'         A = list(A.pred),
#'         effects = list(c(
#'             list(Intercept = 1),
#'             inla.spde.make.index("spatial", mesh$n)
#'         )),
#'         tag = "pred"
#'     )
#' stack <- inla.stack(stack.obs, stack.pred)
#'
#' formula <- y ~ -1 + Intercept + covar + f(spatial, model = spde)
#' result1 <- inla(formula,
#'     data = inla.stack.data(stack.obs, spde = spde),
#'     family = "gaussian",
#'     control.predictor = list(
#'         A = inla.stack.A(stack.obs),
#'         compute = TRUE
#'     )
#' )
#'
#' plot(y, result1$summary.fitted.values[inla.stack.index(stack.obs, "obs")$data, "mean"],
#'     main = "Observations vs posterior predicted values at the data locations"
#' )
#'
#' result2 <- inla(formula,
#'     data = inla.stack.data(stack, spde = spde),
#'     family = "gaussian",
#'     control.predictor = list(
#'         A = inla.stack.A(stack),
#'         compute = TRUE
#'     )
#' )
#'
#' field.pred <- inla.mesh.project(
#'     proj.pred,
#'     result2$summary.fitted.values[inla.stack.index(stack, "pred")$data, "mean"]
#' )
#' field.pred.sd <- inla.mesh.project(
#'     proj.pred,
#'     result2$summary.fitted.values[inla.stack.index(stack, "pred")$data, "sd"]
#' )
#'
#' plot(field, field.pred, main = "True vs predicted field")
#' abline(0, 1)
#' image(inla.mesh.project(mesh,
#'     field = field,
#'     dims = c(200, 200)
#' ),
#' main = "True field"
#' )
#' image(inla.mesh.project(mesh,
#'     field = field.pred,
#'     dims = c(200, 200)
#' ),
#' main = "Posterior field mean"
#' )
#' image(inla.mesh.project(mesh,
#'     field = field.pred.sd,
#'     dims = c(200, 200)
#' ),
#' main = "Prediction standard deviation"
#' )
#' plot(field, (field.pred - field) / 1,
#'     main = "True field vs standardised prediction residuals"
#' )
#' @export inla.stack
inla.stack <- function(..., compress = TRUE, remove.unused = TRUE, multi.family = FALSE) {
    if (all(sapply(list(...), function(x) inherits(x, "inla.data.stack")))) {
        return(do.call(
            inla.stack.join,
            c(list(...),
              compress = compress,
              remove.unused = remove.unused,
              multi.family = multi.family
            )
        ))
    } else {
        return(do.call(
            inla.stack.sum,
            c(list(...),
              compress = compress,
              remove.unused = remove.unused
            )
        ))
    }
}


#' @describeIn inla.stack Create data stack as a sum of predictors
#' @param responses A list of response vectors, matrices, data.frame, or other special
#' response objects, such as [inla.mdata()] and [inla.surv()]. Each list element corresponds to 
#' one response family. In ordinary user-side code, the list has length 1, and longer
#' lists are created by joining stacks with `inla.stack(..., multi.family = TRUE)`.
#'
#' @export
inla.stack.sum <- function(data, A, effects, responses = NULL,
                           tag = "",
                           compress = TRUE,
                           remove.unused = TRUE) {
    
    if (!is.null(responses)) {
        responses <- ensure.canonical.response.objects(responses)
    }
    
    input.nrow <- function(x) {
        return(inla.ifelse(
            is.matrix(x) || is(x, "Matrix"),
            nrow(x),
            inla.ifelse(
                inherits(x, "inla.mdata"),
                stop("inla.mdata objects must be given as 'response' input"),
                inla.ifelse(
                    inherits(x, "inla.surv"),
                    stop("inla.surv objects must be given as 'response' input"),
                    inla.ifelse(
                        is.data.frame(x),
                        rep(nrow(x), ncol(x)),
                        length(x)
                    )
                )
            )
        ))
    }
    input.ncol <- function(x) {
        return(inla.ifelse(
            is.matrix(x) || is(x, "Matrix"),
            ncol(x),
            inla.ifelse(
                inherits(x, "inla.mdata"),
                rep(1L, length(x)),
                inla.ifelse(
                    inherits(x, "inla.surv"),
                    rep(1L, length(x)),
                    inla.ifelse(
                        is.data.frame(x),
                        rep(1L, ncol(x)),
                        1L
                    )
                )
            )
        ))
    }
    
    input.list.nrow <- function(l) {
        if (is.data.frame(l)) {
            return(input.nrow(l))
        }
        return(do.call(c, lapply(l, input.nrow)))
    }
    input.list.ncol <- function(l) {
        if (is.data.frame(l)) {
            return(input.ncol(l))
        }
        return(do.call(c, lapply(l, input.ncol)))
    }
    input.list.names <- function(l) {
        if (is.data.frame(l)) {
            return(colnames(l))
        }
        is.df <- vapply(l, is.data.frame, logical(1))
        name <- vector("list", length(l))
        if (!is.null(names(l))) {
            name[!is.df] <-
                lapply(
                    names(l)[!is.df],
                    function(x) list(x)
                )
        } else {
            name[!is.df] <- ""
        }
        name[is.df] <-
            lapply(
                l[is.df],
                function(x) as.list(colnames(x))
            )
        
        return(do.call(c, name))
    }
    
    
    parse.input.list <- function(l, n.A, error.tag, tag = "", n.A.strict = FALSE) {
        ncol <- input.list.ncol(l)
        nrow <- input.list.nrow(l)
        names <- input.list.names(l)
        if ((n.A > 1) && any(nrow == 1)) {
            for (k in which(nrow == 1)) {
                if (ncol[k] == 1) {
                    l[[k]] <- rep(l[[k]], n.A)
                    nrow[k] <- n.A
                } else {
                    stop(paste(error.tag,
                               "Automatic expansion only available for scalars.",
                               sep = ""
                    ))
                }
            }
        }
        
        if (length(unique(c(names, ""))) < length(c(names, ""))) {
            stop(paste(error.tag,
                       "All variables must have unique names\n",
                       "Names: ('",
                       paste(names, collapse = "', '", sep = ""),
                       "')",
                       sep = ""
            ))
        }
        
        for (k in seq_len(length(names))) {
            if (ncol[k] == 1) {
                names(names)[k] <- names[[k]][[1]]
                names[[k]] <- c(names[[k]][[1]])
            } else {
                names(names)[k] <- names[[k]][[1]]
                names[[k]] <- paste(names[[k]][[1]], ".", seq_len(ncol[k]), sep = "")
            }
        }
        
        names(nrow) <- names(names)
        names(ncol) <- names(names)
        
        ## data = as.data.frame(do.call(cbind, l))
        data <- as.data.frame(l)
        if (!is.null(names)) {
            names(data) <- do.call(c, names)
        }
        nrow <- nrow(data)
        if ((n.A > 1 || n.A.strict) && (nrow != n.A)) {
            stop(paste(error.tag,
                       "Mismatching row sizes: ",
                       paste(nrow, collapse = ",", sep = ""),
                       ", n.A=", n.A,
                       sep = ""
            ))
        }
        
        index <- list(seq_len(nrow))
        if (!is.null(tag)) {
            names(index) <- tag
        }
        
        info <- list(data = data, nrow = nrow, ncol = ncol, names = names, index = index)
        class(info) <- "inla.data.stack.info"
        
        return(info)
    }
    
    if (is.null(tag)) {
        stop("'tag' must not be 'NULL'")
    }
    
    ## Check if only a single block was specified.
    if (!is.list(A)) {
        A <- list(A)
        effects <- list(effects)
    }
    if (length(A) != length(effects)) {
        stop(paste("length(A)=", length(A),
                   " should be equal to length(effects)=", length(effects),
                   sep = ""
        ))
    }
    
    n.effects <- length(effects)
    
    eff <- list()
    for (k in 1:n.effects) {
        if (is.data.frame(effects[[k]])) {
            eff[[k]] <-
                parse.input.list(
                    list(effects[[k]]),
                    input.ncol(A[[k]]),
                    paste("Effect block ", k, ":\n", sep = ""),
                    tag
                )
        } else {
            if (!is.list(effects[[k]])) {
                tmp <-
                    inla.ifelse(
                        is.null(names(effects)[k]),
                        "",
                        names(effects)[k]
                    )
                effects[[k]] <- list(effects[[k]])
                names(effects[[k]]) <- tmp
            }
            eff[[k]] <-
                parse.input.list(
                    effects[[k]],
                    input.ncol(A[[k]]),
                    paste("Effect block ", k, ":\n", sep = ""),
                    tag
                )
        }
    }
    
    for (k in 1:n.effects) {
        if (is.vector(A[[k]])) {
            A[[k]] <- Matrix(A[[k]], input.nrow(A[[k]]), 1)
        }
        if ((input.ncol(A[[k]]) == 1) && (eff[[k]]$nrow > 1)) {
            if (input.nrow(A[[k]]) != 1) {
                stop(paste("ncol(A) does not match nrow(effect) for block ",
                           k, ": ",
                           input.ncol(A[[k]]), " != ", eff[[k]]$nrow,
                           sep = ""
                ))
            }
            A[[k]] <- Diagonal(eff[[k]]$nrow, A[[k]][1, 1])
        } else if (input.ncol(A[[k]]) != eff[[k]]$nrow) {
            stop(paste("ncol(A) does not match nrow(effect) for block ",
                       k, ": ",
                       input.ncol(A[[k]]), " != ", eff[[k]]$nrow,
                       sep = ""
            ))
        }
    }
    if (length(unique(input.list.nrow(A))) > 1) {
        stop(paste("Row count mismatch for A: ",
                   paste(input.list.nrow(A), collapse = ",", sep = ""),
                   sep = ""
        ))
    }
    A.nrow <- nrow(A[[1]])
    A.ncol <- input.list.ncol(A)
    
    data <-
        parse.input.list(inla.ifelse(
            is.data.frame(data),
            as.list(data),
            data
        ),
        A.nrow,
        paste("Data block:\n", sep = ""),
        tag,
        n.A.strict = TRUE
        )
    
    effects <- do.call(rbind.inla.data.stack.info, eff)
    
    A.matrix <- do.call(cbind, A)
    A.nrow <- nrow(A.matrix)
    A.ncol <- ncol(A.matrix)
    
    if (length(unique(c(names(data$names), names(effects$names)))) <
        length(c(names(data$names), names(effects$names)))) {
        stop(paste("Names for data and effects must not coincide.\n",
                   "Data names:   ",
                   paste(names(data$names), collapse = ", ", sep = ""),
                   "\n",
                   "Effect names: ",
                   paste(names(effects$names), collapse = ", ", sep = ""),
                   sep = ""
        ))
    }
    
    stack <- list(
        A = A.matrix,
        data = data,
        effects = effects,
        responses = responses
    )
    class(stack) <- "inla.data.stack"
    
    if (compress) {
        return(inla.stack.compress(stack, remove.unused = remove.unused))
    } else if (remove.unused) {
        return(inla.stack.remove.unused(stack))
    } else {
        return(stack)
    }
}

#' Expand observation vectors/matrices in stacks into to a multicolumn matrix for multiple likelihoods
#' 
#'  Internal helper method for `inla.stack.join`.
#'
#' @aliases inla.stack.mexpand
#' @name inla.stack.mexpand
#' @export
#' @param ... List of stacks that contain vector observations
#'            (existing multilikelihood observation matrices are also permitted)
#' @param old.names A vector of strings with the names of the observation vector/matrix for each stack.
#'        If a single string, this is assumed for all the stacks. (default "response")
#' @param new.name The name to be used for the expanded observation matrix,
#'        possibly the same as an old name. (default "response")
#' @return a list of modified stacks with multicolumn observations
#' @author Fabian E. Bachl \email{f.e.bachl@@bath.ac.uk} and Finn Lindgren \email{finn.lindgren@@gmail.com}
#' @export
#' @keywords internal

inla.stack.mexpand <- function(...,
                               old.names = "response",
                               new.name = "response") {
    stacks <- list(...)
    if (length(old.names) == 1) {
        old.names <- rep(old.names, length(stacks))
    }
    y.cols <- unlist(lapply(
        seq_along(stacks),
        function(x, stacks, old.names) {
            LHS <- INLA::inla.stack.LHS(stacks[[x]])[[old.names[x]]]
            ifelse(is.vector(LHS), 1, NCOL(LHS))
        },
        stacks = stacks, old.names = old.names
    ))
    y.offset <- c(0, cumsum(y.cols))
    y.cols.total <- sum(y.cols)
    for (j in seq_along(stacks)) {
        LHS <- INLA::inla.stack.LHS(stacks[[j]])
        RHS <- INLA::inla.stack.RHS(stacks[[j]])
        A <- INLA::inla.stack.A(stacks[[j]])
        responses <- INLA::inla.stack.response(stacks[[j]], drop = FALSE)
        # Access the raw tag indexing information
        tags <- list(
            data = stacks[[j]]$data$index,
            effects = stacks[[j]]$effects$index
        )
        
        if (!is.null(LHS[[old.names[j]]])) {
            # Expand the observation vector/matrix into a multilikelihood observation matrix:
            y.rows <- NROW(LHS[[old.names[j]]])
            LHS[[new.name]] <-
                cbind(
                    matrix(NA, nrow = y.rows, ncol = y.offset[j]),
                    LHS[[old.names[j]]],
                    matrix(NA, nrow = y.rows, ncol = y.cols.total - y.offset[j + 1])
                )
        }
        
        # Create the modified stack, with model compression disabled to prevent modifications:
        stacks[[j]] <-
            INLA::inla.stack.sum(
                data = LHS,
                A = A,
                effects = RHS,
                compress = FALSE,
                remove.unused = FALSE,
                responses = responses
            )
        # Since the row indexing is unchanged, copy the tag index information:
        stacks[[j]]$data$index <- tags$data
        stacks[[j]]$effects$index <- tags$effects
    }
    stacks
}

#' @describeIn inla.stack Join two or more data stacks
#'
#' @export
inla.stack.join <- function(..., compress = TRUE, remove.unused = TRUE, multi.family = FALSE) {
    if (is.character(multi.family)) {
        # NA expand response variable in the data part as matrix
        S.input <- inla.stack.mexpand(...,
                                      old.names = multi.family,
                                      new.name = multi.family)
        # Nothing more to do for multi-likelihoods
        # For non-multi, need to join the matrices (see below)
        multi.family <- FALSE
    } else {
        S.input <- list(...)
    }
    
    # Join response lists as a single list
    responses <- do.call(c,
                         lapply(S.input,
                                function(x) {
                                    if (is.null(x[["responses"]])) {
                                        list(NULL)
                                    } else {
                                        x[["responses"]]
                                    }
                                }))
    if (isFALSE(multi.family)) {
        responses <- rbind.inla.stack.responses(responses)
    }
    if (isTRUE(multi.family)) {
        responses <- expand.inla.stack.responses(responses)
    }
    
    data <- do.call(
        rbind.inla.data.stack.info,
        lapply(S.input, function(x) x$data)
    )
    effects <- do.call(
        rbind.inla.data.stack.info,
        lapply(S.input, function(x) x$effects)
    )
    ## The .bdiag form of bdiag takes a list as input.
    A <- .bdiag(lapply(S.input, function(x) x$A))
    
    S.output <- list(A = A, data = data, effects = effects, responses = responses)
    class(S.output) <- "inla.data.stack"
    
    if (length(unique(c(names(data$names), names(effects$names)))) <
        length(c(names(data$names), names(effects$names)))) {
        stop(paste("Names for data and effects must not coincide.\n",
                   "Data names:   ",
                   paste(names(data$names), collapse = ", ", sep = ""),
                   "\n",
                   "Effect names: ",
                   paste(names(effects$names), collapse = ", ", sep = ""),
                   sep = ""
        ))
    }
    
    if (compress) {
        return(inla.stack.compress(S.output, remove.unused = remove.unused))
    } else if (remove.unused) {
        return(inla.stack.remove.unused(S.output))
    } else {
        return(S.output)
    }
}





#' @describeIn inla.stack Extract tagged indices
#'
#' @export
inla.stack.index <- function(stack, tag) {
    inla.require.inherits(stack, "inla.data.stack", "'stack'")
    
    if (is.null(tag)) {
        return(list(
            data = as.vector(do.call(c, stack$data$index)),
            effects = (stack$data$nrow +
                           as.vector(do.call(c, stack$effects$index)))
        ))
    } else {
        return(list(
            data = as.vector(do.call(c, stack$data$index[tag])),
            effects = (stack$data$nrow +
                           as.vector(do.call(c, stack$effects$index[tag])))
        ))
    }
}

## Extract an object list from an inla.stack internal structure
inla.stack.do.extract <- function(dat) {
    inla.require.inherits(dat, "inla.data.stack.info", "'dat'")
    
    handle.entry <- function(x) {
        if (dat$ncol[[x]] > 1) {
            return(matrix(
                do.call(
                    c,
                    dat$data[dat$names[[x]]]
                ),
                dat$nrow,
                dat$ncol[[x]]
            ))
        } else if (is.factor(dat$data[[dat$names[[x]]]])) {
            return(dat$data[[dat$names[[x]]]])
        }
        return(as.vector(dat$data[[dat$names[[x]]]]))
    }
    
    out <- lapply(names(dat$names), handle.entry)
    names(out) <- names(dat$names)
    
    return(out)
}


#' @describeIn inla.stack Extract data associated with the "left hand side" of the model
#' (e.g. the data itself, `Ntrials`, `link`, `E`)
#'
#' @export
inla.stack.LHS <- function(stack) {
    inla.require.inherits(stack, "inla.data.stack", "'stack'")
    
    return(inla.stack.do.extract(stack$data))
}

#' @describeIn inla.stack Extract data associated with the "right hand side" of the model
#' (all the covariates/predictors)
#'
#' @export
inla.stack.RHS <- function(stack) {
    inla.require.inherits(stack, "inla.data.stack", "'stack'")
    
    return(inla.stack.do.extract(stack$effects))
}

#' @describeIn inla.stack Extract data for an inla call, and optionally join with other variables
#' @param .response.name The name to assign to the response variable when
#' extracting data from the stack. Default is `NULL`, which skips the
#' response object list.
#'
#' @export
inla.stack.data <- function(stack, ..., .response.name = NULL) {
    inla.require.inherits(stack, "inla.data.stack", "'stack'")
    
    if (is.null(.response.name) || is.null(stack[["responses"]])) {
        return(c(
            inla.stack.do.extract(stack$data),
            inla.stack.do.extract(stack$effects),
            list(...)
        ))
    } else {
        resp <- list()
        resp[[.response.name]] <- inla.stack.response(stack, drop = TRUE)
        return(c(
            inla.stack.do.extract(stack$data),
            inla.stack.do.extract(stack$effects),
            list(...),
            resp
        ))
    }
}

#' @describeIn inla.stack Extract the "A matrix" for control.predictor
#'
#' @export
inla.stack.A <- function(stack) {
    inla.require.inherits(stack, "inla.data.stack", "'stack'")
    return(stack$A)
}

#' @describeIn inla.stack Extract the response variable or list of
#' response objects
#' @param drop logical indicating whether to return the contained object
#' instead of the full list, when the stack responses list has length 1.
#' Default is `TRUE`, as needed for `inla()` single family models.
#' Use `drop = FALSE` to extract the internal response storage, regardless of
#' length.
#'
#' @export
inla.stack.response <- function(stack, drop = TRUE) {
    inla.require.inherits(stack, "inla.data.stack", "'stack'")
    if (drop) {
        if (length(stack[["responses"]]) == 0) {
            return(NULL)
        }
        if (length(stack[["responses"]]) == 1) {
            return(stack[["responses"]][[1]])
        }
    }
    return(stack[["responses"]])
}

#' @describeIn inla.stack Print information about an `inla.data.stack`
#'
#' @param x An `inla.data.stack` object for printing
#' @method print inla.data.stack
#' @export
print.inla.data.stack <- function(x, ...) {
    inla.require.inherits(x, "inla.data.stack", "'stack'")
    LHS <- inla.stack.LHS(x)
    RHS <- inla.stack.RHS(x)
    A <- inla.stack.A(x)
    response <- inla.stack.response(x, drop = FALSE)
    LHS_n <- if (is.data.frame(LHS)) {
        NROW(LHS)
    } else if (length(LHS) > 0) {
        unique(vapply(LHS, NROW, 1L))
    } else {
        0
    }
    RHS_n <- if (is.data.frame(RHS)) {
        NROW(RHS)
    } else if (length(RHS) > 0) {
        unique(vapply(RHS, NROW, 1L))
    } else {
        0
    }
    cat("Data stack with\n  ",
        "data:    ", "(", paste0(names(LHS), collapse = ", "), ")",
        ", size: ", paste0(LHS_n, collapse = ", "), "\n  ",
        "effects: ", "(", paste0(names(RHS), collapse = ", "), ")",
        ", size: ", paste0(RHS_n, collapse = ", "), "\n  ",
        "A:       ", nrow(A), " times ", ncol(A), "\n",
        sep = "")
    if (!is.null(response)) {
        cat("  response: ", length(response), " response objects\n",
            sep = "")
    }
    return(invisible(x))
}



