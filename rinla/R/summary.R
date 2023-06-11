#' Summary for a INLA fit
#' 
#' Takes a fitted `inla` or `surv.inla` object produced by
#' `inla` or `surv.inla` and produces a summary from it.
#' 
#' Posterior mean and standard deviation (together with quantiles or cdf) are
#' printed for the fixed effects in the model.
#' 
#' For the random effects the function `summary()` prints the posterior
#' mean and standard deviations for the hyperparameters
#' 
#' If the option `short.summary` is set to `TRUE` using
#' `inla.setOption`, then a less verbose summary variant will be used,
#' which might be more suitable for Markdown documents.
#' 
#' @aliases summary.inla summary.surv.inla print.summary.inla
#' @param object a fitted `inla` object as produced by `inla`.
#' @param x a `summary.inla` object produced by `summary.inla`
#' @param digits Integer Number of digits
#' @param include.lincomb Logcial Include the summary for the the linear
#' combinations or not
#' @param ...  other arguments.
#' @return `summary.inla` returns an object of class `summary.inla`,
#' a list of components to print.
#' @author Sara Martino and Havard Rue
#' @seealso [inla()]
#' 
#' @method summary inla
#' @rdname summary
#' @export
`summary.inla` <- function(object, digits = 3L, include.lincomb = TRUE, ...) {
    `internal.summary.inla` <- function(object, digits = 3L, include.lincomb = TRUE, ...) {
        inla.eval.dots(...)
        ## provides a summary for a inla object
        ret <- list()

        maxlen <- 2048L
        if (sum(nchar(object$call)) > maxlen) {
            ret <- c(
                ret,
                list(call = paste(
                    substr(inla.paste(object$call), 1L, maxlen),
                    "<<<the rest is not shown>>>"
                ))
            )
        } else {
            ret <- c(ret, list(call = object$call))
        }

        ## might not be if using collect directly
        if (inla.is.element("cpu.used", object)) {
            ret <- c(ret, list(
                cpu.used =
                    paste0(
                        names(object$cpu.used)[1], " = ",
                        format(object$cpu.used[1], digits = digits), ", ",
                        names(object$cpu.used)[2], " = ",
                        format(object$cpu.used[2], digits = digits), ", ",
                        names(object$cpu.used)[3], " = ",
                        format(object$cpu.used[3], digits = digits), ", ",
                        names(object$cpu.used)[4], " = ",
                        format(object$cpu.used[4], digits = digits)
                    )
            ))
        } else {
            ret <- c(ret, list(cpu.used = NA))
        }

        if (!is.null(object$summary.fixed) && length(object$summary.fixed) > 0) {
              ret <- c(ret, list(fixed = round(as.matrix(object$summary.fixed), digits = digits)))
          }

        if (include.lincomb) {
            if (!is.null(object$summary.lincomb) && any(names(object) == "summary.lincomb")
            && (length(object$summary.lincomb) > 0)) {
                  ret <- c(ret, list(lincomb = round(as.matrix(object$summary.lincomb), digits = digits)))
              }
            if (!is.null(object$summary.lincomb.derived) && length(object$summary.lincomb.derived) > 0) {
                  ret <- c(ret, list(lincomb.derived = round(as.matrix(object$summary.lincomb.derived),
                      digits = digits
                  )))
              }
        } else {
            if (!is.null(object$summary.lincomb) && any(names(object) == "summary.lincomb")
            && (length(object$summary.lincomb) > 0)) {
                m <- nrow(as.matrix(object$summary.lincomb))
                ret <- c(ret, list(lincomb = paste("<", m, " lincomb not included>", sep = "")))
            }
            if (!is.null(object$summary.lincomb.derived) && length(object$summary.lincomb.derived) > 0) {
                m <- nrow(as.matrix(object$summary.lincomb.derived))
                ret <- c(ret, list(lincomb.derived = paste("<", m, " lincomb.derived not included>", sep = "")))
            }
        }

        if (!is.null(object$summary.hyperpar) && length(object$summary.hyperpar) > 0) {
              ret <- c(ret, list(hyperpar = round(object$summary.hyperpar, digits = digits)))
          }

        if (!is.null(object$summary.random) && length(object$summary.random) > 0) {
              ret <- c(ret, list(random.names = names(object$summary.random), random.model = object$model.random))
          }

        if (!is.null(object$dic)) {
            ret <- c(ret, list(dic = lapply(object$dic, round, digits = digits)))
        }

        if (!is.null(object$waic)) {
            ret <- c(ret, list(waic = lapply(object$waic, round, digits = digits)))
        }

        if (!is.null(object$mlik)) {
              ret <- c(ret, list(mlik = round(object$mlik, digits = digits)))
          }

        if (!is.null(object$cpo$cpo) && length(object$cpo$cpo) > 0L) {
            ret <- c(ret, list(cpo = lapply(object$cpo, round, digits = digits)))
        }

        if (!is.null(object$gcpo$gcpo) && length(object$gcpo$gcpo) > 0L) {
            ret <- c(ret, list(gcpo = lapply(object$gcpo$gcpo, round, digits = digits)))
        }

        if (!is.null(object$summary.linear.predictor)) {
            ret <- c(ret, list(linear.predictor = round(object$summary.linear.predictor, digits = digits)))
        }

        ret <- c(ret, list(family = object$family))
        class(ret) <- "summary.inla"

        return(ret)
    }

    ## for G & J
    if (!is.null(object$joint.model)) {
        return(summary(object$joint.model, ...))
    }

    ret <- internal.summary.inla(object,
        digits = digits,
        include.lincomb = include.lincomb, ...
    )
    ss <- inla.getOption("short.summary")
    if (!is.null(ss) && ss) {
        new.ret <- list(
            fixed = ret$fixed, hyperpar = ret$hyperpar,
            dic = ret$dic, waic = ret$waic
        )
        class(new.ret) <- class(ret)
        ret <- new.ret
    }
    return(ret)
}


#' @rdname summary
#' @method print summary.inla
#' @export
`print.summary.inla` <- function(x, digits = 3L, ...) 
{
    form <- strwrap(inla.formula2character(x$call))
    if (!is.null(x$call)) {
        cat("\nCall:\n")
        for (i in seq_along(form)) {
            cat("  ", form[i], "\n")
        }
    }

    if (inla.is.element("cpu.used", x)) {
        cat("Time used:\n", "  ", x$cpu.used, "\n")
    }

    if (inla.is.element("fixed", x)) {
        cat("Fixed effects:\n")
        print.default(x$fixed)
        cat("\n")
    }

    if (inla.is.element("lincomb.derived", x)) {
        cat("Linear combinations (derived):\n")
        if (is.character(x$lincomb.derived)) {
            cat("\t", x$lincomb.derived)
        } else {
            print.default(x$lincomb.derived)
        }
        cat("\n")
    }

    if (inla.is.element("lincomb", x)) {
        cat("Linear combinations:\n")
        if (is.character(x$lincomb)) {
            cat("\t", x$lincomb)
        } else {
            print.default(x$lincomb)
        }
        cat("\n")
    }

    if (inla.is.element("random.names", x)) {
        cat("Random effects:\n")
        cat("  Name\t ", "Model\n ")
        for (i in 1:length(x$random.names)) {
            cat("  ", paste0(inla.nameunfix(x$random.names[i]), " ", x$random.model[i], "\n"))
        }
        cat("\n")
    }

    if (inla.is.element("hyperpar", x)) {
        cat("Model hyperparameters:\n")
        print(format(x$hyperpar, digits = digits, nsmall = 2), quote = FALSE)
        cat("\n")
    }

    if (inla.is.element("dic", x)) {
        cat(paste("Deviance Information Criterion (DIC) ...............: ",
                  format(x$dic$dic, digits = digits, nsmall = 2), "\n",
                  "Deviance Information Criterion (DIC, saturated) ....: ",
                  format(x$dic$dic.sat, digits = digits, nsmall = 2), "\n",
                  "Effective number of parameters .....................: ",
                  format(x$dic$p.eff, digits = digits, nsmall = 2), "\n\n",
                  sep = ""
                  ))
    }

    if (inla.is.element("waic", x)) {
        cat(paste("Watanabe-Akaike information criterion (WAIC) ...: ",
                  format(x$waic$waic, digits = digits, nsmall = 2), "\n",
                  "Effective number of parameters .................: ",
                  format(x$waic$p.eff, digits = digits, nsmall = 2), "\n\n",
                  sep = ""
                  ))
    }

    if (inla.is.element("mlik", x)) {
        cat(paste("Marginal log-Likelihood: ", format(x$mlik[2], digits = digits, nsmall = 2), "\n"))
    }

    msg <- ""
    if (inla.is.element("cpo", x)) {
        msg <- paste0(msg, "CPO, PIT")
    }
    if (inla.is.element("gcpo", x)) {
        msg <- paste0(msg, ", GCPO")
    }
    if (inla.is.element("po", x)) {
        msg <- paste0(msg, ", PO")
    }
    msg <- paste0(msg, " is computed")
    cat(msg, "\n")
    
    if (inla.is.element("linear.predictor", x)) {
        cat(sep = "", 
            "Posterior summaries for the linear predictor and the fitted values are computed\n",
            "(Posterior marginals needs also 'control.compute=list(return.marginals.predictor=TRUE)')\n\n")
    }
}
