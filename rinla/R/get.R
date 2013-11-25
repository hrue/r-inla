## Export: inla.get print|inla.get summary|inla.get

##!\name{inla.get}
##!\alias{inla.get}
##!\alias{summary.inla.get}
##!\alias{print.inla.get}
##!\title{Get and view submitted jobs}
##!\description{
##!Functions to get and view submitted jobs
##!}
##!\usage{
##!inla.get(id = NULL, remove = TRUE)
##!\method{summary}{inla.get}(object,...)
##!\method{print}{summary.inla.get}(x,...)
##!}
##!\arguments{
##!  \item{id}{The job-id which is the output from \code{inla} when the job is submitted,  the
##!            job-number or job-name}
##!  \item{remove}{Logical If FALSE, leave the job on the server after getting-it,  otherwise remove it.}
##!  \item{x}{An \code{inla.get}-object,  which is the output from \code{inla.get()}}
##!  \item{object}{An \code{inla.get}-object,  which is the output from \code{inla.get()}}
##!  \item{...}{ other arguments.}
##!}
##!\details{
##!\code{inla.get} returns a list of all current jobs if \code{id=NULL},  otherwise it will get the results
##!from the job with job-id/number/name \code{id}.
##!}
##!\value{
##!  \code{inla.get()} returns an \code{inla.get}-object with information about current jobs.
##!}
##!\author{Havard Rue}
##!\seealso{ \code{\link{inla}} }
##!\examples{
##!\dontrun{
##!r = inla(y~1, data = data.frame(y=rnorm(10)), inla.call="submit")
##!inla.get()
##!r = inla.get(r)
##!}
##!}

`inla.get` = function(id = NULL, remove = TRUE)
{
    ## define some environment variables for remote computing
    inla.eval(paste("Sys.setenv(", "\"INLA_PATH\"", "=\"", system.file("bin", package="INLA"), "\"", ")", sep=""))
    inla.eval(paste("Sys.setenv(", "\"INLA_OS\"", "=\"", inla.os.type() , "\"", ")", sep=""))
    inla.eval(paste("Sys.setenv(", "\"INLA_HGVERSION\"", "=\"", inla.version("hgid") , "\"", ")", sep=""))
    if (inla.os("windows")) {
        inla.eval(paste("Sys.setenv(", "\"INLA_SSH_AUTH_SOCK\"", "=\"", inla.getOption("ssh.auth.sock"), "\"", ")", sep=""))
        inla.eval(paste("Sys.setenv(", "\"INLA_CYGWIN_HOME\"", "=\"", inla.getOption("cygwin.home"), "\"", ")", sep=""))
        inla.eval(paste("Sys.setenv(", "\"INLA_HOME\"", "=\"",
                        inla.cygwin.map.filename(gsub("\\\\", "/", inla.get.HOME())), "\"", ")", sep=""))
    } else {
        inla.eval(paste("Sys.setenv(", "\"INLA_HOME\"", "=\"", inla.get.HOME(), "\"", ")", sep=""))
        ## if SSH_AUTH_SOCK is not set, then we can pass it to the remote computing script
        if (Sys.getenv("SSH_AUTH_SOCK") == "") {
            inla.eval(paste("Sys.setenv(", "\"INLA_SSH_AUTH_SOCK\"", "=\"", inla.getOption("ssh.auth.sock"), "\"", ")", sep=""))
        }
    }

    inla.call = system.file("bin/remote/inla.get", package="INLA")
    if (inla.os("windows")) {
        inla.call = paste(inla.call, ".cygwin", sep="")
    }

    if (is.list(id)) {
        ret = id$ret
        id = id$id
    } else {
        ret = NULL
    }

    if (inla.os("linux") || inla.os("mac")) {
        output = system(paste(shQuote(inla.call), inla.ifelse(is.null(id), "NULL", id), as.numeric(remove)), intern=TRUE)
    } else {
        bat.file = paste(tempfile(), ".BAT",  sep="")
        cat("@echo off\n",  file=bat.file, append=FALSE)
        cat(paste(shQuote(inla.call), inla.ifelse(is.null(id), "NULL", id), as.numeric(remove)))
        output = try(shell(paste("@", shQuote(bat.file), sep=""), wait=TRUE), intern=TRUE)
        unlink(bat.file)
    }

    if (length(grep("^ERROR", output)) > 0) {
        stop(output)
    } else if (length(output) > 0 && length(strsplit(output, " ")[[1]]) == 1) {
        r = inla.collect.results(output)
        rr = c(r, ret)
        class(rr) = class(r)
        if (!is.null(ret) && ret$.args$keep == TRUE && !is.null(ret$.args$working.directory)) {
            ## copy the files to the target directory
            d.fnm = paste(ret$.args$working.directory, "/results.files", sep="")
            inla.dir.create(d.fnm)
            files.to.copy = paste(output, "/", dir(output, recursive=TRUE), sep="")
            file.copy(files.to.copy, d.fnm, recursive=TRUE)
        } else {
            try(unlink(output, recursive = TRUE), silent = TRUE)
        }
        if (!is.null(ret) && ret$.args$keep == FALSE && is.null(ret$.args$working.directory)) {
            try(unlink(ret$misc$inla.dir, recursive=TRUE), silent = TRUE)
        }
        return (rr)
    } else {
        if (length(output) >= 1 && nchar(output[1]) > 0) {
            output = lapply(
                    strsplit(output, " +"),
                    function(a) {
                        names(a) = c("id", "pid", "status")
                        return(as.list(a))
                    })
        } else {
            output = list()
        }
        class(output) = "inla.get"
        return (output)
    }
}

`summary.inla.get` = function(object, ...)
{
    print(object, ...)
}

`print.inla.get` = function(x, ...)
{
    if (length(x) == 0) {
        cat("No jobs available\n")
    } else {
        if (length(x) == 1) {
            cat("You have", length(x), "job\n")
        } else {
            cat("You have", length(x), "jobs\n")
        }
        for(k in seq_along(x)) {
            cat("\t Job:", k, "\tName:", x[[k]]$id, "\tStatus:", x[[k]]$status, "\n")
        }
    }
}
