## Export: inla.qstat inla.qdel inla.qget inla.qnuke inla.qlog
## Export: print!inla.q summary!inla.q

##!\name{inla.qstat}
##!\alias{inla.qstat}
##!\alias{inla.q}
##!\alias{inla.qget}
##!\alias{inla.qdel}
##!\alias{inla.qnuke}
##!\alias{inla.qlog}
##!\alias{summary.inla.q}
##!\alias{print.inla.q}
##!\title{Control and view a remote inla-queue}
##!\description{
##!Control and view a remote inla-queue of submitted jobs
##!}
##!\usage{
##!inla.qget(id, remove = TRUE)
##!inla.qdel(id)
##!inla.qstat(id)
##!inla.qlog(id)
##!inla.qnuke()
##!\method{summary}{inla.q}(object,...)
##!\method{print}{inla.q}(x,...)
##!}
##!\arguments{
##!  \item{id}{The job-id which is the output from \code{inla} when the job is submitted,  the
##!            job-number or job-name. For \code{inla.qstat}, \code{id} is optional and if omitted
##!            all the jobs will be listed.}
##!  \item{remove}{Logical If FALSE, leave the job on the server after retrival,
##!                otherwise remove it (default).}
##!  \item{x}{An \code{inla.q}-object which is the output from \code{inla.qstat}}
##!  \item{object}{An \code{inla.q}-object  which is the output from \code{inla.qstat}}
##!  \item{...}{ other arguments.}
##!}
##!\details{
##!\code{inla.qstat} show job(s) on the server,
##!\code{inla.qget} fetch the results (and by default remove
##!the files on the server),  \code{inla.qdel} removes 
##!a job on the server and \code{inla.qnuke} remove all jobs on the server.
##!\code{inla.qlog} fetches the logfile only.
##!
##!The recommended procedure is to use \code{r=inla(...,
##!inla.call="submit")} and then do \code{r=inla.qget(r)} at a later
##!stage.  If the job is not finished, then \code{r} will not be
##!overwritten and this step can be repeated.  The reason for this
##!procedure, is that some information usually stored in the result
##!object does not go through the remote server, hence have to be
##!appended to the results that are retrieved from the server. Hence
##!doing \code{r=inla(..., inla.call="submit")} and then later retrive
##!it using \code{r=inla.qget(1)}, say, then \code{r} does not contain
##!all the usual information.  All the main results are there, but
##!administrative information which is required to call
##!\code{inla.hyperpar} or \code{inla.rerun} are not there.
##!}
##!\value{
##!  \code{inla.qstat} returns an \code{inla.q}-object with information about current jobs.
##!}
##!\author{Havard Rue}
##!\seealso{ \code{\link{inla}} }
##!\examples{
##!\dontrun{
##!r = inla(y~1, data = data.frame(y=rnorm(10)), inla.call="submit")
##!inla.qstat()
##!r = inla.qget(r, remove=FALSE)
##!inla.qdel(1)
##!inla.qnuke()
##!}
##!}

`summary.inla.q` = function(object, ...)
{
    print(object, ...)
}

`print.inla.q` = function(x, ...)
{
    for(k in seq_along(x)) {
        cat("Job:", x[[k]]$no, "  Id:", x[[k]]$id, "  Size:", x[[k]]$size, "  Status:", x[[k]]$status, "\n")
    }
    return (invisible(x))
}

`inla.qget` = function(id, remove = TRUE)
{
    return (inla.q(cmd = "get", id = id, remove = remove))
}

`inla.qdel` = function(id)
{
    return (inla.q(cmd = "del", id = id))
} 

`inla.qstat` = function(id)
{
    return (inla.q(cmd = "stat", id = id))
}

`inla.qlog` = function(id)
{
    return (inla.q(cmd = "log", id = id))
}

`inla.qnuke` = function()
{
    return (inla.q(cmd = "nuke"))
}

`inla.q` = function(cmd = c("get", "del", "stat", "log", "nuke"), id, remove = TRUE)
{
    cmd = match.arg(cmd)

    ## do a quick return here if possible
    if (cmd == "log" && missing(id))
        return (NULL)

    ## define some environment variables for remote computing
    vars = c(INLA_PATH = system.file("bin", package="INLA"), 
             INLA_OS = inla.os.type(), 
             INLA_HGVERSION = inla.version("hgid"))
    if (inla.os("windows")) {
        vars = c(vars,
                 INLA_SSH_AUTH_SOCK = inla.getOption("ssh.auth.sock"),
                 INLA_CYGWIN_HOME = inla.getOption("cygwin.home"), 
                 INLA_HOME = inla.cygwin.map.filename(gsub("\\\\", "/", inla.get.HOME())))
    } else {
        vars = c(vars,
                 INLA_HOME = inla.get.HOME())
        if (Sys.getenv("SSH_AUTH_SOCK") == "") {
            vars = c(vars,
                     INLA_SSH_AUTH_SOCK = inla.getOption("ssh.auth.sock"))
        }
    }
    do.call("Sys.setenv", as.list(vars))

    inla.call = system.file("bin/remote/inla.q", package="INLA")
    if (inla.os("windows")) {
        inla.call = paste(inla.call, ".cygwin", sep="")
    }

    if (!missing(id) && is.list(id)) {
        ret = id$ret
        id = id$id
    } else {
        ret = NULL
    }

    if (inla.os("linux") || inla.os("mac")) {
        output = system(paste(shQuote(inla.call), cmd, inla.ifelse(missing(id) || is.null(id), "NULL", id), as.numeric(remove)), intern=TRUE)
    } else {
        bat.file = paste(tempfile(), ".BAT",  sep="")
        cat("@echo off\n",  file=bat.file, append=FALSE)
        cat(paste(shQuote(inla.call), cmd, inla.ifelse(missing(id) || is.null(id), "NULL", id), as.numeric(remove)))
        output = try(shell(paste("@", shQuote(bat.file), sep=""), wait=TRUE, intern=TRUE))
        unlink(bat.file)
    }

    if (length(grep("^ERROR", output)) > 0) {
        output = gsub("^ERROR *", "", output)
        stop(output)
    } else if (length(grep("^DELETE", output)) > 0) {
        output = gsub("^DELETE *", "", output)
        cat("Delete job", output, "\n")
    } else if (length(grep("^NUKE", output)) > 0) {
        output = gsub("^NUKE *", "", output)
    } else if (length(grep("^LOG", output)) > 0) {
        output = gsub("^LOG *", "", output)
        if (file.exists(output)) {
            logfile = list(logfile = gsub("\t", "        ", readLines(output)))
            try(unlink(output, recursive = TRUE), silent = TRUE)
        } else {
            logfile = NULL
        }
        return (logfile)
    } else if (length(output) > 0 && length(strsplit(output, " ")[[1]]) == 1) {
        r = inla.collect.results(output, file.log = paste(output, "/results.files/Logfile.txt", sep=""))
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
                        names(a) = c("id", "no", "pid", "status", "size")
                        return(as.list(a))
                    })
        } else {
            output = list()
        }
        class(output) = "inla.q"
        return (output)
    }
}
