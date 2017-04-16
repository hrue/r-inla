## Export: inla.compare.results

##!\name{inla.compare.results}
##!\alias{inla.compare.results}
##!\alias{compare.results}
##!
##!\title{Compare INLA and MCMC results}
##!
##!\description{A small utility to compare INLA and MCMC results (OBSOLETE)}
##!\usage{
##! inla.compare.results(dir.inla = NULL, dir.mcmc = NULL)
##!}
##!
##!\arguments{
##!  \item{dir.inla}{The directory with the INLA results}
##!  \item{dir.mcmc}{The directory with the MCMC results}
##!}
##!
##!\value{%%
##!  Return nothing. This is an interactive function. 
##!
##!  This function is OBSOLETE
##!}
##!%%
##!\author{Havard Rue \email{hrue@r-inla.org}}
##!\examples{
##!## See demo("Tokyo-compare")
##!}


`inla.compare.results` = function(dir.inla = NULL, dir.mcmc = NULL)
{
    ## This function is useful for comparing the results of INLA with
    ## those obtained from a MCMC run (ie running `inla -m mcmc ...')

    if (is.null(dir.inla)) {
        print("Enter directory for INLA results ")
        dir.inla = scan(nmax=1, what = character(0))
    }
    if (!file.exists(dir.inla))
         stop(paste(" *** ERROR *** No such directory ", dir.inla))
    if (file.exists(paste(dir.inla, "/", "results.files", sep="")))
        dir.inla = paste(dir.inla, "/", "results.files", sep="")

    if (is.null(dir.mcmc)) {
        print("Enter directory for MCMC (from ``inla -m mcmc'') results ")
        dir.mcmc = scan(nmax=1, what = character(0))
    }
    if (!file.exists(dir.mcmc))
         stop(paste(" *** ERROR *** No such directory ", dir.mcmc))
    if (file.exists(paste(dir.mcmc, "/", "results.files", sep="")))
        dir.mcmc = paste(dir.mcmc, "/", "results.files", sep="")

    options(device.ask.default=TRUE)
    inla.dev.new()

    dd = c()
    for(d in dir(dir.inla))
        if (file.exists(paste(dir.inla, "/", d, "/marginal-densities.dat", sep="")))
            dd = c(dd, d)

    quit = FALSE
    while(TRUE) {

        cat("\nChose directory to view\n")
        i=1
        for(a in dd) {
            cat("[", i,"] ", a, "\n", sep="")
            i = i+1
        }
        cat("\n(<0: histogram, >0: density: q:quit) : ")

        ok = FALSE
        while(!ok) {
            id = readline("")
            if (id == "q" || id == "quit") {
                ok = TRUE
                quit = TRUE
            } else {
                id = as.integer(id)
                if (!is.na(id))
                    ok = TRUE
                else
                    cat("***Input error.\n")
            }
        }

        if (quit)
            return (invisible())

        if (id < 0) {
            hist = TRUE
            id = abs(id)
        } else {
            hist = FALSE
        }

        d.inla = paste(dir.inla, "/", dd[id], sep="")
        d.mcmc = paste(dir.mcmc, "/", dd[id], sep="")

        fnm = paste(d.inla,"/","marginal-densities.dat", sep="")
        f = paste(d.inla,"/","TAG", sep="")
        TAG = inla.ifelse(file.exists(f), inla.paste(scan(paste(d.inla,"/","TAG", sep=""),
            what=as.character())), "TAG")
        if (file.exists(fnm)) {
            ##print(paste("read marginals from ", fnm))
            marg = read.table(fnm)

            n.marg = dim(marg)[1]
            n.points = dim(marg)[2]

            fnm = paste(d.mcmc,"/","trace.dat", sep="")
            if (file.exists(fnm)) {
                samples = read.table(fnm)
                di = dim(samples)
                samples = samples[-dim(samples)[1],]
                ## with ncol=1, the previous line remove dimension attributes
                samples = as.matrix(samples, nrow=di[1]-1, ncol=di[2])
                n.samples = dim(samples)[1]

                if (dim(samples)[2] == n.marg) {
                    if (n.marg <= 4)
                        par(mfrow=c(n.marg, 1))
                    else
                        par(mfrow=c(4, 3))

                    for(i in 1:n.marg) {
                        x.values = as.double(marg[i, seq(2, n.points, by=2)])
                        f.values = as.double(marg[i, seq(3, n.points, by=2)])
                        plot(x.values, f.values, type="l", lty=1, lwd=2,
                             main=paste("Marginal for ", marg[i, 1]," ", TAG, sep=""),
                             sub = fnm, frame.plot = FALSE)
                        if (hist)
                            hist(samples[, i], n = 50, add = TRUE, probability = TRUE)
                        else
                            lines(density(samples[, i]), lty=3, lwd=2)
                    }
                }
                else
                    print("n.marg != columns(samples) ", n.marg, dim(samples)[2])
            }
            else
                print(paste(" *** ERROR *** file does not exists: ", fnm))

            rm(marg)
            rm(samples)
        }
        else
            print(paste(" *** ERROR *** file does not exits: ", fnm))
    }
}
