

prior = list(theta1 = list(normal=c(1,2)),
             theta2 = list(loggamma=c(3,4)))

read.prior = function(prior) {

    stopifnot(is.list(prior))
    print(paste("number of priors", length(prior)))

    is.listlist = is.list(prior[[1]])

    stopifnot(is.listlist)
    
    for(k in 1:length(prior)) {
        print(paste("prior", k))
        print(paste("prior for", names(prior)[k]))
        print(paste("prior", names(prior[[k]])))
        print(paste(c("param", prior[[k]][[1]])))
    }
}

