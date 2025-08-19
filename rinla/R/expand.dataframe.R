## Nothing to export from here

`inla.expand.dataframe.1` <- function(response, dataframe,
                                      control.hazard = inla.set.control.hazard.default())
{
    n.intervals <- control.hazard$n.intervals
    cutpoints <- control.hazard$cutpoints

    if (!inherits(response, "inla.surv")) {
        stop("Response has to be an object of class `inla.surv'")
    }
    class(response) <- NULL
    event <- response$event
    nn <- length(event)

    time <- numeric(nn)
    time[response$event == 1L] <- response$time[response$event == 1L]
    time[response$event == 0L] <- response$lower[response$event == 0L]

    truncation <- response$truncation

    ## create cutpoints if not provided
    if (is.null(cutpoints)) {
        cutpoints <- seq(0L, max(time), length.out = n.intervals + 1L)
    }

    new.data <- inla.get.poisson.data.1(time = time, truncation = truncation, event = event, cutpoints = cutpoints)
    expand.df <- table(new.data$indicator)

    if (!missing(dataframe) && prod(dim(dataframe)) > 0L) {
        new.dataframe <- as.data.frame(matrix(0.0, length(new.data$y), dim(dataframe)[2L]))
        for (i in 1L:dim(dataframe)[2L]) {
            new.dataframe[, i] <- rep(dataframe[, i], expand.df)
        }

        ## alternative code,  not faster...
        ## new.dataframe = apply(dataframe, 2, rep, times = expand.df)

        new.dataframe <- as.data.frame(new.dataframe)
        ## just give some name that we are able to recognise afterwards
        names(new.dataframe) <- paste("fake.dataframe.names", 1L:dim(new.dataframe)[2L])
    } else {
        new.dataframe <- NULL
    }

    res <- data.frame(
        y..coxph = new.data$y,
        E..coxph = new.data$E,
        expand..coxph = rep(1:nrow(dataframe), expand.df),
        baseline.hazard = cutpoints[new.data$baseline.hazard],
        baseline.hazard.idx = new.data$baseline.hazard,
        baseline.hazard.time = cutpoints[new.data$baseline.hazard],
        baseline.hazard.length = diff(cutpoints)[new.data$baseline.hazard],
        new.dataframe
    )

    names(res)[grep("fake.dataframe.names", names(res))] <- names(dataframe)

    return(list(data = res, data.list = list(baseline.hazard.values = cutpoints)))
}

`inla.get.poisson.data.1` <- function(time, truncation, event, cutpoints)
{
    if (TRUE) {
        ##
        ## this is the new code of (Jan 21st, 2021) that runs better
        ##
        data.new <- numeric(0L)
        nn <- length(event)
        start <- as.numeric(cut(truncation, cutpoints, include.lowest = FALSE))
        end <- as.numeric(cut(time, cutpoints, include.lowest = TRUE))
        ds <- diff(cutpoints)

        res <- lapply(1L:length(time),
                      function(i) ({
                          if (is.na(start[i])) {
                              if (end[i] > 1.0) {
                                  dc <- cbind(ds[1L:(end[i] - 1L)],
                                              rep(0L, (end[i] - 1L)),
                                              rep(i, (end[i] - 1L)),
                                              c(1L:(end[i] - 1L))
                                              )
                              } else {
                                  dc <- numeric(0L)
                              }
                              dc <- rbind(dc, cbind(time[i] - (cutpoints[end[i]]), event[i], i, end[i]))
                          }
                          else {
                              if (start[i] < end[i]) {
                                  dc <- cbind((cutpoints[start[i] + 1L] - truncation[i]), 0L, i, start[i])
                                  if (end[i] > (start[i] + 1L)) {
                                      dc <- rbind(dc,
                                                  cbind(ds[(start[i] + 1L):(end[i] - 1L)],
                                                        rep(0L, (end[i] - start[i] - 1L)),
                                                        rep(i, (end[i] - start[i] - 1L)),
                                                        c((start[i] + 2L):(end[i]) - 1L)))
                                  }
                                  dc <- rbind(dc, cbind(time[i] - (cutpoints[end[i]]), event[i], i, end[i]))
                              } else if (start[i] == end[i]) {
                                  dc <- cbind(time[i] - (cutpoints[end[i]]), event[i], i, end[i])
                              } else {
                                  stop("Truncation cannot be greater than time")
                              }
                          }
                          return(t(dc))
                      }))
        data.new <- matrix(unlist(res), ncol = 4, byrow = TRUE)
        ##
    } else {
        ##
        ## this is the old code that is slow for huge ammount of data
        ##
        data.new <- numeric(0L)
        nn <- length(event)
        start <- as.numeric(cut(truncation, cutpoints, include.lowest = FALSE))
        end <- as.numeric(cut(time, cutpoints, include.lowest = TRUE))
        ds <- diff(cutpoints)

        for (i in 1L:length(time)) {
            if (is.na(start[i])) {
                if (end[i] > 1.0) {
                    dc <- cbind(
                        ds[1L:(end[i] - 1L)], rep(0L, (end[i] - 1L)), rep(i, (end[i] - 1L)),
                        c(1L:(end[i] - 1L))
                    )
                } else {
                    dc <- numeric(0L)
                }
                dc <- rbind(dc, cbind(time[i] - (cutpoints[end[i]]), event[i], i, end[i]))
                data.new <- rbind(data.new, dc)
            }
            else {
                if (start[i] < end[i]) {
                    dc <- cbind((cutpoints[start[i] + 1L] - truncation[i]), 0L, i, start[i])
                    if (end[i] > (start[i] + 1L)) {
                        dc <- rbind(dc, cbind(
                                            ds[(start[i] + 1L):(end[i] - 1L)], rep(0L, (end[i] - start[i] - 1L)),
                                            rep(i, (end[i] - start[i] - 1L)),
                                            c((start[i] + 2L):(end[i]) - 1L)
                                        ))
                    }
                    dc <- rbind(dc, cbind(time[i] - (cutpoints[end[i]]), event[i], i, end[i]))
                    data.new <- rbind(data.new, dc)
                } else if (start[i] == end[i]) {
                    dc <- cbind(time[i] - (cutpoints[end[i]]), event[i], i, end[i])
                    data.new <- rbind(data.new, dc)
                } else {
                    stop("Truncation cannot be greater than time")
                }
            }
        }
    }

    data.new <- data.frame(
        E = data.new[, 1L],
        y = data.new[, 2L],
        indicator = data.new[, 3L],
        baseline.hazard = data.new[, 4L]
    )

    return(data.new)
}

`inla.expand.dataframe.2` <- function (response, dataframe, control.hazard = inla.set.control.hazard.default())
{
    n.intervals <- control.hazard$n.intervals
    cutpoints <- control.hazard$cutpoints
    if (!inherits(response, "inla.surv")) {
        stop("Response has to be an object of class `inla.surv'")
    }
    class(response) <- NULL
    event <- response$event
    nn <- length(event)
    time <- numeric(nn)
    time[response$event == 1L] <- response$time[response$event ==
                                                1L]
    time[response$event == 0L] <- response$lower[response$event ==
                                                 0L]
    truncation <- response$truncation
    subject <- response$subject
    jj <- unique(response$subject)
    subject.first.line <- numeric(length(jj))
    for (i in 1L:length(jj)) {
        rows <- which(subject == i)
        sem <- dataframe[rows, ]
        if (mode(apply(sem, 2L, unique)) == "list") {
            stop("coxph with subject only works for fixed covariates")
        }
        subject.first.line[i] <- rows[1L]
    }
    dataframe.copy <- dataframe[subject.first.line, ]
    if (is.null(cutpoints)) {
        cutpoints <- seq(0, max(time), length.out = n.intervals +
                                           1L)
    }
    new.data <- inla.get.poisson.data.2(time = time, subject = subject,
                                        event = event, cutpoints = cutpoints)
    expand.df <- table(new.data$indicator)
    if (!missing(dataframe.copy) && prod(dim(dataframe.copy)) > 0L) {
        new.dataframe <- as.data.frame(matrix(0, length(new.data$y),
                                              dim(dataframe.copy)[2L]))
        for (i in 1L:dim(dataframe.copy)[2L]) {
            new.dataframe[, i] <- rep(dataframe.copy[, i], expand.df)
        }
        new.dataframe <- as.data.frame(new.dataframe)
        names(new.dataframe) <- paste("fake.dataframe.names",
                                      1L:dim(new.dataframe)[2L])
    }
    else {
        new.dataframe <- NULL
    }
    res <- data.frame(y..coxph = new.data$y, E..coxph = new.data$E,
                      expand..coxph = rep(1:nrow(dataframe.copy), expand.df), baseline.hazard = cutpoints[new.data$baseline.hazard],
                      baseline.hazard.idx = new.data$baseline.hazard, baseline.hazard.time = cutpoints[new.data$baseline.hazard],
                      baseline.hazard.length = diff(cutpoints)[new.data$baseline.hazard],
                      subject = new.data$indicator, new.dataframe)
    names(res)[grep("fake.dataframe.names", names(res))] <- names(dataframe)
    return(list(data = res, data.list = list(baseline.hazard.values = cutpoints)))
}

`inla.get.poisson.data.2` <- function (subject, time, event, cutpoints)
{
    data.new <- numeric(0L)
    nn <- max(subject)
    ds <- diff(cutpoints)
    ris <- matrix(0, max(subject), (length(cutpoints) - 1L))
    dataF <- cbind(subject, time, event)
    end <- integer(0L)
    for (i in 1L:nn) {
        da <- matrix(dataF[dataF[, 1L] == i, ], ncol = 3L)
        rec <- cut(da[, 2L], cutpoints, labels = 1L:(length(cutpoints) -
                                                     1L))
        ris[i, ] <- tapply(da[, 3L], rec, sum)
        ris[i, ][is.na(ris[i, ])] <- 0
        end <- c(end, as.numeric(cut(max(time[subject == i]),
                                     cutpoints, include.lowest = TRUE)))
    }
    totalevent <- 0
    length(totalevent) <- 0
    for (i in 1L:nn) {
        totalevent <- c(totalevent, ris[i, ])
    }
    E <- numeric(0)
    length(E) <- 0
    for (i in 1L:nn) {
        if (end[i] == 1L) {
            E <- c(E, c(max(time[subject == i]) - cutpoints[1L]),
                   rep(0, length(cutpoints) - 1L - end[i]))
        }
        else {
            E <- c(E, ds[1L:end[i] - 1L], max(time[subject ==
                                                   i]) - cutpoints[end[i]], rep(0, length(cutpoints) -
                                                                                   1L - end[i]))
        }
    }
    index <- rep(1L:nn, each = length(cutpoints) - 1L)
    baseline.hazard <- rep(1L:(length(cutpoints) - 1L), nn)
    dc <- cbind(index, E, totalevent, baseline.hazard)
    data.new <- data.frame(indicator = dc[, 1L], E = dc[, 2L],
                           y = dc[, 3L], baseline.hazard = dc[, 4L])
    return(data.new)
}

