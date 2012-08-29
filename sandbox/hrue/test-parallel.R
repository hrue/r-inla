library(multicore)

fifo.file.1 = tempfile(tmpdir = "~/tmp")
fifo.file.2 = tempfile(tmpdir = "~/tmp")

do.a = function()
{
    cat("A\n")
    while (TRUE) {
        if (file.exists(fifo.file.1)) {
            fd1 = fifo(fifo.file.1, "rb", blocking=TRUE)
            break
        } 
    }
    fd2 = fifo(fifo.file.2, "wb", blocking=TRUE)
    cat("AA\n")
    while(TRUE) {
        i = readBin(fd1, what = character())
        cat("\na: get", i, "\n")
        i = paste(i, "a", sep="")
        if (nchar(i) > 50)
            break
        cat("\na: put", i, "\n")
        writeBin(i, fd2)
    }
}

do.b = function()
{
    cat("B\n")
    fd1 = fifo(fifo.file.1, "wb", blocking=TRUE)
    while (TRUE) {
        if (file.exists(fifo.file.2)) {
            fd2 = fifo(fifo.file.2, "rb", blocking=TRUE)
            break
        }
    }
    i = ""
    while(TRUE) {
        i = paste(i, "b", sep="")
        cat("\nb: put", i, "\n")
        writeBin(i, fd1)
        i = readBin(fd2, what = character())
        cat("\nb: get", i, "\n")
        if (nchar(i) > 50)
            break
    }
}

a = parallel(do.a())
b = parallel(do.b())
ab = collect(list(a, b))
