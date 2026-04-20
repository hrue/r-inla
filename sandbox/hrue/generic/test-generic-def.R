library(multicore)

test.read = function(def)
{
    cat("OPEN\n")
    while (TRUE) {
        if (file.exists(def$fifo$R2c)) {
            R2c = fifo(def$fifo$R2c, "rb", blocking=TRUE)
            break
        }
        cat("test.read: sleep\n")
        Sys.sleep(1)
    }
    c2R = fifo(def$fifo$c2R, "wb", blocking=TRUE)
    cat("OPEN DONE\n")

    cat("Put command: initial\n")
    writeBin("initial", c2R)
    ntheta = readBin(R2c, what = integer(), n = 1)
    cat("Got ntheta", ntheta, "\n")
    stopifnot(ntheta >= 0L)
    if (ntheta > 0) {
        initial = readBin(R2c, what = numeric(), n = ntheta)
    } else {
        initial = numeric(0)
    }

    writeBin("graph", c2R)
    len = readBin(R2c, what = integer(), n = 1)
    cat("Got graph len", len, "\n")
    i = readBin(R2c, what = integer(), n = len)
    j = readBin(R2c, what = integer(), n = len)
    print(i)
    print(j)
    n = max(i, j) + 1L
    stopifnot(min(i) == 0 && min(j) == 0 && max(i) == n-1 && max(j) == n-1)
    G = sparseMatrix(i = i+1, j = j+1, x = rep(1, n))
    ##inla.spy(G)
    ##Sys.sleep(2)

    writeBin("Q", c2R)
    initial = rep(1, ntheta)
    writeBin(initial, c2R)
    
    len = readBin(R2c, what = integer(), n = 1)
    cat("Got Q len", len, "\n")
    i = readBin(R2c, what = integer(), n = len)
    j = readBin(R2c, what = integer(), n = len)
    x = readBin(R2c, what = numeric(), n = len)
    stopifnot(min(i) == 0 && min(j) == 0 && max(i) == n-1 && max(j) == n-1)
    Q = sparseMatrix(i = i+1L, j = j+1L, x = x)
    ##inla.display.matrix(Q)
    ##Sys.sleep(2)

    cat("put 'extra'\n")
    writeBin("extra", c2R)
    cat("put theta", initial, "\n")
    writeBin(initial, c2R)
    extra.model = readBin(R2c, what = numeric(), n=1)
    cat("Got extra.model", extra.model, "\n")
    extra.prior = readBin(R2c, what = numeric(), n=1)
    cat("Got extra.prior", extra.prior, "\n")

    cat("send exit\n")
    writeBin("exit", c2R)
}

def = inla.R.generic.define(inla.R.generic.ar1.def, n=10, ntheta=2)
a = parallel(inla.R.generic.loop(def, debug=TRUE))
b = parallel(test.read(def))
ab = collect(list(a, b))
