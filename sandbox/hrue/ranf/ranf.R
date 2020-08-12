ranf = function(n=1, F=2)
{
    if (F <= 1.0) {
        return (rep(1.0, n))
    } else {
        len = F - 1.0 / F;

        nx = rbinom(1, size=n,  prob= len / (len + 2.0 * log(F)))
        ny = n - nx

        if (nx > 0L) {
            x = 1/F + runif(nx, max = len)
        } else {
            x = NULL
        }
        if (ny > 0L) {
            y = F^runif(ny, min=-1, max=1)
        } else {
            y = NULL
        }

        return (c(x, y))
    }
}
