## Export: inla.rerun

##! \name{inla.rerun}
##! \alias{inla.rerun}
##! \alias{rerun}
##!
##! \title{Rerun an analysis}
##! \description{Rerun \code{\link{inla}} on an
##!              inla-object (output from \code{link{inla}})} 
##! \usage{
##! inla.rerun(object, plain=FALSE)
##! }
##! \arguments{
##!
`inla.rerun` <- function(
    ##!\item{object}{An \code{inla}-object, ie the output from an \code{inla}-call} 
    object, 
    ##!\item{plain}{Logical. If \code{FALSE} (default), then make
    ##!changes in \code{object} to improve the performance}
    plain = FALSE) 
    ##!}
    ##!\value{This function will take the result in \code{object},
    ##!       and rerun \code{inla} again.
    ##!       If \code{plain}  is \code{FALSE},  start  the optimization
    ##!       from the mode in \code{object} so that
    ##!       we can obtain an  improvement the mode for the hyperparameters.
    ##!       Otherwise,  start from the same configuration
    ##!       as for \code{object}.
    ##!       The returned value is an \code{inla}-object.}
    ##!\seealso{\code{\link{inla}}}
    ##!\examples{
    ##!r = inla(y ~ 1,  data = data.frame(y=1:10))
    ##!r = inla.rerun(r)
    ##!}
{
    stopifnot(any(inherits(object, "inla")))

    ## need to do this, as if its true it will regenerate the linear combinations
    object$.args$control.fixed$correlation.matrix = FALSE
    
    if (!plain) {
        object$.args$control.mode$result = NULL
        object$.args$control.mode$restart = TRUE
        object$.args$control.mode$theta = object$mode$theta
        object$.args$control.mode$x = object$mode$x
        ## do not want to change this one
        ##object$.args$control.mode$fixed = FALSE
    }
    
    new.obj = do.call("inla",  args = object$.args)
    ## revert this one, as variables are expanded and its long...
    new.obj$call = object$call

    return (new.obj)
}
