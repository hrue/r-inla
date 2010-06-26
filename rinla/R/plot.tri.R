`plot.inla.trimesh` = function(TV, S, color = NULL, ...)
{
    ## Make indices 1 based... should make this safer
    if (min(TV) == 0) {
        TV <- TV + 1
    }


    tTV = t(TV);
    tETV = t(TV[,c(1,2,3,1,NA)]);
    Tx = S[tTV,1]
    Ty = S[tTV,2]
    Tz = S[tTV,3]
    Tcol = color[tTV]
    Ex = S[tETV,1]
    Ey = S[tETV,2]
    Ez = S[tETV,3]
    Ecol = rgb(0.3,0.3,0.3)
    points3d(S, color="black", ...)
    ## need to know if argument "lwd=.." is in the dots or not. In
    ## case it is, do not add lwd=1 argument.
    if (length(grep("lwd", names(match.call(expand.dots=TRUE)))) > 0) {
        lines3d(Ex, Ey, Ez, color=Ecol, ...)
    } else {
        lines3d(Ex, Ey, Ez, color=Ecol, lwd=1, ...)
    }
    ## same issue here...
    if (length(grep("specular", names(match.call(expand.dots=TRUE)))) > 0) {
        triangles3d(Tx, Ty, Tz, color=Tcol, ...)
    } else {
        triangles3d(Tx, Ty, Tz, color=Tcol, specular="black", ...)
    }
}

## library(geometry)
## S = cbind(x=rnorm(30), y=rnorm(30), z=0)
## TV = delaunayn(S[,1:2]) # NOTE: inconsistent triangle orders, only for test.
## trimesh(TV,S)
##
## colors = rgb(runif(30),runif(30),runif(30))
## rgl.viewpoint(0,0,fov=20)
## plot.inla.trimesh(TV,S,colors)

##    Ecol = col2rgb(color)/256
##    Ecol = Ecol*0.5+(1-0.5)*0 # Rescale towards black
##    Ecol = 1-Ecol # Invert
##    Ecol = Ecol[,c(2,3,1)] # Permute
##    Ecol = rgb(Ecol[1,],Ecol[2,],Ecol[3,],maxColorValue = 1)
##    Ecol = Ecol[tETV]

