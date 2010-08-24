`plot.inla.trimesh` = function(TV, S, color = NULL, lwd = 1, specular = "black", ...)
{
    ## Make indices 1 based.  Deprecated and will be deactivated
    if (min(TV) == 0) {
      warning("Zero-based indices in TV are deprecated and will be deactivated in a future version.")
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
    lines3d(Ex, Ey, Ez, color=Ecol, lwd=lwd, ...)
    triangles3d(Tx, Ty, Tz, color=Tcol, specular=specular, ...)

    return (invisible())
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

