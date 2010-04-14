#R code to reproduce FIgure 7 in the inla manual

smooth=read.table("results-0/smoother/summary.dat")
deriv=smooth[513:1024,]

xx=deriv[,2]

split.screen( rbind(c(0,1,0.3,1), c(0,1,0,0.3)))

screen(1)
par( mar=c(2,2,2,2), oma=c(3,3,2,3) )
plot(xx[,2],type="l",ylab="",xlab="",xaxs="i")

screen(2)
par( mar=c(2,2,2,2), oma=c(3,3,2,3) )
image(1:512,1,mm,axes=F,col=gray(seq(0,1,len=3)))
