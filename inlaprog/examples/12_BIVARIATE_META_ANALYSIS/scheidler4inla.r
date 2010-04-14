data <- read.table("scheid.dat", header=T)
test <- rep(data$test, each=2)
lag <- as.numeric(test==1)
ct <- as.numeric(test==2)
mri <- as.numeric(test==3)

id <- 0:(2*nrow(data)-1)
y <-  c(t(cbind(data$tp, data$tn)))
num <- c(t(cbind(data$ndis, data$nndis)))

dat <- data.frame(id,num,y)

intercept_tp <- rep(c(1,0), length(dat$id)/2)
intercept_tn <- rep(c(0,1), length(dat$id)/2)

#write.table(cbind(dat$id, intercept_tp), file="scheidler_intercept_tp.dat", row.names=F, col.names=F, sep=" ")
#write.table(cbind(dat$id, intercept_tn), file="scheidler_intercept_tn.dat", row.names=F, col.names=F, sep=" ")
write.table(cbind(dat$id, lag*intercept_tp),file="covariate_lag_tp.dat",row.names=F,col.names=F,sep=" ")
write.table(cbind(dat$id,lag*intercept_tn),file="covariate_lag_tn.dat",row.names=F,col.names=F,sep=" ")
write.table(cbind(dat$id,ct*intercept_tp),file="covariate_ct_tp.dat",row.names=F,col.names=F,sep=" ")
write.table(cbind(dat$id,mri*intercept_tp),file="covariate_mr_tp.dat",row.names=F,col.names=F,sep=" ")
write.table(cbind(dat$id,ct*intercept_tn),file="covariate_ct_tn.dat",row.names=F,col.names=F,sep=" ")
write.table(cbind(dat$id,mri*intercept_tn),file="covariate_mr_tn.dat",row.names=F,col.names=F,sep=" ")


write.table(cbind(dat$id, dat$id), file="scheidler_cov.dat", row.names=F, col.names=F, sep=" ")
write.table(dat,file="scheidler.dat", row.names=F, col.names=F, sep=" ")
