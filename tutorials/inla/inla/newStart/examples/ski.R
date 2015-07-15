a = read.table("ski.txt", header=F, sep="\t")
year = sort(as.numeric(substr(a[,4], 1,4)), index.return=T)

as = a[year$ix, ]

b=cbind(as[,1], year$x)

yearu = unique(year$x)

m = matrix(0, ncol=2, nrow=length(yearu))

for(i in 1:length(yearu)){
	t =yearu[i] - 1960
	len = max(b[b[,2] ==yearu[i],1])
	m[i,] = c(t,len)
}

pdf("skidata.pdf", width=8, height=6)
par(mar=c(5,5,5, 0.5), cex.lab=1.3, cex.axis=1.3, cex.main=1.4)
plot(m[,1], m[,2], xlab="year", 
	ylab="Length in meters", main="World records in ski jumping, 1961 - 2011",
	pch=19, xaxt="n")
axis(1, at=seq(1, 51, by=5), label=seq(1961, 2011, by=5))
res = lm(m[,2] ~ m[,1])
abline(res, col=2)
coef(res)
legend("topleft", c("Linear model"), col=2, lty=1)
dev.off()

data = data.frame(y=m[,2],x=m[,1])
library(INLA)
formula = y~x
model = inla(formula, family="gaussian", data=data)
plot(model, pdf=TRUE, single=FALSE, prefix="ski-")
