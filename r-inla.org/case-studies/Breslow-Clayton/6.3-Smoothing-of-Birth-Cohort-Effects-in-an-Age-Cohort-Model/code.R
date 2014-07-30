library (INLA)

## Smoothing of Birth Cohort Effects - Iceland

## dataset downloaded from the winbugs example web page

iceland = list(age = c(1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3,
3, 3, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7, 7, 7,
7, 7, 7, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10,
11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13),
year = c(6, 7, 8, 9, 10, 11, 6, 7, 8, 9, 10, 11, 5, 6, 7, 8, 9, 10, 5,
6, 7, 8, 9, 10, 4, 5, 6, 7, 8, 9, 4, 5, 6, 7, 8, 9, 3, 4, 5, 6, 7, 8,
3, 4, 5, 6, 7, 8, 2, 3, 4, 5, 6, 7, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 5,
6, 1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 5), cases = c(2, 0, 1, 1, 1, 2, 0, 2,
1, 1, 5, 5, 1, 1, 3, 7, 12, 10, 6, 11, 9, 14, 20, 14, 7, 14, 22, 25,
29, 37, 21, 11, 29, 33, 57, 24, 15, 8, 22, 27, 38, 52, 10, 15, 22, 26,
47, 31, 8, 11, 17, 23, 31, 38, 8, 10, 24, 30, 53, 26, 5, 3, 10, 18,
22, 30, 1, 7, 11, 26, 32, 17, 5, 8, 17, 32, 31), pyr = c(41380, 43650,
49810, 58105, 57105, 76380, 39615, 42205, 48315, 56785, 55965, 33955,
29150, 38460, 40810, 47490, 55720, 55145, 27950, 37375, 39935, 46895,
54980, 27810, 25055, 27040, 36400, 39355, 46280, 54350, 24040, 26290,
35480, 38725, 45595, 25710, 22890, 23095, 25410, 34420, 37725, 44740,
21415, 21870, 24240, 33175, 36345, 21320, 17450, 19765, 20255, 22760,
31695, 34705, 15350, 17720, 18280, 20850, 29600, 15635, 9965, 12850,
15015, 15725, 18345, 26400, 8175, 11020, 13095, 14050, 16480, 10885,
7425, 10810, 12260, 14780, 13600) )

iceland=data.frame(iceland)
iceland$year2=iceland$year
iceland$year3=iceland$year

x<-cbind(1, 1:11)
extraconstraint<-list(A=matrix((1:11) - mean(1:11), 1, 11), e=matrix(0,1,1))

iceland.inla.fit = inla(cases ~ -1 + as.factor(age) + f(year,
    model="iid", param=c(0.5, 0.00149) ) + year2 + f(year3,
    values=1:11, model="rw2", param=c(0.5, 0.001), constr=T,
    extraconstr=extraconstraint),
        data=iceland, family="poisson",
        offset=I(log(pyr)),
        control.inla= list(int.strategy = "grid",dz=0.25),
        control.predictor=list(compute=TRUE) )
iceland.hyperpar = inla.hyperpar (iceland.inla.fit)
summary(iceland.inla.fit)
inla.emarginal(function(x) 1/x^.5, iceland.hyperpar$marginals[[1]])
