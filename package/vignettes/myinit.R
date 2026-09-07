inla.setOption(num.threads = "2:1:1")
if (FALSE) {
    INLA:::inla.my.update()
    inla.setOption(inla.call="inla.mkl.work")
    inla.setOption(num.threads = "4:1:2")
    ##inla.setOption(smtp = "pardiso")
    ##inla.setOption(pardiso.license="~/sys/licenses/pardiso.lic")
}
