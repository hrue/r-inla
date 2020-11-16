context("test 'latent bym'")

test_that("Case 1", {
    data(Germany)
    g = system.file("demodata/germany.graph", package="INLA")
    Germany = cbind(Germany,region.struct=Germany$region)
    prior.iid = c(1,0.1)
    prior.besag = c(1,0.001)
    formula1 = Y ~ 1 + f(region.struct,
            model="besag",graph=g, constr=TRUE,
            hyper = list(prec = list(param = prior.besag))) +
                f(region,
                  model="iid", hyper = list(prec = list(param = prior.iid)))
    result1  =  inla(formula1,family="poisson",data=Germany,E=E,
            control.predictor = list(compute=TRUE))
    
    formula1.bym = Y ~ 1 + f(region,
            model = "bym", graph = g, constr=TRUE, 
            hyper = list(
                    prec.unstruct = list(param = prior.iid),
                    prec.spatial = list(param = prior.besag)))
    
    result1.bym = inla(formula1.bym,family="poisson",data=Germany,E=E, 
            control.predictor = list(compute=TRUE))
    
    expect_true(mean(abs(result1$summary.linear.predictor$mean -
                         result1.bym$summary.linear.predictor$mean)) < 1e-4)
})

