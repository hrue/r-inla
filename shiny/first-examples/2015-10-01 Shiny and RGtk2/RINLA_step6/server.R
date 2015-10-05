library(shiny)
library(INLA)


shinyServer(function(input, output) {
  
  result = reactive({
    if(input$latent==1){
      if(input$likelihood==1){
        set.seed(1234)
        n = 1000
        x = rnorm(n,  sd = 0.2)
        eta = 1 + x
        E = runif(n, min = 0, max=10)
        mu = E * exp(eta)
        size = 3
        y = rnbinom(n, size=size, mu=mu)
        data = data.frame(y, x, E)
        r = inla(y ~ 1 + x, data = data,family = "nbinomial", E=E)
        return(list(y = y,r = r,ind=1))
      }
      if(input$likelihood==2){
        set.seed(1234)
        n=100
        a=1
        b=1
        z = rnorm(n)
        eta = a + b*z
        tau = 100
        y = rnorm(n, mean = eta, sd = 1/sqrt(tau))
        data = list(y=y, z=z)
        formula = y ~ 1+z
        r = inla(formula, family = "gaussian", data = data)
        return(list(y = y,r = r,ind=2))
      }
    }
  })
  
  
  output$dataPlot <- renderPlot({
    res = result()
    plot(res$y,ylab="y")
  })
  
  
  output$inlafuncs <- renderUI({
    if(isolate(input$latent==1)){
      if(input$likelihood==1){
        return(HTML(paste("r = inla(y ~ 1 + x, data = data.frame(y, x, E),", 
                          "&nbsp &nbsp &nbsp &nbsp &nbsp &nbsp family = \"nbinomial\", E=E)", 
                          sep = '<br/>')))
      }
      if(isolate(input$likelihood==2)){
        return(HTML(paste("r = inla(y ~ 1 + z, data = list(y=y, z=z),", 
                          "&nbsp &nbsp &nbsp &nbsp &nbsp &nbsp family = \"gaussian\", scale=scale, keep=TRUE)", 
                          sep = '<br/>')))
      }
    }
  })
  

  output$inlares <- renderPrint({
    res = result()
    summary(res$r)
  })
  
})



