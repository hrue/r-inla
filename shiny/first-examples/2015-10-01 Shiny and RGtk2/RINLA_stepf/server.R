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
  
  
  
  output$datas <- renderUI({
    if(input$trang_example==1){
      return(HTML(paste("data(SPDEtoy)", 
                 "coords <- as.matrix(SPDEtoy[,1:2])", 
                 "p5 <- coords[1:5,]", 
                 sep = '<br/>')))
    }
    if(input$trang_example==2){
      return(HTML(paste("data(SPDEtoy)", 
                 "coords <- as.matrix(SPDEtoy[,1:2])", 
                 "pl.dom <- cbind(c(0,1,1,0.7,0), c(0,0,0.7,1,1))",
                 sep = '<br/>')))
    }
  })
  
  data(SPDEtoy)
  coords <- as.matrix(SPDEtoy[,1:2])
  pl.dom <- cbind(c(0,1,1,0.7,0), c(0,0,0.7,1,1))
  p5 <- coords[1:5,]
  
  output$funcs1 <- renderText({
    if(input$trang_example==1){
      return(paste("inla.mesh.2d(loc = p5, max.edge=c(",input$trang_maxedge1,",", input$trang_maxedge2,"),
                 offset = c(",input$trang_offset1,",",input$trang_offset2,"), cutoff = ",input$cutoff,")", 
                 sep = ''))
    }
    if(input$trang_example==2){
      return(paste("inla.mesh.2d(loc.domin = pl.dom, max.edge=c(",input$trang_maxedge1,",", input$trang_maxedge2,"),
                 offset = c(",input$trang_offset1,",",input$trang_offset2,"), cutoff = ",input$cutoff,")", 
                   sep = ''))
    }
  })
  
  output$spdeplot <- renderPlot({
    if(input$trang_example==1){
      mi = inla.mesh.2d(loc = p5, max.edge=c(input$trang_maxedge1,input$trang_maxedge2),
                        offset = c(input$trang_offset1,input$trang_offset2), cutoff = input$cutoff)
      plot(pl.dom, type='l', col=3, lwd=2, xlim=c(-0.57,1.57), asp=1, axes=FALSE)
      plot(mi, add=TRUE)
      points(p5, pch=19, col=2)
    }
    if(input$trang_example==2){
      mi = inla.mesh.2d(loc.domain = pl.dom, max.edge=c(input$trang_maxedge1,input$trang_maxedge2),
                        offset = c(input$trang_offset1,input$trang_offset2), cutoff = input$cutoff)
      plot(pl.dom, type='l', col=3, lwd=2, xlim=c(-0.57,1.57), asp=1, axes=FALSE)
      plot(mi, add=TRUE)
      points(p5, pch=19, col=2)
    }
    
  })
  
})



