library(shiny)
library(INLA)

shinyServer(function(input, output) {
  myData <- reactive({
    inFile <- input$data
    if (is.null(inFile)) return(NULL)
    data <- read.csv(inFile$datapath,header=T)
    return(data)
  })
  
  output$datas <- renderUI({
    if(input$trang_example==1){
      HTML(paste("data(SPDEtoy)", 
                 "coords <- as.matrix(SPDEtoy[,1:2])", 
                 "p5 <- coords[1:5,]", 
                 sep = '<br/>'))
    }
    if(input$trang_example==2){
      HTML(paste("data(SPDEtoy)", 
                 "coords <- as.matrix(SPDEtoy[,1:2])", 
                 "pl.dom <- cbind(c(0,1,1,0.7,0), c(0,0,0.7,1,1))",
                 sep = '<br/>'))
    }
  })
  
  output$funcs1 <- renderText({
    if(input$trang_example==1){
      HTML(paste("inla.mesh.2d(loc = p5, max.edge=c(",input$trang_maxedge1, input$trang_maxedge2,"),
                 offset = NULL, n = NULL,", 
                 sep = '<br/>'))
    }
  })
})




HTML('&nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp function (loc = NULL, loc.domain = NULL, offset = NULL, n = NULL,<br>
    &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp boundary = NULL, interior = NULL, max.edge, min.angle = NULL,<br>
     &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp cutoff = 1e-12, plot.delay = NULL)</div>')