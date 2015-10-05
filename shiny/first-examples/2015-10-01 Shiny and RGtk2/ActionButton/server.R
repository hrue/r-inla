library(shiny)

shinyServer(function(input, output) {
  output$distPlot <- renderPlot({
    # Take a dependency on input$goButton
    input$goButton
    
    # Use isolate() to avoid dependency on input$obs
    dist <- isolate(rnorm(input$obs))
    isolate(hist(dist, main=paste("Histogram of ", input$obs," observations")))
  })
}) 