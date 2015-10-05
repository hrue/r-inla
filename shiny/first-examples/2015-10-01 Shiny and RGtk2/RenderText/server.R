library(shiny)

shinyServer(function(input, output) {
  output$distText <- renderText({
    paste("Number of observation is ",input$obs,sep="")
  })
}) 