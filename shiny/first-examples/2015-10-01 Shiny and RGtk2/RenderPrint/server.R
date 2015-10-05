library(shiny)

shinyServer(function(input, output) {
  output$distText <- renderPrint({
    round(summary(rnorm(input$obs)),2)
  })
}) 