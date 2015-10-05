library(shiny)

shinyServer(function(input, output) {
  
  dataInput = reactive({
    x = rnorm(input$obs)
    y = rnorm(input$obs)
    data.frame(x=x,y=y)
  })
  
  output$table <- renderTable({
    head(dataInput(),n=input$header)
  })
}) 