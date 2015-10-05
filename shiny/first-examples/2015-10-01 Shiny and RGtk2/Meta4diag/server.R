library(shiny)

shinyServer(function(input, output) {
  myData <- reactive({
    inFile <- input$data
    if (is.null(inFile)) return(NULL)
    data <- read.csv(inFile$datapath,header=T)
    return(data)
  })
  output$contents <- renderTable({
    myData()
  })
})