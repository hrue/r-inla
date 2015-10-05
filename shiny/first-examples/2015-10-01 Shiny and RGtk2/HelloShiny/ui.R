library(shiny)

shinyUI(fluidPage(
  headerPanel("Hello Shiny"),
  sidebarPanel(
    sliderInput("obs",  "Number of observations:", 0, 1000, 500)
  ),
  mainPanel(
    plotOutput("distPlot")
  )
))