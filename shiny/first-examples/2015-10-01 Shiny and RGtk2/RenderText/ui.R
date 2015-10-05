library(shiny)

shinyUI(pageWithSidebar(
  headerPanel("textOutput    &    renderText"),
  sidebarPanel(
    sliderInput("obs", "Number of observations:",
                min = 0, max = 1000, value = 500)
  ),
  mainPanel(
    wellPanel(textOutput("distText"))
  )
))