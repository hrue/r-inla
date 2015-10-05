library(shiny)

shinyUI(pageWithSidebar(
  headerPanel("tableOutput    &    renderTable"),
  sidebarPanel(
    sliderInput("obs", "Number of observations to simulate:",
                min = 10, max = 1000, value = 50),
    sliderInput("header", "Number of observations to view:",
                min = 1, max = 10, value = 1)
  ),
  mainPanel(
    wellPanel(tableOutput("table"))
  )
))