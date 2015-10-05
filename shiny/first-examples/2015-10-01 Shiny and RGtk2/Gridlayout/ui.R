library(shiny)

shinyUI(fluidPage(
  titlePanel("Grid Layout"),
  fluidRow(column(width = 12, HTML('<p style="color: #ffffff; background-color: #E6E3E3">Fluid 12</p>'),
                  fluidRow(column(width = 6, HTML('<p style="color: #ffffff; background-color: #ABA8A8">Fluid 6</p>')),
                           column(width = 6, HTML('<p style="color: #ffffff; background-color: #646363">Fluid 6</p>')))
           )
    )
))