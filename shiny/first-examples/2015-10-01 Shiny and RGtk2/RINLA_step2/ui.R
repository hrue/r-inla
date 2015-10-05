library(shiny)

shinyUI(navbarPage("A Simple R-INLA tutorial",
                   tabPanel("R-INLA",
                            fluidRow(
                              column(4,
                                     wellPanel("Sidebar Panel")
                              ),
                              
                              column(8,
                                     wellPanel("Main Panel")
                              )
                            )
                            
                            ),
                   tabPanel("SPDE")
)
)