library(shiny)

shinyUI(navbarPage("A Simple R-INLA tutorial",
                   tabPanel("R-INLA", id=3,
                            fluidRow(
                              column(4,
                                     wellPanel(selectInput("latent", label = "Latent models", 
                                                           width='100%',
                                                           choices = list("NULL" = 1, 
                                                                          "iid" = 2,
                                                                          "rw2" = 3), 
                                                           selected = 1)),
                                     wellPanel(selectInput("likelihood", label = "Likelihoods", 
                                                           width='100%',
                                                           choices = list("nbinomial" = 1, 
                                                                          "gaussian" = 2), 
                                                           selected = 1)),
                                     wellPanel(radioButtons("radio", label = "Show simulate data",
                                                            choices = list("No" = 1, 
                                                                           "Yes" = 2),
                                                            selected = 1,
                                                            inline = TRUE))
                                     ),
                              
                              column(8,
                                     conditionalPanel(condition = "input.radio == 2",
                                       wellPanel(plotOutput("dataPlot"))),
                                     wellPanel(htmlOutput("inlafuncs", container = span)),
                                     wellPanel(verbatimTextOutput("inlares"))
                                     
                                     
                                     ) # end of column
                            ) # end of fluidRow
                            
                   ),
                   tabPanel("SPDE",
                            wellPanel("Introduce the SPDE approaches"),
                            fluidRow(
                              
                              column(4,
                                     wellPanel("Select the Examples"),
                                     wellPanel("Numeric Input for max.edge"),
                                     wellPanel("Numeric Input for offset"),
                                     wellPanel("Slide Input for cutoff")
                              ), # end of column
                              
                              column(8,
                                     wellPanel("Show the code to call the data"),
                                     wellPanel("Show the code to call inla.mesh.2d()"),
                                     wellPanel("Show the plot")
                              ) # end of column
                            ) # end o fluidrow
                            
                   )
)
)