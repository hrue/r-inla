library(shiny)

shinyUI(navbarPage("A Simple R-INLA tutorial",
                   tabPanel("R-INLA",
                            fluidRow(
                              column(4,
                                     wellPanel(selectInput("latent", label = "Latent models", 
                                                           width='100%',
                                                           choices = list("NULL" = 1, 
                                                                          "iid" = 2,
                                                                          "rw2" = 3), 
                                                           selected = 1)),
                                     wellPanel("Select Input for Likelihoods"),
                                     wellPanel("Radio to show wether to plot the data")
                                     ),
                              
                              column(8,
                                     wellPanel("If choose to show the data, plot the data"),
                                     wellPanel("Write the INLA code"),
                                     wellPanel("Summary the result from INLA")
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