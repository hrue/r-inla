library(shiny)

shinyUI(navbarPage("A Simple R-INLA tutorial",
                   tabPanel("R-INLA", id=3,
                            fluidRow(
                              column(4,
                                     wellPanel(style = "background-color: #CEF6E3;",
                                       selectInput("latent", label = "Latent models", 
                                                           width='100%',
                                                           choices = list("NULL" = 1, 
                                                                          "iid" = 2,
                                                                          "rw2" = 3), 
                                                           selected = 1)),
                                     wellPanel(style = "background-color: #ECCEF5;",
                                       selectInput("likelihood", label = "Likelihoods", 
                                                           width='100%',
                                                           choices = list("nbinomial" = 1, 
                                                                          "gaussian" = 2), 
                                                           selected = 1)),
                                     wellPanel(style = "background-color: #FBEFF5;",
                                       radioButtons("radio", label = "Show simulate data",
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
                            HTML('For a two dimentional mesh, we have a main function <font color="red">inla.mesh.2d(&nbsp)</font> that is recommended to use for building a mesh. This function creates the Constrained Refined Delaunay Triangulation (CRDT) that we just call mesh. There are a several options:'),
                            br(),
                            HTML('&nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp inla.mesh.2d(loc = NULL, loc.domain = NULL, offset = NULL, n = NULL,<br>
    &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp&nbsp &nbsp &nbsp &nbsp boundary = NULL, interior = NULL, max.edge, min.angle = NULL,<br>
     &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp&nbsp &nbsp &nbsp &nbsp cutoff = 1e-12, plot.delay = NULL)'),
                            br(),
                            fluidRow(
                              
                              column(4,
                                     selectInput("trang_example", label = "Select examples", 
                                                 width='100%',
                                                 choices = list("Toy example 1" = 1, 
                                                                "Toy example 2" = 2), 
                                                 selected = 1),
                                     wellPanel(style = "background-color: #F5BCA9;",
                                       HTML('<p><b>max.edge:</b></p>'),
                                       fluidRow(
                                         column(6, numericInput("trang_maxedge1",label=NULL, value = 0.1, min=0.1, max=5, step=0.01)),
                                         column(6, numericInput("trang_maxedge2",label=NULL, value = 0.1, min=0.1, max=5, step=0.01))
                                       )
                                     ),
                                     wellPanel(
                                       HTML('<p><b>offset:</b></p>'),
                                       fluidRow(
                                         column(6, numericInput("trang_offset1",label=NULL, value = 0.1, min=-10, max=10, step=0.01)),
                                         column(6, numericInput("trang_offset2",label=NULL, value = 0.1, min=-10, max=10, step=0.01))
                                       )
                                     ),
                                     wellPanel(
                                       sliderInput("cutoff", label = "cutoff", min = 0, max = 1, value = 0.1,step=0.1)
                                     )
                              ), # end of column
                              
                              column(8,
                                     wellPanel(htmlOutput("datas", container = span)),
                                     wellPanel(textOutput("funcs1")),
                                     wellPanel(plotOutput("spdeplot"))
                              ) # end of column
                            ) # end o fluidrow
                            
                   )
)
)