library(shiny)

shinyUI(navbarPage("A Simple R-INLA tutorial on SPDE models",
                   tabPanel("The SPDE approach",
                              titlePanel("Prior Settings"),
                              sidebarLayout(
                                sidebarPanel(),
                                mainPanel(plotOutput("plot_prior"))
                              )
                   ),
                   tabPanel("The triangulation",
                            HTML('For a two dimentional mesh, we have a main function <font color="red">inla.mesh.2d(&nbsp)</font> that is recommended to use for building a mesh. This function creates the Constrained Refined Delaunay Triangulation (CRDT) that we just call mesh. There are a several options:'),
                            br(),
                            HTML('&nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp inla.mesh.2d(loc = NULL, loc.domain = NULL, offset = NULL, n = NULL,<br>
    &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp&nbsp &nbsp &nbsp &nbsp boundary = NULL, interior = NULL, max.edge, min.angle = NULL,<br>
     &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp &nbsp&nbsp &nbsp &nbsp &nbsp cutoff = 1e-12, plot.delay = NULL)</div>'),
                            br(),
                            fluidRow(
                              
                              column(4,
                                     selectInput("trang_example", label = "Select examples", 
                                                 width='100%',
                                                 choices = list("Toy example" = 1, 
                                                                "Choice 2" = 2,
                                                                "Choice 3" = 3), 
                                                 selected = 1),
                                     wellPanel(
                                       HTML('<p><b>max.edge:</b></p>'),
                                       fluidRow(
                                         column(6, numericInput("trang_maxedge1",label=NULL, value = 0, min=0, max=100)),
                                         column(6, numericInput("trang_maxedge2",label=NULL, value = 0, min=0, max=100))
                                       )
                                       )
                                     ),
                              
                              column(8,
                                     htmlOutput("datas", container = span),
                                     textOutput("funcs")
                                     )
                              )
                   )
))