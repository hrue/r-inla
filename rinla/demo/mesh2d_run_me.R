library(INLA)
library(shiny)

ui <- fluidPage(
    headerPanel("Two dimensional mesh demostration app"), 
    sidebarLayout(
        sidebarPanel(
            checkboxInput("domain", "Add domain: [0,1]x[0,1]", value=FALSE),
            sliderInput("N", 
                        label="Number of random locatins to drawn", 
                        min=1, max=100, value=30, step=1), 
            sliderInput("offset", 
                        label="The automatic extension distance, 1 or 2 values", 
                        min=0.1, max=2, value=c(0.5, 1), step=0.05), 
            sliderInput("n", 
                        label="Initial nodes in the automatic extensions.",
                        min=3, max=100, value=16, step=1),
            sliderInput("max.edge",
                        label="Largest allowed triangle edge length, 1 or 2 v.",      
                        min=0.05, max=1, value=c(0.3, 0.7), step=0.05),
            sliderInput("cutoff", 
                        label="The minimum allowed distance between points",
                        min=1e-10, max=1, value=0.1, step=0.01), 
            actionButton("Exit", "Exit"), 
            width=4), 
        mainPanel(plotOutput("plot"))
        )
    )

server <- function(input, output) {
    output$plot <- renderPlot({
        loc <- cbind(runif(input$N), runif(input$N))
        loc.dom <- NULL
        if(input$domain)
             loc.dom <- cbind(x=c(0,1,1,0,0), y=c(0,0,1,1,0))
        mesh <- inla.mesh.2d(loc, loc.domain=loc.dom, max.edge=input$max.edge, 
                             offset=input$offset, n=input$n, cutoff=input$cutoff)
        par(mar=c(0,0,1,0), mgp=c(1,0.5,0))
        plot(mesh, asp=1);
        points(loc, pch=19)
    }, width=500, height=500)
    observe({
        if(input$Exit > 0){
            stopApp(NULL)
        }
    })
}

runApp(shinyApp(ui=ui, server=server))
