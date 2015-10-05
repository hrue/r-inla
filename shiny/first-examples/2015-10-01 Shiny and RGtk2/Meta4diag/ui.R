library(shiny)

shinyUI(navbarPage("Meta4diag",
                   tabPanel("Prior",
                              titlePanel("Prior Settings"),
                              sidebarLayout(
                                sidebarPanel(),
                                mainPanel(plotOutput("plot_prior"))
                              )
                   ),
                   tabPanel("Data",
                              titlePanel("Upload Data"),
                              sidebarLayout(
                                sidebarPanel(
                                  fileInput('data', 'Choose A File', accept='.csv')),
                                mainPanel(tableOutput("contents"))
                              )
                   )
))