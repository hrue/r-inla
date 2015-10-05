library(shiny)
library(INLA)

shinyUI(pageWithSidebar(
  
  headerPanel("INLA: Meta-analysis of Diagnostic Test"),
  sidebarPanel(
    checkboxInput("GoldStandard", "Gold-Standard", TRUE),
    br(),
    selectInput("dataset", "Examples: Choose a dataset", 
                choices = c("Telomerase", "ScheidlerLag","NMP22")),
    
    br(),
    br(),
    checkboxInput("prev", "Prevalence", FALSE),
    br(),
    checkboxInput("ppvnpv", "Show PPV&NPV model", FALSE),
    br(),
    fileInput(inputId = "iFile", label = "",
              accept=c('text/csv', 'text/comma-separated-values,text/plain'))
    ),
  mainPanel(
            tabsetPanel(
              tabPanel("Posterior Density",
                       plotOutput('plot') ),
              tabPanel("Model Summary",
                       verbatimTextOutput("summary") ),
              tabPanel("Summary Plot",
                       plotOutput('sumplot') ),
              tabPanel("Data",
                       tableOutput("view"))
              ,id="mainTabUI"))
))


