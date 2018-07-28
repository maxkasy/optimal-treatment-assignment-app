# to upload, run rsconnect::deployApp('webAppDontRandomize')

library(shiny)
rm(list = ls())
source("dontrandomizefunctions.r")



ui <- fluidPage(
#  titlePanel(h1("Optimal treatment assignment given covariates", align = "center")),
  sidebarLayout(
    sidebarPanel(
      width=6,
      fileInput("file1", "Choose CSV File of covariates",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      hr(),
      h3("Design parameters"),
      fluidRow(
        column(6, radioButtons("simplecomparison", "Estimator",
                   choices = c("Difference of means" = TRUE,
                               "Bayes estimator" = FALSE),
                   selected = TRUE)),
        column(6, radioButtons("sqeprior", "Prior",
                   choices = c("Squared exponential" = TRUE,
                               "Linear model" = FALSE),
                   selected = TRUE))),
      hr(),
      
      fluidRow(
        column(6, numericInput("R", 
                   "Re-randomization draws", 
                   value = 10000, min=1)),
        column(6, numericInput("R2", 
                   "Expected R2", 
                   value = .7, min=0, max=1, step=.1))),
      hr(),
      fluidRow(
        column(6, actionButton(inputId = "calcbutton", label = "Calculate optimal design")),
        column(6, downloadButton("downloadData", "Download optimal design")))
      
    ),
    
    mainPanel(
      width=6,
      textOutput("bestMSE"),
      textOutput("meanMSE"),
      textOutput("MSEgain"),
      hr(),
      tableOutput("designtable")
    )
  )  
)



server <- function(input, output, session) {

  v = reactiveValues()
  v$bestMSE=0
  v$meanMSE=.00001
  
  observeEvent(input$calcbutton,{
    req(input$file1)
    
    #loading covariate file
    covariates=read.table(input$file1$datapath,
                          sep=",")
    
    #setting prior covariance matrix
    if (input$sqeprior) C=Csquaredexponential(covariates, R2=input$R2)
      else C=Clinear(covariates, R2=input$R2)
    
    #calculating optimal design
    optimalDesign=maxEMSE(C,
                          R=input$R,
                          simplecomparison=input$simplecomparison,
                          parallel=FALSE)
    
    covariates$Dstar=optimalDesign$Dstar
    
    v$covariates =covariates
    v$bestMSE=optimalDesign$bestMSE
    v$meanMSE=optimalDesign$meanMSE
  })
  
  output$designtable =  renderTable({
    v$covariates
   })
  
  output$bestMSE = renderText({paste("Optimal MSE: ", format(v$bestMSE, digits=3))})
  output$meanMSE = renderText({paste("Average MSE of randomized designs: ",format( v$meanMSE, digits=3)) })
  output$MSEgain = renderText({paste("Percentage gain:", format((1-v$bestMSE/v$meanMSE) * 100, digits=3), "%")})
  
#download optimal design
  output$downloadData <- downloadHandler(
    filename = "optimaltreatmentassignment.csv",
    content = function(file) {
      write.table(v$covariates, file,  row.names = FALSE, sep = ",")
    }
  )
  

 
}

# Run the app ----
shinyApp(ui = ui, server = server)