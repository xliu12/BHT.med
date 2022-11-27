library(parallel)
library(ggplot2)
library(shiny)

source('app_bfmedpower.R')

server <- function(input, output) {
  # power distribution ----
  SSD.reactive = eventReactive(input$update,{
    results=SSD.bfmed(Power.desired=input$power.desired, Nmin=input$Nmin, Nmax=input$Nmax,Nstep=input$Nstep,
                      a=input$std.a,b=input$std.b,cp=input$std.cp, 
                      PriorOdds.a=input$PriorOdds.a, PriorOdds.b=input$PriorOdds.b, PriorOdds.med=input$PriorOdds.med,
                      TypeI=input$TypeI, R=input$R, seed=816051)
    results
  })

  output$outpower.Nseq <- renderDataTable({
    withProgress(message = "In progress", value = 0,
                 detail = NULL, expr = {
      powerdist=SSD.reactive()
    })
    
    df = round( SSD.reactive()$outpower.Nseq[ ,c(1:3, 7, 11,12)] ,3 )
    colnames(df) = c('PriorOdds.a', 'PriorOdds.b', 'PriorOdds.med', 'Planned sample size (N)',
                     'Cutoff for mediation Bayes factor', 'Power')
    df
  }) 

}

  
