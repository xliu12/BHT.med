library(parallel)
library(ggplot2)
library(shiny)

source("functions_BFmed_SSD.R")

server <- function(input, output) {
  # power distribution ----
  SSD.reactive <- eventReactive(input$update, {
    
    results <- BFmed.SSD(
      N = seq(input$Nmin, input$Nmax, by = input$Nstep), 
      
      std.a = input$std.a, std.b = input$std.b, std.cp = input$std.cp, 
      uncertain.effect = input$uncertain_effect,
      sigma.a = input$sigma.a, sigma.b = input$sigma.b, sigma.cp = input$sigma.cp, 
      
      p = input$p, # number of covariates C's
      Rsq.xc = input$Rsq.xc, Rsq.mc = input$Rsq.mc, Rsq.yc = input$Rsq.yc, #  respectively the proportions of the variance explained by all C's for the treatment X, mediator M, and outcome Y, without controlling for any variables.
      
      Design.PriorOdds.a = input$Design.PriorOdds.a, Design.PriorOdds.b = input$Design.PriorOdds.b,
      Analysis.PriorOdds.a = input$Analysis.PriorOdds.a, Analysis.PriorOdds.b = input$Analysis.PriorOdds.b,
      
      cutoff.BF = input$cutoff.BF, absolute.cutoff = TRUE,
      cutoff.FPR = input$cutoff.FPR, relative.cutoff = TRUE,
      R = input$R, seed = 12345
    )
    
    results
  })

  output$outpower.Nseq <- renderDataTable({
    withProgress(message = "In progress", value = 0,
                 detail = NULL, expr = {
      powerdist <- SSD.reactive()
    })
    
    out <- SSD.reactive()
    out <- round(out, 2)
    out
  })
  
  output$plot.TPR <- renderPlot({
    out <- SSD.reactive()
    fig <- plot.true.positive.med(out)
    fig
  })
  
  output$plot.FPR <- renderPlot({
    out <- SSD.reactive()
    fig <- plot.false.positive.med(out)
    fig
  })

}

  
