library(shiny)
library(shinybusy)
# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Sample size determination for testing mediation with the mediation Bayes factor"),
  add_busy_spinner(spin = "fading-circle"),
  br(),
  
  helpText("X: independent variable (e.g., treatment assignment)"),
  helpText("M: mediator"),
  helpText("Y: outcome"),
  
  helpText("a-path: coefficient of independent variable (X) in the mediator (M) model"),
  helpText("b-path: coefficient of mediator (M) in the outcome (Y) model"),
  helpText("The mediation effect is quantified as a*b"),
  br(),
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      numericInput("std.a",
                   "Standardized path coefficient of the a-path if the a-path is truly present (i.e., non-zero):",
                   min = -1,
                   max = 1,
                   value = 0.39),
      numericInput("std.b",
                   "Standardized path coefficient of the b-path if the b-path is truly present (i.e., non-zero):",
                   min = -1,
                   max = 1,
                   value = 0.39),
      numericInput("std.cp",
                   "Standardized path coefficient of the cp-path:",
                   min = -1,
                   max = 1,
                   value = 0),
      
      helpText("If considering the uncertainty about the specified coefficients for a-path, b-path, and cp-path, please check the box below and then provide the standard deviations for these path coefficients."),
      checkboxInput("uncertain_effect",
                    "Considering the uncertainty about the specified coefficients for a-path, b-path, and cp-path?",
                    value = FALSE),
      conditionalPanel(condition = "input.uncertain_effect == true",
                       numericInput("sigma.a", 
                                    "standard deviation of the coefficient of a-path:", min = 0, value = 0),
                       numericInput("sigma.b", 
                                    "standard deviation of the coefficient of b-path:", min = 0, value = 0),
                       numericInput("sigma.cp", 
                                    "standard deviation of the coefficient of cp-path:", min = 0, value = 0)
                       ),
      
      numericInput("p",
                   "The number of baseline covariates:",
                   min = 0,
                   # max = 1,
                   value = 0),
      helpText("If baseline covariates are involved, please provide the proportions of variances explained by the covariates for X, M, and Y (see below)."),
      helpText("Note: These proportions of variances are those obtained without adjustment for any variables (such as from a regression of Y on the covariates)"),
      numericInput("Rsq.xc",
                   "The proportion of variance in the independent variable (X) explained by the baseline covariates:",
                   min = 0,
                   # max = 1,
                   value = 0),
      numericInput("Rsq.mc",
                   "The proportion of variance in the mediator (M) explained by the baseline covariates:",
                   min = 0,
                   # max = 1,
                   value = 0),
      numericInput("Rsq.yc",
                   "The proportion of variance in the outcome (Y) explained by the baseline covariates:",
                   min = 0,
                   # max = 1,
                   value = 0),
      
      helpText("Design prior odds values for a-path and for b-path are needed for generating samples under the null hypothesis of no mediation. Please specify these values below."),
      numericInput("Design.PriorOdds.a", 
                   "Design prior odds for the presence of a-path:", 
                   min = 0, 
                   value = 1),
      numericInput("Design.PriorOdds.b", 
                   "Design prior odds for the presence of b-path:", 
                   min = 0, 
                   value = 1),
      
      helpText("Analysis prior odds values for a-path and for b-path are needed for calculating the mediation Bayes factor (BFmed) with each generated sample. Please specify these values below."),
      numericInput("Analysis.PriorOdds.a", 
                   "Analysis prior odds for the presence of a-path:", 
                   min = 0, 
                   value = 1),
      numericInput("Analysis.PriorOdds.b", 
                   "Analysis prior odds for the presence of b-path:", 
                   min = 0, 
                   value = 1),
      
      helpText("BFmed > cutoff is interpreted as support for the alternative hypothesis that mediation is present. The cutoff can be defined as an absolute value, such as 3; or be definited as a relative value that corresponds to certain false positive rate, such as 5% false positive rate. If a user does not specify the absolute cutoff or the relative cutoff, the defaults are used."),
      numericInput("cutoff.BF", 
                   "Absolute cutoff for interpreting a Bayes factor value as supporting the alternative hypothesis:", 
                   min = 1, 
                   value = 3),
      numericInput("cutoff.FPR", 
                   "False positive rate (used by the relative cutoff):", 
                   min = 0, 
                   value = 0.05),
      
      helpText("Please provide a sequence of candidate sample sizes below."),
      numericInput(inputId = "Nmin",
                  label = "The smallest sample size to be considered:",
                  value = 100,
                  min = 0),
      numericInput(inputId = "Nmax",
                   label = "The largest sample size to be considered:",
                   value = 120),
      numericInput(inputId = "Nstep",
                   label = "a integer specifying the increment in the sequence of sample sizes:",
                   value = 10,
                   min = 1),
      numericInput(inputId = "R",
                   label = "The number of replications:",
                   value = 1000,
                   min = 0),

      actionButton("update", "Go")
      ),
    
    
    mainPanel(
      
      #h1("Sample size determination for testing mediation with the mediation Bayes factor"),
      
      h3("True positive rates and false positive rates:"),
      dataTableOutput("outpower.Nseq"),
      
      # img(src='Untitled.png', align = "center" ,height= 240 ),
      br(),
      # # power/assurance curves
      h3("Plot of the true positive rates"),
      plotOutput("plot.TPR"
                 # , click = "plot_click", dblclick = NULL
                 # , hover = "plot_hover",  inline = FALSE
                 ),
      br(),
      h3("Plot of the false positive rates"),
      plotOutput("plot.FPR"
                 # , click = "plot_click", dblclick = NULL
                 # , hover = "plot_hover",  inline = FALSE
      ),
      br()
    )
    )
  )
