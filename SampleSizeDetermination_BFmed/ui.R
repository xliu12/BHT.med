library(shiny)
library(shinybusy)
# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Sample size determination for testing mediation with the mediation Bayes factor"),
  add_busy_spinner(spin = "fading-circle"),
  br(),
  
  h4("The simple mediation model (standardized) is written as:"),
  p("M = a*X + e_m"),
  p("Y = c'*X + b*M + e_y"),
  
  helpText("a-path: coefficient of the independent variable (X) in the mediator (M) model"),
  helpText("b-path: coefficient of the mediator (M) in the outcome (Y) model"),
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
                   "Standardized path coefficient of the c'-path:",
                   min = -1,
                   max = 1,
                   value = 0),
      radioButtons("whichPriorOdds", "Choice of prior odds specification",
                   choices = c(PriorOdds.a_PriorOdds.b = "a_b",
                               PriorOdds.med_PriorOdds.a = "med_a",
                               PriorOdds.med_PriorOdds.b = "med_b"
                   ),
                   selected = "a_b"),
      conditionalPanel(condition = "input.whichPriorOdds=='a_b' ",
                       numericInput("PriorOdds.a", "Prior odds that the a-path is truly present:", 
                                    min = 0, value = 1),
                       numericInput("PriorOdds.b", "Prior odds that the b-path is truly present:", 
                                    min = 0, value = 1)
                       ),
      conditionalPanel(condition = "input.whichPriorOdds=='med_a' ",
                       numericInput("PriorOdds.med", "Prior odds that the mediation effect is truly present:", 
                                    min = 0, value = 1),
                       numericInput("PriorOdds.a", "Prior odds that the a-path is truly present:", 
                                    min = 0, value = sqrt(0.5)/(1-sqrt(0.5)) ),
                       ),
      conditionalPanel(condition = "input.whichPriorOdds=='med_b' ",
                       numericInput("PriorOdds.med", "Prior odds that the mediation effect is truly present:", 
                                    min = 0, value = 1),
                       numericInput("PriorOdds.b", "Prior odds that the b-path is truly present:", 
                                    min = 0, value = sqrt(0.5)/(1-sqrt(0.5)) ),
      ),
      
    
      
      numericInput(inputId = "TypeI",
                   label = "Type I error rates for determining the mediation Bayes factor cutoff:",
                   value = 0.05,
                   min = 0, max = 1),
      helpText("Note: Type I error rate refers to the probability that the mediation Bayes factor exceeds the cutoff when the null hypothesis [i.e., mediation is absent] is true."),
      numericInput(inputId = "power.desired",
                   label = "Desired power value:",
                   value = 0.8,
                   min = 0, max = 1),
      helpText("Note: Power refers to the probability that the mediation Bayes factor exceeds the cutoff when the alternative hypothesis [i.e., mediation is present] is true."),
      
      numericInput(inputId = "Nmin",
                  label = "The initial sample size to be considered in the sample size determination procedure:",
                  value = 70,
                  min = 0),
      helpText("Note: The sample size determination (SSD) procedure starts from an initial sample size (e.g., 70; a small value). 
      If with this initial sample size, power is lower than the desired power level, 
               the SSD procedure will continue to increase the sample size (by an increment of 'Nstep'), 
               until the sample size is large enough for reaching the desired power level, 
               or unitl the sample size reaches the maximum sample size specified by the user."),
      textInput(inputId = "Nmax",
                   label = "the maximum sample size to be considered in the sample size determination procedure:",
                   value = "500" ),
      numericInput(inputId = "Nstep",
                   label = "a integer specifying the increment in sample size between consecutive rounds of the sample size determination procedure:",
                   value = 10,
                   min = 1),
      numericInput(inputId = "R",
                   label = "the number of simulated samples for constructing the distribution of mediation Bayes factor under each of the competing hypotheses (i.e., the null and the alternative hypotheses):",
                   value = NULL,
                   min = 0),

      actionButton("update", "Go")
      ),
    
    
    mainPanel(
      
      #h1("Sample size determination for testing mediation with the mediation Bayes factor"),
      # power distribution
      # h2("Power distribution for the specific planned sample size"),
      # plotOutput("powerdist.plot"
      #            , width = "100%", height = "400px"
      #            , click = "plot_click", dblclick = NULL
      #            , hover = "plot_hover", 
      #            # hoverDelay = NULL,hoverDelayType = NULL, brush = NULL, clickId = NULL, hoverId = NULL,
      #            inline = FALSE
      # ),
  
      h3("Power and cutoff with the sequence of sample sizes examined in the sample size determination procedure:"),
      dataTableOutput("outpower.Nseq"),
      helpText('In the output table, the column "Cutoff for mediation Bayes factor" lists that with the specified prior odds and sample size in a given row, 
      what is the cutoff value of the mediation Bayes factor associated with the specified Type I error rates.
               That is, with this cutoff value, the probability for the mediation Bayes factor to exceed this cutoff value
               when the null hypothesis is true (i.e., when the mediation effect is absent) is controlled at the specified Type I error rate.'),
      helpText('The column "Power" lists that with the specified prior odds, sample size, and cutoff value in a given row, what is the power of the mediation Bayes factor for detecting a true mediation effect.
              "Power" refers to the probability for the mediation Bayes factor to exceed the cutoff value associated with the specified Type I error rate when the alternative hypothesis is true (i.e., when the mediation effect is present).') ,
      img(src='Untitled.png', align = "center" ,height= 240 ),
      br()
      # # power/assurance curves
      # h2("Power/assurance curves for a range of planned sample sizes"),
      # plotOutput("powercurves"
      #            # , click = "plot_click", dblclick = NULL
      #            # , hover = "plot_hover",  inline = FALSE
      #            )
    )
    )
  )
