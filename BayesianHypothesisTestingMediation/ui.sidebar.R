
library(shiny)
# Define UI for application that draws a histogram
shinyUI(fluidPage(
    # Application title
    titlePanel("Null Hypothesis Bayesian Testing of Mediation"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            helpText("*Note: in this app the numerators of the prior odds, Bayes factor, and posterior odds for an effect are under the alternative hypothesis that the effect is non-zero; and the denominators are under the null hypothesis that the effect is zero."),
            # independent path
            checkboxInput('IndPaths',
                          'Are the paths a and b independent?',
                          value = TRUE),
            conditionalPanel(condition = "input.IndPaths==true",
                             radioButtons("whichPriorOdds", "Choice of prior odds specification",
                                          choices = c(PriorOdds_TwoPaths = "a_b",
                                                      PriorOdds_Med_aPath = "med_a",
                                                      PriorOdds_Med_bPath = "med_b",
                                                      PriorOdds_Med_q = "med_q"
                                                      ),
                                          selected = "a_b"),
                             conditionalPanel(condition = "input.whichPriorOdds=='a_b' ",
                                              numericInput("PriorOdds.a", "PriorOdds_apath:", 
                                                           min = 0, value = NULL),
                                              numericInput("PriorOdds.b", "PriorOdds_bpath:", 
                                                           min = 0, value = NULL)),
                             conditionalPanel(condition = "input.whichPriorOdds=='med_a' ",
                                              numericInput("PriorOdds.amed", "PriorOdds_Med:", 
                                                           min = 0, value = NULL),
                                              numericInput("PriorOdds.aa", "PriorOdds_apath:", 
                                                           min = 0, value = NULL)),
                             conditionalPanel(condition = "input.whichPriorOdds=='med_b' ",
                                              numericInput("PriorOdds.bmed", "PriorOdds_Med:", 
                                                           min = 0, value = NULL),
                                              numericInput("PriorOdds.bb", "PriorOdds_bpath:", 
                                                           min = 0, value = NULL)),
                             conditionalPanel(condition = "input.whichPriorOdds=='med_q' ",
                                              numericInput("PriorOdds.qmed", "PriorOdds_Med:", 
                                                           min = 0, value = NULL),
                                              radioButtons("which_q", "Choose one of the conditional probabilities of the three no mediation scenarios to specify:",
                                                           choices = c(Med0.q10 = "q10",
                                                                       Med0.q01 = "q01",
                                                                       Med0.q00 = "q00"
                                                           ),
                                                           selected = "q10"),
                                              numericInput("value_q", "The specified value:", 
                                                           min = 0, max = 1, value = NULL),
                                              conditionalPanel(condition = "input.which_q=='q00' ",
                                                               radioButtons('which.larger', "which one is larger, Med0.q10 or Med0.q01?",
                                                                            choices = c(Med0.q10 = "10",
                                                                                        Med0.q01 = "01"))
                                                               )
                                              ),
                             numericInput("BF.a",
                                          "Bayes factor for path a:",
                                          min = 0,
                                          value = 1),
                             numericInput("BF.b",
                                           "Bayes factor for path b:",
                                           min = 0,
                                           value = 1)
            ),
            
            # dependent path
            conditionalPanel(condition = "input.IndPaths==false",
                             radioButtons("dep_whichPriorOdds", "Choice of prior odds specification",
                                          choices = c(PriorOdds_Med_abpaths = "med_a_b",
                                                      PriorOdds_Med_qs = "med_qs"
                                          ),
                                          selected = "med_qs"),
                             conditionalPanel(condition = "input.dep_whichPriorOdds=='med_a_b' ",
                                              numericInput("PriorOdds.med", "PriorOdds_Med:", 
                                                           min = 0, value = NULL),
                                              numericInput("PriorOdds.a", "PriorOdds_apath:", 
                                                           min = 0, value = NULL),
                                              numericInput("PriorOdds.b", "PriorOdds_bpath:", 
                                                           min = 0, value = NULL)),
                             conditionalPanel(condition = "input.dep_whichPriorOdds=='med_qs' ",
                                              numericInput("q10", "Med0.q10:", 
                                                           min = 0, max = 1, value = NULL),
                                              numericInput("q01", "Med0.q01:", 
                                                           min = 0, max = 1, value = NULL),
                                              numericInput("q00", "Med0.q00:", 
                                                           min = 0, max = 1, value = NULL)),
                             numericInput("BIC_M11",
                                          "BIC value of the unconstrained model where the paths a and b are both unconstrained:",
                                          min = 0,
                                          value = 1),
                             numericInput("BIC_M01",
                                          "BIC value of the model where the path a is constrained at zero and path b is unconstrained:",
                                          min = 0,
                                          value = 1),
                             numericInput("BIC_M10",
                                          "BIC value of the model where the path a is unconstrained and path b is constrained at zero:",
                                          min = 0,
                                          value = 1),
                             numericInput("BIC_M00",
                                          "BIC value of the model where a and b are both constrained at zero:",
                                          min = 0,
                                          value = 1)
            )
            
        ),

        # Show a plot of the generated distribution
        mainPanel(
            br(),# a line break
            h4("Based on the prior odds specifications,"),
            h4("the mediation Bayes factor and posterior odds are"),
            dataTableOutput("tabout")
        )
    )
))
