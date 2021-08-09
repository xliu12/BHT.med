
# Define UI for application that draws a histogram
shinyUI(fluidPage(

    titlePanel("Null Hypothesis Bayesian Testing of Mediation"),
    
    hr(),
    span(style = "color:black;font-weight:bold;font-size:Large;",
         checkboxInput('IndPaths',
                     'Independent paths a and b?',
                     value = FALSE)),
    # checkboxInput('IndPaths',
    #               'Independent paths a and b?',
    #               value = FALSE),
    helpText("*Path a represents the independent variable to mediator relation; path b represents the mediator to outcome relation."),
    helpText("*Mediation effect is quantified as Med=a*b."),
    
    fluidRow(
        helpText("*Note: in this app the numerators of the prior odds, Bayes factor, and posterior odds are under the alternative hypothesis; and the denominators are under the null hypothesis."),
        
        column(4, align='left',
               h3('Prior belief'),
               h5('prior model odds specifications'),
               conditionalPanel(condition = "input.IndPaths==true",
                                radioButtons("whichPriorOdds", "Choice of prior odds specification",
                                             choices = c(PriorOdds_TwoPaths = "a_b",
                                                         PriorOdds_Med_OnePath = "med_onepath",
                                                         # PriorOdds_Med_bPath = "med_b",
                                                         PriorOdds_Med_q = "med_q"
                                             ),
                                             selected = "a_b"),
                                conditionalPanel(condition = "input.whichPriorOdds=='a_b' ",
                                                 numericInput("PriorOdds.a", "PriorOdds_apath:", 
                                                              min = 0, value = NULL),
                                                 numericInput("PriorOdds.b", "PriorOdds_bpath:", 
                                                              min = 0, value = NULL)),
                                conditionalPanel(condition = "input.whichPriorOdds=='med_onepath' ",
                                                 numericInput("PriorOdds.pathmed", "PriorOdds_Med:", 
                                                              min = 0, value = NULL),
                                                 radioButtons("which_path", "Choose one of the two path prior odds to specify:",
                                                              choices = c(PriorOdds_apath = "PriorOdds.onepatha",
                                                                          PriorOdds_bpath = "PriorOdds.onepathb"
                                                              ),
                                                              selected = "PriorOdds.onepatha"),
                                                 numericInput("value_pathpriorodds", "The specified value:", 
                                                              min = 0, value = NULL)),
                                # conditionalPanel(condition = "input.whichPriorOdds=='med_b' ",
                                #                  numericInput("PriorOdds.bmed", "PriorOdds_Med:", 
                                #                               min = 0, value = NULL),
                                #                  numericInput("PriorOdds.bb", "PriorOdds_bpath:", 
                                #                               min = 0, value = NULL)),
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
                                                 ))
               ),
               conditionalPanel(condition = "input.IndPaths==false",
                                radioButtons("dep_whichPriorOdds", "Choice of prior odds specification",
                                             choices = c(PriorOdds_Med_abpaths = "med_a_b",
                                                         PriorOdds_Med_qs = "med_qs"
                                             ),
                                             selected = "med_qs"),
                                conditionalPanel(condition = "input.dep_whichPriorOdds=='med_a_b' ",
                                                 numericInput("dep_PriorOdds.med", "PriorOdds_Med:", 
                                                              min = 0, value = NULL),
                                                 numericInput("dep_PriorOdds.a", "PriorOdds_apath:", 
                                                              min = 0, value = NULL),
                                                 numericInput("dep_PriorOdds.b", "PriorOdds_bpath:", 
                                                              min = 0, value = NULL)),
                                conditionalPanel(condition = "input.dep_whichPriorOdds=='med_qs' ",
                                                 numericInput("dep_PriorOdds.qmed", "PriorOdds_Med:", 
                                                              min = 0, value = NULL),
                                                 numericInput("dep_q10", "Med0.q10:", 
                                                              min = 0, max = 1, value = NULL),
                                                 numericInput("dep_q01", "Med0.q01:", 
                                                              min = 0, max = 1, value = NULL),
                                                 numericInput("dep_q00", "Med0.q00:", 
                                                              min = 0, max = 1, value = NULL))
               ) ),
        column(4,
               h3('Data evidence'),
               h5('(relative) likelihood of observing the data under the model'),
               conditionalPanel(condition = "input.IndPaths==true",
                                numericInput("BF.a",
                                             "Bayes factor for path a:",
                                             min = 0,
                                             value = 1),
                                numericInput("BF.b",
                                             "Bayes factor for path b:",
                                             min = 0,
                                             value = 1)
               ),
               conditionalPanel(condition = "input.IndPaths==false",
                                numericInput("BIC_M11",
                                             "BIC_M11: BIC value of the unconstrained model where the paths a and b are both unconstrained:",
                                             min = 0,
                                             value = 1),
                                numericInput("BIC_M01",
                                             "BIC_M01: BIC value of the model where the path a is constrained at zero and path b is unconstrained:",
                                             min = 0,
                                             value = 1),
                                numericInput("BIC_M10",
                                             "BIC_M10: BIC value of the model where the path a is unconstrained and path b is constrained at zero:",
                                             min = 0,
                                             value = 1),
                                numericInput("BIC_M00",
                                             "BIC_M00: BIC value of the model where a and b are both constrained at zero:",
                                             min = 0,
                                             value = 1)
               ) )
    ),
    
    br(),# a line break
    h3("Based on the prior odds/probability specifications, the mediation Bayes factor and posterior odds are"
       ,align='left'),
    fluidRow(column(8, align='left',
                    span(style = "color:black;font-weight:bold;font-size:large;",
                         tableOutput("tabout_odds") )
                    )
             ),
    # tableOutput("tabout"),
    # br(),# a line break
    # fluidRow(column(8, align='left',
    #                 span(style = "color:black;font-weight:bold;font-size:large;",
    #                      tableOutput("tabout_prob") )
    # )
    # ),
    
    br(),# a line break
    br()
    # 
    # conditionalPanel(condition = "input.IndPaths==true",
    #                  h4(align='center', "Explore the impact of prior odds specifications"),
    #                  column(6, align='center',
    #                         h5("Impact of PriorOdds.a, PriorOdds.b"),
    #                         plotOutput("plot_ab")),
    #                  column(6, align='center',
    #                         h5("Impact of PriorOdds.med, PriorOdds.a"),
    #                         plotOutput("plot_meda"))
    #                  # ,column(4,
    #                  #        h5("Impact of PriorOdds.med, Med0.q10"),
    #                  #        plotOutput("plot_medq"))
    #                  )
))
