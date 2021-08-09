#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

source('fun_NHBTmed.R') 
source('fun_plot_PriorOdds.R')

library(shiny)
library(BayesFactor)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {

    output$tabout_odds <- renderTable(rownames = T, striped = T, {
        if( input$IndPaths ){
            if( input$whichPriorOdds=='a_b' ){
                out=pathb.a(input$PriorOdds.a, input$PriorOdds.b, input$BF.a, input$BF.b)
            }
            if( input$whichPriorOdds=='med_onepath' ){
                if( input$which_path=='PriorOdds.onepatha' ){
                    out=med.a(input$PriorOdds.pathmed, input$value_pathpriorodds, input$BF.a, input$BF.b)
                }
                if( input$which_path=='PriorOdds.onepathb' ){
                    out=med.b(input$PriorOdds.pathmed, input$value_pathpriorodds, input$BF.a, input$BF.b) 
                }
            }
            # if( input$whichPriorOdds=='med_b' ){
            #     out=med.b(input$PriorOdds.bmed, input$PriorOdds.bb, input$BF.a, input$BF.b)
            # }
            if( input$whichPriorOdds=='med_q' ){
                if( input$which_q=='q10' ){
                    out=med.q10(input$PriorOdds.qmed, input$value_q, input$BF.a, input$BF.b)
                }
                if( input$which_q=='q01' ){
                    out=med.q01(input$PriorOdds.qmed, input$value_q, input$BF.a, input$BF.b)
                }
                if( input$which_q=='q00' ){
                    out=med.q00(input$PriorOdds.qmed, input$value_q, input$which.larger, input$BF.a, input$BF.b)
                }
            }
        }
        
        if( !input$IndPaths ){
            if( input$dep_whichPriorOdds=="med_a_b" ){
                out=med.ab(input$dep_PriorOdds.med, input$dep_PriorOdds.a, input$dep_PriorOdds.b,
                           BIC11=input$BIC_M11, BIC01=input$BIC_M01, BIC10=input$BIC_M10, BIC00=input$BIC_M00 )
            }
            if( input$dep_whichPriorOdds=="med_qs" ){
                out=med.qs(input$dep_PriorOdds.qmed, H0.q10=input$dep_q10, H0.q01=input$dep_q01, H0.q00=input$dep_q00,
                           BIC11=input$BIC_M11, BIC01=input$BIC_M01, BIC10=input$BIC_M10, BIC00=input$BIC_M00 )
            }
        }
        
        out1=matrix(out,3,3, dimnames = list(c('Path a','Path b','Mediation'),c('Prior Odds','Bayes Factor','Posterior Odds')))  
        out1=as.data.frame(out1)
        out1$`Posterior Prob`=(out1$`Posterior Odds`)/(1+out1$`Posterior Odds`) 
        out1
    })
    
    # output$tabout_prob <- renderTable(rownames = T, striped = T,  {
    #     if( input$IndPaths ){
    #         if( input$whichPriorOdds=='a_b' ){
    #             out=pathb.a(input$PriorOdds.a, input$PriorOdds.b, input$BF.a, input$BF.b)
    #         }
    #         if( input$whichPriorOdds=='med_onepath' ){
    #             if( input$which_path=='PriorOdds.onepatha' ){
    #                 out=med.a(input$PriorOdds.pathmed, input$value_pathpriorodds, input$BF.a, input$BF.b)
    #             }
    #             if( input$which_path=='PriorOdds.onepathb' ){
    #                 out=med.b(input$PriorOdds.pathmed, input$value_pathpriorodds, input$BF.a, input$BF.b) 
    #             }
    #         }
    #         if( input$whichPriorOdds=='med_q' ){
    #             if( input$which_q=='q10' ){
    #                 out=med.q10(input$PriorOdds.qmed, input$value_q, input$BF.a, input$BF.b)
    #             }
    #             if( input$which_q=='q01' ){
    #                 out=med.q01(input$PriorOdds.qmed, input$value_q, input$BF.a, input$BF.b)
    #             }
    #             if( input$which_q=='q00' ){
    #                 out=med.q00(input$PriorOdds.qmed, input$value_q, input$which.larger, input$BF.a, input$BF.b)
    #             }
    #         }
    #         out1=matrix(out,3,3, dimnames = list(c('Path a','Path b','Mediation'),c('Prior Odds','Bayes Factor','Posterior Odds')))
    #         posterior_prob_H1med = 1/(1/out1[3,3]+1)
    #         posterior_prob_H0med = 1-posterior_prob_H1med
    #     }
    # 
    #     if( !input$IndPaths ){
    #         if( input$dep_whichPriorOdds=="med_a_b" ){
    #             out=med.ab(input$dep_PriorOdds.med, input$dep_PriorOdds.a, input$dep_PriorOdds.b,
    #                        BIC11=input$BIC_M11, BIC01=input$BIC_M01, BIC10=input$BIC_M10, BIC00=input$BIC_M00 )
    #         }
    #         if( input$dep_whichPriorOdds=="med_qs" ){
    #             out=med.qs(input$dep_PriorOdds.qmed, H0.q10=input$dep_q10, H0.q01=input$dep_q01, H0.q00=input$dep_q00,
    #                        BIC11=input$BIC_M11, BIC01=input$BIC_M01, BIC10=input$BIC_M10, BIC00=input$BIC_M00 )
    #         }
    #         # out1=matrix(out,1,3, dimnames = list(c('Mediation'),c('Prior Odds','Bayes Factor','Posterior Odds')))
    #         out1=matrix(out,3,3, dimnames = list(c('Path a','Path b','Mediation'),c('Prior Odds','Bayes Factor','Posterior Odds')))
    #         posterior_prob_H1med = 1/(1/out1[3,3]+1)
    #         posterior_prob_H0med = 1-posterior_prob_H1med
    #     }
    #     posterior_prob=c(posterior_prob_H1med, posterior_prob_H0med)
    #     out2=matrix(round(posterior_prob,2), 2, 1, 
    #                 dimnames = list(c('non-zero mediation effect','no mediation effect'),
    #                                 c('Posterior Probability')))
    #     data.frame(out2)
    # })
    
    # # explore impact of prior odds independent paths
    # output$plot_ab <- renderPlot(
    #     plot.ab(input$BF.a, input$BF.b)
    # )
    # output$plot_meda <- renderPlot(
    #     plot.meda(input$BF.a, input$BF.b)
    # )
    # # output$plot_medq <- renderPlot(
    # #     plot.medq(input$BF.a, input$BF.b)
    # # )
    
})
