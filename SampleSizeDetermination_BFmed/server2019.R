library(parallel)
library(ggplot2)
library(shiny)

source('fun_powercurve.R')

server <- function(input, output) {
  # power distribution ----
  pd.reactive = eventReactive(input$update,{
    powerdist=pd.jointind(
      Nplan=input$Nplan,
      nrawdata=input$nrawdata,
      a_true=input$std.ahat, 
      b_true=input$std.bhat, 
      cp_true=input$std.cphat,
      alpha=input$alpha,
      Npilot=input$Npilot
      )
    powerdist
  })

  pd.sum=function(powerdist, power_desired){
    pd_sum=c(mean(powerdist)
             ,median(powerdist)
             ,mean(powerdist >= power_desired))
    pd_sum=data.frame(PowerDistribution_SummaryStats = pd_sum)
    return(pd_sum)
  }
  
  output$powerdist.sum <- renderDataTable({
    withProgress(message = "Generating Power Distribution", value = 0,
                 detail = NULL, expr = {
      powerdist=pd.reactive()
    })
    
    output$powerdist.plot <- renderPlot({
      hist(powerdist, freq = F
           ,nclass = 50
           , xlim = c(0,1)
           , main = "Power Distribution", xlab = "Power")
      lines(density(powerdist),col=1)
      abline(v=input$power.desired, col=2 )
      abline(v=mean(powerdist), col=3, lty="dashed")
      abline(v=median(powerdist), col=4, lty="dotted")
      legend("topleft"
             , legend = c("desired power", "mean power", "median power")
             , col = c(2,3,4),lty = c("solid","dashed","dotted"))
      # hist( rnorm(100), freq = F,nclass = 50 )
      
    }) # end of renderPlot
    
    data.frame( MeanPower=mean(powerdist)
                ,MedianPower=median(powerdist)
                ,Assurance=mean(powerdist >= input$power.desired) )  
    
  }) # end of renderDataTable

  # power curves ----
  power.curves <- eventReactive(input$update,{
    rangeNplan= as.numeric(unlist(strsplit(input$rangeNplan,",")))
    Nplan_seq= round(seq(rangeNplan[1], rangeNplan[2], by=input$nstep) )
    
    # our method
    pdjointind=lapply(Nplan_seq, pd.jointind,
                        nrawdata=input$nrawdata,
                        a_true=input$std.ahat,
                        b_true=input$std.bhat,
                        cp_true=input$std.cphat,
                        alpha=input$alpha,
                        Npilot=input$Npilot
                        # nrawdata=1000,
                        # a_true=.35,
                        # b_true=.1,
                        # cp_true=0,
                        # alpha=.05,
                        # Npilot=125
                        # ,mc.cores = 1, mc.preschedule = F
                      )
    # conventional
    res_joint.power.ind=lapply(Nplan_seq, joint.power.ind,
                                 a_true=input$std.ahat, 
                                 b_true=input$std.bhat, 
                                 cp_true=input$std.cphat,
                                 alpha=input$alpha
                                 # , mc.cores = 1, mc.preschedule = F
                                 )
    powercp=data.frame( power=c(unlist(res_joint.power.ind) ), 
                        type="joint.power.ind", desired=input$power.desired, nplan=Nplan_seq )
    
    
    # plot power curves ----
    pdjointind_sum=data.frame(
      c(unlist(lapply(pdjointind,mean,na.rm=T)),
        unlist(lapply(pdjointind,median,na.rm=T)),
        unlist(lapply(pdjointind,function(x) mean(x[!is.na(x)]>=input$power.desired ))) # assurance for reaching the desired power or higher.
      ),
      rep(c("mean.power","median.power","assurance"),each=length(pdjointind)),
      rep(c(input$power.desired, input$power.desired, input$assurance.desired),each=length(pdjointind)),
      rep(Nplan_seq,times=3)
    )
    colnames(pdjointind_sum)=c("power","type", "desired","nplan")
    
    resl=rbind(powercp,pdjointind_sum)
    #the specified rangeNplan may not be able to reach the SSP goal
    # recn=aggregate(nplan~type,data=resl[resl$power>=resl$desired,], min )
    # xaxisb=c(rangeNplan[1], recn$nplan, rangeNplan[2])
    
    power_curves= ggplot(resl, aes(x=nplan,y=power, col=type,linetype=type))+
      geom_path( size=1 )+
      geom_hline(yintercept = input$power.desired, colour="black",size=.5, linetype="dashed")+
      geom_hline(yintercept = input$assurance.desired, colour="black",size=.5, linetype="solid")+
      # geom_segment(data=recn1, x=recn1$nplan, xend=recn1$nplan, y=rep(-0.1,nrow(recn1)), yend=recn1$power
      #              # x=recn$nplan, xend=recn$nplan, y=rep(0,nrow(recn)), yend=recn$power
      #              , inherit.aes = F
      #              , linetype=3, size=0.5) +
      scale_y_continuous(name="Power or Assurance", breaks=c(0, 0.4, 0.7, 0.8,  0.9,  1),limits = c(0,1) )+
      scale_x_continuous(name="Planned sample size (Nplan)"
                         # ,breaks = xaxisb
                         , breaks=round(seq(rangeNplan[1], rangeNplan[2], length.out = 10 ))
                         )+
      scale_linetype_manual(name=""
                            , values = c("assurance"="solid", "joint.power.ind"="dashed", "mean.power"="dotted","median.power"="dotdash" )
                            , labels=c("assurance curve for reaching the desired power level or higher","conventional power curve"
                                       ,"mean power curve","median power curve") ) +
      scale_color_manual(name=""
                         ,  values = c("assurance"="red", "joint.power.ind"="black" ,  "mean.power"="blue","median.power"="purple")                          
                         , labels=c("assurance curve for reaching the desired power level or higher","conventional power curve" ,"mean power curve","median power curve" )
                         )+
      theme_bw()+
      theme(
        panel.grid.major = element_line(colour="grey",linetype = 3,size=0.7)
        # panel.grid.major.y = element_blank()
        # panel.grid = element_blank()
        ,axis.title = element_text(size=20)
        ,axis.text = element_text(size=12)
        ,legend.text = element_text(size=22)
        ,legend.title = element_text(size=20)
        ,legend.key.width = unit(16,"mm")
        ,legend.key.height = unit(1.5,"mm")
        ,legend.position = 'bottom'
        ,legend.direction = 'vertical'
      )
    
    power_curves
  })
  
  output$powercurves <- renderPlot({
    power.curves()
  })

}

  
