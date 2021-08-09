

# independent path scenario
# explore the impact of prior odds specifications

library(tidyverse)
# source('~/Google Drive/Labs/JohnnyLab/BFmed/fun_NHBTmed.R')

# specify PriorOdds.a, PriorOdds.b ----

plot.ab <- function(BF.a=1, BF.b=1){
  vprior_odds_pathb=c(1/3,1,3)
  vprior_odds_a=c(1/100,1/30,1/10,1/3,1,3,10,30,100)
  vprior_prob_a=vprior_odds_a/(vprior_odds_a+1)
  
  pathb_a=as.data.frame(
    expand.grid(BFa=BF.a,BFb=BF.b
                ,prior_odds_pathb=vprior_odds_pathb
                ,prior_odds_a=vprior_odds_a
    ))
  out=mapply(FUN = pathb.a, 
             PriorOdds.b=pathb_a$prior_odds_pathb,
             PriorOdds.a=pathb_a$prior_odds_a,
             BF.a=pathb_a$BFa, BF.b=pathb_a$BFb,
             SIMPLIFY = T)
  pathb_a1=cbind(pathb_a,t(out))
  
  pathb_a1=as_tibble(pathb_a1) %>%
    mutate( prior_prob_a=prior_odds_a/(1+prior_odds_a) ) %>%
    mutate( PriorOdds.b=round(prior_odds_pathb,2) )
  
  lpathb_a1=gather(pathb_a1, key = "BHT",value = "value",
                   BF.med, PosteriorOdds.med) %>%
    mutate( BF.a_BF.b = paste0('BF.a=',round(BFa,2)," BF.b=",round(BFb,2)) )
  
  # plot 
  plot_ab = ggplot(lpathb_a1) +
    geom_line( aes(x=prior_prob_a, y=value )
               , size=1 ) +
    scale_y_continuous('BHT.med'
                       # ,breaks = c(0,0.33,1,3,4,6,9,10)
    ) +
    scale_x_continuous("PriorOdds.a", 
                       limits = c(0, 1)
                       ,breaks = vprior_prob_a[2:8]
                       ,labels = as.character(round(vprior_odds_a[2:8],2))
    ) +
    facet_grid( PriorOdds.b ~ BHT, scales = "free_y"
                ,labeller = labeller(PriorOdds.b=label_both, BHT=label_value) ) +
    theme_bw()  + 
    theme( axis.title = element_text( size = 7 )
           , axis.text = element_text( size = 10 )
           , legend.text = element_text( size = 7 )
           , legend.title = element_text( size = 7 )
           , legend.position = "right"
           , legend.box = "vertical"
           , legend.direction = "vertical"
           , legend.key.width = unit(1.5,"cm")
           , legend.key = element_rect(fill = NULL)
           , strip.text = element_text(size = 7)
           , panel.grid.minor.x = element_blank()
           , panel.grid.minor.y = element_blank()
           , panel.grid.major.x = element_blank()
           # , panel.spacing = unit(4,"mm")
    ) 
  
  plot_ab
  
  return(plot_ab)
}


# specify PriorOdds.med,PriorOdds.a ----

plot.meda <- function(BF.a=1, BF.b=1){
  vprior_odds_med=c(1/3,1,3)
  vprior_odds_a=c(1/100,1/30,1/10,1/3,1,3,10,30,100)
  vprior_prob_a=vprior_odds_a/(vprior_odds_a+1)
  
  med_a=as.data.frame(
    expand.grid(BFa=BF.a,BFb=BF.b
                ,prior_odds_med=vprior_odds_med
                ,prior_odds_a=vprior_odds_a
    ))
  out=mapply(FUN = med.a, 
             PriorOdds.med=med_a$prior_odds_med,
             PriorOdds.a=med_a$prior_odds_a,
             BF.a=med_a$BFa, BF.b=med_a$BFb,
             SIMPLIFY = T)
  med_a1=cbind(med_a,t(out))
  
  med_a1=as_tibble(med_a1) %>%
    mutate(
      BF_med=ifelse(prior_odds_a<prior_odds_med, NA, BF.med) ) %>%
    mutate( 
      BF_med=ifelse(prior_odds_a==prior_odds_med, BFa, BF_med) ) %>%
    mutate(
      PosteriorOdds_med=ifelse(prior_odds_a<prior_odds_med, NA, PosteriorOdds.med) ) %>%
    mutate( 
      PosteriorOdds_med=ifelse(prior_odds_a==prior_odds_med, BFa*prior_odds_a, PosteriorOdds_med) ) %>%
    mutate( BF.med=BF_med) %>%
    mutate( PosteriorOdds.med=PosteriorOdds_med) %>%
    mutate( prior_prob_a=prior_odds_a/(1+prior_odds_a) ) %>%
    mutate( prior_prob_med=prior_odds_med/(1+prior_odds_med) ) %>%
    mutate( PriorOdds.med=round(prior_odds_med,2) )
  
  lmed_a1=gather(med_a1, key = "BHT",value = "value",
                 BF.med, PosteriorOdds.med) %>%
    mutate( BF.a_BF.b = paste0('BF.a=',round(BFa,2)," BF.b=",round(BFb,2)) )
  
  # plot 
  
  plot_meda = ggplot(lmed_a1) +
    geom_line( aes(x=prior_prob_a, y=value)
               , size=1 ) +
    # scale_linetype_manual('',values = lineType) +
    # scale_color_manual('',values = coLor) +
    scale_y_continuous('BHT.med'
                       # ,breaks = c(0,0.33,1,3,4,6,9,10)
    ) +
    scale_x_continuous("PriorOdds.a", 
                       limits = c(0, 1)
                       ,breaks = vprior_prob_a[2:8]
                       ,labels = as.character(round(vprior_odds_a[2:8],2))
    ) +
    facet_grid( PriorOdds.med ~ BHT, scales = "free_y"
                ,labeller = labeller(PriorOdds.med=label_both, BHT=label_value) ) +
    theme_bw()  + 
    theme( axis.title = element_text( size = 7 )
           , axis.text = element_text( size = 10 )
           , legend.text = element_text( size = 7 )
           , legend.title = element_text( size = 7 )
           , legend.position = "right"
           , legend.box = "vertical"
           , legend.direction = "vertical"
           , legend.key.width = unit(1.5,"cm")
           , legend.key = element_rect(fill = NULL)
           , strip.text = element_text(size = 7)
           , panel.grid.minor = element_blank()
           , panel.grid.major.x = element_blank()
           # , panel.spacing = unit(4,"mm")
    ) 
  plot_meda
  return(plot_meda)
}


# specify PriorOdds.med ----
# specify q10 

plot.medq <- function(BF.a=1, BF.b=1){
  vprior_odds_med=c(1/3,1,3)
  vq10=seq(0,1,length.out = 9)
  
  med_q10=as.data.frame(
    expand.grid(BFa=BF.a,BFb=BF.b
                ,prior_odds_med=vprior_odds_med
                ,H0.q10=vq10
    ))
  out=mapply(FUN = med.q10, 
             PriorOdds.med=med_q10$prior_odds_med,
             q10=med_q10$H0.q10,
             BF.a=med_q10$BFa, BF.b=med_q10$BFb,
             SIMPLIFY = T)
  med_q=cbind(med_q10,t(out))
  
  med_q=as_tibble(med_q) %>%
    mutate( PriorOdds.med=round(prior_odds_med,2) )
  
  lmed_q=gather(med_q, key = "BHT",value = "value",
                BF.med, PosteriorOdds.med) %>%
    mutate( BF.a_BF.b = paste0('BF.a=',round(BFa,2)," BF.b=",round(BFb,2)) )
  # plot 
  
  plot_medq = ggplot(lmed_q) +
    geom_line( aes(x=H0.q10, y=value )
               , size=1 ) +
    # scale_linetype_manual('',values = lineType) +
    # scale_color_manual('',values = coLor) +
    scale_y_continuous('BHT.med'
                       # ,breaks = c(0,0.33,1,3,4,6,9,10)
    ) +
    scale_x_continuous("H0.q10", 
                       limits = c(0, 1)
                       # ,breaks = vprior_prob_a[2:8]
                       # ,labels = as.character(round(vprior_odds_a[2:8],2))
    ) +
    facet_grid( PriorOdds.med ~ BHT, scales = "free_y"
                ,labeller = labeller(PriorOdds.med=label_both, BHT=label_value) ) +
    theme_bw()  + 
    theme( axis.title = element_text( size = 7 )
           , axis.text = element_text( size = 10 )
           , legend.text = element_text( size = 7 )
           , legend.title = element_text( size = 7 )
           , legend.position = "right"
           , legend.box = "vertical"
           , legend.direction = "vertical"
           , legend.key.width = unit(1.5,"cm")
           , legend.key = element_rect(fill = NULL)
           , strip.text = element_text(size = 7)
           , panel.grid.minor = element_blank()
           , panel.grid.major.x = element_blank()
           # , panel.spacing = unit(4,"mm")
    ) 
  
  plot_medq
  return(plot_medq)
}





















