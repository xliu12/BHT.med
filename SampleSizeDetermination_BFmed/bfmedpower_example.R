source("~/Library/CloudStorage/Box-Box/Labs/Bayesfactor/BFmed_power/SampleSizeDetermination_BFmed/app_bfmedpower_bfabfb.R")
source("~/Library/CloudStorage/Box-Box/Labs/Bayesfactor/BFmed_power/SampleSizeDetermination_BFmed/app_bfmedpower.R")

source("~/Library/CloudStorage/Box-Box/Labs/Bayesfactor/BFmed_power/SampleSizeDetermination_BFmed/bfmedpower_eff.R")
# illustrative example

# standardized path coeffs 
a<-0.39; b<-0.39;cp<-0.14
corrxm <- a
corrmy <- cp*a+b
corrxy <- cp+a*b

# Sample size determination (SSD) for testing mediation with the mediation Bayes factor
res0<- SSD.bfmed(Power.desired=0.8, Nmin=50,Nmax=100,Nstep=2, 
                 a=0.39,b=0.39,cp=0.14, 
                 PriorOdds.a=1, PriorOdds.b=1, PriorOdds.med=NULL,
                 TypeI=0.05, R=1e4,seed=99881)

tab0<- round(res0$outpower.Nseq[ , c(11, 12, 7)],3); 
tab0
write.csv(tab0, "example_SSD_bfmed.csv")
save(res0, file = "res0_bfmedSSD_example_39_39_14_1_1.RData")

# assessing sensitivity of the SSD results to different prior odds specifications
PriorOdds_combs = rbind(
  c(PriorOdds.a = sqrt(.5)/(1-sqrt(.5)), PriorOdds.b = NA, PriorOdds.med = 1),
  c(PriorOdds.a = 3, PriorOdds.b = 1, PriorOdds.med = NA),
  c(PriorOdds.a = 10, PriorOdds.b = 1, PriorOdds.med = NA)
  )
res1 <- SSD.bfmed.sensPriorOdds(Power.desired=0.8, Nmin=50,Nmax=100,Nstep=2, 
                               a=0.39,b=0.39,cp=0.14, 
                               PriorOdds.combs = PriorOdds_combs,
                               TypeI=0.05, R=1e4,seed=13121)
res1_tab <- res1$outpower.PriorOdds.Nseq
write.csv(res1_tab,"res1-sensPriorOdds.csv")

save(res1, file = "res1-sensPriorOdds_bfmedSSD_example_39_39_14_1_1.RData")

# accounting for uncertainty in the effect size specifications 
# res2 = SSD.bfmed.EffectSizeUncertainty(Assurance =0.8, #Expected.Power = 0.8,
#                                 Nmin=50, Nmax=100, Nstep=2, 
#                                 draws_a_b_cp = 100, est.a_b_cp = c(.39, .39, .14), n_previous = 150, #mean.a_b_cp = NULL, cov.a_b_cp = NULL,  
#                                 Power.desired=0.8, PriorOdds.a=1, PriorOdds.b=1, PriorOdds.med=NULL, TypeI=0.05, R=1e4,seed=816051)

res2 <- SSD.bfmed.uncertain.effectsize(Power.desired = 0.8,
                                       Nmin=50, Nmax=100, Nstep=2, 
                                    est.a_b_cp = c(.39, .39, .14), n_pilot = 150, #mean.a_b_cp = NULL, cov.a_b_cp = NULL,  
                                     PriorOdds.a=1, PriorOdds.b=1, PriorOdds.med=NULL, TypeI=0.05, R=1e4,seed=816051)

# colnames(res2$outpower.Nseq)[c(8:10)] <- c('mean_a','mean_b','mean_cp')
# colnames(res2$out_SSD$outpower)[c(8:10)]<- c('mean_a','mean_b','mean_cp')
# res2$outpower.Nseq[ ,c(8:10)]<- t(replicate(nrow(res2$outpower.Nseq),c(.39, .39, .14)))
# res2_tab<-res2$outpower.Nseq

save.image("~/Library/CloudStorage/Box-Box/Labs/Bayesfactor/BFmed_power/SampleSizeDetermination_BFmed/res2-bfmedpower_eff_example_39_39_14_1_1.RData")
res2

write.csv( res2, "res2-uncertain.effectsize.csv")




