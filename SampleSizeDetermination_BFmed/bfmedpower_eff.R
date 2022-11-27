library(mvtnorm)
library(parallel)
library(BayesFactor)
#################### uncertainty in the effect size specifications ############################

expected.power.bfmed <- function( est.a_b_cp = c(.39, .39, 0.14), n_previous = 150, mean.a_b_cp = NULL, cov.a_b_cp = NULL,  
                                 N=100,PriorOdds.a=1, PriorOdds.b=1, PriorOdds.med=NULL, TypeI=0.05,R=1e4,seed=816051){
  PriorOdds.a=PriorOdds.a; PriorOdds.b=PriorOdds.b; PriorOdds.med=PriorOdds.med;TypeI=TypeI; R=R;seed=seed;
  
  if( (!is.null(cov.a_b_cp)) & (!is.null(mean.a_b_cp)) ){
    mean_a_b_cp = mean.a_b_cp
    cov_a_b_cp = cov.a_b_cp
  }
  # the sampling distribution of the standardized effect size estimates
  hat_a = est.a_b_cp[1]; hat_b = est.a_b_cp[2]; hat_cp = est.a_b_cp[3]
  var_em=1-hat_a^2; var_ey=1-hat_cp^2-hat_b^2-2*hat_a*hat_b*hat_cp
  var_x = 1; cov_xm = matrix(c(1, hat_a, hat_a, 1) ,2,2)
  # hat_a \independent (hat_b, hat_cp) approximately (Wang, 2018)
  var_hat_a = (1/n_previous)*var_em/var_x
  cov_hat_bcp = (1/n_previous)*var_ey*solve(cov_xm)
  cov_hat_a_bcp = cbind(c(var_hat_a, 0, 0), rbind(c(0,0),cov_hat_bcp))
  
  mean_a_b_cp = est.a_b_cp; cov_a_b_cp = cov_hat_a_bcp
  
  # conditional distribution of cp given b=0
  mean_cp_givenb0 = mean_a_b_cp[3] + (cov_a_b_cp[2,3]/cov_a_b_cp[2,2])*( 0 - mean_a_b_cp[2] )
  var_cp_givenb0 = cov_a_b_cp[3,3] - (cov_a_b_cp[2,3])^2/cov_a_b_cp[2,2]
    
  ## prior odds specification 
  if( (!is.null(PriorOdds.a)) & (!is.null(PriorOdds.b))){ 
    # PriorOdds.med = PriorOdds.b*PriorOdds.a / (1+PriorOdds.b+PriorOdds.a) ## using the independent paths assumption
    PriorOdds.med = 1 / (1/(PriorOdds.a*PriorOdds.b) + 1/PriorOdds.a+ 1/PriorOdds.b) ## applicable when PriorOdds.a=Inf
  }
  if( (!is.null(PriorOdds.a)) & (!is.null(PriorOdds.med))){ 
    PriorOdds.b = PriorOdds.med*(1/PriorOdds.a+ 1) / (1 - PriorOdds.med/PriorOdds.a)
  }
  if( (!is.null(PriorOdds.b)) & (!is.null(PriorOdds.med))){ 
    PriorOdds.a = PriorOdds.med*(1/PriorOdds.b+ 1) / (1-PriorOdds.med/PriorOdds.b)
  }
  ## for any priorodds
  q10 = (1-PriorOdds.med/PriorOdds.a )/(1/PriorOdds.a+1) ##applicable when PriorOdds.a=Inf
  q01 = (1-PriorOdds.med/PriorOdds.b)/(1/PriorOdds.b +1)
  q00 = 1-q10-q01
  
  
  # draw from the sampling distribution
  Nsim_H0=Nsim_H1=R # =1e4
  seed_H0=seed #=816051; 
  seed_H1=seed*2
  
  ## number of simulated-datasets under H0
  # Nsim_H0 = 1e4 ## Fu et al. (2021) recommend at least 1e4
  set.seed(seed=seed_H0); seeds_OneData.H0 <- sample(1:1e8, Nsim_H0)
  whichq_H0 = rmultinom(n=Nsim_H0, 1, prob = c(q10,q01,q00)) ## for each data, only one model can be true
  
  OneData.H0 <- function(i=1, N , mean_a_b_cp, cov_a_b_cp){
    set.seed(seed = seeds_OneData.H0[i])
    # a_b_cp = rmvnorm(1, mean = mean_a_b_cp, sigma = cov_a_b_cp)
    # a = a_b_cp[1]; b=a_b_cp[2]; cp=a_b_cp[3]
    # ## error variances in standardized simple mediation model
    # var_em=1-a^2; var_ey=1-cp^2-b^2-2*a*b*cp
    
    seedxmy=3*seeds_OneData.H0[i]; set.seed(seed = seedxmy)
    X=rnorm(N)
    whichq.i=c('q10','q01','q00') [which(whichq_H0[,i]==1)]
    if(whichq.i == 'q01') { ## a=0, b=b
      a=0; 
      b_cp = rmvnorm(1, mean_a_b_cp[2:3], sigma = cov_a_b_cp[2:3,2:3])
      b = b_cp[1]; cp = b_cp[2]
      var_em=1-a^2; var_ey=1-cp^2-b^2-2*a*b*cp
      M=0*X+rnorm(N, 0, sqrt(var_em)) 
      Y=cp*X+b*M+rnorm(N,0, sqrt(var_ey))
    }
    if(whichq.i == 'q10') { ## a=a, b=0
      a=rnorm(1, mean_a_b_cp[1], sd = sqrt(cov_a_b_cp[1,1]) )
      b=0
      cp=rnorm(1, mean_cp_givenb0, sd = sqrt( var_cp_givenb0 ) )
      var_em=1-a^2; var_ey=1-cp^2-b^2-2*a*b*cp
      M=a*X+rnorm(N, 0, sqrt(var_em)) 
      Y=cp*X+0*M+rnorm(N,0, sqrt(var_ey))
    }
    if(whichq.i == 'q00') { ## a=0, b=0
      a=0
      b=0
      cp=rnorm(1, mean_cp_givenb0, sd = sqrt( var_cp_givenb0 ) )
      var_em=1-a^2; var_ey=1-cp^2-b^2-2*a*b*cp
      M=0*X+rnorm(N, 0, sqrt(var_em)) 
      Y=cp*X+0*M+rnorm(N,0, sqrt(var_ey))
    }
    
    ## Bayes factors of paths a and b
    dat1=data.frame(X,M,Y)
    BFa = extractBF( lmBF(M~X, data = dat1) )$bf 
    BFb = extractBF(lmBF(Y~M+X, data = dat1))$bf / extractBF(lmBF(Y~X, data = dat1))$bf
    
    res=c(BFa,BFb)
    names(res)=c('BFa','BFb')
    
    out=res
    return(out)
  }
  
  ptm0=proc.time()[3]
  bfabfb.Nsim_H0= mclapply(1:Nsim_H0, OneData.H0,
                           N=N,mean_a_b_cp=mean_a_b_cp, cov_a_b_cp=cov_a_b_cp,
                           mc.cores = 5,
                           mc.preschedule = F)
  ptm1=proc.time()[3]; #print(ptm1-ptm0)
  
  
  # Nsim_H1=1e4 ### recommend at least 1e4
  set.seed(seed=seed_H1); seeds_OneData.H1 <- sample(1:1e8, Nsim_H1)
  OneData.H1 <- function(i=1, N, mean_a_b_cp, cov_a_b_cp){
    set.seed(seed = seeds_OneData.H1[i])
    a_b_cp = rmvnorm(1, mean = mean_a_b_cp, sigma = cov_a_b_cp)
    a = a_b_cp[1]; b=a_b_cp[2]; cp=a_b_cp[3]
    ## error variances in standardized simple mediation model
    var_em=1-a^2; var_ey=1-cp^2-b^2-2*a*b*cp
   
    seedxmy=3*seeds_OneData.H1[i]; set.seed(seed = seedxmy)
    X=rnorm(N)
    M=a*X+rnorm(N, 0, sqrt(var_em)) 
    Y=cp*X+b*M+rnorm(N,0, sqrt(var_ey))
    
    ## Bayes factors of paths a and b
    dat1=data.frame(X,M,Y)
    BFa = extractBF( lmBF(M~X, data = dat1) )$bf 
    BFb = extractBF(lmBF(Y~M+X, data = dat1))$bf / extractBF(lmBF(Y~X, data = dat1))$bf
    
    res=c(BFa,BFb)
    names(res)=c('BFa','BFb')
    
    out=res
    return(out)
  }
  
  ptm0=proc.time()[3]
  bfabfb.Nsim_H1= mclapply(1:Nsim_H1, OneData.H1,
                           N=N, mean_a_b_cp=mean_a_b_cp, cov_a_b_cp=cov_a_b_cp,
                           mc.cores = 5,
                           mc.preschedule = F)
  ptm1=proc.time()[3]; #print(ptm1-ptm0) ## 90.582 seconds for 1e2 datasets N=200
  
  ## combine the BFa BFb from H0 and H1
  bfabfb_Nsim_H0 = data.frame(do.call(rbind, bfabfb.Nsim_H0), trueH='H0')
  bfabfb_Nsim_H1 = data.frame(do.call(rbind, bfabfb.Nsim_H1), trueH='H1')
  bfabfb_Nsim_H0H1=rbind(bfabfb_Nsim_H0, bfabfb_Nsim_H1) 
  
  ##### calculate BFmed using BFa, BFb ######
  bfmed_Nsim_H0H1 = bfabfb_Nsim_H0H1
  bfmed_Nsim_H0H1$BFmed = bfabfb_Nsim_H0H1$BFa*bfabfb_Nsim_H0H1$BFb / (q00 + bfabfb_Nsim_H0H1$BFb*q01 + bfabfb_Nsim_H0H1$BFa*q10)
  bfmed_Nsim_H0H1$POmed = bfmed_Nsim_H0H1$BFmed * PriorOdds.med
  
  # TypeI = c(.05, 0.1, 0.01) ## this can be modified after obtaining bfmed_Nsim_H0H1
  cutoff_bfmed=quantile(bfmed_Nsim_H0H1$BFmed[bfmed_Nsim_H0H1$trueH=='H0'], probs=c(1-TypeI))
  power_bfmed=sapply(cutoff_bfmed, function(x){mean( bfmed_Nsim_H0H1$BFmed[bfmed_Nsim_H0H1$trueH=='H1'] > x )})
  names(cutoff_bfmed)=paste0(c('cutoff_bfmed_'),TypeI)
  names(power_bfmed)=paste0(c('power_'),TypeI )
  
  ## output 
  priorspeci=c(PriorOdds.a=PriorOdds.a, PriorOdds.b=PriorOdds.b, PriorOdds.med=PriorOdds.med, q10=q10,q01=q01,q00=q00) 
  designspeci=c(N=N,mean_a_b_cp)
  names(designspeci)=c('N','mean_a','mean_b','mean_cp')
  
  outpower=cbind(t(priorspeci),t(designspeci), t(cutoff_bfmed),t(power_bfmed))
  outpower=data.frame(outpower)
  # colnames(outpower)[c(11:16)]=c(paste0(c('cutoff_bfmed_'),TypeI),paste0(c('power_'),TypeI )) 
  
  out=list(outpower=outpower, bfmed_Nsim_H0H1=bfmed_Nsim_H0H1)
  return(out)
}



SSD.bfmed.uncertain.effectsize <- function(Power.desired=0.8, Nmin=50,Nmax=500,Nstep=10,
                                           est.a_b_cp = c(.39, .39, .39), n_pilot = 150,mean.a_b_cp = NULL, cov.a_b_cp = NULL,  
                                           PriorOdds.a=1, PriorOdds.b=1, PriorOdds.med=NULL,TypeI=0.05, R=1e4,seed=816051) {

  Nseq = seq(Nmin, Nmax, by=Nstep)
  outpower.Nseq=NULL
  for(N in Nseq){
    oneN = expected.power.bfmed( est.a_b_cp = est.a_b_cp, n_previous = n_pilot, mean.a_b_cp = mean.a_b_cp, cov.a_b_cp = cov.a_b_cp,  
                               N=N,PriorOdds.a=PriorOdds.a, PriorOdds.b=PriorOdds.b, PriorOdds.med=PriorOdds.med, TypeI=TypeI,
                                             R=R,seed=seed )
    
    print(  oneN$outpower[, paste0(c('power_'),TypeI )] )
    outpower.Nseq = rbind(outpower.Nseq, oneN$outpower)
    
    if (  oneN$outpower[, paste0(c('power_'),TypeI )] >=Power.desired){
      out_SSD=oneN
      break
    }
  }
  # out = list(outpower.Nseq=outpower.Nseq, out_SSD=out_SSD)
  out = outpower.Nseq
  return(out)
}


# power.bfmed.EffectSizeUncertainty <- function( draws_a_b_cp = 500, est.a_b_cp = c(.39, .39, .39), n_previous = 150, mean.a_b_cp = NULL, cov.a_b_cp = NULL,  
#                                                Power.desired=0.8, N=100,PriorOdds.a=NULL, PriorOdds.b=NULL, PriorOdds.med=NULL, TypeI=0.05,R=1e4,seed=816051){
#   PriorOdds.a=PriorOdds.a; PriorOdds.b=PriorOdds.b; PriorOdds.med=PriorOdds.med;TypeI=TypeI; R=R;seed=seed;
#   
#   if( (!is.null(cov.a_b_cp)) & (!is.null(mean.a_b_cp)) ){
#     set.seed(seed = seed)
#     a_b_cp = rmvnorm(draws_a_b_cp, mean = mean.a_b_cp, sigma = cov.a_b_cp)
#   }
#   # the sampling distribution of the standardized effect size estimates
#   hat_a = est.a_b_cp[1]; hat_b = est.a_b_cp[2]; hat_cp = est.a_b_cp[3]
#   var_em=1-hat_a^2; var_ey=1-hat_cp^2-hat_b^2-2*hat_a*hat_b*hat_cp
#   var_x = 1; cov_xm = matrix(c(1, hat_a, hat_a, 1) ,2,2)
#   # hat_a \independent (hat_b, hat_cp) approximately (Wang, 2018)
#   var_hat_a = (1/n_previous)*var_em/var_x
#   cov_hat_bcp = (1/n_previous)*solve(cov_xm)
#   cov_hat_a_bcp = cbind(c(var_hat_a, 0, 0), rbind(c(0,0),cov_hat_bcp))
#   
#   # draw from the sampling distribution
#   set.seed(seed = seed)
#   a_b_cp = rmvnorm(draws_a_b_cp, mean = est.a_b_cp, sigma = cov_hat_a_bcp)
#   
#   set.seed(seed = seed); seeds_abcp = sample(1:9e3, nrow(a_b_cp))
#   outpower_draws_abcp = NULL
#   for(l in 1:nrow(a_b_cp)){
#     a=a_b_cp[l, 1]; b=a_b_cp[l,2]; cp=a_b_cp[l,3]
#     onedraw = bfmed.power(N = N, a=a,b=b,cp = cp, 
#                 PriorOdds.a = PriorOdds.a, PriorOdds.b = PriorOdds.b,PriorOdds.med = PriorOdds.med,
#                 TypeI = TypeI, R = R, seed=seeds_abcp[l]  )
#     outpower_draws_abcp = rbind(outpower_draws_abcp, onedraw$outpower)
#   }
#   outpower_draws_abcp =as.data.frame(outpower_draws_abcp)
#   
#   out = list(expected_power = mean(outpower_draws_abcp[, paste0(c('power_'),TypeI )] ), 
#        assurance=mean(outpower_draws_abcp[, paste0(c('power_'),TypeI )] >= Power.desired),
#        outpower_draws_abcp = outpower_draws_abcp )
#   
#   return(out)
# }


# SSD.bfmed.EffectSizeUncertainty <- function( Assurance =NULL, Expected.Power = 0.8, 
#                                              Nmin=70, Nmax=200,Nstep=10, 
#                                              draws_a_b_cp = 500, est.a_b_cp = c(.39, .39, .39), n_previous = 150,mean.a_b_cp = NULL, cov.a_b_cp = NULL,  
#                                              Power.desired=0.8, PriorOdds.a=NULL, PriorOdds.b=NULL, PriorOdds.med=NULL, TypeI=0.05,R=1e4,seed=816051 ){
#   Nseq = seq(Nmin, Nmax, by=Nstep)
#   SSD_EffectSizeUncertainty = NULL
#   for(N in Nseq){
#     oneN = power.bfmed.EffectSizeUncertainty(draws_a_b_cp = draws_a_b_cp, est.a_b_cp = est.a_b_cp, n_previous = n_previous, mean.a_b_cp = mean.a_b_cp, cov.a_b_cp = cov.a_b_cp,  
#                                              Power.desired=Power.desired, N=N,PriorOdds.a=PriorOdds.a, PriorOdds.b=PriorOdds.b, PriorOdds.med=PriorOdds.med, TypeI=TypeI,
#                                              R=R,seed=seed )
#     oneNres = unlist(c(oneN$outpower_draws_abcp[1,c(1:3,7)], assurance = oneN$assurance, expected_power = oneN$expected_power) )
#     print(oneNres)
#     SSD_EffectSizeUncertainty = rbind(SSD_EffectSizeUncertainty, oneNres)
#     
#     if(!is.null(Assurance)){
#       if( oneN$assurance >= Assurance ){ break }
#     }
#     if(is.null(Assurance)){
#       if( oneN$expected_power >= Expected.Power ){ break }
#     }
#     
#   }
#   return( SSD_EffectSizeUncertainty )
# }




