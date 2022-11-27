
library(BayesFactor)
library(mvtnorm)
library(parallel)

## H0: a*b=0 is true 
##find cut-off using Type I error

## the cutoff for BFmed depends on how the researcher specified H0, i.e., characterized by q00,q10,q01, the conditional prior probs of the three null scenarios under H0 
## Independent paths, BFmed is a function of prior odds/probs and BFa, BFb
## so, only need to emprically approximate the joint sampling distribution of (BFa, BFb), using a large number of simulated datasets from the data generation model
## the data generation model under H0 is a mixture of the 3 null models -- weights given by the qs.
## the data generation model under H1 is just the specified effect sizes

## prior odds specified by the researcher
# PriorOdds.a=1; PriorOdds.b=1 # independent paths
# ## effect sizes in standardized simple mediation model 
# a=.39; b=.39; cp=0
# N=Nplan=100

power.bfmed_bfabfb_pools <- function(N=100, a=0.39,b=0.39,cp=0, PriorOdds.a=NULL, PriorOdds.b=NULL, PriorOdds.med=NULL, TypeI=0.05,R=1e4,seed=816051){
  Nsim_H0=Nsim_H1=R # =1e4
  ## error variances in standardized simple mediation model
  var_em=1-a^2; var_ey=1-cp^2-b^2-2*a*b*cp
  
  bfabfb.samplingdist <- function(N, a1,b1,cp, R=1e4, seed_onecond=816051 ){
    set.seed(seed = seed_onecond)
    seeds_OneData.cond = sample(1:1e8, R)
    
    OneData.cond <- function(i=1, N, a1, b1, cp){
      set.seed(seed = seeds_OneData.cond[i])
      X=rnorm(N)
      M=a1*X+rnorm(N,0,sqrt(var_em)) 
      Y=cp*X+b1*M+rnorm(N,0,sqrt(var_ey))
      
      ## Bayes factors of paths a and b
      dat1=data.frame(X,M,Y)
      BFa = extractBF( lmBF(M~X, data = dat1) )$bf 
      BFb = extractBF(lmBF(Y~M+X, data = dat1))$bf / extractBF(lmBF(Y~X, data = dat1))$bf
      
      res=c(BFa,BFb)
      names(res)=c('BFa','BFb')
      
      out=res
      return(out)
    }
    
    bfabfb.Nsim_cond= mclapply(1:R, OneData.cond,
                               N=N,a1=a1,b1=b1, cp=cp,
                               mc.cores = 5,
                               mc.preschedule = TRUE)
    
    ## output 
    # priorspeci=c(PriorOdds.a=PriorOdds.a, PriorOdds.b=PriorOdds.b, PriorOdds.med=PriorOdds.med, q10=q10,q01=q01,q00=q00) 
    designspeci=c(N=N,a=a1,b=b1,cp=cp)
    bfabfb_out=cbind( t(designspeci),data.frame(do.call(rbind, bfabfb.Nsim_cond) ) )
    out=data.frame(bfabfb_out)
    
    return(out)
  }
  
  ## simulate a pool of BFa, BFb (for given N, and effect sizes) to sample from 
  ##construct the BFmed under H0 by drawing (BFa,BFb) from p(BFa,BFb | each of the null scenarios) according to the conditional prior probabilities of the three scenarios under H0 
  ##this pool of (BFa,BFb) can be used for obtaining BFmed distributions under varied prior odds specifications, so as to assess the sensitivity to the prior odds
  set.seed(seed = seed); abconds=data.frame(expand.grid(a1=c(a, 0),b1=c(b,0)))
  seeds_abconds = sample(1:4e4, nrow(abconds))
  
  bfabfb_pools=NULL
  
  for(i in 1:nrow(abconds)){
    bfabfb_out=bfabfb.samplingdist(N = N,a1 = abconds$a1[i],b1 = abconds$b1[i],cp = cp, R = R,seed_onecond = seeds_abconds[i])
    bfabfb_pools=rbind(bfabfb_pools, bfabfb_out)
  }
  
  ## bfmedpower.givenqs 
  ## prior odds specification 
  
  if( (!is.null(PriorOdds.a)) & (!is.null(PriorOdds.b)) ){ 
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
  
  # bfpower_givenqs <- function(Nplan, asize, bsize,cpsize, Nsim_H0=1e4, Nsim_H1=1e4){
  Nsim_H0=Nsim_H1=R ; Nplan=N; asize=a; bsize=b; cpsize=cp
  set.seed(seed = seed); seeds_sampleqs = sample(3e3, 3)
  ## sample under H0
  aa = bfabfb_pools[ (bfabfb_pools$N==Nplan & bfabfb_pools$a==0 & bfabfb_pools$b==0 & bfabfb_pools$cp==cpsize), ]
  sample_q00 =NULL; set.seed(seed = seeds_sampleqs[1])
  if(nrow(aa)>0){sample_q00 = aa[ sample.int(nrow(aa), round(Nsim_H0*q00), replace = F), ] }
  
  aa = bfabfb_pools[ (bfabfb_pools$N==Nplan & bfabfb_pools$a==asize & bfabfb_pools$b==0 & bfabfb_pools$cp==cpsize), ]
  sample_q10=NULL;   set.seed(seed = seeds_sampleqs[2])
  if(nrow(aa)>0) {sample_q10=aa[ sample.int(nrow(aa), round(Nsim_H0*q10), replace = F), ]  }
  
  aa = bfabfb_pools[ (bfabfb_pools$N==Nplan & bfabfb_pools$a==0 & bfabfb_pools$b==bsize & bfabfb_pools$cp==cpsize), ]
  sample_q01=NULL;  set.seed(seed = seeds_sampleqs[3])
  if(nrow(aa)>0) {sample_q01=aa[ sample.int(nrow(aa), round(Nsim_H0*q01), replace = F), ]  }
  
  sample_H0 = rbind(rbind(sample_q00, sample_q10), sample_q01)
  
  ## sample under H1 # no sampling uncertainty, used all of (bfa, bfb) under H1
  aa = bfabfb_pools[ (bfabfb_pools$N==Nplan & bfabfb_pools$a==asize & bfabfb_pools$b==bsize & bfabfb_pools$cp==cpsize), ]
  sample_H1=NULL;  
  if(nrow(aa)>0) {sample_H1=aa[ sample.int(nrow(aa), round(Nsim_H1*1), replace = F), ]   }
  
  ## combine the BFa BFb from H0 and H1
  bfabfb_Nsim_H0 = data.frame(sample_H0, trueH='H0')
  bfabfb_Nsim_H1 = data.frame(sample_H1, trueH='H1')
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
  designspeci=c(N=N,a=a,b=b,cp=cp)
  
  outpower=cbind(t(priorspeci),t(designspeci), t(cutoff_bfmed),t(power_bfmed))
  outpower=data.frame(outpower)
  # colnames(outpower)[c(11:16)]=c(paste0(c('cutoff_bfmed_'),TypeI),paste0(c('power_'),TypeI )) 
  
  out=list(outpower=outpower,bfabfb_pools=bfabfb_pools, bfmed_Nsim_H0H1=bfmed_Nsim_H0H1)
  return(out)
}

## after obtainning one power with one set of prior odds specification, sensitivity analysis to assess how power changes with alternative prior odds specifications 
## with the same pools of (BFa, BFb), draw under the same design conditions
#bfabfb_pools=out$bfabfb_pools

power.bfmed.sensPriorOdds <- function(bfabfb_pools,N,a,b,cp, PriorOdds.a=1, PriorOdds.b=1, PriorOdds.med=NULL, TypeI=0.05,R=1e4,seed=1111){
  
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
  
  # bfpower_givenqs <- function(Nplan, asize, bsize,cpsize, Nsim_H0=1e4, Nsim_H1=1e4){
  Nsim_H0=Nsim_H1=R ; Nplan=N; asize=a; bsize=b; cpsize=cp
  set.seed(seed = seed); seeds_sampleqs = sample(3e3, 3)
  
  ## sample under H0
  aa = bfabfb_pools[ (bfabfb_pools$N==Nplan & bfabfb_pools$a==0 & bfabfb_pools$b==0 & bfabfb_pools$cp==cpsize), ]
  sample_q00 =NULL; set.seed(seed = seeds_sampleqs[1]) 
  if(nrow(aa)>0){sample_q00 = aa[ sample.int(nrow(aa), round(Nsim_H0*q00), replace = F), ] }
  
  aa = bfabfb_pools[ (bfabfb_pools$N==Nplan & bfabfb_pools$a==asize & bfabfb_pools$b==0 & bfabfb_pools$cp==cpsize), ]
  sample_q10=NULL; set.seed(seed = seeds_sampleqs[2])
  if(nrow(aa)>0) {sample_q10=aa[ sample.int(nrow(aa), round(Nsim_H0*q10), replace = F), ]  }
  
  aa = bfabfb_pools[ (bfabfb_pools$N==Nplan & bfabfb_pools$a==0 & bfabfb_pools$b==bsize & bfabfb_pools$cp==cpsize), ]
  sample_q01=NULL;  set.seed(seed = seeds_sampleqs[3])
  if(nrow(aa)>0) {sample_q01=aa[ sample.int(nrow(aa), round(Nsim_H0*q01), replace = F), ]  }
  
  sample_H0 = rbind(rbind(sample_q00, sample_q10), sample_q01)
  
  ## sample under H1
  aa = bfabfb_pools[ (bfabfb_pools$N==Nplan & bfabfb_pools$a==asize & bfabfb_pools$b==bsize & bfabfb_pools$cp==cpsize), ]
  sample_H1=NULL;  
  if(nrow(aa)>0) {sample_H1=aa[ sample.int(nrow(aa), round(Nsim_H1*1), replace = F), ]   }
  
  ## combine the BFa BFb from H0 and H1
  bfabfb_Nsim_H0 = data.frame(sample_H0, trueH='H0')
  bfabfb_Nsim_H1 = data.frame(sample_H1, trueH='H1')
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
  designspeci=c(N=N,a=a,b=b,cp=cp)
  
  outpower=cbind(t(priorspeci),t(designspeci), t(cutoff_bfmed),t(power_bfmed))
  outpower=data.frame(outpower)
  
  return(outpower)
  
}


SSD.bfmed.sensPriorOdds <- function(Power.desired=0.8, Nmin=70,Nmax=100,Nstep=10, a=0.39,b=0.39,cp=0, 
                                    PriorOdds.combs = data.frame(PriorOdds.a = c(sqrt(.5)/(1-sqrt(.5)),3, 10), PriorOdds.b=c(NA, 1,1 ), PriorOdds.med=c(1, NA,NA)),
                                    TypeI=0.05, R=1e4,seed=816051) {
  a=a;b=b;cp=cp;TypeI=TypeI; R=R;seed=seed; Nseq = seq(Nmin, Nmax, by=Nstep)
  outpower.PriorOdds.Nseq=NULL
  
  for(N in Nseq){
    set.seed(seed = seed); seeds_sens = sample(1:1e4, nrow(PriorOdds.combs))
    
    PriorOdds.a=PriorOdds.combs[1,1]; PriorOdds.b=PriorOdds.combs[1,2]; PriorOdds.med=PriorOdds.combs[1,3];
    if(is.na(PriorOdds.a)){PriorOdds.a=NULL}; if(is.na(PriorOdds.b)){PriorOdds.b=NULL}; if(is.na(PriorOdds.med)){PriorOdds.med=NULL}
    
    out1 = power.bfmed_bfabfb_pools(N = N, a=a,b=b,cp = cp, 
                                    PriorOdds.a = PriorOdds.a, PriorOdds.b = PriorOdds.b,PriorOdds.med = PriorOdds.med,
                       TypeI = TypeI, R = R, seed = seeds_sens[1])
    outpower.PriorOdds = out1$outpower
    if(nrow(PriorOdds.combs)>=2){
      for(j in 2:nrow(PriorOdds.combs)){
        
        PriorOdds.a=PriorOdds.combs[j,1]; PriorOdds.b=PriorOdds.combs[j,2]; PriorOdds.med=PriorOdds.combs[j,3];
        if(is.na(PriorOdds.a)){PriorOdds.a=NULL}; if(is.na(PriorOdds.b)){PriorOdds.b=NULL}; if(is.na(PriorOdds.med)){PriorOdds.med=NULL}
        
        sens1 = power.bfmed.sensPriorOdds(bfabfb_pools =out1$bfabfb_pools, N = N, a=a,b=b,cp = cp, PriorOdds.a = PriorOdds.a, PriorOdds.b = PriorOdds.b,PriorOdds.med = PriorOdds.med,
                                          TypeI = TypeI, R = R, seed = seeds_sens[j] )
        # print( round(sens1[, paste0(c('power_'),TypeI )],2) )
        outpower.PriorOdds = rbind(outpower.PriorOdds, sens1)
      }
    }
    print(outpower.PriorOdds)
    
    outpower.PriorOdds.Nseq = rbind(outpower.PriorOdds.Nseq, outpower.PriorOdds)
    
    if (  min(outpower.PriorOdds[, paste0(c('power_'),TypeI )])  >=Power.desired){

      break # not run next N in Nseq
    }
  }
  out = outpower.PriorOdds.Nseq
  # out = list(outpower.PriorOdds.Nseq=outpower.PriorOdds.Nseq, outpower.PriorOdds=outpower.PriorOdds)
  return(out)
}
  

# aa=SSD.bfmed(Power.desired = .7,Nmin=28, Nmax=100, Nstep=10, R = 120 )





