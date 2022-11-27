
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
# a=.39; b=.39; cp=0.39
# N=Nplan=70

bfmed.power <- function(N=100, a=0.39,b=0.39,cp=0.39, PriorOdds.a=NULL, PriorOdds.b=NULL, PriorOdds.med=NULL, TypeI=0.05,R=1e4,seed=816051){
  Nsim_H0=Nsim_H1=R # =1e4
  seed_H0=seed #=816051; 
  seed_H1=seed*2
  
  ## error variances in standardized simple mediation model
  var_em=1-a^2; var_ey=1-cp^2-b^2-2*a*b*cp
  
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
  

  ## number of simulated-datasets under H0
  # Nsim_H0 = 1e4 ## Fu et al. (2021) recommend at least 1e4
  set.seed(seed=seed_H0); seeds_OneData.H0 <- sample(1:1e8, Nsim_H0)
  whichq_H0 = rmultinom(n=Nsim_H0, 1, prob = c(q10,q01,q00)) ## for each data, only one model can be true
  
  OneData.H0 <- function(i=1, N, a, b, cp){
    set.seed(seed = seeds_OneData.H0[i])
    X=rnorm(N)
    whichq.i=c('q10','q01','q00') [which(whichq_H0[,i]==1)]
    if(whichq.i == 'q01') { ## a=0, b=b
      M=0*X+rnorm(N, 0, sqrt(var_em)) 
      Y=cp*X+b*M+rnorm(N,0, sqrt(var_ey))
    }
    if(whichq.i == 'q10') { ## a=a, b=0
      M=a*X+rnorm(N, 0, sqrt(var_em)) 
      Y=cp*X+0*M+rnorm(N,0, sqrt(var_ey))
    }
    if(whichq.i == 'q00') { ## a=0, b=0
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
                           N=N,a=a,b=b, cp=cp,
                           mc.cores = 5,
                           mc.preschedule = F)
  ptm1=proc.time()[3]; #print(ptm1-ptm0)
  
  
  # Nsim_H1=1e4 ### recommend at least 1e4
  set.seed(seed=seed_H1); seeds_OneData.H1 <- sample(1:1e8, Nsim_H1)
  OneData.H1 <- function(i=1, N, a, b, cp){
    set.seed(seed = seeds_OneData.H1[i])
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
                           N=N,a=a,b=b, cp=cp,
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
  designspeci=c(N=N,a=a,b=b,cp=cp)
  
  outpower=cbind(t(priorspeci),t(designspeci), t(cutoff_bfmed),t(power_bfmed))
  outpower=data.frame(outpower)
  # colnames(outpower)[c(11:16)]=c(paste0(c('cutoff_bfmed_'),TypeI),paste0(c('power_'),TypeI )) 
  
  out=list(outpower=outpower, bfmed_Nsim_H0H1=bfmed_Nsim_H0H1)
  return(out)
}


## to speed up the SSD, start from the required 
SSD.bfmed <- function(Power.desired=0.8, Nmin=50,Nmax=500,Nstep=10, a=0.39,b=0.39,cp=0, PriorOdds.a=1, PriorOdds.b=1, PriorOdds.med=NULL,TypeI=0.05, R=1e4,seed=816051) {
  a=a;b=b;cp=cp;PriorOdds.a=PriorOdds.a; PriorOdds.b=PriorOdds.b; PriorOdds.med=PriorOdds.med;TypeI=TypeI; R=R;seed=seed;
  
  Nseq = seq(Nmin, Nmax, by=Nstep)
  outpower.Nseq=NULL
  for(N in Nseq){
    out1 = bfmed.power(N = N, a=a,b=b,cp = cp, PriorOdds.a = PriorOdds.a, PriorOdds.b = PriorOdds.b,PriorOdds.med = PriorOdds.med,
                       TypeI = TypeI, R = R, seed = seed)
    print(  out1$outpower[, paste0(c('power_'),TypeI )] )
    outpower.Nseq = rbind(outpower.Nseq, out1$outpower)
    
    if (  out1$outpower[, paste0(c('power_'),TypeI )] >=Power.desired){
      out_SSD=out1
      break
    }
  }
  
  # out = list(outpower.Nseq=outpower.Nseq, out_SSD=out_SSD)
  out = outpower.Nseq
  return(out)
}

# ptm0=proc.time()[3]
# aa=SSD.bfmed(Power.desired = .7,Nmin=28, Nmax=100, Nstep=10, R = 40 )
# ptm1=proc.time()[3]; print(ptm1-ptm0)


