library(BayesFactor)
library(tidyverse)
library(mvtnorm)
library(lavaan)
## R functions for conducting Bayesian hypothesis testing of mediation
## source('fun_BHTmed1.R') to source the R functions in this file

# independent paths a and b ----
# BF.a, BF.b need to be computed from the data 
# the function below uses the "BayesFactor" package to compute path Bayes factors for the simple mediation model
BFaBFb <- function(dat
                   ,DepVar='Y', Mediator='M', IndVar='X'
                   ,Covariates.apath=NULL, Covariates.bpath=NULL
){
  Y=dat[,colnames(dat)==DepVar]
  M=dat[,colnames(dat)==Mediator]
  X=dat[,colnames(dat)==IndVar]
  
  # path a
  if( is.null(Covariates.apath) ){
    bfa_full = lmBF(as.formula(
      paste0(Mediator,"~",IndVar)
    ),
    data = dat) ## Bayes factor of full model against intercept only
    
    BF.a = extractBF(bfa_full)$bf 
  }
  if( !is.null(Covariates.apath) ){
    bfa_full = lmBF(as.formula(
      paste0(Mediator,"~",paste(IndVar,Covariates.apath,sep = "+"))
    ),
    data = dat) ## Bayes factor of full model against intercept only
    bfa_redu = lmBF(as.formula(
      paste0(Mediator,"~",paste(Covariates.apath,sep = "+"))
    ),
    data = dat) ## Bayes factor of full model against intercept only
    
    BF.a = extractBF(bfa_full)$bf / extractBF(bfa_redu)$bf
  }
  
  # path b
  if( is.null(Covariates.bpath) ){
    bfb_full = lmBF(as.formula(
      paste0(DepVar,"~",paste(Mediator,IndVar,sep = "+"))
    ),
    data = dat) ## Bayes factor of full model against intercept only
    bfb_redu = lmBF(as.formula(
      paste0(DepVar,"~",paste(IndVar,sep = "+"))
    ),
    data = dat) ## Bayes factor of full model against intercept only
    
    BF.b = extractBF(bfb_full)$bf / extractBF(bfb_redu)$bf
  }
  if( !is.null(Covariates.bpath) ){
    bfb_full = lmBF(as.formula(
      paste0(DepVar,"~",paste(Mediator,IndVar,Covariates.bpath,sep = "+"))
    ),
    data = dat) ## Bayes factor of full model against intercept only
    bfb_redu = lmBF(as.formula(
      paste0(DepVar,"~",paste(IndVar,Covariates.bpath,sep = "+"))
    ),
    data = dat) ## Bayes factor of full model against intercept only
    
    BF.b = extractBF(bfb_full)$bf / extractBF(bfb_redu)$bf
  }
  
  
  # namesout=c('BF.a','BF.b')
  # out=unlist(mget(namesout, envir = environment()))
  out=c(BF.a,BF.b)
  return(out)
}

# (a) specify two path prior odds 
# for example, as PriorOdds.a=PriorOdds.b=1

pathb.a <- function(PriorOdds.a, PriorOdds.b, BF.a, BF.b){
  PriorOdds.med = PriorOdds.b*PriorOdds.a / (1+PriorOdds.b+PriorOdds.a)
  BF.med = BF.a*BF.b * (1+PriorOdds.b+PriorOdds.a) /
    (1 + BF.b*PriorOdds.b + BF.a*PriorOdds.a)
  PosteriorOdds.med = BF.med * PriorOdds.med
  
  PosteriorOdds.a = PriorOdds.a*BF.a
  PosteriorOdds.b = PriorOdds.b*BF.b
  
  namesout=lapply(c('PriorOdds.','BF.','PosteriorOdds.'), paste0, c('a','b','med'))
  out=mget(unlist(namesout), envir = environment())
  out=unlist(out)
  
  out1=matrix(out,3,3, dimnames = list(c('Path a','Path b','Mediation'),c('Prior Odds','Bayes Factor','Posterior Odds')))
  out1=as.data.frame(out1)
  out1$`Posterior Prob`=(out1$`Posterior Odds`)/(1+out1$`Posterior Odds`)
  
  return(out1)
}

# (b1) specify prior odds for the mediation effect and for path a
med.a <- function(PriorOdds.med, PriorOdds.a, BF.a, BF.b){
  
  PriorOdds.b = PriorOdds.med*(1+PriorOdds.a) / (PriorOdds.a-PriorOdds.med)
  BF.med = BF.a*BF.b * (1+PriorOdds.b+PriorOdds.a) /
    (1 + BF.b*PriorOdds.b + BF.a*PriorOdds.a)
  PosteriorOdds.med = BF.med * PriorOdds.med
  
  PosteriorOdds.a = PriorOdds.a*BF.a
  PosteriorOdds.b = PriorOdds.b*BF.b
  
  namesout=lapply(c('PriorOdds.','BF.','PosteriorOdds.'), paste0, c('a','b','med'))
  out=mget(unlist(namesout), envir = environment())
  out=unlist(out)
  
  out1=matrix(out,3,3, dimnames = list(c('Path a','Path b','Mediation'),c('Prior Odds','Bayes Factor','Posterior Odds')))
  out1=as.data.frame(out1)
  out1$`Posterior Prob`=(out1$`Posterior Odds`)/(1+out1$`Posterior Odds`)
  
  return(out1)
}


# (b2) specify prior odds for the mediation effect and for path b
med.b <- function(PriorOdds.med, PriorOdds.b, BF.a, BF.b){
  
  PriorOdds.a = PriorOdds.med*(1+PriorOdds.b) / (PriorOdds.b-PriorOdds.med)
  BF.med = BF.a*BF.b * (1+PriorOdds.b+PriorOdds.a) /
    (1 + BF.b*PriorOdds.b + BF.a*PriorOdds.a)
  PosteriorOdds.med = BF.med * PriorOdds.med
  
  PosteriorOdds.a = PriorOdds.a*BF.a
  PosteriorOdds.b = PriorOdds.b*BF.b
  
  namesout=lapply(c('PriorOdds.','BF.','PosteriorOdds.'), paste0, c('a','b','med'))
  out=mget(unlist(namesout), envir = environment())
  out=unlist(out)
  
  out1=matrix(out,3,3, dimnames = list(c('Path a','Path b','Mediation'),c('Prior Odds','Bayes Factor','Posterior Odds')))
  out1=as.data.frame(out1)
  out1$`Posterior Prob`=(out1$`Posterior Odds`)/(1+out1$`Posterior Odds`)
  
  return(out1)
}



# dependent paths a and b ----

# specify mediation prior odds and the three conditional prior probabilities under the null hypothesis of no mediation 
med.qs <- function(PriorOdds.med, H0.q10, H0.q01, H0.q00, 
         BIC11, BIC01,BIC10,BIC00){
  q10=H0.q10; q01=H0.q01; q00=H0.q00
  BF.med = exp((BIC00-BIC11)/2) / (1*q00 + exp((BIC00-BIC01)/2)*q01 + exp((BIC00-BIC10)/2)*q10)
  PosteriorOdds.med = BF.med * PriorOdds.med
  
  p11=PriorOdds.med/(1+PriorOdds.med)
  p10=(1-p11)*q10
  p01=(1-p11)*q01
  p00=(1-p11)*q00

  PriorOdds.a = (PriorOdds.med+q10)/(1-q10)
  PriorOdds.b = (PriorOdds.med+q01)/(1-q01)
  PosteriorOdds.a = (p11*exp((BIC00-BIC11)/2)+p10*exp((BIC00-BIC10)/2)) / (p01*exp((BIC00-BIC01)/2)+p00*1)
  PosteriorOdds.b = (p11*exp((BIC00-BIC11)/2)+p01*exp((BIC00-BIC01)/2)) / (p10*exp((BIC00-BIC10)/2)+p00*1)
  BF.a=PosteriorOdds.a/PriorOdds.a
  BF.b=PosteriorOdds.b/PriorOdds.b
  
  namesout=lapply(c('PriorOdds.','BF.','PosteriorOdds.'),
                  paste0, c('a','b','med'))
  out=mget(unlist(namesout), envir = environment())
  out=unlist(out)
  
  out1=matrix(out,3,3, dimnames = list(c('Path a','Path b','Mediation'),c('Prior Odds','Bayes Factor','Posterior Odds')))
  out1=as.data.frame(out1)
  out1$`Posterior Prob`=(out1$`Posterior Odds`)/(1+out1$`Posterior Odds`)
  
  return(out1)
}

# specify the mediation prior odds and the prior odds of individual paths a and b
med.ab <- function(PriorOdds.med, PriorOdds.a, PriorOdds.b,  
                   BIC11, BIC01,BIC10,BIC00){
  q10 = (PriorOdds.a-PriorOdds.med)/(1+PriorOdds.a)
  q01 = (PriorOdds.b-PriorOdds.med)/(1+PriorOdds.b)
  q00 = 1-q10-q01
  
  p11=PriorOdds.med/(1+PriorOdds.med)
  p10=(1-p11)*q10
  p01=(1-p11)*q01
  p00=(1-p11)*q00
  
  # PriorOdds.a = (PriorOdds.med+q10)/(1-q10)
  # PriorOdds.b = (PriorOdds.med+q01)/(1-q01)
  PosteriorOdds.a = (p11*exp((BIC00-BIC11)/2)+p10*exp((BIC00-BIC10)/2)) / (p01*exp((BIC00-BIC01)/2)+p00*1)
  PosteriorOdds.b = (p11*exp((BIC00-BIC11)/2)+p01*exp((BIC00-BIC01)/2)) / (p10*exp((BIC00-BIC10)/2)+p00*1)
  BF.a=PosteriorOdds.a/PriorOdds.a
  BF.b=PosteriorOdds.b/PriorOdds.b
  
  # BF.med = BIC11 / (BIC00*q00 + BIC01*q01 + BIC10*q10)
  BF.med = exp((BIC00-BIC11)/2) / (1*q00 + exp((BIC00-BIC01)/2)*q01 + exp((BIC00-BIC10)/2)*q10)
  PosteriorOdds.med = BF.med * PriorOdds.med
  
  namesout=lapply(c('PriorOdds.','BF.','PosteriorOdds.'),
                  paste0, c('a','b','med'))
  out=mget(unlist(namesout), envir = environment())
  out=unlist(out)
  
  out1=matrix(out,3,3, dimnames = list(c('Path a','Path b','Mediation'),c('Prior Odds','Bayes Factor','Posterior Odds')))
  out1=as.data.frame(out1)
  out1$`Posterior Prob`=(out1$`Posterior Odds`)/(1+out1$`Posterior Odds`)
  
  return(out1)
}




