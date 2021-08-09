

# independent paths a and b ----
BF.a=3
BF.b=1/3

# (a) specify two path prior odds 
PriorOdds.a=PriorOdds.b=1

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
  return(out)
}

# (b) specify prior odds for the mediation effect  
PriorOdds.med=1
# specify prior odds for path a
PriorOdds.a=3
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
  
  return(out)
}
# specify prior odds for path b
PriorOdds.b=4
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
  
  return(out)
}



# (c) specify a mediation prior odds 
#and one of the three conditional prior probabilities under the null.

# specify q01
q01=1/2
med.q01<-function(PriorOdds.med, H0.q01, BF.a, BF.b){
  q01=H0.q01
  q00 = (1-q01)*q01 / (PriorOdds.med+q01)
  q10 = (1-q01)*PriorOdds.med / (PriorOdds.med+q01)
  BF.med = BF.a*BF.b / (q00 + BF.b*q01 + BF.a*q10)
  PosteriorOdds.med = BF.med * PriorOdds.med
  
  PriorOdds.a = (PriorOdds.med+q10)/(1-q10)
  PriorOdds.b = (PriorOdds.med+q01)/(1-q01)
  PosteriorOdds.a = PriorOdds.a*BF.a
  PosteriorOdds.b = PriorOdds.b*BF.b
  
  namesout=lapply(c('PriorOdds.','BF.','PosteriorOdds.'), paste0, c('a','b','med'))
  out=mget(unlist(namesout), envir = environment())
  out=unlist(out)

  return(out)
}

# specify q10
q10=1/2
med.q10 <- function(PriorOdds.med, H0.q10, BF.a, BF.b){
  q10=H0.q10
  q00 = (1-q10)*q10 / (PriorOdds.med+q10)
  q01 = (1-q10)*PriorOdds.med / (PriorOdds.med+q10)
  BF.med = BF.a*BF.b / (q00 + BF.b*q01 + BF.a*q10)
  PosteriorOdds.med = BF.med * PriorOdds.med
  
  PriorOdds.a = (PriorOdds.med+q10)/(1-q10)
  PriorOdds.b = (PriorOdds.med+q01)/(1-q01)
  PosteriorOdds.a = PriorOdds.a*BF.a
  PosteriorOdds.b = PriorOdds.b*BF.b
  
  namesout=lapply(c('PriorOdds.','BF.','PosteriorOdds.'), paste0, c('a','b','med'))
  out=mget(unlist(namesout), envir = environment())
  out=unlist(out)
  
  return(out)
}

# specify q00
q00 = 1/8
which.smaller='q10'
med.q00 <- function(PriorOdds.med, H0.q00, which.smaller='q10', BF.a, BF.b){
  q00=H0.q00
  c=1-q00
  d=PriorOdds.med*q00
  if( c^2-4*d>=0 ){
    q.lager = (c+sqrt(c^2-4*d))/2
    q.smaller =(c-sqrt(c^2-4*d))/2
    if (which.smaller=='q10') {q10=q.smaller; q01=q.lager}
    if (which.smaller=='q01') {q01=q.smaller; q10=q.lager}
    BF.med = BF.a*BF.b / (q00 + BF.b*q01 + BF.a*q10)
    PosteriorOdds.med = BF.med * PriorOdds.med
    
    PriorOdds.a = (PriorOdds.med+q10)/(1-q10)
    PriorOdds.b = (PriorOdds.med+q01)/(1-q01)
    PosteriorOdds.a = PriorOdds.a*BF.a
    PosteriorOdds.b = PriorOdds.b*BF.b
    
    namesout=lapply(c('PriorOdds.','BF.','PosteriorOdds.'), paste0, c('a','b','med'))
    out=mget(unlist(namesout), envir = environment())
    out=unlist(out)
    
  }
  if( c^2-4*d<0 ){
    BF.med = NA
    PosteriorOdds.med = NA
    
    PriorOdds.a = NA
    PriorOdds.b = NA
    PosteriorOdds.a = NA
    PosteriorOdds.b = NA
    
    namesout=lapply(c('PriorOdds.','BF.','PosteriorOdds.'), paste0, c('a','b','med'))
    out=mget(unlist(namesout), envir = environment())
    out=unlist(out)
  }
  
  return(out)
}




# dependent paths a and b ----

# specify med and q
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
  
  return(out)
}

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
  
  return(out)
}






