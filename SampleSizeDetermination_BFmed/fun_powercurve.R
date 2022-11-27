library(tidyverse)
library(ggplot2)

# default Nplan.start 
n.sigasigb=function(a_true,b_true,cp_true,alpha,power.desired){
  z_alpha=qnorm(1-alpha/2)
  z_beta=qnorm(power.desired)
  n_siga=( (1-a_true^2)/a_true^2 )*( z_alpha+z_beta )^2
  n_sigb=( (1-b_true^2-cp_true^2-2*a_true*b_true*cp_true)/(b_true^2*(1-a_true^2)) )*( z_alpha+z_beta )^2
  n_sigasigb=max(n_siga,n_sigb)
  return(n_sigasigb)
}
# Nplan.start=ceiling(n.sigasigb(a_true,b_true,cp_true,alpha,power.desired))
# Nplan_seq=round(seq(Nplan.start, 250, by = 2 ))

# conventional power analysis: cp ----
joint.power.ind = function(Nplan,a_true,b_true,cp_true,alpha){
  ncp_a=sqrt(Nplan)*a_true/sqrt(1-a_true^2)
  pr_a=1-pt(q=qt(1-alpha/2,Nplan-2), df=Nplan-2, ncp=ncp_a)+pt(q=-qt(1-alpha/2,Nplan-2), df=Nplan-2, ncp=ncp_a)
  ncp_b=sqrt(Nplan)*b_true/sqrt((1-b_true^2-cp_true^2-2*a_true*b_true*cp_true)/(1-a_true^2))
  pr_b=1-pt(q=qt(1-alpha/2,Nplan-3), df=Nplan-3, ncp=ncp_b)+pt(q=-qt(1-alpha/2,Nplan-3), df=Nplan-3, ncp=ncp_b)
  power_ind=pr_a*pr_b
  return(power_ind)
}

# res_joint.power.ind=mclapply(Nplan_seq, joint.power.ind, 
#                              a_true,b_true,cp_true,alpha
#                              , mc.cores = mccores, mc.preschedule = F)
# powercp=data.frame( power=c(unlist(res_joint.power.ind) ), type="joint.power.ind",nplan=Nplan_seq)


# our pd method: for one Nplan, a distribution of power ----
pd.jointind=function(Nplan,nrawdata,a_true,b_true,cp_true,alpha,Npilot){
  powerdist=NULL
  
  joint.powerind.1pilot=function(a_true,b_true,cp_true,alpha,Npilot,Nplan){
    N=Npilot
    Var.e_M=1-a_true^2
    Var.e_Y=1-b_true^2-cp_true^2-2*a_true*b_true*cp_true
    # simulate a sample with the orginal effect size estimates
    x=rnorm(N)
    m=a_true*x+rnorm(N,0,sqrt(Var.e_M))
    y=b_true*m+cp_true*x+rnorm(N,0,sqrt(Var.e_Y))
    # fit the simple mediation model to the simulated sample and obtain a set of parameter estimates 
    txx=t(x)%*%x
    tmx=t(m)%*%x
    tmm=t(m)%*%m
    ahat=tmx/txx
    bhat=(t(c(txx)*m-c(tmx)*x)%*%y)/(txx*tmm-tmx*tmx) 
    cphat=(t(c(tmm)*x-c(tmx)*m)%*%y)/(txx*tmm-tmx*tmx) 
    Var.e_M.hat=t(m-c(ahat)*x)%*%(m-c(ahat)*x)/(N-1)
    Var.e_Y.hat=t(y-c(cphat)*x-c(bhat)*m)%*%(y-c(cphat)*x-c(bhat)*m)/(N-2)  
    
    # obtain a power value with this set of parameter estimates
    power.ind=joint.power.ind(Nplan=Nplan,a_true=ahat,b_true=bhat,cp_true=cphat,alpha=alpha)
    
    return(power.ind)
  }
  # obtain #nrawdata power values
  for(raw in 1:nrawdata){
    temp=joint.powerind.1pilot(a_true=a_true,b_true=b_true,cp_true=cp_true,alpha=alpha,Npilot=Npilot,Nplan=Nplan)
    powerdist=c(powerdist,temp)
  }
  return(powerdist)
}

# pdjointind=mclapply(Nplan_seq, pd.jointind,
#                     nrawdata=1000,a_true,b_true,cp_true,alpha
#                     ,Npilot=125
#                     ,mc.cores = mccores, mc.preschedule = F)







