
library(BayesFactor)
library(mvtnorm)
library(parallel)



BFaBFb.oneDesign <- function(
    N = 100, 
    std.a = 0.25, std.b = 0.26, std.cp = 0.16, 
    p = 0, # number of covariates C's
    Rsq.xc = 0.2, Rsq.mc = 0.2, Rsq.yc = 0.2, #  respectively the proportions of the variance explained by all C's for the treatment X, mediator M, and outcome Y, without controlling for any variables.
    ## Not used arguments: 
    # Design.PriorOdds.a=1, Design.PriorOdds.b=1,
    # Analysis.PriorOdds.a=1, Analysis.PriorOdds.b=1,
    # cutoff.BF = 3, absolute.cutoff = TRUE,
    # cutoff.FPR = 0.05, relative.cutoff = TRUE,
    R = 1e3,seed = 12345
) {
  
  ## Simulate a pool of BFa, BFb (for given N, and effect sizes) to sample from. 
  ## p(BFmed | H0) is obtained by drawing (BFa, BFb) from p(BFa,BFb | each of the three no-mediation scenarios of H0) according to the conditional prior probabilities of the three scenarios under H0. 
  ## This pool of (BFa, BFb) can be used for obtaining p(BFmed | H0) under varied prior odds specifications, so as to assess the sensitivity to the prior odds.
  
  Nsim_H0 <- Nsim_H1 <- R 
  seed_H0 <- seed
  seed_H1 <- seed*2
  
  a <- std.a
  b <- std.b
  cp <- std.cp
  
  # covariates' coefficients (assuming all C's have the same coefficients)
  if(p > 0) {
    coeff_x.c <- sqrt(Rsq.xc) 
    coeff_x.onec <- coeff_x.c / sqrt(p)
    coeff_m.c <- sqrt(Rsq.mc) - a * sqrt(Rsq.xc)
    coeff_m.onec <- coeff_m.c / sqrt(p)
    coeff_y.c <- sqrt(Rsq.yc) - b * sqrt(Rsq.mc) - cp * sqrt(Rsq.xc)
    coeff_y.onec <- coeff_y.c / sqrt(p)
    # error variances in the standardized model
    var_ex <- 1 - Rsq.xc
    var_em <- 1 - (a^2 * var_ex + Rsq.mc)
    var_ey <- 1 - ((cp + b*a)^2 * var_ex + b^2 * var_em + Rsq.yc)
  }
  if(p == 0) {
    var_ex <- 1
    var_em <- 1 - (a^2 * var_ex)
    var_ey <- 1 - ((cp + b*a)^2 * var_ex + b^2 * var_em)
  }
  
  
  bfabfb.samplingdist <- function(N, a1, b1, cp, R, seed_onecond=816051) {
    set.seed(seed = seed_onecond)
    seeds_OneData.cond = sample(1:1e8, R)
    
    OneData.cond <- function(i=1, N, a1, b1, cp) {
      set.seed(seed = seeds_OneData.cond[i])
      
      if(p > 0) {
        Cvec <- rmvnorm(N, mean = rep(0, p))
        X <- rnorm(N, 0, sqrt(var_ex)) + coeff_x.onec * rowSums(Cvec)
        M <- a1 * X + rnorm(N, 0, sqrt(var_em)) + coeff_m.onec * rowSums(Cvec)
        Y <- cp * X + b1 * M + rnorm(N, 0, sqrt(var_ey)) + coeff_y.onec *  rowSums(Cvec)
        # Bayes factors of paths a and b
        colnames(Cvec) <- paste0("C", 1:p)
        dat1 <- data.frame(X, M, Y, Cvec)
        fM_num <- paste0("M ~ X + ", paste0("C", 1:p, collapse = " + "))
        fM_den <- paste0("M ~ ", paste0("C", 1:p, collapse = " + ")) 
        BFa <- extractBF(lmBF(formula(fM_num), data = dat1))$bf / extractBF(lmBF(formula(fM_den), data = dat1))$bf 
        fY_num <- paste0("Y ~ M + X + ", paste0("C", 1:p, collapse = " + "))
        fY_den <- paste0("Y ~ X + ", paste0("C", 1:p, collapse = " + ")) 
        BFb <- extractBF(lmBF(formula(fY_num), data = dat1))$bf / extractBF(lmBF(formula(fY_den), data = dat1))$bf 
      }
      
      if(p == 0) {
        X <- rnorm(N, 0, sqrt(var_ex)) 
        M <- a1 * X + rnorm(N, 0, sqrt(var_em)) 
        Y <- cp * X + b1 * M + rnorm(N, 0, sqrt(var_ey))
        # Bayes factors of paths a and b
        dat1 <- data.frame(X, M, Y)
        fM_num <- paste0("M ~ X")
        fM_den <- paste0("M ~ 1") 
        BFa <- extractBF(lmBF(formula(fM_num), data = dat1))$bf # Against denominator: Intercept only 
        fY_num <- paste0("Y ~ M + X")
        fY_den <- paste0("Y ~ X") 
        BFb <- extractBF(lmBF(formula(fY_num), data = dat1))$bf / extractBF(lmBF(formula(fY_den), data = dat1))$bf 
      }
      
      res <- c(BFa,BFb)
      names(res) <- c('BFa','BFb')
      
      out <- res
      return(out)
    }
    
    bfabfb.Nsim_cond= mclapply(1:R, OneData.cond,
                               N=N,a1=a1,b1=b1, cp=cp,
                               mc.cores = 5,
                               mc.preschedule = TRUE)
    
    # output 
     
    designspeci <- c(N = N, a = a1, b = b1, cp = cp)
    bfabfb_out <- cbind(t(designspeci), data.frame(do.call(rbind, bfabfb.Nsim_cond)))
    out <- data.frame(bfabfb_out)
    
    return(out)
  }
  
  set.seed(seed = seed)
  abconds <- data.frame(expand.grid(a1 = c(a, 0), b1 = c(b,0)))
  seeds_abconds <- sample(1:4e4, nrow(abconds))
  
  bfabfb_pools <- NULL
  
  for(i in 1:nrow(abconds)) {
    bfabfb_out <- bfabfb.samplingdist(N = N,a1 = abconds$a1[i], b1 = abconds$b1[i], cp = cp, R = R, seed_onecond = seeds_abconds[i])
    bfabfb_pools <- rbind(bfabfb_pools, bfabfb_out)
  }
  
  bfabfb_pools.list <- mget( ls(), envir = environment())
  
  return(bfabfb_pools.list)
}



BFaBFb.oneDesign.random <- function(
    N = 100, 
    std.a = 0.25, std.b = 0.26, std.cp = 0.16, 
    sigma.a = 0.107, sigma.b = 0.111, sigma.cp = 0.110,
    p = 0, # number of covariates C's
    Rsq.xc = 0.2, Rsq.mc = 0.2, Rsq.yc = 0.2, # fixed
    # Design.PriorOdds.a=1, Design.PriorOdds.b=1,
    # Analysis.PriorOdds.a=1, Analysis.PriorOdds.b=1,
    # cutoff.BF = 3, absolute.cutoff = TRUE,
    # cutoff.FPR = 0.05, relative.cutoff = TRUE,
    R = 1e3,seed = 12345
) {
  
  Nsim_H0 <- Nsim_H1 <- R 
  seed_H0 <- seed
  seed_H1 <- seed*2
  
  # consider the standardized model
  a <- std.a
  b <- std.b
  cp <- std.cp
  
  if(p > 0) {
    coeff_x.c <- sqrt(Rsq.xc) 
    coeff_x.onec <- coeff_x.c / sqrt(p)
    coeff_m.c <- sqrt(Rsq.mc) - a * sqrt(Rsq.xc)
    coeff_m.onec <- coeff_m.c / sqrt(p)
    coeff_y.c <- sqrt(Rsq.yc) - b * sqrt(Rsq.mc) - cp * sqrt(Rsq.xc)
    coeff_y.onec <- coeff_y.c / sqrt(p)
    # error variances in the standardized model
    var_ex <- 1 - Rsq.xc
    var_em <- 1 - (a^2 * var_ex + Rsq.mc)
    var_ey <- 1 - ((cp + b*a)^2 * var_ex + b^2 * var_em + Rsq.yc)
  }
  if(p == 0) {
    var_ex <- 1
    var_em <- 1 - (a^2 * var_ex)
    var_ey <- 1 - ((cp + b*a)^2 * var_ex + b^2 * var_em)
  }
  
  # consider uncertainty in the treatment-mediator, mediator-outcome, and treatment-outcome relations, i.e., a, b, cp
  m.a <- std.a
  m.b <- std.b
  m.cp <- std.cp
  
  bfabfb.samplingdist.random <- function(
    N, 
    m.a1, sigma.a1,
    m.b1, sigma.b1, 
    m.cp, sigma.cp,
    R, seed_onecond = 816051 
  ) {
    set.seed(seed = seed_onecond)
    seeds_OneData.cond = sample(1:1e8, R)
    
    OneData.cond.random <- function(
    i=1, N, 
    m.a1, sigma.a1,
    m.b1, sigma.b1, 
    m.cp, sigma.cp
    ) {
      # randomly draw a set of parameter values from the design prior distribution to generate each dataset
      set.seed(seed = seeds_OneData.cond[i])
      a1 = rnorm(1, m.a1, sigma.a1)
      b1 = rnorm(1, m.b1, sigma.b1)
      cp = rnorm(1, m.cp, sigma.cp)
      
      if(p > 0) {
        Cvec <- rmvnorm(N, mean = rep(0, p))
        X <- rnorm(N, 0, sqrt(var_ex)) + coeff_x.onec * rowSums(Cvec)
        M <- a1 * X + rnorm(N, 0, sqrt(var_em)) + coeff_m.onec * rowSums(Cvec)
        Y <- cp * X + b1 * M + rnorm(N, 0, sqrt(var_ey)) + coeff_y.onec *  rowSums(Cvec)
        # Bayes factors of paths a and b
        colnames(Cvec) <- paste0("C", 1:p)
        dat1 <- data.frame(X, M, Y, Cvec)
        fM_num <- paste0("M ~ X + ", paste0("C", 1:p, collapse = " + "))
        fM_den <- paste0("M ~ ", paste0("C", 1:p, collapse = " + ")) 
        BFa <- extractBF(lmBF(formula(fM_num), data = dat1))$bf / extractBF(lmBF(formula(fM_den), data = dat1))$bf 
        fY_num <- paste0("Y ~ M + X + ", paste0("C", 1:p, collapse = " + "))
        fY_den <- paste0("Y ~ X + ", paste0("C", 1:p, collapse = " + ")) 
        BFb <- extractBF(lmBF(formula(fY_num), data = dat1))$bf / extractBF(lmBF(formula(fY_den), data = dat1))$bf 
      }
      
      if(p == 0) {
        X <- rnorm(N, 0, sqrt(var_ex)) 
        M <- a1 * X + rnorm(N, 0, sqrt(var_em)) 
        Y <- cp * X + b1 * M + rnorm(N, 0, sqrt(var_ey))
        # Bayes factors of paths a and b
        dat1 <- data.frame(X, M, Y)
        fM_num <- paste0("M ~ X")
        fM_den <- paste0("M ~ 1") 
        BFa <- extractBF(lmBF(formula(fM_num), data = dat1))$bf # Against denominator: Intercept only 
        fY_num <- paste0("Y ~ M + X")
        fY_den <- paste0("Y ~ X") 
        BFb <- extractBF(lmBF(formula(fY_num), data = dat1))$bf / extractBF(lmBF(formula(fY_den), data = dat1))$bf 
      }
      
      res <- c(BFa,BFb)
      names(res) <- c('BFa','BFb')
      
      out <- res
      return(out)
    }
    
    bfabfb.Nsim_cond <- mclapply(1:R, OneData.cond.random,
                                 N=N, 
                                 m.a1=m.a1, sigma.a1=sigma.a1,
                                 m.b1=m.b1, sigma.b1=sigma.b1, 
                                 m.cp=m.cp, sigma.cp=sigma.cp,
                                 mc.cores = 5,
                                 mc.preschedule = TRUE)
    
    # output 
    
    designspeci <- c(N = N, a = m.a1, b = m.b1, cp = m.cp)
    bfabfb_out <- cbind(t(designspeci), data.frame(do.call(rbind, bfabfb.Nsim_cond)))
    out <- data.frame(bfabfb_out)
    
    return(out)
  }
  
  set.seed(seed = seed); 
  abconds <- data.frame(expand.grid(m.a1=c(m.a, 0), m.b1=c(m.b,0)))
  abconds$sigma.a1 <- (abconds$m.a1 != 0) * sigma.a + (abconds$m.a1 == 0) * 0 
  abconds$sigma.b1 <- (abconds$m.b1 != 0) * sigma.b + (abconds$m.b1 == 0) * 0
  abconds$m.cp <- m.cp
  abconds$sigma.cp <- sigma.cp
  
  seeds_abconds <- sample(1:4e4, nrow(abconds))
  
  bfabfb_pools <- NULL
  
  for(i in 1:nrow(abconds)){
    bfabfb_out <- bfabfb.samplingdist.random(
      N = N,
      m.a1 = abconds$m.a1[i], sigma.a1 = abconds$sigma.a1[i],
      m.b1 = abconds$m.b1[i], sigma.b1 = abconds$sigma.b1[i],
      m.cp = abconds$m.cp[i], sigma.cp = abconds$sigma.cp[i],
      R = R, seed_onecond = seeds_abconds[i] 
    )
    
    bfabfb_pools <- rbind(bfabfb_pools, bfabfb_out)
  }
  a <- m.a
  b <- m.b
  cp <- m.cp
  
  bfabfb_pools.list = mget(ls(), envir = environment())
  
  return(bfabfb_pools.list)
}


BFmed.oneDesign <- function(
    bfabfb_pools.list, 
    Design.PriorOdds.a, 
    Design.PriorOdds.b
) {
  list2env(bfabfb_pools.list, envir = environment())
  # design prior odds specification 
  PriorOdds.a <- Design.PriorOdds.a
  PriorOdds.b <- Design.PriorOdds.b
  
  if((!is.null(PriorOdds.a)) & (!is.null(PriorOdds.b))) { 
    # using the independent paths assumption
    PriorOdds.med <- 1 / (1/(PriorOdds.a*PriorOdds.b) + 1/PriorOdds.a + 1/PriorOdds.b) # applicable when PriorOdds.a=Inf
  }
  if((!is.null(PriorOdds.a)) & (!is.null(PriorOdds.med))) { 
    PriorOdds.b <- PriorOdds.med*(1/PriorOdds.a+ 1) / (1 - PriorOdds.med/PriorOdds.a)
  }
  if((!is.null(PriorOdds.b)) & (!is.null(PriorOdds.med))){ 
    PriorOdds.a <- PriorOdds.med*(1/PriorOdds.b+ 1) / (1-PriorOdds.med/PriorOdds.b)
  }
  # for any priorodds
  Design.q10 <- q10 <- (1 - PriorOdds.med/PriorOdds.a) / (1/PriorOdds.a + 1) ##applicable when PriorOdds.a=Inf
  Design.q01 <- q01 <- (1 - PriorOdds.med/PriorOdds.b) / (1/PriorOdds.b + 1)
  Design.q00 <- q00 <- 1 - q10 - q01
  
  Nsim_H0 <- Nsim_H1 <- R 
  Nplan <- N 
  
  asize= a
  bsize= b
  cpsize= cp
  
  set.seed(seed = seed)
  seeds_sampleqs = sample(3e3, 3)
  
  # sample under H0
  aa <- bfabfb_pools[(bfabfb_pools$N==Nplan & bfabfb_pools$a==0 & bfabfb_pools$b==0 & bfabfb_pools$cp==cpsize), ]
  sample_q00 <- NULL
  set.seed(seed = seeds_sampleqs[1])
  if(nrow(aa) > 0) {
    sample_q00 <- aa[sample.int(nrow(aa), round(Nsim_H0*q00), replace = F), ] 
  }
  
  aa <- bfabfb_pools[(bfabfb_pools$N==Nplan & bfabfb_pools$a==asize & bfabfb_pools$b==0 & bfabfb_pools$cp==cpsize), ]
  sample_q10 <- NULL
  set.seed(seed = seeds_sampleqs[2])
  if(nrow(aa) > 0) {
    sample_q10 <- aa[sample.int(nrow(aa), round(Nsim_H0*q10), replace = F), ]  
  }
  
  aa <- bfabfb_pools[(bfabfb_pools$N==Nplan & bfabfb_pools$a==0 & bfabfb_pools$b==bsize & bfabfb_pools$cp==cpsize), ]
  sample_q01 <- NULL
  set.seed(seed = seeds_sampleqs[3])
  if(nrow(aa) > 0) {
    sample_q01 <- aa[sample.int(nrow(aa), round(Nsim_H0*q01), replace = F), ]  
  }
  
  sample_H0 <- rbind(rbind(sample_q00, sample_q10), sample_q01)
  
  # sample under H1, used all of (bfa, bfb) under H1
  aa <- bfabfb_pools[(bfabfb_pools$N==Nplan & bfabfb_pools$a==asize & bfabfb_pools$b==bsize & bfabfb_pools$cp==cpsize), ]
  sample_H1 <- NULL;  
  if(nrow(aa) > 0) {
    sample_H1 <- aa[sample.int(nrow(aa), round(Nsim_H1*1), replace = F), ]   
  }
  
  # combine the BFa BFb from H0 and H1
  bfabfb_Nsim_H0 <- data.frame(sample_H0, trueH='H0')
  bfabfb_Nsim_H1 <- data.frame(sample_H1, trueH='H1')
  bfabfb_Nsim_H0H1 <- rbind(bfabfb_Nsim_H0, bfabfb_Nsim_H1) 
  
  bfmed_Nsim_H0H1 <- bfabfb_Nsim_H0H1
  
  bfmed_Nsim_H0H1.list <- mget(ls(), envir = environment())
  
  return(bfmed_Nsim_H0H1.list)
}

BFmed.oneAnalysis <- function( 
    bfmed_Nsim_H0H1.list, 
    Analysis.PriorOdds.a, 
    Analysis.PriorOdds.b,
    cutoff.BF, absolute.cutoff,
    cutoff.FPR, relative.cutoff
) {
  
  list2env(bfmed_Nsim_H0H1.list, envir = environment())
  # calculate BFmed using BFa, BFb, and the analysis prior odds 
  if((!is.null(Analysis.PriorOdds.a)) & (!is.null(Analysis.PriorOdds.a))) { 
    # using the independent paths assumption
    Analysis.PriorOdds.med <- 1 / (1/(Analysis.PriorOdds.a*Analysis.PriorOdds.b) + 1/Analysis.PriorOdds.a + 1/Analysis.PriorOdds.b) ## applicable when PriorOdds.a=Inf
  }
   
  ## for any priorodds
  Analysis.q10 <- (1 - Analysis.PriorOdds.med/Analysis.PriorOdds.a) / (1/Analysis.PriorOdds.a + 1) ##applicable when PriorOdds.a=Inf
  Analysis.q01 <- (1 - Analysis.PriorOdds.med/Analysis.PriorOdds.b) / (1/Analysis.PriorOdds.b + 1)
  Analysis.q00 <- 1 - Analysis.q10 - Analysis.q01
  
  bfmed_Nsim_H0H1$BFmed <- with(bfmed_Nsim_H0H1, {
    BFmed <- BFa*BFb / (Analysis.q00 + BFb*Analysis.q01 + BFa*Analysis.q10)
  })
  
  if(relative.cutoff == TRUE){
    # relative to false positive rate
    TypeI <- cutoff.FPR
    cutoff_bfmed <- quantile(bfmed_Nsim_H0H1$BFmed[bfmed_Nsim_H0H1$trueH=='H0'], probs = c(1-TypeI))
    
    TPR_bfmed <- power_bfmed <- sapply(cutoff_bfmed, function(x) {
      mean(bfmed_Nsim_H0H1$BFmed[bfmed_Nsim_H0H1$trueH=='H1'] > x)
    })  # true positive rate
    
  #   # symmetric cutoff for supporting H0
  #   TNR_bfmed <- sapply(cutoff_bfmed, function(x){mean( bfmed_Nsim_H0H1$BFmed[bfmed_Nsim_H0H1$trueH=='H0'] < 1/x )}) # true negative rate
  #   
  #   FNR_bfmed <- sapply(cutoff_bfmed, function(x){mean( bfmed_Nsim_H0H1$BFmed[bfmed_Nsim_H0H1$trueH=='H1'] < 1/x )}) # false negative rate
  # 
  }
  
  if(absolute.cutoff == TRUE){
    cutoff_abs <- cutoff.BF
    ## absolute cutoff for BF
    FPR_bfmed_abs <- TypeI_bfmed_abs <- sapply(cutoff_abs, function(x) {
      mean(bfmed_Nsim_H0H1$BFmed[bfmed_Nsim_H0H1$trueH=='H0'] > x)
    })
    
    TPR_bfmed_abs <- power_bfmed_abs <- sapply(cutoff_abs, function(x) {
      mean(bfmed_Nsim_H0H1$BFmed[bfmed_Nsim_H0H1$trueH=='H1'] > x )
    })
    
    # TNR_bfmed_abs = sapply(cutoff_abs, function(x){mean( 
    #   bfmed_Nsim_H0H1$BFmed[bfmed_Nsim_H0H1$trueH=='H0'] < 1/x )})
    # 
    # FNR_bfmed_abs = sapply(cutoff_abs, function(x){mean( 
    #   bfmed_Nsim_H0H1$BFmed[bfmed_Nsim_H0H1$trueH=='H1'] < 1/x )})
  }
  
  ## output 
  designprior <- c(
    Design.PriorOdds.a=Design.PriorOdds.a, 
    Design.PriorOdds.b=Design.PriorOdds.b, 
    Design.q00=Design.q00, Design.q10=Design.q10, Design.q01=Design.q01
  )
  analysisprior <- c(
    Analysis.PriorOdds.a=Analysis.PriorOdds.a, 
    Analysis.PriorOdds.b=Analysis.PriorOdds.b, 
    Analysis.q00=Analysis.q00, Analysis.q10=Analysis.q10, Analysis.q01=Analysis.q01
  )
  N <- Nplan
  a <- asize
  b <- bsize
  cp <- cpsize
  designspeci <- c(N=N,a=a,b=b,cp=cp)
  
  res <- data.frame(
    N,
    # true_positive_relative.cut = TPR_bfmed,     
    # false_positive_relative.cut = cutoff_bfmed, 
    # 
    # true_positive_absolute.cut = TPR_bfmed_abs, 
    # false_positive_absolute.cut = TypeI_bfmed_abs, 
    
    true_positive_bfmed = TPR_bfmed,     
    false_positive_cutoff = cutoff_bfmed, 
    
    true_positive_bfmed_abs = TPR_bfmed_abs, 
    false_positive_abs = TypeI_bfmed_abs, 
    
    # TNR_bfmed, FNR_bfmed, 
    # TNR_bfmed_abs, FNR_bfmed_abs,
    a, b, cp, t(designprior), t(analysisprior)
  )
  rownames(res) <- NULL
  
  
  out <- res
  return(out)
}




BFmed.SSD <- function(
    N = seq(100, 200, by=50), 
    
    std.a = 0.25, std.b = 0.26, std.cp = 0.16, 
    uncertain.effect = FALSE,
    #sigma.a = 0.107, sigma.b = 0.111, sigma.cp = 0.110, 
    sigma.a = 0, sigma.b = 0, sigma.cp = 0, 
    
    p = 0, # number of covariates C's
    Rsq.xc = 0.2, Rsq.mc = 0.2, Rsq.yc = 0.2, #  respectively the proportions of the variance explained by all C's for the treatment X, mediator M, and outcome Y, without controlling for any variables.
    
    Design.PriorOdds.a = c(1, 3, 20), Design.PriorOdds.b = 1,
    Analysis.PriorOdds.a = 1, Analysis.PriorOdds.b = 1,
    
    cutoff.BF = 3, absolute.cutoff = TRUE,
    cutoff.FPR = 0.05, relative.cutoff = TRUE,
    R=1e3, seed=12345
) {
  
  datacond <- data.frame(expand.grid(std.a = std.a, std.b = std.b, std.cp = std.cp, N = N))
  
  priorspeci <- data.frame(expand.grid(
    Design.PriorOdds.a = Design.PriorOdds.a, Design.PriorOdds.b = Design.PriorOdds.b
  ))
  
  analysis <- data.frame(expand.grid(
    Analysis.PriorOdds.a = Analysis.PriorOdds.a, Analysis.PriorOdds.b=Analysis.PriorOdds.b
  ))
  
  if(uncertain.effect == FALSE){
    mcmapply(BFaBFb.oneDesign, N = datacond$N, 
             std.a = datacond$std.a, std.b = datacond$std.b, std.cp = datacond$std.cp, MoreArgs = list(p=p, Rsq.xc = Rsq.xc, Rsq.mc = Rsq.mc, Rsq.yc = Rsq.yc, R=R, seed=seed), SIMPLIFY = F, mc.preschedule = T, mc.cores = 1) -> bfabfb_pools
  }
  
  if(uncertain.effect == TRUE){
    
    datacond <- data.frame(expand.grid( 
      std.a = std.a, std.b = std.b, std.cp = std.cp,
      sigma.a = sigma.a, sigma.b=sigma.b, sigma.cp = sigma.cp,
      N = N 
    ))
    
    mcmapply(BFaBFb.oneDesign.random, N = datacond$N, 
           std.a = datacond$std.a, std.b = datacond$std.b, std.cp = datacond$std.cp, sigma.a = datacond$sigma.a, sigma.b = datacond$sigma.b, sigma.cp = datacond$sigma.cp, MoreArgs = list(p=p, Rsq.xc = Rsq.xc, Rsq.mc = Rsq.mc, Rsq.yc = Rsq.yc, R=R, seed=seed), SIMPLIFY = F, mc.preschedule = T, mc.cores = 1) -> bfabfb_pools
  }
  
  out.DP <- NULL
  for(d in 1:length(bfabfb_pools)) {
    mcmapply(BFmed.oneDesign, 
             Design.PriorOdds.a=priorspeci$Design.PriorOdds.a, Design.PriorOdds.b=priorspeci$Design.PriorOdds.b, SIMPLIFY = F, mc.preschedule = T,
             mc.cores = 1,
             MoreArgs = list(bfabfb_pools.list = bfabfb_pools[[ d ]])
    ) -> bfmed_Nsim_H0H1.d
    
    for(p in 1:length(bfmed_Nsim_H0H1.d)){
      mcmapply(BFmed.oneAnalysis,
               Analysis.PriorOdds.a = analysis$Analysis.PriorOdds.a, Analysis.PriorOdds.b = analysis$Analysis.PriorOdds.b,
               MoreArgs = list(bfmed_Nsim_H0H1.list = bfmed_Nsim_H0H1.d[[ p ]],
                               cutoff.BF = cutoff.BF, absolute.cutoff = absolute.cutoff,
                               cutoff.FPR = cutoff.FPR, relative.cutoff = relative.cutoff),
               SIMPLIFY = F, mc.preschedule = T,
               mc.cores = 1
      ) -> out.dp
      out.dp <- do.call(rbind, out.dp)
      
      out.DP <- rbind(out.DP, out.dp)
    }
    
  }
  
  out.DP <- out.DP %>% rename( 
    true_positive_relative.cut = true_positive_bfmed,
    true_positive_absolute.cut = true_positive_bfmed_abs,
    relative.cut = false_positive_cutoff,
    false_positive_absolute.cut = false_positive_abs
  ) %>% 
    select(!contains(".q")) 
  SSD_BFmed <- out.DP
  SSD_BFmed$cutoff.BF <- cutoff.BF
  SSD_BFmed$cutoff.FPR <- cutoff.FPR
  
  return(SSD_BFmed)
  
}




# plots -------------------

plot.true.positive.med <- function(SSD_BFmed) {
  
  SSD_BFmed %>%   
    mutate( Design.PriorOdds.b = round(Design.PriorOdds.b, 2), 
            Design.PriorOdds.a = round(Design.PriorOdds.a, 2),
            Analysis.PriorOdds.a = round(Analysis.PriorOdds.a, 2),
            Analysis.PriorOdds.b = round(Analysis.PriorOdds.b),
            Analysis.PriorOdds = paste0("Analysis.PriorOdds.a = ", Analysis.PriorOdds.a, ",\nAnalysis.PriorOdds.b = ", Analysis.PriorOdds.b)
    ) -> plotdf1 
  
  TPR.plot <- ggplot(plotdf1) +
    labs(title = "True positive rate of detecting mediation with the mediation BF") +
    facet_grid( Design.PriorOdds.a + Design.PriorOdds.b ~ a + b,  
                labeller = label_both ) +
    
    geom_point( aes(x=N, y=true_positive_relative.cut, shape=Analysis.PriorOdds), 
                size=2.5 ) +
    geom_line( aes(x=N, y=true_positive_relative.cut, color= "relative to the specified false positive rate,\n cutoff.FPR") ) +
    geom_point( aes(x=N, y=true_positive_absolute.cut, shape=Analysis.PriorOdds),
                size=2.5 ) +
    geom_line( aes(x=N, y=true_positive_absolute.cut, color="the specified absolute cutoff,\n cutoff.BF") ) +
    geom_hline( yintercept = 0.8, size=0.6, color="green4" ) +
    scale_color_manual( name = "Cutoff for BFmed", values = c("blue", "red3") ) +
    # scale_color_discrete(name="Design: \n PriorOdds.a: 1 \n PriorOdds.b:") +
    scale_shape_discrete( name="Analysis" ) +
    
    scale_x_continuous( name = "Sample size N", breaks = unique(plotdf1$N) ) +
    scale_y_continuous( name = "True positive rate of BFmed" ) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 13),
          title = element_text(size = 13),
          legend.text = element_text(size = 12),
          legend.key.height = unit(1.2, "cm")
          # ,legend.spacing.y = unit(10)
    )
  TPR.plot
}


plot.false.positive.med <- function(SSD_BFmed) {
  SSD_BFmed %>%   
    mutate( Design.PriorOdds.b = round(Design.PriorOdds.b, 2), 
            Design.PriorOdds.a = round(Design.PriorOdds.a, 2),
            Analysis.PriorOdds.a = round(Analysis.PriorOdds.a, 2),
            Analysis.PriorOdds.b = round(Analysis.PriorOdds.b),
            Analysis.PriorOdds = paste0("Analysis.PriorOdds.a = ", Analysis.PriorOdds.a, ",\nAnalysis.PriorOdds.b = ", Analysis.PriorOdds.b)
    ) -> plotdf1 
  
  FPR.plot <- ggplot( plotdf1 ) +
    labs(title ="False positive rate of detecting mediation with the mediation BF") +
    facet_grid( Design.PriorOdds.a + Design.PriorOdds.b ~ a + b,  labeller = label_both ) +
    
    geom_point( aes(x=N,y=cutoff.FPR,
                    shape=Analysis.PriorOdds), size=2.5 ) +
    geom_line( aes(x=N, y=cutoff.FPR, color="relative to the specified false positive rate,\n cutoff.FPR") ) +
    geom_point( aes(x=N,y=false_positive_absolute.cut, shape=Analysis.PriorOdds), size=2.5 ) +
    geom_line( aes(x=N, y=false_positive_absolute.cut, color="the specified absolute cutoff,\n cutoff.BF") ) +
    # geom_hline( yintercept = 0.05, size=0.6, color="green4" ) +
    scale_color_manual( name = "Cutoff for BFmed", values = c("blue", "red3") ) +
    # scale_color_discrete(name="Design: \n PriorOdds.a: 1 \n PriorOdds.b:") +
    scale_shape_discrete( name = "Analysis" ) +
    
    scale_x_continuous( name = "Sample size N", breaks = unique(plotdf1$N)  ) +
    scale_y_continuous( name = "False positive rate of BFmed") +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          axis.text = element_text(size = 12),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 13),
          title = element_text(size = 13),
          legend.text = element_text(size = 12),
          legend.key.height = unit(1.2, "cm")
          # ,legend.spacing.y = unit(10)
    )
  
  FPR.plot
}

