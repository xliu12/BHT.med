# Sample size determination (SSD) for testing mediation with the mediation Bayes factor

# illustrative example

# The R functions
source("fun.BFmed.ssd.R")


# Standardized path coefficients in the simple mediation model
# Liu et al. (2019) Example 1
a<-0.25; b<-0.26;cp<-0.16



# Fixed effect sizes ------------
res1 <- BFmed.SSD(
  N=seq(100, 200, by=5), 
  std.a=0.25, std.b=0.26, std.cp=0.16, 
  Design.PriorOdds.a=c(1/10, 1/3, 1), Design.PriorOdds.b=1,
  Analysis.PriorOdds.a=1, Analysis.PriorOdds.b=1,
  cutoff.BF = 3, absolute.cutoff = TRUE,
  cutoff.FPR = 0.05, relative.cutoff = TRUE,
  R=1e4, seed=12345
) 

# Table
library(tidyverse)

res1 %>% 
  mutate(Design.PriorOdds.b = round(Design.PriorOdds.b, 1), 
         Design.PriorOdds.a= round(Design.PriorOdds.a, 1) ) %>%
  select(N, Design.PriorOdds.a, Design.PriorOdds.b,
         true_positive_bfmed, false_positive_cutoff, true_positive_bfmed_abs, false_positive_abs) %>%
  pivot_wider( names_from = c(Design.PriorOdds.a), 
               names_glue = "Design.PriorOdds.a {Design.PriorOdds.a}\n{.value}" ,
               values_from = c(true_positive_bfmed, false_positive_cutoff, true_positive_bfmed_abs, false_positive_abs) ) ->tab1

tab1

# Plot 
SSD_BFmed <- res1
plot.true.positive.med( SSD_BFmed )
plot.false.positive.med( SSD_BFmed )


# Uncertainty in effect sizes --------
ustd.a = .194; 
std.a = .250
# std.a/ustd.a = s_x / s_y
# also holds for the standard error and confidence limits
se_ustd.a = 0.083 # 
se_ustd.a * (std.a/ustd.a) -> se_std.a

ustd.b = 0.089; std.b = .26 # .258 approximately .26
se_ustd.b = 0.038
se_ustd.b * (std.b/ustd.b) -> se_std.b

ustd.cp = 0.042; std.cp=0.16 
se_ustd.cp = 0.029
se_ustd.cp * (std.cp/ustd.cp) -> se_std.cp


ssdres2 <- BFmed.SSD(
  N=seq(100, 400, by=5), 
  std.a=0.25, std.b=0.26, std.cp=0.16, 
  uncertain.effect = TRUE,
  sigma.a = se_std.a, sigma.b = se_std.b, sigma.cp = se_std.cp,
  Design.PriorOdds.a=c(1/10, 1/3, 1), Design.PriorOdds.b=1,
  Analysis.PriorOdds.a=1, Analysis.PriorOdds.b=1,
  cutoff.BF = 3, absolute.cutoff = TRUE,
  cutoff.FPR = 0.05, relative.cutoff = TRUE,
  R=1e4, seed=12345
) 


# Table

ssdres2 %>% 
  mutate(Design.PriorOdds.b = round(Design.PriorOdds.b, 1), 
          Design.PriorOdds.a= round(Design.PriorOdds.a, 1) ) %>%
  select(N, Design.PriorOdds.a, Design.PriorOdds.b,
                true_positive_bfmed, false_positive_cutoff, true_positive_bfmed_abs, false_positive_abs) %>%
  pivot_wider( names_from = c(Design.PriorOdds.a), 
               names_glue = "Design.PriorOdds.a {Design.PriorOdds.a}\n{.value}" ,
               values_from = c(true_positive_bfmed, false_positive_cutoff, true_positive_bfmed_abs, false_positive_abs) ) ->tab2

tab2

# Plot 
SSD_BFmed <- ssdres2
plot.true.positive.med( SSD_BFmed )
plot.false.positive.med( SSD_BFmed )
