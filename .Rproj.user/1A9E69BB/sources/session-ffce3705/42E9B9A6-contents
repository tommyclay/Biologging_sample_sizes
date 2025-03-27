library(lme4)
library(lmtest)
library(multcomp)

in_direx <- "./Data_outputs/"

#### 1. LOAD IN PREDICTIONS ####

# 50% and 95% UDs area occupied predictions for pop size of 100
pred1 <- read.csv(paste0(in_direx, "XX.csv"))

# 50% and 95% UDs for a range of population sizes
pred2 <- read.csv(paste0(in_direx, "XX.csv"))

# predictions from the four model types
pred3 <-  read.csv(paste0(in_direx, "XX.csv"))

# birdlife representativeness predictions
pred_bl <- read.csv(paste0(in_direx, "birdlife_predictions.csv"))
# convert to pred 4 by adding area occupied predictions as well
pred4 <- data.frame()

###### REMOVE THIS BELOW ##########

# comparing asymptotic model types
pred1_50 <- data.frame(species = rep(c("AP", "AP", "AP", "BBA", "BBA", "BP", "BS", "BS", "BS", "GHA", "LMA", "LMA", "NGP", "SGP", "WA", "WA", "WA", "WA", "WCP", "WCP", "WCP"),4),
                      year = rep(c("10", "11", "All", "02", "03", "11", "02", "12", "Both", "99", "04", "05", "00", "00", "04", "08", "09", "All", "03", "13", "Both"), 4),
                      pred = c(7, 18, 12,  8,  7, 3, 16, 19, 19, 20, 11, 8, 9, 4, 30, 25, 18, 22, 21, 18, 26,
                               6, 14, 10, 10, 12, 3, 14, 16, 15, 19, 10, 8, 8, 4, 30, 24, 17, 21, 18, 17, 25,
                               3,  4,  4,  4,  3, 2,  4,  6,  6,  6,  4, 3, 4, 2,  9,  9,  6,  9,  6,  5,  8,
                               3,  5,  5,  4,  4, 2,  5,  7,  7,  8,  6, 4, 5, 3, 10, 11,  9, 12,  6,  6, 10),
                      model = rep(c("MM2", "MM3", "AR2", "AR3"), each = 21))
pred1_95 <- data.frame(species = rep(c("AP", "AP", "AP", "BBA", "BBA", "BP", "BS", "BS", "BS", "GHA", "LMA", "LMA", "NGP", "SGP", "WA", "WA", "WA", "WA", "WCP", "WCP", "WCP"),4),
                      year = rep(c("10", "11", "All", "02", "03",  "11", "02", "12", "Both", "99", "04", "05", "00", "00", "04","08", "09", "All", "03", "13", "Both"), 4),
                      pred = c(12, 13, 17, 32, 25, 29, 14, 28, 32, 29, 18, 14,  9, 19, 27, 31, 34, 33, 21, 15, 20,
                               14, 15, 25, 37, 28, 33, 15, 29, 36, 36, 20, 16, 11, 21, 28, 31, 39, 35, 23, 15, 20,
                               4,  4,  5, 11,  8,  7,  4,  9, 11,  9,  6,  4,  3,  7,  8, 11, 12, 14,  6,  4,  6,
                               6,  5, 10, 16, 12, 10,  5, 11, 15, 15, 10,  7,  7, 10,  8, 15, 18, 21,  8,  6,  8),
                      model = rep(c("MM2", "MM3", "AR2", "AR3"), each = 21))

# comparing for different range of pop sizes
pred2_50 <- data.frame(species = rep(c("AP", "AP", "AP", "BBA", "BBA", "BP", "BS", "BS", "BS", "GHA", "LMA", "LMA", "NGP", "SGP", "WA", "WA", "WA", "WA", "WCP", "WCP", "WCP"),5),
                      year = rep(c("10", "11", "All", "02", "03", "11", "02", "12", "Both", "99", "04", "05", "00", "00", "04", "08", "09", "All", "03", "13", "Both"), 5),
                      pred = c(6, 14, 10, 10, 12, 3, 14, 16, 15, 19, 10, 8, 8, 4, 30, 24, 17, 21, 18, 17, 25,
                               6, 16, 11, 11, 13, 3, 16, 18, 18, 23, 10, 8, 9, 4, 42, 30, 21, 26, 21, 20, 32, 
                               6, 16, 12, 11, 14, 3, 16, 19, 18, 23, 10, 8, 9, 4, 44, 31, 21, 27, 21, 20, 33,
                               6, 16, 12, 11, 14, 3, 16, 19, 18, 23, 10, 8, 9, 4, 44, 31, 21, 27, 21, 21, 33,
                               6, 16, 12, 11, 14, 3, 16, 19, 18, 23, 10, 8, 9, 4, 44, 31, 21, 27, 21, 21, 33),
                      pop = rep(c(100, 1000, 10000, 100000, 1000000), each = 21))
pred2_95 <- data.frame(species = rep(c("AP", "AP", "AP", "BBA", "BBA", "BP", "BS", "BS", "BS", "GHA", "LMA", "LMA", "NGP", "SGP", "WA", "WA", "WA", "WA", "WCP", "WCP", "WCP"),5),
                      year = rep(c("10", "11", "All", "02", "03", "11", "02", "12", "Both", "99", "04", "05", "00", "00", "04", "08", "09", "All", "03", "13", "Both"), 5),
                      pred = c(14, 15, 25, 37, 28, 33, 15, 29, 36, 36, 20, 16, 11, 21, 28, 31, 39, 35, 23, 15, 20,
                               16, 17, 32, 58, 39, 48, 17, 40, 54, 54, 24, 16, 13, 26, 38, 44, 61, 52, 29, 17, 24,
                               16, 17, 33, 61, 40, 51, 17, 42, 57, 57, 25, 19, 13, 27, 40, 46, 65, 54, 30, 17, 24,
                               16, 17, 33, 61, 40, 51, 17, 42, 58, 58, 25, 20, 13, 27, 40, 46, 65, 55, 30, 17, 24,
                               16, 17, 33, 61, 40, 51, 17, 42, 58, 58, 25, 20, 13, 27, 40, 46, 65, 55, 30, 17, 24),
                      pop = rep(c(100, 1000, 10000, 100000, 1000000), each = 21))
# comparing 50% and 95% UD predictions
pred3 <- data.frame(species = rep(c("AP", "AP", "AP", "BBA", "BBA", "BP", "BS", "BS", "BS", "GHA", "LMA", "LMA", "NGP", "SGP", "WA", "WA", "WA", "WA", "WCP", "WCP", "WCP"),2),
                    year = rep(c("10", "11", "All", "02", "03", "11", "02", "12", "Both", "99", "04", "05", "00", "00", "04", "08", "09", "All", "03", "13", "Both"), 2),
                    pred = c( 6, 14, 10, 10,  12, 3, 14, 16, 15, 19, 10,  8,  8,  4, 30, 24, 17, 21, 18, 17, 25,
                              14, 15, 25, 37, 28, 33, 15, 29, 36, 36, 20, 16, 11, 21, 28, 31, 39, 35, 23, 15, 20),
                    iso = rep(c(50, 95), each = 21))
# comparing different methods
pred4 <- data.frame(species = rep(c("AP", "AP", "AP", "BBA", "BBA", "BP", "BS", "BS", "BS", "GHA", "LMA", "LMA", "NGP", "SGP", "WA", "WA", "WA", "WA", "WCP", "WCP", "WCP"),3),
                    year = rep(c("10", "11", "All", "02", "03", "11", "02", "12", "Both", "99", "04", "05", "00", "00", "04", "08", "09", "All", "03", "13", "Both"), 3),
                    pred = c(6,  14, 10, 10, 12,  3, 14, 16, 15, 19, 10, 8,  8,  4, 30, 24, 17, 21, 18, 17, 25,
                               2,  11,  3,  2,  2,  3,  6,  6,  6,  5,  3, 2,  2,  1,  7,  6,  5,  5,  6,  5,  6,
                               11, 68, 22, 12, 12, 32, 28, 43, 46, 37, 22, 9, 12,  5, 55, 47, 36, 41, 44, 30, 44),
                    method = rep(c("Total Area Occupied", "Inclusion70", "Inclusion95"), each = 21))


##### UP TO HERE #####



#### 2. AREA OCCUPIED - COMPARE 50% & 95% UD PREDICTIONS AMONG DIFFERENT ASYMPTOTIC MODEL TYPES ####

##### 2A. 50% UD ######

pred1_50$model <- as.factor(as.character(pred1_50$model))
pred1_50$sp_yr <- as.factor(as.character(paste(pred1_50$species, pred1_50$year, sep = "_")))

mod1_50 <- glmer(pred ~ model + (1|sp_yr), data = pred1_50, family = "poisson")
summary(mod1_50)
mod1_50_2 <- update(mod1_50,~.- model)
lrtest(mod1_50, mod1_50_2)
# Model 1: pred ~ model + (1 | sp_yr)
# Model 2: pred ~ (1 | sp_yr)
# #Df  LogLik Df  Chisq            Pr(>Chisq)    
# 1   5 -205.90                                    
# 2   2 -295.78 -3 179.75 < 0.00000000000000022 ***
plot(mod1_50)

# tukey's tests for adjusted pairwise comparisons
summary(glht(mod1_50, mcp(model="Tukey")))
# Linear Hypotheses:
#                Estimate  Std. Error z value  Pr(>|z|)    
# AR3 - AR2 == 0  0.20997    0.12946   1.622    0.360    
# MM2 - AR2 == 0  1.09861    0.11109   9.889  <0.0001 ***
# MM3 - AR2 == 0  1.03428    0.11201   9.234  <0.0001 ***
# MM2 - AR3 == 0  0.88865    0.10290   8.636  <0.0001 ***
# MM3 - AR3 == 0  0.82432    0.10389   7.934  <0.0001 ***
# MM3 - MM2 == 0 -0.06433    0.07985  -0.806    0.849   

##### 2B. 95% UD ######

pred1_95$model <- as.factor(as.character(pred1_95$model))
pred1_95$sp_yr <- as.factor(as.character(paste(pred1_95$species, pred1_95$year, sep = "_")))

mod1_95 <- glmer(pred ~ model + (1|sp_yr), data = pred1_95, family = "poisson")
summary(mod1_95)
mod1_95_2 <- update(mod1_95,~.- model)
lrtest(mod1_95, mod1_95_2)
# Model 1: pred ~ model + (1 | sp_yr)
# Model 2: pred ~ (1 | sp_yr)
# #Df  LogLik Df  Chisq            Pr(>Chisq)    
# 1   5 -225.31                                    
# 2   2 -379.80 -3 308.96 < 0.00000000000000022 ***
plot(mod1_95)

summary(glht(mod1_95, mcp(model="Tukey")))
# Linear Hypotheses:
# Estimate Std. Error z value Pr(>|z|)    
# AR3 - AR2 == 0  0.37673    0.10467   3.599  0.00176 ** 
# MM2 - AR2 == 0  1.12654    0.09276  12.145  < 0.001 ***
# MM3 - AR2 == 0  1.23676    0.09156  13.507  < 0.001 ***
# MM2 - AR3 == 0  0.74981    0.08102   9.255  < 0.001 ***
# MM3 - AR3 == 0  0.86003    0.07965  10.798  < 0.001 ***
# MM3 - MM2 == 0  0.11022    0.06319   1.744  0.29402    


#### 3. AREA OCCUPIED - COMPARE 50% & 95% UD PREDICTIONS FOR RANGE OF POP SIZES ####

##### 3A. 50% UD ######

pred2_50$pop <- as.factor(as.character(pred2_50$pop))
pred2_50$sp_yr <- as.factor(as.character(paste(pred2_50$species, pred2_50$year, sep = "_")))

mod2_50 <- glmer(pred ~ pop + (1|sp_yr), data = pred2_50, family = "poisson")
summary(mod2_50)
# remove effect of population size and test significance
mod2_50_2 <- update(mod2_50,~.- pop)
lrtest(mod2_50, mod2_50_2)
# Model 1: pred ~ pop + (1 | sp_yr)
# Model 2: pred ~ (1 | sp_yr)
# #Df  LogLik Df  Chisq Pr(>Chisq)  
# 1   6 -282.68                       
# 2   2 -287.52 -4 9.6857    0.04607 *
plot(mod2_50)

# tukey's tests for adjusted pairwise comparisons
summary(glht(mod2_50, mcp(pop="Tukey")))
# Linear Hypotheses:
#                           Estimate   Std. Error z value Pr(>|z|)  
# 1000 - 100 == 0        0.173424053  0.077988386   2.224   0.1708  
# 10000 - 100 == 0       0.195523144  0.077598128   2.520   0.0861 .
# 100000 - 100 == 0      0.198253176  0.077550375   2.556   0.0786 .
# 1000000 - 100 == 0     0.198250214  0.077550431   2.556   0.0786 .
# 10000 - 1000 == 0      0.022099092  0.074130881   0.298   0.9983  
# 100000 - 1000 == 0     0.024829123  0.074080895   0.335   0.9973  
# 1000000 - 1000 == 0    0.024826162  0.074080951   0.335   0.9973  
# 100000 - 10000 == 0    0.002730031  0.073669941   0.037   1.0000  
# 1000000 - 10000 == 0   0.002727070  0.073669996   0.037   1.0000  
# 1000000 - 100000 == 0 -0.000002961  0.073619698   0.000   1.0000  


##### 3A. 95% UD ######

pred2_95$pop <- as.factor(as.character(pred2_95$pop))
pred2_95$sp_yr <- as.factor(as.character(paste(pred2_95$species, pred2_95$year, sep = "_")))

mod2_95 <- glmer(pred ~ pop + (1|sp_yr), data = pred2_95, family = "poisson")
summary(mod2_95)
# remove effect of population size and test significance
mod2_95_2 <- update(mod2_95,~.- pop)
lrtest(mod2_95, mod2_95_2)
# Model 1: pred ~ pop + (1 | sp_yr)
# Model 2: pred ~ (1 | sp_yr)
# #Df  LogLik Df  Chisq Pr(>Chisq)    
# 1   6 -326.29                         
# 2   2 -355.96 -4 59.341  3.989e-12 ***
plot(mod2_95)

# tukey's tests for adjusted pairwise comparisons
summary(glht(mod2_95, mcp(pop="Tukey")))
# Linear Hypotheses:
#                       Estimate    Std. Error   z value Pr(>|z|)    
# 1000 - 100 == 0       0.310665611 0.057263910   5.425 <0.00001 ***
# 10000 - 100 == 0      0.354209109 0.056745570   6.242 <0.00001 ***
# 100000 - 100 == 0     0.359519863 0.056683566   6.343 <0.00001 ***
# 1000000 - 100 == 0    0.359521188 0.056683551   6.343 <0.00001 ***
# 10000 - 1000 == 0     0.043543497 0.052103325   0.836    0.919    
# 100000 - 1000 == 0    0.048854252 0.052035789   0.939    0.881    
# 1000000 - 1000 == 0   0.048855576 0.052035773   0.939    0.881    
# 100000 - 10000 == 0   0.005310755 0.051464821   0.103    1.000    
# 1000000 - 10000 == 0  0.005312079 0.051464804   0.103    1.000    
# 1000000 - 100000 == 0 0.000001325 0.051396429   0.000    1.000    


#### 4. AREA OCCUPIED - COMPARE 50% & 95% UD PREDICTIONS FOR POP SIZE OF 100 ####

pred3$iso <- as.factor(as.character(pred3$iso))
pred3$sp_yr <- as.factor(as.character(paste(pred3$sp, pred3$yr, sep = "_")))

# running glmm
mod3 <- glmer(pred ~ iso + (1|sp_yr), data = pred3, family = "poisson")
summary(mod3)
# remove effect of isopleth and test significance
mod3_2 <- update(mod3,~.- iso)
# lrtest(mod3, mod3_2)
# Model 1: pred ~ iso + (1 | sp_yr)
# Model 2: pred ~ (1 | sp_yr)
# #Df  LogLik Df  Chisq           Pr(>Chisq)    
# 1   3 -146.28                                   
# 2   2 -177.51 -1 62.476 0.000000000000002698 ***
plot(mod3)

# calculating mean and range for each isopleth category

# core
mean(pred3$pred[pred3$iso == "50"]) # 14.33333
sd(pred3$pred[pred3$iso == "50"]) # 7.059273
range(pred3$pred[pred3$iso == "50"]) # 3 30

# general use
mean(pred3$pred[pred3$iso == "95"]) # 25.09524
sd(pred3$pred[pred3$iso == "95"]) #  8.898903
range(pred3$pred[pred3$iso == "95"]) # 11 39


#### 5. COMPARING PREDICTIONS FROM DIFFERENT METHODS - AREA OCCUPIED VS BIRDLIFE REPRESENTATITVENESS ####

pred4$method <- as.factor(as.character(pred4$method))
pred4$sp_yr <- as.factor(as.character(paste(pred4$species, pred4$yr, sep = "_")))

mod4 <- glmer(pred ~ method + (1|sp_yr), data = pred4, family = "poisson")
summary(mod4)
mod4_2 <- update(mod4,~.- method)
lrtest(mod4_2, mod4)
# Model 1: pred ~ (1 | sp_yr)
# Model 2: pred ~ method + (1 | sp_yr)
# #Df  LogLik Df  Chisq            Pr(>Chisq)    
# 1   2 -454.44                                    
# 2   4 -212.30  2 484.29 < 0.00000000000000022 ***
  
# tukey's tests for adjusted pairwise comparisons
summary(glht(mod4, mcp(method ="Tukey")))
# Linear Hypotheses:
#                                         Estimate Std. Error z value            Pr(>|z|)    
# Inclusion95 - Inclusion70 == 0          1.94286    0.11006  17.654 <0.0000000000000002 ***
# Total Area Occupied - Inclusion70 == 0  1.16382    0.11791   9.871 <0.0000000000000002 ***
# Total Area Occupied - Inclusion95 == 0 -0.77904    0.06947 -11.214 <0.0000000000000002 ***

