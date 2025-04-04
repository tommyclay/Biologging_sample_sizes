library(lme4)
library(lmtest)
library(multcomp)
library(emmeans)

in_direx <- "./Data_outputs/"

#### 1. LOAD IN PREDICTIONS ####

# load in area occupied file
pred_50 <- read.csv(paste0(in_direx, "area_occupied_50_predictions.csv"))
pred_95 <- read.csv(paste0(in_direx, "area_occupied_95_predictions.csv"))

# 1. predictions from the four model types
pred1_50 <- subset(pred_50, pop_size == 100,)
pred1_95 <- subset(pred_95, pop_size == 100,)

# 2. 50% and 95% UDs for a range of population sizes
pred2_50 <- subset(pred_50, model == "MM3",)
pred2_95 <- subset(pred_95, model == "MM3",)

# 3. 50% and 95% UDs area occupied predictions for pop size of 100
p3_50 <- subset(pred2_50, pop_size == 100,)
p3_50$iso <- 50
p3_95 <- subset(pred2_95, pop_size == 100,)
p3_95$iso <- 95
pred3 <- rbind(p3_50, p3_95)

# 4. 50% UD area occupued versus representativeness based on inclusion
# incorporating birdlife predictions
pred_bl <- read.csv(paste0(in_direx, "birdlife_predictions.csv"))
pred4 <- data.frame(prediction = c(pred_bl$Rep70, pred_bl$Rep95, p3_50$prediction),
                    method = rep(c("Inclusion70", "Inclusion95", "Area Occupied"), each = length(pred_bl$Rep70)),
                    dataset = rep(pred_bl$Dataset, 3))


#### 2. AREA OCCUPIED - COMPARE 50% & 95% UD PREDICTIONS AMONG DIFFERENT ASYMPTOTIC MODEL TYPES ####

##### 2A. 50% UD ######

pred1_50$model <- as.factor(as.character(pred1_50$model))
pred1_50$dataset <- as.factor(as.character(pred1_50$dataset))

mod1_50 <- glmer(prediction ~ model + (1|dataset), data = pred1_50, family = "poisson")
summary(mod1_50)
mod1_50_2 <- update(mod1_50,~.- model)
lrtest(mod1_50, mod1_50_2)
# Model 1: prediction ~ model + (1 | dataset)
# Model 2: prediction ~ (1 | dataset)
# #Df  LogLik Df  Chisq Pr(>Chisq)    
# 1   5 -205.95                         
# 2   2 -292.58 -3 173.26  < 2.2e-16 ***
plot(mod1_50)

# tukey's tests for adjusted pairwise comparisons
summary(glht(mod1_50, mcp(model="Tukey")))
# Linear Hypotheses:
# Estimate Std. Error z value Pr(>|z|)    
# AR3 - AR2 == 0  0.20997    0.12944   1.622    0.360    
# MM2 - AR2 == 0  1.09861    0.11108   9.891   <1e-04 ***
# MM3 - AR2 == 0  1.00049    0.11250   8.893   <1e-04 ***
# MM2 - AR3 == 0  0.88864    0.10289   8.637   <1e-04 ***
# MM3 - AR3 == 0  0.79052    0.10442   7.571   <1e-04 ***
# MM3 - MM2 == 0 -0.09812    0.08054  -1.218    0.609  

# getting effect sizes
emmeans(mod1_50, pairwise ~model, adjust = "tukey")
# $emmeans
# model emmean    SE  df asymp.LCL asymp.UCL
# AR2     1.51 0.149 Inf      1.22      1.80
# AR3     1.72 0.143 Inf      1.44      2.00
# MM2     2.61 0.126 Inf      2.36      2.86
# MM3     2.51 0.128 Inf      2.26      2.76
exp(1.51); exp(0.149)  # AR2 - 4.526731 +- 1.160673
exp(1.72); exp(0.143)  # AR3 - 5.584528 +- 1.15373
exp(2.61); exp(0.126)  # MM2 - 13.59905 +- 1.134282
exp(2.51); exp(0.128)  # MM3 - 12.30493 +- 1.136553


##### 2B. 95% UD ######

pred1_95$model <- as.factor(as.character(pred1_95$model))
pred1_95$dataset <- as.factor(as.character(pred1_95$dataset))

mod1_95 <- glmer(prediction ~ model + (1|dataset), data = pred1_95, family = "poisson")
summary(mod1_95)
mod1_95_2 <- update(mod1_95,~.- model)
lrtest(mod1_95, mod1_95_2)
# Model 1: prediction ~ model + (1 | dataset)
# Model 2: prediction ~ (1 | dataset)
# #Df  LogLik Df  Chisq Pr(>Chisq)    
# 1   5 -224.93                         
# 2   2 -378.56 -3 307.25  < 2.2e-16 ***
plot(mod1_95)

summary(glht(mod1_95, mcp(model="Tukey")))
# Linear Hypotheses:
# Estimate Std. Error z value Pr(>|z|)    
# AR3 - AR2 == 0  0.38566    0.10448   3.691  0.00133 ** 
# MM2 - AR2 == 0  1.12654    0.09276  12.145  < 0.001 ***
# MM3 - AR2 == 0  1.23676    0.09156  13.507  < 0.001 ***
# MM2 - AR3 == 0  0.74088    0.08077   9.172  < 0.001 ***
# MM3 - AR3 == 0  0.85110    0.07940  10.719  < 0.001 ***
# MM3 - MM2 == 0  0.11022    0.06319   1.744  0.29398    

# getting effect sizes
emmeans(mod1_95, pairwise ~model, adjust = "tukey")
# model emmean     SE  df asymp.LCL asymp.UCL
# AR2     1.93 0.1119 Inf      1.71      2.15
# AR3     2.31 0.1022 Inf      2.11      2.51
# MM2     3.05 0.0902 Inf      2.88      3.23
# MM3     3.16 0.0889 Inf      2.99      3.34
exp(1.93); exp(0.1119)  # AR2 - 6.88951 +- 1.118401
exp(2.31); exp(0.1022)  # AR3 - 10.07442 +- 1.107605
exp(3.05); exp(0.0902)  # MM2 - 21.11534 +- 1.094393
exp(3.16); exp(0.0889)  # MM3 - 23.5706 +- 1.092971


#### 3. AREA OCCUPIED - COMPARE 50% & 95% UD PREDICTIONS FOR RANGE OF POP SIZES ####

##### 3A. 50% UD ######

pred2_50$pop <- as.factor(as.character(pred2_50$pop))
pred2_50$dataset <- as.factor(as.character(pred2_50$dataset))

mod2_50 <- glmer(prediction ~ pop + (1|dataset), data = pred2_50, family = "poisson")
summary(mod2_50)
# remove effect of population size and test significance
mod2_50_2 <- update(mod2_50,~.- pop)
lrtest(mod2_50, mod2_50_2)
# Model 1: prediction ~ pop + (1 | dataset)
# Model 2: prediction ~ (1 | dataset)
# #Df  LogLik Df Chisq Pr(>Chisq)  
# 1   6 -280.66                      
# 2   2 -285.56 -4  9.81    0.04375 *
plot(mod2_50)

# tukey's tests for adjusted pairwise comparisons
summary(glht(mod2_50, mcp(pop="Tukey")))
# Linear Hypotheses:
#                     Estimate Std. Error z value Pr(>|z|)  
# 1000 - 100 == 0     1.789e-01  7.921e-02   2.258   0.1585  
# 10000 - 100 == 0    2.016e-01  7.880e-02   2.558   0.0781 .
# 1e+05 - 100 == 0    2.016e-01  7.880e-02   2.558   0.0781 .
# 1e+06 - 100 == 0    2.016e-01  7.880e-02   2.558   0.0781 .
# 10000 - 1000 == 0   2.273e-02  7.517e-02   0.302   0.9982  
# 1e+05 - 1000 == 0   2.272e-02  7.517e-02   0.302   0.9982  
# 1e+06 - 1000 == 0   2.272e-02  7.517e-02   0.302   0.9982  
# 1e+05 - 10000 == 0 -4.240e-06  7.474e-02   0.000   1.0000  
# 1e+06 - 10000 == 0 -5.333e-06  7.474e-02   0.000   1.0000  
# 1e+06 - 1e+05 == 0 -1.093e-06  7.474e-02   0.000   1.0000  

# getting effect sizes
emmeans(mod2_50, pairwise ~pop, adjust = "tukey")
# $emmeans
# pop   emmean    SE  df asymp.LCL asymp.UCL
# 100     2.43 0.161 Inf      2.11      2.74
# 1000    2.61 0.160 Inf      2.29      2.92
# 10000   2.63 0.159 Inf      2.32      2.94
# 1e+05   2.63 0.159 Inf      2.32      2.94
# 1e+06   2.63 0.159 Inf      2.32      2.94
exp(2.43); exp(0.161)  # pop 100 = 11.35888 +- 1.174685
exp(2.61); exp(0.160)  # pop 1000 = 13.59905 +- 1.173511

# prediction for 1000 pop size
mean(pred2_50$prediction[pred2_50$pop == "1000"]) # 16.57143
range(pred2_50$prediction[pred2_50$pop == "1000"]) #  2 42



##### 3A. 95% UD ######

pred2_95$pop <- as.factor(as.character(pred2_95$pop))
pred2_95$dataset <- as.factor(as.character(pred2_95$dataset))

mod2_95 <- glmer(prediction ~ pop + (1|dataset), data = pred2_95, family = "poisson")
summary(mod2_95)
# remove effect of population size and test significance
mod2_95_2 <- update(mod2_95,~.- pop)
lrtest(mod2_95, mod2_95_2)
# Model 1: prediction ~ pop + (1 | dataset)
# Model 2: prediction ~ (1 | dataset)
# #Df  LogLik Df  Chisq Pr(>Chisq)    
# 1   6 -326.01                         
# 2   2 -355.76 -4 59.492  3.709e-12 ***
plot(mod2_95)

# tukey's tests for adjusted pairwise comparisons
summary(glht(mod2_95, mcp(pop="Tukey")))
# Linear Hypotheses:
# Estimate Std. Error z value Pr(>|z|)    
# 1000 - 100 == 0    0.314822   0.057214   5.503   <1e-04 ***
# 10000 - 100 == 0   0.355534   0.056730   6.267   <1e-04 ***
# 1e+05 - 100 == 0   0.358191   0.056699   6.317   <1e-04 ***
# 1e+06 - 100 == 0   0.359517   0.056684   6.343   <1e-04 ***
# 10000 - 1000 == 0  0.040712   0.052031   0.782    0.936    
# 1e+05 - 1000 == 0  0.043369   0.051997   0.834    0.920    
# 1e+06 - 1000 == 0  0.044695   0.051980   0.860    0.911    
# 1e+05 - 10000 == 0 0.002657   0.051465   0.052    1.000    
# 1e+06 - 10000 == 0 0.003983   0.051448   0.077    1.000    
# 1e+06 - 1e+05 == 0 0.001326   0.051414   0.026    1.000      

# getting effect sizes
emmeans(mod2_95, pairwise ~pop, adjust = "tukey")
# $emmeans
# pop   emmean    SE  df asymp.LCL asymp.UCL
# 100     3.12 0.112 Inf      2.90      3.34
# 1000    3.43 0.109 Inf      3.22      3.65
# 10000   3.47 0.109 Inf      3.26      3.69
# 1e+05   3.48 0.109 Inf      3.26      3.69
# 1e+06   3.48 0.109 Inf      3.26      3.69
exp(3.12); exp(0.112) # pop 100 = 22.64638 +- 1.118513
exp(3.43); exp(0.109) # pop 1000 = 30.87664 +- 1.115162

# prediction for 1000 pop size
mean(pred2_95$prediction[pred2_95$pop == "1000"]) # 34.38095
range(pred2_95$prediction[pred2_95$pop == "1000"]) #  13 61


#### 4. AREA OCCUPIED - COMPARE 50% & 95% UD PREDICTIONS FOR POP SIZE OF 100 ####

pred3$iso <- as.factor(as.character(pred3$iso))
pred3$dataset <- as.factor(as.character(pred3$dataset))

# running glmm
mod3 <- glmer(prediction ~ iso + (1|dataset), data = pred3, family = "poisson")
summary(mod3)
exp(0.59388)
# remove effect of isopleth and test significance
mod3_2 <- update(mod3,~.- iso)
lrtest(mod3, mod3_2)
# Model 1: prediction ~ iso + (1 | dataset)
# Model 2: prediction ~ (1 | dataset)
# #Df  LogLik Df  Chisq Pr(>Chisq)    
# 1   3 -154.90                         
# 2   2 -189.43 -1 69.066  < 2.2e-16 ***
plot(mod3)

# calculating mean and range for each isopleth category

# core
mean(pred3$prediction[pred3$iso == "50"]) # 13.85714
sd(pred3$prediction[pred3$iso == "50"]) # 7.445037
range(pred3$prediction[pred3$iso == "50"]) # 2 30

# general use
mean(pred3$prediction[pred3$iso == "95"]) # 25.09524
sd(pred3$prediction[pred3$iso == "95"]) #  8.898903
range(pred3$prediction[pred3$iso == "95"]) # 11 39

# number of datasets required for less than 10, 11-20, 21-30, 31-40

# 50% UD
nrow(pred3[pred3$prediction <= 10 & pred3$iso == "50",]) 
# 9
nrow(pred3[pred3$prediction <= 20 & pred3$prediction > 10 & pred3$iso == "50",]) 
# 8
nrow(pred3[pred3$prediction <= 30 & pred3$prediction > 20 & pred3$iso == "50",]) 
# 4 

# 95% UD
nrow(pred3[pred3$prediction <= 20 & pred3$prediction > 10 & pred3$iso == "95",]) 
# 8  
nrow(pred3[pred3$prediction <= 30 & pred3$prediction > 20 & pred3$iso == "95",]) 
# 6
nrow(pred3[pred3$prediction <= 40 & pred3$prediction > 30 & pred3$iso == "95",]) 
# 7

# effect sizes 
emmeans(mod3, pairwise ~iso, adjust = "tukey")
# $emmeans
# iso emmean     SE  df asymp.LCL asymp.UCL
# 50    2.59 0.0866 Inf      2.42      2.76
# 95    3.18 0.0773 Inf      3.03      3.33
exp(2.59); exp(0.0866) # 50% =  13.32977 +- 1.09046
exp(3.18); exp(0.0733) # 50% =  24.04675 +- 1.076053


#### 5. COMPARING PREDICTIONS FROM DIFFERENT METHODS - AREA OCCUPIED VS BIRDLIFE REPRESENTATITVENESS ####

pred4$method <- as.factor(as.character(pred4$method))
pred4$dataset <- as.factor(as.character(pred4$dataset))

mod4 <- glmer(prediction ~ method + (1|dataset), data = pred4, family = "poisson")
summary(mod4)
mod4_2 <- update(mod4,~.- method)
lrtest(mod4_2, mod4)
# Model 1: prediction ~ (1 | dataset)
# Model 2: prediction ~ method + (1 | dataset)
# #Df  LogLik Df  Chisq Pr(>Chisq)    
# 1   2 -450.49                         
# 2   4 -186.67  2 527.64  < 2.2e-16 ***
  
# tukey's tests for adjusted pairwise comparisons
summary(glht(mod4, mcp(method ="Tukey")))
# Linear Hypotheses:
#                           Estimate Std. Error z value            Pr(>|z|)    
#   Inclusion70 - Area Occupied == 0 -1.14072    0.11862  -9.617   <2e-16 ***
#   Inclusion95 - Area Occupied == 0  0.85610    0.06968  12.286   <2e-16 ***
#   Inclusion95 - Inclusion70 == 0    1.99682    0.11004  18.146   <2e-16 ***

mean(pred4$prediction[pred4$method == "Area Occupied"]) # 14
range(pred4$prediction[pred4$method == "Area Occupied"]) #  2 30

mean(pred4$prediction[pred4$method == "Inclusion70"]) # 4
range(pred4$prediction[pred4$method == "Inclusion70"]) #  1 10

mean(pred4$prediction[pred4$method == "Inclusion95"]) # 33
range(pred4$prediction[pred4$method == "Inclusion95"]) #  6  79

