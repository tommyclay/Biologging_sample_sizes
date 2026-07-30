library(lme4)
library(lmtest)
library(multcomp)
library(emmeans)

in_direx <- "./Data_outputs/"

#### 1. LOAD IN PREDICTIONS ####

# load in area occupied file
pred_50_gls <- read.csv(paste0(in_direx, "area_occupied_50_predictions.csv"))
pred_50_ptt <- read.csv(paste0(in_direx, "gps+ptt_area_occupied_50_predictions.csv"))
pred_95_gls <- read.csv(paste0(in_direx, "area_occupied_95_predictions.csv"))
pred_95_ptt <- read.csv(paste0(in_direx, "gps+ptt_area_occupied_95_predictions.csv"))
pred_50_both <- rbind(pred_50_gls, pred_50_ptt)
pred_95_both <- rbind(pred_95_gls, pred_95_ptt)

# 1. predictions from the four model types
pred1_50 <- subset(pred_50_both, pop_size == 100,)
pred1_95 <- subset(pred_95_both, pop_size == 100,)

# 2. 50% and 95% UDs for a range of population sizes
pred2_50 <- subset(pred_50_both, model == "MM3",)
pred2_95 <- subset(pred_95_both, model == "MM3",)

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
# give ID for non-breeding versus breeding
pred1_50$season <- as.factor(as.character(ifelse(grepl("GPS|PTT", pred1_50$dataset),"breeding", "non-breeding")))

mod1_50 <- glmer(prediction ~ model + season + (1|dataset), data = pred1_50, family = "poisson")
summary(mod1_50)
mod1_50_2 <- update(mod1_50,~.- model)
lrtest(mod1_50, mod1_50_2)
# Model 1: prediction ~ model + (1 | dataset)
# Model 2: prediction ~ (1 | dataset)
# #Df  LogLik Df  Chisq Pr(>Chisq)    
# 1   5 -354.84                         
# 2   2 -508.58 -3 307.48  < 2.2e-16 ***
mod1_50_3 <- update(mod1_50,~.- season)
lrtest(mod1_50, mod1_50_3)
# Model 1: prediction ~ model + season + (1 | dataset)
# Model 2: prediction ~ model + (1 | dataset)
# #Df  LogLik Df  Chisq Pr(>Chisq)    
# 1   6 -347.64                         
# 2   5 -354.84 -1 14.403  0.0001475 ***
plot(mod1_50)

# tukey's tests for adjusted pairwise comparisons
summary(glht(mod1_50, mcp(model="Tukey")))
# Linear Hypotheses:
# Estimate Std. Error z value Pr(>|z|)    
# AR3 - AR2 == 0 0.308860   0.081847   3.774   <0.001 ***
# MM2 - AR2 == 0 0.967118   0.073014  13.246   <0.001 ***
# MM3 - AR2 == 0 0.973018   0.072955  13.337   <0.001 ***
# MM2 - AR3 == 0 0.658258   0.065610  10.033   <0.001 ***
# MM3 - AR3 == 0 0.664158   0.065545  10.133   <0.001 ***
# MM3 - MM2 == 0 0.005901   0.054114   0.109        1    

# getting effect sizes
emmeans(mod1_50, pairwise ~model, adjust = "tukey")
# $emmeans
# model emmean    SE  df asymp.LCL asymp.UCL
# AR2     1.99 0.118 Inf      1.76      2.22
# AR3     2.30 0.114 Inf      2.07      2.52
# MM2     2.95 0.107 Inf      2.74      3.16
# MM3     2.96 0.107 Inf      2.75      3.17
exp(1.99); exp(0.118)  # AR2 - 7.315534 +- 1.125244
exp(2.30); exp(0.114)  # AR3 - 9.974182 +- 1.120752
exp(2.95); exp(0.107)  # MM2 - 19.10595 +- 1.112934
exp(2.96); exp(0.107)  # MM3 - 19.29797 +- 1.112934

emmeans(mod1_50, pairwise ~season, adjust = "tukey")
# $emmeans
# season       emmean    SE  df asymp.LCL asymp.UCL
# breeding       2.99 0.165 Inf      2.66      3.31
# non-breeding   2.11 0.123 Inf      1.87      2.35
exp(2.99); exp(0.165)  # breeding - 19.88568 +- 1.179393
exp(2.11); exp(0.123)  # non-breeding - 8.248241 +- 1.130884


##### 2B. 95% UD ######

pred1_95$model <- as.factor(as.character(pred1_95$model))
pred1_95$dataset <- as.factor(as.character(pred1_95$dataset))
# give ID for non-breeding versus breeding
pred1_95$season <- as.factor(as.character(ifelse(grepl("GPS|PTT", pred1_95$dataset),"breeding", "non-breeding")))
                             
mod1_95 <- glmer(prediction ~ model + season + (1|dataset), data = pred1_95, family = "poisson")
summary(mod1_95)
mod1_95_2 <- update(mod1_95,~.- model)
lrtest(mod1_95, mod1_95_2)
# Model 1: prediction ~ model + (1 | dataset)
# Model 2: prediction ~ (1 | dataset)
# #Df  LogLik Df  Chisq Pr(>Chisq)    
# 1   6 -377.35                         
# 2   3 -628.63 -3 502.55  < 2.2e-16 ***
plot(mod1_95)

summary(glht(mod1_95, mcp(model="Tukey")))
# Linear Hypotheses:
#    Estimate Std. Error z value Pr(>|z|)    
# AR3 - AR2 == 0  0.34098    0.06440   5.295   <0.001 ***
# MM2 - AR2 == 0  0.92934    0.05814  15.984   <0.001 ***
# MM3 - AR2 == 0  1.03239    0.05733  18.007   <0.001 ***
# MM2 - AR3 == 0  0.58836    0.05177  11.364   <0.001 ***
# MM3 - AR3 == 0  0.69141    0.05086  13.595   <0.001 ***
# MM3 - MM2 == 0  0.10305    0.04266   2.415   0.0725 .  

# getting effect sizes
emmeans(mod1_95, pairwise ~model, adjust = "tukey")
# model emmean     SE  df asymp.LCL asymp.UCL
# AR2     2.56 0.0746 Inf      2.41      2.70
# AR3     2.90 0.0697 Inf      2.76      3.03
# MM2     3.49 0.0640 Inf      3.36      3.61
# MM3     3.59 0.0632 Inf      3.47      3.71
exp(2.56); exp(0.0746)  # AR2 - 12.93582 +- 1.077453
exp(2.90); exp(0.0697)  # AR3 - 18.17415 +- 1.072186
exp(3.49); exp(0.0649)  # MM2 - 32.78595 +- 1.067052
exp(3.59); exp(0.0632)  # MM3 - 36.23408 +- 1.06524

emmeans(mod1_95, pairwise ~season, adjust = "tukey")
# $emmeans
# season       emmean     SE  df asymp.LCL asymp.UCL
# breeding       3.61 0.0937 Inf      3.43      3.79
# non-breeding   2.66 0.0717 Inf      2.52      2.80
exp(3.61); exp(0.0937)  # breeding - 36.96605 +- 1.09823
exp(2.66); exp(0.0717)  # non-breeding - 14.29629 +- 1.074333



#### 3. AREA OCCUPIED - COMPARE 50% & 95% UD PREDICTIONS FOR RANGE OF POP SIZES ####

##### 3A. 50% UD - JUST NON-BREEDING ######

pred2_50$pop <- as.factor(as.character(pred2_50$pop))
pred2_50$dataset <- as.factor(as.character(pred2_50$dataset))
# give ID for non-breeding versus breeding
pred2_50$season <- as.factor(as.character(ifelse(grepl("GPS|PTT", pred2_50$dataset),"breeding", "non-breeding")))
pred2_50 <- subset(pred2_50, season == "non-breeding",)

mod2_50 <- glmer(prediction ~ pop  + (1|dataset), data = pred2_50, family = "poisson")
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
#   Estimate Std. Error z value Pr(>|z|)  
# 1000 - 100 == 0     1.789e-01  7.921e-02   2.258   0.1584  
# 10000 - 100 == 0    2.016e-01  7.880e-02   2.558   0.0781 .
# 1e+05 - 100 == 0    2.016e-01  7.880e-02   2.558   0.0782 .
# 1e+06 - 100 == 0    2.016e-01  7.880e-02   2.558   0.0782 .
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
range(pred2_50$prediction[pred2_50$pop == "1000"]) # 2 42



##### 3A. 95% UD ######

pred2_95$pop <- as.factor(as.character(pred2_95$pop))
pred2_95$dataset <- as.factor(as.character(pred2_95$dataset))
# give ID for non-breeding versus breeding
pred2_95$season <- as.factor(as.character(ifelse(grepl("GPS|PTT", pred2_95$dataset),"breeding", "non-breeding")))
pred2_95 <- subset(pred2_95, season == "non-breeding",)

mod2_95 <- glmer(prediction ~ pop + (1|dataset), data = pred2_95, family = "poisson")
summary(mod2_95)
# remove effect of population size and test significance
mod2_95_2 <- update(mod2_95,~.- pop)
lrtest(mod2_95, mod2_95_2)
# Model 1: prediction ~ pop + (1 | dataset)
# Model 2: prediction ~ (1 | dataset)
# #Df  LogLik Df  Chisq Pr(>Chisq)    
# 1    6 -326.01                         
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
range(pred2_95$prediction[pred2_95$pop == "1000"]) # 13 61


#### 4. AREA OCCUPIED - COMPARE 50% & 95% UD PREDICTIONS FOR POP SIZE OF 100 ####

pred3$iso <- as.factor(as.character(pred3$iso))
pred3$dataset <- as.factor(as.character(pred3$dataset))
# give ID for non-breeding versus breeding
pred3$season <- as.factor(as.character(ifelse(grepl("GPS|PTT", pred3$dataset),"breeding", "non-breeding")))

# running glmm
mod3 <- glmer(prediction ~ iso + season + (1|dataset), data = pred3, family = "poisson")
summary(mod3)
exp(0.52889) # 1.697048
# remove effect of isopleth and test significance
mod3_2 <- update(mod3,~.- iso)
lrtest(mod3, mod3_2)
# Model 1: prediction ~ iso + (1 | dataset)
# Model 2: prediction ~ (1 | dataset)
# #Df  LogLik Df  Chisq Pr(>Chisq)    
# 1   3 -280.81                         
# 2   2 -342.76 -1 123.91  < 2.2e-16 ***
mod3_3 <- update(mod3,~.- season)
lrtest(mod3, mod3_3)
#Df  LogLik Df Chisq Pr(>Chisq)    
# 1   4 -263.34                        
# 2   3 -280.81 -1 34.94    3.4e-09 ***
plot(mod3)

# calculating mean and range for each isopleth category

# core - non-breeding
mean(pred3$prediction[pred3$iso == "50" & pred3$season == "non-breeding"]) # 13.85714
sd(pred3$prediction[pred3$iso == "50" & pred3$season == "non-breeding"]) #  7.445037
range(pred3$prediction[pred3$iso == "50" & pred3$season == "non-breeding"]) # 2 30

# core - breeding
mean(pred3$prediction[pred3$iso == "50" & pred3$season == "breeding"]) # 35.36364
sd(pred3$prediction[pred3$iso == "50" & pred3$season == "breeding"]) # 19.81047
range(pred3$prediction[pred3$iso == "50" & pred3$season == "breeding"]) # 11 67

# general use - non-breeding
mean(pred3$prediction[pred3$iso == "95"& pred3$season == "non-breeding"]) # 25.09524
sd(pred3$prediction[pred3$iso == "95"& pred3$season == "non-breeding"]) #  8.898903
range(pred3$prediction[pred3$iso == "95"& pred3$season == "non-breeding"]) # 11 39

# general use - breeding
mean(pred3$prediction[pred3$iso == "95"& pred3$season == "breeding"]) # 57
sd(pred3$prediction[pred3$iso == "95"& pred3$season == "breeding"]) #  8.318654
range(pred3$prediction[pred3$iso == "95"& pred3$season == "breeding"]) # 46 66


# number of datasets required for less than 10, 11-20, 21-30, 31-40

nrow(pred3) # 64

# 50% UD
nrow(pred3[pred3$iso == "50",]) # 32
nrow(pred3[pred3$prediction <= 10 & pred3$iso == "50",]) 
# 9
nrow(pred3[pred3$prediction <= 20 & pred3$prediction > 10 & pred3$iso == "50",]) 
# 12
nrow(pred3[pred3$prediction <= 30 & pred3$prediction > 20 & pred3$iso == "50",]) 
# 6

# 95% UD
nrow(pred3[pred3$iso == "95",]) # 32
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
# 50    3.07 0.0621 Inf      2.95      3.19
# 95    3.60 0.0570 Inf      3.49      3.71
exp(3.07); exp(0.0621) # 50% =  21.5419 +- 1.064069
exp(3.60); exp(0.0570) # 50% =  36.59823 +- 1.058656

emmeans(mod3, pairwise ~season+iso, adjust = "tukey")
# season       iso emmean     SE  df asymp.LCL asymp.UCL
# breeding     50    3.51 0.0899 Inf      3.33      3.68
# non-breeding 50    2.63 0.0738 Inf      2.49      2.78
# breeding     95    4.04 0.0865 Inf      3.87      4.21
# non-breeding 95    3.16 0.0696 Inf      3.03      3.30
exp(3.51); exp(0.0899) # 50% breeding =  33.44827 +- 1.094065
exp(2.63); exp(0.0738) # 50% non-breeding =  13.87377 +- 1.076591
exp(4.04); exp(0.0865) # 95% breeding =  56.82634 +- 1.090351
exp(3.16); exp(0.0696) # 95% non-breeding =  23.5706 +- 1.072079



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
# 1   2 -919.16                         
# 2   4 -369.02  2 1100.3  < 2.2e-16 ***
  
# tukey's tests for adjusted pairwise comparisons
summary(glht(mod4, mcp(method ="Tukey")))
# Linear Hypotheses:
#                           Estimate Std. Error z value            Pr(>|z|)    
# Inclusion70 - Area Occupied == 0 -1.26460    0.08147  -15.52   <2e-16 ***
# Inclusion95 - Area Occupied == 0  0.75239    0.04637   16.23   <2e-16 ***
# Inclusion95 - Inclusion70 == 0    2.01698    0.07658   26.34   <2e-16 ***

pred4$season <- as.factor(as.character(ifelse(grepl("GPS|PTT", pred4$dataset),"breeding", "non-breeding")))

mean(pred4$prediction[pred4$method == "Area Occupied"]) # 21.25
range(pred4$prediction[pred4$method == "Area Occupied"]) #  2 67

mean(pred4$prediction[pred4$method == "Inclusion70" & pred4$season == "non-breeding"]) # 4.428571
range(pred4$prediction[pred4$method == "Inclusion70" & pred4$season == "non-breeding"]) #  1 10

mean(pred4$prediction[pred4$method == "Inclusion70" & pred4$season == "breeding"]) # 9
range(pred4$prediction[pred4$method == "Inclusion70" & pred4$season == "breeding"]) #  3 22

mean(pred4$prediction[pred4$method == "Inclusion95"& pred4$season == "non-breeding"]) # 32.61905
range(pred4$prediction[pred4$method == "Inclusion95"& pred4$season == "non-breeding"]) #  6 79

mean(pred4$prediction[pred4$method == "Inclusion95"& pred4$season == "breeding"]) # 68.90909
range(pred4$prediction[pred4$method == "Inclusion95"& pred4$season == "breeding"]) # 21 175
