library(ggplot2)
library(plyr)

#### 1. LOAD IN MINIMUM SAMPLE SIZES AND AMONG-INDIVIDUAL OVERLAP SCORES ####

# load in minimum sample sizes for gls data
pred_gls <- read.csv("./Data_outputs/area_occupied_50_predictions.csv")
pred_gls <- subset(pred_gls, pop_size == 100 & model == "MM3",)
pred_gls <- pred_gls[order(pred_gls$dataset),]
# load in overlap scores
overlap_gls <- read.csv("./Data_outputs/ba_overlap_scores.csv")
overlap_gls <- overlap_gls[order(overlap_gls$dataset),]
# load in PTT+GPS sample sizes and overlap scores
pred_br <- read.csv("./Data_inputs/ptt+gps_overlap_sample_predictions.csv")
# combine all datasets
all_preds <- data.frame(prediction = c(pred_gls$prediction, pred_br$prediction), ba = c(overlap_gls$ba, pred_br$ba),
                        species = c(overlap_gls$species, pred_br$species), year = c(overlap_gls$year, rep("ALL", nrow(pred_br))),
                        device = c(rep("GLS", nrow(overlap_gls)), pred_br$device), stage = c(rep("NB", nrow(overlap_gls)), pred_br$stage))


#### 2. EXAMINE RELATIONSHIPS BETWEEN THE TWO ####

# plotting the two
ggplot(all_preds, aes(ba, prediction)) + geom_point()

# using nls (non-linear least squares) to run model
x = all_preds$ba
y = all_preds$prediction
model <- nls(y ~ k/log(1-x), start = list(k = -10)) 
summary(model)
# Parameters:
#     Estimate Std. Error t value Pr(>|t|)    
# k  -3.9926     0.2263  -17.64   <2e-16 ***

# predicting and including on plot
fake_df <- data.frame(x = seq(0.05, 1, 0.005))
fake_df$y <- predict(model, fake_df)
ggplot(all_preds, aes(ba, prediction)) + 
  geom_point() + theme_bw() +   geom_line(data = fake_df, aes(x, y))


#### 3. PREDICTING ONTO BA VALUES FROM OPPEL ET AL. 2018 PAPER ####

# load in values
all_ba <- read.csv("./Data_inputs/oppel_et_al_2018_ba_scores.csv")
hist(all_ba$ba)
range(all_ba$ba) #  0.06 0.95
mean(all_ba$ba)
# 0.4130488

# defining model
pred_fun <- function(x) { y = -3.9926/log(1-x)
      return(y)}
# predicting
all_ba$prediction <- pred_fun(all_ba$ba) 

range(all_ba$prediction) #  1.332763 64.526448
mean(all_ba$prediction) # 12.2497


#### 4. SUMMARIZE PREDICTIONS FOR DIFFERENT TAXA ####

sum_family <- ddply(all_ba,.(family), summarize, mean_pred = mean(prediction), 
                    N = length(family), min_pred = min(prediction), max_pred = max(prediction))
sum_family <- sum_family[order(sum_family$mean_pred),]
sum_family
#               family mean_pred  N  min_pred  max_pred
# 7  Phalacrocoracidae  3.769221  6  1.501395  7.102768
# 5            Laridae  4.820900  8  3.225368  8.077341
# 1            Alcidae  5.081296  8  2.797670  7.567006
# 10           Sulidae  6.859456  5  1.332763 13.878515
# 9       Spheniscidae  8.522902 13  2.636898 24.566984
# 4       Hydrobatidae  9.344078  5  3.601276 18.947326
# 6      Phaethontidae 11.527889  3 10.352576 13.878515
# 8     Procellariidae 17.713615 20  3.803120 55.016699
# 3         Fregatidae 18.240050  3  9.608788 31.232846
# 2        Diomedeidae 24.295508 11  7.815974 64.526448

# proportion datasets with preditions below 10
nrow(all_ba[all_ba$prediction<10,]) # 47
nrow(all_ba[all_ba$prediction<10,])/nrow(all_ba)*100 # 57.31707

# proportion datasets with preditions below 20
nrow(all_ba[all_ba$prediction<20,]) # 69
nrow(all_ba[all_ba$prediction<20,])/nrow(all_ba)*100 # 84.14634

# proportion datasets with preditions below 30
nrow(all_ba[all_ba$prediction<30,]) # 74
nrow(all_ba[all_ba$prediction<30,])/nrow(all_ba)*100 # 90.2439

# proportion above 50
nrow(all_ba[all_ba$prediction>50,]) # 2
nrow(all_ba[all_ba$prediction>50,])/nrow(all_ba)*100 # 2.439024
all_ba[all_ba$prediction>50,]
