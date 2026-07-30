library(ggplot2)
library(plyr)

#### 1. LOAD FUNCTIONS #####

source("./Codes/1_predicting_minimum_sample_sizes/functions.R")

#### 2. LOAD IN MINIMUM SAMPLE SIZES AND AMONG-INDIVIDUAL OVERLAP SCORES ####

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


#### 3. EXAMINE RELATIONSHIPS BETWEEN THE TWO ####

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

# using function from functions.R code to predict
ba_preds <- prediction_from_ba(all_ba$ba) 
all_ba$prediction <- ba_preds$prediction
all_ba$lower <- ba_preds$lower
all_ba$upper <- ba_preds$upper

range(all_ba$prediction) #  1 65
mean(all_ba$prediction) # 12.2439


#### 4. SUMMARIZE PREDICTIONS FOR DIFFERENT TAXA ####

sum_family <- ddply(all_ba,.(family), summarize, mean_pred = round(mean(prediction)), 
                    N = length(family), min_pred = min(prediction), max_pred = max(prediction))
sum_family <- sum_family[order(sum_family$mean_pred),]
sum_family
#                 family mean_pred N min_pred max_pred
# 7  Phalacrocoracidae         4  6        2        7
# 1            Alcidae         5  8        3        8
# 5            Laridae         5  8        3        8
# 10           Sulidae         7  5        1       14
# 9       Spheniscidae         8 13        3       25
# 4       Hydrobatidae        10  5        4       19
# 6      Phaethontidae        11  3       10       14
# 3         Fregatidae        18  3       10       31
# 8     Procellariidae        18 20        4       55
# 2        Diomedeidae        24 11        8       65

# proportion datasets with preditions below 10
nrow(all_ba[all_ba$prediction<10,]) # 43
nrow(all_ba[all_ba$prediction<10,])/nrow(all_ba)*100 # 52.43902

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
