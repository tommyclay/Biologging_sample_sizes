library(ggplot2)
library(plyr)
library(cowplot)

# setting up directories
resample_direx <- "./Data_outputs/resampled_area_outputs/gls_data/"
model_direx <- "./Data_outputs/resampled_area_outputs/gls_models/"
plot_direx <- "./Plots/"

# loading in minimium sample sizes from models
min_sample <- read.csv("./Data_outputs/area_occupied_50_predictions.csv")
min_sample <- subset(min_sample, model == "MM3" & pop_size == 100,)

#### 1. PLOTTING FOR EACH YEAR COMBINATION ####

##### 1A. AP 2010 #####

# loading in data
ap10 <- read.csv(paste0(resample_direx, "AP10_50.csv"))
load(paste0(model_direx, "AP10_50.Rdata"))
# predicting from model
ap10_pred <- data.frame(sample = seq(1, max(ap10$sample), 1))
ap10_pred2 <- predict(mod, ap10_pred, se.fit = TRUE)
ap10_pred$pred <- ap10_pred2[,1]
# getting median values
ap10_med <- ddply(ap10,.(sample), summarize,  tot_area_med = median(tot_area))
# picking out the first 100 iterations
ap10 <- subset(ap10, iteration <= 100,)

ap10_p <- ggplot() + geom_line(data = ap10, aes(x = sample, y = tot_area, group = as.factor(iteration)), linewidth = 0.3, col = "grey80") +
  geom_hline(yintercept=mod$coefficients[2]*0.95, linetype="dashed", linewidth = 0.3)+
  geom_hline(yintercept=mod$coefficients[2], linetype="solid", linewidth = 0.3)+
  geom_line(data = ap10_pred, aes(x= sample, y = pred), colour = "blue") + theme_bw() +
  geom_point(data = ap10_med, aes(x = sample, y = tot_area_med), pch = 21, size = 1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(t = 0.05, r = 0.01,b = -0.02, l = 0), "cm"),
        axis.text = element_text(size = 7), axis.title.y = element_blank(),axis.title.x = element_blank()) +
  ylab(expression("HR area (10"^6~km^2~")"))+
  scale_y_continuous(limits = c(0, 1.45), breaks = seq(0, 1.25, 0.25), expand = c(0, 0))+
  scale_x_continuous(breaks = seq(0, max(ap10$sample), 1)) +
  geom_vline(xintercept=min_sample$prediction[min_sample$dataset == "AP10"], linetype="dotted")
ap10_p


##### 1B. AP 2011 #####

ap11 <- read.csv(paste0(resample_direx, "AP11_50.csv"))
load(paste0(model_direx, "AP11_50.Rdata"))
ap11_pred <- data.frame(sample = seq(1, max(ap11$sample), 1))
ap11_pred2 <- predict(mod, ap11_pred, se.fit = TRUE)
ap11_pred$pred <- ap11_pred2[,1]
ap11_med <- ddply(ap11,.(sample), summarize,  tot_area_med = median(tot_area))
ap11 <- subset(ap11, iteration <= 100,)

ap11_p <- ggplot() + geom_line(data = ap11, aes(x = sample, y = tot_area, group = as.factor(iteration)), linewidth = 0.3, col = "grey80") +
  geom_hline(yintercept=mod$coefficients[2]*0.95, linetype="dashed", size = 0.3)+
  geom_hline(yintercept=mod$coefficients[2], linetype="solid", size = 0.3)+
  geom_line(data = ap11_pred, aes(x= sample, y = pred), colour = "blue") + theme_bw() +
  geom_point(data = ap11_med, aes(x = sample, y = tot_area_med), pch = 21, size = 1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(t = 0.05, r = 0.01,b = -0.02, l = 0), "cm"),
        axis.text = element_text(size = 7), axis.title.y = element_blank(),axis.title.x = element_blank()) +
  ylab(expression("HR area (10"^6~km^2~")"))+
  scale_y_continuous(limits = c(0, 1.45), breaks = seq(0, 1.25, 0.25), expand = c(0, 0))+
  scale_x_continuous(breaks = seq(0, max(ap11$sample), 1)) 
  #geom_vline(xintercept=min_sample$prediction[min_sample$dataset == "AP11"], linetype="dotted")
ap11_p

##### 1C. AP BOTH #####

apboth <- read.csv(paste0(resample_direx, "APBOTH_50.csv"))
load(paste0(model_direx, "APBOTH_50.Rdata"))
apboth_pred <- data.frame(sample = seq(1, max(apboth$sample), 1))
apboth_pred2 <- predict(mod, apboth_pred, se.fit = TRUE)
apboth_pred$pred <- apboth_pred2[,1]
apboth_med <- ddply(apboth,.(sample), summarize,  tot_area_med = median(tot_area))
apboth <- subset(apboth, iteration <= 100,)

apboth_p <- ggplot() + geom_line(data = apboth, aes(x = sample, y = tot_area, group = as.factor(iteration)), linewidth = 0.3, col = "grey80") +
  geom_hline(yintercept=mod$coefficients[2]*0.95, linetype="dashed", size = 0.3)+
  geom_hline(yintercept=mod$coefficients[2], linetype="solid", size = 0.3)+
  geom_line(data = apboth_pred, aes(x= sample, y = pred), colour = "blue") + theme_bw() +
  geom_point(data = apboth_med, aes(x = sample, y = tot_area_med), pch = 21, size = 1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(t = 0.05, r = 0.01,b = -0.02, l = 0), "cm"),
        axis.text = element_text(size = 7), axis.title.y = element_blank(),axis.title.x = element_blank()) +
  ylab(expression("HR area (10"^6~km^2~")"))+
  scale_y_continuous(limits = c(0, 1.45), breaks = seq(0, 1.25, 0.25), expand = c(0, 0))+
  scale_x_continuous(breaks = seq(0, max(apboth$sample), 2)) +
  geom_vline(xintercept=min_sample$prediction[min_sample$dataset == "APBOTH"], linetype="dotted")
apboth_p


##### 1D. BBA 2002 #####

bba02 <- read.csv(paste0(resample_direx, "BBA02_50.csv"))
load(paste0(model_direx, "BBA02_50.Rdata"))
bba02_pred <- data.frame(sample = seq(1, max(bba02$sample), 1))
bba02_pred2 <- predict(mod, bba02_pred, se.fit = TRUE)
bba02_pred$pred <- bba02_pred2[,1]
bba02_med <- ddply(bba02,.(sample), summarize,  tot_area_med = median(tot_area))
bba02 <- subset(bba02, iteration <= 100,)

bba02_p <- ggplot() + geom_line(data = bba02, aes(x = sample, y = tot_area, group = as.factor(iteration)), linewidth = 0.3, col = "grey80") +
  geom_hline(yintercept=mod$coefficients[2]*0.95, linetype="dashed", size = 0.3)+
  geom_hline(yintercept=mod$coefficients[2], linetype="solid", size = 0.3)+
  geom_line(data = bba02_pred, aes(x= sample, y = pred), colour = "blue") + theme_bw() +
  geom_point(data = bba02_med, aes(x = sample, y = tot_area_med), pch = 21, size = 1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(t = 0.05, r = 0.01,b = -0.02, l = 0), "cm"),
        axis.text = element_text(size = 7), axis.title.y = element_blank(),axis.title.x = element_blank()) +
  ylab(expression("HR area (10"^6~km^2~")"))+
  scale_y_continuous(limits = c(0, 3.6), breaks = seq(0, 3.5, 0.5), expand = c(0, 0))+
  scale_x_continuous(breaks = seq(0, max(bba02$sample), 4)) + 
  geom_vline(xintercept=min_sample$prediction[min_sample$dataset == "BBA02"], linetype="dotted")
bba02_p


##### 1E. BBA 2003 #####

bba03 <- read.csv(paste0(resample_direx, "BBA03_50.csv"))
load(paste0(model_direx, "BBA03_50.Rdata"))
bba03_pred <- data.frame(sample = seq(1, max(bba03$sample), 1))
bba03_pred2 <- predict(mod, bba03_pred, se.fit = TRUE)
bba03_pred$pred <- bba03_pred2[,1]
bba03_med <- ddply(bba03,.(sample), summarize,  tot_area_med = median(tot_area))
bba03 <- subset(bba03, iteration <= 100,)

bba03_p <- ggplot() + geom_line(data = bba03, aes(x = sample, y = tot_area, group = as.factor(iteration)), linewidth = 0.3, col = "grey80") +
  geom_hline(yintercept=mod$coefficients[2]*0.95, linetype="dashed", size = 0.3)+
  geom_hline(yintercept=mod$coefficients[2], linetype="solid", size = 0.3)+
  geom_line(data = bba03_pred, aes(x= sample, y = pred), colour = "blue") + theme_bw() +
  geom_point(data = bba03_med, aes(x = sample, y = tot_area_med), pch = 21, size = 1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(t = 0.05, r = 0.01,b = -0.02, l = 0), "cm"),
        axis.text = element_text(size = 7), axis.title.y = element_blank(),axis.title.x = element_blank()) +
  ylab(expression("HR area (10"^6~km^2~")"))+
  scale_y_continuous(limits = c(0, 3.6), breaks = seq(0, 3.5, 0.5), expand = c(0, 0))+
  scale_x_continuous(breaks = seq(0, max(bba03$sample), 4)) + 
  geom_vline(xintercept=min_sample$prediction[min_sample$dataset == "BBA03"], linetype="dotted")
bba03_p


##### 1F. BP 2011 #####

bp11 <- read.csv(paste0(resample_direx, "BP11_50.csv"))
load(paste0(model_direx, "BP11_50.Rdata"))
bp11_pred <- data.frame(sample = seq(1, max(bp11$sample), 1))
bp11_pred2 <- predict(mod, bp11_pred, se.fit = TRUE)
bp11_pred$pred <- bp11_pred2[,1]
bp11_med <- ddply(bp11,.(sample), summarize,  tot_area_med = median(tot_area))
bp11 <- subset(bp11, iteration <= 100,)

bp11_p <- ggplot() + geom_line(data = bp11, aes(x = sample, y = tot_area, group = as.factor(iteration)), linewidth = 0.3, col = "grey80") +
  geom_hline(yintercept=mod$coefficients[2]*0.95, linetype="dashed", size = 0.3)+
  geom_hline(yintercept=mod$coefficients[2], linetype="solid", size = 0.3)+
  geom_line(data = bp11_pred, aes(x= sample, y = pred), colour = "blue") + theme_bw() +
  geom_point(data = bp11_med, aes(x = sample, y = tot_area_med), pch = 21, size = 1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(t = 0.05, r = 0.01,b = -0.02, l = 0), "cm"),
        axis.text = element_text(size = 7), axis.title.y = element_blank(),axis.title.x = element_blank()) +
  ylab(expression("HR area (10"^6~km^2~")"))+
  scale_y_continuous(limits = c(0, 4.7), breaks = seq(0, 4, 1), expand = c(0, 0))+
  scale_x_continuous(breaks = seq(0, max(bp11$sample), 1))  + 
  geom_vline(xintercept=min_sample$prediction[min_sample$dataset == "BP11"], linetype="dotted")
bp11_p


##### 1G. BS 2002 #####

bs02 <- read.csv(paste0(resample_direx, "BS02_50.csv"))
load(paste0(model_direx, "BS02_50.Rdata"))
bs02_pred <- data.frame(sample = seq(1, max(bs02$sample), 1))
bs02_pred2 <- predict(mod, bs02_pred, se.fit = TRUE)
bs02_pred$pred <- bs02_pred2[,1]
bs02_med <- ddply(bs02,.(sample), summarize,  tot_area_med = median(tot_area))
bs02 <- subset(bs02, iteration <= 100,)

bs02_p <- ggplot() + geom_line(data = bs02, aes(x = sample, y = tot_area, group = as.factor(iteration)), linewidth = 0.3, col = "grey80") +
  geom_hline(yintercept=mod$coefficients[2]*0.95, linetype="dashed", size = 0.3)+
  geom_hline(yintercept=mod$coefficients[2], linetype="solid", size = 0.3)+
  geom_line(data = bs02_pred, aes(x= sample, y = pred), colour = "blue") + theme_bw() +
  geom_point(data = bs02_med, aes(x = sample, y = tot_area_med), pch = 21, size = 1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(t = 0.05, r = 0.01,b = -0.02, l = 0), "cm"),
        axis.text = element_text(size = 7), axis.title.y = element_blank(),axis.title.x = element_blank()) +
  ylab(expression("HR area (10"^6~km^2~")"))+
  scale_y_continuous(limits = c(0, 3.85), breaks = seq(0, 3.5, 0.5), expand = c(0, 0))+
  scale_x_continuous(breaks = seq(0, max(bs02$sample), 1)) 
  #geom_vline(xintercept=min_sample$prediction[min_sample$dataset == "BS02"], linetype="dotted")
bs02_p



##### 1H. BS 2012 #####

bs12 <- read.csv(paste0(resample_direx, "BS12_50.csv"))
load(paste0(model_direx, "BS12_50.Rdata"))
bs12_pred <- data.frame(sample = seq(1, max(bs12$sample), 1))
bs12_pred2 <- predict(mod, bs12_pred, se.fit = TRUE)
bs12_pred$pred <- bs12_pred2[,1]
bs12_med <- ddply(bs12,.(sample), summarize,  tot_area_med = median(tot_area))
bs12 <- subset(bs12, iteration <= 100,)

bs12_p <- ggplot() + geom_line(data = bs12, aes(x = sample, y = tot_area, group = as.factor(iteration)), linewidth = 0.3, col = "grey80") +
  geom_hline(yintercept=mod$coefficients[2]*0.95, linetype="dashed", size = 0.3)+
  geom_hline(yintercept=mod$coefficients[2], linetype="solid", size = 0.3)+
  geom_line(data = bs12_pred, aes(x= sample, y = pred), colour = "blue") + theme_bw() +
  geom_point(data = bs12_med, aes(x = sample, y = tot_area_med), pch = 21, size = 1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(t = 0.05, r = 0.01,b = -0.02, l = 0), "cm"),
        axis.text = element_text(size = 7), axis.title.y = element_blank(),axis.title.x = element_blank()) +
  ylab(expression("HR area (10"^6~km^2~")"))+
  scale_y_continuous(limits = c(0, 3.85), breaks = seq(0, 3.5, 0.5), expand = c(0, 0))+
  scale_x_continuous(breaks = seq(0, max(bs12$sample), 2)) + 
  geom_vline(xintercept=min_sample$prediction[min_sample$dataset == "BS12"], linetype="dotted")
bs12_p



##### 1I. BS BOTH #####

bsboth <- read.csv(paste0(resample_direx, "BSBOTH_50.csv"))
load(paste0(model_direx, "BSBOTH_50.Rdata"))
bsboth_pred <- data.frame(sample = seq(1, max(bsboth$sample), 1))
bsboth_pred2 <- predict(mod, bsboth_pred, se.fit = TRUE)
bsboth_pred$pred <- bsboth_pred2[,1]
bsboth_med <- ddply(bsboth,.(sample), summarize,  tot_area_med = median(tot_area))
bsboth <- subset(bsboth, iteration <= 100,)

bsboth_p <- ggplot() + geom_line(data = bsboth, aes(x = sample, y = tot_area, group = as.factor(iteration)), linewidth = 0.3, col = "grey80") +
  geom_hline(yintercept=mod$coefficients[2]*0.95, linetype="dashed", size = 0.3)+
  geom_hline(yintercept=mod$coefficients[2], linetype="solid", size = 0.3)+
  geom_line(data = bsboth_pred, aes(x= sample, y = pred), colour = "blue") + theme_bw() +
  geom_point(data = bsboth_med, aes(x = sample, y = tot_area_med), pch = 21, size = 1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(t = 0.05, r = 0.01,b = -0.02, l = 0), "cm"),
        axis.text = element_text(size = 7), axis.title.y = element_blank(),axis.title.x = element_blank()) +
  ylab(expression("HR area (10"^6~km^2~")"))+
  scale_y_continuous(limits = c(0, 3.85), breaks = seq(0, 3.5, 0.5), expand = c(0, 0))+
  scale_x_continuous(breaks = seq(0, max(bsboth$sample), 4)) + 
  geom_vline(xintercept=min_sample$prediction[min_sample$dataset == "BSBOTH"], linetype="dotted")
bsboth_p



##### 1I. GHA 1999 #####

gha99 <- read.csv(paste0(resample_direx, "GHA99_50.csv"))
load(paste0(model_direx, "GHA99_50.Rdata"))
gha99_pred <- data.frame(sample = seq(1, max(gha99$sample), 1))
gha99_pred2 <- predict(mod, gha99_pred, se.fit = TRUE)
gha99_pred$pred <- gha99_pred2[,1]
gha99_med <- ddply(gha99,.(sample), summarize,  tot_area_med = median(tot_area))
gha99 <- subset(gha99, iteration <= 100,)

gha99_p <- ggplot() + geom_line(data = gha99, aes(x = sample, y = tot_area, group = as.factor(iteration)), linewidth = 0.3, col = "grey80") +
  geom_hline(yintercept=mod$coefficients[2]*0.95, linetype="dashed", size = 0.3)+
  geom_hline(yintercept=mod$coefficients[2], linetype="solid", size = 0.3)+
  geom_line(data = gha99_pred, aes(x= sample, y = pred), colour = "blue") + theme_bw() +
  geom_point(data = gha99_med, aes(x = sample, y = tot_area_med), pch = 21, size = 1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(t = 0.05, r = 0.01,b = -0.02, l = 0), "cm"),
        axis.text = element_text(size = 7), axis.title.y = element_blank(),axis.title.x = element_blank()) +
  ylab(expression("HR area (10"^6~km^2~")"))+
  scale_y_continuous(limits = c(0, 7.6), breaks = seq(0, 7, 1), expand = c(0, 0))+
  scale_x_continuous(breaks = seq(0, max(gha99$sample), 4))  +
  geom_vline(xintercept=min_sample$prediction[min_sample$dataset == "GHA99"], linetype="dotted")
gha99_p



##### 1J. LMA 2004 #####

lma04 <- read.csv(paste0(resample_direx, "LMA04_50.csv"))
load(paste0(model_direx, "LMA04_50.Rdata"))
lma04_pred <- data.frame(sample = seq(1, max(lma04$sample), 1))
lma04_pred2 <- predict(mod, lma04_pred, se.fit = TRUE)
lma04_pred$pred <- lma04_pred2[,1]
lma04_med <- ddply(lma04,.(sample), summarize,  tot_area_med = median(tot_area))
lma04 <- subset(lma04, iteration <= 100,)

lma04_p <- ggplot() + geom_line(data = lma04, aes(x = sample, y = tot_area, group = as.factor(iteration)), linewidth = 0.3, col = "grey80") +
  geom_hline(yintercept=mod$coefficients[2]*0.95, linetype="dashed", size = 0.3)+
  geom_hline(yintercept=mod$coefficients[2], linetype="solid", size = 0.3)+
  geom_line(data = lma04_pred, aes(x= sample, y = pred), colour = "blue") + theme_bw() +
  geom_point(data = lma04_med, aes(x = sample, y = tot_area_med), pch = 21, size = 1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(t = 0.05, r = 0.01,b = -0.02, l = 0), "cm"),
        axis.text = element_text(size = 7), axis.title.y = element_blank(),axis.title.x = element_blank()) +
  ylab(expression("HR area (10"^6~km^2~")"))+
  scale_y_continuous(limits = c(0, 4.85), breaks = seq(0, 4, 1), expand = c(0, 0))+
  scale_x_continuous(breaks = seq(0, max(lma04$sample), 4)) + 
  geom_vline(xintercept=min_sample$prediction[min_sample$dataset == "LMA04"], linetype="dotted")
lma04_p



##### 1K. LMA 2005 #####

lma05 <- read.csv(paste0(resample_direx, "LMA05_50.csv"))
load(paste0(model_direx, "LMA05_50.Rdata"))
lma05_pred <- data.frame(sample = seq(1, max(lma05$sample), 1))
lma05_pred2 <- predict(mod, lma05_pred, se.fit = TRUE)
lma05_pred$pred <- lma05_pred2[,1]
lma05_med <- ddply(lma05,.(sample), summarize,  tot_area_med = median(tot_area))
lma05 <- subset(lma05, iteration <= 100,)

lma05_p <- ggplot() + geom_line(data = lma05, aes(x = sample, y = tot_area, group = as.factor(iteration)), linewidth = 0.3, col = "grey80") +
  geom_hline(yintercept=mod$coefficients[2]*0.95, linetype="dashed", size = 0.3)+
  geom_hline(yintercept=mod$coefficients[2], linetype="solid", size = 0.3)+
  geom_line(data = lma05_pred, aes(x= sample, y = pred), colour = "blue") + theme_bw() +
  geom_point(data = lma05_med, aes(x = sample, y = tot_area_med), pch = 21, size = 1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(t = 0.05, r = 0.01,b = -0.02, l = 0), "cm"),
        axis.text = element_text(size = 7), axis.title.y = element_blank(),axis.title.x = element_blank()) +
  ylab(expression("HR area (10"^6~km^2~")"))+
  scale_y_continuous(limits = c(0, 4.85), breaks = seq(0, 4, 1), expand = c(0, 0))+
  scale_x_continuous(breaks = seq(0, max(lma05$sample), 2)) + 
  geom_vline(xintercept=min_sample$prediction[min_sample$dataset == "LMA05"], linetype="dotted")
lma05_p



##### 1L. NGP 2000 #####

ngp00 <- read.csv(paste0(resample_direx, "NGP00_50.csv"))
load(paste0(model_direx, "NGP00_50.Rdata"))
ngp00_pred <- data.frame(sample = seq(1, max(ngp00$sample), 1))
ngp00_pred2 <- predict(mod, ngp00_pred, se.fit = TRUE)
ngp00_pred$pred <- ngp00_pred2[,1]
ngp00_med <- ddply(ngp00,.(sample), summarize,  tot_area_med = median(tot_area))
ngp00 <- subset(ngp00, iteration <= 100,)

ngp00_p <- ggplot() + geom_line(data = ngp00, aes(x = sample, y = tot_area, group = as.factor(iteration)), linewidth = 0.3, col = "grey80") +
  geom_hline(yintercept=mod$coefficients[2]*0.95, linetype="dashed", size = 0.3)+
  geom_hline(yintercept=mod$coefficients[2], linetype="solid", size = 0.3)+
  geom_line(data = ngp00_pred, aes(x= sample, y = pred), colour = "blue") + theme_bw() +
  geom_point(data = ngp00_med, aes(x = sample, y = tot_area_med), pch = 21, size = 1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(t = 0.05, r = 0.01,b = -0.02, l = 0), "cm"),
        axis.text = element_text(size = 7), axis.title.y = element_blank(),axis.title.x = element_blank()) +
  ylab(expression("HR area (10"^6~km^2~")"))+
  scale_y_continuous(limits = c(0, 2.3), breaks = seq(0, 2.5, 0.5), expand = c(0, 0))+
  scale_x_continuous(breaks = seq(0, max(ngp00$sample), 4)) + 
  geom_vline(xintercept=min_sample$prediction[min_sample$dataset == "NGP00"], linetype="dotted")
ngp00_p



##### 1M. SGP 2000 #####

sgp00 <- read.csv(paste0(resample_direx, "SGP00_50.csv"))
load(paste0(model_direx, "SGP00_50.Rdata"))
sgp00_pred <- data.frame(sample = seq(1, max(sgp00$sample), 1))
sgp00_pred2 <- predict(mod, sgp00_pred, se.fit = TRUE)
sgp00_pred$pred <- sgp00_pred2[,1]
sgp00_med <- ddply(sgp00,.(sample), summarize,  tot_area_med = median(tot_area))
sgp00 <- subset(sgp00, iteration <= 100,)

sgp00_p <- ggplot() + geom_line(data = sgp00, aes(x = sample, y = tot_area, group = as.factor(iteration)), linewidth = 0.3, col = "grey80") +
  geom_hline(yintercept=mod$coefficients[2]*0.95, linetype="dashed", size = 0.3)+
  geom_hline(yintercept=mod$coefficients[2], linetype="solid", size = 0.3)+
  geom_line(data = sgp00_pred, aes(x= sample, y = pred), colour = "blue") + theme_bw() +
  geom_point(data = sgp00_med, aes(x = sample, y = tot_area_med), pch = 21, size = 1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(t = 0.05, r = 0.01,b = -0.02, l = 0), "cm"),
        axis.text = element_text(size = 7), axis.title.y = element_blank(),axis.title.x = element_blank()) +
  ylab(expression("HR area (10"^6~km^2~")"))+
  scale_y_continuous(limits = c(0, 2.8), breaks = seq(0, 3, 0.5), expand = c(0, 0))+
  scale_x_continuous(breaks = seq(0, max(sgp00$sample), 4)) + 
  geom_vline(xintercept=min_sample$prediction[min_sample$dataset == "SGP00"], linetype="dotted")
sgp00_p


##### 1N. WA 2004 #####

wa04 <- read.csv(paste0(resample_direx, "WA04_50.csv"))
load(paste0(model_direx, "WA04_50.Rdata"))
wa04_pred <- data.frame(sample = seq(1, max(wa04$sample), 1))
wa04_pred2 <- predict(mod, wa04_pred, se.fit = TRUE)
wa04_pred$pred <- wa04_pred2[,1]
wa04_med <- ddply(wa04,.(sample), summarize,  tot_area_med = median(tot_area))
wa04 <- subset(wa04, iteration <= 100,)

wa04_p <- ggplot() + geom_line(data = wa04, aes(x = sample, y = tot_area, group = as.factor(iteration)), linewidth = 0.3, col = "grey80") +
  geom_hline(yintercept=mod$coefficients[2]*0.95, linetype="dashed", size = 0.3)+
  geom_hline(yintercept=mod$coefficients[2], linetype="solid", size = 0.3)+
  geom_line(data = wa04_pred, aes(x= sample, y = pred), colour = "blue") + theme_bw() +
  geom_point(data = wa04_med, aes(x = sample, y = tot_area_med), pch = 21, size = 1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(t = 0.05, r = 0.01,b = -0.02, l = 0), "cm"),
        axis.text = element_text(size = 7), axis.title.y = element_blank(),axis.title.x = element_blank()) +
  ylab(expression("HR area (10"^6~km^2~")"))+
  scale_y_continuous(limits = c(0, 9.9), breaks = seq(0, 8, 2), expand = c(0, 0))+
  scale_x_continuous(breaks = seq(0, max(wa04$sample), 2)) 
  #geom_vline(xintercept=min_sample$prediction[min_sample$dataset == "WA04"], linetype="dotted")
wa04_p


##### 1O. WA 2008 #####

wa08 <- read.csv(paste0(resample_direx, "WA08_50.csv"))
load(paste0(model_direx, "WA08_50.Rdata"))
wa08_pred <- data.frame(sample = seq(1, max(wa08$sample), 1))
wa08_pred2 <- predict(mod, wa08_pred, se.fit = TRUE)
wa08_pred$pred <- wa08_pred2[,1]
wa08_med <- ddply(wa08,.(sample), summarize,  tot_area_med = median(tot_area))
wa08 <- subset(wa08, iteration <= 100,)

wa08_p <- ggplot() + geom_line(data = wa08, aes(x = sample, y = tot_area, group = as.factor(iteration)), linewidth = 0.3, col = "grey80") +
  geom_hline(yintercept=mod$coefficients[2]*0.95, linetype="dashed", size = 0.3)+
  geom_hline(yintercept=mod$coefficients[2], linetype="solid", size = 0.3)+
  geom_line(data = wa08_pred, aes(x= sample, y = pred), colour = "blue") + theme_bw() +
  geom_point(data = wa08_med, aes(x = sample, y = tot_area_med), pch = 21, size = 1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(t = 0.05, r = 0.01,b = -0.02, l = 0), "cm"),
        axis.text = element_text(size = 7), axis.title.y = element_blank(),axis.title.x = element_blank()) +
  ylab(expression("HR area (10"^6~km^2~")"))+
  scale_y_continuous(limits = c(0, 9.9), breaks = seq(0, 8, 2), expand = c(0, 0))+
  scale_x_continuous(breaks = seq(0, max(wa08$sample), 4)) + 
  geom_vline(xintercept=min_sample$prediction[min_sample$dataset == "WA08"], linetype="dotted")
wa08_p


##### 1P. WA 2009 #####

wa09 <- read.csv(paste0(resample_direx, "WA09_50.csv"))
load(paste0(model_direx, "WA09_50.Rdata"))
wa09_pred <- data.frame(sample = seq(1, max(wa09$sample), 1))
wa09_pred2 <- predict(mod, wa09_pred, se.fit = TRUE)
wa09_pred$pred <- wa09_pred2[,1]
wa09_med <- ddply(wa09,.(sample), summarize,  tot_area_med = median(tot_area))
wa09 <- subset(wa09, iteration <= 100,)

wa09_p <- ggplot() + geom_line(data = wa09, aes(x = sample, y = tot_area, group = as.factor(iteration)), linewidth = 0.3, col = "grey80") +
  geom_hline(yintercept=mod$coefficients[2]*0.95, linetype="dashed", size = 0.3)+
  geom_hline(yintercept=mod$coefficients[2], linetype="solid", size = 0.3)+
  geom_line(data = wa09_pred, aes(x= sample, y = pred), colour = "blue") + theme_bw() +
  geom_point(data = wa09_med, aes(x = sample, y = tot_area_med), pch = 21, size = 1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(t = 0.05, r = 0.01,b = -0.02, l = 0), "cm"),
        axis.text = element_text(size = 7), axis.title.y = element_blank(),axis.title.x = element_blank()) +
  ylab(expression("HR area (10"^6~km^2~")"))+
  scale_y_continuous(limits = c(0, 9.9), breaks = seq(0, 8, 2), expand = c(0, 0))+
  scale_x_continuous(breaks = seq(0, max(wa09$sample), 4)) + 
  geom_vline(xintercept=min_sample$prediction[min_sample$dataset == "WA09"], linetype="dotted")
wa09_p

##### 1Q. WA ALL #####

waall <- read.csv(paste0(resample_direx, "WAALL_50.csv"))
load(paste0(model_direx, "WAALL_50.Rdata"))
waall_pred <- data.frame(sample = seq(1, max(waall$sample), 1))
waall_pred2 <- predict(mod, waall_pred, se.fit = TRUE)
waall_pred$pred <- waall_pred2[,1]
waall_med <- ddply(waall,.(sample), summarize,  tot_area_med = median(tot_area))
waall <- subset(waall, iteration <= 100,)

waall_p <- ggplot() + geom_line(data = waall, aes(x = sample, y = tot_area, group = as.factor(iteration)), linewidth = 0.3, col = "grey80") +
  geom_hline(yintercept=mod$coefficients[2]*0.95, linetype="dashed", size = 0.3)+
  geom_hline(yintercept=mod$coefficients[2], linetype="solid", size = 0.3)+
  geom_line(data = waall_pred, aes(x= sample, y = pred), colour = "blue") + theme_bw() +
  geom_point(data = waall_med, aes(x = sample, y = tot_area_med), pch = 21, size = 1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(t = 0.05, r = 0.01,b = -0.02, l = 0), "cm"),
        axis.text = element_text(size = 7), axis.title.y = element_blank(),axis.title.x = element_blank()) +
  ylab(expression("HR area (10"^6~km^2~")"))+
  scale_y_continuous(limits = c(0, 9.9), breaks = seq(0, 8, 2), expand = c(0, 0))+
  scale_x_continuous(breaks = seq(0, max(waall$sample), 10)) + 
  geom_vline(xintercept=min_sample$prediction[min_sample$dataset == "WAALL"], linetype="dotted")
waall_p


##### 1R. WCP 2003 #####

wcp03 <- read.csv(paste0(resample_direx, "WCP03_50.csv"))
load(paste0(model_direx, "WCP03_50.Rdata"))
wcp03_pred <- data.frame(sample = seq(1, max(wcp03$sample), 1))
wcp03_pred2 <- predict(mod, wcp03_pred, se.fit = TRUE)
wcp03_pred$pred <- wcp03_pred2[,1]
wcp03_med <- ddply(wcp03,.(sample), summarize,  tot_area_med = median(tot_area))
wcp03 <- subset(wcp03, iteration <= 100,)

wcp03_p <- ggplot() + geom_line(data = wcp03, aes(x = sample, y = tot_area, group = as.factor(iteration)), linewidth = 0.3, col = "grey80") +
  geom_hline(yintercept=mod$coefficients[2]*0.95, linetype="dashed", size = 0.3)+
  geom_hline(yintercept=mod$coefficients[2], linetype="solid", size = 0.3)+
  geom_line(data = wcp03_pred, aes(x= sample, y = pred), colour = "blue") + theme_bw() +
  geom_point(data = wcp03_med, aes(x = sample, y = tot_area_med), pch = 21, size = 1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(t = 0.05, r = 0.01,b = -0.02, l = 0), "cm"),
        axis.text = element_text(size = 7), axis.title.y = element_blank(),axis.title.x = element_blank()) +
  ylab(expression("HR area (10"^6~km^2~")"))+
  scale_y_continuous(limits = c(0, 2.35), breaks = seq(0, 2, 0.5), expand = c(0, 0))+
  scale_x_continuous(breaks = seq(0, max(wcp03$sample), 1)) 
  #geom_vline(xintercept=min_sample$prediction[min_sample$dataset == "WCP03"], linetype="dotted")
wcp03_p


##### 1S. WCP 2013 #####

wcp13 <- read.csv(paste0(resample_direx, "WCP13_50.csv"))
load(paste0(model_direx, "WCP13_50.Rdata"))
wcp13_pred <- data.frame(sample = seq(1, max(wcp13$sample), 1))
wcp13_pred2 <- predict(mod, wcp13_pred, se.fit = TRUE)
wcp13_pred$pred <- wcp13_pred2[,1]
wcp13_med <- ddply(wcp13,.(sample), summarize,  tot_area_med = median(tot_area))
wcp13 <- subset(wcp13, iteration <= 100,)

wcp13_p <- ggplot() + geom_line(data = wcp13, aes(x = sample, y = tot_area, group = as.factor(iteration)), linewidth = 0.3, col = "grey80") +
  geom_hline(yintercept=mod$coefficients[2]*0.95, linetype="dashed", size = 0.3)+
  geom_hline(yintercept=mod$coefficients[2], linetype="solid", size = 0.3)+
  geom_line(data = wcp13_pred, aes(x= sample, y = pred), colour = "blue") + theme_bw() +
  geom_point(data = wcp13_med, aes(x = sample, y = tot_area_med), pch = 21, size = 1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(t = 0.05, r = 0.01,b = -0.02, l = 0), "cm"),
        axis.text = element_text(size = 7), axis.title.y = element_blank(),axis.title.x = element_blank()) +
  ylab(expression("HR area (10"^6~km^2~")"))+
  scale_y_continuous(limits = c(0, 2.35), breaks = seq(0, 2, 0.5), expand = c(0, 0))+
  scale_x_continuous(breaks = seq(0, max(wcp13$sample), 1)) 
  #geom_vline(xintercept=min_sample$prediction[min_sample$dataset == "WCP13"], linetype="dotted")
wcp13_p


##### 1T. WCP BOTH #####

wcpboth <- read.csv(paste0(resample_direx, "WCPBOTH_50.csv"))
load(paste0(model_direx, "WCPBOTH_50.Rdata"))
wcpboth_pred <- data.frame(sample = seq(1, max(wcpboth$sample), 1))
wcpboth_pred2 <- predict(mod, wcpboth_pred, se.fit = TRUE)
wcpboth_pred$pred <- wcpboth_pred2[,1]
wcpboth_med <- ddply(wcpboth,.(sample), summarize,  tot_area_med = median(tot_area))
wcpboth <- subset(wcpboth, iteration <= 100,)

wcpboth_p <- ggplot() + geom_line(data = wcpboth, aes(x = sample, y = tot_area, group = as.factor(iteration)), linewidth = 0.3, col = "grey80") +
  geom_hline(yintercept=mod$coefficients[2]*0.95, linetype="dashed", size = 0.3)+
  geom_hline(yintercept=mod$coefficients[2], linetype="solid", size = 0.3)+
  geom_line(data = wcpboth_pred, aes(x= sample, y = pred), colour = "blue") + theme_bw() +
  geom_point(data = wcpboth_med, aes(x = sample, y = tot_area_med), pch = 21, size = 1) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(t = 0.05, r = 0.01,b = -0.02, l = 0), "cm"),
        axis.text = element_text(size = 7), axis.title.y = element_blank(),axis.title.x = element_blank()) +
  ylab(expression("HR area (10"^6~km^2~")"))+
  scale_y_continuous(limits = c(0, 2.35), breaks = seq(0, 2, 0.5), expand = c(0, 0))+
  scale_x_continuous(breaks = seq(0, max(wcpboth$sample), 2))  
  #geom_vline(xintercept=min_sample$prediction[min_sample$dataset == "WCPBOTH"], linetype="dotted")
wcpboth_p


### 2. JOIN INTO MULTIPANEL PLOT AND PASTE OUT ####


p_all <- plot_grid(ap10_p, ap11_p, apboth_p, bba02_p, 
                  bba03_p, bp11_p, bs02_p, bs12_p, 
                  bsboth_p, gha99_p, lma04_p, lma05_p,
                  ngp00_p, sgp00_p, wa04_p, wa08_p, 
                  wa09_p, waall_p, wcp03_p, wcp13_p, 
                  wcpboth_p,
                  nrow = 6, ncol = 4, align= "v")
print(p_all)

out.name <- paste0(plot_direx, "2_asymptote_plots.png")
ggsave(out.name, dpi = 600, width = 3200, height = 3762, units = "px")
dev.off()

