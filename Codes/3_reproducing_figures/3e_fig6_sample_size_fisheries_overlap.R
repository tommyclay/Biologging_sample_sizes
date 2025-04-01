library(ggplot2)
library(cowplot)
library(plyr)
library(scales)

# setting up directories
in_direx <- "./Data_outputs/resampled_fisheries_overlap_outputs/"
plot_direx <- "./Plots/"


#### 1. CREATING A SEPARATE PANEL FOR EACH SPECIES ####

##### 1A. BBA 2002 #####

bba <- read.csv(paste0(in_direx, "bba02_50.csv"))
bba$o_change_p[is.na(bba$o_change_p) == FALSE & bba$o_change_p<0] <- bba$o_change_p[is.na(bba$o_change_p) == FALSE & bba$o_change_p<0]  *-1
bba_med <- ddply(bba,.(sample), summarize, overlap_med = median(overlap),  o_change_med = median(o_change_p))
bba <- subset(bba, iteration <= 100,)
bba_5 <- subset(bba_med,  o_change_med <= 5,)
bba_5$sample[1] # 6

# 1. plotting overlap
p_bba <- ggplot(bba, aes(sample, overlap)) + geom_line(aes(group = as.factor(iteration)),linewidth = 0.3, col = "grey80")+
  stat_summary(fun=mean, colour="blue", geom="line", size = 0.5) + theme_bw() +
  geom_point(data = bba_med, aes(x = sample, y = overlap_med), pch = 21, size= 0.9) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(t = 0.05, r = 0.05,b = -0.02, l = 0), "cm"),
        axis.text = element_text(size = 7), axis.title.y = element_blank(),axis.title.x = element_blank()) +
  scale_y_continuous(limits = c(0, 1.7), expand = c(0, 0))+scale_x_continuous(breaks = seq(0, 24, 4)) 
p_bba

# 2. plotting change in overlap
p_bba_c <- ggplot(bba, aes(sample, o_change_p)) + geom_line(aes(group = as.factor(iteration)),linewidth = 0.3, col = "grey80")+
  stat_summary(fun=mean, colour="blue", geom="line", size = 0.5) + theme_bw() +
  geom_point(data = bba_med, aes(x = sample, y = o_change_med), pch = 21, size= 0.9) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(t = 0.05, r = 0.05,b = -0.02, l = 0), "cm"),
        axis.text = element_text(size = 7),axis.title.y = element_blank(),axis.title.x = element_blank()) +
  scale_y_continuous(limits = c(0.1, 87050), trans = "log", breaks = c(0.1, 1, 10, 100, 1000, 10000),expand = c(0, 0),
                     label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", "%*%10^", scientific_format()(x)))))})+
  scale_x_continuous(limits = c(0, 25), breaks = seq(0, 24, 4)) +
  geom_hline(yintercept=5, linetype="dashed", size = 0.3)+ geom_vline(xintercept=bba_5$sample[1], linetype="dotted")
p_bba_c 




##### 1B. GHA 1999 #####

gha <- read.csv(paste0(in_direx, "gha99_50.csv"))
gha$o_change_p[is.na(gha$o_change_p) == FALSE & gha$o_change_p<0] <- gha$o_change_p[is.na(gha$o_change_p) == FALSE & gha$o_change_p<0]  *-1
gha_med <- ddply(gha,.(sample), summarize, overlap_med = median(overlap), o_change_med = median(o_change_p))
gha <- subset(gha, iteration <= 100,)
gha_5 <- subset(gha_med,  o_change_med <= 5,)
gha_5$sample[1] # 13

# 1. plotting overlap
p_gha <- ggplot(gha, aes(sample, overlap)) + geom_line(aes(group = as.factor(iteration)),linewidth = 0.3, col = "grey80")+
  stat_summary(fun=mean, colour="blue", geom="line", size = 0.5) + theme_bw() +
  geom_point(data = gha_med, aes(x = sample, y = overlap_med), pch = 21, size= 0.9) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(t = 0.05, r = 0.05,b = -0.02, l = 0), "cm"),
        axis.text = element_text(size = 7),axis.title.y = element_blank(),axis.title.x = element_blank()) +
  scale_y_continuous(limits = c(0, 1.7), expand = c(0, 0))+
  scale_x_continuous(breaks = seq(0, 24, 4)) 
p_gha

# 2. plotting change in overlap
p_gha_c <- ggplot(gha, aes(sample, o_change_p)) + geom_line(aes(group = as.factor(iteration)),linewidth = 0.3, col = "grey80")+
  stat_summary(fun=mean, colour="blue", geom="line", size = 0.5) + theme_bw() +
  geom_point(data = gha_med, aes(x = sample, y = o_change_med), pch = 21, size= 0.9) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(t = 0.05, r = 0.05,b = -0.02, l = 0), "cm"),
        axis.text = element_text(size = 7), axis.title.y = element_blank(),axis.title.x = element_blank()) +
  scale_y_continuous(limits = c(0.1, 87050), trans = "log", breaks = c(0.1, 1, 10, 100, 1000, 10000),expand = c(0, 0),
                     label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", "%*%10^", scientific_format()(x)))))})+
  scale_x_continuous(limits = c(0, 22), breaks = seq(0, 24, 4)) +
  geom_hline(yintercept=5, linetype="dashed", size = 0.3)+ geom_vline(xintercept=gha_5$sample[1], linetype="dotted")
p_gha_c 


##### 1C. WA ALL YEARS #####

wa <- read.csv(paste0(in_direx, "waall_50.csv"))
# making negative %s positive
wa$o_change_p[is.na(wa$o_change_p) == FALSE & wa$o_change_p<0] <- wa$o_change_p[is.na(wa$o_change_p) == FALSE & wa$o_change_p<0]  *-1
# getting median values
wa_med <- ddply(wa,.(sample), summarize, overlap_med = median(overlap), o_change_med = median(o_change_p))
# picking out the first 100 iterations to plot
wa <- subset(wa, iteration <= 100,)
# getting first median value <= 5%
wa_5 <- subset(wa_med,  o_change_med <= 5,)
wa_5$sample[1] # 11

# 1. plotting raw overlap
p_wa <- ggplot(wa, aes(sample, overlap)) + geom_line(aes(group = as.factor(iteration)), linewidth = 0.3, col = "grey80")+
  stat_summary(fun=mean, colour="blue", geom="line", size = 0.5) + theme_bw() +
  geom_point(data = wa_med, aes(x = sample, y = overlap_med), pch = 21, size= 0.9) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(t = 0.05, r = 0.05,b = -0.02, l = 0), "cm"),
        axis.text = element_text(size = 7), axis.title.y = element_blank(),axis.title.x = element_blank()) +
  scale_y_continuous(limits = c(0, 1.7),expand = c(0, 0))+ scale_x_continuous(breaks = seq(0, 80, 10)) 
p_wa

# 2. plotting change in overlap
p_wa_c <- ggplot(wa, aes(sample, o_change_p)) + geom_line(aes(group = as.factor(iteration)),linewidth = 0.3, col = "grey80")+
  stat_summary(fun=mean, colour="blue", geom="line", size = 0.5) + theme_bw() +
  geom_point(data = wa_med, aes(x = sample, y = o_change_med), pch = 21, size= 0.9) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(t = 0.05, r = 0.05,b = -0.02, l = 0), "cm"),
        axis.text = element_text(size = 7), axis.title.y = element_blank(), axis.title.x = element_blank()) +
  scale_y_continuous(limits = c(0.1, 87050), trans = "log", breaks = c(0.1, 1, 10, 100, 1000, 1000, 10000), expand = c(0, 0),
                     label = function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", "%*%10^", scientific_format()(x)))))})+
  scale_x_continuous(limits = c(0, 80), breaks = seq(0, 80, 10)) +
  geom_hline(yintercept=5, linetype="dashed", size = 0.3)+ geom_vline(xintercept=wa_5$sample[1], linetype="dotted")
p_wa_c


##### 1D. WCP BOTH #####

wcp <- read.csv(paste0(in_direx, "wcpboth_50.csv"))
wcp$o_change_p[is.na(wcp$o_change_p) == FALSE & wcp$o_change_p<0] <- wcp$o_change_p[is.na(wcp$o_change_p) == FALSE & wcp$o_change_p<0]  *-1
wcp_med <- ddply(wcp,.(sample), summarize, overlap_med = median(overlap), o_change_med = median(o_change_p))
wcp <- subset(wcp, iteration <= 100,)
wcp_5 <- subset(wcp_med,  o_change_med <= 5,)
wcp_5 # 5% threshold not reached

# 1. plotting overlap
p_wcp <- ggplot(wcp, aes(sample, overlap)) + geom_line(aes(group = as.factor(iteration)), linewidth = 0.3, col = "grey80")+
  stat_summary(fun=mean, colour="blue", geom="line", size = 0.5) + theme_bw() +
  geom_point(data = wcp_med, aes(x = sample, y = overlap_med), pch = 21, size= 0.9) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(t = 0.05, r = 0.05,b = -0.02, l = 0), "cm"),
        axis.text = element_text(size = 7), axis.title.y = element_blank(),axis.title.x = element_blank()) +
  scale_y_continuous(limits = c(0, 1.7), expand = c(0, 0)) + scale_x_continuous(limits = c(0, 20), breaks = seq(0, 20, 4)) 
p_wcp

# 2. plotting change in overlap
p_wcp_c <- ggplot(wcp, aes(sample, o_change_p)) + geom_line(aes(group = as.factor(iteration)), linewidth = 0.3, col = "grey80")+
  stat_summary(fun=mean, colour="blue", geom="line", size = 0.5) + theme_bw() +
  geom_point(data = wcp_med, aes(x = sample, y = o_change_med), pch = 21, size= 0.9) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.margin = unit(c(t = 0.05, r = 0.05,b = -0.02, l = 0), "cm"),
        axis.text = element_text(size = 7), axis.title.y = element_blank(), axis.title.x = element_blank()) +
  scale_y_continuous(limits = c(0.1, 87050), trans = "log", breaks = c(0.1, 1, 10, 100, 1000, 10000), expand = c(0, 0),
                     label= function(x) {ifelse(x==0, "0", parse(text=gsub("[+]", "", gsub("e", "%*%10^", scientific_format()(x)))))})+
  scale_x_continuous(limits = c(0, 20), breaks = seq(0, 24, 4)) +
  geom_hline(yintercept=5, linetype="dashed", size = 0.3) #+geom_vline(xintercept=wcp_5$sample[1], linetype="dotted")
p_wcp_c 



#### 2. MERGING TOGETHER AND PASTING OUT ####

plot_all <- plot_grid(p_bba, p_gha, p_wa, p_wcp,
                   p_bba_c, p_gha_c, p_wa_c, p_wcp_c,
                   nrow = 2, ncol = 4, align= "v")
print(plot_all)

out.name <- paste0(plot_direx, "6_fisheries_overlap.png")
ggsave(out.name, dpi = 600, width = 3800, height = 1500, units = "px")
dev.off()
