
library(ggplot2)
library(ggrepel)

#### 1. LOAD IN DATA ####

# this first section is copied from 2d (and edited slightly)
pred_red <- read.csv("./Data_outputs/sensitivity_analysis_predictions.csv")
pred_red$type <-"subset"
pred_og <- read.csv("./Data_outputs/area_occupied_50_predictions.csv")
pred_og <- subset(pred_og, model == "MM3",)
pred_og <- subset(pred_og, pop_size == 100,)
pred_og$type <-"final"
# pasting in original sample sizes
pred_og$n_ind <- c(9, 6, 14, 25, 24, 9, 6, 19, 25, 22, 25, 11, 25, 29, 15, 37, 33, 78, 10, 9, 19)
# combine into one dataframe 
preds <- data.frame(prediction = c(pred_og$prediction, pred_red$prediction), dataset = c(pred_og$dataset, pred_red$dataset), 
                    sample_og = c(pred_og$n_ind, pred_red$n_ind), type = c(pred_og$type, pred_red$type))
preds$Species <- gsub("[0-9]", "", preds$dataset)
preds$yr <-gsub("[a-zA-Z]", "", preds$dataset)
preds$sp_yr <- as.factor(as.character(paste(preds$Species, preds$yr, sep = "_")))
# remove combined years of data
preds <- subset(preds, dataset != "APBOTH",)
preds <- subset(preds, dataset != "WCPBOTH",)
preds <- subset(preds, dataset != "WAALL",)
preds <- subset(preds, dataset != "BSBOTH",)

preds$Species <- as.factor(as.character(preds$Species))
preds$yr <- as.factor(as.character(preds$yr))


#### 2. PLOTTING ####

# colourblind friendly palette - colours chosen from - https://clauswilke.com/dataviz/color-pitfalls.html
cols <- c("#E69F00", # orange,
          "#56B4E9", # sky blue
          "#009E73", # bluish green
          "#F0E442", # yellow
          "#0072B2", # blue
          "#D55E00", # vermillon
          "#CC79A7", # reddish purple
          "#000000", # black
          "grey80", # light grey
          "grey30") #  dark grey)
# specify labels
sp_lab <- subset(preds, type == "final",)
sp_lab$sp_yr <- gsub("_", "", sp_lab$sp_yr)

p1 <- ggplot(data = preds, aes(x = sample_og, y = prediction, colour = Species, group = sp_yr)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dotted")+
  geom_point(aes(shape = type), size = 1.7) + geom_line(linewidth = 0.7) +
  scale_colour_manual(values = c("AP" = "#D55E00", "BBA" = "#E69F00", "BP" = "#F0E442",
                                 "BS" = "#009E73", "GHA" = "#56B4E9", "LMA" = "#0072B2",
                                 "NGP" = "#CC79A7", "SGP" = "#882255",
                                 "WA" = "#000000", "WCP" = "grey60")) +
  scale_shape_manual(values =c("final" = 17, "subset" = 21)) +
  theme_bw()+guides(shape = "none") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.title = element_text(size = 13), axis.text = element_text(size = 11),
        legend.title = element_text(size = 13), legend.text = element_text(size = 10))+
  xlab("Original sample size") + ylab("Minimum sample size") +
  scale_x_continuous(limits =c(5,37), breaks = seq(5, 40, 5))+
  scale_y_continuous(limits=c(0, 40), breaks = seq(0, 40, 5)) +
  geom_text_repel(data = sp_lab, aes(x = sample_og, y = prediction, label = sp_yr), size = 2.2, nudge_x = 0.5, nudge_y = 0.5)
p1

plot_direx <- "./Plots/"
out.name <- paste0(plot_direx, "3_og_sample_predictions", ".png")
ggsave(out.name, dpi = 600, width = 3075, height = 2140, units = "px")
dev.off()
