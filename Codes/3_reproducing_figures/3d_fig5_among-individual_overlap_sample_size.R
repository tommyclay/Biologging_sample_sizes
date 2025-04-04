library(ggplot2)

plot_direx <- "./Plots/"

#### 1. LOAD IN AND PREPARE DATA ####

# much of this section is copied from 2f.
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

# setting up function and predicting onto dummy dataset
# editing function from function codes to not round up numbers to predict to - see functions.R code
prediction_from_ba_woround <- function(ba) {
  x <- ba$ba
  # Compute predictions
  y <- -3.9926/log(1-x)
  # Compute the gradient manually: derivative w.r.t k
  grad_matrix <- matrix(1 / log(1 - x), ncol = 1)  # Ensure it's a matrix
  # Compute standard error of predictions from vcov of model
  se_fit <- sqrt(grad_matrix^2 * 0.05120739)
  # Compute 95% confidence intervals
  alpha <- 0.05
  z_score <- qnorm(1 - alpha / 2)  # 1.96 for 95% CI
  upper <- y + z_score * se_fit
  lower <- y - z_score * se_fit
  # Output dataframe as whole numbers
  out_df <- data.frame(ba = x, prediction = y, lower = lower, upper = upper)
  return(out_df)
}
fake_df <- data.frame(ba = seq(0.05, 1, 0.005))
pred_df <- prediction_from_ba_woround(fake_df)

# load in oppel 2018 values
oppel <- read.csv("./Data_inputs/oppel_et_al_2018_ba_scores.csv")
oppel_pred <- prediction_from_ba(oppel) 
oppel$prediction <- oppel_pred$prediction


#### 2. PLOTTING #####


##### PANEL A. RELATIONSHIP BETWEEN BA AND MINIMUM SAMPLE SIZES  #####

all_preds$Species <- all_preds$species
all_preds$Device <- all_preds$device

p1 <- ggplot() + geom_point(data = all_preds, aes(x = ba, y = prediction, colour = Species, shape = Device)) + 
  geom_ribbon(data = pred_df, aes(x = ba, ymin = lower, ymax = upper), alpha = 0.3, fill = "grey50") +
  geom_line(data = pred_df, aes(x = ba, y = prediction)) + theme_bw() + xlab("BA overlap") + ylab("Minimum sample size") +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1))+
  scale_y_continuous(limits = c(1, 87), breaks = seq(0, 80, 10))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        axis.title = element_text(size = 13), axis.text = element_text(size = 11),
        legend.title = element_text(size = 11), legend.text = element_text(size = 9),
        legend.margin = margin(0,0,0,0, unit="cm"),
        legend.key.size = unit(0.45, "cm")) +
  scale_colour_manual(values = c("AP" = "#D55E00", "BBA" = "#E69F00", "BP" = "#F0E442",
                                 "BS" = "#009E73", "GHA" = "#56B4E9", "LMA" = "#0072B2",
                                 "NGP" = "#CC79A7", "SGP" = "#882255",
                                 "WA" = "#000000", "WCP" = "grey60")) +
  geom_hline(yintercept=1, linetype="dashed")
p1

out.name <- paste0(plot_direx, "5a_ba_sample_relationship.png")
ggsave(out.name, dpi = 600, width = 2700, height = 2200, units = "px")
dev.off()


##### PANEL B. HISTOGRAMS OF BA #####

oppel$family <- factor(oppel$family, levels = c("Phalacrocoracidae", "Laridae",  "Alcidae", "Sulidae", "Spheniscidae", "Hydrobatidae",
                                          "Phaethontidae", "Procellariidae", "Fregatidae", "Diomedeidae"))

p2 <- ggplot(oppel, aes(family, ba)) + geom_boxplot() + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  axis.text.x = element_text(angle = 90, size = 14, vjust = 0.6),
  axis.text.y = element_text(size = 14)) +
  scale_y_continuous(limits = c(0,1), breaks= seq(0,1, 0.2))
p2

out.name <- paste0(plot_direx, "5b_oppel2018_ba.png")
ggsave(out.name, dpi = 600, width = 2000, height = 2500, units = "px")
dev.off()


##### PANEL C. HISTOGRAMS OF PREDICTED SAMPLE SIZE #####

p3 <- ggplot(oppel, aes(family, prediction)) + geom_boxplot() + theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 90, size = 14, vjust = 0.6),
        axis.text.y = element_text(size = 14)) +
  scale_y_continuous(limits = c(0,65), breaks= seq(0,60, 10))
p3

out.name <- paste0(plot_direx, "5c_oppel2018_predictions.png")
ggsave(out.name, dpi = 600, width = 2000, height = 2500, units = "px")
dev.off()
