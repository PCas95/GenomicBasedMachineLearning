#### heatmappe.R ####
#### R script to produce heatmaps and correlation plots from metrics of machine learning performance for the research article ####
#### "Source attribution of Listeria monocytogenes through a versatile supervised machine learning workflow based on genomic data" ####
#### performed with R version 4.2.3 (2023-03-15) and RStudio 2022.12.0+353 on Ubuntu 20.04 5 LTS "Focal Fossa" ####


# clean environment
rm(list=ls())

# clean graphical device
graphics.off()

# set working directory
setwd("/home/IZSNT/p.castelli/Documenti/GitHub_repos/heatmappeR/demo")

# install regular packages
#install.packages("ggplot2")
#install.packages("plyr")
#install.packages("ggpmisc")
#install.packages("reshape2")
#install.packages("lubridate")

# call library
library(ggplot2)
library(plyr)
library(ggpmisc)
library(reshape2)
library(lubridate)


# read dataframe
data_raw = read.table("fake_Metrics.tsv", dec = ".", header=TRUE, sep = "\t", quote = "")
## check dimension
dim(data_raw)
## check first 20 lines
head(data_raw, 20)

# check nature of variables (integer or factor)
str(data_raw)
dim(data_raw)

# transform scientific numbers into decimal numbers with R
## format numbers in columns
new_col05 <- as.numeric(format(data_raw$average_accuracy_training, digits=3))
new_col06 <- format(data_raw$CI95_accuracy_training, digits=3)
new_col07 <- as.numeric(format(data_raw$average_accuracy_testing, digits=3))
new_col08 <- format(data_raw$CI95_accuracy_testing, digits=3)
new_col09 <- as.numeric(format(data_raw$F.score, digits=3))
new_col10 <- as.numeric(format(data_raw$Cohen_kappa_training, digits=3))
new_col11 <- as.numeric(format(data_raw$Cohen_kappa_testing, digits=3))
new_col12 <- as.numeric(format(data_raw$ROC_AUC, digits=3))
new_col13 <- as.numeric(format(data_raw$PR_AUC, digits=3))
new_col14 <- as.numeric(format(data_raw$PRG_AUC, digits=3))
new_col15 <- as.numeric(format(data_raw$seconds, digits=5))
## create new dataframe with desired value format
data_long = data.frame(data_raw$mutation,data_raw$preprocessing,data_raw$splitting,data_raw$model,new_col05,new_col06,new_col07,new_col08,new_col09,new_col10,new_col11,new_col12,new_col13,new_col14,new_col15)
## change column names of new dataframe back to original names
colnames(data_long)[1] <- "mutation"
colnames(data_long)[2] <- "preprocessing"
colnames(data_long)[3] <- "splitting"
colnames(data_long)[4] <- "model"
colnames(data_long)[5] <- "average_accuracy_training"
colnames(data_long)[6] <- "CI95_accuracy_training"
colnames(data_long)[7] <- "average_accuracy_testing"
colnames(data_long)[8] <- "CI95_accuracy_testing"
colnames(data_long)[9] <- "F1_score"
colnames(data_long)[10] <- "Cohen_kappa_training"
colnames(data_long)[11] <- "Cohen_kappa_testing"
colnames(data_long)[12] <- "ROC_AUC"
colnames(data_long)[13] <- "PR_AUC"
colnames(data_long)[14] <- "PRG_AUC"
colnames(data_long)[15] <- "seconds"

# check dimension
dim(data_long)

# check 10 first lines
head(data_long, 10)
tail(data_long, 10)

# check nature of variables (integer or factor)
str(data_long)

# express accuracy and AUC in percentage
data_long$average_accuracy_training <- data_long$average_accuracy_training*100
data_long$average_accuracy_testing <- data_long$average_accuracy_testing*100
data_long$ROC_AUC <- data_long$ROC_AUC*100
data_long$PR_AUC <- data_long$PR_AUC*100
data_long$PRG_AUC <- data_long$PRG_AUC*100

# rename levels of a factor
data_long$mutation <- as.factor(data_long$mutation)
levels(data_long$mutation)
data_long$mutation <- revalue(data_long$mutation, 
                              c(
                                "mutation_A"="genes",
                                "mutation_B"="SNPs"
                              ))
levels(data_long$mutation)

# reorganize levels of variables
levels(data_long$mutation)
data_long$mutation <- factor(data_long$mutation, levels=c("genes", "SNPs"))
levels(data_long$mutation)

# reorganize levels of variables
data_long$preprocessing <- as.factor(data_long$preprocessing)
levels(data_long$preprocessing)
data_long$preprocessing <- factor(data_long$preprocessing, levels=c("no", "yes"))
levels(data_long$preprocessing)

# reorganize levels of variables
data_long$model <- as.factor(data_long$model)
levels(data_long$model)
data_long$model <- factor(data_long$model, levels=c("mod_E", "mod_D", "mod_C", "mod_B", "mod_A"))
levels(data_long$model)

# create a different vector for average_accuracy_testing heatmap labels
newCI <- c("lower_CI95_accuracy_testing", "upper_CI95_accuracy_testing")
dfCIte <- colsplit(data_long$CI95_accuracy_testing, "-", newCI)
dfCIte$lower_CI95_accuracy_testing <- dfCIte$lower_CI95_accuracy_testing*100
dfCIte$upper_CI95_accuracy_testing <- dfCIte$upper_CI95_accuracy_testing*100
dfCIte$lower_CI95_accuracy_testing <- format(round(dfCIte$lower_CI95_accuracy_testing, 1), nsmall = 1)
dfCIte$upper_CI95_accuracy_testing <- format(round(dfCIte$upper_CI95_accuracy_testing, 1), nsmall = 1)
CI95_accuracy_testing <- paste(dfCIte$lower_CI95_accuracy_testing, dfCIte$upper_CI95_accuracy_testing, sep = "-")
average_accuracy_testing <- format(round(data_long$average_accuracy_testing, 1), nsmall = 1)
labels01 <- paste(average_accuracy_testing, CI95_accuracy_testing, sep = "\n")

# create a different vector for Cohen Kappa heatmap labels
labels02 <- as.character(sprintf("%.3f", data_long$Cohen_kappa_testing))

# create a different vector for F1 score heatmap labels
labels03 <- as.character(sprintf("%.3f", data_long$F1_score))

# create a different vector for AUC heatmap labels 
labelsROC <- sprintf("%.1f", data_long$ROC_AUC)
labelsPR <- sprintf("%.1f", data_long$PR_AUC)
labelsPRG <- sprintf("%.1f", data_long$PRG_AUC)

# create a different vector for execution time heatmap labels
seconds_int <- as.integer(data_long$seconds)
seconds_dur <- as.duration(seconds_int)
exec_t01 <- gsub(".*~", "", seconds_dur, perl=TRUE)
exec_t02 <- gsub("minutes\\)$", "min", exec_t01, perl=TRUE)
exec_t03 <- gsub("hours\\)$", "h", exec_t02, perl=TRUE)
exec_t04 <- gsub("days\\)$", "d", exec_t03, perl=TRUE)
exec_t05 <- gsub("s$", " sec", exec_t04, perl=TRUE)


# plot heatmaps

## for average_accuracy_testing
p = ggplot(data = data_long, aes(x = splitting, y = model)) + 
  theme_light(base_size = 16) + 
  geom_tile(aes(fill = average_accuracy_testing)) + 
  geom_text(size=2.8, aes(label = labels01)) + 
  scale_fill_gradient(
    low = "#ff0000", 
    high = "#008000", 
    name = "average accuracy \nfrom the testing dataset \n(%)",
    limits = c(round(min(data_long$average_accuracy_testing), digits=2)-0.1, round(max(data_long$average_accuracy_testing), digits=2)+0.1), 
    breaks = c(round(min(data_long$average_accuracy_testing), digits=1), 
               round((max(data_long$average_accuracy_testing)-min(data_long$average_accuracy_testing))/2+min(data_long$average_accuracy_testing), digits=1), 
               round(max(data_long$average_accuracy_testing), digits=1))
  ) + 
  theme(
    legend.position = "top", 
    legend.title = element_text(colour="#000000", size=12, face="bold"), 
    legend.text = element_text(size=12), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.text.x = element_text (color = "#000000", size = 12, angle = 0),
    axis.text.y = element_text (color = "#000000", size = 12, angle = 0), 
    strip.text.x = element_text(size = 10, face = "bold", color = "#000000"),
    strip.text.y = element_text(size = 10, face = "bold", color = "#000000"),
    strip.background = element_rect(colour="#000000", fill="#A9A9A9"), 
    panel.border = element_rect(colour="#FFFFFF"),
    axis.ticks = element_line(colour="#FFFFFF"), 
    plot.title = element_text(size = 15, face = "bold")
  ) +
  scale_x_continuous(
    name = "splitting proportion of the training dataset (%)", 
    limits = c(25, 75), breaks = c(30,40,50,60,70)
  ) +
  scale_y_discrete(
    name = "machine learning models", 
    position = "right"
  ) +
  facet_grid(mutation ~ preprocessing, switch = 'y') +
  ggtitle("A")
p
plot(p)
ggsave("Demo-heatmap-accuracy.tiff",device="tiff",width=21,height=29,units="cm",dpi="retina")
ggsave("Demo-heatmap-accuracy.pdf",device="pdf",width=21,height=29,units="cm",dpi="retina")
ggsave("Demo-heatmap-accuracy.png",device="png",width=21,height=29,units="cm",dpi="retina")
dev.off()

## for F1 score
p = ggplot(data = data_long, aes(x = splitting, y = model)) + 
  theme_light(base_size = 16) + 
  geom_tile(aes(fill = F1_score)) + 
  geom_text(aes(label = labels03), size = 3) + 
  scale_fill_gradient(
    low = "#ff0000", 
    high = "#008000", 
    name = "F1 score", 
    limits = c(round(min(data_long$F1_score), digits=2)-0.01, round(max(data_long$F1_score), digits=2)+0.01), 
    breaks = c(round(min(data_long$F1_score), digits=2), 
               round((max(data_long$F1_score)-min(data_long$F1_score))/2+min(data_long$F1_score), digits=2), 
               round(max(data_long$F1_score), digits=2))
  ) + 
  theme(
    legend.position = "top", 
    legend.title = element_text(colour="#000000", size=12, face="bold"), 
    legend.text = element_text(size=12), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.text.x = element_text (color = "#000000", size = 12, angle = 0),
    axis.text.y = element_text (color = "#000000", size = 12, angle = 0), 
    strip.text.x = element_text(size = 10, face = "bold", color = "#000000"),
    strip.text.y = element_text(size = 10, face = "bold", color = "#000000"),
    strip.background = element_rect(colour="#000000", fill="#A9A9A9"), 
    panel.border = element_rect(colour="#FFFFFF"),
    axis.ticks = element_line(colour="#FFFFFF"), 
    plot.title = element_text(size = 15, face = "bold")
  ) +
  scale_x_continuous(
    name = "splitting proportion of the training dataset (%)", 
    limits = c(25, 75), breaks = c(30,40,50,60,70)
  ) +
  scale_y_discrete(
    name = "machine learning models", 
    position = "right"
  ) +
  facet_grid(mutation ~ preprocessing, switch = 'y') + 
  ggtitle("B")
p
plot(p)
ggsave("Demo-heatmap-F1score.tiff",device="tiff",width=21,height=29,units="cm",dpi="retina")
ggsave("Demo-heatmap-F1score.pdf",device="pdf",width=21,height=29,units="cm",dpi="retina")
ggsave("Demo-heatmap-F1score.png",device="png",width=21,height=29,units="cm",dpi="retina")
dev.off()

## for Cohen_kappa_testing
p = ggplot(data = data_long, aes(x = splitting, y = model)) + 
  theme_light(base_size = 16) + 
  geom_tile(aes(fill = Cohen_kappa_testing)) + 
  geom_text(aes(label = labels02), size = 3) + 
  scale_fill_gradient(
    low = "#ff0000", 
    high = "#008000", 
    name = "average Cohen kappa from the testing dataset", 
    limits = c(round(min(data_long$Cohen_kappa_testing), digits=2)-0.01, round(max(data_long$Cohen_kappa_testing), digits=2)+0.01), 
    breaks = c(round(min(data_long$Cohen_kappa_testing), digits=2), 
               round((max(data_long$Cohen_kappa_testing)-min(data_long$Cohen_kappa_testing))/2+min(data_long$Cohen_kappa_testing), digits=2), 
               round(max(data_long$Cohen_kappa_testing), digits=2))
  ) + 
  theme(
    legend.position = "top", 
    legend.title = element_text(colour="#000000", size=12, face="bold"), 
    legend.text = element_text(size=12), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.text.x = element_text (color = "#000000", size = 12, angle = 0),
    axis.text.y = element_text (color = "#000000", size = 12, angle = 0), 
    strip.text.x = element_text(size = 10, face = "bold", color = "#000000"),
    strip.text.y = element_text(size = 10, face = "bold", color = "#000000"),
    strip.background = element_rect(colour="#000000", fill="#A9A9A9"), 
    panel.border = element_rect(colour="#FFFFFF"),
    axis.ticks = element_line(colour="#FFFFFF"), 
    plot.title = element_text(size = 15, face = "bold")
  ) +
  scale_x_continuous(
    name = "splitting proportion of the training dataset (%)", 
    limits = c(25, 75), breaks = c(30,40,50,60,70)
  ) +
  scale_y_discrete(
    name = "machine learning models", 
    position = "right"
  ) +
  facet_grid(mutation ~ preprocessing, switch = 'y') + 
  ggtitle("C")
p
plot(p)
ggsave("Demo-heatmap-kappa.tiff",device="tiff",width=21,height=29,units="cm",dpi="retina")
ggsave("Demo-heatmap-kappa.pdf",device="pdf",width=21,height=29,units="cm",dpi="retina")
ggsave("Demo-heatmap-kappa.png",device="png",width=21,height=29,units="cm",dpi="retina")
dev.off()

## for ROC_AUC
p = ggplot(data = data_long, aes(x = splitting, y = model)) + 
  theme_light(base_size = 16) + 
  geom_tile(aes(fill = ROC_AUC)) + 
  geom_text(aes(label = labelsROC), size = 3) + 
  scale_fill_gradient(
    low = "#ff0000", 
    high = "#008000", 
    name = "area under the curve of the \nreceiver operating characteristic curve \n(%)", 
    limits = c(round(min(data_long$ROC_AUC), digits=2)-0.01, round(max(data_long$ROC_AUC), digits=2)+0.01),
    breaks = c(round(min(data_long$ROC_AUC), digits=2), 
               round((max(data_long$ROC_AUC)-min(data_long$ROC_AUC))/2+min(data_long$ROC_AUC), digits=2), 
               round(max(data_long$ROC_AUC), digits=2))
  ) + 
  theme(
    legend.position = "top", 
    legend.title = element_text(colour="#000000", size=12, face="bold"), 
    legend.text = element_text(size=12), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.text.x = element_text (color = "#000000", size = 12, angle = 0),
    axis.text.y = element_text (color = "#000000", size = 12, angle = 0), 
    strip.text.x = element_text(size = 10, face = "bold", color = "#000000"),
    strip.text.y = element_text(size = 10, face = "bold", color = "#000000"),
    strip.background = element_rect(colour="#000000", fill="#A9A9A9"), 
    panel.border = element_rect(colour="#FFFFFF"),
    axis.ticks = element_line(colour="#FFFFFF"), 
    plot.title = element_text(size = 15, face = "bold")
  ) +
  scale_x_continuous(
    name = "splitting proportion of the training dataset (%)", 
    limits = c(25, 75), breaks = c(30,40,50,60,70)
  ) +
  scale_y_discrete(
    name = "machine learning models", 
    position = "right"
  ) +
  facet_grid(mutation ~ preprocessing, switch = 'y') + 
  ggtitle("A")
p
plot(p)
ggsave("Demo-heatmap-ROC_AUC.tiff",device="tiff",width=21,height=29,units="cm",dpi="retina")
ggsave("Demo-heatmap-ROC_AUC.pdf",device="pdf",width=21,height=29,units="cm",dpi="retina")
ggsave("Demo-heatmap-ROC_AUC.png",device="png",width=21,height=29,units="cm",dpi="retina")
dev.off()

## for PR_AUC
p = ggplot(data = data_long, aes(x = splitting, y = model)) + 
  theme_light(base_size = 16) + 
  geom_tile(aes(fill = PR_AUC)) + 
  geom_text(aes(label = labelsPR), size = 3) + 
  scale_fill_gradient(
    low = "#ff0000", 
    high = "#008000", 
    name = "area under the curve of the \nprecision recall curve \n(%)", 
    limits = c(round(min(data_long$PR_AUC), digits=2)-0.01, round(max(data_long$PR_AUC), digits=2)+0.01),
    breaks = c(round(min(data_long$PR_AUC), digits=2), 
               round((max(data_long$PR_AUC)-min(data_long$PR_AUC))/2+min(data_long$PR_AUC), digits=2), 
               round(max(data_long$PR_AUC), digits=2))
  ) + 
  theme(
    legend.position = "top", 
    legend.title = element_text(colour="#000000", size=12, face="bold"), 
    legend.text = element_text(size=12), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.text.x = element_text (color = "#000000", size = 12, angle = 0),
    axis.text.y = element_text (color = "#000000", size = 12, angle = 0), 
    strip.text.x = element_text(size = 10, face = "bold", color = "#000000"),
    strip.text.y = element_text(size = 10, face = "bold", color = "#000000"),
    strip.background = element_rect(colour="#000000", fill="#A9A9A9"), 
    panel.border = element_rect(colour="#FFFFFF"),
    axis.ticks = element_line(colour="#FFFFFF"), 
    plot.title = element_text(size = 15, face = "bold")
  ) +
  scale_x_continuous(
    name = "splitting proportion of the training dataset (%)", 
    limits = c(25, 75), breaks = c(30,40,50,60,70)
  ) +
  scale_y_discrete(
    name = "machine learning models", 
    position = "right"
  ) +
  facet_grid(mutation ~ preprocessing, switch = 'y') + 
  ggtitle("B")
p
plot(p)
ggsave("Demo-heatmap-PR_AUC.tiff",device="tiff",width=21,height=29,units="cm",dpi="retina")
ggsave("Demo-heatmap-PR_AUC.pdf",device="pdf",width=21,height=29,units="cm",dpi="retina")
ggsave("Demo-heatmap-PR_AUC.png",device="png",width=21,height=29,units="cm",dpi="retina")
dev.off()

## for PRG_AUC
p = ggplot(data = data_long, aes(x = splitting, y = model)) + 
  theme_light(base_size = 16) + 
  geom_tile(aes(fill = PRG_AUC)) + 
  geom_text(aes(label = labelsPRG), size = 3) + 
  scale_fill_gradient(
    low = "#ff0000", 
    high = "#008000", 
    name = "area under the curve of the \nprecision recall gain curve \n(%)", 
    limits = c(round(min(data_long$PRG_AUC), digits=2)-0.01, round(max(data_long$PRG_AUC), digits=2)+0.01),
    breaks = c(round(min(data_long$PRG_AUC), digits=2), 
               round((max(data_long$PRG_AUC)-min(data_long$PRG_AUC))/2+min(data_long$PRG_AUC), digits=2), 
               round(max(data_long$PRG_AUC), digits=2))
  ) + 
  theme(
    legend.position = "top", 
    legend.title = element_text(colour="#000000", size=12, face="bold"), 
    legend.text = element_text(size=12), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.text.x = element_text (color = "#000000", size = 12, angle = 0),
    axis.text.y = element_text (color = "#000000", size = 12, angle = 0), 
    strip.text.x = element_text(size = 10, face = "bold", color = "#000000"),
    strip.text.y = element_text(size = 10, face = "bold", color = "#000000"),
    strip.background = element_rect(colour="#000000", fill="#A9A9A9"), 
    panel.border = element_rect(colour="#FFFFFF"),
    axis.ticks = element_line(colour="#FFFFFF"), 
    plot.title = element_text(size = 15, face = "bold")
  ) +
  scale_x_continuous(
    name = "splitting proportion of the training dataset (%)", 
    limits = c(25, 75), breaks = c(30,40,50,60,70)
  ) +
  scale_y_discrete(
    name = "machine learning models", 
    position = "right"
  ) +
  facet_grid(mutation ~ preprocessing, switch = 'y') + 
  ggtitle("C")
p
plot(p)
ggsave("Demo-heatmap-PRG_AUC.tiff",device="tiff",width=21,height=29,units="cm",dpi="retina")
ggsave("Demo-heatmap-PRG_AUC.pdf",device="pdf",width=21,height=29,units="cm",dpi="retina")
ggsave("Demo-heatmap-PRG_AUC.png",device="png",width=21,height=29,units="cm",dpi="retina")
dev.off()

## for seconds
p = ggplot(data = data_long, aes(x = splitting, y = model)) + 
  theme_light(base_size = 16) + 
  geom_tile(aes(fill = seconds)) + 
  geom_text(aes(label = exec_t04), size = 3) + 
  scale_fill_gradient(
    high = "#ff0000", 
    low = "#008000", 
    name = "execution time \n(seconds)", 
    limits = c(round(min(data_long$seconds), digits=0)-5, round(max(data_long$seconds), digits=0))+1, 
    breaks = c(round(min(data_long$seconds), digits=0), round((max(data_long$seconds)-min(data_long$seconds))/2+min(data_long$seconds), digits=0), round(max(data_long$seconds), digits=0))
  ) + 
  theme(
    legend.position = "top", 
    legend.title = element_text(colour="#000000", size=12, face="bold"), 
    legend.text = element_text(size=12), 
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.text.x = element_text (color = "#000000", size = 12, angle = 0),
    axis.text.y = element_text (color = "#000000", size = 12, angle = 0), 
    strip.text.x = element_text(size = 10, face = "bold", color = "#000000"),
    strip.text.y = element_text(size = 10, face = "bold", color = "#000000"),
    strip.background = element_rect(colour="#000000", fill="#A9A9A9"), 
    panel.border = element_rect(colour="#FFFFFF"),
    axis.ticks = element_line(colour="#FFFFFF"), 
    plot.title = element_text(size = 15, face = "bold")
  ) +
  scale_x_continuous(
    name = "splitting proportion of the training dataset (%)", 
    limits = c(25, 75), breaks = c(30,40,50,60,70)
  ) +
  scale_y_discrete(
    name = "machine learning models", 
    position = "right"
  ) +
  facet_grid(mutation ~ preprocessing, switch = 'y') + 
  ggtitle("D")
p
plot(p)
ggsave("Demo-heatmap-time.tiff",device="tiff",width=21,height=29,units="cm",dpi="retina")
ggsave("Demo-heatmap-time.pdf",device="pdf",width=21,height=29,units="cm",dpi="retina")
ggsave("Demo-heatmap-time.png",device="png",width=21,height=29,units="cm",dpi="retina")
dev.off()


# Pearson correlation tests

## average accuracy training vs kappa training
### R2
cor(x=data_long$average_accuracy_training, y=data_long$Cohen_kappa_training, method="pearson")
### p-value
cor.test(x=data_long$average_accuracy_training, y=data_long$Cohen_kappa_training, method="pearson")

## average accuracy training vs F1 score
### R2
cor(x=data_long$average_accuracy_training, y=data_long$F1_score, method="pearson")
### p-value
cor.test(x=data_long$average_accuracy_training, y=data_long$F1_score, method="pearson")

## average accuracy training vs ROC-AUC
### R2
cor(x=data_long$average_accuracy_training, y=data_long$ROC_AUC, method="pearson")
### p-value
cor.test(x=data_long$average_accuracy_training, y=data_long$ROC_AUC, method="pearson")

## average accuracy training vs PR-AUC
### R2
cor(x=data_long$average_accuracy_training, y=data_long$PR_AUC, method="pearson")
### p-value
cor.test(x=data_long$average_accuracy_training, y=data_long$PR_AUC, method="pearson")

## average accuracy training vs PRG-AUC
### R2
cor(x=data_long$average_accuracy_training, y=data_long$PRG_AUC, method="pearson")
### p-value
cor.test(x=data_long$average_accuracy_training, y=data_long$PRG_AUC, method="pearson")

## average accuracy testing vs kappa testing
### R2
cor(x=data_long$average_accuracy_testing, y=data_long$Cohen_kappa_testing, method="pearson")
### p-value
cor.test(x=data_long$average_accuracy_testing, y=data_long$Cohen_kappa_testing, method="pearson")


# correlation plots

## reorganize levels of variables
levels(data_long$model)
data_long$model <- factor(data_long$model, levels=c("mod_A", "mod_B", "mod_C", "mod_D", "mod_E"))
levels(data_long$model)

## average_accuracy_training versus average_accuracy_testing
my.formula <- y ~ x
p = ggplot(data = data_long, aes(x = average_accuracy_training, y = average_accuracy_testing, group=model, color=model)) +
  theme_light(base_size = 16) +
  geom_point(aes(shape = model, color = model, size = model)) +
  geom_smooth(aes(shape = model, color = model), size = 0.5, fill = "#A9A9A9", method=lm, linetype="dashed", se=FALSE, formula = my.formula) +
  scale_shape_manual(values=c(20, 20, 20, 20, 20, 20)) +
  scale_color_manual(values=c('#FF0000','#0000FF','#008000','#EE82EE','#FFA500','#582900')) +
  scale_size_manual(values=c(3,3,3,3,3,3)) +
  theme(legend.position="right", 
        plot.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15)) + 
  scale_x_continuous(name = "Average accuracy \n of the training dataset (%)", limits = c(0, 105), breaks = c(0,20,40,60,80,100)) +
  scale_y_continuous(name = "Average accuracy \n of the testing dataset (%)", limits = c(0, 80), breaks = c(0,20,40,60,80)) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE,label.y = "top", label.x = "left", rr.digits = 3, coef.digits = 3, f.digits = 6) + 
  ggtitle("A")
p
plot(p)
ggsave("Demo-correlation-accuracy_training_vs_accuracy_testing.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("Demo-correlation-accuracy_training_vs_accuracy_testing.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
ggsave("Demo-correlation-accuracy_training_vs_accuracy_testing.png",device="png",width=17,height=17,units="cm",dpi="retina")
dev.off()

## average_accuracy_training versus Cohen_kappa_training
my.formula <- y ~ x
p = ggplot(data = data_long, aes(x = average_accuracy_training, y = Cohen_kappa_training, group=model, color=model)) +
  theme_light(base_size = 16) +
  geom_point(aes(shape = model, color = model, size = model)) +
  geom_smooth(aes(shape = model, color = model), size = 0.5, fill = "#A9A9A9", method=lm, linetype="dashed", se=FALSE, formula = my.formula) +
  scale_shape_manual(values=c(20, 20, 20, 20, 20, 20)) +
  scale_color_manual(values=c('#FF0000','#0000FF','#008000','#EE82EE','#FFA500','#582900')) +
  scale_size_manual(values=c(3,3,3,3,3,3)) +
  theme(legend.position="right", 
        plot.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15)) + 
  scale_x_continuous(name = "Average accuracy \n of the training dataset (%)", limits = c(0, 110), breaks = c(0,20,40,60,80,100)) +
  scale_y_continuous(name = "Cohen kappa \n of the training dataset", limits = c(0, 1.0), breaks = c(0.0,0.2,0.4,0.6,0.8,1)) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE,label.y = "bottom", label.x = "right", rr.digits = 3, coef.digits = 3, f.digits = 6) + 
  ggtitle("B")
p
plot(p)
ggsave("Demo-correlation-accuracy_training_vs_kappa_training.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("Demo-correlation-accuracy_training_vs_kappa_training.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
ggsave("Demo-correlation-accuracy_training_vs_kappa_training.png",device="png",width=17,height=17,units="cm",dpi="retina")
dev.off()

## average_accuracy_training versus F1_score
my.formula <- y ~ x
p = ggplot(data = data_long, aes(x = average_accuracy_training, y = F1_score, group=model, color=model)) +
  theme_light(base_size = 16) +
  geom_point(aes(shape = model, color = model, size = model)) +
  geom_smooth(aes(shape = model, color = model), size = 0.5, fill = "#A9A9A9", method=lm, linetype="dashed", se=FALSE, formula = my.formula) +
  scale_shape_manual(values=c(20, 20, 20, 20, 20, 20)) +
  scale_color_manual(values=c('#FF0000','#0000FF','#008000','#EE82EE','#FFA500','#582900')) +
  scale_size_manual(values=c(3,3,3,3,3,3)) +
  theme(legend.position="right", 
        plot.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15)) + 
  scale_x_continuous(name = "Average accuracy \n of the training dataset (%)", limits = c(40, 100), breaks = c(40,60,80,100)) +
  scale_y_continuous(name = "F1-score", limits = c(0.2, 1.0), breaks = c(0.2,0.4,0.6,0.8,1.0)) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE,label.y = "bottom", label.x = "right", rr.digits = 3, coef.digits = 3, f.digits = 6) + 
  ggtitle("C")
p
plot(p)
ggsave("Demo-correlation-accuracy_training_vs_F1_score.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("Demo-correlation-accuracy_training_vs_F1_score.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
ggsave("Demo-correlation-accuracy_training_vs_F1_score.png",device="png",width=17,height=17,units="cm",dpi="retina")
dev.off()

## average_accuracy_training versus ROC-AUC
my.formula <- y ~ x
p = ggplot(data = data_long, aes(x = average_accuracy_training, y = ROC_AUC, group=model, color=model)) +
  theme_light(base_size = 16) +
  geom_point(aes(shape = model, color = model, size = model)) +
  geom_smooth(aes(shape = model, color = model), size = 0.5, fill = "#A9A9A9", method=lm, linetype="dashed", se=FALSE, formula = my.formula) +
  scale_shape_manual(values=c(20, 20, 20, 20, 20, 20)) +
  scale_color_manual(values=c('#FF0000','#0000FF','#008000','#EE82EE','#FFA500','#582900')) +
  scale_size_manual(values=c(3,3,3,3,3,3)) +
  theme(legend.position="right", 
        plot.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15)) + 
  scale_x_continuous(name = "Average accuracy \n of the training dataset (%)", limits = c(40, 100), breaks = c(40,60,80,100)) +
  scale_y_continuous(name = "ROC-AUC (%)", limits = c(10, 100), breaks = c(20,40,60,80,100)) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE,label.y = "bottom", label.x = "right", rr.digits = 3, coef.digits = 3, f.digits = 6) + 
  ggtitle("D")
p
plot(p)
ggsave("Demo-correlation-accuracy_training_vs_ROC-AUC.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("Demo-correlation-accuracy_training_vs_ROC-AUC.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
ggsave("Demo-correlation-accuracy_training_vs_ROC-AUC.png",device="png",width=17,height=17,units="cm",dpi="retina")
dev.off()

## average_accuracy_training versus PR-AUC
my.formula <- y ~ x
p = ggplot(data = data_long, aes(x = average_accuracy_training, y = PR_AUC, group=model, color=model)) +
  theme_light(base_size = 16) +
  geom_point(aes(shape = model, color = model, size = model)) +
  geom_smooth(aes(shape = model, color = model), size = 0.5, fill = "#A9A9A9", method=lm, linetype="dashed", se=FALSE, formula = my.formula) +
  scale_shape_manual(values=c(20, 20, 20, 20, 20, 20)) +
  scale_color_manual(values=c('#FF0000','#0000FF','#008000','#EE82EE','#FFA500','#582900')) +
  scale_size_manual(values=c(3,3,3,3,3,3)) +
  theme(legend.position="right", 
        plot.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15)) + 
  scale_x_continuous(name = "Average accuracy \n of the training dataset (%)", limits = c(40, 100), breaks = c(40,60,80,100)) +
  scale_y_continuous(name = "PR-AUC (%)", limits = c(10, 100), breaks = c(20,40,60,80,100)) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE,label.y = "bottom", label.x = "right", rr.digits = 3, coef.digits = 3, f.digits = 6) + 
  ggtitle("E")
p
plot(p)
ggsave("Demo-correlation-accuracy_training_vs_PR-AUC.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("Demo-correlation-accuracy_training_vs_PR-AUC.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
ggsave("Demo-correlation-accuracy_training_vs_PR-AUC.png",device="png",width=17,height=17,units="cm",dpi="retina")
dev.off()

## average_accuracy_training versus PRG-AUC
my.formula <- y ~ x
p = ggplot(data = data_long, aes(x = average_accuracy_training, y = PRG_AUC, group=model, color=model)) +
  theme_light(base_size = 16) +
  geom_point(aes(shape = model, color = model, size = model)) +
  geom_smooth(aes(shape = model, color = model), size = 0.5, fill = "#A9A9A9", method=lm, linetype="dashed", se=FALSE, formula = my.formula) +
  scale_shape_manual(values=c(20, 20, 20, 20, 20, 20)) +
  scale_color_manual(values=c('#FF0000','#0000FF','#008000','#EE82EE','#FFA500','#582900')) +
  scale_size_manual(values=c(3,3,3,3,3,3)) +
  theme(legend.position="right", 
        plot.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15)) + 
  scale_x_continuous(name = "Average accuracy \n of the training dataset (%)", limits = c(40, 100), breaks = c(40,60,80,100)) +
  scale_y_continuous(name = "PRG-AUC (%)", limits = c(10, 100), breaks = c(20,40,60,80,100)) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE,label.y = "top", label.x = "left", rr.digits = 3, coef.digits = 3, f.digits = 6) + 
  ggtitle("F")
p
plot(p)
ggsave("Demo-correlation-accuracy_training_vs_PRG-AUC.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("Demo-correlation-accuracy_training_vs_PRG-AUC.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
ggsave("Demo-correlation-accuracy_training_vs_PRG-AUC.png",device="png",width=17,height=17,units="cm",dpi="retina")
dev.off()

## average_accuracy_testing versus Cohen_kappa_testing
my.formula <- y ~ x
p = ggplot(data = data_long, aes(x = average_accuracy_testing, y = Cohen_kappa_testing, group=model, color=model)) +
  theme_light(base_size = 16) +
  geom_point(aes(shape = model, color = model, size = model)) +
  geom_smooth(aes(shape = model, color = model), size = 0.5, fill = "#A9A9A9", method=lm, linetype="dashed", se=FALSE, formula = my.formula) +
  scale_shape_manual(values=c(20, 20, 20, 20, 20, 20)) +
  scale_color_manual(values=c('#FF0000','#0000FF','#008000','#EE82EE','#FFA500','#582900')) +
  scale_size_manual(values=c(3,3,3,3,3,3)) +
  theme(legend.position="right", 
        plot.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15)) + 
  scale_x_continuous(name = "Average accuracy \n of the testing dataset (%)", limits = c(40, 100), breaks = c(40,60,80,100)) +
  scale_y_continuous(name = "Cohen kappa \n of the testing dataset", limits = c(0.2, 1.0), breaks = c(0.2,0.4,0.6,0.8,1.0)) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE,label.y = "bottom", label.x = "right", rr.digits = 3, coef.digits = 3, f.digits = 6) + 
  ggtitle("G")
p
plot(p)
ggsave("Demo-correlation-average_accuracy_testing_vs_Cohen_kappa_testing.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("Demo-correlation-average_accuracy_testing_ vs_Cohen_kappa_testing.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
ggsave("Demo-correlation-average_accuracy_testing_ vs_Cohen_kappa_testing.png",device="png",width=17,height=17,units="cm",dpi="retina")
dev.off()

## average_accuracy_testing versus F1_score
my.formula <- y ~ x
p = ggplot(data = data_long, aes(x = average_accuracy_testing, y = F1_score, group=model, color=model)) +
  theme_light(base_size = 16) +
  geom_point(aes(shape = model, color = model, size = model)) +
  geom_smooth(aes(shape = model, color = model), size = 0.5, fill = "#A9A9A9", method=lm, linetype="dashed", se=FALSE, formula = my.formula) +
  scale_shape_manual(values=c(20, 20, 20, 20, 20, 20)) +
  scale_color_manual(values=c('#FF0000','#0000FF','#008000','#EE82EE','#FFA500','#582900')) +
  scale_size_manual(values=c(3,3,3,3,3,3)) +
  theme(legend.position="right", 
        plot.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15)) + 
  scale_x_continuous(name = "Average accuracy \n of the training dataset (%)", limits = c(40, 100), breaks = c(40,60,80,100)) +
  scale_y_continuous(name = "F1-score", limits = c(0.2, 1.0), breaks = c(0.2,0.4,0.6,0.8,1.0)) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE,label.y = "bottom", label.x = "right", rr.digits = 3, coef.digits = 3, f.digits = 6) + 
  ggtitle("H")
p
plot(p)
ggsave("Demo-correlation-average_accuracy_testing_vs_F1_score.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("Demo-correlation-average_accuracy_testing_vs_F1_score.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
ggsave("Demo-correlation-average_accuracy_testing_vs_F1_score.png",device="png",width=17,height=17,units="cm",dpi="retina")
dev.off()

## average_accuracy_testing versus ROC-AUC
my.formula <- y ~ x
p = ggplot(data = data_long, aes(x = average_accuracy_testing, y = ROC_AUC, group=model, color=model)) +
  theme_light(base_size = 16) +
  geom_point(aes(shape = model, color = model, size = model)) +
  geom_smooth(aes(shape = model, color = model), size = 0.5, fill = "#A9A9A9", method=lm, linetype="dashed", se=FALSE, formula = my.formula) +
  scale_shape_manual(values=c(20, 20, 20, 20, 20, 20)) +
  scale_color_manual(values=c('#FF0000','#0000FF','#008000','#EE82EE','#FFA500','#582900')) +
  scale_size_manual(values=c(3,3,3,3,3,3)) +
  theme(legend.position="right", 
        plot.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15)) + 
  scale_x_continuous(name = "Average accuracy \n of the testing dataset (%)", limits = c(40, 100), breaks = c(40,60,80,100)) +
  scale_y_continuous(name = "ROC-AUC (%)", limits = c(10, 100), breaks = c(20,40,60,80,100)) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE,label.y = "bottom", label.x = "right", rr.digits = 3, coef.digits = 3, f.digits = 6) + 
  ggtitle("I")
p
plot(p)
ggsave("Demo-correlation-accuracy_testing_vs_ROC-AUC.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("Demo-correlation-accuracy_testing_vs_ROC-AUC.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
ggsave("Demo-correlation-accuracy_testing_vs_ROC-AUC.png",device="png",width=17,height=17,units="cm",dpi="retina")
dev.off()

## average_accuracy_testing versus PR-AUC
my.formula <- y ~ x
p = ggplot(data = data_long, aes(x = average_accuracy_testing, y = PR_AUC, group=model, color=model)) +
  theme_light(base_size = 16) +
  geom_point(aes(shape = model, color = model, size = model)) +
  geom_smooth(aes(shape = model, color = model), size = 0.5, fill = "#A9A9A9", method=lm, linetype="dashed", se=FALSE, formula = my.formula) +
  scale_shape_manual(values=c(20, 20, 20, 20, 20, 20)) +
  scale_color_manual(values=c('#FF0000','#0000FF','#008000','#EE82EE','#FFA500','#582900')) +
  scale_size_manual(values=c(3,3,3,3,3,3)) +
  theme(legend.position="right", 
        plot.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15)) + 
  scale_x_continuous(name = "Average accuracy \n of the testing dataset (%)", limits = c(40, 100), breaks = c(40,60,80,100)) +
  scale_y_continuous(name = "PR-AUC (%)", limits = c(10, 100), breaks = c(20,40,60,80,100)) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE,label.y = "bottom", label.x = "right", rr.digits = 3, coef.digits = 3, f.digits = 6) + 
  ggtitle("J")
p
plot(p)
ggsave("Demo-correlation-accuracy_testing_vs_PR-AUC.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("Demo-correlation-accuracy_testing_vs_PR-AUC.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
ggsave("Demo-correlation-accuracy_testing_vs_PR-AUC.png",device="png",width=17,height=17,units="cm",dpi="retina")
dev.off()

## average_accuracy_testing versus PRG-AUC
my.formula <- y ~ x
p = ggplot(data = data_long, aes(x = average_accuracy_testing, y = PRG_AUC, group=model, color=model)) +
  theme_light(base_size = 16) +
  geom_point(aes(shape = model, color = model, size = model)) +
  geom_smooth(aes(shape = model, color = model), size = 0.5, fill = "#A9A9A9", method=lm, linetype="dashed", se=FALSE, formula = my.formula) +
  scale_shape_manual(values=c(20, 20, 20, 20, 20, 20)) +
  scale_color_manual(values=c('#FF0000','#0000FF','#008000','#EE82EE','#FFA500','#582900')) +
  scale_size_manual(values=c(3,3,3,3,3,3)) +
  theme(legend.position="right", 
        plot.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15)) + 
  scale_x_continuous(name = "Average accuracy \n of the testing dataset (%)", limits = c(40, 100), breaks = c(40,60,80,100)) +
  scale_y_continuous(name = "PRG-AUC (%)", limits = c(10, 100), breaks = c(20,40,60,80,100)) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")), 
               parse = TRUE,label.y = "top", label.x = "left", rr.digits = 3, coef.digits = 3, f.digits = 6) + 
  ggtitle("K")
p
plot(p)
ggsave("Demo-correlation-accuracy_testing_vs_PRG-AUC.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("Demo-correlation-accuracy_testing_vs_PRG-AUC.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")
ggsave("Demo-correlation-accuracy_testing_vs_PRG-AUC.png",device="png",width=17,height=17,units="cm",dpi="retina")
dev.off()