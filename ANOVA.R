#### ANOVA.R ####
#### R script to perform statistical analyses on prediction results for the research article ####
#### "Source attribution of Listeria monocytogenes through a versatile supervised machine learning workflow based on genomic data" ####
#### performed with R version 4.2.3 (2023-03-15) and RStudio 2022.12.0+353 on Ubuntu 20.04 5 LTS "Focal Fossa" ####

# clean environment
rm(list=ls())

# clean graphical device
graphics.off()

# set working directory
setwd("/home/IZSNT/p.castelli/Documenti/work/ML-Pyseer/R-prediction-res")

# call library
library(ggplot2) # ‘3.4.2’
library(plyr) # ‘1.8.8’
library(ggpmisc) # ‘0.5.2’
library(reshape2) #  ‘1.4.4’
library(lubridate)  # ‘1.9.2’

# read dataframe
data_raw = read.table("Metrics.tsv", dec = ".", header=TRUE, sep = "\t", quote = "")
data_raw = read.table("fake_Metrics.tsv", dec = ".", header=TRUE, sep = "\t", quote = "")
## check dimension
dim(data_raw)
#[1] 270  18
## check first 20 lines
head(data_raw, 20)

# check nature of variables (integer or factor)
str(data_raw)

# transform scientific numbers into decimal numbers and fix column names
## format numbers in columns
new_col05 <- as.numeric(format(data_raw$average_accuracy_training, digits=3))
new_col06 <- format(data_raw$CI95_accuracy_training, digits=3)
new_col07 <- as.numeric(format(data_raw$average_accuracy_testing, digits=3))
new_col08 <- format(data_raw$CI95_accuracy_testing, digits=3)
new_col09 <- as.numeric(format(data_raw$F1.score, digits=3))
new_col10 <- as.numeric(format(data_raw$Cohen_kappa_training, digits=3))
new_col11 <- as.numeric(format(data_raw$Cohen_kappa_testing, digits=3))
new_col12 <- as.numeric(format(data_raw$ROC_AUC, digits=3))
new_col13 <- as.numeric(format(data_raw$PR_AUC, digits=3))
new_col14 <- as.numeric(format(data_raw$PRG_AUC, digits=3))
new_col15 <- as.numeric(format(data_raw$seconds, digits=5))
## create new dataframe with desired value format
data_long = data.frame(data_raw$mutation,data_raw$preprocessing,data_raw$splitting,data_raw$model,new_col05,new_col06,new_col07,new_col08,new_col09,new_col10,new_col11,new_col12,new_col13,new_col14,new_col15,data_raw$after_constant_descr_removal,data_raw$after_near_zero.variance_descr_removal,data_raw$selected_predictors)
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
colnames(data_long)[16] <- "after_constant_descr_removal"
colnames(data_long)[17] <- "after_near_zero-variance_descr_removal"
colnames(data_long)[18] <- "selected_predictors"
## save table as file
#write.table(df, file = "Metrics_2.tsv", quote=FALSE, sep='\t', row.names = FALSE)

# check dimension
dim(data_long)
#[1] 270  18

# check first 10 lines and last 10 lines
head(data_long, 10)
tail(data_long,10)

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
                                "MLST"="7_locus_alleles",
                                "coregenome_alleles"="core_alleles",
                                "accessory_genes"="acc_genes",
                                "coregenome_variants"="core_SNPs",
                                "pangenome_kmer"="acc_kmers"
                              ))
levels(data_long$mutation)

# reorganize levels of variables (mutations)
levels(data_long$mutation)
data_long$mutation <- factor(data_long$mutation, levels=c("7_locus_alleles", "core_alleles", "acc_genes", "core_SNPs", "acc_kmers"))
levels(data_long$mutation)

# reorganize levels of variables (pre-processing)
data_long$preprocessing <- as.factor(data_long$preprocessing)
levels(data_long$preprocessing)
data_long$preprocessing <- factor(data_long$preprocessing, levels=c("no_preprocessing", "near-zero_variances"))
levels(data_long$preprocessing)

# reorganize levels of variables (models)
data_long$model <- as.factor(data_long$model)
levels(data_long$model)
data_long$model <- factor(data_long$model, levels=c("XGB", "SVM", "SGB", "RF", "ERT", "BLR"))
levels(data_long$model)

# convert splitting to factor
data_long$splitting <- as.factor(data_long$splitting)
str(data_long)


# calculate accuracy mean and sd
library(dplyr) # ‘1.1.2’
# calculate accuracy mean and sd through mutation
df.mutation <- group_by(data_long, mutation) %>%
  summarise(
    count = n(),
    mean = mean(average_accuracy_testing, na.rm = TRUE),
    sd = sd(average_accuracy_testing, na.rm = TRUE)) %>%
  mutate(se = sd / sqrt(count),
         lower.ci = mean - qt(1 - (0.05 / 2), count - 1) * se,
         upper.ci = mean + qt(1 - (0.05 / 2), count - 1) * se)
df.mutation <- as.data.frame(df.mutation)
#mutation count     mean       sd        se lower.ci upper.ci
#1 7_locus_alleles    30 53.68000 4.919588 0.8981898 51.84300 55.51700
#2    core_alleles    60 65.72667 9.277709 1.1977471 63.32998 68.12335
#3       acc_genes    60 68.88667 4.534892 0.5854520 67.71518 70.05815
#4       core_SNPs    60 59.99000 9.624857 1.2425636 57.50364 62.47636
#5       acc_kmers    60 67.32500 6.096140 0.7870083 65.75020 68.89980
# calculate accuracy mean and sd through splitting
df.splitting <- group_by(data_long, splitting) %>%
  summarise(
    count = n(),
    mean = mean(average_accuracy_testing, na.rm = TRUE),
    sd = sd(average_accuracy_testing, na.rm = TRUE)) %>%
  mutate(se = sd / sqrt(count),
         lower.ci = mean - qt(1 - (0.05 / 2), count - 1) * se,
         upper.ci = mean + qt(1 - (0.05 / 2), count - 1) * se)
df.splitting <- as.data.frame(df.splitting)
#splitting count     mean        sd       se lower.ci upper.ci
#1        50    54 60.50370  7.585213 1.032217 58.43334 62.57407
#2        60    54 64.38519  7.653891 1.041563 62.29608 66.47430
#3        70    54 65.48889  7.609742 1.035555 63.41183 67.56595
#4        80    54 67.64074  9.165985 1.247333 65.13891 70.14257
#5        90    54 62.83519 10.491491 1.427711 59.97156 65.69881
# calculate accuracy mean and sd through preprocessing
df.preprocessing <- group_by(data_long, preprocessing) %>%
  summarise(
    count = n(),
    mean = mean(average_accuracy_testing, na.rm = TRUE),
    sd = sd(average_accuracy_testing, na.rm = TRUE)) %>%
  mutate(se = sd / sqrt(count),
         lower.ci = mean - qt(1 - (0.05 / 2), count - 1) * se,
         upper.ci = mean + qt(1 - (0.05 / 2), count - 1) * se)
df.preprocessing <- as.data.frame(df.preprocessing)
#preprocessing count     mean       sd        se lower.ci upper.ci
#1    no_preprocessing   150 63.67267 9.472963 0.7734642 62.14429 65.20104
#2 near-zero_variances   120 64.79333 8.004408 0.7306991 63.34648 66.24019
# calculate accuracy mean and sd through model
df.model <- group_by(data_long, model) %>%
  summarise(
    count = n(),
    mean = mean(average_accuracy_testing, na.rm = TRUE),
    sd = sd(average_accuracy_testing, na.rm = TRUE)) %>%
  mutate(se = sd / sqrt(count),
         lower.ci = mean - qt(1 - (0.05 / 2), count - 1) * se,
         upper.ci = mean + qt(1 - (0.05 / 2), count - 1) * se)
df.model <- as.data.frame(df.model)
#model count     mean       sd        se lower.ci upper.ci
#1   XGB    45 67.83778 7.227765 1.0774516 65.66632 70.00924
#2   SVM    45 71.21111 7.198460 1.0730831 69.04845 73.37377
#3   SGB    45 64.03556 6.365573 0.9489237 62.12313 65.94799
#4    RF    45 56.46444 9.437261 1.4068238 53.62918 59.29971
#5   ERT    45 60.74444 8.684169 1.2945595 58.13543 63.35346
#6   BLR    45 64.73111 5.699315 0.8496037 63.01885 66.44337


# calculate accuracy mean and sd from preprocessing of SNPs
sbs.SNPs <- subset(data_long, data_long$mutation == "core_SNPs")

df.SNPs.preprocessing <- group_by(sbs.SNPs, preprocessing) %>%
  summarise(
    count = n(),
    mean = mean(average_accuracy_testing, na.rm = TRUE),
    sd = sd(average_accuracy_testing, na.rm = TRUE)) %>%
  mutate(se = sd / sqrt(count),
         lower.ci = mean - qt(1 - (0.05 / 2), count - 1) * se,
         upper.ci = mean + qt(1 - (0.05 / 2), count - 1) * se)
df.SNPs.preprocessing <- as.data.frame(df.SNPs.preprocessing)
#preprocessing count     mean        sd       se lower.ci upper.ci
#1    no_preprocessing    30 62.18333 11.399670 2.081285 57.92663 66.44004
#2 near-zero_variances    30 57.79667  6.968821 1.272327 55.19447 60.39887


### --------------------------------------
df.global <- group_by(data_long, splitting, preprocessing, model) %>%
  summarise(
    count = n(),
    mean = mean(average_accuracy_testing, na.rm = TRUE),
    sd = sd(average_accuracy_testing, na.rm = TRUE)) %>%
  mutate(se = sd / sqrt(count),
         lower.ci = mean - qt(1 - (0.05 / 2), count - 1) * se,
         upper.ci = mean + qt(1 - (0.05 / 2), count - 1) * se)
df.global <- as.data.frame(df.global)
### --------------------------------------


df <- subset(data_long, data_long$mutation != "7_locus_alleles")
dim(data_long)
dim(df)

# Levene’s test: Compare the variances of k samples, where k can be more than two samples. It’s an alternative to the Bartlett’s test that is less sensitive to departures from normality.
library(car) #  ‘3.1.2’
# Levene test
Levene.mutation <- leveneTest(average_accuracy_testing ~ mutation, df)
#Levene's Test for Homogeneity of Variance (center = median)
#       Df F value    Pr(>F)    
#group   3  14.437 1.148e-08 ***
#      236 
Levene.splitting <- leveneTest(average_accuracy_testing ~ splitting, df)
#Levene's Test for Homogeneity of Variance (center = median)
#       Df F value Pr(>F)
#group   4  1.8407 0.1218
#      235  
Levene.preprocessing <- leveneTest(average_accuracy_testing ~ preprocessing, df)
#Levene's Test for Homogeneity of Variance (center = median)
#       Df F value Pr(>F)
#group   1   0.048 0.8268
#      238 
Levene.model <- leveneTest(average_accuracy_testing ~ model, df)
#Levene's Test for Homogeneity of Variance (center = median)
#       Df F value    Pr(>F)    
#group   5  15.335 4.885e-13 ***
#      234  


levels(data_long$mutation)
dim(df)

# perform multi-way ANOVA test requiring assumptions of normal distribution of residuals and homogeneity of variance in all groups
mwANOVA.global <- aov(average_accuracy_testing ~ mutation + preprocessing + splitting + model, data = df)
results.mwANOVA.global <- summary(mwANOVA.global)
#Df Sum Sq Mean Sq F value   Pr(>F)    
#mutation        3   2713   904.2  33.565  < 2e-16 ***
#  preprocessing   1    114   113.9   4.226    0.041 *  
#  splitting       4   1450   362.5  13.457 7.45e-10 ***
#  model           5   6298  1259.6  46.756  < 2e-16 ***
#  Residuals     226   6088    26.9  

df.core_alleles <- subset(data_long, data_long$mutation == "core_alleles")
mwANOVA.core_alleles <- aov(average_accuracy_testing ~ preprocessing + splitting + model, data = df.core_alleles)
results.mwANOVA.core_alleles <- summary(mwANOVA.core_alleles)
#Df Sum Sq Mean Sq F value   Pr(>F)    
#preprocessing  1      2     1.9   0.306    0.583    
#splitting      4    345    86.2  14.073 1.02e-07 ***
#  model          5   4432   886.3 144.699  < 2e-16 ***
#  Residuals     49    300     6.1          
df.acc_genes <- subset(data_long, data_long$mutation == "acc_genes")
mwANOVA.acc_genes <- aov(average_accuracy_testing ~ preprocessing + splitting + model, data = df.acc_genes)
results.mwANOVA.acc_genes <- summary(mwANOVA.acc_genes)
#Df Sum Sq Mean Sq F value   Pr(>F)    
#preprocessing  1    0.0    0.00   0.001    0.982    
#splitting      4  366.8   91.71  18.113 3.43e-09 ***
#  model          5  598.4  119.68  23.636 5.23e-12 ***
#  Residuals     49  248.1    5.06 
df.core_SNPs <- subset(data_long, data_long$mutation == "core_SNPs")
mwANOVA.core_SNPs <- aov(average_accuracy_testing ~ preprocessing + splitting + model, data = df.core_SNPs)
results.mwANOVA.core_SNPs <- summary(mwANOVA.core_SNPs)
#Df Sum Sq Mean Sq F value   Pr(>F)    
#preprocessing  1    289   288.6  20.390 3.98e-05 ***
#  splitting      4    340    84.9   5.998  0.00052 ***
#  model          5   4144   828.7  58.543  < 2e-16 ***
#  Residuals     49    694    14.2   
df.acc_kmers <- subset(data_long, data_long$mutation == "acc_kmers")
mwANOVA.acc_kmers <- aov(average_accuracy_testing ~ preprocessing + splitting + model, data = df.acc_kmers)
results.mwANOVA.acc_kmers <- summary(mwANOVA.acc_kmers)
#Df Sum Sq Mean Sq F value   Pr(>F)    
#preprocessing  1   33.3   33.30   5.729   0.0206 *  
#  splitting      4  590.5  147.62  25.398 1.98e-11 ***
#  model          5 1284.0  256.80  44.182  < 2e-16 ***
#  Residuals     49  284.8    5.81   

# check normality assumption with Shapiro-Wilk test
# global
mwANOVA.global.residuals <- residuals(object = mwANOVA.global)
Shapiro.mwANOVA.global <- shapiro.test(x = mwANOVA.global.residuals)
#Shapiro-Wilk normality test
#data:  mwANOVA.global.residuals
#W = 0.99174, p-value = 0.1961

# core_alleles
mwANOVA.core_alleles.residuals <- residuals(object = mwANOVA.core_alleles)
Shapiro.mwANOVA.core_alleles <- shapiro.test(x = mwANOVA.core_alleles.residuals)
#Shapiro-Wilk normality test
#data:  mwANOVA.core_alleles.residuals
#W = 0.97062, p-value = 0.1565

# acc_genes
mwANOVA.acc_genes.residuals <- residuals(object = mwANOVA.acc_genes)
Shapiro.mwANOVA.acc_genes <- shapiro.test(x = mwANOVA.acc_genes.residuals)
#Shapiro-Wilk normality test
#data:  mwANOVA.acc_genes.residuals
#W = 0.98342, p-value = 0.5886

# core_SNPs
mwANOVA.core_SNPs.residuals <- residuals(object = mwANOVA.core_SNPs)
Shapiro.mwANOVA.core_SNPs <- shapiro.test(x = mwANOVA.core_SNPs.residuals)
#Shapiro-Wilk normality test
#data:  mwANOVA.core_SNPs.residuals
#W = 0.98279, p-value = 0.5571

# acc_kmers
mwANOVA.acc_kmers.residuals <- residuals(object = mwANOVA.acc_kmers)
Shapiro.mwANOVA.acc_kmers <- shapiro.test(x = mwANOVA.acc_kmers.residuals)
#Shapiro-Wilk normality test
#data:  mwANOVA.acc_kmers.residuals
#W = 0.95975, p-value = 0.04573

# perform one-way ANOVA test with no assumption of equal variances
## global
owTST.global.mutation <- oneway.test(average_accuracy_testing ~ mutation, data = df)
#One-way analysis of means (not assuming equal variances)
#data:  average_accuracy_testing and mutation
#F = 14.332, num df = 3.00, denom df = 125.63, p-value = 4.339e-08
owTST.global.splitting <- oneway.test(average_accuracy_testing ~ splitting, data = df)
#One-way analysis of means (not assuming equal variances)
#data:  average_accuracy_testing and splitting
#F = 5.9386, num df = 4.00, denom df = 117.08, p-value = 0.0002188
owTST.global.preprocessing <- oneway.test(average_accuracy_testing ~ preprocessing, data = df)
#One-way analysis of means (not assuming equal variances)
#data:  average_accuracy_testing and preprocessing
#F = 1.6374, num df = 1.00, denom df = 236.54, p-value = 0.2019
owTST.global.model <- oneway.test(average_accuracy_testing ~ model, data = df)
#One-way analysis of means (not assuming equal variances)
#data:  average_accuracy_testing and model
#F = 27.273, num df = 5.00, denom df = 108.24, p-value < 2.2e-16

## core_alleles
owTST.core_alleles.splitting <- oneway.test(average_accuracy_testing ~ splitting, data = df.core_alleles)
#One-way analysis of means (not assuming equal variances)
#data:  average_accuracy_testing and splitting
#F = 1.1815, num df = 4.000, denom df = 27.322, p-value = 0.3409
owTST.core_alleles.preprocessing <- oneway.test(average_accuracy_testing ~ preprocessing, data = df.core_alleles)
#One-way analysis of means (not assuming equal variances)
#data:  average_accuracy_testing and preprocessing
#F = 0.021395, num df = 1.000, denom df = 57.968, p-value = 0.8842
owTST.core_alleles.model <- oneway.test(average_accuracy_testing ~ model, data = df.core_alleles)
#One-way analysis of means (not assuming equal variances)
#data:  average_accuracy_testing and model
#F = 91.978, num df = 5.000, denom df = 25.095, p-value = 2.606e-15

## acc_genes
owTST.acc_genes.splitting <- oneway.test(average_accuracy_testing ~ splitting, data = df.acc_genes)
#One-way analysis of means (not assuming equal variances)
#data:  average_accuracy_testing and splitting
#F = 5.6048, num df = 4.000, denom df = 27.448, p-value = 0.001981
owTST.acc_genes.preprocessing <- oneway.test(average_accuracy_testing ~ preprocessing, data = df.acc_genes)
#One-way analysis of means (not assuming equal variances)
#data:  average_accuracy_testing and preprocessing
#F = 0.00012747, num df = 1.000, denom df = 57.945, p-value = 0.991
owTST.acc_genes.model <- oneway.test(average_accuracy_testing ~ model, data = df.acc_genes)
#One-way analysis of means (not assuming equal variances)
#data:  average_accuracy_testing and model
#F = 12.911, num df = 5.000, denom df = 25.061, p-value = 2.856e-06

## core_SNPs
owTST.core_SNPs.splitting <- oneway.test(average_accuracy_testing ~ splitting, data = df.core_SNPs)
#One-way analysis of means (not assuming equal variances)
#data:  average_accuracy_testing and splitting
#F = 0.94215, num df = 4.000, denom df = 27.407, p-value = 0.4545
owTST.core_SNPs.preprocessing <- oneway.test(average_accuracy_testing ~ preprocessing, data = df.core_SNPs)
#One-way analysis of means (not assuming equal variances)
#data:  average_accuracy_testing and preprocessing
#F = 3.2338, num df = 1.000, denom df = 48.019, p-value = 0.07842
owTST.core_SNPs.model <- oneway.test(average_accuracy_testing ~ model, data = df.core_SNPs)
#One-way analysis of means (not assuming equal variances)
#data:  average_accuracy_testing and model
#F = 42.766, num df = 5.000, denom df = 24.877, p-value = 1.992e-11

## acc_kmers
owTST.acc_kmers.splitting <- oneway.test(average_accuracy_testing ~ splitting, data = df.acc_kmers)
#One-way analysis of means (not assuming equal variances)
#data:  average_accuracy_testing and splitting
#F = 5.5868, num df = 4.00, denom df = 26.99, p-value = 0.002074
owTST.acc_kmers.preprocessing <- oneway.test(average_accuracy_testing ~ preprocessing, data = df.acc_kmers)
#One-way analysis of means (not assuming equal variances)
#data:  average_accuracy_testing and preprocessing
#F = 0.89449, num df = 1.000, denom df = 57.936, p-value = 0.3482
owTST.acc_kmers.model <- oneway.test(average_accuracy_testing ~ model, data = df.acc_kmers)
#One-way analysis of means (not assuming equal variances)
#data:  average_accuracy_testing and model
#F = 15.495, num df = 5.000, denom df = 25.125, p-value = 5.532e-07




# perform Tukey multiple pairwise-comparisons
results.Tukey.mwANOVA.global <- TukeyHSD(mwANOVA.global)
'''
Tukey multiple comparisons of means
95% family-wise confidence level
Fit: aov(formula = average_accuracy_testing ~ mutation + preprocessing + splitting + model, data = df)
$mutation
diff         lwr        upr     p adj
acc_genes-core_alleles  3.160000   0.7073410  5.6126590 0.0054754
core_SNPs-core_alleles -5.736667  -8.1893257 -3.2840076 0.0000000
acc_kmers-core_alleles  1.598333  -0.8543257  4.0509924 0.3330877
core_SNPs-acc_genes    -8.896667 -11.3493257 -6.4440076 0.0000000
acc_kmers-acc_genes    -1.561667  -4.0143257  0.8909924 0.3539373
acc_kmers-core_SNPs     7.335000   4.8823410  9.7876590 0.0000000
$preprocessing
diff       lwr         upr     p adj
near-zero_variances-no_preprocessing -1.3775 -2.697872 -0.05712781 0.0409535
$splitting
diff        lwr        upr     p adj
60-50  3.777083  0.8636329  6.6905338 0.0040217
70-50  5.018750  2.1052996  7.9322004 0.0000375
80-50  7.422917  4.5094662 10.3363671 0.0000000
90-50  2.722917 -0.1905338  5.6363671 0.0794461
70-60  1.241667 -1.6717838  4.1551171 0.7673715
80-60  3.645833  0.7323829  6.5592838 0.0061511
90-60 -1.054167 -3.9676171  1.8592838 0.8574171
80-70  2.404167 -0.5092838  5.3176171 0.1588032
90-70 -2.295833 -5.2092838  0.6176171 0.1961228
90-80 -4.700000 -7.6134504 -1.7865496 0.0001385
$model
diff          lwr         upr     p adj
SVM-XGB   3.3250  -0.01090908   6.6609091 0.0513004
SGB-XGB  -4.1450  -7.48090908  -0.8090909 0.0057293
RF-XGB  -12.3600 -15.69590908  -9.0240909 0.0000000
ERT-XGB  -8.2725 -11.60840908  -4.9365909 0.0000000
BLR-XGB  -4.5300  -7.86590908  -1.1940909 0.0017280
SGB-SVM  -7.4700 -10.80590908  -4.1340909 0.0000000
RF-SVM  -15.6850 -19.02090908 -12.3490909 0.0000000
ERT-SVM -11.5975 -14.93340908  -8.2615909 0.0000000
BLR-SVM  -7.8550 -11.19090908  -4.5190909 0.0000000
RF-SGB   -8.2150 -11.55090908  -4.8790909 0.0000000
ERT-SGB  -4.1275  -7.46340908  -0.7915909 0.0060354
BLR-SGB  -0.3850  -3.72090908   2.9509091 0.9994626
ERT-RF    4.0875   0.75159092   7.4234091 0.0067922
BLR-RF    7.8300   4.49409092  11.1659091 0.0000000
BLR-ERT   3.7425   0.40659092   7.0784091 0.0179245
'''
results.Tukey.mwANOVA.core_alleles <- TukeyHSD(mwANOVA.core_alleles)
'''
Tukey multiple comparisons of means
95% family-wise confidence level
Fit: aov(formula = average_accuracy_testing ~ preprocessing + splitting + model, data = df.core_alleles)
$preprocessing
diff        lwr      upr     p adj
near-zero_variances-no_preprocessing 0.3533333 -0.9308403 1.637507 0.5828296
$splitting
diff        lwr      upr     p adj
60-50  4.2750000  1.4136322 7.136368 0.0009269
70-50  5.8083333  2.9469655 8.669701 0.0000055
80-50  6.9500000  4.0886322 9.811368 0.0000001
90-50  5.3500000  2.4886322 8.211368 0.0000268
70-60  1.5333333 -1.3280345 4.394701 0.5563168
80-60  2.6750000 -0.1863678 5.536368 0.0772499
90-60  1.0750000 -1.7863678 3.936368 0.8238765
80-70  1.1416667 -1.7197011 4.003034 0.7899780
90-70 -0.4583333 -3.3197011 2.403034 0.9909908
90-80 -1.6000000 -4.4613678 1.261368 0.5148173
$model
diff         lwr         upr     p adj
SVM-XGB   3.95   0.6677844   7.2322156 0.0099824
SGB-XGB  -3.04  -6.3222156   0.2422156 0.0841403
RF-XGB  -20.08 -23.3622156 -16.7977844 0.0000000
ERT-XGB -15.88 -19.1622156 -12.5977844 0.0000000
BLR-XGB  -3.79  -7.0722156  -0.5077844 0.0149960
SGB-SVM  -6.99 -10.2722156  -3.7077844 0.0000011
RF-SVM  -24.03 -27.3122156 -20.7477844 0.0000000
ERT-SVM -19.83 -23.1122156 -16.5477844 0.0000000
BLR-SVM  -7.74 -11.0222156  -4.4577844 0.0000001
RF-SGB  -17.04 -20.3222156 -13.7577844 0.0000000
ERT-SGB -12.84 -16.1222156  -9.5577844 0.0000000
BLR-SGB  -0.75  -4.0322156   2.5322156 0.9836245
ERT-RF    4.20   0.9177844   7.4822156 0.0051669
BLR-RF   16.29  13.0077844  19.5722156 0.0000000
BLR-ERT  12.09   8.8077844  15.3722156 0.0000000
'''
results.Tukey.mwANOVA.acc_genes <- TukeyHSD(mwANOVA.acc_genes)
'''
Tukey multiple comparisons of means
95% family-wise confidence level
Fit: aov(formula = average_accuracy_testing ~ preprocessing + splitting + model, data = df.acc_genes)
$preprocessing
                                           diff       lwr      upr     p adj
near-zero_variances-no_preprocessing 0.01333333 -1.154231 1.180897 0.9817842
$splitting
           diff        lwr       upr     p adj
60-50  3.308333  0.7067926  5.909874 0.0063364
70-50  3.033333  0.4317926  5.634874 0.0147793
80-50  7.516667  4.9151259 10.118207 0.0000000
90-50  1.908333 -0.6932074  4.509874 0.2462367
70-60 -0.275000 -2.8765407  2.326541 0.9981897
80-60  4.208333  1.6067926  6.809874 0.0002985
90-60 -1.400000 -4.0015407  1.201541 0.5522585
80-70  4.483333  1.8817926  7.084874 0.0001102
90-70 -1.125000 -3.7265407  1.476541 0.7371903
90-80 -5.608333 -8.2098741 -3.006793 0.0000016
$model
         diff         lwr         upr     p adj
SVM-XGB  0.82  -2.1641734  3.80417336 0.9634211
SGB-XGB -7.23 -10.2141734 -4.24582664 0.0000001
RF-XGB  -2.21  -5.1941734  0.77417336 0.2581568
ERT-XGB  1.62  -1.3641734  4.60417336 0.5961890
BLR-XGB -4.72  -7.7041734 -1.73582664 0.0003061
SGB-SVM -8.05 -11.0341734 -5.06582664 0.0000000
RF-SVM  -3.03  -6.0141734 -0.04582664 0.0446578
ERT-SVM  0.80  -2.1841734  3.78417336 0.9670711
BLR-SVM -5.54  -8.5241734 -2.55582664 0.0000192
RF-SGB   5.02   2.0358266  8.00417336 0.0001129
ERT-SGB  8.85   5.8658266 11.83417336 0.0000000
BLR-SGB  2.51  -0.4741734  5.49417336 0.1458075
ERT-RF   3.83   0.8458266  6.81417336 0.0049959
BLR-RF  -2.51  -5.4941734  0.47417336 0.1458075
BLR-ERT -6.34  -9.3241734 -3.35582664 0.0000012
'''
results.Tukey.mwANOVA.core_SNPs <- TukeyHSD(mwANOVA.core_SNPs)
'''
Tukey multiple comparisons of means
95% family-wise confidence level
Fit: aov(formula = average_accuracy_testing ~ preprocessing + splitting + model, data = df.core_SNPs)
$preprocessing
                                          diff       lwr       upr    p adj
near-zero_variances-no_preprocessing -4.386667 -6.338889 -2.434445 3.98e-05
$splitting
             diff        lwr        upr     p adj
60-50  3.60833333 -0.7415656  7.9582323 0.1470840
70-50  5.91666667  1.5667677 10.2665656 0.0030096
80-50  5.85833333  1.5084344 10.2082323 0.0033754
90-50  1.31666667 -3.0332323  5.6665656 0.9109438
70-60  2.30833333 -2.0415656  6.6582323 0.5656472
80-60  2.25000000 -2.0998990  6.5998990 0.5896837
90-60 -2.29166667 -6.6415656  2.0582323 0.5725133
80-70 -0.05833333 -4.4082323  4.2915656 0.9999995
90-70 -4.60000000 -8.9498990 -0.2501010 0.0333402
90-80 -4.54166667 -8.8915656 -0.1917677 0.0367076
$model
          diff         lwr         upr     p adj
SVM-XGB   4.55  -0.4396788   9.5396788 0.0926607
SGB-XGB  -2.32  -7.3096788   2.6696788 0.7392394
RF-XGB  -18.91 -23.8996788 -13.9203212 0.0000000
ERT-XGB -13.62 -18.6096788  -8.6303212 0.0000000
BLR-XGB  -0.18  -5.1696788   4.8096788 0.9999979
SGB-SVM  -6.87 -11.8596788  -1.8803212 0.0021517
RF-SVM  -23.46 -28.4496788 -18.4703212 0.0000000
ERT-SVM -18.17 -23.1596788 -13.1803212 0.0000000
BLR-SVM  -4.73  -9.7196788   0.2596788 0.0724628
RF-SGB  -16.59 -21.5796788 -11.6003212 0.0000000
ERT-SGB -11.30 -16.2896788  -6.3103212 0.0000003
BLR-SGB   2.14  -2.8496788   7.1296788 0.7986267
ERT-RF    5.29   0.3003212  10.2796788 0.0318332
BLR-RF   18.73  13.7403212  23.7196788 0.0000000
BLR-ERT  13.44   8.4503212  18.4296788 0.0000000
'''
results.Tukey.mwANOVA.acc_kmers <- TukeyHSD(mwANOVA.acc_kmers)
'''
Tukey multiple comparisons of means
95% family-wise confidence level
Fit: aov(formula = average_accuracy_testing ~ preprocessing + splitting + model, data = df.acc_kmers)
$preprocessing
                                      diff       lwr        upr    p adj
near-zero_variances-no_preprocessing -1.49 -2.740938 -0.2390619 0.020557
$splitting
           diff        lwr       upr     p adj
60-50  3.916667  1.1293537  6.703980 0.0020389
70-50  5.316667  2.5293537  8.103980 0.0000185
80-50  9.366667  6.5793537 12.153980 0.0000000
90-50  2.316667 -0.4706463  5.103980 0.1457169
70-60  1.400000 -1.3873130  4.187313 0.6164832
80-60  5.450000  2.6626870  8.237313 0.0000116
90-60 -1.600000 -4.3873130  1.187313 0.4886778
80-70  4.050000  1.2626870  6.837313 0.0013379
90-70 -3.000000 -5.7873130 -0.212687 0.0290822
90-80 -7.050000 -9.8373130 -4.262687 0.0000000
$model
          diff         lwr         upr     p adj
SVM-XGB   3.98   0.7827311   7.1772689 0.0070042
SGB-XGB  -3.99  -7.1872689  -0.7927311 0.0068169
RF-XGB   -8.24 -11.4372689  -5.0427311 0.0000000
ERT-XGB  -5.21  -8.4072689  -2.0127311 0.0001910
BLR-XGB  -9.43 -12.6272689  -6.2327311 0.0000000
SGB-SVM  -7.97 -11.1672689  -4.7727311 0.0000000
RF-SVM  -12.22 -15.4172689  -9.0227311 0.0000000
ERT-SVM  -9.19 -12.3872689  -5.9927311 0.0000000
BLR-SVM -13.41 -16.6072689 -10.2127311 0.0000000
RF-SGB   -4.25  -7.4472689  -1.0527311 0.0033189
ERT-SGB  -1.22  -4.4172689   1.9772689 0.8658342
BLR-SGB  -5.44  -8.6372689  -2.2427311 0.0000931
ERT-RF    3.03  -0.1672689   6.2272689 0.0726003
BLR-RF   -1.19  -4.3872689   2.0072689 0.8774958
BLR-ERT  -4.22  -7.4172689  -1.0227311 0.0036113
'''

# perform multi-way ANOVA test requiring assumptions of normal distribution of residuals and homogeneity of variance in all groups
mwANOVA.global.interaction <- aov(average_accuracy_testing ~ preprocessing * splitting * model, data = df)
results.mwANOVA.global.interaction <- summary(mwANOVA.global.interaction)
#Df Sum Sq Mean Sq F value   Pr(>F)    
#preprocessing                   1    114   113.9   2.542    0.113    
#splitting                       4   1450   362.5   8.096 5.03e-06 ***
#  model                           5   6298  1259.6  28.127  < 2e-16 ***
#  preprocessing:splitting         4     19     4.9   0.109    0.979    
#preprocessing:model             5    224    44.8   0.999    0.420    
#splitting:model                20    457    22.8   0.510    0.960    
#preprocessing:splitting:model  20     40     2.0   0.045    1.000    
#Residuals                     180   8061    44.8
mwANOVA.core_alleles.interaction <- aov(average_accuracy_testing ~ splitting * model, data = df.core_alleles)
results.mwANOVA.core_alleles.interaction <- summary(mwANOVA.core_alleles.interaction)
#Df Sum Sq Mean Sq F value   Pr(>F)    
#splitting        4    345    86.2   67.79 1.37e-14 ***
#  model            5   4432   886.3  696.98  < 2e-16 ***
#  splitting:model 20    264    13.2   10.38 1.40e-08 ***
#  Residuals       30     38     1.3
mwANOVA.acc_genes.interaction <- aov(average_accuracy_testing ~ splitting * model, data = df.acc_genes)
results.mwANOVA.acc_genes.interaction <- summary(mwANOVA.acc_genes.interaction)
#Df Sum Sq Mean Sq F value   Pr(>F)    
#splitting        4  366.8   91.71  20.529 3.09e-08 ***
#  model            5  598.4  119.68  26.790 3.21e-10 ***
#  splitting:model 20  114.1    5.70   1.277    0.266    
#Residuals       30  134.0    4.47 
mwANOVA.core_SNPs.interaction <- aov(average_accuracy_testing ~ splitting * model, data = df.core_SNPs)
results.mwANOVA.core_SNPs.interaction <- summary(mwANOVA.core_SNPs.interaction)
#splitting        4    340    84.9   3.432    0.020 *  
#  model            5   4144   828.7  33.499 2.07e-11 ***
#  splitting:model 20    240    12.0   0.485    0.952    
#Residuals       30    742    24.7    
mwANOVA.acc_kmers.interaction <- aov(average_accuracy_testing ~ splitting * model, data = df.acc_kmers)
results.mwANOVA.acc_kmers.interaction <- summary(mwANOVA.acc_kmers.interaction)
#Df Sum Sq Mean Sq F value   Pr(>F)    
#splitting        4  590.5  147.62  36.723 3.72e-11 ***
#  model            5 1284.0  256.80  63.884 4.41e-15 ***
#  splitting:model 20  197.5    9.88   2.457   0.0127 *  
#  Residuals       30  120.6    4.02  


# check normality assumption with Shapiro-Wilk test
# global
mwANOVA.global.residuals.interaction <- residuals(object = mwANOVA.global.interaction)
Shapiro.mwANOVA.global.interaction <- shapiro.test(x = mwANOVA.global.residuals.interaction)
#Shapiro-Wilk normality test
#data:  mwANOVA.global.residuals.interaction
#W = 0.98679, p-value = 0.0259

# core_alleles
mwANOVA.core_alleles.residuals.interaction <- residuals(object = mwANOVA.core_alleles.interaction)
Shapiro.mwANOVA.core_alleles.interaction <- shapiro.test(x = mwANOVA.core_alleles.residuals.interaction)
#Shapiro-Wilk normality test
#data:  mwANOVA.core_alleles.residuals.interaction
#W = 0.97405, p-value = 0.2291

# acc_genes
mwANOVA.acc_genes.residuals.interaction <- residuals(object = mwANOVA.acc_genes.interaction)
Shapiro.mwANOVA.acc_genes.interaction <- shapiro.test(x = mwANOVA.acc_genes.residuals.interaction)
#Shapiro-Wilk normality test
#data:  mwANOVA.acc_genes.residuals.interaction
#W = 0.97514, p-value = 0.2581

# core_SNPs
mwANOVA.core_SNPs.residuals.interaction <- residuals(object = mwANOVA.core_SNPs.interaction)
Shapiro.mwANOVA.core_SNPs.interaction <- shapiro.test(x = mwANOVA.core_SNPs.residuals.interaction)
#Shapiro-Wilk normality test
#data:  mwANOVA.core_SNPs.residuals.interaction
#W = 0.98601, p-value = 0.7225

# acc_kmers
mwANOVA.acc_kmers.residuals.interaction <- residuals(object = mwANOVA.acc_kmers.interaction)
Shapiro.mwANOVA.acc_kmers.interaction <- shapiro.test(x = mwANOVA.acc_kmers.residuals.interaction)
#Shapiro-Wilk normality test
#data:  mwANOVA.acc_kmers.residuals.interaction
#W = 0.93863, p-value = 0.004707


# perform pairwise t-tests with no assumption of equal variances and with assumption of normality
# global
results.t.test.global.mutation <- pairwise.t.test(df$average_accuracy_testing, df$mutation, p.adjust.method = "BH", pool.sd = FALSE)
#Pairwise comparisons using t tests with non-pooled SD 
#data:  df$average_accuracy_testing and df$mutation 
#core_alleles acc_genes core_SNPs
#acc_genes 0.0300       -         -        
#  core_SNPs 0.0024       3.6e-08   -        
#  acc_kmers 0.2674       0.1371    7.8e-06  
#P value adjustment method: BH
results.t.test.global.splitting <- pairwise.t.test(df$average_accuracy_testing, df$splitting, p.adjust.method = "BH", pool.sd = FALSE)
#Pairwise comparisons using t tests with non-pooled SD 
#data:  df$average_accuracy_testing and df$splitting 
#50      60      70      80     
#60 0.03664 -       -       -      
#  70 0.00402 0.43851 -       -      
#  80 0.00012 0.05150 0.19335 -      
#  90 0.19335 0.55295 0.24394 0.03664
#P value adjustment method: BH
results.t.test.global.preprocessing <- pairwise.t.test(df$average_accuracy_testing, df$preprocessing, p.adjust.method = "BH", pool.sd = FALSE)
#Pairwise comparisons using t tests with non-pooled SD 
#data:  df$average_accuracy_testing and df$preprocessing 
#no_preprocessing
#near-zero_variances 0.2             
#P value adjustment method: BH 
results.t.test.global.model <- pairwise.t.test(df$average_accuracy_testing, df$model, p.adjust.method = "BH", pool.sd = FALSE)
#Pairwise comparisons using t tests with non-pooled SD 
#data:  df$average_accuracy_testing and df$model 
#XGB     SVM     SGB     RF      ERT    
#SVM 0.00347 -       -       -       -      
#  SGB 0.00023 3.1e-09 -       -       -      
#  RF  3.1e-09 6.0e-12 1.6e-05 -       -      
#  ERT 5.9e-06 3.1e-09 0.01385 0.05332 -      
#  BLR 0.00029 8.9e-09 0.73671 5.7e-05 0.03119
#P value adjustment method: BH 

# core_alleles
results.t.test.core_alleles.splitting <- pairwise.t.test(df.core_alleles$average_accuracy_testing, df.core_alleles$splitting, p.adjust.method = "BH", pool.sd = FALSE)
#Pairwise comparisons using t tests with non-pooled SD 
#data:  df.core_alleles$average_accuracy_testing and df.core_alleles$splitting 
#50   60   70   80  
#60 0.56 -    -    -   
#  70 0.46 0.89 -    -   
#  80 0.46 0.89 0.89 -   
#  90 0.56 0.89 0.92 0.89
#P value adjustment method: BH 
results.t.test.core_alleles.preprocessing <- pairwise.t.test(df.core_alleles$average_accuracy_testing, df.core_alleles$preprocessing, p.adjust.method = "BH", pool.sd = FALSE)
#Pairwise comparisons using t tests with non-pooled SD 
#data:  df.core_alleles$average_accuracy_testing and df.core_alleles$preprocessing 
#no_preprocessing
#near-zero_variances 0.88            
#P value adjustment method: BH 
results.t.test.core_alleles.model <- pairwise.t.test(df.core_alleles$average_accuracy_testing, df.core_alleles$model, p.adjust.method = "BH", pool.sd = FALSE)
#Pairwise comparisons using t tests with non-pooled SD 
#data:  df.core_alleles$average_accuracy_testing and df.core_alleles$model 
#XGB     SVM     SGB     RF      ERT    
#SVM 0.01924 -       -       -       -      
#  SGB 0.09765 0.00039 -       -       -      
#  RF  9.0e-10 2.6e-12 1.5e-08 -       -      
#  ERT 2.2e-08 9.6e-11 5.4e-07 0.00723 -      
#  BLR 0.05531 0.00032 0.68456 6.9e-08 2.8e-06
#P value adjustment method: BH 

# acc_genes
results.t.test.acc_genes.splitting <- pairwise.t.test(df.acc_genes$average_accuracy_testing, df.acc_genes$splitting, p.adjust.method = "BH", pool.sd = FALSE)
#Pairwise comparisons using t tests with non-pooled SD 
#data:  df.acc_genes$average_accuracy_testing and df.acc_genes$splitting 
#50     60     70     80    
#60 0.1343 -      -      -     
#  70 0.1343 0.8578 -      -     
#  80 0.0021 0.0312 0.0178 -     
#  90 0.4056 0.4926 0.5188 0.0076
#P value adjustment method: BH 
results.t.test.acc_genes.preprocessing <- pairwise.t.test(df.acc_genes$average_accuracy_testing, df.acc_genes$preprocessing, p.adjust.method = "BH", pool.sd = FALSE)
#Pairwise comparisons using t tests with non-pooled SD 
#data:  df.acc_genes$average_accuracy_testing and df.acc_genes$preprocessing 
#no_preprocessing
#near-zero_variances 0.99            
#P value adjustment method: BH 
results.t.test.acc_genes.model <- pairwise.t.test(df.acc_genes$average_accuracy_testing, df.acc_genes$model, p.adjust.method = "BH", pool.sd = FALSE)
#Pairwise comparisons using t tests with non-pooled SD 
#data:  df.acc_genes$average_accuracy_testing and df.acc_genes$model 
#XGB     SVM     SGB     RF      ERT    
#SVM 0.53688 -       -       -       -      
#  SGB 0.00018 2.1e-05 -       -       -      
#  RF  0.23347 0.08039 0.00934 -       -      
#  ERT 0.30945 0.53688 2.1e-05 0.04387 -      
#  BLR 0.03576 0.01188 0.23347 0.24560 0.00902
#P value adjustment method: BH 

# core_SNPs
results.t.test.core_SNPs.splitting <- pairwise.t.test(df.core_SNPs$average_accuracy_testing, df.core_SNPs$splitting, p.adjust.method = "BH", pool.sd = FALSE)
#Pairwise comparisons using t tests with non-pooled SD 
#data:  df.core_SNPs$average_accuracy_testing and df.core_SNPs$splitting 
#50   60   70   80  
#60 0.62 -    -    -   
#  70 0.62 0.73 -    -   
#  80 0.62 0.73 0.99 -   
#  90 0.83 0.73 0.62 0.62
#P value adjustment method: BH 
results.t.test.core_SNPs.preprocessing <- pairwise.t.test(df.core_SNPs$average_accuracy_testing, df.core_SNPs$preprocessing, p.adjust.method = "BH", pool.sd = FALSE)
#Pairwise comparisons using t tests with non-pooled SD 
#data:  df.core_SNPs$average_accuracy_testing and df.core_SNPs$preprocessing 
#no_preprocessing
#near-zero_variances 0.078           
#P value adjustment method: BH 
results.t.test.core_SNPs.model <- pairwise.t.test(df.core_SNPs$average_accuracy_testing, df.core_SNPs$model, p.adjust.method = "BH", pool.sd = FALSE)
#Pairwise comparisons using t tests with non-pooled SD 
#data:  df.core_SNPs$average_accuracy_testing and df.core_SNPs$model 
#XGB     SVM     SGB     RF      ERT    
#SVM 0.1333  -       -       -       -      
#  SGB 0.3287  0.0279  -       -       -      
#  RF  1.5e-07 1.1e-06 9.6e-07 -       -      
#  ERT 2.6e-06 9.9e-06 3.1e-05 0.0024  -      
#  BLR 0.9413  0.1472  0.4236  2.6e-06 5.1e-05
#P value adjustment method: BH

# acc_kmers
results.t.test.acc_kmers.splitting <- pairwise.t.test(df.acc_kmers$average_accuracy_testing, df.acc_kmers$splitting, p.adjust.method = "BH", pool.sd = FALSE)
#Pairwise comparisons using t tests with non-pooled SD 
#data:  df.acc_kmers$average_accuracy_testing and df.acc_kmers$splitting 
#50     60     70     80    
#60 0.1218 -      -      -     
#  70 0.0332 0.4823 -      -     
#  80 0.0011 0.0382 0.0468 -     
#  90 0.4823 0.5600 0.3418 0.0382
#P value adjustment method: BH 
results.t.test.acc_kmers.preprocessing <- pairwise.t.test(df.acc_kmers$average_accuracy_testing, df.acc_kmers$preprocessing, p.adjust.method = "BH", pool.sd = FALSE)
#Pairwise comparisons using t tests with non-pooled SD 
#data:  df.acc_kmers$average_accuracy_testing and df.acc_kmers$preprocessing 
#no_preprocessing
#near-zero_variances 0.35            
#P value adjustment method: BH 
results.t.test.acc_kmers.model <- pairwise.t.test(df.acc_kmers$average_accuracy_testing, df.acc_kmers$model, p.adjust.method = "BH", pool.sd = FALSE)
#Pairwise comparisons using t tests with non-pooled SD 
#data:  df.acc_kmers$average_accuracy_testing and df.acc_kmers$model 
#XGB     SVM     SGB     RF      ERT    
#SVM 0.04117 -       -       -       -      
#  SGB 0.04590 0.00039 -       -       -      
#  RF  0.00040 5.5e-06 0.03598 -       -      
#  ERT 0.01359 8.3e-05 0.49893 0.08569 -      
#  BLR 0.00159 8.3e-05 0.04117 0.59025 0.08569
#P value adjustment method: BH 



# perform Kruskal-Wallis with no assumptions of equal variances and normality
# global
results.kruskal.global.mutation <- kruskal.test(average_accuracy_testing ~ mutation, data = df)
#Kruskal-Wallis rank sum test
#data:  average_accuracy_testing by mutation
#Kruskal-Wallis chi-squared = 29.833, df = 3, p-value = 1.496e-06
results.kruskal.global.splitting <- kruskal.test(average_accuracy_testing ~ splitting, data = df)
#Kruskal-Wallis rank sum test
#data:  average_accuracy_testing by splitting
#Kruskal-Wallis chi-squared = 26.443, df = 4, p-value = 2.576e-05
results.kruskal.global.preprocessing <- kruskal.test(average_accuracy_testing ~ preprocessing, data = df)
#Kruskal-Wallis rank sum test
#data:  average_accuracy_testing by preprocessing
#Kruskal-Wallis chi-squared = 2.7361, df = 1, p-value = 0.0981
results.kruskal.global.model <- kruskal.test(average_accuracy_testing ~ model, data = df)
#Kruskal-Wallis rank sum test
#data:  average_accuracy_testing by model
#Kruskal-Wallis chi-squared = 88.225, df = 5, p-value < 2.2e-16

# core_alleles
results.kruskal.core_alleles.splitting <- kruskal.test(average_accuracy_testing ~ splitting, data = df.core_alleles)
#Kruskal-Wallis rank sum test
#data:  average_accuracy_testing by splitting
#Kruskal-Wallis chi-squared = 5.8374, df = 4, p-value = 0.2116
results.kruskal.core_alleles.preprocessing <- kruskal.test(average_accuracy_testing ~ preprocessing, data = df.core_alleles)
#Kruskal-Wallis rank sum test
#data:  average_accuracy_testing by preprocessing
#Kruskal-Wallis chi-squared = 0.015802, df = 1, p-value = 0.9
results.kruskal.core_alleles.model <- kruskal.test(average_accuracy_testing ~ model, data = df.core_alleles)
#Kruskal-Wallis rank sum test
#data:  average_accuracy_testing by model
#Kruskal-Wallis chi-squared = 48.514, df = 5, p-value = 2.79e-09

# acc_genes
results.kruskal.acc_genes.splitting <- kruskal.test(average_accuracy_testing ~ splitting, data = df.acc_genes)
#Kruskal-Wallis rank sum test
#data:  average_accuracy_testing by splitting
#Kruskal-Wallis chi-squared = 17.741, df = 4, p-value = 0.001386
results.kruskal.acc_genes.preprocessing <- kruskal.test(average_accuracy_testing ~ preprocessing, data = df.acc_genes)
#Kruskal-Wallis rank sum test
#data:  average_accuracy_testing by preprocessing
#Kruskal-Wallis chi-squared = 0.031487, df = 1, p-value = 0.8592
results.kruskal.acc_genes.model <- kruskal.test(average_accuracy_testing ~ model, data = df.acc_genes)
#Kruskal-Wallis rank sum test
#data:  average_accuracy_testing by model
#Kruskal-Wallis chi-squared = 28.29, df = 5, p-value = 3.195e-05

# core_SNPs
results.kruskal.core_SNPs.splitting <- kruskal.test(average_accuracy_testing ~ splitting, data = df.core_SNPs)
#Kruskal-Wallis rank sum test
#data:  average_accuracy_testing by splitting
#Kruskal-Wallis chi-squared = 3.3017, df = 4, p-value = 0.5087
results.kruskal.core_SNPs.preprocessing <- kruskal.test(average_accuracy_testing ~ preprocessing, data = df.core_SNPs)
#Kruskal-Wallis rank sum test
#data:  average_accuracy_testing by preprocessing
#Kruskal-Wallis chi-squared = 4.0434, df = 1, p-value = 0.04434
results.kruskal.core_SNPs.model <- kruskal.test(average_accuracy_testing ~ model, data = df.core_SNPs)
#Kruskal-Wallis rank sum test
#data:  average_accuracy_testing by model
#Kruskal-Wallis chi-squared = 42.973, df = 5, p-value = 3.742e-08

# acc_kmers
results.kruskal.acc_kmers.splitting <- kruskal.test(average_accuracy_testing ~ splitting, data = df.acc_kmers)
#Kruskal-Wallis rank sum test
#data:  average_accuracy_testing by splitting
#Kruskal-Wallis chi-squared = 15.554, df = 4, p-value = 0.003679
results.kruskal.acc_kmers.preprocessing <- kruskal.test(average_accuracy_testing ~ preprocessing, data = df.acc_kmers)
#Kruskal-Wallis rank sum test
#data:  average_accuracy_testing by preprocessing
#Kruskal-Wallis chi-squared = 0.67348, df = 1, p-value = 0.4118
results.kruskal.acc_kmers.model <- kruskal.test(average_accuracy_testing ~ model, data = df.acc_kmers)
#Kruskal-Wallis rank sum test
#data:  average_accuracy_testing by model
#Kruskal-Wallis chi-squared = 34.568, df = 5, p-value = 1.835e-06




packageVersion("plyr")
#[1] ‘1.8.8’
packageVersion("ggplot2")
#[1] ‘3.4.2’
packageVersion("ggpmisc")
#[1] ‘0.5.2’
packageVersion("reshape2")
#[1] ‘1.4.4’
packageVersion("lubridate")
#[1] ‘1.9.2’
packageVersion("dplyr")
#[1] ‘1.1.2’
packageVersion("car")
#[1] ‘3.1.2’
library(ggpubr)
packageVersion("ggpubr")
#[1] ‘0.6.0’
