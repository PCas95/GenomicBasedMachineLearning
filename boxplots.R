#### boxplots.R ####
#### R script to produce graphical representations of ML for for the research article ####
#### "Source attribution of Listeria monocytogenes through a versatile supervised machine learning workflow based on genomic data" ####
#### performed with R version 4.2.3 (2023-03-15) and RStudio 2022.12.0+353 on Ubuntu 20.04 5 LTS "Focal Fossa" ####
#### input: dataframe_complete.tsv from Bash the workflow used to filter samples of ML for FSA ####


# clean environment
rm(list=ls())

# set working directory
setwd("/home/IZSNT/p.castelli/Documenti/GitHub_repos/heatmappeR/demo")

# install regular packages
install.packages("ggplot2")
install.packages("dplyr")
install.packages("ape")
install.packages("ggprism")
install.packages("reshape2")

# call libraries
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(ggprism))
suppressPackageStartupMessages(library(reshape2))

# read dataframe
data_long <- read.table("dataframe_complete.tsv", dec = ".", header=TRUE, sep = "\t", quote = "")
# check dataframe
head(data_long, 20)
str(data_long)
# check dimension
dim(data_long)


# perform subset to filter on contamination + bad quality and assembly values
## subset 0: filter for contamination (number of contaminated SNVs > 10)
sbs0 <- subset(data_long, data_long$NumContamSNVs < 10)
### check dataframe
head(sbs0, 20)
dim(sbs0)
## subset 1: filter for Q30 (phred score higher than 30) lower than 90%
sbs1 <- subset(sbs0, sbs0$higher_than_Q30 >= 90.00)
### check dataframe
head(sbs1, 20)
dim(sbs1)
## subset 2: filter for depth of coverage values below 35.00%
sbs2 <- subset(sbs1, sbs1$Depth_of_coverage >= 35.000)
### check dataframe
head(sbs2, 20)
dim(sbs2)
## subset 3: filter for breadth of coverage values below 85.00%
sbs3 <- subset(sbs2, sbs2$Breadth_of_coverage >= 85.000)
### check dataframe
head(sbs3, 20)
dim(sbs3)
## subset 4: filter for samples with a number of contigs higher than 400
sbs4 <- subset(sbs3, sbs3$Num_contigs <= 400)
### check dataframe
head(sbs4, 20)
dim(sbs4)
## subset 5: filter for samples with total length of contigs greater than 3.5Mbp
sbs5 <- subset(sbs4, sbs4$Total_length < 3600000)
### check dataframe
head(sbs5, 20)
dim(sbs5)

### save subset dataframe
write.table(sbs5, file = "dataframe_complete-filtered.tsv", quote=FALSE, sep='\t', col.names = NA)


### phenotypes table ###

# create table
ph = table(sbs5$Phenotypes)
# write table as a tsv file with "phenotype" and "samples" as column names
cln <- c("phenotypes", "samples")
write.table(ph, file = "samples_x_phenotypes.tsv", quote=FALSE, sep='\t', col.names = cln, row.names = FALSE)


### boxplot for breadth of coverage ###

# create a short dataframe for breadth of coverage from long dataframe
# (sample vs phenotype, table with breadth values)
short_brcov <- dcast(sbs5, formula = GenPat_ID~Phenotypes, value.var = "Breadth_of_coverage")
# check dataframe
head(short_brcov, 20)
str(short_brcov)
# check dimension
dim(short_brcov)


# Kolmogorov Smirnov test for p-value ks(d$x, d$y)
## p-value dairy-fruit -- Kolmogorov
p_breadth_dairyfruit <- ks.test(short_brcov$dairy, short_brcov$fruit)
v_breadth_dairyfruit <- format(p_breadth_dairyfruit$p.value, digits=7)
cv_breadth_dairyfruit <- as.character(v_breadth_dairyfruit)

## p-value dairy-leafy_greens -- Kolmogorov
p_breadth_dairyleafy <- ks.test(short_brcov$dairy, short_brcov$leafy_greens)
v_breadth_dairyleafy <- format(p_breadth_dairyleafy$p.value, digits=7)
cv_breadth_dairyleafy <- as.character(v_breadth_dairyleafy)

## p-value dairy-meat -- Kolmogorov
p_breadth_dairymeat <- ks.test(short_brcov$dairy, short_brcov$meat)
v_breadth_dairymeat <- format(p_breadth_dairymeat$p.value, digits=7)
cv_breadth_dairymeat <- as.character(v_breadth_dairymeat)

## p-value dairy-poultry -- Kolmogorov
p_breadth_dairypoultry <- ks.test(short_brcov$dairy, short_brcov$poultry)
v_breadth_dairypoultry <- format(p_breadth_dairypoultry$p.value, digits=7)
cv_breadth_dairypoultry <- as.character(v_breadth_dairypoultry)

## p-value dairy-seafood -- Kolmogorov
p_breadth_dairyseafood <- ks.test(short_brcov$dairy, short_brcov$seafood)
v_breadth_dairyseafood <- format(p_breadth_dairyseafood$p.value, digits=7)
cv_breadth_dairyseafood <- as.character(v_breadth_dairyseafood)

## p-value dairy-vegetables -- Kolmogorov
p_breadth_dairyvegeta <- ks.test(short_brcov$dairy, short_brcov$vegetables)
v_breadth_dairyvegeta <- format(p_breadth_dairyvegeta$p.value, digits=7)
cv_breadth_dairyvegeta <- as.character(v_breadth_dairyvegeta)

# create dataframes for 'add_pvalue()' argument in ggplot
## dairy-fruit
pvalue1 <- tibble::tribble(
    ~group1, ~group2, ~p,     ~y.position,
    "dairy",    "fruit",    cv_breadth_dairyfruit, 110
)
## dairy-leafy_greens
pvalue2 <- tibble::tribble(
  ~group1, ~group2, ~p,     ~y.position,
  "dairy",    "leafy_greens",    cv_breadth_dairyleafy, 108
)
## dairy-meat
pvalue3 <- tibble::tribble(
    ~group1, ~group2, ~p,     ~y.position,
    "dairy",    "meat",    cv_breadth_dairymeat, 106
)
## dairy-poultry
pvalue4 <- tibble::tribble(
    ~group1, ~group2, ~p,     ~y.position,
    "dairy",    "poultry",    cv_breadth_dairypoultry, 104
)
## dairy-seafood
pvalue5 <- tibble::tribble(
    ~group1, ~group2, ~p,     ~y.position,
    "dairy",    "seafood",    cv_breadth_dairyseafood, 102
)
## dairy-vegetables
pvalue6 <- tibble::tribble(
    ~group1, ~group2, ~p,     ~y.position,
    "dairy",    "vegetables",    cv_breadth_dairyvegeta, 100
)

# graphical representation of breadth of coverage
b <- seq(85, 100, 5)
p <- ggplot(data = sbs5, aes(x = Phenotypes, y = Breadth_of_coverage)) +
  theme_light(base_size = 16) +
  geom_boxplot(color = "#000000", fill = "#A9A9A9", coef = 6, outlier.colour = "white", outlier.shape = 0, outlier.size = 0) +
  geom_point(position = position_jitter(width = 0.35), size = 1, color = "#000000", alpha = 0.7, shape = "x") +
  theme(axis.text.x = element_text (color = "#000000", size = 12, angle = 0, vjust = 0.5)) +
  scale_y_continuous(name = "Breadth of coverage (%)", limits = c(85.0, 110.0), breaks = b) +
  scale_x_discrete(name = "Phenotypes") +
  add_pvalue(pvalue1) +
  add_pvalue(pvalue2) +
  add_pvalue(pvalue3) +
  add_pvalue(pvalue4) +
  add_pvalue(pvalue5) +
  add_pvalue(pvalue6) +
  theme(plot.title = element_text(color="black", size=12, face="bold.italic"),
        strip.text.x = element_text(size=16, face = "bold"),
        strip.text.y = element_text(size=16, face="bold"),
        strip.background = element_rect(colour="black", fill="#A9A9A9"))
p
plot(p)
ggsave("BofC-ML4FSA.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("BofC-ML4FSA.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")


### boxplot for depth of coverage ###

# create a short dataframe for depth of coverage from long dataframe
# (sample vs phenotype, table with depth values)
short_dpcov <- dcast(sbs5, formula = GenPat_ID~Phenotypes, value.var = "Depth_of_coverage")
# check dataframe
head(short_dpcov, 20)
str(short_dpcov)
# check dimension
dim(short_dpcov)

# Kolmogorov Smirnov test for p-value ks(d$x, d$y)
## p-value dairy-fruit -- Kolmogorov
p_depth_dairyfruit <- ks.test(short_dpcov$dairy, short_dpcov$fruit)
v_depth_dairyfruit <- format(p_depth_dairyfruit$p.value, digits=7)
cv_depth_dairyfruit <- as.character(v_depth_dairyfruit)

## p-value dairy-leafy_greens -- Kolmogorov
p_depth_dairyleafy <- ks.test(short_dpcov$dairy, short_dpcov$leafy_greens)
v_depth_dairyleafy <- format(p_depth_dairyleafy$p.value, digits=7)
cv_depth_dairyleafy <- as.character(v_depth_dairyleafy)

## p-value dairy-meat -- Kolmogorov
p_depth_dairymeat <- ks.test(short_dpcov$dairy, short_dpcov$meat)
v_depth_dairymeat <- format(p_depth_dairymeat$p.value, digits=7)
cv_depth_dairymeat <- as.character(v_depth_dairymeat)

## p-value dairy-poultry -- Kolmogorov
p_depth_dairypoultry <- ks.test(short_dpcov$dairy, short_dpcov$poultry)
v_depth_dairypoultry <- format(p_depth_dairypoultry$p.value, digits=7)
cv_depth_dairypoultry <- as.character(v_depth_dairypoultry)

## p-value dairy-seafood -- Kolmogorov
p_depth_dairyseafood <- ks.test(short_dpcov$dairy, short_dpcov$seafood)
v_depth_dairyseafood <- format(p_depth_dairyseafood$p.value, digits=7)
cv_depth_dairyseafood <- as.character(v_depth_dairyseafood)

## p-value dairy-vegetables -- Kolmogorov
p_depth_dairyvegeta <- ks.test(short_dpcov$dairy, short_dpcov$vegetables)
v_depth_dairyvegeta <- format(p_depth_dairyvegeta$p.value, digits=7)
cv_depth_dairyvegeta <- as.character(v_depth_dairyvegeta)

# dataframes for 'add_pvalue()' argument in ggplot
## dairy-fruit
pvalue7 <- tibble::tribble(
    ~group1, ~group2, ~p,     ~y.position,
    "dairy",    "fruit",    cv_depth_dairyfruit, 95
)
## dairy-leafy_greens
pvalue8 <- tibble::tribble(
    ~group1, ~group2, ~p,     ~y.position,
    "dairy",    "leafy_greens",    cv_depth_dairyleafy, 90
)
## dairy-meat
pvalue9 <- tibble::tribble(
    ~group1, ~group2, ~p,     ~y.position,
    "dairy",    "meat",    cv_depth_dairymeat, 85
)
## dairy-poultry
pvalue10 <- tibble::tribble(
    ~group1, ~group2, ~p,     ~y.position,
    "dairy",    "poultry",    cv_depth_dairypoultry, 80
)
## dairy-seafood
pvalue11 <- tibble::tribble(
    ~group1, ~group2, ~p,     ~y.position,
    "dairy",    "seafood",    cv_depth_dairyseafood, 75
)
## dairy-vegetables
pvalue12 <- tibble::tribble(
    ~group1, ~group2, ~p,     ~y.position,
    "dairy",    "vegetables",    cv_depth_dairyvegeta, 70
)


# graphical representation of depth of coverage
b2 <- seq(30, 100, 10)
p <- ggplot(data = sbs5, aes(x = Phenotypes, y = Depth_of_coverage)) +
    theme_light(base_size = 16) +
    geom_boxplot(color = "#000000", fill = "#A9A9A9", coef = 6, outlier.colour = "white", outlier.shape = 0, outlier.size = 0) +
    geom_point(position = position_jitter(width = 0.35), size = 1, color = "#000000", alpha = 0.7, shape = "x") +
    theme(axis.text.x = element_text (color = "#000000", size = 12, angle = 0, vjust = 0.5)) +
    scale_y_continuous(name = "Depth of coverage (%)", limits = c(30.0, 100.0), breaks = b2) +
    scale_x_discrete(name = "Phenotypes") +
    add_pvalue(pvalue7) +
    add_pvalue(pvalue8) +
    add_pvalue(pvalue9) +
    add_pvalue(pvalue10) +
    add_pvalue(pvalue11) +
    add_pvalue(pvalue12) +
    theme(plot.title = element_text(color="black", size=12, face="bold.italic"),
          strip.text.x = element_text(size=16, face = "bold"),
          strip.text.y = element_text(size=16, face="bold"),
          strip.background = element_rect(colour="black", fill="#A9A9A9"))
p
plot(p)
ggsave("DofC-ML4FSA.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("DofC-ML4FSA.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")



### boxplot for total length of contigs ###

# create a short dataframe for total length of contigs from long dataframe
# (sample vs phenotype, table with total length of contigs)
short_totl <- dcast(sbs5, formula = GenPat_ID~Phenotypes, value.var = "Total_length")
# check dataframe
head(short_totl, 20)
str(short_totl)
# check dimension
dim(short_totl)

# Kolmogorov Smirnov test for p-value ks(d$x, d$y)
## p-value dairy-fruit -- Kolmogorov
p_totl_dairyfruit <- ks.test(short_totl$dairy, short_totl$fruit)
v_totl_dairyfruit <- format(p_totl_dairyfruit$p.value, digits=7)
cv_totl_dairyfruit <- as.character(v_totl_dairyfruit)

## p-value dairy-leafy_greens -- Kolmogorov
p_totl_dairyleafy <- ks.test(short_totl$dairy, short_totl$leafy_greens)
v_totl_dairyleafy <- format(p_totl_dairyleafy$p.value, digits=7)
cv_totl_dairyleafy <- as.character(v_totl_dairyleafy)

## p-value dairy-meat -- Kolmogorov
p_totl_dairymeat <- ks.test(short_totl$dairy, short_totl$meat)
v_totl_dairymeat <- format(p_totl_dairymeat$p.value, digits=7)
cv_totl_dairymeat <- as.character(v_totl_dairymeat)

## p-value dairy-poultry -- Kolmogorov
p_totl_dairypoultry <- ks.test(short_totl$dairy, short_totl$poultry)
v_totl_dairypoultry <- format(p_totl_dairypoultry$p.value, digits=7)
cv_totl_dairypoultry <- as.character(v_totl_dairypoultry)

## p-value dairy-seafood -- Kolmogorov
p_totl_dairyseafood <- ks.test(short_totl$dairy, short_totl$seafood)
v_totl_dairyseafood <- format(p_totl_dairyseafood$p.value, digits=7)
cv_totl_dairyseafood <- as.character(v_totl_dairyseafood)

## p-value dairy-vegetables -- Kolmogorov
p_totl_dairyvegeta <- ks.test(short_totl$dairy, short_totl$vegetables)
v_totl_dairyvegeta <- format(p_totl_dairyvegeta$p.value, digits=7)
cv_totl_dairyvegeta <- as.character(v_totl_dairyvegeta)

# dataframes for 'add_pvalue()' argument in ggplot
## dairy-fruit
pvalue13 <- tibble::tribble(
    ~group1, ~group2, ~p,     ~y.position,
    "dairy",    "fruit",    cv_totl_dairyfruit, 3500000
)
## dairy-leafy_greens
pvalue14 <- tibble::tribble(
    ~group1, ~group2, ~p,     ~y.position,
    "dairy",    "leafy_greens",    cv_totl_dairyleafy, 3450000
)
## dairy-meat
pvalue15 <- tibble::tribble(
    ~group1, ~group2, ~p,     ~y.position,
    "dairy",    "meat",    cv_totl_dairymeat, 3400000
)
## dairy-poultry
pvalue16 <- tibble::tribble(
    ~group1, ~group2, ~p,     ~y.position,
    "dairy",    "poultry",    cv_totl_dairypoultry, 3350000
)
## dairy-seafood
pvalue17 <- tibble::tribble(
    ~group1, ~group2, ~p,     ~y.position,
    "dairy",    "seafood",    cv_totl_dairyseafood, 3300000
)
## dairy-vegetables
pvalue18 <- tibble::tribble(
    ~group1, ~group2, ~p,     ~y.position,
    "dairy",    "vegetables",    cv_totl_dairyvegeta, 3250000
)

# graphical representation of total length of contigs
b3 <- seq(2700000, 3300000, 100000)
p <- ggplot(data = sbs5, aes(x = Phenotypes, y = Total_length)) +
    theme_light(base_size = 16) +
    geom_boxplot(color = "#000000", fill = "#A9A9A9", coef = 6, outlier.colour = "white", outlier.shape = 0, outlier.size = 0) +
    geom_point(position = position_jitter(width = 0.35), size = 1, color = "#000000", alpha = 0.7, shape = "x") +
    theme(axis.text.x = element_text (color = "#000000", size = 12, angle = 0, vjust = 0.5)) +
    scale_y_continuous(name = "Total length (bp)", limits = c(2700000.0, 3500000.0), breaks = b3) +
    scale_x_discrete(name = "Phenotypes") +
    add_pvalue(pvalue13) +
    add_pvalue(pvalue14) +
    add_pvalue(pvalue15) +
    add_pvalue(pvalue16) +
    add_pvalue(pvalue17) +
    add_pvalue(pvalue18) +
    theme(plot.title = element_text(color="black", size=12, face="bold.italic"),
          strip.text.x = element_text(size=16, face = "bold"),
          strip.text.y = element_text(size=16, face="bold"),
          strip.background = element_rect(colour="black", fill="#A9A9A9"))
p
plot(p)
ggsave("totl-ML4FSA.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("totl-ML4FSA.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")


### boxplot for total number of contigs ###

# create a short dataframe for total number of contigs from long dataframe
# (sample vs phenotype, table with total number of contigs)
short_totn <- dcast(sbs5, formula = GenPat_ID~Phenotypes, value.var = "Num_contigs")
# check dataframe
head(short_totn, 20)
str(short_totn)
# check dimension
dim(short_totn)

# Kolmogorov Smirnov test for p-value ks(d$x, d$y)
## p-value dairy-fruit -- Kolmogorov
p_totn_dairyfruit <- ks.test(short_totn$dairy, short_totn$fruit)
v_totn_dairyfruit <- format(p_totn_dairyfruit$p.value, digits=7)
cv_totn_dairyfruit <- as.character(v_totn_dairyfruit)

## p-value dairy-leafy_greens -- Kolmogorov
p_totn_dairyleafy <- ks.test(short_totn$dairy, short_totn$leafy_greens)
v_totn_dairyleafy <- format(p_totn_dairyleafy$p.value, digits=7)
cv_totn_dairyleafy <- as.character(v_totn_dairyleafy)

## p-value dairy-meat -- Kolmogorov
p_totn_dairymeat <- ks.test(short_totn$dairy, short_totn$meat)
v_totn_dairymeat <- format(p_totn_dairymeat$p.value, digits=7)
cv_totn_dairymeat <- as.character(v_totn_dairymeat)

## p-value dairy-poultry -- Kolmogorov
p_totn_dairypoultry <- ks.test(short_totn$dairy, short_totn$poultry)
v_totn_dairypoultry <- format(p_totn_dairypoultry$p.value, digits=7)
cv_totn_dairypoultry <- as.character(v_totn_dairypoultry)

## p-value dairy-seafood -- Kolmogorov
p_totn_dairyseafood <- ks.test(short_totn$dairy, short_totn$seafood)
v_totn_dairyseafood <- format(p_totn_dairyseafood$p.value, digits=7)
cv_totn_dairyseafood <- as.character(v_totn_dairyseafood)

## p-value dairy-vegetables -- Kolmogorov
p_totn_dairyvegeta <- ks.test(short_totn$dairy, short_totn$vegetables)
v_totn_dairyvegeta <- format(p_totn_dairyvegeta$p.value, digits=7)
cv_totn_dairyvegeta <- as.character(v_totn_dairyvegeta)

# dataframes for 'add_pvalue()' argument in ggplot
## dairy-fruit
pvalue19 <- tibble::tribble(
    ~group1, ~group2, ~p,     ~y.position,
    "dairy",    "fruit",    cv_totn_dairyfruit, 630
)
## dairy-leafy_greens
pvalue20 <- tibble::tribble(
    ~group1, ~group2, ~p,     ~y.position,
    "dairy",    "leafy_greens",    cv_totn_dairyleafy, 590
)
## dairy-meat
pvalue21 <- tibble::tribble(
    ~group1, ~group2, ~p,     ~y.position,
    "dairy",    "meat",    cv_totn_dairymeat, 550
)
## dairy-poultry
pvalue22 <- tibble::tribble(
    ~group1, ~group2, ~p,     ~y.position,
    "dairy",    "poultry",    cv_totn_dairypoultry, 510
)
## dairy-seafood
pvalue23 <- tibble::tribble(
    ~group1, ~group2, ~p,     ~y.position,
    "dairy",    "seafood",    cv_totn_dairyseafood, 470
)
## dairy-vegetables
pvalue24 <- tibble::tribble(
    ~group1, ~group2, ~p,     ~y.position,
    "dairy",    "vegetables",    cv_totn_dairyvegeta, 430
)

# graphical representation of total length of contigs
b4 <- seq(0, 500, 100)
p <- ggplot(data = sbs5, aes(x = Phenotypes, y = Num_contigs)) +
    theme_light(base_size = 16) +
    geom_boxplot(color = "#000000", fill = "#A9A9A9", coef = 6, outlier.colour = "white", outlier.shape = 0, outlier.size = 0) +
    geom_point(position = position_jitter(width = 0.35), size = 1, color = "#000000", alpha = 0.7, shape = "x") +
    theme(axis.text.x = element_text (color = "#000000", size = 12, angle = 0, vjust = 0.5)) +
    scale_y_continuous(name = "Contigs", limits = c(0.0, 650.0), breaks = b4) +
    scale_x_discrete(name = "Phenotypes") +
    add_pvalue(pvalue19) +
    add_pvalue(pvalue20) +
    add_pvalue(pvalue21) +
    add_pvalue(pvalue22) +
    add_pvalue(pvalue23) +
    add_pvalue(pvalue24) +
    theme(plot.title = element_text(color="black", size=12, face="bold.italic"),
          strip.text.x = element_text(size=16, face = "bold"),
          strip.text.y = element_text(size=16, face="bold"),
          strip.background = element_rect(colour="black", fill="#A9A9A9"))
p
plot(p)
ggsave("totn-ML4FSA.tiff",device="tiff",width=17,height=17,units="cm",dpi="retina")
ggsave("totn-ML4FSA.pdf",device="pdf",width=17,height=17,units="cm",dpi="retina")


### table "class phenotype vs Clonal Complex" ###

# create a contingency table of Phenotypes_x_CC
ccxph <- table(sbs5$CC, sbs5$Phenotypes)
# save contingency table as tsv file
write.table(ccxph, file = "phenotypes_x_CC.tsv", quote=FALSE, sep='\t', col.names = NA)

