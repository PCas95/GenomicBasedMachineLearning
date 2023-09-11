https://lukaswallrich.github.io/GoldCoreQuants/chi-squared-tests-for-associations-between-categorical-variables.html

#clean environment
rm(list=ls())

# install regular packages
install.packages("dplyr")

# call library
library(dplyr)

# set working directory
setwd("/home/IZSNT/n.radomski/Downloads")

# read dataframe
data_long = read.table("sources-CCs.tsv", dec = ".", header=TRUE, sep = "\t", quote = "")
str(data_long)

# check table
array <- table(data_long$CC, data_long$Phenotypes)

dairy fruit leafy_greens meat poultry seafood vegetables
CC1       10    11            3    8       3       0          2
CC101      8     0            0    0       0       4          0
CC1056     0     1            0    0       0       0          0
CC1070     0     3            0    0       0       0          0
CC11       3     0            1    1       1       0         16
CC1127     0     0            1    0       0       0          0
CC121      0     0            0    4      11      13          2
CC124      0     0            0    0       0       0         12

# association test
chisq.test(data_long$CC, data_long$Phenotypes, correct = FALSE)

Pearsons Chi-squared test
data:  data_long$CC and data_long$Phenotypes
X-squared = 2161.4, df = 546, p-value < 2.2e-16

# association test for small number
## check expectations < 7
chisq.test(data_long$CC, data_long$Phenotypes, correct = FALSE)$expected

# Pearson's Chi-squared test with Yates' continuity correction (based on 2000 Monte Carlo simulations)
chisq.test(data_long$CC, data_long$Phenotypes, simulate.p.value = TRUE)

Pearsons Chi-squared test with simulated p-value (based on 2000 replicates)
data:  data_long$CC and data_long$Phenotypes
X-squared = 2161.4, df = NA, p-value = 0.0004998

=> poor association between Sources and CCs

# Post-hoc tests
## check prop array
prop.array <- prop.table(array, margin = 1) %>% #Share of each cell as part of the row
  round(2) #Round to 2 decimal places

dairy fruit leafy_greens meat poultry seafood vegetables
CC1     0.27  0.30         0.08 0.22    0.08    0.00       0.05
CC101   0.67  0.00         0.00 0.00    0.00    0.33       0.00
CC1056  0.00  1.00         0.00 0.00    0.00    0.00       0.00
CC1070  0.00  1.00         0.00 0.00    0.00    0.00       0.00
CC11    0.14  0.00         0.05 0.05    0.05    0.00       0.73
CC1127  0.00  0.00         1.00 0.00    0.00    0.00       0.00
CC121   0.00  0.00         0.00 0.13    0.37    0.43       0.07
CC124   0.00  0.00         0.00 0.00    0.00    0.00       1.00

# compare each cell to the expected value
sudo apt-get install libfreetype6-dev
install.packages("devtools")
devtools::install_github("ebbertd/chisq.posthoc.test", force = TRUE)
install.packages("pacman")
pacman::p_load(chisq.posthoc.test)
getNamespaceVersion("pacman")
getNamespaceVersion("stats")

chisq.posthoc.test.results <- table(data_long$CC, data_long$Phenotypes) %>% chisq.posthoc.test(simulate.p.value = TRUE, method = "bonferroni")

Dimension     Value               dairy              fruit       leafy_greens               meat            poultry            seafood          vegetables
1         CC1 Residuals     1.6332913100691  0.960124643703037 -0.116365027016053   2.45450619855182 -0.372977597116192  -2.38227553098023   -2.08142568369914
2         CC1  p values                   1                  1                  1                  1                  1                  1                   1
3       CC101 Residuals    4.58731136239676  -1.91343737905119  -1.07090584681755  -1.14928717285655   -1.1551780807182   2.12170125371394   -1.65700689417222
4       CC101  p values             0.0029*                  1                  1                  1                  1                  1                   1
5      CC1056 Residuals  -0.454233333998195   1.82119211205805 -0.307592872712426 -0.330106091138668  -0.33179811782564 -0.385175729722786  -0.475936808261314
6      CC1056  p values                   1                  1                  1                  1                  1                  1                   1
7      CC1070 Residuals  -0.787472074366889   3.15727143510051 -0.533251920116208  -0.57228148815504 -0.575214834665181 -0.667751809870008  -0.825097846453035
8      CC1070  p values                   1                  1                  1                  1                  1                  1                   1
9        CC11 Residuals  -0.434816038914397  -2.60279544666528 -0.690026838383968 -0.839573229541886 -0.850550762000105  -1.82414652507077    6.62870044173472
10       CC11  p values                   1                  1                  1                  1                  1                  1                  0*

# export
## array
write.table(array, 
            file = "array.tsv", 
            append = FALSE, 
            quote = TRUE, 
            sep = "\t", 
            dec = ".", 
            row.names = TRUE, 
            col.names = NA)
## prop.array
write.table(prop.array, 
            file = "prop.array.tsv", 
            append = FALSE, 
            quote = TRUE, 
            sep = "\t", 
            dec = ".", 
            row.names = TRUE, 
            col.names = NA)
## chisq.posthoc.test.results
write.table(chisq.posthoc.test.results, 
            file = "chisq.posthoc.test.results.tsv", 
            append = FALSE, 
            quote = TRUE, 
            sep = "\t", 
            dec = ".", 
            row.names = FALSE, 
            col.names = TRUE)
