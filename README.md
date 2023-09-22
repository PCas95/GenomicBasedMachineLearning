# Introduction to GenomicBasedMachineLearning

The repository "GenomicBasedMachineLearning" aims at making available the 3 scripts used for data management and production of Figures 2 and 3 plus Additional Files 3, 4, 7 and 8 from the research article "Harmonization of supervised machine learning practices for efficient source attribution of Listeria monocytogenes based on genomic data" (<https://doi.org/10.1186/s12864-023-09667-w>).

The 3 scripts available in this repository are:

- [**`heatmappe.R`**](#heatmapper), which was used to create heatmaps (Figure 3 and Additional File 8 of the research article) and correlation plots (Figure 2 and Additional File 7 of the research article) from metrics of machine learning performance;
- [**`boxplots.R`**](#boxplotsr), which was used to produce boxplot‑based distributions of the samples' collection (Additional Files 3 and 4 of the research article);
- [**`ANOVA.R`**](#anovar), used to perform statistical analyses on prediction results (multi-way and one-way Analysis Of Variance - ANOVA, Tukey multiple pairwise-comparisons, pairwise t-tests, Kruskal-Wallis rank sum test, Shapiro-Wilk test, Levene’s test and Wilcoxon test).

**None of the scripts is designed to be automatic**, because they are highly dependent on the dataset and thus on highly variable premises. We encourage step-by-step execution using RStudio.

# Dependencies

The R scripts heatmappe.R, boxplots.R and ANOVA.R were developed and tested with R version 4.2.3 (2023-03-15) and RStudio 2022.12.0+353.

- For heatmappe.R:
    * `library(ggplot2)`
    * `library(plyr)`
    * `library(ggpmisc)`
    * `library(reshape2)`
    * `library(lubridate)`
- For boxplots.R:
    * `library(ggplot2)`
    * `library(dplyr)`
    * `library(ape)`
    * `library(ggprism)`
    * `library(reshape2)`
- For ANOVA.R the same libraries of heatmappe.R will be necessary, plus the following:
    * `library(dplyr)`
    * `library(car)`
    * `library(ggpubr)`

> **Note:** older versions of R and RStudio suffered from unavailability of the `ggpmisc` package. We recommend you use versions equal or higher to R v4.2.3 and RStudio v2022.12.0 to avoid library availability issues.

# Installation (Ubuntu 20.04 / Linux Mint 20.3)

## 1/ Update software repository and update system

```
sudo apt-get update && sudo apt-get upgrade
```
## 2/ Install R

### 2.1/ Check latest R version available

```
sudo apt-cache showpkg r-base
```

### 2.2/ Install latest R version

```
sudo apt-get install r-base
```

### 2.3/ Check installed R version

```
R --version
```

> If manual installation of a specific version is needed or in case of Windows / Mac Operative Systems, source code and installation packages for R latest and older versions can be reached at the official R Archive Network: <https://cran.r-project.org/>

## 3/ Install RStudio

### 3.1/ Check latest RStudio version available

```
sudo apt-cache showpkg rstudio
```

### 3.2/ Install latest RStudio version

```
sudo apt-get install rstudio
```

### 3.3/ Check installed RStudio version

```
rstudio --version
```

> In case manual installation of a specific version is needed or in case of Windows / Mac Operative Systems, Installers and Tarballs of RStudio latest and older versions are available at the official download web page: <https://posit.co/download/rstudio-desktop/>

## 4/ Install R libraries with RStudio

Open RStudio and execute the following lines of code:

```R
install.packages("ggplot2")
install.packages("plyr")
install.packages("ggpmisc")
install.packages("reshape2")
install.packages("lubridate")
```

> Note: The example above only lists libraries related to heatmappe.R. 

# `heatmappe.R`

The script heatmappe.R is a tool for data representation through heatmaps, starting from a dataframe of metrics.

A fully automated version of heatmappe.R may be developed in the future.

The heatmappe.R script is divided into the following 3 sections:
- Preparation:
   - Cleaning environment, path to working directory, import dataframe and library calls;
   - Data formatting (management of decimal digits and percentages, levels reorganization and names);
   - Creation of labels for heatmaps;
- Plot heatmaps;
- Pearson correlation tests & Correlation plots.

With the provided demo dataframe, heatmaps will look like those in the following pictures:

![](./demo/Demo-heatmaps-A-D.png)
![](./demo/Demo-heatmaps-AUC.png)

Likewise, examples of correlation plots from the artificial dataset can be found below:

![](./demo/Demo-correlation-accuracy_training_vs_ROC-AUC.png)
![](./demo/Demo-correlation-accuracy_training_vs_PRG-AUC.png)

For the user's ease, some notes on the script and a general reference to used functions from the `ggplot2` package will be given in the [heatmapper-guide.md](./heatmapper-guide.md) of this repository.

# `boxplots.R`

The `boxplots.R` script is provided for peer-review reproducibility of dataset graphical representations.

The script is divided in 3 sections:
- the first one performs the last dataset filtration steps for samples used in the "Source attribution of Listeria monocytogenes through a versatile supervised machine learning workflow based on genomic data" research article, (*i.e.* dataframe subset based on quality, Phred score, contigs, depth and breadth of coverage);
- the second part executes Kolmogorov-Smirnov tests;
- the third and last part performs the graphical representation.

Below is an example of boxplot produced with the script:

![](./demo/Demo-boxplot.png)

# `ANOVA.R`

The ANOVA.R script is provided for peer-review reproducibility of statistical analyses.

- The first part reads and prepares the same input dataframe used for heatmappe.R. It's used to ensure the dataset is being processed in the same way as the heatmappe.R script.
- The second part prepares to and performs the statistical analyses on machine learning prediction results (mwANOVA and owANOVA, Tukey multiple pairwise-comparisons, pairwise t-tests, Kruskal-Wallis rank sum test, Shapiro-Wilk test, Levene’s test and Wilcoxon test).

# References

- Castelli P., De Ruvo A., Bucciacchio A., D’Alterio N., Cammà C., Di Pasquale A. and Radomski N. "Harmonization of supervised machine learning practices for efficient source attribution of Listeria monocytogenes based on genomic data". BMC Genomics 24, 560 (2023), <https://doi.org/10.1186/s12864-023-09667-w>.

# Acknowledgements

My colleague and mentor Nicolas Radomski for directions, guide and collaboration in writing the three R scripts provided in this repository.

# Authors

Pierluigi Castelli and Nicolas Radomski 
