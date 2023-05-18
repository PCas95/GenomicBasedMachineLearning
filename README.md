# Introduction to GenomicBasedMachineLearning

The repository "GenomicBasedMachineLearning" aims at making available the 3 scripts used for data management and production of Figures 2 and 3 plus Additional Files 3, 4, 7 and 8 from the research article "Source attribution of Listeria monocytogenes through a versatile supervised machine learning workflow based on genomic data".

The 3 scripts available in this repository are:

- [**`heatmappe.R`**](#heatmapper), which was used to create heatmaps (Figure 3 and Additional File 8 of the research article) and correlation plots (Figure 2 and Additional File 7 of the research article) from metrics of machine learning performance;
- [**`boxplots.R`**](#boxplotsr), which was used to produce boxplot‑based distributions of the samples' collection (Additional Files 3 and 4 of the research article);
- [**`ANOVA.R`**](#anovar), used to perform statistical analyses on prediction results (multi-way and one-way Analysis Of Variance - ANOVA, Tukey multiple pairwise-comparisons, pairwise t-tests, Kruskal-Wallis rank sum test, Shapiro-Wilk test, Levene’s test and Wilcoxon test).

In this 'readme' will be also fully explained how to use `heatmappe.R`, while `boxplots.R` and `ANOVA.R` will be made available to grant reproducibility, but their execution will not be described, because they highly depend on the dataset and thus on highly variable premises. 

**For the reasons explained above, none of the scripts is designed to be automatic.**

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

> The libraries above only refer to heatmappe.R. Use the same function on libraries used by the other 2 scripts as needed.

The needed libraries have now been installed on your system. The `install.packages()` R function calls are also present as commented lines in the heatmappe.R, boxplots.R and ANOVA.R source code.

# `heatmappe.R`

The GenomicBasedMachineLearning repository provides access to the R script heatmappe.R, a tool to build heatmaps from a dataframe of metrics.

The heatmappe.R has been developed and updated as a tool for data representation for the research article "Source attribution of Listeria monocytogenes through a versatile supervised machine learning workflow based on genomic data".

Since heatmappe.R was developed for research purpuses, with no need or desire for automation, and since the user's research will most certainly have different data and goals, the tool is not designed to be automatic, just as the other 2 scripts. We encourage step-by-step execution using RStudio.

Below you can find a full guide to usage with RStudio. A fully automated version of heatmappe.R may be developed in the future.

## Usage with RStudio

Since heatmappe.R is not designed to be automatic we encourage manual execution using RStudio, in order to tailor all the function to the used dataset.

Open the heatmappe.R Rscript file with RStudio, then proceed to the following steps. Each part of the heatmappe.R script is commented to allow easier identification of each step of the guide.

## Preparation

The preparation section includes cleaning the environment from variables, objects and graphs (**lines 8 and 11**), setting the path to working directory (**line 14**, all input files should reside in the same directory; all output files will be saved in there) and library calls.

The necessary libraries should already be installed on the system ([step 4 of the "Install RStudio" section](#4-install-r-libraries-with-rstudio)), so we just need to call them (*e.g.* `library(ggplot2)`, **lines 24-28**). If you didn't install the libraries before, you can do that now, *e.g.* with `install.packages("ggplot2")` (**lines 17-21**).

## Read dataframe

Read your dataframe with the `read.table` function (**line 32**). Name, decimal separator and column separator should be edited as needed.

> **Our example dataframe has 15 columns. Data extraction  froma the dataframe is managed through header name of each column, rather than number. This also means that the reference example names hardcoded in the script should be changed to match those of the dataframe that is being used.**

**Lines 34, 36, 39 and 40** are optional checks to let the user inspect the dataframe, to ensure it has been imported correctly.

## Data formatting

### Decimal numbers

**Lines 44 to 54** are used to ensure data in the dataframe will follow the formatting and number of decimals needed. In our case 3 digits after decimal separator for numeric values.

**Line 56** creates a new dataframe with the updated columns. Take extra care in listing all columns (both unchanged and updated with the previous lines) to ensure all data will be retrieved for heatmap construction.

### Column names

**Lines 58 to 72** serve the purpose of changing names to columns, in order to label elements of heatmaps produced by `ggplot()` functions (**in this case column order matter**s).

> We suggest you refrain from using dashes ("-") in your column names, to avoid R forcing them into dots.

### Percentages

For decimal values that should be expressed as percentages, **lines 85 to 89** deliver examples of arithmetic transformation, and variable assignment to do that. Skip if your dataframe's values are alreday expressed as percentages.

## Level management

### Level names

**Lines 92 to 99** are used to rename levels of mutations. Add as many levels to rename as the number of levels (mutations) in your dataframe.

### Level reorganization

**Lines 102 to 116** reorder levels of mutations (**102-104**), preprocessing (**107-110**) and models (**113-116**). Level reorganization is important for the order in which all variables will be displayed in the graphs.

## Creation of labels

Before the actual plot construction, **lines 119 to 147** provide examples for the creation of new vectors for heatmap labels. Some data may be difficult to display or may not be shown properly in the heatmap's cells. Using a different vector for labels will allow to improve graphic layout without impairing data representation.

In this example the vectors created for labels are for:
- accuracy (**lines 119 to 127**), which incorporates CI95 range;
- Cohen's kappa (**line 130**), F1-score (**lines 133**) and AUC from ROC, PR and PRG (**lines 136 to 138**), in which the number of decimal digits is  adjusted;
- execution time (**lines 141 to 147**), which will now be displayed using the most appropriate time unit.

## Plot heatmaps

From **line 153** starts the actual heatmap creation.

Picture names, dataset names, spacing and dimensions should be adapted to personal needs and formats. This step is the reason why the script is highly dependent on the characteristics of the dataframe and its contents. If our settings do not fit your data, explore the `ggplot` function and test what it's best for you. 

For the user's ease, a general reference to used functions from the `ggplot2` package will be given below: the same structure and functions are shared among all `ggplot` blocks of code.

* `ggplot(data, aes(x, y))` : defines the dataframe from which to retrieve data and what vectors will be represented on the 2 axes.
* `theme_light()` : enables light theme and defines base text size.
* `geom_tile(aes(fill))` : defines the vector with which to populate the heatmap's tiles (or cells).
* `geom_text(size, aes(label))` : determines the vector whose values will be displayed as text labels in the tiles and the text size. If you have a different vector for labels, put it here, otherwise  the vector should be the same used for `geom_tile(aes(fill))`.
* `scale_fill_gradient(low, high, name, limits, breaks)` : controls the characteristics of the gradient gauge bar. 
    * `low` and `high` : define the colors to be assigned to the upper and lower limit of the bar, respectively, thus determining the color gradient with which the heatmap tiles will be colored, based on the value in them. For both of them, the argument that should be passed is the desired color's hexadecimal code. 
    * `name` : allows to display test as a title.
    * `limits` : determines the values for the upper and lower limits of the gauge bar.
    * `breaks` : controls the values at which breaks and corresponding labels should be displayed for the gauge bar.
* `theme()` : contains all available controls (color, position and size, each where applicable) for elements of the plot (legend, borders, title, background and grid).

* After plotting, each graph is saved as pdf, tiff and png, through the function `ggsave()`.

The ggplot blocks of code in the script will produce the following graphs in this order:

* **Heatmap for average accuracy of prediction on the testing group**

* **Heatmap for F1-score (global)**

* **Heatmap for Cohen kappa of the testing group**

* **Heatmap for ROC-AUC (Receiver Operating Characteristic Area Under the Curve)**

* **Heatmap for PR-AUC (Precision-Recall Area Under the Curve)**

* **Heatmap for PRG-AUC (Precision-Recall-Gain Area Under the Curve)**

* **Heatmap for execution time**

With the provided demo dataframe, heatmaps will look like those in the following pictures:

![](./demo/Demo-heatmaps-A-D.png)

![](./demo/Demo-heatmaps-AUC.png)

## Correlation plots

Heatmappe.R also provides code to perform Pearson correlation tests (**lines 473-509**), which calculate R2 and p-value, and to produce correlation plots between two metrics (**lines 512-803**).

**Line 516** reverses the order of model levels, which is necessary to keep the same order as the heatmaps in the correlation plots, due to a different element layout in the latter.

Correlation plots are produced with `ggplot()` just like heatmaps, so please refer to the [previous section](#plot-heatmaps) and to R's help page for ggplot2 (`help("ggplot")`) for reference.

Below are a couple of examples for correlation plots. Bear in mind that the dataset used in these esamples is still artificial.

![](./demo/Demo-correlation-accuracy_training_vs_ROC-AUC.png)
![](./demo/Demo-correlation-accuracy_training_vs_PRG-AUC.png)

# `boxplots.R`

The `boxplots.R` script is provided for peer-review reproducibility of dataset graphical representations.

The script is divided in 3 sections:
- the first one performs the last dataset filtration steps for samples used in the "Source attribution of Listeria monocytogenes through a versatile supervised machine learning workflow based on genomic data" research article, (*i.e.* dataframe subset based on quality, Phred score, contigs, depth and breadth of coverage);
- the second part executes Kolmogorov-Smirnov tests;
- the third and last part performs the graphical representation.

# `ANOVA.R`

The ANOVA.R script is provided for peer-review reproducibility of statistical analyses. The first part reads and prepares the same input dataframe used for heatmappe.R. It's used to ensure the dataset is being processed in the same way as the heatmappe.R script.

The second part prepares to and performs the statistical analyses on machine learning prediction results (mwANOVA and owANOVA, Tukey multiple pairwise-comparisons, pairwise t-tests, Kruskal-Wallis rank sum test, Shapiro-Wilk test, Levene’s test and Wilcoxon test).

# References

- P. Castelli, A. De Ruvo, C. Cammà, A. Di Pasquale and N. Radomski. "Supervised machine learning workflow based on genomic data : Case study of Listeria monocytogenes source attribution". 2023, XXXXX, X(XXXX): X-XX, doi.org/XXXXXXXXXX

# Acknowledgements

My colleague and mentor Nicolas Radomski for directions, guide and collaboration in writing the three R scripts provided in this repository.

# Authors

Pierluigi Castelli and Nicolas Radomski 