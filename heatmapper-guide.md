#  heatmappe.R guide

## Notes on preparation

> **Our example dataframe has 15 columns. Data extraction from the dataframe is managed through header name of each column, which is hardcoded and should be changed as needed.**

### Column names and levels

> **Lines 58 to 72** handle columns, in order to label elements of heatmaps.
>
> We suggest you refrain from using dashes ("-") in your column names, to avoid R forcing them into dots.

**Lines 92 to 99** are used to rename levels of mutations.
**Lines 102 to 116** reorder levels of mutations, preprocessing and models. Level reorganization is important for the order in which all variables will be displayed in the graphs.

### Creation of labels

**Lines 119 to 147** create new vectors for heatmap labels. Using a different vector for labels will allow to improve graphic layout without impairing data representation.

New vectors are created for labels of:
- accuracy, to incorporate CI95 range;
- Cohen's kappa, F1-score and AUC from ROC, PR and PRG, to adjust the number of decimal digits;
- execution time, to display it using the most appropriate time unit.

## Plots with `ggplot`

A general reference to used functions from the `ggplot2` package will be given below. Code structure and functions are shared among `ggplot` blocks of code (heatmaps and correlation plots).

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

After plotting, each graph is saved as pdf, tiff and png, through the function `ggsave()`.