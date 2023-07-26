# ComBatFamQC

The **ComBatFamQC** package is a powerful tool designed to streamline statistical analysis and interactive visualization for harmonization quality control needs. This package is sepcifically tailored for evaluating batch effects before and after applying ComBatFamily harmonization method. In terms of the final delivery, it will provide an interactive visualization through R shiny.

## Features

The ComBatFamQC package offers the following two key functionalities:

-   **Statistical Analysis Results of Batch Effects for Visualization**: ComBatFamQC simplifies the process of performing statistical analyses to detect potential batch effects. It provides you with all relevant statistical test results for batch effect visualization and evaluation.

-   **Interactive Visualization Tool with R shiny**: The ComBatFamQC package comes with an interactive visualization tool built on R shiny, providing an intuitive user interface to explore and evaluate batch effects. The output is organized into multiple tabs, which includes:

    -   **Summary**: Sample Size and Covariate Distribution
    -   **Residual Plot**: Additive and Multiplicative Batch Effect
    -   **Residual Dimensionality Reduction**: PCA & T-SNE
    -   **Empirical Distribution**: Location and Scale Paramaters Distribution
    -   **Statistical Test**:
        -   *Batch Effect Test*: MDMR, Kenward-Roger (liner mix model), ANOVA, Kruskal-Wallis
        -   *Equality of Variance Test*: Fligner-Killeen, Levene's Test, Bartlett's Test

## Installation

```{r}
library(devtools)

devtools::install_github("Zheng206/ComBatFam_Pipeline/ComBatFamQC")

```