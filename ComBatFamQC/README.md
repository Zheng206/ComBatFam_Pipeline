# ComBatFamQC

The **ComBatFamQC** package is a powerful tool designed to streamline statistical analysis and interactive visualization for harmonization quality control needs. This package is sepcifically tailored for evaluating batch effects before and after applying ComBatFamily harmonization method, and providing life span age trends of brain structures if needed. In terms of the final delivery, it will provide an interactive visualization through R shiny.

## Package Features

The ComBatFamQC package offers the following three key functionalities:

-   **Statistical Analysis and Harmonization**: ComBatFamQC simplifies the process of performing statistical analyses to detect potential batch effects. It provides you with all relevant statistical test results for batch effect visualization and evaluation. This step also includes the harmonization process and a dataset after harmonization will be returned. 

-   **Interactive Visualization through R shiny**: The ComBatFamQC package comes with an interactive visualization tool built on R shiny, providing an intuitive user interface to explore and evaluate batch effects. The output is organized into multiple tabs, which includes:

    -   **Summary**: Sample Size and Covariate Distribution
    -   **Residual Plot**: Additive and Multiplicative Batch Effect
    -   **Residual Dimensionality Reduction**: PCA & T-SNE
    -   **Empirical Bayes Assumption Check**: Location and Scale Paramaters Distribution (check whether EB priors' distribution overlaps well with empirical distributions)
    -   **Statistical Test**:
        -   *Batch Effect Test*: MDMR, Kenward-Roger (liner mix model), ANOVA, Kruskal-Wallis
        -   *Equality of Variance Test*: Fligner-Killeen, Levene's Test, Bartlett's Test
    
-   **Age Trajectory** \
    Sex adjustment and customized lower bounds and upper bounds of quantiles to look at are enabled as well.
    -  **Age Trend Plots**
    -  **Age Trend Table** 

## Installation

```{r}
library(devtools)

devtools::install_github("Zheng206/ComBatFam_Pipeline/ComBatFamQC")

```