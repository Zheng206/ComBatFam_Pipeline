# ComBatFamQC

The **ComBatFamQC** package is a powerful tool designed to streamline statistical analysis and interactive visualization for harmonization quality control needs. This package is sepcifically tailored for evaluating batch effects before and after applying ComBatFamily harmonization method, providing life span age trends of brain structures and residual data sets eliminating specific covariates' effects if needed. In terms of the final delivery, it will provide an interactive visualization through R shiny for batch effect and age trend visualization. Additionaly, it itegrated harmonization process and can provide harmonized data set, fitted combat model, residual data set, fitted regression model etc.

## Diagram
![ComBatFamQC Diagram](/figure/diagram.png)

## Package Features

The ComBatFamQC package offers the following four key functionalities:

1. <u>**Harmonization and Visualization**</u>

-   **Statistical Analysis and Harmonization**: ComBatFamQC simplifies the process of performing statistical analyses to detect potential batch effects. It provides you with all relevant statistical test results for batch effect visualization and evaluation. This step also includes the harmonization process and a dataset after harmonization will be returned. 

-   **Interactive Visualization through R shiny**: The ComBatFamQC package comes with an interactive visualization tool built on R shiny, providing an intuitive user interface to explore and evaluate batch effects. The output is organized into multiple tabs, which includes:

    -   **Summary**: Sample Size and Covariate Distribution
    -   **Residual Plot**: Additive and Multiplicative Batch Effect
    -   **Residual Dimensionality Reduction**: PCA & T-SNE
    -   **Empirical Bayes Assumption Check**: Location and Scale Paramaters Distribution (check whether EB priors' distribution overlaps well with empirical distributions)
    -   **Statistical Test**:
        -   *Batch Effect Test*: MDMR, Kenward-Roger (liner mix model), ANOVA, Kruskal-Wallis
        -   *Equality of Variance Test*: Fligner-Killeen, Levene's Test, Bartlett's Test
2. <u>**Post-Harmonization Downstream Analysis**</u>

-   **Age Trajectory** \
    Generate age trend of each brain structure (roi), adjusting sex and ICV. Customized centiles are enabled as well.
    -  **Age Trend Plots**
    -  **Age Trend Table** 

-   **Residual Generation** \
    Generate residual data set, removing specific covariates' effetcs.


## Installation

```{r}
library(devtools)

devtools::install_github("Zheng206/ComBatFam_Pipeline/ComBatFamQC", build_vignettes = TRUE)

```

## Tutorial

```{r}
vignette("ComBatQC")
vignette("Post-Harmonization")
```