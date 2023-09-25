# ComBatFamily Pipeline

## Summary

The ComBatFamily Pipeline includes both ComBat Family harmonization methods and a visualization tool that provides interactive diagnostics and statistical analysis to evaluate how batch effects impact the data.

To simplify the harmonization process, ComBatFamily Pipeline seamlessly integrates both the harmonization methods and the visualization tool into a single, unified command line interface (CLI). This integration enables users to execute both harmonization and visualization steps in a smooth and efficient manner. Beyond, the ComBatFamily Pipleline also includes a visualization option for investigating life span age trend of brain structures.

In a short summary, the following two command line interfaces are included in the ComBatFamily pipeline:

-   **Harmonization & Visualization**: combatfam_CLI.R
    1) Batch effect visualization
    2) Data Harmonization
-   **Age Trajectory**: age_CLI.R
    1) Life span age trend of brain structures visualization

## Usage

Use the following command to access the Harmonization & Visualization pipeline:

```
Rscript combatfam_CLI.R --help
```

Use the following command to access the Harmonization & Visualization pipeline:

```
Rscript age_CLI.R --help
```

This step will display the available commands and options for both pipelines.

## Package Installation

-   ComBatFamily

```{r}
library(devtools)

devtools::install_github("Zheng206/ComBatFam_Pipeline/ComBatFamily")

```

-   ComBatFamQC

```{r}
library(devtools)

devtools::install_github("Zheng206/ComBatFam_Pipeline/ComBatFamQC")
```