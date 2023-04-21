# PSI_final_models
Final LASSO regression models for predicting smoking status from placental DNA methylation data.

## Use Case and Background
This repository includes source functions and reference data enabling users to predict cigarette/nicotine smoking status from DNA methylation data. Specifically, these models were trained on Illumina EPIC array data with DNA methylation being signified by beta values.

Users can choose from three different models which differ in the cpg probes used during training. These include:

- **Model 1:** Probes from the Illumina EPIC array
- **Model 2:** Probes from the EPIC array which overlap the 450k array
- **Model 3:** Probes from the EPIC array which overlap a meta analysis of smoking associated CpGs (PACE)

## User Guide
In order to use these models:

1. Clone this repository to have access to source functions and reference tables.
2. Load your DNA methylation beta matrix into an R session as a dataframe with CpGs as rows and samples as columns.
3. Load the reference table from one of the 3 models available in the data folder.
4. Use the "psi_predict" function to predict smoking status and calculate a Placenta Smoking Index (PSI) score for each sample.
