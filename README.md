# CPMM: Longitudinal Plasma Proteomics Analysis Package

**CPMM** is designed for **longitudinal plasma proteomics analysis**, with a special focus on modeling disease onset and progression (e.g., Alzheimer's disease). It integrates preprocessing, change-point mixed models, statistical tests, survival analysis, and visualization.

## Table of Contents

-   [Installation](#installation)
-   [Data Requirements](#data-requirements)
-   [Quick Start](#quick-start)
-   [Data Preprocessing](#data-preprocessing)
-   [Statistical Modeling](#statistical-modeling)
-   [Visualization](#visualization)
-   [Pathway Analysis](#pathway-analysis)
-   [Survival Analysis](#survival-analysis)
-   [Predictive Modelling](#predictive-modelling)
-   [Function Reference](#function-reference)

## Installation

``` r
# Install from GitHub 
devtools::install_github("xiaoqinghuanglab/CPMM")
library(CPMM)
```

The package can also be installed using the tar.gz file from the zip folder
```r
install.packages("path/to/CPMM_0.0.0.9000.tar.gz", repos = NULL, type = "source")
library(CPMM)
```
Installing via RStudio GUI
- Go to Tools → Install Packages

- Set "Install from" to Package Archive File (.zip, .tar.gz)

- Locate and select your .tar.gz file

- Click Install


## Data Requirements

### Input Data Structure

Most functions expect longitudinal data frames with the following columns:

-   `SUBID` – subject identifier
-   `PROCEDURE_AGE` – age at each visit
-   `ONSET_AGE` – disease onset age
-   `SEX` – biological sex (factor)
-   `BASELINE_AGE` – age at baseline visit
-   `CATEGORY` – diagnosis category (Normal, SCD, MCI, AD Dementia, FTD Dementia)
-   Protein columns – numeric expression values for 100s–1000s of proteins

### Required CSV Files

-   `df_all` - CSV file with all the data
-   `df_normal_only` - CSV file with normal only patients
-   `df_abnormal_only` - CSV file with abnormal only patients
-   `df_status_change` - CSV file with patients who converted from Normal to Abnormal

All the data frames must contain same numbers of proteins

## Quick Start

This Code shows how to get started with fitting the model for all the proteins/genes, visualize them and run Wald test on the results for a toy dataet

``` r
library(RPackage)

# Load your data
df_all <- read.csv("./data/df_all_toy.csv")
df_normal_only <- read.csv("./data/df_normal_only_toy.csv")
df_abnormal_only <- read.csv("./data/df_abnormal_only_toy.csv")
df_status_change <- read.csv("./data/df_status_change_toy.csv")
df_pathway <- read.csv("./data/toy_pathways.csv")
df_mci <- read.csv("./data/df_mci_toy.csv")

# Select the proteins
proteins <- colnames(df_status_change)[9:27]

# Fit change-point mixed models across proteins
results <- fit_cpmm_all_proteins(
  df_status_change = df_status_change,
  protein_list = proteins,
  covariates = c("SEX", "BASELINE_AGE"),
  subject_id_col = "SUBID",                     
  years_since_onset_col = "years_since_onset"
)

# Visualize CPMM trajectory for a single protein
plot_cpmm(
  df_status_change,
  protein = "P1",
  covariates = c("SEX","BASELINE_AGE"),
  years_since_onset_col = "years_since_onset"
)

# Perform Wald test for significant slope changes
wald <- compute_wald_test(
  results_df = results,
  adjust_p = TRUE,
  alpha = 0.05             # User defined threshold
)
```

## Statistical Modeling

### Change-Point Mixed Models (CPMM)

The `fit_cpmm_all_proteins()` function estimates slopes before and after onset for status-change subjects:

``` r
results <- fit_cpmm_all_proteins(
  df_status_change = df_status_change,
  protein_list = proteins,
  covariates = c("SEX","BASELINE_AGE"),
  subject_id_col = "SUBID",
  years_since_onset_col = "years_since_onset"
)
```

#### Result Columns:

-   **Beta 1, SE Beta 1** = pre-onset slope & SE
-   **Beta 2, SE Beta 2** = post-onset slope & SE
-   Model metrics (AIC, BIC, MSE) per group

### Wald Test

Find proteins with significant slope changes:

``` r
wald <- compute_wald_test(
  results_df = results,
  adjust_p = TRUE,    # BH FDR correction
  alpha = 0.05        # User defined threshold
)
```

### Mann-Whitney U Test

Compare distributions within categories:

``` r
# Prepare combined expression data frame
expr_df <- prepare_combined_expression(
  inputs = list(
    list(df = df_normal_only, label = "Normal_only"),
    list(df = df_status_change, label = "Status_change"),
    list(df = df_abnormal_only, label = "Abnormal_only"),
    list(
      df = df_all,
      label = "MCI",
      category_col = "CATEGORY",
      categories = "MCI"
    )
  ),
  subset_genes = proteins
)

# Perform Mann-Whitney U test
mw <- compare_groups_mannwhitney(
  combined_expr = expr_df,
  gene_list = c("P1","P2","P3"),
  group1 = "Normal_only",
  group2 = "MCI"
)
```

## Visualization

### CPMM Trajectory Plots

Visualize protein trajectories for proteins of interest:

``` r
plot_cpmm(
  df_status_change,
  protein = "P6",                                     # protein/gene of inerest
  covariates = c("SEX","BASELINE_AGE"),
  subject_id_col = "SUBID",
  years_since_onset_col = "years_since_onset"
)
```
![Protein Trajectories](./assets/cpmm_plot.png)

### Volcano Plot

``` r
plot_wald_volcano(
  wald_df = wald,
  pval_col = "Adjusted P-value",
  annotate = TRUE
)
```
![Volcano Plot](./assets/volcano.png)

### Quadrant Plot

``` r
plot_quadrant_beta(
  wald_df = wald,
  beta_x_col = "Beta 1",
  beta_y_col = "Beta 2",
  fdr_col = "Adjusted P-value",
  annotate = TRUE
)
```
![Quadrant Plot](./assets/quadrant.png)

### Expression Boxplots

``` r
plot_expression_boxplot(
  expr_df,
  gene_order = proteins,  # genes/proteins of interest
  hue_col = "Source",
  expression_col = "Expression",
  gene_col = "Gene"
)
```
![Expression Plot](./assets/expr.png)

## Pathway Analysis

### Data Requirements

Pathway analysis requires a data frame with the following columns:

-   `Pathway` - Pathway Names
-   `Gene` - Gene Enriched by the Pathway
-   `Source` - The tool used for enrichment (e.g., DAVID, Metascape)
-   `CategoryGroup` - Category of the Pathway provided by the tool
-   `LogQValues` - Log Q or Log P values of the pathways
-   `Cleaned Pathway` - Cleaned Pathway names for better readability
-   `Category` - Category Name to categorize the pathways (manually defined)

### Pathway Visualization

#### Bubble Plot

``` r
plot_pathway_bubble(
  df = df_pathway,
  pathway_col = "Cleaned_Pathway",
  category_col = "Category",
  source_col = "Source",
  logq_col = "LogQValue",
  gene_col = "Gene",
  title = "Pathway Enrichment by Source",
  size_scale = 15  # maximum bubble size
)
```
![Bubble Plot](./assets/bubble.svg)

#### Heatmap

``` r
plot_pathway_gene_heatmap(
  df = df_pathway,
  pathway_col = "Cleaned_Pathway",
  category_col = "Category",
  gene_col = "Gene",
  title = "Pathway–Gene Membership Heatmap"
)
```
![Heatmap Plot](./assets/heatmap.svg)

#### Top Pathways Bar Plot

``` r
plot_top_pathways_bar(
  df = df_pathway,
  pathway_col = "Cleaned_Pathway",
  gene_col = "Gene",
  logq_col = "LogQValue",         # logq (or logp) column
  category_col = "Category",
  top_n = 6,                      # number of pathways in the plot
  annotate = TRUE                 # annotate with logq (or logp) values if available
)
```
![Bar Plot](./assets/bar.svg)

## Survival Analysis

Compute time-to-threshold events per subject and plot Kaplan-Meier curves:

``` r
plot_km_with_threshold(
  biomarker_name = "P1",
  threshold = 3,
  groups = list(
    Normal = list(df = df_normal_only),
    Status_change = list(df = df_status_change),
    SCD = list(
      df = df_status_change,
      filter = CATEGORY == "SCD"
    )
  ),
  time_points = seq(-6, 6, by = 2)
)

```
![Survival Plot](./assets/survival.png)

## Predictive Modelling

LASSO-Selected Logistic Regression with CV ROC

```r
# Select the proteins for the LASSO
prot <- paste0("P", 1:6)

roc_res <- plot_lasso_cv_roc(
  df = df_mci,
  protein_features = prot,
  covariates = c("PROCEDURE_AGE","Sex_bin"),
  label_col = "label",
  n_folds = 5,
  title = "Cross-validated ROC"
)

roc_res$mean_auc
roc_res$final_features  # Shows the final selected features by the LASSO
roc_res$plot            # If the LASSO didn't select any features a Warning appears: "No protein features were selected by LASSO; model uses covariates only."
```
![AUC Plot](./assets/auc.png)

## Data Preprocessing

These functions can be used to get `years_since_onset`, setting onset age for normal patients, correcting onset age and status to make a unidirectional flow and to separate status change patients from the data. It is not an requirement to run these function when using the toy dataset as the required columns are alredy present in the dataset.

### Calculate Years Since Onset

Compute time relative to onset: age - onset_age

``` r
calculate_years_since_onset(
  df, 
  age_col = "age", 
  onset_age_col = "onset_age", 
  new_col = "years_since_onset"
)
```

### Set Onset Age For Normal Subjects

For subjects always Normal, set onset to their max observed age (so all timepoints are pre-onset):

``` r
set_onset_age_for_normals(
  df,
  subject_id_col = "SUBID", 
  status_col = "status_raw",
  age_col = "age", 
  mutated_col = "DECAGE"
)
```

### Correct Status by Onset

Fix label inconsistencies relative to onset (post-onset "Normal" → "Abnormal", etc.):

``` r
correct_status_by_onset(
  df, 
  age_col = "AGE", 
  onset_age_col = "onset_age", 
  status_col = "status_cleaned"
)
```

### Enforce Unidirectional Status Change

"Once Abnormal, always Abnormal" (applies forward in time per subject):

``` r
enforce_unidirectional_status_change(
  df,
  subject_id_col = "SUBID", 
  status_col = "status_cleaned", 
  date_col = "procedure_date"
)
```

### Identify Status Change Patients

Keep only subjects who change from Normal → Abnormal:

``` r
identify_status_change_subjects(
  df,
  subject_id_col = "SUBID",
  status_col = "status_cleaned",
  date_col = "procedure_date"
)
```

## Function Reference {#function-reference}

### Data Requirements by Function Type

-   **Longitudinal frames** (df_status_change, df_normal_only, df_abnormal_only):
    -   Columns: SUBID, years_since_onset, covariates (e.g., SEX, BASELINE_AGE), and protein columns (numeric)
    -   Optional: PROCEDURE_AGE, ONSET_AGE if you need to compute years_since_onset
-   **CPMM results** for Wald (results):
    -   Columns: Protein, Beta 1, SE Beta 1, Beta 2, SE Beta 2
-   **Expression long table** (expr_df):
    -   Columns: Gene, Expression, Source
-   **Pathway table** (path_df):
    -   Columns: Cleaned_Pathway, BioCategory_Manual, Source, LogQValue, Gene

### Help Documentation

For detailed function documentation, use:

``` r
?fit_cpmm_all_proteins
?plot_wald_volcano
?compute_wald_test
# ... and other function names
```








