# miRNA qPCR Analysis Platform: Comprehensive Methods Description

## Overview

The miRNA qPCR Analysis Platform is a web-based bioinformatics tool developed using R Shiny for comprehensive analysis of microRNA (miRNA) expression data obtained from quantitative real-time PCR (qPCR) experiments. The platform integrates multiple analytical methods to provide fold change calculations, reference gene stability assessment, statistical comparisons, and predictive modeling capabilities for differential miRNA expression between normal and tumor samples.

## Software Implementation

**Platform:** R Shiny (Interactive web application)  
**Programming Language:** R (version ≥4.0.0)  
**Key Dependencies:** 
- Data processing: dplyr, tidyr, readxl
- Visualization: ggplot2, plotly, DT
- Statistical analysis: pROC, glmnet
- Reference gene analysis: geNorm algorithm implementation

## Input Data Requirements

### Required Columns
- **Well**: Well position identifier
- **Fluor**: Fluorophore used
- **Target**: miRNA or reference gene name
- **Sample**: Sample identifier
- **Cq/Ct**: Quantification cycle (threshold cycle) value
- **Cq Mean**: Mean Cq value across technical replicates
- **Type**: Sample classification (Normal/Tumor/Control)
- **Sd**: Standard deviation of Cq values

### Accepted Reference Genes
The platform recognizes the following commonly used reference genes for normalization:
- **U6** (U6 snRNA)
- **UniSP6** (Universal spike-in control)

## Analytical Methods

### 1. Data Normalization and Fold Change Calculation

#### 1.1 ΔCt Method (Delta Ct)
The ΔCt method normalizes target miRNA expression to reference gene(s):

```
ΔCt = Ct(target miRNA) - Ct(reference gene)
```

Where multiple reference genes are present, the geometric mean of their Ct values is used:

```
Ct(reference) = (Ct(ref1) + Ct(ref2) + ... + Ct(refN)) / N
```

#### 1.2 ΔΔCt Method (Comparative Ct Method)
The ΔΔCt method calculates relative expression compared to a calibrator (normal samples):

```
ΔΔCt = ΔCt(sample) - mean(ΔCt(normal samples))
```

#### 1.3 Fold Change Calculation
Fold change is calculated using the 2^(-ΔΔCt) formula:

```
Fold Change (FC) = 2^(-ΔΔCt)
```

**Interpretation:**
- FC = 1: No change in expression
- FC > 1: Upregulation in tumor samples
- FC < 1: Downregulation in tumor samples
- FC = 2: 2-fold increase (100% increase)
- FC = 0.5: 2-fold decrease (50% of normal expression)

### 2. Reference Gene Stability Analysis (geNorm Algorithm)

The platform implements the geNorm algorithm (Vandesompele et al., 2002) to assess reference gene expression stability across all samples.

#### 2.1 geNorm M-value Calculation

For each candidate reference gene, the algorithm:

1. **Calculates pairwise variations:** For each gene pair (j, k), computes the standard deviation of log2-transformed expression ratios across all samples:

```
SDjk = SD[log2(Ctj / Ctk)]
```

2. **Computes gene stability (M-value):** For each gene j, the M-value is the arithmetic mean of all pairwise standard deviations:

```
Mj = mean(SDjk) for all k ≠ j
```

#### 2.2 M-value Interpretation

- **M < 0.5:** Excellent stability (ideal for normalization)
- **M < 1.0:** Good stability (acceptable for normalization)
- **M < 1.5:** Acceptable stability (use with caution)
- **M ≥ 1.5:** Poor stability (not recommended for normalization)

#### 2.3 Additional Stability Metrics

The platform also calculates:

- **Coefficient of Variation (CV):**
```
CV = (SD / Mean) × 100
```

- **Mean Ct value:** Average expression level across all samples
- **Standard Deviation (SD):** Measure of variation in Ct values

### 3. Statistical Testing

#### 3.1 Wilcoxon Rank-Sum Test (Mann-Whitney U Test)

For each miRNA, the platform performs non-parametric comparison between Normal and Tumor groups:

- **Null Hypothesis (H₀):** No difference in miRNA expression between groups
- **Alternative Hypothesis (H₁):** Significant difference exists
- **Test Statistic:** Wilcoxon W statistic
- **P-value calculation:** Exact or asymptotic method depending on sample size

**Significance Levels:**
- ***p < 0.001:** Highly significant (displayed as ***)
- **p < 0.01:** Very significant (displayed as **)
- **p < 0.05:** Significant (displayed as *)
- **p ≥ 0.05:** Not significant (displayed as ns)

**Rationale for Non-parametric Testing:**
- Does not assume normal distribution of fold changes
- Robust to outliers
- Appropriate for small sample sizes (n < 30)
- Handles skewed distributions common in qPCR data

### 4. Receiver Operating Characteristic (ROC) Analysis

#### 4.1 Individual miRNA ROC Curves

For each miRNA, ROC analysis assesses diagnostic/predictive performance:

- **True Positive Rate (Sensitivity):**
```
TPR = TP / (TP + FN)
```

- **False Positive Rate (1 - Specificity):**
```
FPR = FP / (FP + TN)
```

- **Area Under the Curve (AUC):** Measures overall discriminatory ability
  - AUC = 0.5: No discriminatory power (random classifier)
  - 0.5 < AUC < 0.7: Poor discrimination
  - 0.7 ≤ AUC < 0.8: Acceptable discrimination
  - 0.8 ≤ AUC < 0.9: Excellent discrimination
  - AUC ≥ 0.9: Outstanding discrimination

#### 4.2 Combined Predictive Model

The platform employs **Regularized Logistic Regression** with Elastic Net penalty to build multi-miRNA classifiers.

##### 4.2.1 Model Formulation

The logistic regression model estimates the probability of tumor classification:

```
P(Y = Tumor | X) = 1 / (1 + e^(-(β₀ + β₁X₁ + β₂X₂ + ... + βₙXₙ)))
```

Where:
- Y: Binary outcome (Normal = 0, Tumor = 1)
- X₁, X₂, ..., Xₙ: Fold change values for n miRNAs
- β₀: Intercept
- β₁, β₂, ..., βₙ: Regression coefficients

##### 4.2.2 Elastic Net Regularization

To prevent overfitting and perform feature selection, the model uses Elastic Net penalty:

```
Loss = -LogLikelihood + λ[α||β||₁ + (1-α)||β||₂²]
```

Where:
- **λ (lambda):** Regularization strength parameter
- **α (alpha):** Mixing parameter (set to 0.5 for equal L1/L2 penalty)
  - α = 1: LASSO (L1) - promotes sparse solutions
  - α = 0: Ridge (L2) - shrinks coefficients
  - α = 0.5: Elastic Net - balances both properties

**Model Selection:**
- Lambda is automatically selected from the middle of the regularization path
- This approach balances model complexity and predictive performance
- Coefficients shrunk to zero indicate non-informative miRNAs

##### 4.2.3 Model Evaluation

The combined model generates:
- **Combined ROC curve:** Overall predictive performance
- **Combined AUC:** Multi-miRNA signature discrimination ability
- **Non-zero coefficients:** Selected miRNAs contributing to classification
- **Coefficient magnitudes:** Relative importance of each miRNA

### 5. Data Visualization

#### 5.1 Fold Change Dot Plots
- Box plots showing distribution by group
- Individual data points with jitter to avoid overlap
- Mean values marked with diamond symbols
- Log-scale y-axis for better visualization of fold changes
- Reference line at FC = 1 (no change)
- Statistical significance annotations

#### 5.2 Reference Gene Quality Control
- Grouped bar charts showing Ct values
- Faceted by sample type (Normal/Tumor)
- Color-coded by reference gene
- Facilitates visual inspection of normalization consistency

#### 5.3 ROC Curves
- Individual curves for each miRNA
- Combined curve for multi-miRNA model
- AUC values displayed
- Diagonal reference line (AUC = 0.5)

### 6. Summary Statistics and Reporting

For each miRNA, the platform calculates:

- **Sample size (n)** per group
- **Mean fold change** with standard error of the mean (SEM)
- **Standard deviation of fold changes**
- **Mean M-values** with standard deviation
- **Expression direction** (Upregulated/Downregulated/No change)
- **Statistical significance** (p-value)

## Quality Control Features

### 1. Reference Gene Monitoring
- Visualization of reference gene Ct values across samples
- Detection of outliers or inconsistent normalization
- Coefficient of variation assessment

### 2. Data Validation
- Automatic detection of missing values
- Filtering of negative control (NTC) wells
- Validation of required column presence
- Sample size requirements for statistical tests

### 3. Model Diagnostics
- Minimum sample size enforcement (n ≥ 10 for combined modeling)
- Presence of both sample groups (Normal and Tumor)
- Complete case analysis for multi-miRNA models
- Regularization parameter reporting

## Output and Export Capabilities

### 1. Interactive Tables
- Sortable and searchable data tables
- Per-miRNA statistics
- Reference gene stability rankings
- Color-coded stability classifications

### 2. High-Resolution Graphics
- Publication-quality plots (PNG format)
- Interactive plotly visualizations
- Downloadable in multiple formats

### 3. Comprehensive Reports
- Excel workbooks with multiple sheets:
  - Wide-format data with all calculations
  - Per-miRNA detailed results
  - Summary statistics
  - Reference gene metrics
  - Model coefficients

## Statistical Software and Packages

**Core Analysis:**
- glmnet (v4.1+): Elastic Net regularized logistic regression
- pROC (v1.18+): ROC curve analysis and AUC calculation
- stats: Wilcoxon rank-sum test

**Data Manipulation:**
- dplyr (v1.1+): Data transformation and summarization
- tidyr (v1.3+): Data reshaping

**Visualization:**
- ggplot2 (v3.4+): Grammar of graphics plotting
- plotly (v4.10+): Interactive visualizations

## Computational Requirements

**Minimum Requirements:**
- R version ≥ 4.0.0
- 4 GB RAM
- Modern web browser (Chrome, Firefox, Safari, Edge)

**Recommended:**
- R version ≥ 4.3.0
- 8 GB RAM
- Multi-core processor for faster computation

## Limitations and Considerations

1. **Sample Size:** 
   - Minimum n ≥ 3 per group for basic statistics
   - Minimum n ≥ 10 total for combined modeling
   - Small sample sizes reduce statistical power

2. **Reference Gene Selection:**
   - At least one stable reference gene required
   - M-value < 1.0 recommended for reliable normalization
   - Multiple reference genes preferred (n ≥ 2)

3. **Assumptions:**
   - Amplification efficiency near 100% (E ≈ 2)
   - Technical replicates properly averaged
   - No systematic batch effects between groups

4. **Model Interpretation:**
   - Regularization may shrink coefficients to zero
   - Reduced models better for small sample sizes
   - Cross-validation not performed (single dataset analysis)

## Data Privacy and Security

- All analyses performed locally in user's browser session
- No data uploaded to external servers
- Session-based computation (data not retained)
- HIPAA/GDPR compliant for clinical data

## Citation and References

**Key Methodological References:**

1. **ΔΔCt Method:**
   - Livak KJ, Schmittgen TD. Analysis of relative gene expression data using real-time quantitative PCR and the 2^(-ΔΔCt) Method. *Methods*. 2001;25(4):402-408.

2. **geNorm Algorithm:**
   - Vandesompele J, De Preter K, Pattyn F, et al. Accurate normalization of real-time quantitative RT-PCR data by geometric averaging of multiple internal control genes. *Genome Biol*. 2002;3(7):RESEARCH0034.

3. **ROC Analysis:**
   - Hanley JA, McNeil BJ. The meaning and use of the area under a receiver operating characteristic (ROC) curve. *Radiology*. 1982;143(1):29-36.

4. **Elastic Net Regularization:**
   - Zou H, Hastie T. Regularization and variable selection via the elastic net. *J R Stat Soc Series B Stat Methodol*. 2005;67(2):301-320.

5. **Wilcoxon Test:**
   - Wilcoxon F. Individual comparisons by ranking methods. *Biometrics Bulletin*. 1945;1(6):80-83.

## Example Methods Section for Publication

*"miRNA expression analysis was performed using a custom R Shiny-based bioinformatics platform. Quantification cycle (Cq) values were normalized to reference genes (U6 and/or UniSP6) using the ΔCt method. Relative expression levels between tumor and normal samples were calculated using the 2^(-ΔΔCt) method (Livak and Schmittgen, 2001). Reference gene stability was assessed using the geNorm algorithm, with M-values < 1.0 considered acceptable (Vandesompele et al., 2002). Statistical comparisons between groups were performed using the Wilcoxon rank-sum test, with p < 0.05 considered significant. Receiver operating characteristic (ROC) curve analysis was conducted to evaluate the diagnostic performance of individual miRNAs, with area under the curve (AUC) values calculated. A combined predictive model was built using elastic net-regularized logistic regression (α = 0.5) incorporating multiple miRNAs, with performance assessed by ROC-AUC analysis. All statistical analyses were performed in R version 4.3.0 using the glmnet, pROC, dplyr, and ggplot2 packages."*

## Troubleshooting and Support

**Common Issues:**

1. **Missing reference genes:** Ensure U6 or UniSP6 are included in Target column
2. **Low M-values:** Consider using multiple reference genes for better stability
3. **Small sample warnings:** Combine technical replicates or increase biological replicates
4. **Model convergence:** Reduce number of miRNAs or increase sample size

## Version Information

**Current Version:** 1.0.0  
**Last Updated:** 2024  
**License:** Open source (specify license type)  
**Availability:** Web-based platform accessible via modern browsers

## Acknowledgments

This platform integrates established statistical methods and algorithms from the quantitative PCR and bioinformatics communities. The implementation follows best practices for qPCR data analysis as outlined by the MIQE guidelines (Minimum Information for Publication of Quantitative Real-Time PCR Experiments).
