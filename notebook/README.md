# GTEx Gene Expression Research

## Overview

This project analyzes gene expression across human tissues using data from the **Genotype-Tissue Expression (GTEx) Project**.

The main goal is to study the relationship between genes and tissues and determine which tissues are most strongly associated with a selected group of genes.

The project uses median TPM gene expression values from GTEx and converts them into relative expression values so that genes with very high overall expression do not automatically have more influence than genes with lower overall expression.

The analysis is being developed in Python using Jupyter Notebook and VS Code.

---

## Dataset

The dataset used in this project is:

`GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct.gz`

The data contains:

- Ensembl gene IDs
- Gene names
- 68 tissue and cell-type columns
- Median TPM expression values for each gene-tissue pair

The structure of the data is approximately:

```text
Gene ID | Gene Name | Tissue 1 | Tissue 2 | Tissue 3 | ...
```

Each row represents one gene, and each tissue column contains the median TPM expression value for that gene.

The raw GTEx data is not uploaded to this repository because the `data/` folder is ignored by Git.

---

## What is TPM?

TPM stands for **Transcripts Per Million**.

TPM is a normalized measurement used to represent gene expression. In simple terms, a larger TPM value means that a gene has a higher level of expression in that tissue.

For example:

```text
Gene A

Brain = 10 TPM
Liver = 30 TPM
Heart = 10 TPM
```

This tells us that Gene A has higher expression in the liver compared with the brain and heart.

---

## Current Analysis

### 1. Load the GTEx Data

The GTEx `.gct.gz` file is loaded using Pandas.

Because GCT files contain two metadata rows before the real table header, the first two rows are skipped.

```python
df = pd.read_csv(
    file_name,
    sep="\t",
    skiprows=2
)
```

---

### 2. Identify Gene and Tissue Columns

The first two columns contain gene information:

```text
Name
Description
```

The remaining columns contain tissue expression values.

```python
tissue_cols = df.columns[2:]
```

This gives the list of tissues that will be used in the analysis.

---

### 3. Create the Gene-Tissue Table

A new table is created containing the gene ID, gene name, and tissue expression values.

```python
gene_tissue_table = df[
    ["Name", "Description"] + list(tissue_cols)
].copy()
```

The table has this structure:

```text
Name | Description | Brain | Liver | Heart | ...
```

This table represents the gene-to-tissue expression network that will be used for the rest of the project.

---

### 4. Check the Data

Before doing calculations, I checked the data types and looked for missing values.

The gene ID and gene name columns contain strings, while the tissue expression columns contain numeric values.

I also checked for missing values using:

```python
gene_tissue_table.isnull().sum()
```

No missing values were found in the dataset.

---

### 5. Find Genes With Zero Expression

Some genes have a TPM value of 0 in every tissue.

To identify these genes, I first created a table containing only the tissue expression values:

```python
expression_table = gene_tissue_table[tissue_cols]
```

Then I added the tissue values across each gene:

```python
zero_genes = (expression_table.sum(axis=1) == 0)
```

If the total expression for a gene is 0, this means that the gene has no expression in any of the tissues in this dataset.

These genes are removed before normalization because dividing their expression by their total expression would require dividing by 0.

They also do not provide information about tissue preference.

---

## Normalization

### Why Normalize the Data?

Different genes can have very different overall TPM values.

For example:

```text
Gene A

Brain = 1000
Liver = 3000
Heart = 1000
```

and:

```text
Gene B

Brain = 1
Liver = 3
Heart = 1
```

Gene A has much larger raw values, but both genes actually have the same expression pattern across the tissues.

If raw TPM values were used directly, Gene A could have much more influence on the analysis simply because its overall expression is higher.

To make the comparison more fair, each gene is normalized across all tissues.

The method used in this project follows the research writeup, which states that gene-to-tissue expression values should be converted into relative expression values so that the edges leaving each gene sum to 1.

For each gene:

```text
Normalized expression =
TPM in one tissue / Total TPM across all tissues
```

For example:

```text
Brain = 10
Liver = 30
Heart = 10

Total = 50
```

After normalization:

```text
Brain = 10 / 50 = 0.20
Liver = 30 / 50 = 0.60
Heart = 10 / 50 = 0.20
```

The values now add up to:

```text
0.20 + 0.60 + 0.20 = 1.00
```

The normalized values can also be understood as the relative distribution of that gene's expression across tissues.

---

### Normalization Code

```python
expression_table = gene_tissue_table[tissue_cols]

gene_total = expression_table.sum(axis=1)

normalized_expression = expression_table.div(
    gene_total,
    axis=0
)
```

The gene information is then added back to the normalized expression table:

```python
normalized_gene_tissue_table = pd.concat(
    [
        gene_tissue_table[["Name", "Description"]].reset_index(drop=True),
        normalized_expression.reset_index(drop=True)
    ],
    axis=1
)
```

The normalization is checked using:

```python
normalized_gene_tissue_table[tissue_cols].sum(axis=1).head()
```

Each gene should have a total close to:

```text
1.0
```

---

## Tissue Score Calculation

After normalization, the next goal is to determine which tissues are most strongly associated with a selected group of genes.

The research writeup defines the tissue score as:

```text
S(T) = Σ wi × e(gi,T)
```

where:

- `T` is a tissue
- `gi` is one gene from the input gene set
- `wi` is the weight assigned to that gene
- `e(gi,T)` is the normalized expression of that gene in tissue `T`

If no special gene weights are provided, the weight can default to:

```text
wi = 1
```

In that case, the tissue score is simply the sum of the normalized expression values for the selected genes.

---

### Simple Example

Assume there are two genes:

```text
Gene A

Brain = 0.20
Liver = 0.60
Heart = 0.20
```

and:

```text
Gene B

Brain = 0.20
Liver = 0.10
Heart = 0.70
```

If both genes have weight 1, then:

```text
Brain score
= 0.20 + 0.20
= 0.40
```

```text
Liver score
= 0.60 + 0.10
= 0.70
```

```text
Heart score
= 0.20 + 0.70
= 0.90
```

The result would be:

```text
Heart = 0.90
Liver = 0.70
Brain = 0.40
```

This means that this group of genes has the strongest combined relative expression in the heart.

---

### Tissue Score Code

After selecting an input gene set:

```python
selected_genes = normalized_gene_tissue_table[
    normalized_gene_tissue_table["Description"].isin(input_genes)
]
```

the tissue scores can be calculated using:

```python
tissue_scores = selected_genes[tissue_cols].sum(axis=0)
```

The tissues are then ranked from highest score to lowest:

```python
tissue_scores = tissue_scores.sort_values(
    ascending=False
)
```

Finally, the scores are converted into a table:

```python
tissue_score_table = tissue_scores.reset_index()

tissue_score_table.columns = [
    "Tissue",
    "Score"
]
```

A higher tissue score means that the selected gene set has a stronger combined relative expression in that tissue.

---

## Why Use a Selected Gene Set?

The goal is not only to find which tissue has the highest expression across every human gene.

Instead, the goal is to take a biologically meaningful group of genes and ask:

> Which tissues are these genes most strongly associated with?

The input genes could come from:

- a GWAS study
- disease-associated genes
- genes from an experiment
- differentially expressed genes
- a biological pathway
- a curated gene list

This makes the tissue score useful for interpreting gene sets that may contain a tissue-specific signal.

---

## Statistical Testing

A high tissue score by itself does not tell us whether the score is unusual.

For example, a tissue might naturally receive high scores from many different gene sets.

The next step of the method is therefore to compare the real gene set against random gene sets.

The research writeup proposes generating:

```text
1000 random gene sets
```

with the same number of genes as the real input set.

For every random gene set:

1. Select random genes.
2. Calculate tissue scores using the same method.
3. Store the random scores.
4. Compare the real tissue score with the random distribution.

This allows an empirical p-value to be calculated.

Conceptually:

```text
Real gene set
      ↓
Real tissue score
      ↓
Compare with
      ↓
1000 random gene sets
      ↓
Random tissue scores
      ↓
Empirical p-value
```

A small empirical p-value would indicate that obtaining a tissue score that large by chance is uncommon.

This statistical-testing section is the next stage of the project.

---

## Relationship to Enrichr

This research is related to the idea behind **gene-set enrichment analysis**.

Enrichr takes a list of genes and compares that list with existing biological gene-set libraries to determine which terms, pathways, tissues, transcription factors, diseases, or other biological categories are enriched.

The original Enrichr method uses several calculations, including:

- Fisher's exact test
- a rank-based z-score
- a combined enrichment score

The combined score described in the original Enrichr paper is:

```text
Combined Score = log(p-value) × z-score
```

The purpose is to combine information about statistical overlap with how unusually high a biological term ranks.

The 2016 update to Enrichr also introduced the idea of **fuzzy gene sets**.

Traditional gene sets use:

```text
0 = gene is not in the set
1 = gene is in the set
```

A fuzzy gene set allows values between 0 and 1:

```text
0.0 = no membership
0.4 = partial membership
1.0 = full membership
```

This is relevant to this project because the normalized GTEx values also represent gene-to-tissue relationships using values between 0 and 1.

However, the exact tissue-score calculation implemented in this project follows the research writeup provided for this project rather than directly reproducing the Enrichr scoring algorithm.

---

## Research Workflow

The current project workflow is:

```text
GTEx median TPM data
        ↓
Load and inspect data
        ↓
Create gene-tissue table
        ↓
Check missing values
        ↓
Identify genes with zero expression
        ↓
Remove all-zero genes
        ↓
Normalize each gene across tissues
        ↓
Verify each gene sums to 1
        ↓
Create normalized gene-tissue table
        ↓
Choose an input gene set
        ↓
Calculate tissue scores
        ↓
Rank tissues by score
        ↓
Generate random gene sets
        ↓
Calculate empirical p-values
        ↓
Gene subnetwork analysis
```

---

## Current Project Status

Completed:

- Loaded GTEx median TPM data
- Identified gene and tissue columns
- Created the gene-tissue table
- Checked data types
- Checked for missing values
- Identified genes with zero expression across all tissues
- Removed all-zero genes
- Normalized gene expression across tissues
- Verified that each normalized gene sums to 1
- Created the normalized gene-tissue table
- Implemented the first tissue-score calculation
- Created a ranked tissue-score table

Next steps:

- Test tissue scoring using meaningful input gene sets
- Generate random gene sets
- Calculate empirical p-values
- Study the distribution of tissue scores
- Continue toward gene subnetwork analysis

---

## Tools

- Python
- Pandas
- NumPy
- Jupyter Notebook
- VS Code
- Git
- GitHub

---

## References

### Enrichr — Original Paper

Chen, E. Y., Tan, C. M., Kou, Y., et al.  
**Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool.**  
BMC Bioinformatics, 2013, 14:128.

This paper describes the original Enrichr gene-set enrichment method, including Fisher's exact test, expected-rank correction, z-scores, combined enrichment scores, and gene-set libraries.

### Enrichr — 2016 Update

Kuleshov, M. V., Jones, M. R., Rouillard, A. D., et al.  
**Enrichr: a comprehensive gene set enrichment analysis web server 2016 update.**  
Nucleic Acids Research, 2016, 44.

This update expands Enrichr's gene-set libraries, adds GTEx-derived expression libraries, and introduces fuzzy gene-set enrichment where genes can have membership values between 0 and 1.

### Research Method Writeup

The tissue normalization, tissue-score calculation, randomization strategy, empirical p-value calculation, and planned gene subnetwork analysis in this repository follow the research writeup provided for this project.

---

## Author

**El Mehdi Bensouda**  
Data Science Student  
Duquesne University

---

## Project Status

This project is currently **in progress**.