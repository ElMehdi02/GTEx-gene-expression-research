# GTEx Data Analysis Process

This file explains how I am analyzing the GTEx gene expression dataset and how the research papers I am reading are helping me build the analysis.

The main notebook for this work is:

`data_exploration1.ipynb`

The goal of this analysis is to understand how genes are connected to different human tissues using GTEx gene expression data.

---

## Research Papers I Am Using

I am using two Enrichr research papers to understand how gene-set enrichment analysis works and how gene expression can be connected to biological terms such as tissues.

### 1. Enrichr: Interactive and Collaborative HTML5 Gene List Enrichment Analysis Tool

Chen et al., 2013.

This paper explains how Enrichr takes a list of genes and compares it with existing gene-set libraries.

The paper uses different methods to measure enrichment, including:

- Fisher exact test
- expected-rank correction
- Z-score
- combined score

The main idea I am using from this paper is that a group of genes can be analyzed together to find biological signals that are more represented than expected.

---

### 2. Enrichr: A Comprehensive Gene Set Enrichment Analysis Web Server 2016 Update

Kuleshov et al., 2016.

This paper expands the Enrichr method and adds more gene-set libraries, including gene expression information from GTEx.

Another important idea in this paper is **fuzzy membership**.

Instead of only saying:

```text
0 = gene is not connected
1 = gene is connected
```

a gene can have a value between 0 and 1.

For example:

```text
Gene A = 0.80
Gene B = 0.30
Gene C = 0.05
```

This allows the strength of the relationship to also be included in the analysis.

This idea is useful for my project because I am also converting the GTEx gene expression values into relative values between 0 and 1.

---

# My Analysis Process

## 1. Load the GTEx Dataset

I started by loading the GTEx median TPM dataset using Pandas.

The file contains:

- gene IDs
- gene names
- tissue names
- median TPM values

Each row represents a gene and each tissue column contains the expression of that gene in that tissue.

The first two rows of the GCT file contain information about the file, so I skip them when loading the data.

```python
df = pd.read_csv(
    file_name,
    sep="\t",
    skiprows=2
)
```

---

## 2. Identify the Tissue Columns

The first two columns contain information about the genes:

```text
Name
Description
```

The remaining columns are the tissues.

I selected the tissue columns using:

```python
tissue_cols = df.columns[2:]
```

This allows me to work only with the expression values when I need to perform calculations.

---

## 3. Create a Gene-Tissue Table

I created a table that contains:

```text
Gene ID | Gene Name | Tissue 1 | Tissue 2 | Tissue 3 | ...
```

Each value in the tissue columns is the median TPM expression for that gene.

```python
gene_tissue_table = df[
    ["Name", "Description"] + list(tissue_cols)
].copy()
```

This table is the main table I use for the analysis.

---

## 4. Check the Dataset

Before doing calculations, I checked the structure of the data.

I checked:

- number of genes
- number of tissues
- column names
- data types
- missing values

For missing values I used:

```python
gene_tissue_table.isnull().sum()
```

The dataset did not contain missing values.

---

## 5. Find Genes With Zero Expression

Even though there were no missing values, some genes had a TPM value of 0 in every tissue.

I checked this using:

```python
expression_table = gene_tissue_table[tissue_cols]

zero_genes = (
    expression_table.sum(axis=1) == 0
)
```

This adds all tissue values for each gene.

If the total is 0, the gene has no expression in any tissue in the dataset.

These genes are removed because they do not give information about tissue expression and they cannot be normalized by dividing by their total expression.

---

## 6. Normalize Gene Expression

This is one of the most important parts of the analysis.

Some genes have very large TPM values while other genes have much smaller values.

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

The values are very different, but both genes have the same expression pattern.

If I only used the original TPM values, Gene A could have much more influence just because its numbers are larger.

To make the comparison more fair, I convert the expression into relative values.

For each gene I calculate:

```text
relative expression =
tissue TPM / total TPM for that gene
```

Example:

```text
Brain = 10
Liver = 30
Heart = 10
```

The total is:

```text
10 + 30 + 10 = 50
```

The normalized values become:

```text
Brain = 10 / 50 = 0.20
Liver = 30 / 50 = 0.60
Heart = 10 / 50 = 0.20
```

The values now add up to:

```text
1.00
```

This means that the normalized values show how the expression of one gene is distributed across the tissues.

I calculated this using:

```python
expression_table = gene_tissue_table[tissue_cols]

gene_total = expression_table.sum(axis=1)

normalized_expression = expression_table.div(
    gene_total,
    axis=0
)
```

---

## 7. Check the Normalization

After normalization, I checked that the tissue values for every gene add up to 1.

```python
normalized_expression.sum(axis=1).head()
```

The result should be close to:

```text
1.0
1.0
1.0
1.0
1.0
```

This check confirms that the normalization worked correctly.

---

## 8. Create the Normalized Gene-Tissue Table

After the normalization, I added the gene ID and gene name back to the table.

```python
normalized_gene_tissue_table = pd.concat(
    [
        gene_tissue_table[
            ["Name", "Description"]
        ].reset_index(drop=True),

        normalized_expression.reset_index(
            drop=True
        )
    ],
    axis=1
)
```

The final table looks like:

```text
Gene | Brain | Liver | Heart | ...
```

but now the tissue values represent relative expression instead of raw TPM.

---

# Tissue Score

The next part of my analysis is to calculate a tissue score for a group of genes.

The idea is to take a selected gene set and see which tissues have the strongest combined relative expression.

For example:

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

For brain:

```text
0.20 + 0.20 = 0.40
```

For liver:

```text
0.60 + 0.10 = 0.70
```

For heart:

```text
0.20 + 0.70 = 0.90
```

The result is:

```text
Heart = 0.90
Liver = 0.70
Brain = 0.40
```

This means that this group of genes has the strongest combined relative expression in the heart.

The code I use is:

```python
tissue_scores = selected_genes[
    tissue_cols
].sum(axis=0)
```

Then I sort the scores:

```python
tissue_scores = tissue_scores.sort_values(
    ascending=False
)
```

and create a table:

```python
tissue_score_table = (
    tissue_scores.reset_index()
)

tissue_score_table.columns = [
    "Tissue",
    "Score"
]
```

---

# Connection With the Research Papers

My analysis is not an exact reproduction of the Enrichr algorithm.

Instead, I am using the papers to understand the ideas behind gene-set enrichment and then applying similar ideas to GTEx tissue expression data.

The 2013 Enrichr paper shows that a list of genes can be analyzed together and compared with biological gene-set libraries.

It also explains why statistical testing is needed to determine if an observed biological signal is stronger than what would normally happen by chance.

The 2016 paper is especially useful because it includes GTEx gene-expression libraries and introduces fuzzy gene sets, where membership can have values between 0 and 1 instead of only 0 or 1.

This is related to my normalized gene-tissue table because each gene has a relative expression value between 0 and 1 for each tissue.

---

# Next Steps

The next steps of the analysis are:

```text
Normalized gene-tissue table
        ↓
Select a gene set
        ↓
Calculate tissue scores
        ↓
Rank tissues
        ↓
Generate random gene sets
        ↓
Calculate random tissue scores
        ↓
Compare real scores with random scores
        ↓
Calculate statistical significance
```

The goal of the randomization step is to determine whether a high tissue score represents a real signal or something that could happen by chance.

---

# Current Progress

Completed:

- Load GTEx data
- Identify gene and tissue columns
- Create gene-tissue table
- Check missing values
- Find genes with zero expression
- Remove all-zero genes
- Normalize gene expression
- Check that normalized values sum to 1
- Create normalized gene-tissue table
- Calculate tissue scores
- Create ranked tissue-score table

Still working on:

- selecting meaningful gene sets
- random gene-set testing
- statistical significance
- additional enrichment analysis
- comparing the results with methods described in the research papers

---

# References

Chen, E. Y., Tan, C. M., Kou, Y., et al. (2013).  
**Enrichr: interactive and collaborative HTML5 gene list enrichment analysis tool.**  
BMC Bioinformatics, 14:128.

Kuleshov, M. V., Jones, M. R., Rouillard, A. D., et al. (2016).  
**Enrichr: a comprehensive gene set enrichment analysis web server 2016 update.**  
Nucleic Acids Research, 44, W90-W97.