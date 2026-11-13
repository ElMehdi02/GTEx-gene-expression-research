# GTEx Data Analysis Process

This file explains how we are analyzing the GTEx gene expression dataset and how the research papers we are reading are helping us build the analysis.

The main notebook for this work is:

`data_exploration1.ipynb`

The goal of this analysis is to understand how genes are connected to different human tissues using GTEx gene expression data.

---

## Research Papers We Are Using

We are using two Enrichr research papers to understand how gene-set enrichment analysis works and how gene expression can be connected to biological terms such as tissues.

### 1. Enrichr: Interactive and Collaborative HTML5 Gene List Enrichment Analysis Tool

Chen et al., 2013.

This paper explains how Enrichr takes a list of genes and compares it with existing gene-set libraries.

The paper uses different methods to measure enrichment, including:

- Fisher exact test
- expected-rank correction
- Z-score
- combined score

The main idea we are using from this paper is that a group of genes can be analyzed together to find biological signals that are more represented than expected.

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

This idea is useful for our project because we are also working with relative gene-tissue values between 0 and 1.

---

# Our Analysis Process

## 1. Load the GTEx Dataset

We started by loading the GTEx median TPM dataset using Pandas.

The file contains:

- gene IDs
- gene names
- tissue names
- median TPM values

Each row represents a gene and each tissue column contains the expression of that gene in that tissue.

The first two rows of the GCT file contain information about the file, so we skip them when loading the data.

## 2. Identify the Tissue Columns

The first two columns contain information about the genes:

```text
Name
Description
```

The remaining columns are the tissues.

We selected the tissue columns so we can work only with the expression values when performing calculations.

## 3. Create a Gene-Tissue Table

We created a table that contains:

```text
Gene ID | Gene Name | Tissue 1 | Tissue 2 | Tissue 3 | ...
```

Each value in the tissue columns is the median TPM expression for that gene.

This table is the main table we use for the analysis.

## 4. Check the Dataset

Before doing calculations, we checked the structure of the data.

We checked:

- number of genes
- number of tissues
- column names
- data types
- missing values

The dataset did not contain missing values.

## 5. Find Genes With Zero Expression

Even though there were no missing values, some genes had a TPM value of 0 in every tissue.

We checked for these genes by adding the tissue values for each gene.

If the total is 0, the gene has no expression in any tissue in the dataset.

These genes are removed because they do not give information about tissue expression and they cannot be normalized by dividing by their total expression.

## 6. Normalize Gene Expression

This is one of the most important parts of our analysis.

Some genes have very large TPM values while other genes have much smaller values.

If we only used the original TPM values, genes with larger numbers could have more influence just because their overall expression is higher.

To make the comparison more fair, we convert the expression into relative values.

For each gene we calculate:

```text
relative expression =
tissue TPM / total TPM for that gene
```

The normalized values for each gene add up to:

```text
1.00
```

This means that the normalized values show how the expression of one gene is distributed across the tissues.

## 7. Check the Normalization

After normalization, we checked that the tissue values for every gene add up to 1.

This confirms that the normalization worked correctly.

## 8. Create the Normalized Gene-Tissue Table

After the normalization, we added the gene ID and gene name back to the table.

The final table contains gene information together with the normalized relative-expression values.

---

# Tissue Score

The next part of our analysis is to calculate a tissue score for a group of genes.

The idea is to take a selected gene set and see which tissues have the strongest combined relative expression.

A higher tissue score means that the selected group of genes has a stronger combined signal in that tissue.

---

# Connection With the Research Papers

Our analysis is not an exact reproduction of the Enrichr algorithm.

Instead, we are using the papers to understand the ideas behind gene-set enrichment and applying similar ideas to GTEx tissue-expression data.

The 2013 Enrichr paper shows that a list of genes can be analyzed together and compared with biological gene-set libraries.

The 2016 paper is especially useful because it includes GTEx gene-expression libraries and introduces fuzzy gene sets, where membership can have values between 0 and 1 instead of only 0 or 1.

This is related to our normalized gene-tissue table because each gene can have different relative expression values for different tissues.

---

# Next Steps

The next steps of our analysis are:

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