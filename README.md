# GTEx V11 Gene Expression Analysis

This project explores human gene expression using data from the
Genotype-Tissue Expression (GTEx) project.

The dataset contains median gene expression values across different
human tissues. The goal of this project is to explore how gene
expression changes between tissues and use data science methods to
identify patterns in the data.

## Dataset

The dataset used in this project is:

`GTEx_Analysis_2025-08-22_v11_RNASeQCv2.4.3_gene_median_tpm.gct.gz`

This is a GTEx V11 gene-level median TPM dataset.

The dataset contains:

- 74,628 genes
- 68 tissue types
- Ensembl gene IDs
- Gene symbols
- Median TPM expression values for each tissue

Some of the tissues included in the dataset are:

- Brain
- Heart
- Liver
- Lung
- Pancreas
- Kidney
- Skin
- Skeletal Muscle
- Thyroid
- Whole Blood

## What is TPM?

TPM stands for **Transcripts Per Million**.

TPM is a normalized measurement used to represent gene expression.
A higher TPM value usually means that a gene has higher expression
in that tissue.

For example:

| Gene | Brain | Liver | Muscle |
|------|------:|------:|-------:|
| Gene A | 20.5 | 2.1 | 0.5 |
| Gene B | 1.2 | 50.3 | 4.8 |

In this example, Gene A has higher expression in the brain while
Gene B has higher expression in the liver.

The values in this dataset are median TPM values, meaning that the
expression values from multiple samples of the same tissue were
combined using the median.

## Research Goal

The main goal of this project is to study how gene expression differs
between human tissues.

Some questions I want to explore are:

- Which genes have the highest expression across tissues?
- Which genes are expressed in almost every tissue?
- Which genes are specific to certain tissues?
- Which tissues have similar gene expression patterns?
- Which tissues are very different from each other?
- Can tissues be grouped based on their gene expression?
- Which genes contribute the most to differences between tissues?

## Tools

This project will use:

- Python
- VS Code
- Jupyter Notebook
- Pandas
- NumPy
- Matplotlib
- Scikit-learn
- Git
- GitHub

