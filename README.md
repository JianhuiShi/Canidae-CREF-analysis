Copyright © 2024, Jianhui Shi & Lei M. Li. Academy of Mathematics and Systems Science, Chinese Academy of Sciences, Beijing 100190, China

# Canidae-CREF-analysis

## 1. Introduction

Canidae-CREF-analysis is a comprehensive and comparative analysis of genome-wide cis-regulatory element frequency (CREF) for five canids: dog, dingo, red fox, dhole, and wolf.

### 1.1 Pipeline description

1. Input the species-specific CREF matrices of five canids.
2. Perform the robust singular value decomposition (SVD) to stratify each CREF matrix into multiple dual eigen-modules at frequency levels from high to low. Each module is comprised of the singular value and the pair of gene- and motif-eigenvectors.
3. Polarize the gene- and motif-eigenvectors by sorting their loadings.
4. Evaluate the correlation between the motif-eigenvectors.
5. Analyze the rotation between motif-eigenvectors by projecting the fourth and fifth motif-eigenvectors of the other four canids onto the 2-D eigenspaces of dogs.
6. Compute the rank of all motifs according to their loadings in polarized motif-eigenvectors.
7. Perform the enrichment analysis on the polarized gene-eigenvectors to identify the biological pathways significantly enriched at the two poles.

### 1.2 Major results

1. The top three eigen-modules are highly conserved while a phase transition occurred between the fourth and fifth ones in dogs.

2. The dogs underwent a rotation in its motif-eigenvectors with reference to the wolf. The motif-eigenvectors of dingo and dhole are highly correlated with that of the dog.

3. Long-term memory, myelination, and cochlear development are significantly enhanced at level four in dogs.

4. The red fox is closest to the singularity or fusion point that characterizes the phase transition onset.

**Note**: This pipeline provides a systematic framework for analyzing the genome-wide CREF and is applicable to the comparison between other species.

## 2. Install

You can download the codes and data by the command:
```bash
git clone https://github.com/JianhuiShi/Canidae-CREF-analysis.git
```
or click the **Download Zip** button and decompress the package.

## 3. Dependencies

1. `R` with the following packages: `tidyverse`, `rlist`, `openxlsx`, and `snowfall`.

2. `MATLAB`

## 4. Usage

```bash
cd code
bash run.sh
```
All computational results can be found in the `results` directory generated.

## 5. Directory description

- ### `data` directory

  * **`CREF_matrix`** directory
  contains the CREF matrices of five canids: dog, dingo, red fox, dhole, and wolf.

  * **`geneName`** directory
  contains the gene names files of five canids. Gene names represent the row names of CREF matrix.

  * **`motif_symbol.txt`**
  contains all symbols of 1403 motifs used in this study. Motif symbols represent the column names of CREF matirx.

  * **`pathway`** directory
  contains the pathway data of [GO](https://geneontology.org/) (GO.BP, GO.CC, and GO.MF), [KEGG](https://www.kegg.jp/), and [Reactome](https://reactome.org/).

- ### `code` directory

  * **`run.sh`**
  This is the script used for running other scripts all in sequence.

  * **`config.yaml`**
  This is the configuration file that contains the basic information of the species used in computation.

  * **`motifRank.R`**
  This script computes the rank of all motifs according to their loadings in polarized motif-eigenvectors. The results are saved in the file `results/motifRank.xlsx`.

  * **`PCC.R`**
  This script computes the Pearson correlation coefficients between motif-eigenvectors across speices. The results are saved in the file `results/PCC.xlsx`.

  * **`projection.R`**
  This script computes the projections of the 4th and 5th motif-eigenvectors of one species onto
the 2-D eigen space spanned by the 4th and 5th motif-eigenvectors of another species. The results are saved in the file `results/projection_4-5.xlsx`.

  * #### `SVD` directory

    * **`robustSVD.m`**
    This script performs the robust singular value decomposition (robust SVD) on the CREF matrices of five canids. The singular values, gene-eigenvectors, and motif-eigenvectors are saved in the directory `results/SVD`. The polarized eigenvectors are saved in the directory `results/loadings`.

    * **`scripts`** directory
    contains matlab functions we defined for robust SVD.

    * **`yamlmatlab`** directory
    is the [yamlmatlab](https://github.com/ewiger/yamlmatlab) package.

    * **`inexact_alm_rpca`** directory
    is the matlab package for the [inexact augmented Lagrange multipliers (IALM)](https://arxiv.org/abs/1009.5055) method.

  * #### `enrichmentAnalysis` directory

    * **`gene_enrichment_analysis.R`**
    This script performs the gene enrichment analysis by the Wilcoxon scoring method on the polarized gene eigenvectors. The results are saved in the directory `results/enrichment/RData`.

    * **`format_enrichment_result.R`**
    This script formats the raw enrichment results (.RData) in `results/enrichment/RData` into readable files (.xlsx). The results are saved in the directory `results/enrichment/xlsx`.

    * **`enrichment_utils.R`**
    This script contains functions we defined for enrichment analysis.

## 6. Contact

Please contact shi.jianhui@foxmail.com for any questions.

## 7. License

**GNU General Public License v3.0 only**

For details, please read `Canidae-CREF-analysis/LICENSE`.
