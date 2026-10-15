# Targeted DNA damage sequencing and multi-omic analysis reveal repair kinetics and genomic context shaping UV-induced mutagenesis

This repository contains the analysis code used for the manuscript:

**Targeted DNA damage sequencing and multi-omic analysis reveal repair kinetics and genomic context shaping UV-induced mutagenesis**

The scripts are provided as used for the analyses in the paper. Shell scripts were run on the University of Georgia Sapelo2 high-performance computing cluster (SLURM), and file paths reflect our original project directories.

## Repository structure


## Data_analysis

| Script | Description |
|---|---|
| `Dinu_seq_analysis.sh` | Dinucleotide sequence analysis at damage sites |
| `Each_dinucleotide_percentage.py` | Calculates the percentage of each dinucleotide |
| `SVNs_analysis.sh` | Single-nucleotide variant (SNV) analysis |
| `XR_analysis.sh` | XR-seq (excision repair) analysis |
| `RNA_analysis.sh` | RNA-seq (transcription) analysis |
| `Histone_analysis.sh` | Histone modification and CTCF analysis |
| `Chrom3D_analysis.sh` | 3D chromatin organization analysis |
| `build_all_combine.sh` | Combines all features into the final tables used for analysis and modeling |

## Figure scripts

| Folder | Script | Analysis |
|---|---|---|
| `Figure1_Sfig1` | `Pearson_Cor_2_repli.R` | Pearson correlation between biological replicates |
| | `Nucl_frequency.Rmd` | Nucleotide frequency |
| | `Dinu_frequency.Rmd` | Dinucleotide frequency |
| | `Dinu_count_barplot.R` | Dinucleotide count bar plots |
| `Figure2_Fig3` | `01.SNV_4group.Rmd` | SNV analysis across four gene groups |
| | `02.Dinucleotide_4_group.Rmd` | Dinucleotide analysis across four gene groups |
| | `127gene_cor_MELA_CPD.Rmd` | Correlation between CPD damage and melanoma mutations across 127 genes |
| | `MELA_CPD_group.Rmd` | Gene grouping by CPD damage and melanoma mutation |
| `Figure4-5_Sfig3567` | `03.XR_Overall_Plot.Rmd` | Repair (XR-seq) analysis |
| | `04.Overall_RNA.Rmd` | Transcription analysis |
| | `05.Overall_Histone_CTCF.Rmd` | Histone modification and CTCF analysis |
| | `06.Group4_distance.Rmd` | Distance analysis across four gene groups |
| | `Model_analysis.Rmd` | Modeling of genomic features associated with mutagenesis |

## Software

- **Linux:** bash, BEDTools, SLURM
- **Python:** Python 3
- **R:** R (v4.x) with dplyr, ggplot2, reshape2, and other packages loaded in each script

## Data availability

Sequencing data generated in this study are available at GEO under accession **GSE315195**.

## License

This code is released under the MIT License.
