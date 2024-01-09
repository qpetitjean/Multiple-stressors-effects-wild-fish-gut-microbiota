# Multiples stressors effects on gut microbiota composition in wild fish populations.
## General
![GitHub](https://img.shields.io/github/license/qpetitjean/Multiple-stressors-effects-wild-fish-gut-microbiota)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/qpetitjean/Multiple-stressors-effects-wild-fish-gut-microbiota)
![GitHub last commit](https://img.shields.io/github/last-commit/qpetitjean/Multiple-stressors-effects-wild-fish-gut-microbiota)
![GitHub repo size](https://img.shields.io/github/repo-size/qpetitjean/Multiple-stressors-effects-wild-fish-gut-microbiota)

This repository contains the scripts used to manage the data and test the effects of multiples stressors (i.e., metal contamination and immune challenge) on the gut microbiota of wild fish populations (n=5, _Gobio Occitaniae_). These results are discussed in the following manuscript: 

PETITJEAN, Q., GRANADA, M., JEAN, S., MANZI S., VEYSSIERE C., PERRAULT, A., COUSSEAU, M., LAFFAILLE P., WHITE J., JACQUIN, L., 2023. Environmentally relevant metal contamination levels affect gut microbiota composition in wild fish populations (_Gobio Occitaniae_).

## step by step procedure to reproduce the analyses:
### Data preparation
0- download the raw dataset\* from the figshare repository: 
1- clean the dataset using Stats_MetabaRLab_MergedRep.R
2- construct the phylogenetic tree using PhyloTree.R 

### Analyses conducted on MOTUs
3- compute and test (LMM) alpha diversity metrics on fish gut microbiota using AlphaDivLMM_MOTUs.R and water microbiota using AlphaDivLMM_Water.R
4- compute and test (LMM & envfit) beta diversity metrics on fish gut microbiota using BetaDivLMM_MOTUs.R and water microbiota using BetaDivLM_Water.R
5- compute and test (betadisper) variance homogeneity of beta diversity metrics on fish gut microbiota  using BetaDisper_Tax.R
6 (optional)- check whether results are consistent across normalization, ordination methods, and indices using ExtractNormRes.R
7- Test treatments and covariates effects on Taxonomic differential abundance using Linda (LMM) using LINDA_TaxComp_Family.R and LINDA_TaxComp_Phylum.R to perform analysis at the family and phylum level, respectively.


8- ratios....................

### Analyses conducted on Functions
<ins>9- perform functional inferences:</ins>
  9a- prepare the file needed to convert the dataset to Biom format using ConvertToBiom.R 
  9b- perform functional inferences (on Linux distribution) using the annotated code in FuncInferences_Picrust2.txt

10- compute and test (LMM) alpha diversity metrics using AlphaDivLMM_Function.R
11- compute and test (LMM) beta diversity metrics using BetaDivLMM_Functions.R

12 TODO- compute and test (betadisper) variance homogeneity of beta diversity metrics on fish gut microbiota using BetaDisper_Func.R
9f- Test treatments and covariates effects on Taxonomic differential abundance using Linda (LMM) using LINDA_FuncComp.R

10- Corr host traits................


\* The raw dataset available in the figshare repository has already been processed trough a bioinformatic pipeline including pair-end assembly, demultiplexing, filtration of short, poorly aligned, duplicated and chimeric sequences as well as sequences containing non attributed nucleotides. Also, as sequencing has been conducted on triplicate samples, triplicates were merged based on a 97% similarity treshold and taxonomic annotation was performed using the SILVAngs database for small sequences (16S and 18S; v138.1 released on August 27th, 2020).