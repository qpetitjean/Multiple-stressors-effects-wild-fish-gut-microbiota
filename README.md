# Multiples stressors effects on gut microbiota composition in wild fish populations.
## General
![GitHub](https://img.shields.io/github/license/qpetitjean/Multiple-stressors-effects-wild-fish-gut-microbiota)
![GitHub release (latest by date)](https://img.shields.io/github/v/release/qpetitjean/Multiple-stressors-effects-wild-fish-gut-microbiota)
![GitHub last commit](https://img.shields.io/github/last-commit/qpetitjean/Multiple-stressors-effects-wild-fish-gut-microbiota)
![GitHub repo size](https://img.shields.io/github/repo-size/qpetitjean/Multiple-stressors-effects-wild-fish-gut-microbiota)

Datasets: <br /><a href="https://doi.org/10.5281/zenodo.14989875"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.14989875.svg" alt="DOI"></a><br />

Permanent copy of this code repository: <br /><a href="https://doi.org/10.5281/zenodo.14990007"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.14990007.svg" alt="DOI"></a><br />

This repository contains:
<ul>
<li>Some R useful R function to allow some computations, model selection and reproduce the visualization of the results within the <code>R_Func</code> directory</li>
<li>The R script used to preprocess, clean, and analyse the data within the <code>Statistical_analyses</code> directory</li>
</ul><br />
Files within this code repository, Datas from [Zenodo repository](https://doi.org/10.5281/zenodo.14989875) and the step by step procedure detailed below allow to fully reproduce the results of the paper by testing the effects of multiples stressors (i.e., metal contamination and immune challenge) on the gut microbiota of wild fish populations (n=5, <i>Gobio Occitaniae</i>). These results are discussed in the following manuscript: <br /><br />

>Petitjean, Q., Jean, S., Granada, M., Manzi, S., Veyssiere, C., Perrault, A., Cousseau, M., Laffaille, P., Jacquin, L., White, J., 2025. Experimental metal contamination reduces gut microbiota diversity and alters its composition and function in wild-caught fish. https://doi.org/10.1101/2025.03.21.644596 <br /><br />
Preprint is available here: <br /> 
https://www.biorxiv.org/content/10.1101/2025.03.21.644596v1
<br />

## Step by step procedure to reproduce the analyses:

<u>NB:</u> Before skimming trough the following scripts, check the list of R packages to install below. 

### Data preparation
<ol start="1">
<li>Download the raw and cleaned dataset &ast; from the Zenodo repository (<code>Data</code> directory) </li>
https://doi.org/10.5281/zenodo.14989874 <br />

<li>Clean the dataset using <code>Statistical_analyses/Data_pre-processing_MetabaR/Stats_MetabaRLab_MergedRep.R</code> <br />
<u>NB:</u> Optional, the cleaned dataset is available in <code>Data/CleanedData</code> as <code>fguts_Bact_agg_MergedRep.RDS</code>) </li><br />

<li>Construct the phylogenetic tree using <code>Statistical_analyses/Phylogenetic_Tree/PhyloTree.R</code> <br />
<u>NB:</u> Optional, the phylogenetic tree is available in <code>Data/PhyloTree</code></li><br />
</ol>

### Analyses conducted on MOTUs
<ol start="4">
<li>Compute and test (LMM) alpha diversity metrics on:
<ul>
 <li>Fish gut microbiota using <code>Statistical_analyses/Alpha_Diversity/AlphaDivLMM_MOTUs.R</code></li>
 <li>Water microbiota using <code>Statistical_analyses/Alpha_Diversity/AlphaDivLMM_Water.R</code> </li>
  </li>
 </ul><br />
 
<li>Compute and test (LMM & envfit) beta diversity metrics on:
<ul>
 <li>Fish gut microbiota using <code>Statistical_analyses/Beta_Diversity/BetaDivLMM_MOTUs.R</code></li>
 <li>Water microbiota using <code>Statistical_analyses/Beta_Diversity/BetaDivLM_Water.R</code></li>
  </li>
 </ul><br />

<li>Compute and test (betadisper) variance homogeneity of beta diversity metrics on fish gut microbiota using <code>Statistical_analyses/Beta_Diversity/BetaDisper_Tax.R</code> </li><br />

<li>Check whether results are consistent across normalization, ordination methods, and indices using <code>Statistical_analyses/Check_Normalization/ExtractNormRes.R</code> <br />
<u>NB:</u> Optional, this step is not mandatory to perform the next steps</li><br />

<li>Test treatments and covariates effects on Taxonomic differential abundance using Linda (LMM) on:
<ul>
 <li>Fish gut microbiota at:</li> 
 <ul>
  <li>The family level using <code>Statistical_analyses/Taxonomic_analysis/LINDA_TaxComp_Family.R</code></li> 
  <li>The phylum level using <code>Statistical_analyses/Taxonomic_analysis/LINDA_TaxComp_Phylum.R</code></li>
 </ul>
 <li>Water microbiota at:</li> 
 <ul>
  <li>The family level using <code>Statistical_analyses/Taxonomic_analysis/LINDA_TaxComp_Family_water.R</code></li> 
  <li>The phylum level using <code>Statistical_analyses/Taxonomic_analysis/LINDA_TaxComp_Phylum_water.R</code></li>
 </ul>
  </li>
 </ul><br />

<li>Compute and test common taxonomic levels ratios indicating dysbiosis such as Bacteroidota/Proteobacteria and Firmicutes/Bacteroidota ratio using <code>Statistical_analyses/Taxonomic_analysis/Dysbiosis_Ratio.R</code>  </li> 
</ol><br />

### Analyses conducted on Functions
<ol start="10">
<li>Perform functional inferences:  </li>
<ul>
  <li>Prepare the file needed to convert the dataset to Biom format using <code>Statistical_analyses/Functional_analysis/ConvertToBiom.R</code><br />
<u>NB:</u> Optional, the dataset converted to BIOM format is available in <code>Data/FunctionalInferences/BiomOutput</code> </li>
  
  <li>Perform functional inferences (on Linux distribution) using the annotated code in <code>Statistical_analyses/Functional_analysis/FuncInferences_Picrust2.txt</code><br />
<u>NB:</u> Optional, the results of the functional inferences are available in <code>Data/FunctionalInferences/Picrust2Output</code></li>
 </ul><br />

<li>Compute and test (LMM) alpha diversity metrics on inferred functions using <code>Statistical_analyses/Alpha_Diversity/AlphaDivLMM_Function.R</code> </li><br />

<li>Compute and test (LMM) beta diversity metrics on inferred functions using <code>Statistical_analyses/Beta_Diversity/BetaDivLMM_Functions.R</code> </li><br />

<li>Compute and test (betadisper) variance homogeneity of beta diversity metrics on inferred functions using <code>Statistical_analyses/Beta_Diversity/BetaDisper_Functions.R</code> </li><br />

<li>Test treatments and covariates effects on functions differential abundance using Linda (LMM) using <code>Statistical_analyses/Functional_analysis/LINDA_FuncComp.R</code> </li><br />
</ol>

&ast; The raw dataset available in the [Zenodo repository](https://doi.org/10.5281/zenodo.14989875) has already been processed trough a bioinformatic pipeline including pair-end assembly, demultiplexing, filtration of short, poorly aligned, duplicated and chimeric sequences as well as sequences containing non attributed nucleotides. Also, as sequencing has been conducted on triplicate samples, triplicates were merged based on a 97% similarity treshold and taxonomic annotation was performed using the [SILVAngs database](https://www.arb-silva.de/) for small sequences (16S and 18S; v138.1 released on August 27th, 2020).

## List of R packages needed 

```{r} 
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
if (!require("remotes", quietly = TRUE))
    install.packages("remotes")
if (!require("ape", quietly = TRUE))
    install.packages("ape")
if (!require("car", quietly = TRUE))
    install.packages("car")
if (!require("cowplot", quietly = TRUE))
    install.packages("cowplot")
if (!require("basicPlotteR", quietly = TRUE))
    devtools::install_github("JosephCrispell/basicPlotteR")
if (!require("biomformat", quietly = TRUE))
    BiocManager::install("biomformat")
if (!require("Biostrings", quietly = TRUE))
    BiocManager::install("Biostrings")
if (!require("DECIPHER", quietly = TRUE))
    BiocManager::install("DECIPHER")
if (!require("Deseq2", quietly = TRUE))
    BiocManager::install("Deseq2")
if (!require("edgeR", quietly = TRUE))
    BiocManager::install("edgeR")
if (!require("ggplot2", quietly = TRUE))
    install.packages("ggplot2")
if (!require("ggpubr", quietly = TRUE))
    install.packages("ggpubr")
if (!require("ggtree", quietly = TRUE))
    BiocManager::install("ggtree")
if (!require("graphics", quietly = TRUE))
    install.packages("graphics")
if (!require("grdevices", quietly = TRUE))
    install.packages("grdevices")
if (!require("kableExtra", quietly = TRUE))
    install.packages("kableExtra")
if (!require("lme4", quietly = TRUE))
    install.packages("lme4")
if (!require("lme4", quietly = TRUE))
    install.packages("lme4")
if (!require("metabaR", quietly = TRUE))
    remotes::install_github("metabaRfactory/metabaR")
if (!require("metagenomeSeq", quietly = TRUE))
    BiocManager::install("metagenomeSeq")
if (!require("MicrobiomeStat", quietly = TRUE))
    devtools::install_github("cafferychen777/MicrobiomeStat")
if (!require("MuMIn", quietly = TRUE))
    install.packages("MuMIn")
if (!require("performance", quietly = TRUE))
    install.packages("performance")
if (!require("phangorn", quietly = TRUE))
    install.packages("phangorn")
if (!require("phyloseq", quietly = TRUE))
    BiocManager::install("phyloseq")
if (!require("picante", quietly = TRUE))
    install.packages("picante")
if (!require("plyr", quietly = TRUE))
    install.packages("plyr")
if (!require("progress", quietly = TRUE))
    install.packages("progress")
if (!require("RcolorBrewer", quietly = TRUE))
    install.packages("RcolorBrewer")
if (!require("reshape2", quietly = TRUE))
    install.packages("reshape2")
if (!require("seqinr", quietly = TRUE))
    install.packages("seqinr")
if (!require("stats", quietly = TRUE))
    install.packages("stats")
if (!require("vegan", quietly = TRUE))
    install.packages("vegan")
```

## Citation

Article:<br />Petitjean, Q., Jean, S., Granada, M., Manzi, S., Veyssiere, C., Perrault, A., Cousseau, M., Laffaille, P., Jacquin, L., White, J., 2025. Experimental metal contamination reduces gut microbiota diversity and alters its composition and function in wild-caught fish. https://doi.org/10.1101/2025.03.21.644596 <br />

Datasets:<br />PETITJEAN, Q. (2025). Dataset from: Petitjean et al. (Submitted) Experimental metal contamination reduces gut microbiota diversity and alters its composition and function in wild-caught fish. [Data set]. Zenodo. https://doi.org/10.5281/zenodo.14989875 <br />

Code:<br />Quentin PETITJEAN. (2025). Code from: Petitjean et al. (Submitted) Experimental metal contamination reduces gut microbiota diversity and alters its composition and function in wild-caught fish (V1.0). Zenodo. https://doi.org/10.5281/zenodo.14990007 <br />