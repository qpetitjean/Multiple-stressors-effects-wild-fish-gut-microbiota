# Script Title: Conduct LINDA Analysis for Differential Abundance Testing of Taxa 
#               at the Phylum Level in Water Samples using LMM
#
# Author: Quentin PETITJEAN
# Date Created: 03/2023
# Last Modified: 05/05/2024
# ==============================================================================
# Requirements:
# - R version 4.2.3
# - Packages:
#   - metabaR v1.0.0: For processing and subsetting metabarcoding datasets.
#   - phyloseq v1.42.0: For converting metabaR objects to phyloseq objects.
#   - MicrobiomeStat v1.1: For conducting LINDA analysis using linear mixed models.
#   - reshape2 v1.4.4: For reshaping and aggregating data.
#   - ggplot2 v3.5.1: For generating volcano and relative abundance plots.
#   - kableExtra v1.3.4: For creating enhanced HTML summary tables.
# - Source Files:
#   - volcanoPlot.R: Custom function to generate volcano plots from LINDA analysis results.
# ==============================================================================
# Script Overview:
# This script performs differential abundance analysis at the phylum level for water samples 
# using the LINDA method as described Zhou et al., 2022. The analysis tests the effects 
# of treatments (e.g., Contamination, Immune Challenge, and Origin) on the relative 
# abundance of taxa using a linear mixed model (LMM) framework.
#
# The workflow includes:
# 1. Installing and loading required packages.
# 2. Specifying directory paths for input and output.
# 3. Importing custom functions (volcanoPlot.R).
# 4. Importing and subsetting a cleaned metabaR object for water samples.
# 5. Retrieving OTU taxonomy and aggregating read counts at the phylum level.
# 6. Converting the aggregated data into a phyloseq object.
# 7. Merging sample metadata and releveling treatment factors.
# 8. Filtering out taxa detected in less than 5 samples.
# 9. Running LINDA analysis using an LMM approach via the MicrobiomeStat package.
# 10. Visualizing LINDA results with volcano plots and storing significant results in an HTML table.
# 11. Generating relative abundance bar plots for both top abundant taxa and for taxa with 
#     significant changes.
#
# Usage:
# 1. Update the 'savingDir' variable to the correct directory containing your input data 
#    (e.g., DesignData.csv, fguts_Bact_agg_MergedRep.RDS) and desired output location.
# 2. Ensure all required packages and source files are installed and accessible.
# 3. Run the script in an R environment to perform the LINDA analysis, generate volcano plots,
#    and create relative abundance plots.
# 4. The script outputs include:
#    - Volcano plots (TIF format) summarizing differential abundance results.
#    - Relative abundance plots (SVG format) at the phylum level.
#    - An HTML table summarizing significant LINDA analysis results.
# ==============================================================================
# References:
# Zhou H, He K, Chen J, Zhang X. LinDA: linear models for differential abundance analysis of microbiome compositional data. Genome Biol 2022; 23: 1–23. 
# ==============================================================================

#######################
# specify some path   #
#######################
# the path to the directory from where file should be both imported 
savingDir <- "W:/POSTDOC_INP_GOLFECH_2023/Outputs"

# directory where the output files should be saved
DirSave <- file.path(savingDir, "Normalized_Data")

##############################################
#   Import some custom functions             #
##############################################
source("https://raw.githubusercontent.com/qpetitjean/Multiple-stressors-effects-wild-fish-gut-microbiota/main/R_Func/volcanoPlot.R") # import a function used to display a vulcano plot from the result of Linda analysis

#######################
# Import the data     #
#######################

# import the cleaned dataset (metabaR object) and select only gut samples
labWater <- readRDS(file.path(savingDir, "Data/CleanedData", "fguts_Bact_agg_MergedRep.RDS"))
labWater <- metabaR::subset_metabarlist(labWater,
                                        table = "pcrs",
                                        indices = labWater$pcrs$matrix == "water")

names(labWater$samples)[which(names(labWater$samples) == "treatment_contam")] <- "Contam"
names(labWater$samples)[which(names(labWater$samples) == "pop")] <- "Pop"
names(labWater$samples)[which(names(labWater$samples) == "treatment_inj")] <- "Inj"

# retrieve OTUs taxonomy
OTUsTax <- data.frame(MOTU = colnames(labWater$reads),
                      phylum = labWater$motus$phylum, 
                      class = labWater$motus$class,
                      order = labWater$motus$order, 
                      family = labWater$motus$family, 
                      genus = labWater$motus$genus,
                      species = labWater$motus$species)

# aggregating the reads to the specified taxonomic level - Phylum
taxLevel <- "phylum"
labWaterFam <- metabaR::aggregate_motus(labWater, groups = labWater$motus[[taxLevel]])

# Convert to phyloseq object 
## makes otu table for phyloseq
read.phylo <-
  phyloseq::otu_table(t(labWaterFam[["reads"]]), taxa_are_rows = T)

## make sample table for phyloseq
sample.phylo <- phyloseq::sample_data(labWaterFam[["samples"]])

### merge phyloseq objects
PhyloDat  <- 
  phyloseq::phyloseq(read.phylo, sample.phylo)

##################################################################################################
# Test treatment effect on taxa (Phylum level) relative abundance
# using LINDA - https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02655-5 
#################################################################################################
PhyloDat@sam_data$ContPop <- ifelse(PhyloDat@sam_data$Pop %in% c("AUSCOR", "RIOU"), "C", "NC")

# test differential abundance of taxa levels among treatment using LMM approach
meta <- data.frame(Contam = factor(PhyloDat@sam_data$Contam),
                   Inj = factor(PhyloDat@sam_data$Inj),
                   Or = factor(PhyloDat@sam_data$ContPop),
                   Pop = factor(PhyloDat@sam_data$Pop))
rownames(meta) <- rownames(PhyloDat@sam_data)

# relevel the treatment
meta$Contam <- factor(meta$Contam, levels = c("NC", "C"))
meta$Inj <- factor(meta$Inj, levels = c("PBS", "AMIX"))
meta$Or <- factor(meta$Or, levels = c("NC", "C"))

toRemove <- which(apply(PhyloDat@otu_table, 1, function(x) length(which(x != 0 ))) <= 5)
if(length(toRemove) != 0){
  PhyloDat@otu_table <- PhyloDat@otu_table[-toRemove,]
}

linda.obj <- MicrobiomeStat::linda(feature.dat = as.data.frame(PhyloDat@otu_table), 
                                   meta.dat = meta,
                                   feature.dat.type = "count",
                                   formula = "~ Contam + Inj + Or", 
                                   alpha = 0.05, p.adj.method = "fdr",
                                   is.winsor = F,
                                   adaptive = T, 
                                   prev.filter = 0)

# check the results 
MicrobiomeStat::linda.plot(
  linda.obj,
  c('ContamC', "InjLPS", "OrC"),
  titles = c('Contam: NC vs C', 'Inj: PBS vs AMIX', 'Origin: NC vs C'),
  alpha = 0.05,
  lfc.cut = 0.5,
  legend = TRUE,
  directory = NULL,
  width = 11,
  height = 8
)

# store the significant results in a table
## keep only the significant results to store them in a table
SignifRes <- lapply(linda.obj$output, function(x) {
  x[which(x$log2FoldChange <= -0.5 &
            x$padj < 0.05 |
            x$log2FoldChange >= 0.5 &
            x$padj < 0.05),]
})

## remove the empty list element (non-signif. effects)
SignifRes <- Filter(function(x) nrow(x) > 0, SignifRes)

## rename variable
names(SignifRes)[names(SignifRes) == "ContamC"] <- "Contamination (NC vs. C)"
names(SignifRes)[names(SignifRes) == "InjAMIX"] <- "Imm. Chall. (PBS vs. AMIX)" 
names(SignifRes)[names(SignifRes) == "OrC"] <- "Origin (LP vs HP)"       

## round numeric value to 2 signbificant digits
SignifRes <- lapply(SignifRes, function(df) {
  numeric_columns <- sapply(df, is.numeric)  
  df[, numeric_columns] <- signif(df[, numeric_columns], digits = 2) 
  return(df)
})

## keep only desired columns
SignifRes <- lapply(SignifRes, function(x){
  x[c("baseMean", "log2FoldChange", "lfcSE", "stat", "df", "pvalue", "padj")]
})

## merge the list elements
table_data <- do.call(rbind, lapply(names(SignifRes), function(name) {
  df <- SignifRes[[name]]
  empty_row <- setNames(as.data.frame(matrix(NA, ncol = ncol(df), nrow = 1)), names(df))
  rbind(empty_row, df)  
}))
table_data <- table_data[-which(is.na(table_data$padj)), ]

## Create the HTML table using kable and kableExtra
html_table <-
  kableExtra::kbl(
    table_data,
    caption = "<span style='font-size:12px; font-weight: bold; font-style: italic'>Differential abundance of Taxa at the Phylum level - LINDA's Results </span>",
    format = "html",
    align = "c",
    escape = FALSE
  ) |>
  kableExtra::kable_classic_2(
    bootstrap_options = c("striped", "hover", "condensed", "responsive"),
    full_width = F,
    html_font = "Arial",
    font_size = 10
  )
start_row <- 1
for(name in names(SignifRes)) {
  end_row <- start_row + nrow(SignifRes[[name]]) - 1
  html_table <-
    kableExtra::pack_rows(html_table, name, start_row, end_row, label_row_css = "background-color: #fff; color: #000; border-bottom: 1px solid; border-top: 1px solid; font-style: italic;")
  start_row <-
    end_row + 1
}

html_string <- as.character(html_table)
modified_html_string <- gsub('(<td[^>]*>\\s*\\b[\\w\\s]+)\\d+(\\s*</td>)', '\\1\\2', html_string, perl = TRUE)

## save the table as html format
if (length(list.dirs(file.path(savingDir, "Normalized_Data"))) == 0) {
  dir.create(file.path(savingDir, "Normalized_Data"))
}
DirSave <- file.path(savingDir, "Normalized_Data")
if (length(list.dirs(file.path(DirSave, "Signif_Effects"))) == 0) {
  dir.create(file.path(DirSave, "Signif_Effects"))
}
if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots"))) == 0) {
  dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots"))
} 
if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Water"))) == 0) {
  dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Water"))
} 
if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Water", "TaxoLindaWater"))) == 0) {
  dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Water", "TaxoLindaWater"))
} 

kableExtra::save_kable(modified_html_string, 
                       file.path(DirSave, "Signif_Effects", "Signif_Plots", "Water", "TaxoLindaWater", 
                                 paste0("TaxPhylumLindaWaterTabSignif", ".HTML")))

# specify p value and log2 fold change threshold
FcTresh <- c(-0.5,0.5)
pTresh <- -log10(0.05)

# specify the colors to represent significance thresholds
colors <- c(adjustcolor("#808080", alpha.f = 0.5),
            adjustcolor("#420693", alpha.f = 0.5))
names(colors) <- c("NS", "p-adj. & Log2FC")

# specify a color vecotr to display each taxa 
colorsTax <- c(
  "#4F7387", # Darker Light Blue
  "#6FA176", # Darker Light Green
  "#D35F5F", # Richer Dark Pink
  "#8E6E8D", # Darker Lavender
  "#E57342", # Vibrant Peach
  "#B8A832", # Darker Yellow-Green
  "#8A5A3D", # Richer Beige
  "#3F8667", # Darker Forest Green
  "#C3708E", # Darker Pastel Pink
  "#E4B3B3", # Warm Beige (replacing Darker Soft White)
  "#4A76A8", # Darker Mid-Blue
  "#B76082", # Darker Rose
  "#88A865", # Darker Lime Green
  "#D1A84D", # Vibrant Light Cream
  "#A66DAA", # Darker Lavender
  "#607FA1", # Darker Sky Blue
  "#73A679", # Darker Olive Green
  "#6B6196", # Deep Lilac
  "#C69665", # Darker Apricot
  "#8A5F9A", # Darker Purple
  "#B3727D", # Darker Coral Pink
  "#D05078"  # Deep Coral Red
)

# Customized volcano Plot
## for contamination 
i = "ContamC"
linda.obj$output[[i]]$VolcanoGroups <- rep(NA, nrow(linda.obj$output[[i]]))

# specify non significant OTU
linda.obj$output[[i]]$VolcanoGroups[which(
  linda.obj$output[[i]]$log2FoldChange > FcTresh[1] &
    linda.obj$output[[i]]$log2FoldChange < FcTresh[2] &
    -log10(linda.obj$output[[i]]$padj) < pTresh
)] <- "NS"

# specify OTU above pvalue
linda.obj$output[[i]]$VolcanoGroups[which(
  -log10(linda.obj$output[[i]]$padj) > pTresh
)] <- "p-adj."

# specify OTU above Log2 FC
linda.obj$output[[i]]$VolcanoGroups[which(
  linda.obj$output[[i]]$log2FoldChange <= FcTresh[1] |
    linda.obj$output[[i]]$log2FoldChange >= FcTresh[2]
)] <- "Log2FC"

# specify significant OTU
linda.obj$output[[i]]$VolcanoGroups[which(
  linda.obj$output[[i]]$log2FoldChange <= FcTresh[1] &
    -log10(linda.obj$output[[i]]$padj) > pTresh |
    linda.obj$output[[i]]$log2FoldChange >= FcTresh[2] &
    -log10(linda.obj$output[[i]]$padj) > pTresh
)] <- "p-adj. & Log2FC"

# append specified colors to significance thresholds
linda.obj$output[[i]]$colors = colors[as.character(linda.obj$output[[i]]$VolcanoGroups)]  
linda.obj$output[[i]]$colors[which(is.na(linda.obj$output[[i]]$colors))] <- colors[["NS"]]

# instead replace the color of significantly less or over abundant taxa (p-adj. & Log2FC) by selecting 1 color per taxa
taxL <- length(linda.obj$output[[i]]$colors[which(linda.obj$output[[i]]$VolcanoGroups == "p-adj. & Log2FC")])
linda.obj$output[[i]]$colors[which(linda.obj$output[[i]]$VolcanoGroups == "p-adj. & Log2FC")] <- sample(colorsTax, taxL, replace = FALSE)

# retrieve total number of reads per tax level and rescale them to specify dot size in the plot
linda.obj$output[[i]][["readsCount"]] <- apply(PhyloDat@otu_table, 1, function(x) sum(x, na.rm = T))
DotSize <- c(1, 6)
maxVal <- max(linda.obj$output[[i]][["readsCount"]])
minVal <- min(linda.obj$output[[i]][["readsCount"]])
normVal <- (linda.obj$output[[i]][["readsCount"]] - 
              minVal) / (maxVal - minVal)
cexPTS <- DotSize[1] + normVal * (DotSize[2] - DotSize[1])

# draw the volcano plot
if (length(list.dirs(file.path(savingDir, "Normalized_Data"))) == 0) {
  dir.create(file.path(savingDir, "Normalized_Data"))
}
DirSave <- file.path(savingDir, "Normalized_Data")
if (length(list.dirs(file.path(DirSave, "Signif_Effects"))) == 0) {
  dir.create(file.path(DirSave, "Signif_Effects"))
}
if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots"))) == 0) {
  dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots"))
} 
if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Water"))) == 0) {
  dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Water"))
} 
if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Water", "TaxoLindaWater"))) == 0) {
  dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Water", "TaxoLindaWater"))
} 
toSave <- file.path(savingDir, "Normalized_Data", "Signif_Effects", "Signif_Plots", "Water", "TaxoLindaWater", "Phylum_Water_Contam.tif")
if(file.exists(toSave)){
  unlink(toSave)
}
grDevices::tiff(filename = toSave,
                width = 1800,
                height = 1200,
                res = 180,
                compression = "lzw",
                pointsize = 14)
volcanoPlot(
  x = linda.obj$output[[i]]$log2FoldChange,
  y = -log10(linda.obj$output[[i]]$padj),
  group = linda.obj$output[[i]]$VolcanoGroups,
  col = linda.obj$output[[i]]$colors,
  labels = rownames(linda.obj$output[[i]]),
  selLabels = rownames(linda.obj$output[[i]][which(
    linda.obj$output[[i]]$VolcanoGroups != "NS" &
      linda.obj$output[[i]]$VolcanoGroups != "Log2FC" &
      linda.obj$output[[i]]$VolcanoGroups != "p-adj." 
  ), ]),
  ylim = c(-5, 25),
  xlim = c(-10, 10),
  ylab = bquote( ~ -Log[10] ~ 'p-adj.'),
  xlab = bquote( ~ Log[2] ~ 'fold change'),
  vlines = FcTresh,
  hlines = pTresh,
  cex.lab = 1.2,
  leg.order = c("NS", "p-adj. & Log2FC"),
  cex.pts = cexPTS
)
# add legend for dot size (number of reads)
leg <- c(min(cexPTS), max(cexPTS)/2, max(cexPTS))
legValues <- ((leg - DotSize[1]) / (DotSize[2] - DotSize[1])) * 
  (maxVal - minVal) + minVal
text(0.45, 24.65,"# of reads:", cex = 1.0)
legend(
  x = 1.15,
  y = 25.85,
  legend = round(legValues),
  horiz = T,
  pch = 16,
  col = adjustcolor("#000000", alpha.f = 0.8),
  pt.bg =  adjustcolor("#000000", alpha.f = 0.8),
  pt.cex = leg,
  #y.intersp = c(1, 2, 4),
  x.intersp = c(1, 1, 2),
  text.width = c(1, 1.5, 2),
  bty = "n"
)
dev.off()

# draw relative abundance plot 
## for contamination and for the 15 top families (there is only 13 here)
# Aggregate per treatment
NFunc <- 15 # number of level to display

tmp <-
  reshape2::melt(stats::aggregate(t(PhyloDat@otu_table), by = list(
    meta$Contam), sum))

# Sorting and selecting top Ntax abundant MOTUS
summedDat <- stats::aggregate(value ~ variable, data = tmp, FUN = sum)
Top <- summedDat[order(summedDat$value, decreasing = T), ][1:NFunc, ]

# keep only the top Ntax abundant MOTUS
tmpTop <- tmp[which(tmp[["variable"]] %in% Top[["variable"]]),]

# Splitting the dataframe
tmp_C <- tmpTop[tmpTop$Group.1 == "C", ]
tmp_NC <- tmpTop[tmpTop$Group.1 == "NC", ]

# compute relative abundances
tmp_C$Relvalue <- tmp_C$value/sum(tmp_C$value)
tmp_NC$Relvalue <- tmp_NC$value/sum(tmp_NC$value)

# Merging the dataframes
tmp_merged <- rbind(tmp_C, tmp_NC)

# Plot 
newOrder <- c("NC", "C")
tmp_merged[["ContamOrd"]] <- factor(tmp_merged[["Group.1"]], levels = newOrder)

nb.cols <- NFunc
RelaAB_Plot_Contam_Top <-
  ggplot2::ggplot(tmp_merged,
                  ggplot2::aes(x = ContamOrd, y = Relvalue, fill = variable)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::labs(x = NULL, y = "Relative abundance", fill = "Taxa")  + ggplot2::theme_bw() +
  ggplot2::scale_fill_manual(values = colorsTax) +
  ggplot2::theme_classic() +
  ggplot2::theme(text = ggplot2::element_text(size=20)) +
  ggplot2::ggtitle("")

# save the plot
toSave <- file.path(savingDir, "Normalized_Data", "Signif_Effects", "Signif_Plots", 
                    "Water", "TaxoLindaWater", "EC_Contam_Water_RelAb_Phylum_TOP.svg")
if(file.exists(toSave)){
  unlink(toSave)
} 
ggplot2::ggsave(
  file = toSave,
  plot = RelaAB_Plot_Contam_Top,
  width = 9,
  height = 6,
  dpi = 300
)

## for contamination and for the taxa with significant changes in relative abundance 
# Aggregate per treatment
tmp <-
  reshape2::melt(stats::aggregate(t(PhyloDat@otu_table), by = list(
    meta$Contam), sum))

# Splitting the dataframe
tmp_C <- tmp[tmp$Group.1 == "C", ]
tmp_NC <- tmp[tmp$Group.1 == "NC", ]

# compute relative abundances
tmp_C$Relvalue <- tmp_C$value/sum(tmp_C$value)
tmp_NC$Relvalue <- tmp_NC$value/sum(tmp_NC$value)

# Merging the dataframes
tmp_merged <- rbind(tmp_C, tmp_NC)

# Selecting only taxa with significantly altered relative abundance 
taxSignif <- linda.obj$output[[i]][which(linda.obj$output[[i]]$VolcanoGroups == "p-adj. & Log2FC"),]
taxSignif$taxa <- rownames(taxSignif)
tmp_merged_signif <- tmp_merged[which(tmp_merged$variable %in% taxSignif$taxa),]
tmp_merged_Nonsignif <- tmp_merged[which(!tmp_merged$variable %in% taxSignif$taxa),]
tmp_merged_Nonsignif <- stats::aggregate(
  . ~ Group.1,
  data = tmp_merged_Nonsignif[, c("Group.1", "value", "Relvalue")],
  sum,
  na.rm = TRUE
)
tmp_merged_Nonsignif$variable <- "Background Taxa"
tmp_merged_Nonsignif <- tmp_merged_Nonsignif[, c("Group.1", "variable", "value", "Relvalue")]
tmp_merged <- rbind(tmp_merged_signif, tmp_merged_Nonsignif)
tmp_merged <- merge(tmp_merged,
                    taxSignif[, c("taxa", "colors")],
                    by.x = "variable",
                    by.y = "taxa",
                    all.x = T)
tmp_merged$colors[which(is.na(tmp_merged$colors))] <- adjustcolor("#808080", alpha.f = 0.5)

# Plot 
newOrder <- c("NC", "C")
tmp_merged[["ContamOrd"]] <- factor(tmp_merged[["Group.1"]], levels = newOrder)

RelaAB_Plot_Contam_signif <-
  ggplot2::ggplot(tmp_merged,
                  ggplot2::aes(x = ContamOrd, y = Relvalue, fill = variable)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::labs(x = NULL, y = "Relative abundance", fill = "Taxa")  + ggplot2::theme_bw() +
  ggplot2::scale_fill_manual(values = unique(tmp_merged$colors)) +
  ggplot2::theme_classic() +
  ggplot2::theme(text = ggplot2::element_text(size=20)) +
  ggplot2::ggtitle("")

# save the plot
toSave <- file.path(savingDir, "Normalized_Data", "Signif_Effects", "Signif_Plots", 
                    "Water", "TaxoLindaWater", "EC_Contam_Water_RelAb_Phylum_Signif.svg")
if(file.exists(toSave)){
  unlink(toSave)
} 
ggplot2::ggsave(
  file = toSave,
  plot = RelaAB_Plot_Contam_signif,
  width = 9,
  height = 6,
  dpi = 300
)
