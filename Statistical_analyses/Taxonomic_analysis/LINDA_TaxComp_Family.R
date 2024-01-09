##############################################################################################################
#  Conduct LINDA analysis (effects of envir. variable. on relative abundance among taxa with LMM)            #
#############################################################################################################
source("D:/POSTDOC_INP_GOLFECH_2023/R scripts/GitRepository/R_Func/volcanoPlot.R") # import a function used to import and filter cleaned dataset for further use

#######################
# specify some path   #
#######################
# the path to the directory from where file should be both imported 
savingDir <- "D:/POSTDOC_INP_GOLFECH_2023/Outputs"

# directory where the output files should be saved
DirSave <- "D:/POSTDOC_INP_GOLFECH_2023/Outputs/Normalized_Data"

###############################################################
# data preparation
###############################################################

# import the full dataset from the experiment (treatment & behavioral & physiological measures)
FullDat <-
  read.csv2(file.path(savingDir, "Varia_contam_fullv5.csv"), dec = ".", sep = ";")

# import the cleaned dataset (metabaR object) and select only gut samples
labFish <- readRDS(file.path(savingDir, "Preprocessing-Metabar/fguts_Bact_agg_MergedRep.RDS"))
labFish <- metabaR::subset_metabarlist(labFish,
                                       table = "pcrs",
                                       indices = labFish$pcrs$matrix == "fish_gut")

# retrieve OTUs taxonomy
OTUsTax <- data.frame(MOTU = colnames(labFish$reads),
                      phylum = labFish$motus$phylum, 
                      class = labFish$motus$class,
                      order = labFish$motus$order, 
                      family = labFish$motus$family, 
                      genus = labFish$motus$genus,
                      species = labFish$motus$species)

# aggregating the reads to the specified taxonomic level - family
taxLevel <- "family"
labFishFam <- metabaR::aggregate_motus(labFish, groups = labFish$motus[[taxLevel]])

# merge FullDat to the metabar sample object
labFishFam[["samples"]]["RowNames"] <- rownames(labFishFam[["samples"]])
labFishFam[["samples"]] <- merge(
  x = FullDat,
  y = labFishFam[["samples"]],
  by.x = "Ind",
  by.y = "Num_prlvt_Euth",
  all.y = TRUE
)
rownames(labFishFam[["samples"]]) <- labFishFam[["samples"]]$Ind
rownames(labFishFam[["reads"]]) <- labFishFam[["samples"]][["Ind"]][match(rownames(labFishFam[["reads"]]), labFishFam[["samples"]]$RowNames)]

# Convert to phyloseq object 
## makes otu table for phyloseq
read.phylo <-
  phyloseq::otu_table(t(labFishFam[["reads"]]), taxa_are_rows = T)

## make sample table for phyloseq
sample.phylo <- phyloseq::sample_data(labFishFam[["samples"]])

### merge phyloseq objects
PhyloDat  <- 
  phyloseq::phyloseq(read.phylo, sample.phylo)

##################################################################################################
# Test treatment effect on taxa (family level) relative abundance
# using LINDA - https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02655-5 
#################################################################################################

# test differential abundance of taxa levels among treatment using LMM approach
meta <- data.frame(Contam = factor(PhyloDat@sam_data$Contam),
                   Inj = factor(PhyloDat@sam_data$Inj),
                   Or = factor(PhyloDat@sam_data$ContPop),
                   Pop = factor(PhyloDat@sam_data$Pop),
                   Bac = factor(PhyloDat@sam_data$Bac),
                   Size = PhyloDat@sam_data$Taille_mm_M,
                   Sex = factor(PhyloDat@sam_data$Sexe),
                   Session = factor(PhyloDat@sam_data$Session))
rownames(meta) <- rownames(PhyloDat@sam_data)

# relevel the treatment
meta$Contam <- factor(meta$Contam, levels = c("NC", "C"))
meta$Inj <- factor(meta$Inj, levels = c("PBS", "AMIX"))
meta$Or <- factor(meta$Or, levels = c("NC", "C"))

# remove the motus detected in less than 5 samples
## it remove 17 families: 
# "Acidothermaceae"     "Armatimonadaceae"    "Caedibacteraceae"    "Chthoniobacteraceae"
# "Fimbriimonadaceae"   "Gallionellaceae"     "Gemmatimonadaceae"   "Hydrogenophilaceae" 
# "Hymenobacteraceae"   "Lachnospiraceae"     "Methyloligellaceae"  "Micavibrionaceae"   
# "Nitrosomonadaceae"   "Pedosphaeraceae"     "Sandaracinaceae"     "Solimonadaceae"     
# "Sulfuricellaceae"
toRemove <- which(apply(PhyloDat@otu_table, 1, function(x) length(which(x != 0 ))) <= 5)
if(length(toRemove) != 0){
  PhyloDat@otu_table <- PhyloDat@otu_table[-toRemove,]
}

linda.obj <- MicrobiomeStat::linda(feature.dat = as.data.frame(PhyloDat@otu_table), 
                                   meta.dat = meta,
                                   feature.dat.type = "count",
                                   formula = "~ Contam + Inj + Or + Sex + Size + (1|Session/Bac) + (1|Pop)", 
                                   alpha = 0.05, p.adj.method = "fdr",
                                   is.winsor = F,
                                   adaptive = T, 
                                   prev.filter = 0)

# check the results 
MicrobiomeStat::linda.plot(
  linda.obj,
  c('ContamC', "InjLPS", "OrC", "SexI", "SexM"),
  titles = c('Contam: NC vs C', 'Inj: PBS vs AMIX', 'Origin: NC vs C', 'Sex: F vs. I', 'Sex: F vs. M'),
  alpha = 0.05,
  lfc.cut = 0.5,
  legend = TRUE,
  directory = NULL,
  width = 11,
  height = 8
)

# specify p value and log2 fold change threshold
FcTresh <- c(-0.5,0.5)
pTresh <- -log10(0.05)

# specify the colors to represent significance thresholds
colors <- c(adjustcolor("#808080", alpha.f = 0.5),
            adjustcolor("#420693", alpha.f = 0.5))
names(colors) <- c("NS", "p-adj. & Log2FC")

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
if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "TaxoLinda"))) == 0) {
  dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "TaxoLinda"))
} 
toSave <- file.path(savingDir, "Normalized_Data", "Signif_Effects", "Signif_Plots", "TaxoLinda", "Family_Contam.tif")
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
  ylim = c(-5, 20),
  xlim = c(-5, 5),
  ylab = bquote( ~ -Log[10] ~ 'p-adj.'),
  xlab = bquote( ~ Log[2] ~ 'fold change'),
  vlines = FcTresh,
  hlines = pTresh,
  cex.lab = 1.2,
  leg.order = c("NS", "p-adj.","Log2FC", "p-adj. & Log2FC"),
  cex.pts = cexPTS
)
# add legend for dot size (number of reads)
leg <- c(min(cexPTS), max(cexPTS)/2, max(cexPTS))
legValues <- ((leg - DotSize[1]) / (DotSize[2] - DotSize[1])) * 
  (maxVal - minVal) + minVal
text(0.45, 19.65,"# of reads:", cex = 1.0)
legend(
  x = 1.15,
  y = 20.85,
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
## for contamination and for the 15 top phylum (there is only 13 here)
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
mycolors <-
  colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(nb.cols)
RelaAB_Plot_Contam_Top <-
  ggplot2::ggplot(tmp_merged,
                  ggplot2::aes(x = ContamOrd, y = Relvalue, fill = variable)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::labs(x = NULL, y = "Relative abundance", fill = "Taxa")  + ggplot2::theme_bw() +
  ggplot2::scale_fill_manual(values = mycolors) +
  ggplot2::theme_classic() +
  ggplot2::theme(text = ggplot2::element_text(size=20)) +
  ggplot2::ggtitle("")

# save the plot
toSave <- file.path(savingDir, "Normalized_Data", "Signif_Effects", "Signif_Plots", 
                    "TaxoLinda", "EC_Contam_RelAb_Family_TOP.svg")
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

## for Sex (only between F and M since there is no difference with I)
i = "SexM"
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
if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "TaxoLinda"))) == 0) {
  dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "TaxoLinda"))
} 
toSave <- file.path(savingDir, "Normalized_Data", "Signif_Effects", "Signif_Plots", "TaxoLinda", "Family_SexFvsM.tif")
if(file.exists(toSave)){
  unlink(toSave)
}
grDevices::tiff(filename = toSave,
                width = 1200,
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
  ylim = c(-1, 5),
  xlim = c(-1.5, 1.5),
  ylab = bquote( ~ -Log[10] ~ 'p-adj.'),
  xlab = bquote( ~ Log[2] ~ 'fold change'),
  vlines = FcTresh,
  hlines = pTresh,
  cex.lab = 1.2,
  leg.order = c("NS", "p-adj.","Log2FC", "p-adj. & Log2FC"),
  cex.pts = cexPTS
)
# add legend for dot size (number of reads)
leg <- c(min(cexPTS), max(cexPTS)/2, max(cexPTS))
legValues <- ((leg - DotSize[1]) / (DotSize[2] - DotSize[1])) * 
  (maxVal - minVal) + minVal
text(-1.2, 4.35,"# of reads:", cex = 1.0)
legend(
  x = -0.9,
  y = 4.65,
  legend = round(legValues),
  horiz = T,
  pch = 16,
  col = adjustcolor("#000000", alpha.f = 0.8),
  pt.bg =  adjustcolor("#000000", alpha.f = 0.8),
  pt.cex = leg,
  #y.intersp = c(1, 2, 4),
  x.intersp = c(1, 1, 2),
  text.width = c(0.2, 0.6, 0.6),
  bty = "n"
)
dev.off()

## for Sex (Male vs female) and for the 15 top phylum (there is only 13 here)
# Aggregate per treatment
NFunc <- 15 # number of level to display

tmp <-
  reshape2::melt(stats::aggregate(t(PhyloDat@otu_table), by = list(
    meta$Sex), sum))

# Sorting and selecting top Ntax abundant MOTUS
summedDat <- stats::aggregate(value ~ variable, data = tmp, FUN = sum)
Top <- summedDat[order(summedDat$value, decreasing = T), ][1:NFunc, ]

# keep only the top Ntax abundant MOTUS
tmpTop <- tmp[which(tmp[["variable"]] %in% Top[["variable"]]),]

# Splitting the dataframe
tmp_F <- tmpTop[tmpTop$Group.1 == "F", ]
tmp_M <- tmpTop[tmpTop$Group.1 == "M", ]

# compute relative abundances
tmp_F$Relvalue <- tmp_F$value/sum(tmp_F$value)
tmp_M$Relvalue <- tmp_M$value/sum(tmp_M$value)

# Merging the dataframes
tmp_merged <- rbind(tmp_F, tmp_M)

# Plot 
newOrder <- c("F", "M")
tmp_merged[["SexOrd"]] <- factor(tmp_merged[["Group.1"]], levels = newOrder)

nb.cols <- NFunc
mycolors <-
  colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(nb.cols)
RelaAB_Plot_Sex_Top <-
  ggplot2::ggplot(tmp_merged,
                  ggplot2::aes(x = SexOrd, y = Relvalue, fill = variable)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::labs(x = NULL, y = "Relative abundance", fill = "Taxa")  + ggplot2::theme_bw() +
  ggplot2::scale_fill_manual(values = mycolors) +
  ggplot2::theme_classic() +
  ggplot2::theme(text = ggplot2::element_text(size=20)) +
  ggplot2::ggtitle("")

# save the plot
toSave <- file.path(savingDir, "Normalized_Data", "Signif_Effects", 
                    "Signif_Plots", "TaxoLinda", "EC_Sex_RelAb_Family_TOP.svg")

if(file.exists(toSave)){
  unlink(toSave)
} 
ggplot2::ggsave(
  file = toSave,
  plot = RelaAB_Plot_Sex_Top,
  width = 9,
  height = 6,
  dpi = 300
)
