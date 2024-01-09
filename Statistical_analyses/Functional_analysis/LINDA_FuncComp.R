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

# merge FullDat to the metabar sample object
labFish[["samples"]]["RowNames"] <- rownames(labFish[["samples"]])
labFish[["samples"]] <- merge(
  x = FullDat,
  y = labFish[["samples"]],
  by.x = "Ind",
  by.y = "Num_prlvt_Euth",
  all.y = TRUE
)
rownames(labFish[["samples"]]) <- labFish[["samples"]]$Ind
rownames(labFish[["reads"]]) <- labFish[["samples"]][["Ind"]][match(rownames(labFish[["reads"]]), labFish[["samples"]]$RowNames)]

# Convert to phyloseq object 
## makes otu table for phyloseq
read.phylo <-
  phyloseq::otu_table(t(labFish[["reads"]]), taxa_are_rows = T)

## make sample table for phyloseq
sample.phylo <- phyloseq::sample_data(labFish[["samples"]])

### merge phyloseq objects
PhyloDat  <- 
  phyloseq::phyloseq(read.phylo, sample.phylo)


###############################################################
# import Picrust2 results
###############################################################
## import EC data including description column (tibble)
ECPath <-  as.data.frame(readr::read_tsv("D:/POSTDOC_INP_GOLFECH_2023/Outputs/Picrust2Output/EC_pathways_out/pred_metagenome_unstrat_descrip.tsv.gz"))

# manual change of the PWY-5182 and ARGORNPROST-PWY pathways which are not found in the Metacyc db 
# but found trough manual search using description given by picrust2 
ECPath[grep("PWY-5182", ECPath[,1]), 1] <- "TOLUENE-DEG-3-OH-PWY"
ECPath[grep("ARGORNPROST-PWY", ECPath[,1]), 1] <- "ARGDEGRAD-PWY"

######################################################################
# prepare the database containing different levels of class/ pathways
######################################################################

# Save the pathway in a csv file to import it into the metacyc smart table (https://metacyc.org/smarttables)
# briefly, in the "operation tab", select "SmartTable of Object from Text Entry" paste the pathways from the file created below
# and in the "Add transform column" tab select "ontology - direct parent of entity", do it again by selecting the created column 
# to retrieve the upper level of the pathway until having the full pathway for each entry.
# write.table(ECPath[["pathway"]], file = file.path(savingDir, "MetacycPathways.csv"), sep = ";") 
# then export the smart table to spreadsheet file using frame Id 
# no need to export ECPath[["pathway"]] and make smart table on metacyc website since we provide it in the repository (see "SelectedPathway-of-MetaCyc_FrameId.txt")
# So then import reference db from metacyc (smart table saved on 13-Dec-2023) - for the sack of simplicity and reproducibility, 
# the file is available in the git hub repository within the "functional analysis" directory
MetaCycDb <- read.delim(file.path(savingDir, "SelectedPathway-of-MetaCyc_FrameId.txt"),
                        quote = '', sep = "\t", na.strings = c("", "NA"))

# split pathways with several entities (separated by //)
# in case a sub pathway is not found, retrieve the upper one
MetacycDf <- data.frame(MetaCycDb[c(1,2)])
for (i in 3:ncol(MetaCycDb)) {
  Tmp <- lapply(MetaCycDb[, i], function(x)
    strsplit(x, "//"))
  DfnCol <- max(unlist(lapply(Tmp, function(x)
    length(unlist(
      x
    )))))
  DfTemp <- data.frame(matrix(ncol = DfnCol, nrow = nrow(MetaCycDb)))
  for (r in seq_along(Tmp)) {
    input <- unlist(Tmp[[r]])
    if(length(input) == 0 | length(input) == 1 && is.na(input)){
      input <- "Unknown"
    }
    if(i == 3 && length(input) == 1 && input == "Super-Pathways"){
      input <- MetaCycDb[r, "Common.Name"]
    }
    for (c in seq_along(input)) {
      DfTemp[r, c] <- input[c]
    }
    NoNa <- which(!is.na(DfTemp[r,]))
    if (length(NoNa) > 1) {
      NoNa <- max(NoNa)
    }
    DfTemp[r, which(is.na(DfTemp[r,]))] <- DfTemp[r, NoNa]
  }
  MetacycDf <- cbind(MetacycDf, DfTemp)
}

# remove potential artifactual space in pathways (trim leading and trailing white space)
MetacycDf <- apply(MetacycDf, 2, function(x) trimws(x, which = "both"))

# also remove uninformative strings such as "Pathways", "FRAMES"
# and replace it with the first informative level of pathway

toSub <-
  c("",
    "NA",
    "Unknown",
    "Pathways",
    "FRAMES",
    "THINGS",
    "GeneralizedReactions",
    "Generalized-Reactions",
    "Superpathways",
    "Super-Pathways",
    "All Pathways and Reactions")
for (i in 2:ncol(MetacycDf)) {
  for(j in toSub){
    MetacycDf[which(MetacycDf[, i] == j), i] <-
      MetacycDf[which(MetacycDf[, i] == j), i - 1]
  }
}

MetacycDf <- stats::setNames(as.data.frame(MetacycDf),
                             c(tolower(colnames(MetacycDf)[c(1,2)]), paste("Path_level", seq(1:(ncol(MetacycDf)-2)), sep = "_")))

# manually rename the pathway (GLYCOLYSIS // GLYCOLYSIS-VARIANTS) returned by the metacyc df to match EC data
MetacycDf[which(MetacycDf[["pathway"]] == "GLYCOLYSIS // GLYCOLYSIS-VARIANTS"), "pathway"] <- "GLYCOLYSIS"

# check the number of pathway levels created and select the desired one
pathLev <- grep("Path_level", names(MetacycDf), value = T)
pathLev
Tofind <- "Path_level_5"

# make a quick function to homogeneize path level
homogenPathLev <- function(df) {
  # Number of path level columns
  pathLev <- ncol(df) - 2  # Assuming the first two columns are 'common.name' and 'pathway'
  
  for (i in 3:(pathLev+1)) {
    # Iterate through each element in the current column
    for (j in 1:nrow(df)) {
      # Get the corresponding value in the next column
      nextCol <- df[j, i+1]
      
      # Check if this value is present anywhere in the current column
      if (nextCol %in% df[, i]) {
        # Replace the current value with the value from the next column
        df[j, i] <- nextCol
      }
      # If not present, do not replace
    }
  }
  # Remove duplicated columns
  df <- df[!duplicated(as.list(df))]
  pathLevFinal <- setdiff(names(df), c("common.name", "pathway"))
  newName <- paste0("Path_level_", seq_along(pathLevFinal))
  names(df)[names(df) %in% pathLevFinal] <- newName
  return(df)
}

# Apply the function to the data frame
MetacycDf <- homogenPathLev(MetacycDf)

# merge the pathways with the current dataset to retrieve only useful pathways
EC <- merge(ECPath, MetacycDf, all.x = T, by = "pathway")

# Remove pathway to transform the dataset in matrix
PathCol <- grep("Path_level", names(EC))
ECMat <- EC[,-c(1, 2, PathCol, min(PathCol)-1)]

## remove character columns (class) and use it as rownames and keep only the selected pathway level
ECMat <- as.matrix(ECMat)
rownames(ECMat) <- EC[, which(names(EC) == Tofind)]
class(ECMat) <- "numeric"

## order the function (pathways class) alphabetically
ECMat <- ECMat[order(rownames(ECMat), decreasing = TRUE), ]

## sum duplicates pathways' reads
ECMat <- stats::aggregate(ECMat, list(row.names(ECMat)), sum)
rownames(ECMat) <- ECMat[[1]]
ECMat <- ECMat[-1]

# homogeneize pathways appearence (only first letter capitalized, removing hyphen as word separation)
## make a quick function to do it
HomogenStrings <- function(x) {
  # lower lettering
  lower <- tolower(x)
  # Replace hyphens that are not after a single letter word
  noHyph <- gsub("(?<![[:space:]]\\b\\w)\\-", " ", lower, perl = TRUE)
  # Replace abbreviations
  NoAbb <- gsub("\\bdeg\\b", "degradation", noHyph)
  NoAbb <- gsub("\\bsyn\\b", "synthesis", NoAbb)
  NoAbb <- gsub("\\bbiosyn\\b", "biosynthesis", NoAbb)
  # capitalize the first letter of each word
  words <- strsplit(NoAbb, " ")[[1]]
  capWords <- sapply(words, function(word) {
    if (nchar(word) > 0) {
      return(paste0(toupper(substring(word, 1, 1)), substring(word, 2)))
    } else {
      return("")
    }
  })
  Res <- paste(capWords, collapse = " ")
  return(Res)
}
rownames(ECMat) <- unlist(lapply(rownames(ECMat), HomogenStrings))

##################################################################################################
# Test treatment effect on functions relative abundance
# using LINDA - https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02655-5 
#################################################################################################

# prepare metadata for lmm testing
meta <- data.frame(Contam = factor(PhyloDat@sam_data$Contam),
                   Inj = factor(PhyloDat@sam_data$Inj),
                   Or = factor(PhyloDat@sam_data$ContPop),
                   Pop = factor(PhyloDat@sam_data$Pop),
                   Bac = factor(PhyloDat@sam_data$Bac),
                   Size = PhyloDat@sam_data$Taille_mm_M,
                   Sex = factor(PhyloDat@sam_data$Sexe),
                   Session = factor(PhyloDat@sam_data$Session))
rownames(meta) <- PhyloDat@sam_data$Ind

# relevel the treatment
meta$Contam <- factor(meta$Contam, levels = c("NC", "C"))
meta$Inj <- factor(meta$Inj, levels = c("PBS", "AMIX"))
meta$Or <- factor(meta$Or, levels = c("NC", "C"))

# remove the pathways detected in less than 5 samples (there is no pathway to remove in this case)
toRemove <- which(apply(ECMat, 1, function(x) length(which(x != 0 ))) <= 5)
if(length(toRemove) != 0){
  ECMat <- ECMat[-toRemove,]
}

# run the model
linda.obj <- MicrobiomeStat::linda(feature.dat = as.data.frame(ECMat), 
                                   meta.dat = meta,
                                   feature.dat.type = "count",
                                   formula = "~ Contam + Inj + Or + Sex + Size + (1|Session/Bac) + (1|Pop)", 
                                   alpha = 0.05, p.adj.method = "fdr",
                                   is.winsor = F,
                                   adaptive = T)

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
linda.obj$output[[i]][["readsCount"]] <- apply(ECMat, 1, function(x) sum(x, na.rm = T))
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
if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Functions"))) == 0) {
  dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Functions"))
} 
if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Functions", "FuncLinda"))) == 0) {
  dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Functions", "FuncLinda"))
} 
toSave <- file.path(savingDir, "Normalized_Data", "Signif_Effects", "Signif_Plots", "Functions", "FuncLinda", "EC_Contam.tif")
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
text(-0.05, 4.95,"# of reads:", cex = 1.0)
legend(
  x = 0.1,
  y = 5.2,
  legend = round(legValues),
  horiz = T,
  pch = 16,
  col = adjustcolor("#000000", alpha.f = 0.8),
  pt.bg =  adjustcolor("#000000", alpha.f = 0.8),
  pt.cex = leg,
  #y.intersp = c(1, 2, 4),
  x.intersp = c(1, 1, 2),
  text.width = c(0.2, 0.4, 0.6),
  bty = "n"
)
dev.off()

# draw relative abundance plot 
## for contamination and for the 15 top functions
# Aggregate per treatment
NFunc <- 15 # number of level to display

tmp <-
  reshape2::melt(stats::aggregate(t(ECMat), by = list(
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
RelaAB_Plot_Top <-
  ggplot2::ggplot(tmp_merged,
                  ggplot2::aes(x = ContamOrd, y = Relvalue, fill = variable)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::labs(x = NULL, y = "Relative abundance", fill = "Function")  + ggplot2::theme_bw() +
  ggplot2::scale_fill_manual(values = mycolors) +
  ggplot2::theme_classic() +
  ggplot2::theme(text = ggplot2::element_text(size=20)) +
  ggplot2::ggtitle("")

# save the plot
toSave <- file.path(savingDir, "Normalized_Data", "Signif_Effects", "Signif_Plots", "Functions", "FuncLinda", "EC_Contam_RelAb_TOP.svg")

if(file.exists(toSave)){
  unlink(toSave)
} 
ggplot2::ggsave(
  file = toSave,
  plot = RelaAB_Plot_Top,
  width = 9,
  height = 6,
  dpi = 300
)

## for contamination and for the 15 bottom functions
# Aggregate per treatment
NFunc <- 15 # number of level to display

tmp <-
  reshape2::melt(stats::aggregate(t(ECMat), by = list(
    meta$Contam), sum))

# Sorting and selecting Bottom Ntax abundant MOTUS
summedDat <- stats::aggregate(value ~ variable, data = tmp, FUN = sum)
Bottom <- summedDat[order(summedDat$value, decreasing = F), ][1:NFunc, ]

# keep only the Bottom Ntax abundant MOTUS
tmpBottom <- tmp[which(tmp[["variable"]] %in% Bottom[["variable"]]),]

# Splitting the dataframe
tmp_C <- tmpBottom[tmpBottom$Group.1 == "C", ]
tmp_NC <- tmpBottom[tmpBottom$Group.1 == "NC", ]

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
RelaAB_Plot_Bottom <-
  ggplot2::ggplot(tmp_merged,
                  ggplot2::aes(x = ContamOrd, y = Relvalue, fill = variable)) +
  ggplot2::geom_bar(stat = "identity") +
  ggplot2::labs(x = NULL, y = "Relative abundance", fill = "Function")  + ggplot2::theme_bw() +
  ggplot2::scale_fill_manual(values = mycolors) +
  ggplot2::theme_classic() +
  ggplot2::theme(text = ggplot2::element_text(size=20)) +
  ggplot2::ggtitle("")

# save the plot
toSave <- file.path(savingDir, "Normalized_Data", "Signif_Effects", "Signif_Plots", "Functions", "FuncLinda", "EC_Contam_RelAb_Bottom.svg")

if(file.exists(toSave)){
  unlink(toSave)
} 
ggplot2::ggsave(
  file = toSave,
  plot = RelaAB_Plot_Bottom,
  width = 9,
  height = 6,
  dpi = 300
)

# list the over or lower abundant functions
rownames(linda.obj$output[[i]][which(linda.obj$output[[i]]$VolcanoGroups == "p-adj. & Log2FC"), ])
