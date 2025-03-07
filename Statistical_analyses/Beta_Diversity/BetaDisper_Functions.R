# Script Title: Compute Beta Diversity Indices and Test Treatment Effects using LMM and Dredge in fish gut microbiota
#
# Author: Quentin PETITJEAN
# Date Created: 02/02/2024
# Last Modified: 02/02/2024
# ==============================================================================
# Requirements:
# - R version 4.2.3
# - Packages:
#   - vegan v2.6-4: For computing betadisper, running permutation tests, and environmental fitting.
#   - kableExtra v1.3.4: For creating enhanced HTML summary tables of statistical results.
# - Source Files:
#   - multiNorm.R: Custom function to compute multiple normalization methods.
#   - DredgedLMM.R: Custom function to perform model dredging (summary, ANOVA table, R² computation).
#   - multiDist.R: Custom function to compute multiple beta diversity distance metrics.
#   - BxpltFunc.R: Custom function to generate customized boxplots.
# ==============================================================================
# Script Overview:
# This script computes beta diversity indices from normalized functional inference data 
# (derived from PICRUSt2 enzyme classification pathway predictions) and tests the effects 
# of experimental treatments (e.g., Contamination, Injection, Population) on beta diversity 
# using linear mixed models (LMM) with a dredging approach (AICc model selection). The workflow 
# includes:
#
# 1. Specifying file paths for input data and output directories.
# 2. Importing the experimental design data (DesignData.csv), 
#    and the cleaned fish gut microbiota dataset (as a metabaR object).
# 3. Importing Picrust2 enzyme classification pathway predictions (with manual adjustment of specific
#    pathway IDs) and merging these with the sample metadata.
# 4. Normalizing the functional inference data using multiple normalization methods.
# 5. Computing beta diversity distance matrices using multiple metrics (e.g., Bray, Jaccard, 
#    Robust.Aitchison, Hellinger).
# 6. Deriving ordination coordinates via PCOA and NMDS from the distance matrices.
# 7. Testing treatment effects on beta diversity through LMM and model dredging, and evaluating 
#    the dispersion of sample groups using betadisper and permutation tests.
# 8. Generating and saving diagnostic plots (TIFF format) and summarizing the statistical results 
#    in HTML tables.
#
# ==============================================================================
# Usage:
# 1. Ensure that the required input files (DesignData.csv, and the metabaR object) 
#    are available in the specified directory (savingDir).
# 2. Update the 'savingDir' variable to reflect your local directory structure.
# 3. Install and load all required R packages and source files before running the script.
# 4. Run the script in an R environment; outputs (plots, HTML tables, and RDS files) will be saved 
#    in the designated output directories.
# ==============================================================================


##############################################
#       	Install needed packages            #
##############################################

if(!require(kableExtra)){
  install.packages("kableExtra")
}
if(!require(vegan)){
  install.packages("vegan")
}

#######################
# specify some path   #
#######################
# the path to the directory from where file should be both imported and saved
savingDir <- "w:/POSTDOC_INP_GOLFECH_2023/Outputs"

# directory where the output files should be saved
DirSave <- file.path(savingDir, "Normalized_Data")

##############################################
#   Import some custom functions             #
##############################################
source("https://raw.githubusercontent.com/qpetitjean/Multiple-stressors-effects-wild-fish-gut-microbiota/main/R_Func/BxpltFunc.R") # import a function compute draw customized boxplots

#######################
# Import the data     #
#######################

# import the full dataset from the experiment (treatment & behavioral & physiological measures)
FullDat <-
  read.csv2(file.path(savingDir, "Data/DesignData", "DesignData.csv"), dec = ".", sep = ";")

# import the distance matrix according to the beta diversity indices selected 
DistList <- readRDS(file.path(DirSave, "DistMatrices_Functions.rds"))

# import the scaled coordinates of the samples in a 2d space (PCOA or NMDS), derived from the distance matrix
BetaDivCoords <- readRDS(file.path(DirSave, "BetaDivCoords_Functions.rds"))

# import the results of the LMM approach on beta diversity indices
LMMRes <- readRDS(file.path(DirSave, "Signif_Effects", "Functions", "LMMRes_Functions.rds"))

#########################################################################################
# compute betadisper and run permutation test over the normalized data                  #
#########################################################################################

# initialize an empty list to store the results of betadisper and permutation test
betadispRes <- list()

for(i in names(LMMRes)){
  for(k in names(LMMRes[[i]])){
    for(l in names(LMMRes[[i]][[k]])){
      # merge the BetaDiv coords to the full dataset of the experiment
      TempDist <- as.data.frame(BetaDivCoords[[l]][[k]][[i]][[ifelse(l == "PCOA", "vectors", "points")]])
      TempDist[["Samples"]] <- rownames(TempDist)
      TempDist <- merge(
        x = FullDat,
        y = TempDist,
        by.x = "Ind",
        by.y = "Samples",
        all = TRUE
      )
      
      # retrieve the significant effect from LMM results
      toKeep <-
        unique(c(
          rownames(LMMRes[[i]][[k]][[l]][[ifelse(l == "PCOA", "Axis.1", "MDS1")]]$ModAnov)[-1],
          rownames(LMMRes[[i]][[k]][[l]][[ifelse(l == "PCOA", "Axis.2", "MDS2")]]$ModAnov)[-1]
        ))
      
      # iterate trough significant effect highlighted in the LMM approach
      for(j in toKeep) {
        if (is.character(TempDist[[j]]) | is.factor(TempDist[[j]])) {
          groups <-
            TempDist[match(labels(DistList[[i]][[k]]), TempDist[["Ind"]]), ]
          rownames(groups) <- groups[["Ind"]]

          betadispRes[[i]][[k]][[j]][["betaDisper"]] <-
            vegan::betadisper(
              DistList[[i]][[k]],
              groups[[j]],
              type = "median",
              add = "lingoes",
              bias.adjust = TRUE
            )
          
          betadispRes[[i]][[k]][[j]][["PermTest"]] <-
            vegan::permutest(
              betadispRes[[i]][[k]][[j]][["betaDisper"]],
              pairwise = TRUE,
              permutations = 999,
              by = "margin"
            )
          
          betadispRes[[i]][[k]][[j]][["FDRadjustedP"]] <-
            p.adjust(betadispRes[[i]][[k]][[j]][["PermTest"]]$pairwise$observed, method = "fdr")
          
          TempDistOrd <- TempDist[order(TempDist[[j]]),]
          
          tempDf <-
            data.frame(
              gr = betadispRes[[i]][[k]][[j]][["betaDisper"]][["group"]],
              dist = betadispRes[[i]][[k]][[j]][["betaDisper"]][["distances"]],
              Ind = rownames(betadispRes[[i]][[k]][[j]][["betaDisper"]][["vectors"]])
            )

          tempDfPlot <- merge(
            x = tempDf,
            y = TempDistOrd,
            by.x = "Ind",
            by.y = "Ind",
            all = TRUE
          )
          
          if(l == "PCOA") {
            VarAxis1 <- (BetaDivCoords[[l]][[k]][[i]][["values"]][["Eigenvalues"]][[1]]/
                           sum(BetaDivCoords[[l]][[k]][[i]][["values"]][["Eigenvalues"]])) * 100
            VarAxis2 <- (BetaDivCoords[[l]][[k]][[i]][["values"]][["Eigenvalues"]][[2]]/
                           sum(BetaDivCoords[[l]][[k]][[i]][["values"]][["Eigenvalues"]])) * 100
          }
          ordLab <- c("NC", "C")
          ordCol <- as.character(TempDist$Contam_col[match(ordLab, TempDist$Contam)])
          ggplot2::ggplot(tempDfPlot,
                          ggplot2::aes(x = .data[[ifelse(l == "PCOA", "Axis.1", "MDS1")]], 
                                       y = .data[[ifelse(l == "PCOA", "Axis.2", "MDS2")]], 
                                       colour = .data[["Contam_col"]], 
                                       fill = .data[["Contam_col"]])
          ) +
            ggplot2::geom_point(size = 4) +
            #ggplot2::annotate("text", x = tempDfPlot[[ifelse(l == "PCOA", "Axis.1", "MDS1")]], 
            #         y = tempDfPlot[[ifelse(l == "PCOA", "Axis.2", "MDS2")]], 
            #         label = round(tempDfPlot[["dist"]]), 2) +
            ggplot2::ggtitle(ifelse(l == "PCOA", i, paste0(
              i,
              " (stress = ",
              signif(BetaDivCoords[[l]][[k]][[i]][["stress"]], digits = 2),
              ")"
            ))) +
            ggplot2::xlab(ifelse(l == "PCOA", paste0("Axis.1 (", round(VarAxis1, digits = 1), "%", ")"), "MDS1")) +
            ggplot2::ylab(ifelse(l == "PCOA", paste0("Axis.2 (", round(VarAxis2, digits = 1), "%", ")"), "MDS2")) +
            ggplot2::stat_ellipse(geom = "polygon",
                                  ggplot2::aes(fill = .data[["Contam_col"]], colour = .data[["Contam_col"]]), 
                                  alpha = 0.25,
                                  type = "t") + 
            ggplot2::scale_colour_identity(guide = "legend", 
                                           breaks = ordCol, 
                                           labels = ordLab, 
                                           name = "Contamination") + 
            ggplot2::scale_fill_identity(guide = "legend", 
                                         breaks = ordCol, 
                                         labels = ordLab, 
                                         name = "Contamination")+
            ggplot2::theme_classic()
          

          if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Functions", "BetadisperPlots"))) == 0) {
            dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Functions", "BetadisperPlots"))
          }
          if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Functions", "BetadisperPlots", i))) == 0) {
            dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Functions", "BetadisperPlots", i))
          }
          if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Functions", "BetadisperPlots", i, k))) == 0) {
            dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Functions", "BetadisperPlots", i, k))
          }
          if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Functions", "BetadisperPlots", i, k, l))) == 0) {
            dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Functions", "BetadisperPlots", i, k, l))
          }
          ## save the figs representing distance to spatial median accroding to the treatments
          toSave <- file.path(DirSave, "Signif_Effects", "Signif_Plots", "Functions", "BetadisperPlots", i, k, l, paste0(paste("Betadisper-Plots", i, k, l, j, sep = "_"), ".tif" ))
          if(file.exists(toSave)){
            unlink(toSave)
          } 
          if(j == "Contam" | j == "ContPop"){
            newOrder <- c("NC", "C")
          } else if(j == "Inj") {
            newOrder <- c("PBS", "AMIX")
          } else {
            newOrder <- levels(factor(tempDfPlot[[j]]))
          } 
          
          tempDfPlot[["Ordered"]] <- factor(tempDfPlot[[j]], levels = newOrder)
          grDevices::tiff(filename = toSave,
                          width = 1200,
                          height = 1200,
                          res = 300,
                          compression = "lzw",
                          pointsize = 12)
          
          Bxplt(
            x = tempDfPlot[["Ordered"]],
            y =  tempDfPlot[["dist"]],
            fill = tempDfPlot[match(levels(tempDfPlot[["Ordered"]]), tempDfPlot[[j]]), paste(j, "col", sep = "_")],
            colpts = tempDfPlot[[paste(j, "col", sep = "_")]],
            las = ifelse(length(unique(tempDfPlot[[j]])) > 4, 2, 1),
            boxwex = 0.5,
            ylab = "Distance to spatial median",
            cex.axis = ifelse(length(unique(tempDfPlot[[j]])) > 4, 0.7, 1.1),
            main = paste("BetaDisper", paste(i, l, sep = " "), sep = " : "),
            yscale = 6)
          
          dev.off()
          
        } else if (is.numeric(TempDistOrd[[j]])) {
          next()
        }
      }
    }
  }
}

saveRDS(betadispRes, file = file.path(DirSave, "Signif_Effects", "Functions", "BetaDispRes_Functions.rds"), compress = TRUE)
#betadispRes <- readRDS(file.path(DirSave, "Signif_Effects", "Functions", "BetaDispRes_Functions.rds"))

################################################################
# summarize the results in html table for each normalization   #
################################################################

BestMods <- data.frame(
  mods = names(betadispRes),
  modsId = c(
    "Original",
    "Original (log-transformed)",
    "Proportion",
    "Proportion (log-transformed)",
    "CSS",
    "CSS (log-transformed)",
    "Rarefied",
    "Rarefied (log-transformed)",
    "DESEQ",
    "DESEQ (log-transformed)",
    "TMM",
    "TMM (log-transformed)"
  )
)

# retrieve information from model summary and anova table
options(scipen = 999)
VRList <- list()
distances <- names(betadispRes[[1]])

for(Ord in c("PCOA", "NMDS")){
  Val <- setNames(lapply(distances, function(i) {
    ValuesRes <- data.frame()
    for (r in seq(nrow(BestMods))) {
      VarTest <- betadispRes[[BestMods[r, "mods"]]][[i]]
      effects <- list()
      if(is.null(names(VarTest))){
        next
      }
      for(v in names(VarTest)){
        effects[[v]] <- betadispRes[[BestMods[r, "mods"]]][[i]][[v]]$PermTest$tab
        rownames(effects[[v]])[1] <- v
      }
      Values <- do.call(rbind, effects)
      Values <- signif(Values, digits = 3)
      Values <- cbind(var = rownames(Values), Values)
      rownames(Values) <- NULL
      splitVar <- strsplit(Values$var, "\\.")
      Values["var"] <- sapply(splitVar, "[", 2)
      Values["mod"] <- sapply(splitVar, "[", 1)
      Values <- cbind(mods = rep(BestMods[r, "mods"], nrow(Values)), Values)
      ValuesRes <- rbind(ValuesRes, Values, stringsAsFactors = FALSE)
      rownames(ValuesRes) <- NULL
    }
    # leave the two first column without names
    names(ValuesRes)[c(1, 2)] <- c("", "")
    
    # specify the new predictor names (to display in the table)
    PredictorsNames <- data.frame(toreplace = unique(ValuesRes[[2]]), replacement = NA)
    PredictorsNames$replacement[PredictorsNames$toreplace == "Inj"] <- "Imm. Chall."
    PredictorsNames$replacement[PredictorsNames$toreplace == "Sexe"] <- "Sex"
    PredictorsNames$replacement[PredictorsNames$toreplace == "Taille_mm_M"] <- "Size (cm)"
    PredictorsNames$replacement[PredictorsNames$toreplace == "Contam"] <- "Contamination"
    PredictorsNames$replacement[PredictorsNames$toreplace == "ContPop"] <- "Origin"
    PredictorsNames$replacement[PredictorsNames$toreplace == "Contam:ContPop"] <- "Contam.:Origin"
    PredictorsNames$replacement[PredictorsNames$toreplace == "Contam:Inj"] <- "Contam.:Imm. Chall."
    PredictorsNames <- PredictorsNames[!is.na(PredictorsNames$replacement),]
    
    # rename the predictors
    for (pn in PredictorsNames[["toreplace"]]) {
      ValuesRes[[2]][which(ValuesRes[[2]] == pn)] <-
        PredictorsNames[["replacement"]][which(PredictorsNames[["toreplace"]] == pn)]
    }
    
    # replace p.value with large decimal by approximations
    ValuesRes[["Pr(>F)"]][which(ValuesRes[["Pr(>F)"]] < 0.0001)] <-
      "<0.0001"
    ValuesRes[["Pr(>F)"]][which(ValuesRes[["Pr(>F)"]] > 0.0001 &
                                  ValuesRes[["Pr(>F)"]] < 0.01)] <- "<0.01"
    ValuesRes[["Pr(>F)"]][which(ValuesRes[["Pr(>F)"]] > 0.01 &
                                  ValuesRes[["Pr(>F)"]] < 0.05)] <- "<0.05"
    
    return(ValuesRes)
  }), distances)
  VRList[[Ord]] <- Val
}
options(scipen = 0)

# generate the table (here one table per metric - e.g., richness, shannon, ... - including all normalisation methods)
if (length(list.dirs(file.path(DirSave, "Signif_Effects", "BetaDispFunctionsBestModTab"))) == 0) {
  dir.create(file.path(DirSave, "Signif_Effects", "BetaDispFunctionsBestModTab"))
} 


for(Ord in c("PCOA", "NMDS")){
  for(m in distances){
    Values <- VRList[[Ord]]
    ValIndex <- which(names(VRList[[Ord]]) == m)
    toTab <- Values[[ValIndex]]
    toTab[is.na(toTab)] <- ""
    l <- vector("list", nrow(unique(toTab[1])))
    l[[1]] <- kableExtra::kable_classic_2(
      kableExtra::kbl(toTab[-c(1, ncol(toTab))], 
                      align = "c", 
                      caption = paste0("<span style='font-size:12px; font-weight: bold; font-style: italic'>", 
                                       m, 
                                       " ",
                                       paste(Ord, sep = ": "), 
                                       "</span>")),
      bootstrap_options = c("striped", "hover", "condensed", "responsive"),
      full_width = F,
      html_font = "arial",
      font_size = 10
    )
    
    BestModsTemp <- BestMods[which(BestMods$mods %in% unique(toTab[[1]])),]
    for (i in seq(nrow(unique(toTab[1])))) {
      l[[i + 1]] <- kableExtra::pack_rows(l[[i]],
                                          BestModsTemp[i, "modsId"],
                                          min(which(toTab[-ncol(toTab)][[1]] ==  BestModsTemp[i, "mods"])),
                                          max(which(toTab[-ncol(toTab)][[1]] ==  BestModsTemp[i, "mods"])),
                                          label_row_css = "background-color: #fff; color: #000; border-bottom: 1px solid; border-top: 1px solid; font-style: italic;")
    }
    l[[length(l)]]
    kableExtra::save_kable(l[[length(l)]], 
                           file.path(DirSave, "Signif_Effects", "BetaDispFunctionsBestModTab", 
                                     paste0(paste("BetaDispFunctionsBestModTab", m, Ord, sep = "_"), 
                                            ".HTML")))
  }
}

# generate the table including all metric for the original data log transformed (Table ...)
## retrieve the summary from the distances of the original Log transformed data in PCOA ordination (for both axes)
dataNorm <- "original.log"
ordination <- "PCOA"

ValuesPart <- VRList[[ordination]]
OrLog <- lapply(ValuesPart, function(x)
  x[which(x[[1]] == dataNorm), ])
OrLogDf <- do.call(rbind,OrLog)
OrLogDf$index <- gsub("\\.\\d+", "", rownames(OrLogDf))
rownames(OrLogDf) <- NULL
OrLogDf[is.na(OrLogDf)] <- ""

l <- vector("list", length(unique(OrLogDf[["index"]])))
l[[1]] <- kableExtra::kable_classic_2(
  kableExtra::kbl(
    OrLogDf[-c(1, length(OrLogDf),length(OrLogDf)-1 )],
    align = "c",
    caption = paste0(
      "<span style='font-size:12px; font-weight: bold; font-style: italic'>",
      paste(dataNorm,
            ordination, sep = " "),
      "</span>"
    )
  ),
  bootstrap_options = c("striped", "hover", "condensed", "responsive"),
  full_width = F,
  html_font = "arial",
  font_size = 10
)

for (i in seq_along(distances)) {
  l[[i + 1]] <- kableExtra::pack_rows(l[[i]],
                                      distances[i], 
                                      min(which(OrLogDf[["index"]] == distances[i])),
                                      max(which(OrLogDf[["index"]] == distances[i])),
                                      label_row_css = "background-color: #fff; color: #000; border-bottom: 1px solid; border-top: 1px solid; font-style: italic;")
}
kableExtra::save_kable(l[[length(l)]], 
                       file.path(DirSave, "Signif_Effects", "BetaDispFunctionsBestModTab", 
                                 paste0(paste("BetaDispFunctionsBestModTab", dataNorm, ordination, sep = "_"), 
                                        ".HTML")))

# compute mean distance to median and se for some indices that give significant differences 
## For "Robust.Aitchison" in contamination treatment
df <- data.frame(dist = betadispRes[[dataNorm]][["Robust.Aitchison"]]$Contam$betaDisper$distances,
                 group = betadispRes[[dataNorm]][["Robust.Aitchison"]]$Contam$betaDisper$group)

meanRA <- stats::aggregate(df$dist, list(df$group), FUN = function(x) mean(x, na.rm = T))
seRA <- stats::aggregate(df$dist, list(df$group), FUN = function(x) sd(x, na.rm = T))$x /
  stats::aggregate(df$dist, list(df$group), function(x)
    sqrt(sum(!is.na(x))))$x

## For "Robust.Aitchison" according to sex
df <- data.frame(dist = betadispRes[[dataNorm]][["Robust.Aitchison"]]$Sexe$betaDisper$distances,
                 group = betadispRes[[dataNorm]][["Robust.Aitchison"]]$Sexe$betaDisper$group)

meanRA <- stats::aggregate(df$dist, list(df$group), FUN = function(x) mean(x, na.rm = T))
seRA <- stats::aggregate(df$dist, list(df$group), FUN = function(x) sd(x, na.rm = T))$x /
  stats::aggregate(df$dist, list(df$group), function(x)
    sqrt(sum(!is.na(x))))$x

p <- betadispRes[[dataNorm]][["Robust.Aitchison"]]$Sexe$FDRadjustedP


## For "Bray" according to origin
df <- data.frame(dist = betadispRes[[dataNorm]][["Bray"]]$ContPop$betaDisper$distances,
                 group = betadispRes[[dataNorm]][["Bray"]]$ContPop$betaDisper$group)

meanRA <- stats::aggregate(df$dist, list(df$group), FUN = function(x) mean(x, na.rm = T))
seRA <- stats::aggregate(df$dist, list(df$group), FUN = function(x) sd(x, na.rm = T))$x /
  stats::aggregate(df$dist, list(df$group), function(x)
    sqrt(sum(!is.na(x))))$x

p <- betadispRes[[dataNorm]][["Bray"]]$ContPop$FDRadjustedP
p

## For "jaccard" according to origin
df <- data.frame(dist = betadispRes[[dataNorm]][["Jaccard"]]$ContPop$betaDisper$distances,
                 group = betadispRes[[dataNorm]][["Jaccard"]]$ContPop$betaDisper$group)

meanRA <- stats::aggregate(df$dist, list(df$group), FUN = function(x) mean(x, na.rm = T))
seRA <- stats::aggregate(df$dist, list(df$group), FUN = function(x) sd(x, na.rm = T))$x /
  stats::aggregate(df$dist, list(df$group), function(x)
    sqrt(sum(!is.na(x))))$x

p <- betadispRes[[dataNorm]][["Jaccard"]]$ContPop$FDRadjustedP
p
