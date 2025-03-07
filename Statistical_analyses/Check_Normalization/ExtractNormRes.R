# Script Title: Consistency Check among Normalization Methods and Ordination Approaches (PCOA vs NMDS) 
#              for Diversity Analyses in Fish Gut Microbiota
#
# Author: Quentin PETITJEAN
# Date Created: 03/2023
# Last Modified: 05/05/2024
# ==============================================================================
# Requirements:
# - R version 4.2.3
# - Packages:
#   - ggplot2 v3.4.0: For generating publication-quality plots.
# ==============================================================================
# Script Overview:
# This script checks whether the results are consistent among different normalization 
# methods and between ordination approaches (PCOA vs NMDS) applied to the fish gut microbiota 
# dataset. The analysis is conducted over both alpha and beta diversity measures using 
# linear mixed models (LMM), Envfit, and betadisper frameworks. The workflow includes:
#
# 1. Specifying file paths and importing experimental design data and precomputed analysis 
#    outputs (LMM, Envfit, and Betadisper results).
# 2. Extracting and grouping results from LMM analyses for alpha diversity and beta diversity 
#    (from both PCOA and NMDS ordinations) into summary data frames.
# 3. Plotting treatment and covariate effects (with significance indicators) on diversity 
#    measures using bar plots.
# 4. Comparing the outcomes across normalization methods and between the two ordination 
#    techniques (PCOA vs NMDS) to assess result consistency.
# 5. Generating and saving effect plots (in SVG format) for further 
#    interpretation and reporting.
#
# ==============================================================================
# Usage:
# 1. Update the 'savingDir' variable to point to the directory containing both the input 
#    and output files.
# 2. Ensure that all required input files (e.g., LMMAlphaDivRes.rds, LMMRes.rds, EnvfitRes.rds, 
#    BetaDispRes.rds) are present in the specified directories.
# 3. Install and load the required R packages if they are not already installed.
# 4. Run the script in an R environment to generate effect plots.
# 5. The outputs include SVG plots displaying effect sizes and significance for various diversity indices.
# ==============================================================================

##############################################
#       	Install needed packages            #
##############################################
if(!require(ggplot2)){
  install.packages("ggplot2")
}

#######################
# specify some path   #
#######################
# the path to the directory from where file should be both imported and saved
savingDir <- "W:/POSTDOC_INP_GOLFECH_2023/Outputs/Normalized_Data/Signif_Effects"

#############################################################################################################
# Effect of treatments and covariates on alpha diversity tested within LMM framework                        #
#############################################################################################################

## extract the results of LMM conducted on PCOA and NMDS axes (beta-diversity)
AlphaLMMRes <- readRDS(file.path(savingDir, "LMMAlphaDivRes.rds"))

# set an empty dataframe to retrieve the desired results
res <- stats::setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("Chisq", "p", "treat", "norm", "Index"))

for(i in names(AlphaLMMRes)){
  tempDat <- AlphaLMMRes[[i]]
  for(j in names(tempDat)){
    tempDat2 <- tempDat[[j]]
      
      # group the results in a dataframe
      if(!length(tempDat2$ModAnov[["Chisq"]][-1]) == 0){
        tempRes <- stats::setNames(data.frame(tempDat2$ModAnov[["Chisq"]][-1],
                                              tempDat2$ModAnov[["Pr(>Chisq)"]][-1],
                                              rownames(tempDat2$ModAnov)[-1],
                                              i,
                                              j
        ), names(res))
      } else { 
        tempRes <- stats::setNames(data.frame(NA,
                                              NA,
                                              NA,
                                              i,
                                              j
        ), names(res))
      }
      
      for(m in c("Contam", "Inj", "ContPop", "Contam:Inj", "Contam:ContPop", "ContPop:Inj", "Contam:ContPop:Inj", "Sexe", "Taille_mm_M")){
        if(!m %in% tempRes[[3]]){
          tempRes <- rbind(tempRes, stats::setNames(data.frame(NA, NA, m, tempRes[1, c(4:5)]), names(res)))
          rownames(tempRes) <- NULL
        }
      }
    if(is.na(tempRes[["treat"]][1])){
      tempRes <- tempRes[-is.na(tempRes[["treat"]]),]
    }
    res <- rbind(res, tempRes)
  }
}

# plot the effects 
library(gridExtra)
for(i in seq_along(unique(res[["Index"]]))){
  plotList <-list()
    Temp <-
      res[which(res[["Index"]] == unique(res[["Index"]])[i]), ]
    if(nrow(Temp) == 0){
      next
    }
    # add significance stars according to pvalue and remove the NA for the plot
    #Temp <- Temp[!is.na(Temp$r), ]
    Temp$signif <- ifelse(Temp$p < 0.05, "*", "")
    norm_positions <- as.numeric(as.factor(Temp$norm))
    # Calculate mid-points
    mid_points <- unique(norm_positions) + 0.5
    # modify the legend
    leg <- unique(Temp$treat)[!is.na(unique(Temp$treat))]
    leg <- leg[order(leg)]
    
    PredictorsNames <- data.frame(toreplace = leg, replacement = NA)
    PredictorsNames$replacement[PredictorsNames$toreplace == "Inj"] <- "Imm. Chall."
    PredictorsNames$replacement[PredictorsNames$toreplace == "Sexe"] <- "Sex"
    PredictorsNames$replacement[PredictorsNames$toreplace == "Taille_mm_M"] <- "Size (cm)"
    PredictorsNames$replacement[PredictorsNames$toreplace == "Contam"] <- "Contamination"
    PredictorsNames$replacement[PredictorsNames$toreplace == "ContPop"] <- "Origin"
    PredictorsNames$replacement[PredictorsNames$toreplace == "Contam:ContPop"] <- "Contam.:Origin"
    PredictorsNames$replacement[PredictorsNames$toreplace == "Contam:Inj"] <- "Origin:Imm. Chall."
    PredictorsNames$replacement[PredictorsNames$toreplace == "ContPop:Inj"] <- "Contam.:Imm. Chall."
    PredictorsNames$replacement[PredictorsNames$toreplace == "Contam:ContPop:Inj"] <- "Contam.:Origin:Imm. Chall."
    PredictorsNames <- PredictorsNames[!is.na(PredictorsNames$replacement),]
    
    
      p <- ggplot2::ggplot(Temp, ggplot2::aes(fill=treat, y=Chisq, x=norm)) + 
      ggplot2::geom_bar(position=ggplot2::position_dodge(width = 0.7), stat="identity", width = 0.7, colour="black") + 
      ggplot2::ggtitle(paste(unique(res[["dist"]])[i], unique(res[["coordSys"]])[j], sep = " ")) +
      ggplot2::xlab("") +
      ggplot2::theme_classic()+ 
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::scale_fill_brewer(name="Treatment/Covariate", palette="Dark2", breaks = PredictorsNames$toreplace, labels = PredictorsNames$replacement) +
      ggplot2::geom_text(ggplot2::aes(label = signif, y = Chisq + Chisq/10),  # Adjust 0.5 if you need more space between bars and labels
                         vjust = -0.5, position = ggplot2::position_dodge(width = 0.7))+
      ggplot2::geom_vline(xintercept = mid_points, linetype="dashed", color = "grey30")
    
  # create output folder to save figs
  if (length(list.dirs(file.path(savingDir, "Signif_Plots"))) == 0) {
    dir.create(file.path(savingDir, "Signif_Plots"))
  }
  if (length(list.dirs(file.path(savingDir, "Signif_Plots", "EffectPlot"))) == 0) {
    dir.create(file.path(savingDir, "Signif_Plots", "EffectPlot"))
  }
  if (length(list.dirs(file.path(savingDir, "Signif_Plots", "EffectPlot", "AlphaDivLMM"))) == 0) {
    dir.create(file.path(savingDir, "Signif_Plots", "EffectPlot", "AlphaDivLMM"))
  }
  
  # save the plots' panel
  ggplot2::ggsave(
    file = paste0(file.path(savingDir,"Signif_Plots", "EffectPlot", "AlphaDivLMM", unique(res[["Index"]])[i]), ".svg"),
    plot = p,
    width = 9,
    height = 6
  )
}

####################################################################################################
# Effect of treatments and covariates on Beta diversity tested within LMM framework                #
####################################################################################################

### extract the results of LMM conducted on PCOA and NMDS axes (beta-diversity)
NormRes <- readRDS(file.path(savingDir, "LMMRes.rds"))

# set an empty dataframe to retrieve the desired results
res <- stats::setNames(data.frame(matrix(ncol = 7, nrow = 0)), c("Chisq", "p", "treat", "norm", "dist", "coordSys", "axis"))

for(i in names(NormRes)){
  tempDat <- NormRes[[i]]
  for(j in names(tempDat)){
      tempDat2 <- tempDat[[j]]
    for(k in names(tempDat2)){
        tempDat3 <- tempDat2[[k]]
      for(l in names(tempDat3)){
        tempDat4 <- tempDat3[[l]]
        
        # group the results in a dataframe
        if(!length(tempDat4[["ModAnov"]]$Chisq[-1]) == 0){
        tempRes <- stats::setNames(data.frame(tempDat4[["ModAnov"]]$Chisq[-1],
                     tempDat4[["ModAnov"]]$'Pr(>Chisq)'[-1],
                     rownames(tempDat4[["ModAnov"]])[-1],
                     i,
                     j,
                     k,
                     l
        ), names(res))
        } else { 
          tempRes <- stats::setNames(data.frame(NA,
                                                NA,
                                                NA,
                                                i,
                                                j,
                                                k,
                                                l
          ), names(res))
          }

        for(m in c("Contam", "Inj", "ContPop", "Contam:Inj", "Contam:ContPop", "ContPop:Inj", "Contam:ContPop:Inj", "Sexe", "Taille_mm_M")){
          if(!m %in% tempRes[[3]]){
            tempRes <- rbind(tempRes, stats::setNames(data.frame(NA, NA, m, tempRes[1, c(4:7)]), names(res)))
            
          }
        }
        res <- rbind(res, tempRes)
      }
    }
  }
}

# plot the effects (chisq are plotted when p is significant)
library(gridExtra)
for(i in seq_along(unique(res[["dist"]]))){
  plotList <-list()
  for(j in seq_along(unique(res[["coordSys"]]))){
    for(k in seq_along(unique(res[["axis"]]))){
      Temp <-
        res[which(res[["dist"]] == unique(res[["dist"]])[i] &
                    res[["coordSys"]] == unique(res[["coordSys"]])[j] &
                    res[["axis"]] == unique(res[["axis"]])[k]), ]
      if(nrow(Temp) == 0){
        next
      }
      # add significance stars according to pvalue and remove the NA for the plot
      #Temp <- Temp[!is.na(Temp$Chisq), ]
      Temp$signif <- ifelse(Temp$p < 0.05, "*", "")
      norm_positions <- as.numeric(as.factor(Temp$norm))
      # Calculate mid-points
      mid_points <- unique(norm_positions) + 0.5
      leg <- unique(Temp$treat)[!is.na(unique(Temp$treat))]
      leg <- leg[order(leg)]
      
      PredictorsNames <- data.frame(toreplace = leg, replacement = NA)
      PredictorsNames$replacement[PredictorsNames$toreplace == "Inj"] <- "Imm. Chall."
      PredictorsNames$replacement[PredictorsNames$toreplace == "Sexe"] <- "Sex"
      PredictorsNames$replacement[PredictorsNames$toreplace == "Taille_mm_M"] <- "Size (cm)"
      PredictorsNames$replacement[PredictorsNames$toreplace == "Contam"] <- "Contamination"
      PredictorsNames$replacement[PredictorsNames$toreplace == "ContPop"] <- "Origin"
      PredictorsNames$replacement[PredictorsNames$toreplace == "Contam:ContPop"] <- "Contam.:Origin"
      PredictorsNames$replacement[PredictorsNames$toreplace == "Contam:Inj"] <- "Origin:Imm. Chall."
      PredictorsNames$replacement[PredictorsNames$toreplace == "ContPop:Inj"] <- "Contam.:Imm. Chall."
      PredictorsNames$replacement[PredictorsNames$toreplace == "Contam:ContPop:Inj"] <- "Contam.:Origin:Imm. Chall."
      PredictorsNames <- PredictorsNames[!is.na(PredictorsNames$replacement),]
      
      plotList[[k]] <- ggplot2::ggplot(Temp, ggplot2::aes(fill=treat, y=Chisq, x=norm)) + 
        ggplot2::geom_bar(position=ggplot2::position_dodge(width = 0.7), stat="identity", width = 0.7, colour="black") + 
        ggplot2::ggtitle(paste(unique(res[["dist"]])[i], unique(res[["coordSys"]])[j], unique(res[["axis"]])[k], sep = " ")) +
        ggplot2::xlab("") +
        ggplot2::theme_classic()+ 
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::scale_fill_brewer(name="Treatment/Covariate", palette="Dark2", breaks = PredictorsNames$toreplace, labels = PredictorsNames$replacement) +
        ggplot2::geom_text(ggplot2::aes(label = signif, y = Chisq + 0.5),  # Adjust 0.5 if you need more space between bars and labels
                           vjust = -0.5, position = ggplot2::position_dodge(width = 0.7)) +
        ggplot2::geom_vline(xintercept = mid_points, linetype="dashed", color = "grey30")
      }
  }
  n <- length(plotList)
  nCol <- floor(sqrt(n))
  nRow <- ifelse(nCol == 1 & n > 1, 2, nCol)
  panelPlot <- do.call(ggpubr::ggarrange, c(plotList, ncol = nCol, nrow=nRow, common.legend = TRUE, legend="bottom"))
  
  # create output folder to save figs
  if (length(list.dirs(file.path(savingDir, "Signif_Plots"))) == 0) {
    dir.create(file.path(savingDir, "Signif_Plots"))
  }
  if (length(list.dirs(file.path(savingDir, "Signif_Plots", "EffectPlot"))) == 0) {
    dir.create(file.path(savingDir, "Signif_Plots", "EffectPlot"))
  }
  if (length(list.dirs(file.path(savingDir, "Signif_Plots", "EffectPlot", "BetaDivLMM"))) == 0) {
    dir.create(file.path(savingDir, "Signif_Plots", "EffectPlot", "BetaDivLMM"))
  }
  
  # save the plots' panel
  ggplot2::ggsave(
    file = paste0(file.path(savingDir,"Signif_Plots", "EffectPlot", "BetaDivLMM", unique(res[["dist"]])[i]), ".svg"),
    plot = panelPlot,
    width = 16,
    height = 12
  )
}

#######################################################################################################
# Effect of treatments and covariates on Beta diversity tested within Envfit framework                #
#######################################################################################################

### do the same for Envfit (correlation between ordination score and treatment and covariates)
## extract the results of LMM conducted on PCOA and NMDS axes (beta-diversity)
EnvFitRes <- readRDS(file.path(savingDir, "EnvfitRes.rds"))

# set an empty dataframe to retrieve the desired results
resEnvfit <- stats::setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("r", "p", "treat", "norm", "dist", "coordSys"))

for(i in names(EnvFitRes)){
  tempDat <- EnvFitRes[[i]]
  for(j in names(tempDat)){
    tempDat2 <- tempDat[[j]]
    for(k in names(tempDat2)){
      tempDat3 <- tempDat2[[k]]
      
        # group the results in a dataframe
      if(!length(tempDat3[["vectors"]]$r) == 0){
        tempRes <- stats::setNames(data.frame(c(tempDat3[["vectors"]]$r, tempDat3[["factors"]]$r),
                                              c(tempDat3[["vectors"]]$pvals, tempDat3[["factors"]]$pvals),
                                              c(names(tempDat3[["vectors"]]$r), names(tempDat3[["factors"]]$r)),
                                              i,
                                              j,
                                              k
        ), names(resEnvfit))
        } else { 
          tempRes <- stats::setNames(data.frame(NA,
                                                NA,
                                                NA,
                                                i,
                                                j,
                                                k
          ), names(res))
        }
        
        for(m in c("Contam", "Inj", "ContPop", "Contam:Inj", "Contam:ContPop", "ContPop:Inj", "Contam:ContPop:Inj", "Sexe", "Taille_mm_M")){
          if(!m %in% tempRes[[3]]){
            tempRes <- rbind(tempRes, stats::setNames(data.frame(NA, NA, m, tempRes[1, c(4:6)]), names(resEnvfit)))
            rownames(tempRes) <- NULL
          }
        }
        resEnvfit <- rbind(resEnvfit, tempRes)
    }
  }
}

# plot the effects (chisq are plotted when p is significant)
library(gridExtra)
for(i in seq_along(unique(resEnvfit[["dist"]]))){
  plotList <-list()
  for(j in seq_along(unique(resEnvfit[["coordSys"]]))){
      Temp <-
        resEnvfit[which(resEnvfit[["dist"]] == unique(resEnvfit[["dist"]])[i] &
                    resEnvfit[["coordSys"]] == unique(resEnvfit[["coordSys"]])[j]), ]
      if(nrow(Temp) == 0){
        next
      }
      # add significance stars according to pvalue and remove the NA for the plot
      #Temp <- Temp[!is.na(Temp$r), ]
      Temp$signif <- ifelse(Temp$p < 0.05, "*", "")
      norm_positions <- as.numeric(as.factor(Temp$norm))
      # Calculate mid-points
      mid_points <- unique(norm_positions) + 0.5
      # modify the legend
      leg <- unique(Temp$treat)[!is.na(unique(Temp$treat))]
      leg <- leg[order(leg)]
      
      PredictorsNames <- data.frame(toreplace = leg, replacement = NA)
      PredictorsNames$replacement[PredictorsNames$toreplace == "Inj"] <- "Imm. Chall."
      PredictorsNames$replacement[PredictorsNames$toreplace == "Sexe"] <- "Sex"
      PredictorsNames$replacement[PredictorsNames$toreplace == "Taille_mm_M"] <- "Size (cm)"
      PredictorsNames$replacement[PredictorsNames$toreplace == "Contam"] <- "Contamination"
      PredictorsNames$replacement[PredictorsNames$toreplace == "ContPop"] <- "Origin"
      PredictorsNames$replacement[PredictorsNames$toreplace == "Contam:ContPop"] <- "Contam.:Origin"
      PredictorsNames$replacement[PredictorsNames$toreplace == "Contam:Inj"] <- "Origin:Imm. Chall."
      PredictorsNames$replacement[PredictorsNames$toreplace == "ContPop:Inj"] <- "Contam.:Imm. Chall."
      PredictorsNames$replacement[PredictorsNames$toreplace == "Contam:ContPop:Inj"] <- "Contam.:Origin:Imm. Chall."
      PredictorsNames <- PredictorsNames[!is.na(PredictorsNames$replacement),]
      
      
      plotList[[j]] <- ggplot2::ggplot(Temp, ggplot2::aes(fill=treat, y=r, x=norm)) + 
        ggplot2::geom_bar(position=ggplot2::position_dodge(width = 0.7), stat="identity", width = 0.7, colour="black") + 
        ggplot2::ggtitle(paste(unique(resEnvfit[["dist"]])[i], unique(resEnvfit[["coordSys"]])[j], sep = " ")) +
        ggplot2::xlab("") +
        ggplot2::theme_classic()+ 
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::scale_fill_brewer(name="Treatment/Covariate", palette="Dark2", breaks = PredictorsNames$toreplace, labels = PredictorsNames$replacement) +
        ggplot2::geom_text(ggplot2::aes(label = signif, y = r + r/10),  # Adjust 0.5 if you need more space between bars and labels
                           vjust = -0.5, position = ggplot2::position_dodge(width = 0.7))+
        ggplot2::geom_vline(xintercept = mid_points, linetype="dashed", color = "grey30")
  }
  n <- length(plotList)
  nCol <- floor(sqrt(n))
  nRow <- ifelse(nCol == 1 & n > 1, 2, nCol)
  panelPlot <- do.call(ggpubr::ggarrange, c(plotList, ncol = nCol, nrow=nRow, common.legend = TRUE, legend="bottom"))
  
  # create output folder to save figs
  if (length(list.dirs(file.path(savingDir, "Signif_Plots"))) == 0) {
    dir.create(file.path(savingDir, "Signif_Plots"))
  }
  if (length(list.dirs(file.path(savingDir, "Signif_Plots", "EffectPlot"))) == 0) {
    dir.create(file.path(savingDir, "Signif_Plots", "EffectPlot"))
  }
  if (length(list.dirs(file.path(savingDir, "Signif_Plots", "EffectPlot", "EnvFit"))) == 0) {
    dir.create(file.path(savingDir, "Signif_Plots", "EffectPlot", "EnvFit"))
  }
  
  # save the plots' panel
  ggplot2::ggsave(
    file = paste0(file.path(savingDir,"Signif_Plots", "EffectPlot", "EnvFit", unique(resEnvfit[["dist"]])[i]), ".svg"),
    plot = panelPlot,
    width = 16,
    height = 12
  )
}

#############################################################################################################
# Effect of treatments and covariates on Beta diversity tested within Betadisper framework                  #
#############################################################################################################

### extract the results of BetaDISPER conducted on PCOA and NMDS axes
BetadispRes <- readRDS("D:/POSTDOC_INP_GOLFECH_2023/Outputs/Normalized_Data/Signif_Effects/BetaDispRes.rds")

res <- stats::setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("F", "Pr(>F)", "treat", "norm", "dist"))

for(i in names(BetadispRes)){
  tempDat <- BetadispRes[[i]]
  for(j in names(tempDat)){
    tempDat2 <- tempDat[[j]]
        FVal <- do.call("rbind", lapply(tempDat2, function(x) x[["PermTest"]]$tab["Groups", "F"]))
        pVal <- do.call("rbind", lapply(tempDat2, function(x) x[["PermTest"]]$tab["Groups", "Pr(>F)"]))

        # group the results in a dataframe
        if(!length(tempDat2) == 0){
        tempRes <- stats::setNames(data.frame(FVal[,1],
                                              pVal[,1],
                                              rownames(FVal),
                                              i,
                                              j
        ), names(res))  
        } else { 
          tempRes <- stats::setNames(data.frame(NA,
                                                NA,
                                                NA,
                                                i,
                                                j
          ), names(res))
        }
        
        
        for(m in c("Contam", "Inj", "Pop", "Contam:Inj", "Contam:Pop", "Inj:Pop", "Sexe", "Taille_mm_M")){
          if(!m %in% tempRes[[3]]){
            tempRes <- rbind(tempRes, stats::setNames(data.frame(NA, NA, m, tempRes[1, c(4:5)]), names(res)))
            rownames(tempRes) <- NULL
          }
        }
        res <- rbind(res, tempRes)
  }
}

# plot the effects (chisq are plotted when p is significant)
library(gridExtra)
for(i in seq_along(unique(res[["dist"]]))){
    plotList <-list()
      Temp <-
        res[which(res[["dist"]] == unique(res[["dist"]])[i]), ]
      if(nrow(Temp) == 0){
        next
      }
      
      # add significance stars according to pvalue and remove the NA for the plot
      #Temp <- Temp[!is.na(Temp$r), ]
      Temp$signif <- ifelse(Temp[["Pr(>F)"]] < 0.05, "*", "")
      norm_positions <- as.numeric(as.factor(Temp$norm))
      # Calculate mid-points
      mid_points <- unique(norm_positions) + 0.5
      # modify the legend
      leg <- unique(Temp$treat)[!is.na(unique(Temp$treat))]
      leg <- leg[order(leg)]
      
      PredictorsNames <- data.frame(toreplace = leg, replacement = NA)
      PredictorsNames$replacement[PredictorsNames$toreplace == "Inj"] <- "Imm. Chall."
      PredictorsNames$replacement[PredictorsNames$toreplace == "Sexe"] <- "Sex"
      PredictorsNames$replacement[PredictorsNames$toreplace == "Taille_mm_M"] <- "Size (cm)"
      PredictorsNames$replacement[PredictorsNames$toreplace == "Contam"] <- "Contamination"
      PredictorsNames$replacement[PredictorsNames$toreplace == "ContPop"] <- "Origin"
      PredictorsNames$replacement[PredictorsNames$toreplace == "Contam:ContPop"] <- "Contam.:Origin"
      PredictorsNames$replacement[PredictorsNames$toreplace == "Contam:Inj"] <- "Origin:Imm. Chall."
      PredictorsNames$replacement[PredictorsNames$toreplace == "ContPop:Inj"] <- "Contam.:Imm. Chall."
      PredictorsNames$replacement[PredictorsNames$toreplace == "Contam:ContPop:Inj"] <- "Contam.:Origin:Imm. Chall."
      PredictorsNames <- PredictorsNames[!is.na(PredictorsNames$replacement),]
 
      p <-  ggplot2::ggplot(Temp, ggplot2::aes(fill=treat, y=F, x=norm)) + 
        ggplot2::geom_bar(position=ggplot2::position_dodge(width = 0.7), stat="identity", width = 0.7, colour="black") + 
        ggplot2::ggtitle(unique(res[["dist"]])[i]) +
        ggplot2::xlab("") +
        ggplot2::theme_classic()+ 
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
        ggplot2::scale_fill_brewer(name="Treatment/Covariate", palette="Dark2", breaks = PredictorsNames$toreplace, labels = PredictorsNames$replacement) +
        ggplot2::geom_text(ggplot2::aes(label = signif, y = F + F/10),  # Adjust 0.5 if you need more space between bars and labels
                           vjust = -0.5, position = ggplot2::position_dodge(width = 0.7))+
        ggplot2::geom_vline(xintercept = mid_points, linetype="dashed", color = "grey30")

      
      # create output folder to save figs
      if (length(list.dirs(file.path(savingDir, "Signif_Plots"))) == 0) {
        dir.create(file.path(savingDir, "Signif_Plots"))
      }
      if (length(list.dirs(file.path(savingDir, "Signif_Plots", "EffectPlot"))) == 0) {
        dir.create(file.path(savingDir, "Signif_Plots", "EffectPlot"))
      }
      if (length(list.dirs(file.path(savingDir, "Signif_Plots", "EffectPlot", "BetadisperPermut"))) == 0) {
        dir.create(file.path(savingDir, "Signif_Plots", "EffectPlot", "BetadisperPermut"))
      }
      
      # save the plots' panel
      ggplot2::ggsave(
        file = file.path(savingDir, "Signif_Plots", "EffectPlot", "BetadisperPermut", paste0(unique(res[["dist"]])[i], ".svg")),
        plot = p,
        width = 9,
        height = 6
      )
}
