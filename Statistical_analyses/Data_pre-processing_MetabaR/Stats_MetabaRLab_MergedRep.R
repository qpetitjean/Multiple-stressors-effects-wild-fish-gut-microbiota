##############################################
#       	Data pre-processing                #
##############################################

#######################
# specify some path   #
#######################
# the path to the directory from where file should be both imported and saved
savingDir <- "D:/POSTDOC_INP_GOLFECH_2023/Outputs/Preprocessing-Metabar"

#######################
# Import data         #
#######################
# import the dataset
labData <- read.table(file.path(savingDir, "Phypat1_2_3_R1R2_concatenated_ngsfilt_n50_ali_uniq_noS_cl97_silva_138.1_ssuparc_full_tagged_merged_fullTaxo_LabOnly.tab"), 
                          sep="\t", dec= ".")

####################################################################
# Data curation according to Bokulich et al 2013                   #
####################################################################

  ## Remove MOTUs with relative abundance below 0.00005
  totreads <- sum(labData$count)
  labData$MOTUrelative_abund <- labData$count/totreads
  labDataClean <- labData[-which(labData$MOTUrelative_abund < 0.00005),]
  
  ## Remove MOTUs detected in less than 3 samples
  labDataClean$MOTUoccurence <- rowSums(labDataClean[,1:(grep("count", colnames(labDataClean))-1)]>0)
  if(length(which(labDataClean$MOTUoccurence < 3)) > 0){
  labDataClean <- labDataClean[-which(labDataClean$MOTUoccurence < 3),]
  }
  labDataClean <- labDataClean[order(labDataClean$MOTUoccurence, decreasing = TRUE),]

# generate new phylum, class ... columns accordingly
taxo <- lapply(labDataClean$path, function(x) {
  matchPos <- gregexpr("(?<=:)(.*?)(?=@)", x, perl = TRUE)
  results <- unlist(regmatches(x, matchPos))
  matchPosCol <- gregexpr("(?<=@)(.*?)(?=:)|(?<=@)(.*?)(?=$)", x, perl = TRUE)
  col <- unlist(regmatches(x, matchPosCol))
  col <- col[-which(col == "root")]
  df <- stats::setNames(as.data.frame(t(results)), col)
  return(df)
})
taxodf <- do.call(plyr::rbind.fill, taxo)

# merge them to the full data
labDataClean <- cbind(labDataClean, taxodf)

# prepare the motus table for metabaR (retrieve only MOTUs information without read counts)
charcol <- labDataClean[,grep("count", colnames(labDataClean)):ncol(labDataClean)]

# prepare the reads table for metabaR
labDataCleanFinal <- as.data.frame(t(labDataClean[,-c(grep("count", colnames(labDataClean)):ncol(labDataClean))]))
rownames(labDataCleanFinal) <-  gsub("^X", "", rownames(labDataCleanFinal))

# save them on the hardrive
write.table(labDataCleanFinal, file.path(savingDir, "file_motus_P1P2P3_BokulFilt_MergedRep.txt"), sep="\t", dec= ".")
write.table(charcol, file.path(savingDir, "file_motus_P1P2P3_BokulFiltTaxo_MergedRep.txt"), sep="\t", dec= ".")

##########################
# Import metabar list    #
##########################

# Import libraries needed for the following part
library(metabaR)
library(ggplot2)
library(reshape2)
library(cowplot)
library(RColorBrewer)

fguts_Bact <- metabaR::tabfiles_to_metabarlist(file_reads = file.path(savingDir, "file_motus_P1P2P3_BokulFilt_MergedRep.txt"),
                                               file_motus = file.path(savingDir, "file_motus_P1P2P3_BokulFiltTaxo_MergedRep.txt"),
                                               file_pcrs =  file.path(savingDir, "fishgut_Bact_pcrs_final.txt"),
                                               file_samples = file.path(savingDir, "fishgut_Bact_samples_final.txt"),
                                               files_sep = "\t")

metabaR::summary_metabarlist(fguts_Bact)

# some sample have 0 reads after the first filtration step (Bokulich)
## identify and remove them
### Compute the number of reads per pcr
fguts_Bact$pcrs$nb_reads <- rowSums(fguts_Bact$reads)

### Compute the number of motus per pcr
fguts_Bact$pcrs$nb_motus <- rowSums(fguts_Bact$reads>0)

### identify the pcrs with 0 reads (excepted the blank - BLC)
toRemove <- rownames(fguts_Bact$pcrs[fguts_Bact$pcrs$nb_motus == 0 & !grepl("BLC", rownames(fguts_Bact$pcrs)),])

### remove those pcrs from the dataset
fguts_Bact <- subset_metabarlist(fguts_Bact, 
                                 table = "pcrs", 
                                 indices = !rownames(fguts_Bact$pcrs) %in% toRemove)

metabaR::summary_metabarlist(fguts_Bact)

###############################
#       Diagnostic plot       #
###############################

# Create an input table (named check1) for ggplot of 3 columns: 
#  (i) control type 
#  (ii) a vector indicated whether it corresponds to nb_reads or nb_motus, 
#  (iii) the corresponding values.

check1 <- reshape2::melt(fguts_Bact$pcrs[,c("control_type", "nb_reads", "nb_motus")])
head(check1)
ggplot(data <- check1, aes(x=control_type, y=value, color=control_type)) + 
  geom_boxplot() + theme_bw() + 
  geom_jitter(alpha=0.2) + 
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  facet_wrap(~variable, scales = "free_y") + 
  theme(axis.text.x = element_text(angle=45, h=1))

# Using the nb_reads and nb_motus defined previously in the fguts_Bact$pcrs table
ggplot(fguts_Bact$pcrs, aes(x=nb_reads, y=nb_motus, color = control_type)) + 
  geom_point() + theme_bw() + 
  scale_y_log10() + scale_x_log10() + 
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey")
ggpcrplate(fguts_Bact, FUN = function(m){rowSums(m$reads)}, legend_title = "# of reads per PCR")

# Here the list of all tag/indices used in the experiment 
# is available in the column "tag_rev" of the fguts_Bact$pcrs table
tag.list <- as.character(unique(fguts_Bact$pcrs$tag_rev))
ggpcrtag(fguts_Bact, 
         legend_title = "# of reads per PCR", 
         FUN = function(m) {rowSums(m$reads)},
         taglist = tag.list) 

# performed rarefaction curves for the whole dataset (excluding control sample)
fguts_Bact_sample_Only <- subset_metabarlist(fguts_Bact, 
                                             table = "pcrs", 
                                             indices = fguts_Bact$pcrs$type == "sample")

fguts_Bact.raref = hill_rarefaction(fguts_Bact_sample_Only, nboot = 20, nsteps = 10)
tail(fguts_Bact.raref$hill_table)
dim(fguts_Bact.raref$hill_table)
gghill_rarefaction(fguts_Bact.raref) 

#### performed rarefactions curves for each sample type
# Define a vector containing the Material info for each pcrs 
fguts_Bact_sample_Only$samples$Material=paste(fguts_Bact_sample_Only$samples$experiment, fguts_Bact_sample_Only$samples$matrix, sep="_")
material <- fguts_Bact_sample_Only$samples$Material[match(fguts_Bact_sample_Only$pcrs$sample_id,
                                                          rownames(fguts_Bact_sample_Only$samples))]

# Use of gghill_rarefaction requires a vector with named pcrs
material <- setNames(material,rownames(fguts_Bact_sample_Only$pcrs))

p <- gghill_rarefaction(fguts_Bact.raref, group=material)
p + scale_fill_manual(values = c("goldenrod4", "blue4")) +
  scale_color_manual(values = c("goldenrod4", "blue4")) +
  labs(color="Material type")

# no surprise, water samples tend to be more diverse than fish gut samples

###############################
#  Flagging spurious signal   #
###############################

### ... during extraction step

# Identifying extraction contaminants
fguts_Bact <- contaslayer(fguts_Bact, 
                          control_types = "extraction",
                          output_col = "not_an_extraction_conta",
                          method = "all")

# indicating whether the MOTU is a genuine MOTU (TRUE) or a contaminant (FALSE).
table(fguts_Bact$motus$not_an_extraction_conta) 

#get the list of MOTU considered as contaminants in during extraction
Contam_extract_MOTU <- subset(fguts_Bact$motus, not_an_extraction_conta =="FALSE")  

# Identify the most common contaminant
# get contaminant ids
fguts_Bact$motus$count=colSums(fguts_Bact$reads)
id <- !fguts_Bact$motus$not_an_extraction_conta
max.conta <- rownames(fguts_Bact$motus[id,])[which.max(fguts_Bact$motus[id, "count"])]
max.conta

#... and its distribution and relative abundance in each pcr
ggpcrplate(fguts_Bact, legend_title = "#reads of most \nabundant contaminant",
           FUN = function(m) {m$reads[, max.conta]/rowSums(m$reads)})

# do the same for the other contaminants 
contamId <- rownames(fguts_Bact$motus[id,])
MeanRelAb <- c()
MeanRelAb <- lapply(contamId, function(x){
  p <- ggpcrplate(fguts_Bact, legend_title = "#reads of most \nabundant contaminant",
                  FUN = function(m) {m$reads[, x]/rowSums(m$reads)})
  temp <- mean(fguts_Bact$reads[, x]/rowSums(fguts_Bact$reads), na.rm = T)
  MeanRelAb <- c(MeanRelAb, temp)
  p <- p + scale_x_continuous(sec.axis = sec_axis(~ . , name = paste(x, temp, sep = " "), breaks = NULL, labels = NULL))
  print(p)
  return(MeanRelAb)
})

MotusRelAbContam <- data.frame(contamId, MeanRelAb = unlist(MeanRelAb))

Contam_extract_MOTU_final <- subset(fguts_Bact$motus, not_an_extraction_conta =="FALSE")  
write.table(table(Contam_extract_MOTU$path), file.path(savingDir, "Extract_Contam.txt"), sep=";")

# Compute relative abundance of all pcr contaminants together 
ContamRelAb <- data.frame(conta.relab = rowSums(fguts_Bact$reads[,!fguts_Bact$motus$not_an_extraction_conta]) / 
                            rowSums(fguts_Bact$reads))

# Add information on control types
ContamRelAb$control_type <- fguts_Bact$pcrs$control_type[match(rownames(ContamRelAb), rownames(fguts_Bact$pcrs))]

ggplot(ContamRelAb, aes(x=control_type, y=conta.relab, color=control_type)) + 
  geom_boxplot() + geom_jitter(alpha=0.5) +
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  labs(x=NULL, y="Prop. Reads (log10)") + 
  theme_bw() + 
  scale_y_log10()

# flag pcrs with total contaminant relative abundance > 10% of reads)
fguts_Bact$pcrs$low_contamination_level <- 
  ifelse(ContamRelAb$conta.relab[match(rownames(fguts_Bact$pcrs), rownames(ContamRelAb))]>1e-1,  F, T)
fguts_Bact$pcrs$low_contamination_level 
# Proportion of potentially functional (TRUE) vs. dysfunctional (FALSE) pcrs
# (controls included) based on this criterion
table(fguts_Bact$pcrs$low_contamination_level) / nrow(fguts_Bact$pcrs)

### ... during sequencing step

# Identifying sequencing contaminants
fguts_Bact <- contaslayer(fguts_Bact, 
                          control_types = "sequencing",
                          output_col = "not_a_sequencing_conta", 
                          method = "all")

table(fguts_Bact$motus$not_a_sequencing_conta)

#get the list of MOTU considered as contaminants in extraction
Contam_sequencing_MOTU <- subset(fguts_Bact$motus, not_a_sequencing_conta=="FALSE" )                                                                            
Contam_sequencing_MOTU

# Identify the most common contaminant
# get contaminant ids
names(fguts_Bact$motus)
id <- !fguts_Bact$motus$not_a_sequencing_conta

# there is no sequencing contaminant detected, most contamination seems to have occured during extraction step

# Intersection with extraction contaminant flags (not contaminant = T)
table(fguts_Bact$motus$not_an_extraction_conta,
      fguts_Bact$motus$not_a_sequencing_conta)

### identify MOTUs whose sequence is too dissimilar from references (here we used the minimum identification score for a given cluster)

# Plot the unweighted distribution of MOTUs similarity scores 
a <- 
  ggplot( fguts_Bact$motus, aes(x=best_identity.silva_min)) + 
  geom_histogram(color="grey", fill="white", bins=20) + 
  geom_vline(xintercept = 0.8, col="orange", lty=2) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  labs(x="% similarity against best match", y="# OTUs")

# Same for the weighted distribution
b <- 
  ggplot( fguts_Bact$motus, 
          aes(x=best_identity.silva_min, y = ..count.., weight = count)) + 
  geom_histogram(color="grey", fill="white", bins=20) + 
  geom_vline(xintercept = 0.8, col="orange", lty=2) + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  labs(x="% similarity against best match", y="# Reads")

# Combine plots into one
cowplot::ggdraw() + 
  cowplot::draw_plot(a, x=0, y=0, width = 0.5) + 
  cowplot::draw_plot(b, x=0.5, y=0, width = 0.5)

# Flag not degraded (TRUE) vs. potentially degraded sequences (FALSE)
fguts_Bact$motus$not_degraded <-
  ifelse(fguts_Bact$motus$best_identity.silva_min < 0.8, F, T)

# Proportion of each of these over total number of MOTUs
table(fguts_Bact$motus$not_degraded) / nrow(fguts_Bact$motus)

# there is not remaining degraded MOTUs

# Intersection with other flags
table(fguts_Bact$motus$not_an_extraction_conta, 
      fguts_Bact$motus$not_a_sequencing_conta,
      fguts_Bact$motus$not_degraded)

# main contamination come from the extraction step with 6 MOTUs detected as contaminants

#############################
#  Detecting PCR outliers   #
#############################

# sequencing depth (treshold value = 1000) 
ggplot(fguts_Bact$pcrs, aes(nb_reads)) +
  geom_histogram(bins=40, color="grey", fill="white") + 
  geom_vline(xintercept = 1e3, lty=2, color="orange") + # threshold
  scale_x_log10() + 
  labs(x="# Reads (with all OTUs and PCRs)", 
       y="# PCRs") +
  theme_bw() + 
  theme(panel.grid = element_blank())

# Flag pcrs with an acceptable sequencing depth (treshold value = 1000) (TRUE) or inacceptable one (FALSE)
fguts_Bact$pcrs$seqdepth_ok <- ifelse(fguts_Bact$pcrs$nb_reads < 1e3, F, T)

# Proportion of each of these over total number of pcrs, control excluded
table(fguts_Bact$pcrs$seqdepth_ok[fguts_Bact$pcrs$type=="sample"]) /
  nrow(fguts_Bact$pcrs[fguts_Bact$pcrs$type=="sample",])

###  assess PCR quality trough their reproducibility (among replicates)
# Subsetting the metabarlist to remove pcrs yielding no reads, as well as negative controls
fguts_Bact_sub <- subset_metabarlist(fguts_Bact, 
                                     table = "pcrs", 
                                     indices = fguts_Bact$pcrs$nb_reads>0 & (
                                       is.na(fguts_Bact$pcrs$control_type) |
                                         fguts_Bact$pcrs$control_type=="positive"))

# First visualization
comp1 = pcr_within_between(fguts_Bact_sub)
check_pcr_thresh(comp1)

# Flagging
fguts_Bact_sub <- pcrslayer(fguts_Bact_sub, output_col = "replicating_pcr", plot = F)

# Proportion of replicating pcrs (TRUE)
table(fguts_Bact_sub$pcrs$replicating_pcr) /
  nrow(fguts_Bact_sub$pcrs)

# Intersection with the sequencing depth criterion
table(fguts_Bact_sub$pcrs$seqdepth_ok, 
      fguts_Bact_sub$pcrs$replicating_pcr)

fguts_Bact_sub$pcrs[fguts_Bact_sub$pcrs$seqdepth_ok,]

# Distinguish between pcrs obtained from samples from positive controls
mds = check_pcr_repl(fguts_Bact_sub, 
                     groups = fguts_Bact_sub$pcrs$type, 
                     funcpcr = fguts_Bact_sub$pcrs$replicating_pcr)
mds + labs(color = "pcr type") + scale_color_manual(values = c("cyan4", "gray"))

# Now report the flagging in the initial metabarlist
fguts_Bact$pcrs$replicating_pcr <- NA
fguts_Bact$pcrs[rownames(fguts_Bact_sub$pcrs),"replicating_pcr"] <- fguts_Bact_sub$pcrs$replicating_pcr

#############################
#    Lowering tag-jumps     #
#############################

# Define a vector of thresholds to test
thresholds <- c(0,1e-4,1e-3, 1e-2, 3e-2, 5e-2) 

# Run the tests and stores the results in a list
tests <- lapply(thresholds, function(x) tagjumpslayer(fguts_Bact,x))
names(tests) <- paste("t_", thresholds, sep="")

# Format the data for ggplot with amount of reads at each threshold
tmp <- reshape2::melt(as.matrix(do.call("rbind", lapply(tests, function(x) rowSums(x$reads)))))
colnames(tmp) <- c("threshold", "sample", "abundance")

# Add richness in MOTUs at each threshold
tmp$richness <-
  reshape2::melt(as.matrix(do.call("rbind", lapply(tests, function(x) {
    rowSums(x$reads > 0)
  }))))$value

# Add control type information on pcrs and make data curation threshold numeric
tmp$controls <- fguts_Bact$pcrs$control_type[match(tmp$sample, rownames(fguts_Bact$pcrs))]
tmp$threshold <- as.numeric(gsub("t_", "", tmp$threshold))

# New table formatting for ggplot
tmp2 <- reshape2::melt(tmp, id.vars=colnames(tmp)[-grep("abundance|richness", colnames(tmp))])

ggplot(tmp2, aes(x=as.factor(threshold), y=value)) + 
  geom_boxplot(color="grey40") + 
  geom_vline(xintercept = which(levels(as.factor(tmp2$threshold)) == "0.01"), col="orange", lty=2) + 
  geom_jitter(aes(color=controls), width = 0.2, alpha=0.5) + 
  scale_color_manual(values = c("brown", "red", "cyan4","pink"), na.value = "darkgrey") +
  facet_wrap(~variable+controls, scale="free_y", ncol=5) + 
  theme_bw() + 
  scale_y_log10() +
  labs(x="MOTU pcr : total abundance filtering threshold", y="# Reads/MOTUs") + 
  theme(panel.grid = element_blank(), 
        strip.background = element_blank(), 
        axis.text.x = element_text(angle=40, h=1), 
        legend.position = "none")


#################################################
#     Summarizing the noise in the dataset      #
#################################################

### MOTUs artifact viz

# Create a table of MOTUs quality criteria 
# noise is identified as FALSE in fguts_Bact, the "!" transforms it to TRUE
motus.qual <- !fguts_Bact$motus[,c("not_an_extraction_conta","not_a_sequencing_conta", "not_degraded")]
colnames(motus.qual) <- c("extraction_conta", "sequencing_conta", "degraded_seq")

# Proportion of MOTUs potentially artifactual (TRUE) based on the criteria used
prop.table(table(apply(motus.qual, 1, sum) > 0))

# Corresponding proportion of artifactual reads (TRUE)
prop.table(xtabs(fguts_Bact$motus$count~apply(motus.qual, 1, sum) > 0))

# Proportion of MOTUs and reads potentially artifactual for each criterion
apply(motus.qual, 2, sum) / nrow(motus.qual)

apply(motus.qual, 2, function(x) sum(fguts_Bact$motus$count[x])/sum(fguts_Bact$motus$count))

tmp.motus <- 
  apply(sapply(1:ncol(motus.qual), function(x) {
    ifelse(motus.qual[,x]==T, colnames(motus.qual)[x], NA)}), 1, function(x) {
      paste(sort(unique(x)), collapse = "|")
    })
tmp.motus <- as.data.frame(gsub("^$", "not_artefactual", tmp.motus))
colnames(tmp.motus) <-  "artefact_type"

ggplot(tmp.motus, aes(x=1, fill=artefact_type)) +
  geom_bar() +  xlim(0, 2) +
  labs(fill="Artifact type") + 
  coord_polar(theta="y") + theme_void() + 
  scale_fill_brewer(palette = "Set3") + 
  theme(legend.direction = "vertical") + 
  ggtitle("MOTUs artefacts overview")


### pcrs artifacts viz

# Create a table of pcrs quality criteria 
# noise is identified as FALSE in fguts_Bact, the "!" transforms it to TRUE
pcrs.qual <- !fguts_Bact$pcrs[,c("low_contamination_level", "seqdepth_ok", "replicating_pcr")]
colnames(pcrs.qual) <- c("high_contamination_level", "low_seqdepth", "outliers")

# Proportion of pcrs potentially artifactual (TRUE) based on the criteria used excluding controls
prop.table(table(apply(pcrs.qual[fguts_Bact$pcrs$type=="sample",], 1, sum) > 0))

# Proportion of MOTUs and reads potentially artifactual for each criterion
apply(pcrs.qual[fguts_Bact$pcrs$type=="sample",], 2, sum) / nrow(pcrs.qual[fguts_Bact$pcrs$type=="sample",])

tmp.pcrs <- 
  apply(sapply(1:ncol(pcrs.qual), function(x) {
    ifelse(pcrs.qual[fguts_Bact$pcrs$type=="sample",x]==T, 
           colnames(pcrs.qual)[x], NA)}), 1, function(x) {
             paste(sort(unique(x)), collapse = "|")
           })
tmp.pcrs <- as.data.frame(gsub("^$", "not_artefactual", tmp.pcrs))

colnames(tmp.pcrs) <- "artefact_type"

ggplot(tmp.pcrs, aes(x=1, fill=artefact_type)) +
  geom_bar() +  xlim(0, 2) +
  labs(fill="Artifact type") + 
  coord_polar(theta="y") + theme_void() + 
  scale_fill_brewer(palette = "Set3") + 
  theme(legend.direction = "vertical") +
  ggtitle("PCR artefacts overview")

#################################################
#        Data cleaning and aggregation          #
#################################################

#### Removing spurious signal in MOTUs

# Use tag-jump corrected metabarlist with the threshold identified above
tmp <- tests[["t_0.01"]]

# Subset on MOTUs: we keep motus that are defined as TRUE following the 
# three criteria below (sum of four TRUE is equal to 4 with the rowSums function)
tmp <- subset_metabarlist(tmp, "motus", 
                          indices = rowSums(tmp$motus[,c("not_an_extraction_conta","not_a_sequencing_conta" ,"not_degraded")]) == 3)
summary_metabarlist(tmp)
summary_metabarlist(fguts_Bact)

#### Removing spurious signal in pcrs
# Subset on pcrs and keep only controls 
## identify the wrong replicates
WrongRep <- gsub("_.*", "", rownames(tmp$pcrs[!tmp$pcrs$replicating_pcr,]))
WrongRep <- unique(WrongRep[duplicated(gsub("_.*", "", WrongRep))])

## as the sequencing is performed in triplicate check whether we could remove one replicate to avoid removing the whole sample
Test <- tmp$reads[grep(paste(WrongRep, collapse="|"), rownames(tmp$reads)),]

### Distinguish between pcrs obtained from samples from positive controls
#### after removing contaminants, some replicate have 0 reads, remove them and plot the remaining one to check their distances\
tmp2 <- subset_metabarlist(tmp, "pcrs", 
                   indices = !rownames(tmp$pcrs) %in% names(which(rowSums(Test) == 0)))

### consider remaining replicate as Ok since the discrepancy among replicates mainly came from sample without reads
tmp2$pcrs$replicating_pcr[grep(paste(names(which(rowSums(Test) != 0)), collapse="|"), rownames(tmp2$pcrs))] <- TRUE

# Subset on pcrs and keep only controls 
fguts_Bact_clean <- subset_metabarlist(tmp2, "pcrs", 
                                       indices = tmp2$pcrs$replicating_pcr == TRUE & tmp2$pcrs$type == "sample")

summary_metabarlist(fguts_Bact_clean)

#Now check if previous subsetting leads to any empty pcrs or MOTUs
if(sum(colSums(fguts_Bact_clean$reads)==0)>0){print("empty motus present")}
if(sum(rowSums(fguts_Bact_clean$reads)==0)>0){print("empty pcrs present")}

# as there is still some pcrs with 0 reads count, we are removing them
fguts_Bact_clean <- subset_metabarlist(fguts_Bact_clean, "pcrs", 
                                       indices = !rownames(fguts_Bact_clean$pcrs) %in% names(which(rowSums(fguts_Bact_clean$reads) == 0)))

summary_metabarlist(fguts_Bact_clean)

# look at final sample size
# remove water sample from the dataset before looking at sample size
smplSize <- fguts_Bact_clean$samples[-grep("Eau", rownames(fguts_Bact_clean$samples)),]

# display sample size 
table(smplSize$treatment_contam, smplSize$treatment_inj, smplSize$pop)

# Since we have now removed certain OTUs or reduced their read counts, 
# we need to update some parameters in the metabarlist (e.g. counts, etc.)
fguts_Bact_clean$motus$count = colSums(fguts_Bact_clean$reads)
fguts_Bact_clean$pcrs$nb_reads_postmetabaR = rowSums(fguts_Bact_clean$reads)
fguts_Bact_clean$pcrs$nb_motus_postmetabaR = rowSums(ifelse(fguts_Bact_clean$reads>0, T, F))

# We can now compare some basic characteristics of the dataset before and after data curation with metabaR
check <- reshape2::melt(fguts_Bact_clean$pcrs[,c("nb_reads", "nb_reads_postmetabaR", 
                                                 "nb_motus", "nb_motus_postmetabaR")])
check$type <- ifelse(grepl("motus", check$variable), "richness", "abundance")

ggplot(data = check, aes(x = variable, y = value)) +
  geom_boxplot( color = "darkgrey") +
  geom_jitter(alpha=0.1, color = "darkgrey") +
  theme_bw() +
  facet_wrap(~type, scales = "free", ncol = 5) +
  theme(axis.text.x = element_text(angle=45, h=1))

# Let's now see if the signal is changed in terms of beta diversity

# Get row data only for samples
tmp <- subset_metabarlist(fguts_Bact, table = "pcrs",
                          indices = fguts_Bact$pcrs$type == "sample")

# Add sample biological information for checks
tmp$pcrs$Material <- tmp$samples$Material[match(tmp$pcrs$sample_id, rownames(tmp$samples))]
tmp$pcrs$Habitat <- tmp$samples$Habitat[match(tmp$pcrs$sample_id, rownames(tmp$samples))]

fguts_Bact_clean$pcrs$Material <-
  fguts_Bact_clean$samples$Material[match(fguts_Bact_clean$pcrs$sample_id,
                                          rownames(fguts_Bact_clean$samples))]
fguts_Bact_clean$pcrs$Habitat <-
  fguts_Bact_clean$samples$Habitat[match(fguts_Bact_clean$pcrs$sample_id,
                                         rownames(fguts_Bact_clean$samples))]

# Build PCoA ordinations 
mds1 <- check_pcr_repl(tmp,
                       groups = tmp$pcrs$matrix)

mds2 <- check_pcr_repl(fguts_Bact_clean,
                       groups = fguts_Bact_clean$pcrs$matrix)

# Custom colors
a <- mds1 + labs(color = "matrix") + 
  scale_color_manual(values = c("goldenrod4", "blue4")) +
  theme(legend.position = "none") + 
  ggtitle("Raw data")
b <- mds2 + labs(color = "matrix") +
  scale_color_manual(values = c("goldenrod4", "blue4")) + 
  ggtitle("Clean data")

# at the first sight we can see that water samples are very close together contrary to gut samples

# Assemble plots
leg <- get_legend(b + guides(shape=F) + 
                    theme(legend.position = "right", 
                          legend.direction = "vertical"))
cowplot::ggdraw() +
  cowplot::draw_plot(a, x=0, y=0, width = 0.4, height = 1) + 
  cowplot::draw_plot(b + guides(color=F, shape=F), x=0.42, y=0, width = 0.4, height = 1) +
  cowplot::draw_grob(leg, x=0.4, y=0)

################################################
#        Data (replicates) aggregation         #
################################################
fguts_Bact_agg <- aggregate_pcrs(fguts_Bact_clean, FUN = FUN_agg_pcrs_mean)
summary_metabarlist(fguts_Bact_agg)

#######################################
#         some visualizations         #
#######################################

# Aggregate the metabarlist at the phylum level
fish_bacteria_phy <-
  aggregate_motus(fguts_Bact_agg, groups = fguts_Bact_agg$motus$phylum)

# Aggregate per Habitat/material
tmp <-
  reshape2::melt(aggregate(fish_bacteria_phy$reads, by = list(
    paste(fish_bacteria_phy$samples$experiment,
          fish_bacteria_phy$samples$matrix,
          sep = " | ")), sum))

# Plot 
nb.cols <- 18
mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(nb.cols)
ggplot(tmp, aes(x=Group.1, y=value, fill=variable)) +
  geom_bar(stat="identity") + 
  labs(x=NULL, y="#reads", fill="Phyla") + 
  coord_flip() + theme_bw() + 
  scale_fill_manual(values = mycolors) + 
  theme(legend.position = "bottom") +
  ggtitle("Bacteria phyla")


#######################################
#        Rarefaction curves           #
#######################################

# lets take a final look at the diversity coverage of samples with rarefaction curves 
#(rather than our previous look at the pcrs level)

fguts_Bact_agg.raref = hill_rarefaction(fguts_Bact_agg, nboot = 20, nsteps = 10)
material <- paste(fguts_Bact_agg$samples$experiment, fguts_Bact_agg$samples$matrix)
material <- setNames(material,rownames(fguts_Bact_agg$samples))
options(max.print=99999)
fguts_Bact_agg.raref$hill_table
which(is.na(fguts_Bact_agg.raref$hill_table))

#plot
p <- gghill_rarefaction(fguts_Bact_agg.raref, group=material)
p + scale_fill_manual(values = c("goldenrod4", "blue4")) +
  scale_color_manual(values = c("goldenrod4", "blue4")) +
  labs(color="matrix type")

saveRDS(fguts_Bact_agg, file.path(savingDir, "fguts_Bact_agg_MergedRep.RDS"))
