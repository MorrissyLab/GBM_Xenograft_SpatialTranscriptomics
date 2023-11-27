## CellChat 
library(dplyr)
library(tidyr)
library(stringr)
library(Seurat)
library(CellChat)
library(nichenetr)
library(UpSetR)
library(tibble)
library(patchwork)
library(openxlsx)
library(reshape2)
library(ggplot2)
library(wesanderson)
library(ggsci)
library(ggtree)
library(ggraph)
library(igraph)
library(tidyverse)
library(RColorBrewer)
library(circlize)
library(chorddiag)
library(vegan)
library(gridExtra)
library(cowplot)
options(stringsAsFactors = FALSE)


###### 1. CellChat - D1 : D4 
## Prepare data
metadata_all <- read.csv("Spots_metadata.csv", row.names = 1, stringsAsFactor=FALSE)
metadata_all <- metadata_all %>% select(Sample, Line, Mouse, Replicate, Day, Site, Patient, type, Time_point, Hs_admix)
metadata_sub <- metadata_all %>% filter(Density_group != "")

data.orig <- Read10X_h5("/work/morrissy_lab/vthoppey/SpatialData/Analysis/SpatialAnalysis/morrissy-chan/pdx_merge_all/filtered_feature_bc_matrix.h5")
rownames(data.orig) <- gsub("GRCh38_","GRCH38.",rownames(data.orig)); rownames(data.orig) <- gsub("mm10___","MM10.",rownames(data.orig))
data.orig_sub <- data.orig[,rownames(metadata_sub)]

visium.brain <- CreateSeuratObject(counts = data.orig_sub, project = "visium", meta.data = metadata_sub) 
visium.brain <- NormalizeData(visium.brain, normalization.method = "LogNormalize", scale.factor = 10000)

## Run CellChat 
data.input = visium.brain@assays$RNA@data # normalized data matrix
meta = data.frame(labels = visium.brain@meta.data$Density_group, row.names = rownames(visium.brain@meta.data))
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents))

CellChatDB <- readRDS("CellChatDB/CellChatDB_Hs_Mm_updated.rds")
CellChatDB[[1]] <- CellChatDB[[1]][-grep("H2-BI", CellChatDB[[1]]$ligand),]
CellChatDB[[1]] <- CellChatDB[[1]][-grep("H2-Ea-ps", CellChatDB[[1]]$ligand),]
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

df.net <- subsetCommunication(cellchat)
df.net <- df.net %>% separate(interaction_name_2, c("ligand_2","receptor_2"),sep = " - ")
df.net$ligand_species <- "" ; df.net$receptor_species <- ""
df.net$ligand_2 <- gsub(" ","",df.net$ligand_2); df.net$receptor_2 <- gsub(" ","",df.net$receptor_2)
df.net$ligand_species[grep("GRCH38.", df.net$ligand_2)] <- "Human"; df.net$receptor_species[grep("GRCH38.", df.net$receptor_2)] <- "Human"
df.net$ligand_species[grep("MM10.", df.net$ligand_2)] <- "Mouse"; df.net$receptor_species[grep("MM10.", df.net$receptor_2)] <- "Mouse"
df.net$interaction_name_2 <- paste(df.net$ligand_2, df.net$receptor_2, sep=" - ")
df.net$category <- paste(df.net$ligand_species, df.net$receptor_species, sep="->")
df.net$ST <- paste(df.net$source, df.net$target, sep="->")
categories <- unique(df.net$category)
cellchat@LR$LRsig$interaction_name_2 <- gsub("  - "," - ",cellchat@LR$LRsig$interaction_name_2)
for(i in 1:nrow(df.net)){
    df.net$interaction_name[i] <- cellchat@LR$LRsig$interaction_name[which(cellchat@LR$LRsig$interaction_name_2 == df.net$interaction_name_2[i])]
}
write.table(df.net, "cellchatNET_byTumorDensity.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


###### 1. CellChat - WM, CP, Vasculature (m30, m52, m59, m86)
## Prepare data
metadata_all <- read.csv("Spots_metadata.csv", row.names = 1, stringsAsFactor=FALSE)
metadata_all <- metadata_all %>% select(Sample, Line, Mouse, Replicate, Day, Site, Patient, type, Time_point, Hs_admix)
metadata_sub <- metadata_all %>% filter(D1_annot != "")

data.orig <- Read10X_h5("/work/morrissy_lab/vthoppey/SpatialData/Analysis/SpatialAnalysis/morrissy-chan/pdx_merge_all/filtered_feature_bc_matrix.h5")
rownames(data.orig) <- gsub("GRCh38_","GRCH38.",rownames(data.orig)); rownames(data.orig) <- gsub("mm10___","MM10.",rownames(data.orig))
data.orig_sub <- data.orig[,rownames(metadata_sub)]

visium.brain <- CreateSeuratObject(counts = data.orig_sub, project = "visium", meta.data = metadata_sub) 
visium.brain <- NormalizeData(visium.brain, normalization.method = "LogNormalize", scale.factor = 10000)

## Run CellChat 
data.input = visium.brain@assays$RNA@data # normalized data matrix
meta = data.frame(labels = visium.brain@meta.data$D1_annot, row.names = rownames(visium.brain@meta.data))
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents))

CellChatDB <- readRDS("/work/morrissy_lab/vthoppey/SpatialData/Analysis/CellChat2/CellChatDB/CellChatDB_Hs_Mm_updated.rds")
CellChatDB[[1]] <- CellChatDB[[1]][-grep("H2-BI", CellChatDB[[1]]$ligand),]
CellChatDB[[1]] <- CellChatDB[[1]][-grep("H2-Ea-ps", CellChatDB[[1]]$ligand),]
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.05)
#cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = 1)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

df.net <- subsetCommunication(cellchat)
df.net <- df.net %>% separate(interaction_name_2, c("ligand_2","receptor_2"),sep = " - ")
df.net$ligand_species <- "" ; df.net$receptor_species <- ""
df.net$ligand_2 <- gsub(" ","",df.net$ligand_2); df.net$receptor_2 <- gsub(" ","",df.net$receptor_2)
df.net$ligand_species[grep("GRCH38.", df.net$ligand_2)] <- "Human"; df.net$receptor_species[grep("GRCH38.", df.net$receptor_2)] <- "Human"
df.net$ligand_species[grep("MM10.", df.net$ligand_2)] <- "Mouse"; df.net$receptor_species[grep("MM10.", df.net$receptor_2)] <- "Mouse"
df.net$interaction_name_2 <- paste(df.net$ligand_2, df.net$receptor_2, sep=" - ")
df.net$category <- paste(df.net$ligand_species, df.net$receptor_species, sep="->")
df.net$ST <- paste(df.net$source, df.net$target, sep="->")
categories <- unique(df.net$category)
cellchat@LR$LRsig$interaction_name_2 <- gsub("  - "," - ",cellchat@LR$LRsig$interaction_name_2)
for(i in 1:nrow(df.net)){
    df.net$interaction_name[i] <- cellchat@LR$LRsig$interaction_name[which(cellchat@LR$LRsig$interaction_name_2 == df.net$interaction_name_2[i])]
}
write.table(df.net, "cellchatNET_Invasion.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

