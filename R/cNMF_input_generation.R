library(dplyr)
library(tidyr)
library(stringr)
library(Seurat)

# 1. Data normalizeation and calculation of admixture
GBM.data <- Read10X_h5("filtered_feature_bc_matrix.h5")

# Initialize the Seurat object with the raw (non-normalized data).
GBM_xeno <- CreateSeuratObject(counts = GBM.data, project = "GBM_xenografts", min.cells = 3)

# Calculate percent mito for GRCh38 + mm10
mito_genes <- c( rownames(GBM_xeno)[grep("^GRCh38-MT-", rownames(GBM_xeno))], rownames(GBM_xeno)[grep("^mm10---mt-", rownames(GBM_xeno))] )
GBM_xeno[["percent.mt"]] <- PercentageFeatureSet(GBM_xeno, features = mito_genes)

# Filter spots
GBM_xeno <- subset(GBM_xeno, subset = nFeature_RNA > 200 & nCount_RNA > 500 & percent.mt < 50)

# Normalize
GBM_xeno <- NormalizeData(GBM_xeno, normalization.method = "LogNormalize", scale.factor = 10000)

# Admixture calculation
GRCh38_index <- grep("GRCh38-",rownames(GBM_xeno@assays$RNA@counts))
mm10_index <- grep("mm10---",rownames(GBM_xeno@assays$RNA@counts))
for(i in 1:nrow(GBM_xeno@meta.data)){
    GBM_xeno@meta.data$GRCh38_nCount_RNA[i] <- sum(GBM_xeno@assays$RNA@counts[,i][GRCh38_index])
    GBM_xeno@meta.data$mm10_nCount_RNA[i] <- sum(GBM_xeno@assays$RNA@counts[,i][mm10_index])
    GBM_xeno@meta.data$admixture[i] <- sum(GBM_xeno@assays$RNA@counts[,i][GRCh38_index]) / (sum(GBM_xeno@assays$RNA@counts[,i][GRCh38_index]) + sum(GBM_xeno@assays$RNA@counts[,i][mm10_index]))
}

#### Selection of OD genes
## GRCh38 
GRCh38_raw <- GBM_xeno@assays$RNA@counts[GRCh38_index,]
GRCh38_corpus <- restrictCorpus(GRCh38_raw,
                             removeAbove = 1.0,
                             removeBelow = 0.001,
                             alpha = 0.05,
                             plot = TRUE,
                             verbose = TRUE,
                             nTopOD = NA)
write.table(rownames(GRCh38_corpus), "GRCh38_ODgenes.txt",sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# Filter out cells with zero expression OD genes
spots_to_remove <- colnames(GRCh38_corpus)[which(colSums(GRCh38_corpus) == 0)]
GRCh38_raw_filt <- GRCh38_raw[-which(colnames(GRCh38_raw_filt) %in% spots_to_remove)]
GRCh38_norm_filt <- GBM_xeno@assays$RNA@data[GRCh38_index,][,colnames(GRCh38_raw_filt)]
write.table(t(GRCh38_raw_filt), "GRCh38_GBM_Xeno_raw_filt.txt", sep="\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(t(GRCh38_norm_filt), "GRCh38_GBM_Xeno_norm_filt.txt", sep="\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
## mm10
mm10_raw <- GBM_xeno@assays$RNA@counts[mm10_index,]
mm10_corpus <- restrictCorpus(mm10_raw,
                             removeAbove = 1.0,
                             removeBelow = 0.001,
                             alpha = 0.05,
                             plot = TRUE,
                             verbose = TRUE,
                             nTopOD = NA)
write.table(rownames(mm10_corpus), "mm10_ODgenes.txt",sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
# Filter out cells with zero expression OD genes
spots_to_remove <- colnames(mm10_corpus)[which(colSums(mm10_corpus) == 0)]
mm10_raw_filt <- mm10_raw[-which(colnames(mm10_raw_filt) %in% spots_to_remove)]
mm10_norm_filt <- GBM_xeno@assays$RNA@data[mm10_index,][,colnames(mm10_raw_filt)]
write.table(t(mm10_raw_filt), "mm10_GBM_Xeno_raw_filt.txt", sep="\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(t(mm10_norm_filt), "mm10_GBM_Xeno_norm_filt.txt", sep="\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

