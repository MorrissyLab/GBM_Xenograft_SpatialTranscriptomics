library(AUCell)

### Load data
expr_file <- "Hs_all.gene_spectra_score.k_15.dt_0_10.txt"  ### Mm_all.gene_spectra_score.k_15.dt_0_10.txt for mm10

exprMat <- read.table(expr_file, sep="\t", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
rownames(exprMat) <- paste("GEP_",seq(1,nrow(exprMat)),sep="")
colnames(exprMat) <- gsub("GRCh38_","",colnames(exprMat))
colnames(exprMat) <- gsub("mm10___","",colnames(exprMat))
exprMat <- t(exprMat)
exprMat <- as.matrix(exprMat)
print("Loaded exprMat")

### Initialize settings
library(SCENIC)
org <- "hgnc" # or hgnc, or dmel
dbDir <- "cisTarget_databases" # RcisTarget databases location
myDatasetTitle <- "SCENIC on visium pdx - mouse" # choose a name for your analysis
dbs <- c("hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather","hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather") # mm10 - c("mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather","mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather")
names(dbs) <- c("500bp","10kb")
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=20) 
saveRDS(scenicOptions, file="int/scenicOptions.Rds")
print("Initialized scenicOptions")


### Co-expression network
print("Calculating correlation")
#runCorrelation(exprMat, scenicOptions)
print("runCorrelation complete")

print("Starting runGenie3")
#runGenie3(exprMat, scenicOptions)
print("runGenie3 complete")

### Build and score the GRN
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123


print("runSCENIC_1")
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)

print("runSCENIC_2")
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions)

print("runSCENIC_3")
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat)

saveRDS(scenicOptions, file="int/scenicOptions.Rds")
print("saved scenicOptions")

