library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)
library(tibble)
library(ggplot2)
library(survival)
library(survminer)
library(DESeq2)
library(fgsea)
library(RColorBrewer)
set.seed(123456)


## Load expression data
# TCGA - NEW
# expMat <- read.table("TCGA/TcgaTargetGtex_RSEM_Hugo_norm_count", sep="\t", header = TRUE, stringsAsFactor=FALSE) ; colnames(expMat) <- gsub("[.]","-",colnames(expMat))
# category <- read.table("TCGA/TCGA_GTEX_category.txt", sep="\t", header = TRUE, stringsAsFactor=FALSE) ; category_GBM <- category %>% filter(TCGA_GTEX_main_category == "TCGA Glioblastoma Multiforme") ; rownames(category_GBM) <- category_GBM$sample
# survival <- read.table("TCGA/TCGA_survival_data", sep="\t", header = TRUE, stringsAsFactor=FALSE) ; survival_GBM <- survival %>% filter(sample %in% category_GBM$sample) ; rownames(survival_GBM) <- survival_GBM$sample ; survival_GBM <- survival_GBM[category_GBM$sample,]
# phenotype <- read.table("TCGA/TcgaTargetGTEX_phenotype.txt", sep="\t", header = TRUE, stringsAsFactor=FALSE) ; phenotype_GBM <- phenotype %>% filter(sample %in% category_GBM$sample) ; rownames(phenotype_GBM) <- phenotype_GBM$sample; phenotype_GBM <- phenotype_GBM[category_GBM$sample,]
# rownames(expMat) <- expMat[,1]
# expMat_GBM <- expMat[,category_GBM$sample]
# metadata <- cbind(category_GBM, survival_GBM, phenotype_GBM)
# metadata <- metadata[colnames(expMat_GBM),]
# metadata <- metadata[,-c(3,12)] 
# metadata$patient <- metadata$sample
# metadata <- metadata %>% separate(patient, c("X1","X2","X3","X4"), sep="-") 
# metadata$patient <- paste(metadata$X1, metadata$X2, metadata$X3, sep="-") 
# saveRDS(expMat_GBM, "TCGA/TCGA_GBM.rds")
#  write.table(metadata, "TCGA/TCGA_metadata.txt", sep="\t", quote= FALSE, row.names = FALSE, col.names = TRUE)

# Load gene-sets
pathways <- readRDS("c400_TS_and_hub_genelist.rds")

# Load TCGA data
expMat <- readRDS("TCGA_GBM.rds")
expMat <- as.matrix(expMat)
dataset <- "TCGA"
CAT <- data.frame(V1 = c("allSamples","primary","recurrent","recurrentIDHwt","primaryIDHwtClassical","primaryIDHwtMesenchymal","primaryIDHwtNeural","primaryIDHwtProneural"))

# Create output directories
dir.create(dataset) ; dir.create(file.path(dataset,"page_rank"))
output.dir <- file.path(dataset, "page_rank", "c400")
dir.create(output.dir)
prefix <- paste(dataset,"__", "page_rank", "__", "c400", sep="")

        
# Run fgsea
FGSEA_RES <- data.frame()
for(p in 1:length(pathways)){
    pathways[[p]] <- intersect(pathways[[p]], rownames(expMat))
}

if(!file.exists(paste(output.dir, "/", prefix, "_fgseaRes.txt",sep=""))){
    for(s in 1:ncol(expMat)){
        print(s)
        fgseaRes <- fgsea(pathways = pathways, 
                stats    = sort(expMat[,s], decreasing = TRUE),
                minSize  = 1,
                maxSize  = 3000) #,
                #nperm    = 1000)
        fgseaRes <- data.frame(fgseaRes)
        for(f in 1:nrow(fgseaRes)){
            fgseaRes$leadingEdge[f] <- paste(unlist(fgseaRes$leadingEdge[f]), collapse = "|")
        }
        fgseaRes$leadingEdge <- unlist(fgseaRes$leadingEdge)
        fgseaRes <- fgseaRes %>% add_column("sample" = colnames(expMat)[s], .before = "pathway")
        FGSEA_RES <- rbind(FGSEA_RES, fgseaRes)
    }
    write.table(FGSEA_RES, paste(output.dir, "/", prefix, "_fgseaRes.txt",sep=""),sep="\t", quote =FALSE, row.names= FALSE, col.names = TRUE)
}


    



# Create empty data.frame for results 
KM_Surv_table <- data.frame() 

### Read Data - TCGA
expMat <- expMat <- readRDS("TCGA_GBM.rds")
expMat <- as.matrix(expMat)
dataset <- "TCGA"
metadata <- read.table("TCGA_metadata.txt", sep="\t", header = TRUE, stringsAsFactor=FALSE)
     

### Run survival analysis for each category
# In each category
    # For 25, 33, 50 pct thresholds
        # For each GEP
        
for(C in 1:nrow(CAT)){

    if(C == 1){
        CI <- metadata %>% select(sample, patient_id, OS, OS.time) 
    }
    if(C == 2){
        CI <- metadata %>% filter(X_sample_type == "Primary Tumor") %>% select(sample, patient_id, OS, OS.time) 
    }
    if(C == 3){
        CI <- metadata %>% filter(X_sample_type == "Recurrent Tumor") %>% select(sample, patient_id, OS, OS.time) 
    }
    if(C == 4){
        CI <- metadata %>% filter(X_sample_type == "Primary Tumor") %>% filter(IDH_status == "WT") %>% select(sample, patient_id, OS, OS.time) 
    }
    if(C == 5){
        CI <- metadata %>% filter(X_sample_type == "Recurrent Tumor") %>% filter(IDH_status == "WT") %>% select(sample, patient_id, OS, OS.time)
    }
    if(C == 6){
        CI <- metadata %>% filter(X_sample_type == "Primary Tumor") %>% filter(IDH_status == "WT") %>% filter(Subtype == "Classical") %>% select(sample, patient_id, OS, OS.time)
    }
    if(C == 7){
        CI <- metadata %>% filter(X_sample_type == "Primary Tumor") %>% filter(IDH_status == "WT") %>% filter(Subtype == "Mesenchymal") %>% select(sample, patient_id, OS, OS.time)
    }
    if(C == 8){
        CI <- metadata %>% filter(X_sample_type == "Primary Tumor") %>% filter(IDH_status == "WT") %>% filter(Subtype == "Neural") %>% select(sample, patient_id, OS, OS.time)
    }
    if(C == 9){
        CI <- metadata %>% filter(X_sample_type == "Primary Tumor") %>% filter(IDH_status == "WT") %>% filter(Subtype == "Proneural") %>% select(sample, patient_id, OS, OS.time)
    }


    colnames(CI) <- c("sample_id","patient_id","status","time")
    rownames(CI) <- CI$sample_id


    visium_gep <- names(pathways)
    FGSEA_RES <- read.table(paste(output.dir, "/", prefix, "_fgseaRes.txt",sep=""), sep="\t", header = TRUE, stringsAsFactors=FALSE)
    samples_int <- intersect(rownames(CI), FGSEA_RES$sample)
    FGSEA_RES <- FGSEA_RES %>% filter(sample %in% samples_int)
    CI_sub <- CI[samples_int,]
    
    KM_Surv_table1 <- data.frame() 
    
    PCT <- c("25pctsamples", "33pctsamples", "50pctsamples")
    
    for(P in 1:length(PCT)){
        pdf(paste(output.dir, "/", prefix, "_OS_SurvivalCurves_", CAT$V1[C], "_", PCT[P], ".pdf",sep=""))
        
        for(g in 1:length(visium_gep)){
            tab <- FGSEA_RES %>% filter(pathway == visium_gep[g])
            #tab <- tab %>% filter(padj < 0.05)
            tab <- tab[order(tab$NES, decreasing = TRUE),]
            rownames(tab) <- NULL 
            tab$category <- ""

            if(P == 1){
                tab$category[1:round(nrow(tab)/4)] <- "top25pct"
                tab$category[((nrow(tab) - round(nrow(tab)/4) + 1)):nrow(tab)] <- "bottom25pct"
            }
            if(P == 2){
                tab$category[1:round(nrow(tab)/3)] <- "top33"
                tab$category[seq(1, nrow(tab), nrow(tab)/3)[3]:nrow(tab)] <- "bottom33"
            }
            if(P == 3){
                tab$category[1:round(nrow(tab)/2)] <- "top50pct"
                tab$category[((nrow(tab) - round(nrow(tab)/2) + 1)):nrow(tab)] <- "bottom50pct"
            }

            tab <- tab %>% filter(category != "")
            CI_gep <- CI[tab$sample,]
            tab$status <- CI_gep$status
            tab$time <- CI_gep$time
            fit <- survfit(Surv(time, status) ~ category, data = tab)

            fit_P1 <- surv_pvalue(fit)
            KM_Surv_table_sub1 <- data.frame(Dataset = dataset, SURV = "OS", CAT = CAT$V1[C], PCT = PCT[P], GEP = visium_gep[g], p.value = fit_P1$pval)

            p1 <- ggsurvplot(fit, data = tab, 
                #palette = c("#E7B800", "#2E9FDF"),
                palette = c("blue", "red"),
                conf.int = FALSE,  # Add confidence interval
                risk.table = TRUE,
                risk.table.col = "strata",
                pval = TRUE) + 
                ggtitle(visium_gep[g])
            print(p1)
            KM_Surv_table1 <- rbind(KM_Surv_table1, KM_Surv_table_sub1)
        }
        dev.off()         
    }
    KM_Surv_table <- rbind(KM_Surv_table, KM_Surv_table1)
}
write.table(KM_Surv_table, paste(output.dir, "/", prefix, "_OS_SurvivalCurves.txt",sep=""))




### T.test
Ttest_violin_df <- data.frame() 

for(C in 1:nrow(CAT)){
    
    expMat <- expMat <- readRDS("TCGA_GBM.rds")
    expMat <- as.matrix(expMat)
    dataset <- "TCGA"
    metadata <- read.table("TCGA_metadata.txt", sep="\t", header = TRUE, stringsAsFactor=FALSE)
    metadata <- metadata[order(metadata$OS.time, decreasing = TRUE),]
    

    if(C == 1){
        metadata <- metadata 
    }
    if(C == 2){
        metadata <- metadata %>% filter(X_sample_type == "Primary Tumor") 
    }
    if(C == 3){
        metadata <- metadata %>% filter(X_sample_type == "Recurrent Tumor") 
    }
    if(C == 4){
        metadata <- metadata %>% filter(X_sample_type == "Primary Tumor") %>% filter(IDH_status == "WT")
    }
    if(C == 5){
        metadata <- metadata %>% filter(X_sample_type == "Recurrent Tumor") %>% filter(IDH_status == "WT")
    }
    if(C == 6){
        metadata <- metadata %>% filter(X_sample_type == "Primary Tumor") %>% filter(IDH_status == "WT") %>% filter(Subtype == "Classical") %>% select(sample, patient_id, OS, OS.time)
    }
    if(C == 7){
        metadata <- metadata %>% filter(X_sample_type == "Primary Tumor") %>% filter(IDH_status == "WT") %>% filter(Subtype == "Mesenchymal") %>% select(sample, patient_id, OS, OS.time)
    }
    if(C == 8){
        metadata <- metadata %>% filter(X_sample_type == "Primary Tumor") %>% filter(IDH_status == "WT") %>% filter(Subtype == "Neural") %>% select(sample, patient_id, OS, OS.time)
    }
    if(C == 9){
        metadata <- metadata %>% filter(X_sample_type == "Primary Tumor") %>% filter(IDH_status == "WT") %>% filter(Subtype == "Proneural") %>% select(sample, patient_id, OS, OS.time)
    }


    metadata$category50 <- ""
    metadata$category50[1:round(nrow(metadata)/2)] <- "top50"
    metadata$category50[(round(nrow(metadata)/2) + 1):nrow(metadata)] <- "bottom50"
    metadata$category25 <- ""
    metadata$category25[1:round(nrow(metadata)/4)] <- "top25"
    metadata$category25[((nrow(metadata) - round(nrow(metadata)/4) + 1)):nrow(metadata)] <- "bottom25"
    metadata$category33 <- ""
    metadata$category33[1:round(nrow(metadata)/3)] <- "top33"
    metadata$category33[seq(1, nrow(metadata), nrow(metadata)/3)[3]:nrow(metadata)] <- "bottom33"


    visium_gep <- names(pathways)
    FGSEA_RES <- read.table(paste(output.dir, "/", prefix, "_fgseaRes.txt",sep=""), sep="\t", header = TRUE, stringsAsFactors=FALSE)
    samples_int <- intersect(metadata$sample, FGSEA_RES$sample)
    FGSEA_RES <- FGSEA_RES %>% filter(sample %in% samples_int)

    Ttest_res <- data.frame()

    high50 <- metadata$sample[which(metadata$category50 == "top50")]
    low50 <- metadata$sample[which(metadata$category50 == "bottom50")]
    high25 <- metadata$sample[which(metadata$category25 == "top25")]
    low25 <- metadata$sample[which(metadata$category25 == "bottom25")]
    high33 <- metadata$sample[which(metadata$category33 == "top33")]
    low33 <- metadata$sample[which(metadata$category33 == "bottom33")]
    
    

    for(g in 1:length(visium_gep)){
        
        res_sub <- FGSEA_RES %>% filter(pathway == visium_gep[g])
        
        NES_high50 <- res_sub$NES[which(res_sub$sample %in% high50)]
        NES_low50 <- res_sub$NES[which(res_sub$sample %in% low50)]
        plot50 <- data.frame(NES = c(NES_high50, NES_low50), category1 = c(rep("NES_high50", length(NES_high50)), rep("NES_low50", length(NES_low50))))
        plot50$Patient_percent <- "pct50"
        test50 <- t.test(NES_high50, NES_low50)
        test_df50 <- data.frame(category = "high_low_50pct",
                            pathway = visium_gep[g], 
                            high_mean = as.numeric(test50$estimate[1]),
                            low_mean = as.numeric(test50$estimate[2]),
                            statistic = as.numeric(test50$statistic),
                            p.value = as.numeric(test50$p.value))
        
        NES_high25 <- res_sub$NES[which(res_sub$sample %in% high25)]
        NES_low25 <- res_sub$NES[which(res_sub$sample %in% low25)]
        plot25 <- data.frame(NES = c(NES_high25, NES_low25), category1 = c(rep("NES_high25", length(NES_high25)), rep("NES_low25", length(NES_low25))))
        plot25$Patient_percent <- "pct25"
        test25 <- t.test(NES_high25, NES_low25)
        test_df25 <- data.frame(category = "high_low_25pct",
                            pathway = visium_gep[g], 
                            high_mean = as.numeric(test25$estimate[1]),
                            low_mean = as.numeric(test25$estimate[2]),
                            statistic = as.numeric(test25$statistic),
                            p.value = as.numeric(test25$p.value))

        NES_high33 <- res_sub$NES[which(res_sub$sample %in% high33)]
        NES_low33 <- res_sub$NES[which(res_sub$sample %in% low33)]
        plot33 <- data.frame(NES = c(NES_high33, NES_low33), category1 = c(rep("NES_high33", length(NES_high33)), rep("NES_low33", length(NES_low33))))
        plot33$Patient_percent <- "pct33"
        test33 <- t.test(NES_high33, NES_low33)
        test_df33 <- data.frame(category = "high_low_33pct",
                            pathway = visium_gep[g], 
                            high_mean = as.numeric(test33$estimate[1]),
                            low_mean = as.numeric(test33$estimate[2]),
                            statistic = as.numeric(test33$statistic),
                            p.value = as.numeric(test33$p.value))

        test_df <- rbind(test_df50, test_df25, test_df33)
        Ttest_res <- rbind(Ttest_res, test_df)

        plot_df <- rbind(plot25, plot33, plot50)   
        plot_df$Survival_categ <- "high"
        plot_df$Survival_categ[grep("low", plot_df$category1)] <- "low"
        plot_df$Survival <- "OS"
        plot_df$Disease_category <- CAT$V1[C]
        plot_df$GEP <- visium_gep[g]

        Ttest_violin_df <- rbind(Ttest_violin_df, plot_df)
    }
    # dev.off()
    
    Ttest_res <- Ttest_res[order(Ttest_res$category),]
    write.table(Ttest_res, paste(output.dir, "/", prefix, "_", "OS", "_", CAT$V1[C], "__tTest_results.txt",sep=""), sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

    Ttest_OS <- Ttest_violin_df %>% filter(Survival == "OS")
    pdf(paste(output.dir, "/", prefix, "_", "OS", "__tTest_results_violin.pdf",sep=""), width = 12, height = 10)
    for(g in 1:length(visium_gep)){
        sub1 <- Ttest_OS %>% filter(GEP == visium_gep[g])
        T1 <- ggplot(sub1, aes(x = Survival_categ, y=NES, fill = Survival_categ)) + 
                    geom_violin(trim = TRUE) + geom_boxplot(width=0.1) + 
                    theme_classic(base_size = 13) +
                    stat_compare_means(comparisons = list( c("low", "high") ), method = "t.test", label = "p.signif") +
                    scale_fill_manual(values=brewer.pal(n = 3, name = "RdBu")[c(1,3)]) +
                    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
                    theme(legend.position = "bottom") + 
                    facet_grid(Patient_percent ~ Disease_category) + 
                    labs(title = visium_gep[g])
        print(T1)
    }
    dev.off()
}







# library(tibble) 
# files <- list.files("TCGA/")[grep("_tTest_results", list.files("TCGA/"))][1:4]
# all <- data.frame()
# for(f in 1:length(files)){
#     f1 <- read.table(file.path("TCGA", files[f]), sep="\t", header = TRUE, stringsAsFactor=FALSE)
#     name <- gsub("SurvivalAnalysis_","",files[f])
#     name <- gsub("_tTest_results.txt","", name)
#     f1 <- f1 %>% add_column(sample_type = name, .before = "category")
#     all <- rbind(all, f1)
# }
# write.table(all, "TCGA/TCGA_PFI_SurvivalAnalyis_tTest_results.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


