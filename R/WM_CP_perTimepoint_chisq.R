library(dplyr)
library(tidyr)
library(stringr)
library(ComplexHeatmap)
library(ggplot2)
library(tibble)
library(reshape2)
library(openxlsx)
library(ggpubr)
library(gridExtra)
library(cowplot)
library(ggprism)

MD <- read.csv('Spots_metadata.csv', stringsAsFactors=FALSE) 
MD$Sample2 <- paste(MD$Patient,MD$Line, MD$Mouse, "_", MD$Day, "d_", MD$Time_point, sep="") #_", MD$Site,sep="")
MD$Time_point[which(MD$Sample2 == "BT143_x2_48d")] <- "mid"

MD <- MD %>% select(Barcode, Sample, Sample2, Patient, Line, Mouse, Day, Site, Hs_admix, Time_point) 
MD$Density_group <- ""
MD$Density_group[which(MD$Hs_admix > 0 & MD$Hs_admix <= 0.05)] <- "D0"
MD$Density_group[which(MD$Hs_admix > 0.05 & MD$Hs_admix <= 0.2)] <- "D1"
MD$Density_group[which(MD$Hs_admix > 0.2 & MD$Hs_admix <= 0.5)] <- "D2"
MD$Density_group[which(MD$Hs_admix > 0.5 & MD$Hs_admix <= 0.8)] <- "D3"
MD$Density_group[which(MD$Hs_admix > 0.8 & MD$Hs_admix <= 1)] <- "D4"

# MD <- MD #%>% filter(Density_group %in% c("D1","D2"))
MD$Sample_RepMerged <- paste(MD$Patient,"_", MD$Line, MD$Mouse,"_", MD$Day,"d",sep="") 
MD$Sample_RepMerged_SideSep <- paste(MD$Patient,"_", MD$Line, MD$Mouse,"_", MD$Day,"d","_", MD$Site, sep="") 
 

mCon <- read.table(cNMF_dir$consensus_file[1], sep="\t", header = TRUE, stringsAsFactor=FALSE)
colnames(mCon) <- gsub("Usage_","mGEP",colnames(mCon))


### Per time-point
TP <- c("early","mid","late")
cont_tab <- vector("list",3)
names(cont_tab) <- TP
TAB_all <- data.frame()

for(s in 1:length(TP)){
    df1 <- data.frame(matrix(0, nrow = 6, ncol = 4, dimnames = list(c("Brain","Tumor","Tumor_D1","Tumor_D2","Tumor_D3","Tumor_D4"), c("nr_WM","nr_nonWM", "pct_WM","pct_nonWM"))))
    sample_sub <- MDsub %>% filter(Time_point == TP[s])
    
    mCon_sub <- mCon[sample_sub$Barcode,]
    mGEP19_spots <- rownames(mCon_sub)[which(mCon_sub[,"mGEP19"] > 0.1)]   ##mGEP1 > 0.1 for CP 
    mGEP71_spots <- rownames(mCon_sub)[which(mCon_sub[,"mGEP71"] > 0.1)]   ##mGEP73 > 0.1 for CP
    
    #Brain
    WM_spots <- unique(c(mGEP19_spots,mGEP71_spots))
    nonWM_spots <- setdiff(sample_sub$Barcode, WM_spots)
    df1$nr_WM[1] <- length(WM_spots)
    df1$nr_nonWM[1] <- length(nonWM_spots)
    df1$pct_WM[1] <- length(WM_spots)/length(sample_sub$Barcode)
    df1$pct_nonWM[1] <- length(nonWM_spots)/length(sample_sub$Barcode)

    #Tumor
    Tumor_spots <- sample_sub$Barcode[which(sample_sub$Hs_admix > 0.05)]
    Tumor_spots_WM <- intersect(Tumor_spots, WM_spots)
    Tumor_spots_nonWM <- setdiff(Tumor_spots, Tumor_spots_WM)
    df1$nr_WM[2] <- length(Tumor_spots_WM)
    df1$nr_nonWM[2] <- length(Tumor_spots_nonWM)
    df1$pct_WM[2] <- length(Tumor_spots_WM)/length(Tumor_spots)
    df1$pct_nonWM[2] <- length(Tumor_spots_nonWM)/length(Tumor_spots)

    #Tumor_D1
    TumorD1_spots <- sample_sub$Barcode[which(sample_sub$Hs_admix > 0.05 & sample_sub$Density_group == "D1")]
    TumorD1_spots_WM <- intersect(TumorD1_spots, WM_spots)
    TumorD1_spots_nonWM <- setdiff(TumorD1_spots, TumorD1_spots_WM)
    df1$nr_WM[3] <- length(TumorD1_spots_WM)
    df1$nr_nonWM[3] <- length(TumorD1_spots_nonWM)
    df1$pct_WM[3] <- length(TumorD1_spots_WM)/length(TumorD1_spots)
    df1$pct_nonWM[3] <- length(TumorD1_spots_nonWM)/length(TumorD1_spots)

    #Tumor_D2
    TumorD2_spots <- sample_sub$Barcode[which(sample_sub$Hs_admix > 0.05 & sample_sub$Density_group == "D2")]
    TumorD2_spots_WM <- intersect(TumorD2_spots, WM_spots)
    TumorD2_spots_nonWM <- setdiff(TumorD2_spots, TumorD2_spots_WM)    
    df1$nr_WM[4] <- length(TumorD2_spots_WM)
    df1$nr_nonWM[4] <- length(TumorD2_spots_nonWM)
    df1$pct_WM[4] <- length(TumorD2_spots_WM)/length(TumorD2_spots)
    df1$pct_nonWM[4] <- length(TumorD2_spots_nonWM)/length(TumorD2_spots)

    #Tumor_D3
    TumorD3_spots <- sample_sub$Barcode[which(sample_sub$Hs_admix > 0.05 & sample_sub$Density_group == "D3")]
    TumorD3_spots_WM <- intersect(TumorD3_spots, WM_spots)
    TumorD3_spots_nonWM <- setdiff(TumorD3_spots, TumorD3_spots_WM) 
    df1$nr_WM[5] <- length(TumorD3_spots_WM)
    df1$nr_nonWM[5] <- length(TumorD3_spots_nonWM)
    df1$pct_WM[5] <- length(TumorD3_spots_WM)/length(TumorD3_spots)
    df1$pct_nonWM[5] <- length(TumorD3_spots_nonWM)/length(TumorD3_spots)

    #Tumor_D4
    TumorD4_spots <- sample_sub$Barcode[which(sample_sub$Hs_admix > 0.05 & sample_sub$Density_group == "D4")]
    TumorD4_spots_WM <- intersect(TumorD4_spots, WM_spots)
    TumorD4_spots_nonWM <- setdiff(TumorD4_spots, TumorD4_spots_WM) 
    df1$nr_WM[6] <- length(TumorD4_spots_WM)
    df1$nr_nonWM[6] <- length(TumorD4_spots_nonWM)
    df1$pct_WM[6] <- length(TumorD4_spots_WM)/length(TumorD4_spots)
    df1$pct_nonWM[6] <- length(TumorD4_spots_nonWM)/length(TumorD4_spots)

    cont_tab[[s]] <- df1
}
#saveRDS(cont_tab, "/work/morrissy_lab/vthoppey/SpatialData/Analysis/Overlap_with_WM/WM_nrSpots_perTimepoint.rds")

#cont_tab <- readRDS("/work/morrissy_lab/vthoppey/SpatialData/Analysis/Overlap_with_WM/WM_nrSpots_perTimepoint.rds")

pdf("WM_chisq_and_fisherexact_perTimpoint.pdf", width = 30, height = 5)
for(s in 1:length(cont_tab)){
    df1 <- cont_tab[[s]]
    df1 <- df1[,1:2]
    
    # Chi-square test
    chsq_DFlist <- vector("list",3)
    if(length(which(rowSums(df1) < 5)) > 0){
        chsq_DFlist[[1]] <- df1[3:6,]
        chsq_DFlist[[1]] <- chsq_DFlist[[1]][-which(rowSums(chsq_DFlist[[1]]) < 5),]
    }else{
        chsq_DFlist[[1]] <- df1[3:6,]
    }
    chsq_DFlist[[2]] <- df1[-which(rowSums(df1) <= 5 | rownames(df1) == "Brain"),]
    if(length(which(rowSums(df1) <= 5)) > 0){
        chsq_DFlist[[3]] <- df1[-which(rowSums(df1) <= 5),]
    }else{
        chsq_DFlist[[3]] <- df1
    }
    index = 3
    chsq_Plotlist <- vector("list",6)
    for(l in index){
        print(l)
        if(l == 1){
            sample_name = names(cont_tab)[s]
        }else{
            sample_name = ""
        }
        chisqTest1 <- chisq.test(chsq_DFlist[[which(index == l)]])
        obs1 <- data.frame(chisqTest1$observed) %>% add_column(value_type = "observed", group = rownames(chisqTest1$observed), .before = "nr_WM")
        exp1 <- data.frame(chisqTest1$expected) %>% add_column(value_type = "expected", group = rownames(chisqTest1$expected), .before = "nr_WM")
        v1 <- rbind(obs1, exp1) 
        v1_melt <- melt(v1, id.vars = c("value_type","group")) %>% filter(variable == "nr_WM")
        chsq_Plotlist[[l]] <- ggplot(v1_melt, aes(x = group, y = value, fill = value_type)) + 
            geom_bar(stat =  "identity", position = "dodge") + 
            theme_classic(base_size = 13) + 
            theme(legend.position = "bottom") + 
            ylab("count") + 
            labs(title = paste(sample_name), subtitle = paste("chi-sq. pval = ", signif(chisqTest1$p.value, digits = 3))) +  
            scale_fill_manual(values = c("#999999", "#E69F00"))  +
            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
        ## Chi-square test - v1 - residuals
        chisqTest1_res <- data.frame(chisqTest1$residuals) %>% add_column(group = rownames(chisqTest1$residuals), .before = "nr_WM")
        chisqTest1_res_melt <- melt(chisqTest1_res)
        chsq_Plotlist[[l+1]] <- ggplot(chisqTest1_res_melt, aes(x = group, y = value, fill = variable)) + 
                geom_bar(stat =  "identity", position = "dodge") + 
                theme_classic(base_size = 13) + 
                theme(legend.position = "bottom") + 
                ylab("residuals") + 
                labs(title = "", subtitle = paste("chi-sq. pval = ", signif(chisqTest1$p.value, digits = 3))) + 
                scale_fill_manual(values = c("red", "blue")) +
                theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
    }
    
    g1 <- ggarrange(plotlist = chsq_Plotlist, ncol = 6)
    print(g1)

    TAB_sub1 <- data.frame(Structure = "WM", Time_point = names(cont_tab)[s],
                            category = rownames(chsq_DFlist[[3]]), 
                            nr_WM_obs = v1$nr_WM[which(v1$value_type == "observed")], 
                            nr_nonWM_obs = v1$nr_nonWM[which(v1$value_type == "observed")], 
                            nr_WM_expected = v1$nr_WM[which(v1$value_type == "expected")],
                            nr_nonWM_expected = v1$nr_nonWM[which(v1$value_type == "expected")],
                            nr_WM_residuals = chisqTest1_res$nr_WM,
                            nr_nonWM_residuals = chisqTest1_res$nr_nonWM,
                            chisq_test_pvalue = p_value)
    TAB_all <- rbind(TAB_all, TAB_sub1)
}
dev.off()
write.table(TAB_all, "WM_chisq_all_samples_per_timepoint.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)