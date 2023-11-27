## GBM States (Neftel et al.,)
library(dplyr)
library(tidyr)
library(stringr)
library(pheatmap)
library(ggplot2)
library(pals)
library(reshape2)

####### OVERLAP 
THRESHOLDS <- data.frame(V1 = c(0.95, 0.8, 0.5, 0.2, 0.05, 0, 0.05, 0.1, 0.025), V2 = c(1.0, 0.95, 0.8, 0.5, 0.2, 0.05, 1, 1, 1))

MD <- read.csv("Spots_metadata.csv", stringsAsFactors=FALSE) 
rownames(MD) <- MD$Barcode
annot <- MD %>% select("MES", "AC", "OPC", "NPC")
MD <- MD[,-which(colnames(MD) %in% colnames(annot))]


SAMPLES <- unique(MD$Sample)
UMAT <- read.table("GRCh38_consensus_norm_scaled.txt",sep="\t", header = TRUE, stringsAsFactor=FALSE) 
colnames(UMAT) <- gsub("Usage","GEP",colnames(UMAT))
MD <- MD[rownames(UMAT),]
annot <- annot[rownames(UMAT),]

colFun <- rev(RColorBrewer::brewer.pal(6, "RdBu"))

for(S in 1:length(SAMPLES)){
    for(i in 1:nrow(THRESHOLDS)){
        
        spots <- MD %>% filter(Hs_admix > THRESHOLDS$V1[i]) %>% filter(Hs_admix <= THRESHOLDS$V2[i]) %>% filter(Sample == SAMPLES[S])
        print( paste(THRESHOLDS$V1[i]," < HsAdmix <= ", THRESHOLDS$V2[i], " - ", SAMPLES[S], " - ", length(spots$Barcode), " spots", sep=""))
        spots <- spots$Barcode
        if(length(spots) > 10){
            annot_sub <- annot[spots,]
            MD2 <- MD[spots,]
            UMAT_sample <- UMAT[spots,]

            OVERLAP <- data.frame(matrix(0,nrow= ncol(UMAT_sample), ncol = 4))
            colnames(OVERLAP) <- colnames(annot_sub)
            rownames(OVERLAP) <- colnames(UMAT_sample)

            UMAT_sample <- as.matrix(UMAT_sample)
            annot_sub <- as.matrix(annot_sub)
            for(j in 1:nrow(OVERLAP)){
                df1 <- UMAT_sample[,rownames(OVERLAP)[j]]
                df2 <- df1[which(df1 > 0.1)]
                spots1 <- names(df2) 
                for(k in 1:ncol(OVERLAP)){
                    df3 <- annot_sub[,colnames(OVERLAP)[k]]
                    df4 <- df3[which(df3 > 0.1)]
                    spots2 <- names(df4)
                    intersect <- length(intersect(spots1, spots2))
                    intersect_prop <- intersect/length(df2)
                    OVERLAP[j,k] <- intersect_prop
                }
            }

            for(r in 1:nrow(OVERLAP)){
                OVERLAP[r,][which(is.na(OVERLAP[r,]))] <- 0
            }

            if(sum(rowSums(OVERLAP)) != 0){
                dir.create(file.path("OVERLAP_PER_SAMPLE",SAMPLES[S]))
                pdf(file.path("OVERLAP_PER_SAMPLE", SAMPLES[S], paste("OVERLAP_HsAdmix_", THRESHOLDS$V1[i],"__", THRESHOLDS$V2[i], "_", SAMPLES[S], ".pdf", sep="")), width = 4, height = 6)
                p1 <- pheatmap(OVERLAP,cluster_rows = TRUE, cluster_cols = TRUE, main = paste(THRESHOLDS$V1[i]," < HsAdmix <= ", THRESHOLDS$V2[i], " \n ", SAMPLES[S],"(",length(spots)," spots)", sep=""),
                               fontsize = 7, fontsize_row = 7, fontsize_col = 7)
                #p1 <- Heatmap(matrix = OVERLAP, col = colFun, column_title = paste(THRESHOLDS$V1[i]," < HsAdmix <= ", THRESHOLDS$V2[i], " - ", SAMPLES[S],"(",length(spots)," spots)", sep=""))
                print(p1)
                dev.off()

                pdf(file.path("OVERLAP_PER_SAMPLE", SAMPLES[S], paste("OVERLAP_HsAdmix_", THRESHOLDS$V1[i],"__", THRESHOLDS$V2[i], "_", SAMPLES[S], "_unclustered.pdf", sep="")),  width = 3, height = 4)
                p1_2 <- pheatmap(OVERLAP,cluster_rows = FALSE, cluster_cols = FALSE, main = paste(THRESHOLDS$V1[i]," < HsAdmix <= ", THRESHOLDS$V2[i], " \n ", SAMPLES[S],"(",length(spots)," spots)", sep=""),
                                 fontsize = 4, fontsize_row = 7, fontsize_col = 7)
                #p1_2 <- Heatmap(matrix = OVERLAP, col = colFun,  column_title = paste(THRESHOLDS$V1[i]," < HsAdmix <= ", THRESHOLDS$V2[i], " - ", SAMPLES[S],"(",length(spots)," spots)", sep=""),
                #            cluster_columns=FALSE, cluster_rows = FALSE)
                print(p1_2)
                dev.off()
            }
        }
    }
}


#### Summary 
library(dplyr)
library(tidyr)
library(stringr)
library(pheatmap)
library(ggplot2)
library(pals)
library(reshape2)
library(openxlsx)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

# Load thresholds
THRESHOLDS <- data.frame(V1 = c(0.8, 0.5, 0.2, 0.05), V2 = c(1.0, 0.8, 0.5, 0.2, 0.05))

# Step1: Load usage matrix and assign GEP to each spot based on max value
UMAT <- read.table("GRCh38_consensus_norm_scaled.txt",sep="\t", header = TRUE, stringsAsFactors=FALSE) 
colnames(UMAT) <- gsub("Usage","GEP",colnames(UMAT)) 
UMAT$GEP_annot <- ""
for(i in 1:nrow(UMAT)){
    UMAT$GEP_annot[i] <- colnames(UMAT)[which.max(UMAT[i,])]
}

# Step2: Load file with baseline state for each sample
BL <- read.table("BaselineState.txt",sep="\t", header =TRUE, stringsAsFactor=FALSE)
BL$Baseline <- gsub(" ","", BL$Baseline)
BL_states <- setdiff(unique(BL$Baseline),"None")
#BL_states <- c("AC","OPC","NPC")

# Step3: Load metadata with state signature scores
MD <- read.csv("Spots_metadata.csv", stringsAsFactors=FALSE) 
rownames(MD) <- MD$Barcode
MD_sub <- MD %>% select("Barcode", "MES", "AC", "OPC", "NPC","Hs_admix", "Sample","Time_point")

# Step4: Subset metadata to desired Admix category
for(t in 1:nrow(THRESHOLDS)){
    
    MD_sub1 <- MD_sub %>% filter(Hs_admix > THRESHOLDS$V1[t]) %>% filter(Hs_admix <= THRESHOLDS$V2[t]) 

    for(BLS in BL_states){
        # Step5: Subset further to desired baseline state 
        MD_sub2 <- MD_sub1 %>% filter(Sample %in% BL$Sample[which(BL$Baseline == BLS)]) 

        # Step6: Make sure the dimensions of usage matrix and metadata match
        spots <- intersect(rownames(MD_sub2), rownames(UMAT))
        annot <- MD_sub2 %>% select("MES","AC","OPC","NPC","Time_point")
        annot_sub <- annot[spots,]
        UMAT_sub <- UMAT[,1:15][spots,]

        # Step7: For the given baseline, calculate overlap with GEPs per time point 
        MATRIX_LIST <- vector("list",3)
        names(MATRIX_LIST) <- c("early","mid","late")

        BLS_Transition <- vector("list",3)
        names(BLS_Transition) <- c("early","mid","late")

        for(TP in c("early","mid","late")){
            # Create overlap matrix layout
            OVERLAP <- data.frame(matrix(0,nrow= cNMF_dir$selected_k[8], ncol = 4))
            colnames(OVERLAP) <- colnames(annot)[1:4]
            rownames(OVERLAP) <- colnames(UMAT_sub)
    
            annot_TP <- annot_sub[which(annot_sub$Time_point == TP),]
            annot_TP <- as.matrix(annot_TP[,1:4])
            UMAT_TP <- as.matrix(UMAT_sub)[rownames(annot_TP),]
        
            # Calculate overlap
            for(j in 1:nrow(OVERLAP)){
                df1 <- UMAT_TP[,rownames(OVERLAP)[j]]
                df2 <- df1[which(df1 > 0.1)]
                spots1 <- names(df2) 
                for(k in 1:ncol(OVERLAP)){
                    df3 <- annot_TP[,colnames(OVERLAP)[k]]
                    df4 <- df3[which(df3 > 0.1)]
                    spots2 <- names(df4)
                    intersect <- length(intersect(spots1, spots2))
                    intersect_prop <- intersect/length(df2)
                    OVERLAP[j,k] <- intersect_prop
                }
            }
            for(r in 1:nrow(OVERLAP)){
                OVERLAP[r,][which(is.na(OVERLAP[r,]))] <- 0
                OVERLAP[r,] <- OVERLAP[r,]/max(OVERLAP[r,])
                OVERLAP[r,][which(is.na(OVERLAP[r,]))] <- 0
            }
            MATRIX_LIST[[TP]] <- OVERLAP
        
            ## Calculate transition for a given baseline state
            if(BLS %in% c("AC","OPC","NPC")){
                BLS_DF <- OVERLAP[,-which(colnames(OVERLAP) == BLS)]
                for(c in 1:ncol(BLS_DF)){
                    BLS_DF[,c] <- BLS_DF[,c]/OVERLAP[[BLS]]
                }
                colnames(BLS_DF) <- paste(BLS,"__",colnames(BLS_DF),sep="")
            }else{
                OVERLAP_merge <- OVERLAP
                OVERLAP_merge$AC_OPC <- (OVERLAP_merge$AC + OVERLAP_merge$OPC)/2
                BLS_DF <- OVERLAP_merge[,-which(colnames(OVERLAP_merge) %in% c("AC","OPC","AC_OPC"))]
                for(c in 1:ncol(BLS_DF)){
                    BLS_DF[,c] <- BLS_DF[,c]/OVERLAP_merge[["AC_OPC"]]
                }
                colnames(BLS_DF) <- paste("AC_OPC__",colnames(BLS_DF),sep="")
            }
            for(r2 in 1:nrow(BLS_DF)){
                BLS_DF[r2,][which(is.na(BLS_DF[r2,]))] <- 0
            }
            for(c2 in 1:ncol(BLS_DF)){
                BLS_DF[,c2][which(BLS_DF[,c2] == "Inf")] <- OVERLAP[[gsub(paste(BLS,"__",sep=""),"",colnames(BLS_DF)[c2])]][which(BLS_DF[,c2] == "Inf")]
            }
        
            BLS_Transition[[TP]] <- BLS_DF
        } 
        saveRDS(BLS_Transition, paste("T0.1__",BLS, "___HsAdmix", THRESHOLDS$V1[t],"_", THRESHOLDS$V2[t],"__TransitionMatrix.rds",sep=""))
    }
}



#### Transitions
annot <- read.table("/work/morrissy_lab/vthoppey/SpatialData/Analysis/cNMF_new/pdx_merge_all/Mm_Hs_GEP_Annotations_final.txt",sep="\t",header = TRUE, stringsAsFactors=FALSE) 
annot <- annot %>% filter(Annotation1 == "cell_state")
annot$cNMF_GEP <- gsub("Hs_","", annot$cNMF_GEP)
annot$label <- paste(annot$cNMF_GEP,annot$Annotation2,sep="__")

THRESHOLDS <- data.frame(V1 = c(0.8, 0.5, 0.2, 0.05), V2 = c(1.0, 0.8, 0.5, 0.2, 0.05))
THRESHOLDS$label <- paste("T",seq(1,nrow(THRESHOLDS)),sep="")

BL <- read.table("BaselineState.txt",sep="\t", header =TRUE, stringsAsFactor=FALSE)
BL$Baseline <- gsub(" ","", BL$Baseline)
BL_states <- setdiff(unique(BL$Baseline),"None")
DF <- data.frame(GEP = paste("GEP_",seq(1,15),sep=""))

for(t in 1:nrow(THRESHOLDS)){
    for(BLS in BL_states){
        list1 <- readRDS(paste("T0.05__", BLS, "___HsAdmix", THRESHOLDS$V1[t],"_", THRESHOLDS$V2[t],"__TransitionMatrix.rds",sep=""))
        NTrans = ncol(list1[[1]])
        for(a in 1:length(list1)){
            colnames(list1[[a]]) <- paste(THRESHOLDS$label[t],names(list1)[a], colnames(list1[[a]]),sep="_")
        }
        for(b in 1:NTrans){
            for(c in 1:length(list1)){
                DF <- cbind(DF, list1[[c]][b])
            }
        }
    }
}

Transitions <- expand.grid(BL_states, c("OPC","AC","NPC","MES")) 
Transitions$Var1 <- as.character(Transitions$Var1) ; Transitions$Var2 <- as.character(Transitions$Var2)
Transitions <- Transitions[-which(Transitions$Var1 == Transitions$Var2),] 
Transitions$label <- paste(Transitions$Var1, Transitions$Var2, sep="__")
Transitions <- Transitions[-which(Transitions$label %in% c("AC_OPC__AC","AC_OPC__OPC")),]

Transitions <- Transitions %>% filter(Var1 == "AC")

Transitions_TP <- data.frame(V1 = c( paste( c("early"), Transitions$label,sep="_"),
                               paste( c("mid"), Transitions$label,sep="_"),
                               paste( c("late"), Transitions$label,sep="_")),
                          V2 = rep(Transitions$label,3))
Transitions_TP_THRESH <- c( paste("T1", Transitions_TP$V1, sep="_"),
                            paste("T2", Transitions_TP$V1, sep="_"),
                            paste("T3", Transitions_TP$V1, sep="_"),
                            paste("T4", Transitions_TP$V1, sep="_"))
DF_sub <- DF[,Transitions_TP_THRESH]
rownames(DF_sub) <- annot$label 

metadata <- data.frame(V1 = colnames(DF_sub), V2 = colnames(DF_sub))
metadata <- metadata %>% separate(V2, c("X1","TransitionTo"),sep="__")
metadata$X1 <- gsub("AC_OPC","ACOPC",metadata$X1)
metadata <- metadata %>% separate(X1, c("Cellularity","Time_point","Baseline"),sep="_")
metadata$Baseline <- gsub("ACOPC","AC_OPC",metadata$Baseline)
rownames(metadata) <- metadata$V1
metadata$V1 <- NULL
metadata <- metadata[colnames(DF_sub),]
metadata <- metadata[order(match(metadata$Baseline, sort(metadata$Baseline) )),] 
metadata <- metadata[order(match(metadata$Baseline, sort(metadata$TransitionTo) )),] 
DF_sub <- DF_sub[,rownames(metadata)]



#### Barplots for transitions
DF1 <- t(DF_sub)
DF2 <- DF1 #[,c(4,11)]
DF2 <- data.frame(DF2) %>% add_column("V1" = rownames(DF2), .before = "GEP_1__OC_OP")
DF2 <- DF2 %>% separate(V1, c("V2","EndState"),sep="__")
DF2$V2 <- gsub("AC_OPC","AC|OPC",DF2$V2)
DF2 <- DF2 %>% separate(V2, c("Density","TimePoint","Baseline"),sep="_")
DF2$Baseline[which(DF2$Baseline == "AC|OPC")] <- "AC_OPC"

#GEP = "GEP_4__AC3_PAN_Ribosome_Translation"
Den <- unique(DF2$Density)

DF3 <- DF2 #%>% filter(Density == Den[d]) 
DF3$xaxis <- paste(DF3$Baseline,DF3$EndState,sep="->")
DF3$EndState <- NULL
DF4 <- melt(DF3) #%>% filter(variable == GEP)

DF4$variable <- as.character(DF4$variable)
DF4 <- DF4 %>% separate(variable, c("X1","X2"), sep="__")
DF4$X2 <- NULL; DF4$X1 <- gsub("GEP_","hGEP",DF4$X1)
colnames(DF4)[which(colnames(DF4) == "X1")] <- "GEP"

DF4$GEP <- factor(DF4$GEP, levels = paste("hGEP",seq(1,15),sep=""))
DF4$TimePoint <- factor(DF4$TimePoint, levels = c("early","mid","late"))
DF4$value2 <- DF4$value ; DF4$value[which(DF4$value > 2)] <- 2
DF4$xaxis <- paste(DF4$Density, DF4$xaxis,sep="__")

#DF4_sub <- DF4 %>% filter(Baseline == "OPC")
p1 <- ggplot(DF4, aes(x = xaxis, y = value, fill = TimePoint)) + 
        geom_bar(position="dodge", stat="identity") + 
        theme_classic(base_size = 13) +
        facet_grid(GEP~Baseline+Density, scales = "free_x") + 
        scale_fill_manual(values = c("early" = "#FF0000", "mid" = "#00A08A", "late" = "#F2AD00")) + 
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))  + 
        theme(panel.spacing.x = unit(2, "lines"))
pdf("Barplot_TransitionScore.pdf", width = 18, height = 24)
print(p1)
dev.off()
 
# pdf("test3.pdf", width = 15, height = 18)
# ggarrange(p1,p2,p3,ncol=3,common.legend = TRUE)
# dev.off()

