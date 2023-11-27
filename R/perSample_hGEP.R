library(dplyr)
library(tidyr)
library(stringr)
library(ComplexHeatmap)
library(ggplot2)
library(tibble)
library(reshape2)
library(openxlsx)
library(graphics)
library(vcd)
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

hCon <- read.table("GRCh38_consensus_norm_scaled.txt", sep="\t", header = TRUE, stringsAsFactor=FALSE)
colnames(hCon) <- gsub("Usage_","hGEP",colnames(hCon))

### All Samples combined
GEPs <- colnames(hCon)
cont_tab <- vector("list",length(GEPs))
names(cont_tab) <- GEPs


### ALL GEPs - per sample
GEPs <- colnames(hCon)
for(NRSPOT in 50){
    for(thresh in 0.01){

        TAB_all <- data.frame()

        for(g in 1:length(GEPs)){

            GEP_spots <- names(which(as.matrix(hCon)[,GEPs[g]] > thresh)) ## USAGE THRESHOLD
            barcode_sub <- MD %>% filter(Barcode %in% GEP_spots)
            
            sample_freq <- data.frame(table(barcode_sub$Sample_RepMerged ))
            sample_freq_sub <- sample_freq %>% filter(Freq > NRSPOT) ## NR SPOT THRESHOLD
            samples <- as.character(sample_freq_sub$Var1)

            MD_subset <- MD %>% filter(Sample_RepMerged %in% samples) 

            cont_tab <- vector("list",length(samples))
            names(cont_tab) <- samples
            
            samples <- names(cont_tab)
            for(s in 1:length(samples)){
                
                df1 <- data.frame(matrix(0, nrow = 4, ncol = 4, dimnames = list(c("Tumor_D1","Tumor_D2","Tumor_D3","Tumor_D4"), c("nr_GEP","nr_nonGEP", "pct_GEP","pct_nonGEP"))))
                ss_sample <- MD %>% filter(Sample_RepMerged == samples[s])
                
                hCon_sub <- hCon[intersect(ss_sample$Barcode,rownames(hCon)),]
                GEP_spots <- rownames(hCon_sub)[which(hCon_sub[,GEPs[g]] > thresh)] ## USAGE THRESHOLD

                #Tumor_D1
                TumorD1_spots <- ss_sample$Barcode[which(ss_sample$Density_group == "D1")]
                TumorD1_spots_GEP <- intersect(TumorD1_spots, GEP_spots)
                TumorD1_spots_nonGEP <- setdiff(TumorD1_spots, TumorD1_spots_GEP)
                df1$nr_GEP[1] <- length(TumorD1_spots_GEP)
                df1$nr_nonGEP[1] <- length(TumorD1_spots_nonGEP)
                df1$pct_GEP[1] <- length(TumorD1_spots_GEP)/length(TumorD1_spots)
                df1$pct_nonGEP[1] <- length(TumorD1_spots_nonGEP)/length(TumorD1_spots)

                #Tumor_D2
                TumorD2_spots <- ss_sample$Barcode[which(ss_sample$Density_group == "D2")]
                TumorD2_spots_GEP <- intersect(TumorD2_spots, GEP_spots)
                TumorD2_spots_nonGEP <- setdiff(TumorD2_spots, TumorD2_spots_GEP)    
                df1$nr_GEP[2] <- length(TumorD2_spots_GEP)
                df1$nr_nonGEP[2] <- length(TumorD2_spots_nonGEP)
                df1$pct_GEP[2] <- length(TumorD2_spots_GEP)/length(TumorD2_spots)
                df1$pct_nonGEP[2] <- length(TumorD2_spots_nonGEP)/length(TumorD2_spots)

                #Tumor_D3
                TumorD3_spots <- ss_sample$Barcode[which(ss_sample$Density_group == "D3")]
                TumorD3_spots_GEP <- intersect(TumorD3_spots, GEP_spots)
                TumorD3_spots_nonGEP <- setdiff(TumorD3_spots, TumorD3_spots_GEP) 
                df1$nr_GEP[3] <- length(TumorD3_spots_GEP)
                df1$nr_nonGEP[3] <- length(TumorD3_spots_nonGEP)
                df1$pct_GEP[3] <- length(TumorD3_spots_GEP)/length(TumorD3_spots)
                df1$pct_nonGEP[3] <- length(TumorD3_spots_nonGEP)/length(TumorD3_spots)

                #Tumor_D4
                TumorD4_spots <- ss_sample$Barcode[which(ss_sample$Density_group == "D4")]
                TumorD4_spots_GEP <- intersect(TumorD4_spots, GEP_spots)
                TumorD4_spots_nonGEP <- setdiff(TumorD4_spots, TumorD4_spots_GEP) 
                df1$nr_GEP[4] <- length(TumorD4_spots_GEP)
                df1$nr_nonGEP[4] <- length(TumorD4_spots_nonGEP)
                df1$pct_GEP[4] <- length(TumorD4_spots_GEP)/length(TumorD4_spots)
                df1$pct_nonGEP[4] <- length(TumorD4_spots_nonGEP)/length(TumorD4_spots)

                cont_tab[[s]] <- df1
            }
            #saveRDS(cont_tab, paste(GEPs[g],"_nrSpots_perSample_usagethreshold", thresh, "_nrSpotThresh", NRSPOT, ".rds",sep=""))

            pdf(paste(GEPs[g],"_",GEP_name,"_chisq_and_fisherexact_perSample_usagethreshold", thresh, "_nrSpotThresh", NRSPOT, ".pdf",sep=""), width = 15, height = 5)
            for(s in 1:length(cont_tab)){
                df1 <- cont_tab[[s]]
                df1 <- df1[,1:2]
                
                # Chi-square test
                chsq_DFlist <- vector("list",1)
                zero_rows <- c()
                for(r in 1:nrow(df1)){
                    if(df1[r,][1] == 0 | df1[r,][2] == 0){
                        zero_rows <- c(zero_rows, r)
                    }
                }
                if(length(zero_rows > 0)){
                    df1 <- df1[-zero_rows,]
                }else{
                    df1 <- df1
                }

                if(nrow(df1) > 0){
                    chsq_DFlist[[1]] <- df1

                    index = c(1)
                    chsq_Plotlist <- vector("list",2)
                    for(l in index){
                        print(l)
                        if(l == 1){
                            sample_name = names(cont_tab)[s]
                        }else{
                            sample_name = ""
                        }
                        chisqTest1 <- chisq.test(chsq_DFlist[[which(index == l)]])
                        if(nrow(chsq_DFlist[[which(index == l)]]) == 1){
                            obs1 <- chsq_DFlist[[which(index == l)]] %>% add_column(value_type = "observed", group = rownames(chsq_DFlist[[which(index == l)]]), .before = "nr_GEP")
                            exp1 <- chsq_DFlist[[which(index == l)]] %>% add_column(value_type = "expected", group = rownames(chsq_DFlist[[which(index == l)]]), .before = "nr_GEP")
                        }else{
                            obs1 <- data.frame(chisqTest1$observed) %>% add_column(value_type = "observed", group = rownames(chisqTest1$observed), .before = "nr_GEP")
                            exp1 <- data.frame(chisqTest1$expected) %>% add_column(value_type = "expected", group = rownames(chisqTest1$expected), .before = "nr_GEP")
                        }
                        v1 <- rbind(obs1, exp1) 
                        v1_melt <- melt(v1, id.vars = c("value_type","group")) %>% filter(variable == "nr_GEP")
                        chsq_Plotlist[[l]] <- ggplot(v1_melt, aes(x = group, y = value, fill = value_type)) + 
                            geom_bar(stat =  "identity", position = "dodge") + 
                            theme_classic(base_size = 13) + 
                            theme(legend.position = "bottom") + 
                            ylab("count") + 
                            labs(title = paste(sample_name), subtitle = paste("chi-sq. pval = ", signif(chisqTest1$p.value, digits = 3))) +  
                            scale_fill_manual(values = c("#999999", "#E69F00"))  +
                            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
                        ## Chi-square test - v1 - residuals
                        if(nrow(chsq_DFlist[[which(index == l)]]) == 1){
                            chisqTest1_res <- chsq_DFlist[[which(index == l)]]
                            chisqTest1_res[1] <- chisqTest1$residuals[1]; chisqTest1_res[2] <- chisqTest1$residuals[2]
                            chisqTest1_res <- chisqTest1_res %>% add_column(group = rownames(chsq_DFlist[[which(index == l)]]), .before = "nr_GEP")
                        }else{
                            chisqTest1_res <- data.frame(chisqTest1$residuals) %>% add_column(group = rownames(chisqTest1$residuals), .before = "nr_GEP")
                        }
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
                }else{
                    full <- cont_tab[[s]][,1:2]
                    full <- full %>% add_column(value_type = "observed", group = rownames(full) ,before = "nrGEP")
                    full_melt <- melt(full)
                    chsq_Plotlist[[l]] <- ggplot(full_melt, aes(x = group, y = value, fill = variable)) + 
                            geom_bar(stat =  "identity", position = "dodge") + 
                            theme_classic(base_size = 13) + 
                            theme(legend.position = "bottom") + 
                            ylab("nrSpots") + 
                            labs(title = "", subtitle = paste("chi-sq. not performed",sep="")) + 
                            theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
                    chsq_Plotlist[[l+1]] <- ggplot() + theme_void()
                }

                
                g1 <- ggarrange(plotlist = chsq_Plotlist, ncol = 2)
                print(g1)

                TAB_sub1 <- data.frame(GEP = GEPs[g], sample = names(cont_tab)[s],
                            category = rownames(df1), 
                            nr_GEP_obs = v1$nr_GEP[which(v1$value_type == "observed")], 
                            nr_nonGEP_obs = v1$nr_nonGEP[which(v1$value_type == "observed")], 
                            nr_GEP_expected = v1$nr_GEP[which(v1$value_type == "expected")],
                            nr_nonGEP_expected = v1$nr_nonGEP[which(v1$value_type == "expected")],
                            nr_GEP_residuals = chisqTest1_res$nr_GEP,
                            nr_nonGEP_residuals = chisqTest1_res$nr_nonGEP,
                            chisq_test_pvalue = p_value)
                TAB_all <- rbind(TAB_all, TAB_sub1) 
            }
            dev.off()
        }
    }
}
write.table(TAB_all, "hGEP_chisq_individual_samples_per_program.txt", sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

