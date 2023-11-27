# Load libraries
library(dplyr)
library(tidyr)
library(stringr)
library(mclust)
library(ggplot2)
library(gridExtra)


# Function 
TopScoringGenes <- function(gene_spectra_score = gene_spectra_score, rank = rank, output_dir = output_dir, file_prefix = file_prefix){
    
    GS <- read.table(gene_spectra_score,sep="\t", header = TRUE, stringsAsFactors=FALSE, row.names = 1)
    rownames(GS) <- paste("GEP_",seq(1,nrow(GS)),sep="")
    GS <- t(GS)
    for(y in 1:ncol(GS)){
        GS[,y][which(GS[,y] < 0)] <- 0  ## Set negative values to zero 
    }
    GS <- as.matrix(GS)
    K = rank 

    GENE_LIST <- vector("list",K)
    PLOT_LIST <- vector("list",K)

    for(i in 1:K){
        mod <- densityMclust(GS[,i][which(GS[,i] != 0)], plot = FALSE)
        df <- data.frame(Usage = GS[,i][which(GS[,i] != 0)], Density = mod$density) 

        df2 <- data.frame(data = mod$data, classification = mod$classification)
        divisions <- c()
        for(d in sort(unique(mod$classification))){
            df3 <- df2 %>% filter(classification == d)
            max_val <- max(df3$data)
            divisions <- c(divisions, max_val)
        }
        perc95 <- quantile(df$Usage, 0.95)
        p <- ggplot(df,aes(x = Usage, y = Density)) + 
                    geom_histogram(aes(y=..density..), colour="white", fill="grey",bins = 50) +
                    geom_density(aes(x = Usage, y = Density), stat = 'identity') + theme_classic() +
                    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size=2)) + 
                    xlab(colnames(GS)[i]) + ylab("density") + 
                    geom_vline(xintercept = divisions, linetype="dotted", size = 0.3, color = "red") + 
                    geom_vline(xintercept = perc95, linetype="solid", size = 0.3, color = "blue") + 
                    geom_text(data = data.frame(x = c(round(divisions,7),round(perc95,7)), y = rep(max(df$Density), length(divisions) + 1)), 
                    mapping = aes(x = x, y = y, label = x),
                    inherit.aes = FALSE,
                    vjust = 1.5, angle=90, size = 2)
        
        PLOT_LIST[[i]] <- p 
        GENE_LIST[[i]] <- rownames(df)[which(df$Usage > sort(divisions)[length(divisions) - 1])]
    }
        
    if(K <= 5){
        ncol = K; nrow = 1
    }else{
        ncol = 5; nrow = ceiling(K/5)
    }
    
    pdf(file.path(output_dir, paste(file_prefix,"_topscoringgenes_mclustplot.pdf",sep="")), width = 5*ncol + 5, height = 5*nrow)
    G <- grid.arrange(grobs = PLOT_LIST, top = paste(file_prefix,"_",K,sep=""), ncol = ncol, nrow = nrow)
    print(G)
    dev.off()
    
    names(GENE_LIST) <- colnames(GS)
    saveRDS(GENE_LIST, file.path(output_dir, paste(file_prefix,"_topscoringgenes.rds",sep="")))
}

# Run function
TopScoringGenes(gene_spectra_score = gene_spectra_score, rank = rank, output_dir = output_dir, file_prefix = file_prefix)


