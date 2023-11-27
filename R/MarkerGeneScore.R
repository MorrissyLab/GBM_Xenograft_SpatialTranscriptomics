## Marker gene score

calculate_marker_gene_score <- function(gene_spectra_score_file, geneset, scale = TRUE){

    ## Get top 50 genes per GEP
    gene_score <- read.table(gene_spectra_score_file, sep="\t", header = TRUE, stringsAsFactors=FALSE, row.names = 1) 
    rownames(gene_score) <- paste("GEP_",seq(1,nrow(gene_score)),sep="")
    gene_score <- t(gene_score)
    top50 <- data.frame(V1 = seq(1, 50))
    for(i in 1:ncol(gene_score)){
        df <- data.frame(genes = rownames(gene_score), value = gene_score[,i]) 
        df <- df[order(match(df$value, sort(df$value, decreasing=TRUE) )),]
        genes <- gsub("GRCh38_","",genes) 
        genes <- gsub("mm10___","",genes)
        df2 <- data.frame(V1 = genes)
        colnames(df2) <- paste("GEP_",i,sep="")
        top50 <- cbind(top50, df2)
    }
    top50[,1] <- NULL
    
    ## Read gmt
    geneset <- GSEABase::getGmt(geneset_file)
    
    ## Format gene score 
    rownames(gene_score) <- paste("Usage_",seq(1, nrow(gene_score)),sep="")
    colnames(gene_score) <- gsub("GRCh38_","",colnames(gene_score))
    colnames(gene_score) <- gsub("mm10___","",colnames(gene_score))
    gene_score <- t(gene_score)
    gene_score <- data.frame(gene_score)
    

    geneset_names <- names(geneset)
    for(G in 1:length(geneset_names)){
        GS_name <- geneset_names[G]
        marker_genes <- gsub(" ","",unique( geneset[[GS_name]]@geneIds ))
        df1 <- data.frame(cell_type = GS_name, n_markers = length(marker_genes))
        
        SUB1 <- data.frame(Index = 1)
        for(U in 1:ncol(gene_score)){
            int_N <- length(intersect(top50[,U], marker_genes ))
            ranks <- 1/(which(top50[,U] %in% intersect(top50[,U], marker_genes )))
            ranks_sum <- sum(ranks)
            marker_score <- int_N*ranks_sum 
            marker_score <- marker_score/round(log10(length(marker_genes)))

            SUB2 <- data.frame(topN_marker_score = marker_score)
            colnames(SUB2) <- paste("GEP",U,"_",colnames(SUB2), sep="")
            colnames(SUB2) <- gsub("topN",paste("top",50,sep=""), colnames(SUB2))
            DF_SUB1 <- cbind(DF_SUB1,SUB2)
    }
    DF_SUB1 <- cbind(df1, DF_SUB1)
    FINAL_DF <- rbind(FINAL_DF, DF_SUB1)
    if(scale){
        if(length(which(rowSums(FINAL_DF) == 0)) > 0){
            FINAL_DF2 <- FINAL_DF[-which(rowSums(FINAL_DF) == 0),]
        }else{
            FINAL_DF2 <- FINAL_DF
        }
        FINAL_DF_scaled <- scale(t(FINAL_DF2))
        return(FINAL_DF_scaled)
    }else{
        return(FINAL_DF)
    }
    ## Plot
    #pheatmap(FINAL_DF, fontsize = 15)
    #pheatmap(FINAL_DF_scaled, fontsize = 15)
}
