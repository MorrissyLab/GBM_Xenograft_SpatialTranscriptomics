## Spatial concordance / overlap

calculate_spatial_concordance <- function(usage_matrix_1 = path_to_usage_matrix_1, usage_matrix_2 = NULL, threshold_1 = threshold_1, threshold_2 = NULL){
    
    matrix_1 <- usage_matrix_1
    
    if(is.null(usage_matrix_2)){
        matrix_2 <- matrix_1
    }else{
        matrix_2 <- usage_matrix_2
    }

    if(is.null(threshold_2)){
        threshold_2 <- threshold_1
    }else{
        threshold_2 <- threshold_2
    }

    matrix_1 <- as.matrix(matrix_1) ; matrix_2 <- as.matrix(matrix_2)
    OVERLAP <- data.frame(matrix(0,nrow=ncol(matrix_1), ncol = ncol(matrix_2), dimnames = list(colnames(matrix_1), colnames(matrix_2))))

    for(i in 1:nrow(OVERLAP)){
        df1 <- matrix_1[,rownames(OVERLAP)[i]]
        df2 <- df1[which(df1 > threshold_1)]
        spots1 <- names(df2) 
        for(j in 1:ncol(OVERLAP)){
            df3 <- matrix_2[,colnames(OVERLAP)[j]]
            df4 <- df3[which(df3 > threshold_2)]
            spots2 <- names(df4)
            intersect <- length(intersect(spots1, spots2))
            intersect_prop <- intersect/length(df2)
            OVERLAP[i,j] <- intersect_prop
        }
    }

    return(OVERLAP)
}

# GRCh38_overlap_matrix <- calculate_usage_overlap(usage_matrix = "GRCh38_consensus_norm_scaled.txt", usage_matrix_2 = NULL, threshold_1 = 0.1, threshold_2 = NULL )    
# mm10_overlap_matrix <- calculate_usage_overlap(usage_matrix_1 = "mm10_consensus_norm_scaled.txt", usage_matrix_2 = NULL, threshold_1 = 0.05, threshold_2 = NULL )       
# GRCh38_mm10_overlap_matrix <- calculate_usage_overlap(usage_matrix = "GRCh38_consensus_norm_scaled.txt", usage_matrix_2 = "mm10_consensus_norm_scaled.txt", threshold_1 = 0.1, threshold_2 = 0.05)   
                                          