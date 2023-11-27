library(dplyr)
library(tidyr)
library(stringr)
library(ComplexHeatmap)
library(tibble)

# Tumor
# TME
# CellChat 


csv_file <- "gProfiler_Human_LR_allterms.csv" 
gProf_res <- read.csv(csv_file, stringsAsFactors = FALSE)

## Get columns 
gProf_res <- gProf_res[,c(1,2,4,grep("adjusted_p_value__", colnames(gProf_res)))]
colnames(gProf_res) <- gsub("adjusted_p_value__","",colnames(gProf_res))
gProf_res <- gProf_res %>% add_column(category = "", .after = "term_size") %>% filter(term_size > 10) %>% filter(term_size < 2000)

## Convert to -log10(adj_pval)
for(i in 5:ncol(gProf_res)){
    gProf_res[,i] <- -log10(gProf_res[,i])
}

## Select GO terms
gProf_res <- gProf_res %>% filter(source == "GO:BP")
gProf_res_sub <- gProf_res[,c(2, 5:ncol(gProf_res))]
rownames(gProf_res_sub) <- make.unique(gProf_res_sub$term_name)
gProf_res_sub$term_name <- NULL 


## Cap
gProf_res_sub <- as.matrix(gProf_res_sub)
gProf_res_sub[which(gProf_res_sub > 10)] <- 10

## Assign a group to each term bases on max -log10(adj_pval)
gProf_res_sub2 <- gProf_res_sub
gProf_res_sub2 <- data.frame(gProf_res_sub2) %>% add_column("group" = "", .before = colnames(gProf_res_sub2)[1])
for(i in 1:nrow(gProf_res_sub2)){
    gProf_res_sub2$LRgroup[i] <- names(which.max(gProf_res_sub2[i,][2:ncol(gProf_res_sub2)]))
}

## For each group, 
   # Select terms assigned to the group
   # Order (decreasing) by -log10(adj_pval)
LRgroups <- paste("LR",seq(1,8),sep="")
Terms_order <- data.frame() 
for(i in 1:length(LRgroups)){
    gProf_res_sub3 <- gProf_res_sub2 %>% filter(LRgroup == LRgroups[i]) 
    gProf_res_sub3 <- gProf_res_sub3[which(colnames(gProf_res_sub3) %in% c("LRgroup", LRgroups[i]))]
    gProf_res_sub3 <- gProf_res_sub3[order(gProf_res_sub3[,2], decreasing = TRUE),]
    Terms_order <- rbind(Terms_order, data.frame(V1 = rownames(gProf_res_sub3), V2 = LRgroups[i]))
}

## Order the original matrix by "term_order" obtained above
gProf_res_sub4 <- gProf_res_sub[Terms_order$V1,]


## Plot
col_fun <- RColorBrewer::brewer.pal(name = "Blues", n = 9)

pdf("Heatmap.pdf", height = 10, width = 6)
Heatmap(gProf_res_sub4, col = col_fun, cluster_rows = FALSE, show_row_names=FALSE, cluster_columns = FALSE)
dev.off()



### Hub gene pathway enrichment
library(clusterProfiler)
library(enrichplot)
library(DOSE) # needed to convert to enrichResult object
library(gprofiler2)
library(rutils)

gp_mod <- read.table("hub_gene_pathwayenrichment.txt", sep="\t", header = TRUE, stringsAsFactor=FALSE)
gp_mod <- gp_mod %>% filter(source == "GO:BP")

gp_mod$GeneRatio = gp_mod$intersection_size / gp_mod$query_size
gp_mod$BgRatio = paste0(gp_mod$term_size, "/", gp_mod$effective_domain_size)
names(gp_mod) = c("Cluster", "Category", "ID", "Description", "p.adjust", 
                    "query_size", "Count", "term_size", "effective_domain_size", 
                    "GeneRatio", "BgRatio")

# define as enrichResult object
gp_mod_enrich  = new("enrichResult", result = gp_mod)

plot_df_sub <- gp_mod_enrich@result %>% filter(term_size > 5) %>% filter(term_size < 500) %>% dplyr::select(Cluster, ID, term_size, Description, p.adjust, GeneRatio)

thresholds <- read.table("go_reduce_thresholds.txt", sep="\t", header = TRUE, stringsAsFactor=FALSE)
GEP <- paste("GEP_",seq(1,15),sep="")
final_res <- data.frame()
for(i in 1:length(GEP)){
    sub1 <- plot_df_sub %>% filter(Cluster == GEP[[i]])
    sub2 <- data.frame(go_id = sub1$ID, go_type = "BP")
    res <- go_reduce(sub2, orgdb="org.Hs.eg.db", threshold = thresholds$threshold[which(thresholds$GEP == GEP[[i]])]) 
    sub1 <- sub1[order(match(sub1$ID, res$go_id)),]
    res$p.adjust <- sub1$p.adjust
    res$Description <- sub1$Description
    res <- res %>% add_column(GEP = GEP[[i]], .before = "go_id")
    final_res <- rbind(final_res, res)
}

plot1 <- final_res %>% dplyr::select(GEP, parent_term, p.adjust) 
plot1 <- plot1 %>% count(GEP, parent_term)
colnames(plot1)[3] <- "nr_terms"
plot1 <- plot1[-which(is.na(plot1$parent_term)),]
plot1$parent_term <- factor(plot1$parent_term, levels = unique(plot1$parent_term))
plot1$GEP <- factor(plot1$GEP, levels = names(hub_genes))
colnames(plot1)[3] <- "nr_term"
plot1$nr_term <- log2(plot1$nr_term + 1)
pal <- wes_palette("Zissou1", 4, type = "continuous")
sc <- scale_colour_gradientn(colours = as.character(pal), limits=c(min(plot1$nr_term),max(plot1$nr_term)))

pdf("Dotplot.pdf", width = 10, height = 12)
p1 <- ggplot(plot1, aes(x = GEP, y = parent_term, size = nr_term, color = nr_term )) + 
      geom_point() + theme_minimal(base_size = 15) + sc + 
      theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))
print(p1)
dev.off()

