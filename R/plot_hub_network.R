library(dplyr)
library(tidyr)
library(stringr)
library(STRINGdb)
library(igraph)
library(ggplot2)

# for(i in 1:15){
#     plot_hub_network(gene_list = 'Hs_all_15_GeneSpectraScore_95perc_list.rds', 
#                      GEP = i, 
#                      corr_file = "Hs_all_GEP_Corr_with_Expr.rds", 
#                      confidence = 400, 
#                      int_type = "all")
# }


plot_hub_network <- function(gene_list, GEP, corr_file, confidence , combine, int_type = c("known","all")){

    int_list  <- list("known" = c("experimentally_determined_interaction", "database_annotated"),
                 "all" = c("experimentally_determined_interaction", "database_annotated", "known_int", "neighborhood_on_chromosome", "gene_fusion", "phylogenetic_cooccurrence", "homology", "coexpression", "automated_textmining"))
    
    ### Step1: Create a STRINGdb object of genes in a GEP with
    ###           confidence value C and a given interaction type (Known OR ALL)
    #Step1a : Gene genes in a GEP
    GENE_LIST <- readRDS(gene_list)
    GENE_LIST <- GENE_LIST[GEP]
    
    #Step1b : Select interaction type
    type <- int_list[[int_type]]
    
    #Step1c : Create a STRINGdb object - this will include all interaction types
    string_db <- STRINGdb$new( version="11.5", species=9606, score_threshold=confidence, network_type="full", input_directory="")
    
    #Step1d : Create a dataframe of all genes in a GEP 
    tab_gene <- data.frame(gene = GENE_LIST) ## Creating a object to map information
    colnames(tab_gene) <- "gene"
    
    #Step1e : Select genes of a given interaction type (Known OR ALL)
    interac_tab <- read.table(paste(gsub("_","", names(GENE_LIST)), "/interaction_files/", gsub("_","", names(GENE_LIST)), "_string_interactions_c",confidence,"_",int_type,".tsv",sep="") ,sep="\t", header = TRUE, stringsAsFactors=FALSE, comment.char = "")
    
    colnames(interac_tab)[1] <- "node1"
    interac_tab_subset <- interac_tab[which(rowSums(interac_tab[, which( colnames(interac_tab) %in% type )]) != 0),] %>% filter(combined_score >= confidence/1000)
    tab_gene_sub <- tab_gene %>% filter(gene %in% unique(c(interac_tab_subset$node1, interac_tab_subset$node2)))
    
    #Step1f : Map data to create a STRINGdb network of given confidence in the GEP
    data_mapped <- string_db$map( tab_gene_sub, "gene", removeUnmappedRows = TRUE )
    STRING_g <- string_db$get_subnetwork(data_mapped$STRING_id)   ## STRING_g has the network to work with. 
    clustersList <- string_db$get_clusters(data_mapped$STRING_id, algorithm = "fastgreedy")
    for(cl1 in 1:length(clustersList)){
        sublist <- clustersList[[cl1]]
        for(cl2 in 1:length(sublist)){
            clustersList[[cl1]][cl2] <- data_mapped$gene[which(data_mapped$STRING_id == sublist[cl2])]
        }
    }
    clustersList <- clustersList[order(sapply(clustersList,length), decreasing=TRUE)]
    #names(clustersList) <- paste("C",seq(1,length(clustersList)),sep="")
    names(clustersList) <- paste("CC",sprintf("%02d",seq(1,length(clustersList))),sep="_")
    
    #Step1g: Start a dataframe to include info about genes (nodes)
    DF <- data.frame()
    for(d in 1:length(clustersList)){
        df1 <- data.frame(community = names(clustersList)[d], genes = clustersList[[d]], community.size = length(clustersList[[d]]))
        DF <- rbind(DF, df1)
    }


    ### Step2: Create a igraph network with gene names (instead of ensembl ids), subsetting to interactions of the required type (Known OR All)
    #Step2a : Edges
    e <- data.frame(get.edgelist(STRING_g))
    for(r1 in 1:nrow(e)){
        e$X1[r1] <- data_mapped$gene[which(data_mapped$STRING_id == e$X1[r1])]
        e$X2[r1] <- data_mapped$gene[which(data_mapped$STRING_id == e$X2[r1])]
    }
    colnames(e) <- c("From","To")
    #Step2b : Nodes
    v <- data.frame(Name = V(STRING_g)$name)
    for(r2 in 1:nrow(v)){
        v$Name[r2] <- data_mapped$gene[which(data_mapped$STRING_id == v$Name[r2])]
    }
    #Step2d : Create an igraph object
    g <- graph.empty()
    g <- add.vertices(g, nrow(v), 
	         name=as.character(v[,1]))
    g <- add.edges(g,t(e))
    g <- as.undirected(g)


    ### Step3: Select hub genes in the network based on hub_score() and add to dataframe DF
    #Step3a: Select hub genes in the network based on page_rank()
    hub_score <- page_rank(g, vids = V(g),directed = FALSE) ; hub_score <- sort(hub_score$vector, decreasing = TRUE) 
    hub_genes <- names( hub_score[which(hub_score > quantile(hub_score, 0.90))] ) 
    hub_score <- sort(degree(g, v = V(g), mode = "all", loops = TRUE, normalized = TRUE), decreasing = TRUE)
    hub_genes <- names(hub_score[which(hub_score > quantile(hub_score, 0.90))])
    
    #Step3b: add to dataframe DF
    DF$hub_score <- hub_score[DF$genes]
    DF$degree <- degree(g)[DF$genes]
    DF$is.hub <- ""
    DF$is.hub[which(DF$genes %in% hub_genes)] <- "YES"
    if(combine == "FALSE"){     ## default : write file
        write.table(DF, paste("/work/morrissy_lab/vthoppey/SpatialData/Analysis/Figure5_DruggableLandscape/Density_based/Human/iGraph_networks/",names(GENE_LIST),"/",names(GENE_LIST), "__", int_type,"_int_c",confidence,"__nodeinfo.txt",sep=""), sep="\t", quote = FALSE, col.names = TRUE, row.names = FALSE) 
    }else{
        write.table(DF, paste("/work/morrissy_lab/vthoppey/SpatialData/Analysis/Figure5_DruggableLandscape/Density_based/Human/iGraph_networks/",names(GENE_LIST),"/",names(GENE_LIST), "__", int_type,"_int_c",confidence,"_combinedhub__nodeinfo.txt",sep=""), sep="\t", quote = FALSE, col.names = TRUE, row.names = FALSE) 
    }
    ### Step4: Create attributes dataframes of the edges and nodes
    #Step4a : Create a dataframe of edge width (interaction scores)
    #edgeSize <- data.frame(get.edgelist(g))
    #edgeSize$score <- 0.001
    #for(r3 in 1:nrow(edgeSize)){
    #    sub_df <- interac_tab_subset %>% filter(node1 == edgeSize$X1[r3] | node2 == edgeSize$X1[r3]) %>% filter(node2 == edgeSize$X2[r3] | node1 == edgeSize$X2[r3])
    #    if(nrow(sub_df) > 0){
    #        edgeSize$score[r3] <- unique(sub_df$combined_score)
    #    }
    #}

    #Step4b : Creat a dataframe for nodes with size, color, shape, label, community information 
    corr <- readRDS(corr_file)
    corr <- cbind(data.frame(gene_name = rownames(corr)), corr)

    corr_sub <- corr[,which(colnames(corr) %in% c("gene_name", names(GENE_LIST)))]
    rownames(corr_sub) <- corr$gene_name
    v$size <- round((corr_sub[v$Name,][,2])*10)
    v$shape <- "circle"
    community_colors <- rainbow(max(length(clustersList)))
    names(community_colors) <- names(clustersList)
    v$community <- ""
    v$color_community <- ""
    for(cl3 in 1:length(clustersList)){
        v$community[which(v$Name %in% clustersList[[cl3]])] <- names(clustersList)[cl3]
        v$color_community[which(v$Name %in% clustersList[[cl3]])] <- community_colors[[cl3]]
    }
    v$color_hub <- "gold"
    v$color_hub[which(v$Name %in% hub_genes)] <- "coral1"
    v$label1 <- v$Name
    v$label2 <- v$Name
    v$label2[which(v$color_hub != "coral1")] <- ""
    v$label3 <- ""
 
    
    ### Step5: Create the igraph object 
    g <- graph.empty()
    #g1
    g1 <- add.vertices(g, nrow(v), 
	         name=as.character(v$Name), 
	         color=as.character(v$color_community), 
	         shape=as.character(v$shape), 
	         size=as.character(v$size),
             label = as.character(v$label1) )
    g1 <- add.edges(g1,t(e))
    g1 <- as.undirected(g1)
    #E(g1)$weight <- (edgeSize$score)*2
    #g2
    g2 <- add.vertices(g, nrow(v), 
	         name=as.character(v$Name), 
	         color=as.character(v$color_hub), 
	         shape=as.character(v$shape), 
	         size=as.character(v$size),
             label = as.character(v$label1) )
    g2 <- add.edges(g2,t(e))
    g2 <- as.undirected(g2)
    #E(g2)$weight <- (edgeSize$score)*2
    
    
    ### Step6: Plotting
    #Step6a : Other parameters
    set.seed(123)
    lay1 <- layout.fruchterman.reingold(g1)
    lay2 <- layout.fruchterman.reingold(g2)
    sizeCut <- round(seq(min(v$size), max(v$size), max(v$size)/4))
    comm_legend <- unique(v[,c("community","color_community")])
    comm_legend <- comm_legend[order(match(comm_legend$community, names(clustersList))),]
    rownames(comm_legend) <- NULL 
    #Step6b: Plot
    if(combine == "FALSE"){
        pdf(paste(names(GENE_LIST),"/", names(GENE_LIST), "__", int_type,"_int_c",confidence,".pdf",sep=""), width = 15, height = 15) 
    }else{
        pdf(paste(names(GENE_LIST),"/", names(GENE_LIST), "__", int_type,"_int_c",confidence,"_combinedhub.pdf",sep=""), width = 15, height = 15)
    }

    # Only hub gene labeled
    plot(g2, 
	    layout=lay1, 
        #edge.width = E(g2)$weight,
	    vertex.label=v$label2, 
        vertx.label.dist=0.2,
	    vertex.label.color="black", 
        vertex.frame.color = "white",
	    vertex.color=V(g2)$color, 
	    vertex.shape=V(g2)$shape, 
	    vertex.label.cex=1, 
	    vertex.size=as.numeric(V(g2)$size) ) 
    #legend('topright',legend=unique(sizeCut),pt.cex= sizeCut ,col='black')
    a <- legend('topright',legend=unique(sizeCut),pt.cex=sizeCut/200,col='white',
            pch=21, pt.bg='white', cex = 2)
    x <- (a$text$x + a$rect$left) / 2 
    y <- a$text$y
    symbols(x,y,circles=sizeCut/200,inches=FALSE,add=TRUE,bg='black')

    dev.off()
} 