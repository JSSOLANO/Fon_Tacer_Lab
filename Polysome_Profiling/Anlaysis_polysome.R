

BiocManager::install("org.Mm.eg.db")

library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(igraph)
library(ggplot2)
library(tibble)
library(ggrepel)
library(DESeq2)
library(biomaRt)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(ReactomePA)
library(msigdbr)

#c("pwt1inp","pwt2inp","pwt2inp","pko1inp","pko2inp","pko3inp")
#c("pwt1mon","pwt2mon","pwt3mon","pko1mon","pko2mon","pko3mon")
#c("pwt1lp","pwt2lp","pwt3lp","pko1lp","pko2lp","pko3lp")
#c("pwt1hp","pwt2hp","pwt3hp","pko1hp","pko2hp","pko3hp")

#DATA

counting_matriz <- read.csv(file="/Users/juansola/Documents/TC3R/research/polysome_profiling/240902_polysome_profiling/240522_polysome_profiling_input/Counting/Outputs/rawcounts.csv")

samples <- subset(counting_matriz, select = X25:X30)
head(samples)

#plot a PCA of the samples

colnames(samples) <- c("pwt1hp","pwt2hp","pwt3hp","pko1hp","pko2hp","pko3hp")


pca_samples <- samples
#samples$ensembl <- counting_matriz$X

pca_samples <- pca_samples[apply(pca_samples, 1, var) != 0, ]

#scale the data
data_log <- log1p(pca_samples)
data_scaled <- t(data_log) 


pca_result <- prcomp(data_scaled, scale. = TRUE)

pca_scores <- as.data.frame(pca_result$x)

pca_df <- pca_scores %>%
  mutate(sample = rownames(pca_scores)) %>%
  mutate(group = ifelse(grepl("wt", sample), "wt", "ko"))

# Plotting the PCA results using ggplot2

ggplot(pca_df, aes(x = PC1, y = PC2, label = sample, color = group)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3) +
  labs(title = "PCA of Samples",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme(plot.title = element_text(hjust=0.5))

ggsave("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240524_Pituitary_Heavy_polysomes/pheavy_polysome.pdf")




#DEA with deseq2


condition <- data.frame("samples" = c("pwt1hp","pwt2hp","pwt3hp","pko1hp","pko2hp","pko3hp"),
"condition" =c("control","control","control","treatment","treatment","treatment"))

samples_expression <- samples

rownames(samples_expression) <- counting_matriz$X


dea <- DESeqDataSetFromMatrix(countData =samples_expression , colData = condition,design = ~ condition)
dea$condition <- relevel(dea$condition, ref = "control")
dea <- DESeq(dea)
res <- results(dea)
resOrdered <- res[order(res$padj),]
resOrdered$gene <- row.names(resOrdered)
rownames(resOrdered) <- NULL
sublog <-resOrdered[!is.na(resOrdered$log2FoldChange) & abs(resOrdered$log2FoldChange) >= 0.5,]
resSig <- subset(sublog, pvalue < 0.05) #select genes with a p-value less than 0.05
resSig$gene <- sapply(strsplit(resSig$gene, "[.]"), `[`, 1) # delet decimal in ensmbl code
resOrdered$gene <- sapply(strsplit(resOrdered$gene, "[.]"), `[`, 1) # delet decimal in ensmbl code

write.csv(resSig, "/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240524_log005_and_pvalue005/240624_Pituitary_Heavy_polysomes/Deseq/DEG_pvalue005_log05.csv")
write.csv(resOrdered, "/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240524_log005_and_pvalue005/240624_Pituitary_Heavy_polysomes/Deseq/DEG_not_cutoff.csv")


############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################

#GO enrichment


not_na_pvalue <- resSig


positives <- not_na_pvalue[not_na_pvalue$log2FoldChange > 0 ,]
negatives <- not_na_pvalue[not_na_pvalue$log2FoldChange < 0 ,]


positives <- total_input[total_input$foldchange > 0 ,]
negatives <- total_input[total_input$foldchange < 0 ,]


ego_cutoff_genes <- enrichGO(gene = total_input$genes,
                   keyType = "ENSEMBL",
                   OrgDb = org.Mm.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)


ego_cutoff_genes_positives <- enrichGO(gene = positives$genes,
                   keyType = "ENSEMBL",
                   OrgDb = org.Mm.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)


ego_cutoff_genes_negatives <- enrichGO(gene = negatives$genes,
                   keyType = "ENSEMBL",
                   OrgDb = org.Mm.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)


go_results_positives <- data.frame(ego_cutoff_genes_positives)
go_results_negatives <- data.frame(ego_cutoff_genes_negatives)
go_results_allgenes <- data.frame(ego_cutoff_genes)

write.csv(go_results_positives, "/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240725_hypothalamus/240725_hypothalamus_input_polysomes/GO/upregulated_GO.csv")
write.csv(go_results_negatives, "/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240725_hypothalamus/240725_hypothalamus_input_polysomes/GO/downregulated_GO.csv")
write.csv(go_results_allgenes, "/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240725_hypothalamus/240725_hypothalamus_input_polysomes/GO/all_genes_GO.csv")


if (nrow(go_results_allgenes) < 1 ) { 
  print ("there are not GO categories")
} else {
  plot_go_up <-dotplot(ego_cutoff_genes, showCategory=20, font.size = 14)
  plot_go_p=paste("GO","PLOT_all_gens.png",sep="_")
  #png(filename=paste(GO_folder,plot_go_p,sep="/"), width=400, height=700)
  plot(plot_go_up)
  #dev.off()
}



ggsave("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240725_hypothalamus/240725_hypothalamus_input_polysomes/GO/all_genes_GO.pdf")



if (nrow(go_results_positives) < 1 ) { 
  print ("there are not GO categories")
} else {
  plot_go_up <-dotplot(ego_cutoff_genes_positives, showCategory=20, font.size = 14)
  plot_go_p=paste("GO","PLOT_positives.png",sep="_")
  #png(filename=paste(GO_folder,plot_go_p,sep="/"), width=400, height=700)
  plot(plot_go_up)
  #dev.off()
}
ggsave("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240725_hypothalamus/240725_hypothalamus_input_polysomes/GO/upregulated_GO.pdf")




if (nrow (go_results_negatives) < 1) {
  print ("there are not GO categories")
  
}else{
  plot_go_down <-dotplot(ego_cutoff_genes_negatives, showCategory=20, font.size = 14)
  plot_go_d=paste("GO","PLOT_negatives.png",sep="_")
  #png(filename=paste(GO_folder,plot_go_d,sep="/"), width=400, height=700)
  plot(plot_go_down)
  #dev.off() 
}

ggsave("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240725_hypothalamus/240725_hypothalamus_light_polysomes/GO/downregulated_GO.pdf")

############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################






#GSEA analysis
Gsea <- resOrdered[!is.na(resOrdered$log2FoldChange=="NA"),]


foldchanges <- Gsea$log2FoldChange
names(foldchanges) <- Gsea$gene

foldchanges <- sort(foldchanges, decreasing = TRUE)

gseaenrichment_result <- gseGO(geneList= gsea_input,
                               OrgDb = org.Mm.eg.db,
                               ont = "BP",
                               keyType = "ENSEMBL",
                               minGSSize    = 10,
                               maxGSSize    = 500,
                               pvalueCutoff = 0.05,
                               verbose      = FALSE)

require(DOSE)

gesaplot <- dotplot(gseaenrichment_result, showCategory=10, split=".sign",font.size = 8) + facet_grid(.~.sign)
plot(gesaplot)
ggsave("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240725_hypothalamus/240725_hypothalamus_input_polysomes/GSEA/gsea_h_imput_polysomes_polysomes.pdf")


gsea_final_result <- gseaenrichment_result@result
write.csv(gsea_final_result, "/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240725_hypothalamus/240725_hypothalamus_input_polysomes/GSEA/gsea_h_imput_polysomes_polysomes.csv")




############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################




#Network construction

#implemented database = https://www.inetbio.org/mousenet/downloadnetwork.php

database <- read.csv(file = "/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_input/MouseNetV2_entrez.txt",sep = '\t', header=FALSE)
colnames(database) <- c("gene_1","gene_2","score")

annotated_genes <- AnnotationDbi::select(org.Mm.eg.db, # database
                                            keys = total_input$genes,  # data to use for retrieval
                                            columns = c("ENSEMBL","GENENAME","ENTREZID","SYMBOL"), # information to retreive for given data
                                            keytype = "ENSEMBL")


entrez_numer <- as.list(annotated_genes[!is.na(annotated_genes$ENTREZID=="NA"),])

entrez_numer <- entrez_numer$ENTREZID
df_network <- data.frame()

for (i in 1:nrow(database)){
    current_row <- database[i, ]
    id_1 <- current_row[,1]
    id_2 <- current_row[,2]
    if ((id_1 %in% entrez_numer) & (id_2 %in% entrez_numer)){
        df_network <- rbind(df_network, current_row)
    }
}

#anotate the new dataframe with the symbol of each gene
symbol_1 <- AnnotationDbi::select(org.Mm.eg.db, # database
                                         keys = as.character(df_network$gene_1),  # data to use for retrieval
                                         columns = c("SYMBOL"), # information to retreive for given data
                                         keytype = "ENTREZID")
symbol_2 <- AnnotationDbi::select(org.Mm.eg.db, # database
                                         keys = as.character(df_network$gene_2),  # data to use for retrieval
                                         columns = c("SYMBOL"), # information to retreive for given data
                                         keytype = "ENTREZID")


#Select symbosl infromation and generate the complete network
df_network["symbol_1"] <- symbol_1$SYMBOL
df_network["symbol_2"] <- symbol_2$SYMBOL


#generate the complete network
graph_node <- df_network[c("symbol_1","symbol_2")]
weights <- c(df_network$score)

g <- graph_from_data_frame(graph_node, directed=FALSE)
E(g)$weight <- weights # Use directed=TRUE for directed graphs

# Use Fruchterman-Reingold layout with adjusted parameters
layout_fr <- layout_with_fr(g, niter = 500, area = vcount(g)^2)

pdf("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240725_hypothalamus/240725_hypothalamus_input_polysomes/Network/hypothalamus_input_polysomes_Network.pdf", width = 15, height = 20)

plot_net <- plot(g, layout = layout_fr, vertex.size = 8, vertex.label.cex = 0.7)

title(main = "Hypothalamus input Polysomes Network", cex.main = 4)

dev.off()



############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################



monosomes_cutoff <- read.csv("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240524_log005_and_pvalue005/240624_Pituitary_monosomes_polysomes/Deseq/DEG_pvalue005_log05.csv")


monosomes_not_cutoff<- read.csv("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240524_log005_and_pvalue005/240624_Pituitary_monosomes_polysomes/Deseq/DEG_not_cutoff.csv")


annotated_cutoof <- AnnotationDbi::select(org.Mm.eg.db, # database
                                            keys = monosomes_cutoff$gene,  # data to use for retrieval
                                            columns = c("ENSEMBL","SYMBOL"), # information to retreive for given data
                                            keytype = "ENSEMBL")

annotated_not_cutoof <- AnnotationDbi::select(org.Mm.eg.db, # database
                                            keys = monosomes_not_cutoff$gene,  # data to use for retrieval
                                            columns = c("ENSEMBL","SYMBOL"), # information to retreive for given data
                                            keytype = "ENSEMBL")



annotated_cutoof <- annotated_cutoof[!is.na(annotated_cutoof$SYMBOL=="NA"),]

annotated_not_cutoof <- annotated_not_cutoof[!is.na(annotated_not_cutoof$SYMBOL=="NA"),]





DEG_pvalue005_log05 <- merge(annotated_cutoof, monosomes_cutoff, by.x = "ENSEMBL", by.y = "gene", all = FALSE)
DEG_not_cutoff<- merge(annotated_not_cutoof, monosomes_not_cutoff, by.x = "ENSEMBL", by.y = "gene", all = FALSE)



write.csv(DEG_not_cutoff, "/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240524_log005_and_pvalue005/240624_Pituitary_monosomes_polysomes/Deseq/Annotated_DEG_not_cutoff.csv")

write.csv(DEG_pvalue005_log05, "/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240524_log005_and_pvalue005/240624_Pituitary_monosomes_polysomes/Deseq/Annotated_DEG_pvalue005_log05.csv")



############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################


ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")


gene_info <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'transcript_length'),
                   mart = ensembl)


# Calculate gene lengths in kilobases (kb)

#how to calcualte TPM: https://bioinformatics.ccr.cancer.gov/btep/questions/what-is-the-difference-between-rpkm-fpkm-and-tpm

gene_info$gene_length_kb <- gene_info$transcript_length / 1000


# Get the average transcript length for each gene
average_transcripts <- gene_info %>%
  group_by(ensembl_gene_id) %>%
  summarise(average_gene_length_kb = mean(gene_length_kb))

#inner merge the polysome dataframe with the average length dataframe

merged_data <- inner_join(counting_matriz, average_transcripts,   by = c("X" = "ensembl_gene_id"))

#NORMILZE data to TPM

heatmap_data <- subset(merged_data, select = X1:X32)
divided_by_length <- heatmap_data/ merged_data$average_gene_length_kb

scaling_factor <-  colSums(divided_by_length)

divided_by_scalin_factor <- sweep(divided_by_length, 2, scaling_factor, FUN = "/")

TPM_results <- divided_by_scalin_factor * 10^6

TPM_results$genes_id <- merged_data$X


annotated_for_tpms <- AnnotationDbi::select(org.Mm.eg.db, # database
                                            keys = TPM_results$genes_id,  # data to use for retrieval
                                            columns = c("ENSEMBL","SYMBOL","GENENAME"), # information to retreive for given data
                                            keytype = "ENSEMBL")


annotated_for_tpms <- annotated_for_tpms[!duplicated(annotated_for_tpms$ENSEMBL), ]

TPM_results <- merge(TPM_results, annotated_for_tpms, by.x = "genes_id", by.y = "ENSEMBL", all.x = FALSE)

write.csv(TPM_results,"/Users/juansola/Desktop/tpm_values.csv",row.names = FALSE)
TPM_results$genes_id
selection_vector <- c("ENSMUSG00000033585","ENSMUSG00000045232","ENSMUSG00000059824","ENSMUSG00000056972","ENSMUSG00000100826","ENSMUSG00000005268")

filter_tpms <- TPM_results[TPM_results$genes_id %in% selection_vector, ]
rownames(TPM_results)


###################################################
###################################################
###################################################
###################################################
#violin plot ######################################


ensemblgeneid <- "ENSMUSG00000100826"
#fractions 
h_WT_Pituitary <-   melt(filter_tpms[filter_tpms$genes_id %in% ensemblgeneid,c("X25","X26","X27")], variable.name = "Replicate", value.name = "Value")
h_KO_Pituitary <-  melt(filter_tpms[filter_tpms$genes_id %in% ensemblgeneid,c("X28","X29","X30")], variable.name = "Replicate", value.name = "Value") 

l_WT_Pituitary <- melt(filter_tpms[filter_tpms$genes_id %in% ensemblgeneid,c("X17","X18","X19")], variable.name = "Replicate", value.name = "Value") 
l_KO_Pituitary <- melt(filter_tpms[filter_tpms$genes_id %in% ensemblgeneid,c("X20","X21","X22")], variable.name = "Replicate", value.name = "Value") 

m_WT_Pituitary <- melt(filter_tpms[filter_tpms$genes_id %in% ensemblgeneid,c("X9","X10","X11")], variable.name = "Replicate", value.name = "Value")
m_KO_Pituitary <- melt(filter_tpms[filter_tpms$genes_id %in% ensemblgeneid,c("X12","X13","X14")], variable.name = "Replicate", value.name = "Value") 

i_WT_Pituitary <- melt(filter_tpms[filter_tpms$genes_id %in% ensemblgeneid,c("X1","X2","X3")], variable.name = "Replicate", value.name = "Value")
i_KO_Pituitary <- melt(filter_tpms[filter_tpms$genes_id %in% ensemblgeneid,c("X4","X5","X6")], variable.name = "Replicate", value.name = "Value") 


h_WT_Pituitary$Fraction <- "H_WT_Pituitary"
h_KO_Pituitary$Fraction <- "H_KO_Pituitary"
l_WT_Pituitary$Fraction <- "L_WT_Pituitary"
l_KO_Pituitary$Fraction <- "L_KO_Pituitary"
m_WT_Pituitary$Fraction <- "M_WT_Pituitary"
m_KO_Pituitary$Fraction <- "M_KO_Pituitary"
i_WT_Pituitary$Fraction <- "I_WT_Pituitary"
i_KO_Pituitary$Fraction <- "I_KO_Pituitary"


combined_df <- rbind(h_WT_Pituitary, h_KO_Pituitary, l_WT_Pituitary, l_KO_Pituitary, m_WT_Pituitary, m_KO_Pituitary,i_WT_Pituitary,i_KO_Pituitary)


violin_plot <- ggplot(combined_df, aes(x = Fraction, y = Value, fill = Fraction)) +
  geom_violin(trim = FALSE) +
  labs(title = "TPM Prlr Expression", x = "Fraction", y = "TPM") +
  theme(plot.title = element_text(hjust=0.5))+
  theme(plot.title = element_text(size = 20))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave("/Users/juansola/Desktop/Pituitary_Prlr.png", plot = violin_plot, width = 8, height = 6, dpi = 300)












selection_vector <- c("ENSMUSG00000033585","ENSMUSG00000045232","ENSMUSG00000059824","ENSMUSG00000056972","ENSMUSG00000100826","ENSMUSG00000005268")

ensemblgeneid <- "ENSMUSG00000005268"

h_WT_hypothalamus <-   filter_tpms[filter_tpms$genes_id %in% ensemblgeneid,c("X31","X32","X23","X24","X15","X16","X7","X8")]

colnames(h_WT_hypothalamus) <- c('H_WT_hypothalamus','H_KO_hypothalamus','L_WT_hypothalamus','L_KO_hypothalamus','M_WT_hypothalamus','M_KO_hypothalamus','I_WT_hypothalamus','I_KO_hypothalamus')

row_df <- data.frame(
  Category = colnames(h_WT_hypothalamus),  # Column names as categories
  Values = as.numeric(h_WT_hypothalamus)   # Values in the row as numeric
)

violin_plot <- ggplot(row_df, aes(x = Category, y = Values, fill = Category)) +
  geom_bar(stat = "identity") +
  labs(title = "TPM Prlr Expression", x = "Fraction", y = "TPM") +
  theme(plot.title = element_text(hjust=0.5))+
  theme(plot.title = element_text(size = 20))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave("/Users/juansola/Desktop/hypothalamus_Prlr.png", plot = violin_plot, width = 8, height = 6, dpi = 300)


#create hierarchical clustering

log2_tpm_data <- as.matrix(log2(TPM_results + 1)) # Adding 1 to avoid log(0)

constant_rows <- apply(log2_tpm_data, 1, function(x) sd(x) == 0)

#Select a subset to create the clustering
subset <- log2_tpm_data[!constant_rows, ]

subset <- subset[1:20000,]

# Step 2: Calculate the distance matrix using Euclidean distance
treatment <- c("P.WT1.Inp","P.WT2.Inp","P.WT3.Inp","P.KO1.Inp","P.KO2.Inp","P.KO3.Inp","H.WT3.Inp", "H.KO3.Inp","P.WT1.Mon","P.WT2.Mon","P.WT3.Mon",
"P.KO1.Mon","P.KO2.Mon","P.KO3.Mon","H.WT3.Mon","H.KO3.Mon","P.WT1.LP","P.WT2.LP","P.WT3.LP","P.KO1.LP","P.KO2.LP","P.KO3.LP","H.WT3.LP","H.KO3.LP",
"P.WT1.HP","P.WT2.HP","P.WT3.HP","P.KO1.HP","P.KO2.HP","P.KO3.HP","H.WT3.HP","H.KO3.HP")

colnames(subset) <- treatment


annotation_col <- data.frame(
    Group = treatment
)
rownames(annotation_col) <- colnames(subset)

group_colors <- colorRampPalette(brewer.pal(9, "Set1"))(32)

annotation_colors <- list(
    Group = setNames(group_colors, unique(treatment))
)


custom_legend_labels <- list(title = "FoldChange")
pdf("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/heatmap_polysomes_comparisons.pdf", width = 15, height = 20)

map <- pheatmap(
    subset,
    annotation_col = annotation_col,
    annotation_colors = annotation_colors,
    clustering_distance_rows = "euclidean",
    clustering_method = "complete",
    legend_labels = custom_legend_labels,
    color = colorRampPalette(c("#B10909", "white", "#2c772c"))(60),
    breaks = seq(-2, 2, length.out = 60),
    scale = "row",
    show_rownames = FALSE,
    main = "Sample Clustering"
)

dev.off()



############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################



#analysis of Hypothalamus

#extrac information of each group, then, rest the TPMs of KO samples by the TPMs of the wildtype

input_Hypothalamus <- TPM_results[c("X7","X8")]
input_Hypothalamus <-input_Hypothalamus+1
input_Hypothalamus$"genes" <- merged_data$X

input_Hypothalamus$"foldchange" <-  log2(input_Hypothalamus$X8/input_Hypothalamus$X7)

input_Hypothalamus <- input_Hypothalamus[order(-input_Hypothalamus$foldchange), ]


total_input <- input_Hypothalamus[abs(input_Hypothalamus$foldchange) >= 0.5, ]




annotated_input <- AnnotationDbi::select(org.Mm.eg.db, # database
                                            keys = total_input$genes,  # data to use for retrieval
                                            columns = c("ENSEMBL","GENENAME","ENTREZID","SYMBOL"), # information to retreive for given data
                                            keytype = "ENSEMBL")

annotated_input_genes <- merge(total_input, annotated_input, by.x = "genes", by.y = "ENSEMBL", all = FALSE)


write.csv(input_Hypothalamus, "/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240725_hypothalamus/240725_hypothalamus_input_polysomes/foldchange/fold_change_input.csv")


write.csv(annotated_input_genes, "/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240725_hypothalamus/240725_hypothalamus_input_polysomes/foldchange/annotated_input.csv")


gsea_input <- total_input$foldchange
names(gsea_input) <- total_input$genes

gsea_input <- sort(gsea_input, decreasing = TRUE)





Monosomes_Hypothalamus <- TPM_results[c("X15","X16")]

Monosomes_Hypothalamus <-Monosomes_Hypothalamus+1
Monosomes_Hypothalamus$"genes" <- merged_data$X

Monosomes_Hypothalamus$"foldchange" <-  log2(Monosomes_Hypothalamus$X16/Monosomes_Hypothalamus$X15)

Monosomes_Hypothalamus <- Monosomes_Hypothalamus[order(-Monosomes_Hypothalamus$foldchange), ]



total_monosomes <- Monosomes_Hypothalamus[abs(Monosomes_Hypothalamus$foldchange) >= 0.5, ]

annotated_monosomes <- AnnotationDbi::select(org.Mm.eg.db, # database
                                            keys = total_monosomes$genes,  # data to use for retrieval
                                            columns = c("ENSEMBL","GENENAME","ENTREZID","SYMBOL"), # information to retreive for given data
                                            keytype = "ENSEMBL")

annotated_monosomes_genes <- merge(total_monosomes, annotated_monosomes, by.x = "genes", by.y = "ENSEMBL", all = FALSE)

write.csv(annotated_monosomes_genes, "/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240725_hypothalamus/240725_hypothalamus_monosomes/foldchange/annotated_monosomes.csv")

write.csv(total_monosomes, "/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240725_hypothalamus/240725_hypothalamus_monosomes/foldchange/fold_change_monosomes.csv")



gsea_monosomes <- total_monosomes$foldchange
names(gsea_monosomes) <- total_monosomes$genes

gsea_monosomes <- sort(gsea_monosomes, decreasing = TRUE)





light_polysomes_Hypothalamus <- TPM_results[c("X23","X24")]

light_polysomes_Hypothalamus <-light_polysomes_Hypothalamus+1
light_polysomes_Hypothalamus$"genes" <- merged_data$X

light_polysomes_Hypothalamus$"foldchange" <-  log2(light_polysomes_Hypothalamus$X24/light_polysomes_Hypothalamus$X23)

light_polysomes_Hypothalamus <- light_polysomes_Hypothalamus[order(-light_polysomes_Hypothalamus$foldchange), ]



total_light_polysomes <- light_polysomes_Hypothalamus[abs(light_polysomes_Hypothalamus$foldchange) >= 0.5, ]


annotated_light_polysomes <- AnnotationDbi::select(org.Mm.eg.db, # database
                                            keys = total_light_polysomes$genes,  # data to use for retrieval
                                            columns = c("ENSEMBL","GENENAME","ENTREZID","SYMBOL"), # information to retreive for given data
                                            keytype = "ENSEMBL")

annotated_light_polysomes_genes <- merge(total_light_polysomes, annotated_light_polysomes, by.x = "genes", by.y = "ENSEMBL", all = FALSE)

write.csv(annotated_light_polysomes_genes, "/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240725_hypothalamus/240725_hypothalamus_light_polysomes/foldchange/annotated_light_polysomes.csv")

write.csv(total_light_polysomes, "/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240725_hypothalamus/240725_hypothalamus_light_polysomes/foldchange/fold_change_light_polysomes.csv")


gsea_light_polysomes <- total_light_polysomes$foldchange
names(gsea_light_polysomes) <- total_light_polysomes$genes

gsea_light_polysomes <- sort(gsea_light_polysomes, decreasing = TRUE)







heavy_polysomes_Hypothalamus <- TPM_results[c("X31","X32")]

heavy_polysomes_Hypothalamus <-heavy_polysomes_Hypothalamus+1
heavy_polysomes_Hypothalamus$"genes" <- merged_data$X

heavy_polysomes_Hypothalamus$"foldchange" <-  log2(heavy_polysomes_Hypothalamus$X32/heavy_polysomes_Hypothalamus$X31)

heavy_polysomes_Hypothalamus <- heavy_polysomes_Hypothalamus[order(-heavy_polysomes_Hypothalamus$foldchange), ]


total_heavy_polysomes <- heavy_polysomes_Hypothalamus[abs(heavy_polysomes_Hypothalamus$foldchange) >= 0.5, ]



annotated_heavy_polysomes <- AnnotationDbi::select(org.Mm.eg.db, # database
                                            keys = total_heavy_polysomes$genes,  # data to use for retrieval
                                            columns = c("ENSEMBL","GENENAME","ENTREZID","SYMBOL"), # information to retreive for given data
                                            keytype = "ENSEMBL")

annotated_heavy_polysomes_genes <- merge(total_heavy_polysomes, annotated_heavy_polysomes, by.x = "genes", by.y = "ENSEMBL", all = FALSE)


write.csv(heavy_polysomes_Hypothalamus, "/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240725_hypothalamus/240725_hypothalamus_heavy_polysomes/foldchange/fold_change_heavy_polysomes.csv")

write.csv(annotated_heavy_polysomes_genes, "/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240725_hypothalamus/240725_hypothalamus_heavy_polysomes/foldchange/annotated_heavy_polysomes.csv")



gsea_heavy_polysomes <- total_heavy_polysomes$foldchange
names(gsea_heavy_polysomes) <- total_heavy_polysomes$genes

gsea_heavy_polysomes <- sort(gsea_heavy_polysomes, decreasing = TRUE)

############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################



#dot plot figure

#dot plot figure hypothalamus

hypo_heavy_polysome_log2 <- read.csv("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240725_hypothalamus/240725_hypothalamus_heavy_polysomes/foldchange/fold_change_heavy_polysomes.csv")

hypo_heavy_log2 <- hypo_heavy_polysome_log2[,c("genes","foldchange")]
colnames(hypo_heavy_log2)[2] <- "heavy_polysome"

filter_heavy <- hypo_heavy_log2[abs(hypo_heavy_log2$heavy_polysome) >= 0.5, ]
nrow(filter_heavy)


hypo_input_polysome_log2 <- read.csv("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240725_hypothalamus/240725_hypothalamus_input_polysomes/foldchange/fold_change_input.csv")

hypo_input_log2 <- hypo_input_polysome_log2[,c("genes","foldchange")]
colnames(hypo_input_log2)[2] <- "input_polysome"

filter_input <- hypo_input_log2[abs(hypo_input_log2$input_polysome) >= 0.5, ]







#share and unique genes 
share_genes <- merge(filter_heavy, filter_input, by="genes", all = FALSE)
only_heavy <- setdiff(filter_heavy$genes,filter_input$genes)
only_input <- setdiff(filter_input$genes, filter_heavy$genes)

#dataframe to create the graph
graph_hypo <- merge(hypo_heavy_log2, hypo_input_log2, by="genes", all = FALSE)

#assign a different color to the shared genes

graph_hypo$color <- ifelse(graph_hypo$genes %in% share_genes$genes, "#800e0e", ifelse(graph_hypo$genes %in% only_heavy, "#055a05", ifelse(graph_hypo$genes %in% only_input, "#28285d", "#d0cbcb")))


graph_hypo$label <- ifelse(graph_hypo$genes %in% share_genes$genes, "Shared (n 391)", ifelse(graph_hypo$genes %in% only_heavy, "Heavy Polysome (n 1167)", ifelse(graph_hypo$genes %in% only_input, "Input Polysome (n 1356)", "Others")))

plot_hypo <- ggplot(graph_hypo, aes(x = input_polysome, y = heavy_polysome,color = label)) +
  geom_point(size = 2) +  # Set point size for better visibility
  scale_color_manual(values = c("Shared (n 391)" = "#800e0e", "Heavy Polysome (n 1167)" = "#055a05", "Input Polysome (n 1356)" = "#28285d", "Others" = "#d0cbcb")) +
  geom_point() +
  labs(title = "Hypothalamus Expression", x = "Log2 Fold Change Input Polysome", y = "Log2 Fold Change Heavy Polysome") +
  theme_classic()+
  theme(text=element_text(size=10),plot.title = element_text(size = 20,hjust = 0.8))



ggsave("/Users/juansola/Desktop/Hypothalamus Expression.tif", plot = plot_hypo, width = 8, height = 10, dpi = 300)







#dot plot figure pituitary

pitu_heavy_polysome <- read.csv("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240524_pituitary/240524_log005_and_pvalue005/240624_Pituitary_Heavy_polysomes/Deseq/DEG_not_cutoff.csv")

pitu_heavy_polysome <- pitu_heavy_polysome[,c("gene","log2FoldChange","pvalue")]

filtered_pitu_heavy_polysome <- pitu_heavy_polysome %>% filter(abs(log2FoldChange) >= 0.5 & pvalue < 0.05)

pitu_heavy_log2 <- pitu_heavy_polysome[,c("gene","log2FoldChange")]

colnames(pitu_heavy_log2)[2] <- "heavy_polysome"



pitu_input_polysome <- read.csv("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240524_pituitary/240524_log005_and_pvalue005/240624_Pituitary_input_polysomes/Deseq/DEG_not_cutoff.csv")

pitu_input_polysome <- pitu_input_polysome[,c("gene","log2FoldChange","pvalue")]


filtered_pitu_input_polysome <- pitu_input_polysome %>% filter(abs(log2FoldChange) >= 0.5 & pvalue < 0.05)

pitu_input_log2 <- pitu_input_polysome[,c("gene","log2FoldChange")]


colnames(pitu_input_log2)[2] <- "input_polysome"

#share and unique genes 
share_genes <- merge(filtered_pitu_heavy_polysome, filtered_pitu_input_polysome, by="gene", all = FALSE)
only_heavy_pituitary <- setdiff(filtered_pitu_heavy_polysome$gene,filtered_pitu_input_polysome$gene)
only_input_pituitary <- setdiff(filtered_pitu_input_polysome$gene, filtered_pitu_heavy_polysome$gene)

nrow(share_genes)
length(only_heavy_pituitary)
length(only_input_pituitary)
nrow(filtered_pitu_heavy_polysome)

#dataframe to create the graph
pitu_heavy_log2 <- na.omit(pitu_heavy_log2)
pitu_input_log2 <- na.omit(pitu_input_log2)

graph_hypo <- merge(pitu_heavy_log2, pitu_input_log2, by="gene", all = FALSE)

#assign a different color to the shared genes

graph_hypo$color <- ifelse(graph_hypo$gene %in% share_genes$gene, "#800e0e", ifelse(graph_hypo$gene %in% only_heavy_pituitary, "#055a05", ifelse(graph_hypo$gene %in% only_input_pituitary, "#28285d", "#f8f8f6")))


graph_hypo$label <- ifelse(graph_hypo$gene %in% share_genes$gene, "Shared (n 42)", ifelse(graph_hypo$gene %in% only_heavy_pituitary, "Heavy Polysome (n 394)", ifelse(graph_hypo$gene %in% only_input_pituitary, "Input Polysome (n 240)", "Others")))

graph_hypo <- graph_hypo[graph_hypo$label != "Others", ]


plot_pitu <- ggplot(graph_hypo, aes(x =input_polysome , y = heavy_polysome,color = label)) +
  geom_point(size = 2) +  # Set point size for better visibility
  scale_color_manual(values = c("Shared (n 42)" = "#800e0e", "Heavy Polysome (n 394)" = "#055a05", "Input Polysome (n 240)" = "#28285d", "Others" = "#f8f6f686")) +
  geom_point() +
  labs(title = "Pituitary Expression", x = "Log2 Fold Change Input Polysome", y = "Log2 Fold Change Heavy Polysome") +
  theme_classic()+
  theme(text=element_text(size=10),plot.title = element_text(size = 20,hjust = 0.5))

plot(plot_pitu)

ggsave("/Users/juansola/Desktop/Pituitary Expression.tif", plot = plot_pitu, width = 8, height = 10, dpi = 300)





############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################
############################################################




#neirhcment analysisi of translatome, the DEG in Heavy polysome but not in Input fraction


go_enrichment <-  function(datafarame,path){
  positives <- datafarame[datafarame$log2FoldChange > 0 ,]
  negatives <- datafarame[datafarame$log2FoldChange < 0 ,]
  
  ego_cutoff_genes <- enrichGO(gene = datafarame$ENSEMBL,
                   keyType = "ENSEMBL",
                   OrgDb = org.Mm.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)

  ego_cutoff_genes_positives <- enrichGO(gene = positives$ENSEMBL,
                   keyType = "ENSEMBL",
                   OrgDb = org.Mm.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)


  ego_cutoff_genes_negatives <- enrichGO(gene = negatives$ENSEMBL,
                   keyType = "ENSEMBL",
                   OrgDb = org.Mm.eg.db, 
                   ont = "BP", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)

  
  go_results_allgenes <- data.frame(ego_cutoff_genes)
  go_results_positives <- data.frame(ego_cutoff_genes_positives)
  go_results_negatives <- data.frame(ego_cutoff_genes_negatives)
  

  output_allgenes_file <- paste(path,"GO/GO_all-genes.csv",sep = "")
  output_go_positives_file <- paste(path,"GO/GO_positives.csv",sep = "")
  output_go_negatives_file <- paste(path,"GO/GO_genatives.csv",sep = "")

  
  write.csv(go_results_allgenes,output_allgenes_file)
  write.csv(go_results_positives,output_go_positives_file)
  write.csv(go_results_negatives,output_go_negatives_file)
  

  if (nrow(go_results_allgenes) < 1 ) { 
  print ("there are not GO categories")
} else {
  plot_go_all <-dotplot(ego_cutoff_genes, showCategory=20, font.size = 14)
  output_allgenes_plot <- paste(path,"GO/all_genes_GO.pdf",sep = "")
  ggsave(output_allgenes_plot, plot=plot_go_all)
}

  if (nrow(go_results_positives) < 1 ) { 
  print ("there are not GO categories")
} else {
  plot_go_up <-dotplot(ego_cutoff_genes_positives, showCategory=20, font.size = 14)
  output_up_plot <- paste(path,"GO/positives_GO.pdf",sep = "")
  ggsave(output_up_plot, plot=plot_go_up)


  
}

  if (nrow(go_results_negatives) < 1 ) { 
  print ("there are not GO categories")
} else {
  plot_down_up <-dotplot(ego_cutoff_genes_negatives, showCategory=20, font.size = 14)
  output_down_plot <- paste(path, "GO/negatives_GO.pdf",sep = "")
  ggsave(output_down_plot, plot=plot_down_up)


  
}




}


###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################

#translatome enrihcmene in Pituitary 

Pituitary_translatome <- pitu_heavy_polysome [pitu_heavy_polysome$ENSEMBL %in% only_heavy_pituitary,]


Pituitary_translatome_symbol  <- AnnotationDbi::select(org.Mm.eg.db, # database
                                            keys = Pituitary_translatome$ENSEMBL,  # data to use for retrieval
                                            columns = c("ENSEMBL","GENENAME","ENTREZID","SYMBOL"), # information to retreive for given data
                                            keytype = "ENSEMBL")


Pituitary_translatome_annotated <- merge(Pituitary_translatome, Pituitary_translatome_symbol, by = "ENSEMBL")





write.csv(Pituitary_translatome_annotated,"/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240812_translatome/Pituitary/Pituitary_translatome.csv",row.names=FALSE)

######################################################

go_enrichment(filtered_pitu_heavy_polysome,"/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240812_translatome/prueba/")

######################################################


gsea_Pituitary_translatome <- Pituitary_translatome[!is.na(Pituitary_translatome$log2FoldChange=="NA"),]


PT_foldchanges <- gsea_Pituitary_translatome$log2FoldChange
names(PT_foldchanges) <- gsea_Pituitary_translatome$ENSEMBL

PT_foldchanges <- sort(PT_foldchanges, decreasing = TRUE)


######################################################
######################################################
######################################################
######################################################
######################################################
######################################################


#translatome enrihcmene in Hypothalamus 
Hypothalamus_translatome <- hypo_heavy_log2 [hypo_heavy_log2$genes %in% only_heavy,]


colnames(Hypothalamus_translatome) <- c('ENSEMBL','log2FoldChange')

Hypothalamus_translatome_symbol  <- AnnotationDbi::select(org.Mm.eg.db, # database
                                            keys = Hypothalamus_translatome$ENSEMBL,  # data to use for retrieval
                                            columns = c("ENSEMBL","GENENAME","ENTREZID","SYMBOL"), # information to retreive for given data
                                            keytype = "ENSEMBL")


Hypothalamus_translatome_annotated <- merge(Hypothalamus_translatome, Hypothalamus_translatome_symbol, by = "ENSEMBL")

write.csv(Hypothalamus_translatome_annotated,"/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240812_translatome/Hypothalamus/Hypothalamus_translatome.csv",row.names=FALSE)



#################################################################

go_enrichment(Hypothalamus_translatome,"/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240812_translatome/Hypothalamus/")



#################################################################


gsea_Hypothalamus_translatome <- Hypothalamus_translatome[!is.na(Hypothalamus_translatome$log2FoldChange=="NA"),]


HT_foldchanges <- gsea_Hypothalamus_translatome$log2FoldChange
names(HT_foldchanges) <- gsea_Hypothalamus_translatome$ENSEMBL

HT_foldchanges <- sort(HT_foldchanges, decreasing = TRUE)


###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################


#GSEA

gseaenrichment_result <- gseGO(geneList= HT_foldchanges,
                               OrgDb = org.Mm.eg.db,
                               ont = "BP",
                               keyType = "ENSEMBL",
                               minGSSize    = 10,
                               maxGSSize    = 500,
                               pvalueCutoff = 0.05,
                               verbose      = FALSE)

require(DOSE)

gesaplot <- dotplot(gseaenrichment_result, showCategory=10, split=".sign",font.size = 8) + facet_grid(.~.sign)
plot(gesaplot)
ggsave("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240725_hypothalamus/240725_hypothalamus_input_polysomes/GSEA/gsea_h_imput_polysomes_polysomes.pdf")


gsea_final_result <- gseaenrichment_result@result
write.csv(gsea_final_result, "/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240725_hypothalamus/240725_hypothalamus_input_polysomes/GSEA/gsea_h_imput_polysomes_polysomes.csv")

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################


#Reactome
enrich_reactome <- enrichPathway(gene=entrez_numer,organism = "mouse",minGSSize = 10, maxGSSize = 500)

entrez_ids <- unique(unlist(strsplit(enrich_reactome$geneID, "/")))

gene_mapping <- bitr(entrez_ids, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Mm.eg.db)


enrich_reactome$geneID <- sapply(enrich_reactome$geneID, function(ids) {
  symbols <- gene_mapping$SYMBOL[match(unlist(strsplit(ids, "/")), gene_mapping$ENTREZID)]
  paste(symbols, collapse = "/")
})

barplot(enrich_reactome, showCategory=20)

cnet <- cnetplot(enrich_reactome, foldChange=enrich_reactome)
cnet <- cnet + theme_classic()

ggsave("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240812_translatome/Hypothalamus/reactome/Hypothalamus_nodes.png", plot = cnet, width = 8, height = 6, dpi = 300)




###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################



#clusterprofile analysis using whole MSigDB database




pitu_heavy <- read.csv("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240524_pituitary/240624_log005_and_pvalue005/240624_Pituitary_Heavy_polysomes/Deseq/DEG_pvalue005_log05.csv")
pitu_input <- read.csv("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240524_pituitary/240624_log005_and_pvalue005/240624_Pituitary_input_polysomes/Deseq/DEG_pvalue005_log05.csv")
pitu_light <- read.csv("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240524_pituitary/240624_log005_and_pvalue005/240624_Pituitary_light_polysomes/Deseq/DEG_pvalue005_log05.csv")
pitu_mono <- read.csv("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240524_pituitary/240624_log005_and_pvalue005/240624_Pituitary_monosomes_polysomes/Deseq/DEG_pvalue005_log05.csv")



list_piyu_deseq = list(heavy=pitu_heavy,input=pitu_input,light=pitu_light,mono=pitu_mono)
  

for (position in seq_along(list_piyu_deseq)){
  file=list_piyu_deseq[[position]]
  file_name= names(list_piyu_deseq)[position]

  enrich_result <- enricher(gene = file$gene,
                          TERM2GENE = all_gene_sets_msigdb[, c("gs_name", "ensembl_gene")])

  enrich_dataframe <- as.data.frame(enrich_result)

  csv_file_name= paste(file_name,".csv",sep="")
  write.csv(enrich_dataframe,paste("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240828_GSEA_msigdb_MH_pitu",csv_file_name,sep="/"))


  if (nrow(enrich_dataframe) < 1 ) { 
  print ("there are not GO categories")
} else {
  plot_msigdb  <- barplot(enrich_result, showCategory = 20)
  #Create name to save the file
  plot_go_p=paste(file_name,".png",sep="")
  #opne the png file
  png(filename=paste("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240828_GSEA_msigdb_MH_pitu",plot_go_p,sep="/"), width=400, height=700)
  plot(plot_msigdb)
  dev.off()
}
}


hypo_heavy <- read.csv("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240808_hypothalamus/240808_hypothalamus_heavy_polysomes/foldchange/fold_change_heavy_polysomes.csv")
hypo_input <- read.csv("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240808_hypothalamus/240808_hypothalamus_input_polysomes/foldchange/fold_change_input.csv")
hypo_light <- read.csv("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240808_hypothalamus/240808_hypothalamus_light_polysomes/foldchange/fold_change_light_polysomes.csv")
hypo_mono <- read.csv("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240808_hypothalamus/240808_hypothalamus_monosomes/foldchange/fold_change_monosomes.csv")


list_hypo_deseq = list(heavy=hypo_heavy,input=hypo_input,light=hypo_light,mono=hypo_mono)


for (position in seq_along(list_hypo_deseq)){

  file=list_hypo_deseq[[position]]
  filter_file <- file[abs(file$foldchange) >= 0.5, ]

  file_name= names(list_hypo_deseq)[position]

  enrich_result <- enricher(gene = filter_file$genes,
                          TERM2GENE = all_gene_sets_msigdb[, c("gs_name", "ensembl_gene")])

  enrich_dataframe <- as.data.frame(enrich_result)

  csv_file_name= paste(file_name,".csv",sep="")
  write.csv(enrich_dataframe,paste("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240827_enrichment_msigdn_hypo",csv_file_name,sep="/"))


  if (nrow(enrich_dataframe) < 1 ) { 
  print ("there are not GO categories")
} else {
  plot_msigdb  <- barplot(enrich_result, showCategory = 20)
  #Create name to save the file
  plot_go_p=paste(file_name,".png",sep="")
  #opne the png file
  png(filename=paste("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240827_enrichment_msigdn_hypo",plot_go_p,sep="/"), width=400, height=700)
  plot(plot_msigdb)
  dev.off()
}
}




##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################


#GSEA analysis using whole MSigDB database

collection <- c('H','C1','C2','C3','C4','C5','C6','C7','C8')

pitu_heavy <- read.csv("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240524_pituitary/240624_log005_and_pvalue005/240624_Pituitary_Heavy_polysomes/Deseq/DEG_not_cutoff.csv")
pitu_input <- read.csv("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240524_pituitary/240624_log005_and_pvalue005/240624_Pituitary_input_polysomes/Deseq/DEG_not_cutoff.csv")
pitu_light <- read.csv("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240524_pituitary/240624_log005_and_pvalue005/240624_Pituitary_light_polysomes/Deseq/DEG_not_cutoff.csv")
pitu_mono <- read.csv("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240524_pituitary/240624_log005_and_pvalue005/240624_Pituitary_monosomes_polysomes/Deseq/DEG_not_cutoff.csv")

list_piyu_deseq = list(heavy=pitu_heavy,input=pitu_input,light=pitu_light,mono=pitu_mono)



for (position in seq_along(list_piyu_deseq)){
  file=list_piyu_deseq[[position]]

  geneList <- file$log2FoldChange
  names(geneList) <- file$gene
  geneList <- sort(geneList, decreasing = TRUE)
  file_name= names(list_piyu_deseq)[position]


  for (cole in collection){
  all_gene_sets_msigdb <- msigdbr(species = "Mus musculus",category = cole)
  
  gsea_results <- GSEA(geneList = geneList, TERM2GENE = all_gene_sets_msigdb[, c("gs_name", "ensembl_gene")])

  enrich_dataframe <- as.data.frame(gsea_results)

  msigdb_dataframe <- as.data.frame(all_gene_sets_msigdb)
  msigdb_dataframe <- msigdb_dataframe[,c("gs_name","gs_description")]



  enrich_dataframe <- merge(x=enrich_dataframe,y=msigdb_dataframe,  by.x= "ID",  by.y= "gs_name")

  enrich_dataframe <- enrich_dataframe[!duplicated(enrich_dataframe$ID), ]

  csv_file_name= paste(file_name,cole,".csv",sep="")
  write.csv(enrich_dataframe,paste("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240828_GSEA_msigdb_collections_pitu",csv_file_name,sep="/"))


  if (nrow(enrich_dataframe) < 1 ) { 
  print ("there are not GO categories")
} else {
  plot_msigdb  <- dotplot(gsea_results, showCategory=10, split=".sign",font.size = 8) + facet_grid(.~.sign)
  #Create name to save the file
  plot_go_p=paste(file_name,cole,".png",sep="")
  #opne the png file
  png(filename=paste("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240828_GSEA_msigdb_collections_pitu",plot_go_p,sep="/"), width=400, height=700)
  plot(plot_msigdb)
  dev.off()
}
  }
  
}







hypo_heavy <- read.csv("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240725_hypothalamus/240808_hypothalamus/240808_hypothalamus_heavy_polysomes/foldchange/fold_change_heavy_polysomes.csv")
hypo_input <- read.csv("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240725_hypothalamus/240808_hypothalamus/240808_hypothalamus_input_polysomes/foldchange/fold_change_input.csv")
hypo_light <- read.csv("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240725_hypothalamus/240808_hypothalamus/240808_hypothalamus_light_polysomes/foldchange/fold_change_light_polysomes.csv")
hypo_mono <- read.csv("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240725_hypothalamus/240808_hypothalamus/240808_hypothalamus_monosomes/foldchange/fold_change_monosomes.csv")


list_hypo_deseq = list(heavy=hypo_heavy,input=hypo_input,light=hypo_light,mono=hypo_mono)



for (position in seq_along(list_hypo_deseq)){

  file=list_hypo_deseq[[position]]
  geneList <- file$foldchange
  names(geneList) <- file$genes
  file_name= names(list_hypo_deseq)[position]
  for (cole in collection){

  all_gene_sets_msigdb <- msigdbr(species = "Mus musculus",category = cole)
  gsea_results <- GSEA(geneList = geneList, TERM2GENE = all_gene_sets_msigdb[, c("gs_name", "ensembl_gene")])


  enrich_dataframe <- as.data.frame(gsea_results)

  msigdb_dataframe <- as.data.frame(all_gene_sets_msigdb)
  msigdb_dataframe <- msigdb_dataframe[,c("gs_name","gs_description")]



  enrich_dataframe <- merge(x=enrich_dataframe,y=msigdb_dataframe,  by.x= "ID",  by.y= "gs_name")

  enrich_dataframe <- enrich_dataframe[!duplicated(enrich_dataframe$ID), ]


  csv_file_name= paste(file_name,cole,".csv",sep="")
  write.csv(enrich_dataframe,paste("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240828_GSEA_msigdb_collections_hypo",csv_file_name,sep="/"))


  if (nrow(enrich_dataframe) < 1 ) { 
  print ("there are not GO categories")
} else {
  plot_msigdb <- dotplot(gsea_results, showCategory=10, split=".sign",font.size = 8) + facet_grid(.~.sign)
  #Create name to save the file
  plot_go_p=paste(file_name,cole,".png",sep="")
  #opne the png file
  png(filename=paste("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240828_GSEA_msigdb_collections_hypo",plot_go_p,sep="/"), width=400, height=700)
  plot(plot_msigdb)
  dev.off()
}
  }

}


##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################
##############################################################


#Network

pitu_heavy_cutoff <- read.csv("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240524_pituitary/240624_log005_and_pvalue005/240624_Pituitary_Heavy_polysomes/Deseq/DEG_pvalue005_log05.csv")
pitu_input_cutoff <- read.csv("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240524_pituitary/240624_log005_and_pvalue005/240624_Pituitary_input_polysomes/Deseq/DEG_pvalue005_log05.csv")
pitu_light_cutoff <- read.csv("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240524_pituitary/240624_log005_and_pvalue005/240624_Pituitary_light_polysomes/Deseq/DEG_pvalue005_log05.csv")
pitu_mono_cutoff <- read.csv("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240524_pituitary/240624_log005_and_pvalue005/240624_Pituitary_monosomes_polysomes/Deseq/DEG_pvalue005_log05.csv")


list_pitu_deseq_cutoff = list(Heavy=pitu_heavy_cutoff,Input=pitu_input_cutoff,Light=pitu_light_cutoff,Monosome=pitu_mono_cutoff)



database <- read.csv(file = "/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_input/MouseNetV2_entrez.txt",sep = '\t', header=FALSE)


net_output_folder="/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/pitu_network"

for (position in seq_along(list_pitu_deseq_cutoff)) {

file=list_pitu_deseq_cutoff[[position]]
file_name= names(list_pitu_deseq_cutoff)[position]

collection_networ = paste(net_output_folder,file_name,sep = "/")

dir.create(collection_networ)


colnames(database) <- c("gene_1","gene_2","score")

annotated_genes <- AnnotationDbi::select(org.Mm.eg.db, # database
                                            keys = file$gene,  # data to use for retrieval
                                            columns = c("ENSEMBL","GENENAME","ENTREZID","SYMBOL"), # information to retreive for given data
                                            keytype = "ENSEMBL")


entrez_numer <- annotated_genes[!is.na(annotated_genes$ENTREZID=="NA"),]
entrez_numer <- entrez_numer$ENTREZID

df_network <- data.frame()

for (i in 1:nrow(database)){
    current_row <- database[i, ]
    id_1 <- current_row[,1]
    id_2 <- current_row[,2]
    if ((id_1 %in% entrez_numer) & (id_2 %in% entrez_numer)){
        df_network <- rbind(df_network, current_row)
    }
}



#anotate the new dataframe with the symbol of each gene
symbol_1 <- AnnotationDbi::select(org.Mm.eg.db, # database
                                         keys = as.character(df_network$gene_1),  # data to use for retrieval
                                         columns = c("SYMBOL"), # information to retreive for given data
                                         keytype = "ENTREZID")
symbol_2 <- AnnotationDbi::select(org.Mm.eg.db, # database
                                         keys = as.character(df_network$gene_2),  # data to use for retrieval
                                         columns = c("SYMBOL"), # information to retreive for given data
                                         keytype = "ENTREZID")


#Select symbosl infromation and generate the complete network
df_network["symbol_1"] <- symbol_1$SYMBOL
df_network["symbol_2"] <- symbol_2$SYMBOL

#save the interaction table
write.csv(df_network,paste(collection_networ,"network.csv",sep="/"))



#generate the complete network
graph_node <- df_network[c("symbol_1","symbol_2")]
weights <- c(df_network$score)

g <- graph_from_data_frame(graph_node, directed=FALSE)
E(g)$weight <- weights # Use directed=TRUE for directed graphs

layout_fr <- layout_with_fr(g, niter = 500,)
nameplot=paste(file_name,"_netwro.png",sep = "")
png(paste(collection_networ,nameplot,sep = "/"), width = 5000, height = 5000, units = "px", res = 300)

# Plot the graph
par(mar=c(4, 3, 3, 1) + 0.1)
plot(g,  
    layout=layout_fr,               
     vertex.size=12,                       
     vertex.label.cex=1.2,
     edge.width=3,
     repel=TRUE)                          


title(main=paste(file_name,"Network",sep = " "), cex.main=5)
# Close the PNG device to save the file
dev.off()

}







hypo_heavy <- read.csv("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240725_hypothalamus/240808_hypothalamus/240808_hypothalamus_heavy_polysomes/foldchange/fold_change_heavy_polysomes.csv")
hypo_input <- read.csv("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240725_hypothalamus/240808_hypothalamus/240808_hypothalamus_input_polysomes/foldchange/fold_change_input.csv")
hypo_light <- read.csv("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240725_hypothalamus/240808_hypothalamus/240808_hypothalamus_light_polysomes/foldchange/fold_change_light_polysomes.csv")
hypo_mono <- read.csv("/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/240725_hypothalamus/240808_hypothalamus/240808_hypothalamus_monosomes/foldchange/fold_change_monosomes.csv")

list_hypo_deseq = list(heavy=hypo_heavy,input=hypo_input,light=hypo_light,mono=hypo_mono)

net_output_folder="/Users/juansola/Documents/TC3R/research/240522_polysome_profiling/240522_polysome_profiling_output/hypo_network"

for (position in seq_along(list_hypo_deseq)) {

file=list_hypo_deseq[[position]]
filter_file <- file[abs(file$foldchange) >= 0.5, ]
file_name= names(list_hypo_deseq)[position]

collection_networ = paste(net_output_folder,file_name,sep = "/")

dir.create(collection_networ)


colnames(database) <- c("gene_1","gene_2","score")

annotated_genes <- AnnotationDbi::select(org.Mm.eg.db, # database
                                            keys = filter_file$genes,  # data to use for retrieval
                                            columns = c("ENSEMBL","GENENAME","ENTREZID","SYMBOL"), # information to retreive for given data
                                            keytype = "ENSEMBL")


entrez_numer <- annotated_genes[!is.na(annotated_genes$ENTREZID=="NA"),]
entrez_numer <- entrez_numer$ENTREZID

df_network <- data.frame()

for (i in 1:nrow(database)){
    current_row <- database[i, ]
    id_1 <- current_row[,1]
    id_2 <- current_row[,2]
    if ((id_1 %in% entrez_numer) & (id_2 %in% entrez_numer)){
        df_network <- rbind(df_network, current_row)
    }
}



#anotate the new dataframe with the symbol of each gene
symbol_1 <- AnnotationDbi::select(org.Mm.eg.db, # database
                                         keys = as.character(df_network$gene_1),  # data to use for retrieval
                                         columns = c("SYMBOL"), # information to retreive for given data
                                         keytype = "ENTREZID")
symbol_2 <- AnnotationDbi::select(org.Mm.eg.db, # database
                                         keys = as.character(df_network$gene_2),  # data to use for retrieval
                                         columns = c("SYMBOL"), # information to retreive for given data
                                         keytype = "ENTREZID")


#Select symbosl infromation and generate the complete network
df_network["symbol_1"] <- symbol_1$SYMBOL
df_network["symbol_2"] <- symbol_2$SYMBOL

#save the interaction table
write.csv(df_network,paste(collection_networ,"network.csv",sep="/"))



#generate the complete network
graph_node <- df_network[c("symbol_1","symbol_2")]
weights <- c(df_network$score)

g <- graph_from_data_frame(graph_node, directed=FALSE)
E(g)$weight <- weights # Use directed=TRUE for directed graphs

layout_fr <- layout_with_fr(g, niter = 500,)
nameplot=paste(file_name,"_netwro.png",sep = "")
png(paste(collection_networ,nameplot,sep = "/"), width = 5000, height = 5000, units = "px", res = 300)

# Plot the graph
par(mar=c(4, 3, 3, 1) + 0.1)
plot(g,  
    layout=layout_fr,               
     vertex.size=12,                       
     vertex.label.cex=1.2,
     edge.width=3,
     repel=TRUE)                          


title(main=paste(file_name,"Network",sep = " "), cex.main=5)
# Close the PNG device to save the file
dev.off()

}
