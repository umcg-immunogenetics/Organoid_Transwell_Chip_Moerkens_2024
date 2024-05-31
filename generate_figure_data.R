###############################################################################
# Generate figure data  
# Marijn Berg
# m.berg@umcg.nl
###############################################################################
##########    LIBRARIES     ###################################################
library(DESeq2)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(stringr)
library(pheatmap)
library(Hmisc)
library(corrplot)
library(ReactomePA)
library(enrichplot)
library(msigdbr)
library(cowplot)
library(RColorBrewer)
###############################################################################
L2FC.filter <- 1.0
p.adj.filter <- 0.01

###############################################################################
##########  FUNCTIONS  ########################################################
###############################################################################
#FUN: ENSEMBLToEntrezSYMBOL, using ENSEMBL, gets Entrez and SYMBOL
ENSEMBLToEntrezSYMBOL <- function(f) {
  df <- f
  gene <- df$ENSEMBL
  gene.df <- bitr(gene, fromType = "ENSEMBL",
                  toType = c("ENTREZID", "SYMBOL"),
                  OrgDb = org.Hs.eg.db)
  df <- dplyr::full_join(df, gene.df, by = "ENSEMBL")
  df <- df[!duplicated(df$ENSEMBL),]
  return(df)
}

#FUN: rownames2col, changes row names using ENSEMBL column
rownames2col <- function(inDF, RowName = "ENSEMBL") {
  temp <- data.frame(rownames(inDF), inDF, row.names = NULL)
  names(temp)[1] <- RowName
  temp
}

usefulresults <- function(resultdf){
  resultdf <- resultdf[order(resultdf$padj), ]
  resultdf <- as.data.frame(resultdf)
  resultdf <- rownames2col(resultdf)
  rownames(resultdf) <- resultdf$ENSEMBL
  resultdf <- ENSEMBLToEntrezSYMBOL(resultdf)
  resultdf <- resultdf[!is.na(resultdf$padj),]
  return(resultdf)
}

#FUN: add.down_or_up.column, to add extra column as downregulated or upregulated
add.down_or_up.column <- function(df){
  df$Direction <- "Upregulated"
  df$Direction[df$log2FoldChange < 0] <- "Downregulated"
  df$Direction <- as.factor(df$Direction)
  # df$Direction <- factor(df$Direction, levels = c("Upregulated", "Downregulated"))
  return(df)
}


do_enrichments <- function(deg_res_table, res_name){
  #REACTOME enrichment
  DEG.Pathway <- compareCluster(ENTREZID~Comparison, data=deg_res_table, fun="enrichPathway", readable=T)
  DEG.Pathway@compareClusterResult$minuslog10p.adjust <- -log10(DEG.Pathway@compareClusterResult$p.adjust)
  DEG.Pathway.pw <- pairwise_termsim(DEG.Pathway)
  
  #GO: Biological process
  DEG.bp <- compareCluster(ENTREZID~Comparison, data=deg_res_table, fun="enrichGO", OrgDb = org.Hs.eg.db, ont = "BP", readable=T) 
  DEG.bp@compareClusterResult$minuslog10p.adjust <- -log10(DEG.bp@compareClusterResult$p.adjust)
  DEG.bp.pw <- pairwise_termsim(DEG.bp)
  DEG.bp.pw.simp <- simplify(DEG.bp.pw, cutoff=0.7, by="p.adjust", select_fun=min)

  dir.create("./Results/objects/enrichments/DEG/", showWarnings = FALSE, recursive = TRUE)
  #Saving enrichment objects
  saveRDS(DEG.Pathway, paste0("./Results/objects/enrichments/DEG/", res_name, "Pathway.RDS"))
  saveRDS(DEG.bp, paste0("./Results/objects/enrichments/DEG/", res_name, "bp.RDS"))
}

###############################################################################
geneSymbol <- read.table("Data/features.tsv", sep = "\t")[,1:2]
colnames(geneSymbol) <- c("ensembl", "symbol")
rownames(geneSymbol) <- geneSymbol$ensembl

h_gene_set <- msigdbr(species = "Homo sapiens", category = "H")
msigdbr_t2g <- h_gene_set %>% 
  dplyr::distinct(gs_name, entrez_gene) %>% 
  as.data.frame()

###############################################################################
# Different model systems and stimulations
ms_counts <- as.matrix(read.csv(file = "Data/COT_MC_Final_Iris-Moerkens.expression.genelevel.v75.htseq.txt.table", sep="\t", row.names = "probe"))

metadata <- read.csv(file = "Data/MetaData_COTMCFinal_RNAseq_csv.csv",sep=";",stringsAsFactors = TRUE)
# make a single coded name for the samples containing useful info
metadata$codename <- paste0(str_sub(as.character(metadata$Model_system),1,1),  # first character of the model system
                            "_",
                            as.character(metadata$Medium_Bottom_Top),
                            "_",
                            str_sub(as.character(metadata$Cell_line), 6, 7))
metadata$codename
rownames(metadata) <- metadata$codename

metadata$conc.condition <- factor(paste0(str_sub(as.character(metadata$Model_system),1,1),  # first character of the model system
                                         "_",
                                         as.character(metadata$Medium_Bottom_Top)))

metadata <- metadata[,c('Cell_line', 'conc.condition', 'Model_system', 'Medium_Bottom_Top', 'codename', 'Sample')]

ms_counts <- ms_counts[,as.character(metadata$Sample)]
colnames(ms_counts) <- metadata$codename
dds <- DESeqDataSetFromMatrix(countData = ms_counts,
                              colData = metadata,
                              design= ~ Cell_line + Medium_Bottom_Top)

keep <- rowSums(counts(dds) >= 10) >= 3

dds <- dds[keep,]
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind=FALSE)

###############################################################################
###################v Organoid EMc-cDM comparison     ##########################
###############################################################################
dds <- DESeq(dds)
resultsNames(dds)

dds_org <- dds[, dds$Model_system %in% c("Organoid")]
dds_org$Model_system <- droplevels(dds_org$Model_system)
dds_org$Medium_Bottom_Top <- droplevels(dds_org$Medium_Bottom_Top)
dds_org <- DESeq(dds_org)

resultsNames(dds_org)
org.res <- usefulresults(results(dds_org, contrast =c('Medium_Bottom_Top', 'DM_DM', 'EM_EM')))

org_counts <- getVarianceStabilizedData(dds_org)
write.table(org_counts, file = "Results/figure_data/organoid_normalized_counts.csv")


org.res <- add.down_or_up.column(org.res)
org.res$Comparison <- "EE_DD"
org.res$EE_avg_expr <- rowMeans(org_counts[org.res$ENSEMBL ,grep("EM_EM", colnames(org_counts))])
org.res$DD_avg_expr <- rowMeans(org_counts[org.res$ENSEMBL ,grep("DM_DM", colnames(org_counts))])
write.table(org.res, file = "Results/figure_data/organoid_EE_DD_DEG_results.csv")

###############################################################################
###################v Transwell all conditions DGE     ##########################
###############################################################################
dds_twell <- dds[, dds$Model_system %in% c("Transwell")]
dds_twell$Model_system <- droplevels(dds_twell$Model_system)
dds_twell$Medium_Bottom_Top <- droplevels(dds_twell$Medium_Bottom_Top)
dds_twell <- DESeq(dds_twell)

###############################################################################
resultsNames(dds_twell)

teed.res <- usefulresults(results(dds_twell, contrast = c('Medium_Bottom_Top', 'EM_DM', 'EM_EM')))
teed.res <- add.down_or_up.column(teed.res)
teed.res$Comparison <- "EE_ED"

###############################################################################
tedd.res <- usefulresults(results(dds_twell, contrast = c('Medium_Bottom_Top', 'DM_DM', 'EM_EM')))
tedd.res <- add.down_or_up.column(tedd.res)
tedd.res$Comparison <- "EE_DD"

###############################################################################
tddd.res <- usefulresults(results(dds_twell, contrast = c('Medium_Bottom_Top', 'DM_DM', 'EM_DM')))
tddd.res <- add.down_or_up.column(tddd.res)
tddd.res$Comparison <- "ED_DD"

###############################################################################
###############################################################################
org_counts <- getVarianceStabilizedData(dds_twell)
write.table(org_counts, file = "Results/figure_data/transwell_normalized_counts.csv")

DEG.merged <- rbind(rbind(teed.res, tedd.res), tddd.res)
DEG.merged$EE_avg_expr <- rowMeans(org_counts[DEG.merged$ENSEMBL ,grep("EM_EM", colnames(org_counts))])
DEG.merged$ED_avg_expr <- rowMeans(org_counts[DEG.merged$ENSEMBL ,grep("EM_DM", colnames(org_counts))])
DEG.merged$DD_avg_expr <- rowMeans(org_counts[DEG.merged$ENSEMBL ,grep("DM_DM", colnames(org_counts))])

write.table(DEG.merged, file = "Results/figure_data/transwell_DEG_results.csv")

###############################################################################
###################  Model system comparison  #################################
###############################################################################
dds_test <- dds[, (dds$Model_system %in% c("Transwell") & dds$Medium_Bottom_Top %in% c("EM_DM")) | 
                   (dds$Model_system %in% c("Organoid") & dds$Medium_Bottom_Top %in% c("DM_DM")) |
                   (dds$Model_system %in% c("Chip") & dds$Medium_Bottom_Top %in% c("EM_DM"))]

design(dds_test) <- ~ Cell_line + Model_system
dds_test$Model_system <- droplevels(dds_test$Model_system)
dds_test$Medium_Bottom_Top <- droplevels(dds_test$Medium_Bottom_Top)
dds_test <- DESeq(dds_test)

odted.res <- usefulresults(results(dds_test, contrast = c('Model_system', 'Transwell', 'Organoid')))
odted.res <- add.down_or_up.column(odted.res)
odted.res$Comparison <- "Organoid_DM_Transwell_EMDM"

cedted.res <- usefulresults(results(dds_test, contrast = c('Model_system', 'Transwell', 'Chip')))
cedted.res <- add.down_or_up.column(cedted.res)
cedted.res$Comparison <- "Chip_EMDM_Transwell_EMDM"

odced.res <- usefulresults(results(dds_test, contrast = c('Model_system', 'Chip', 'Organoid')))
odced.res <- add.down_or_up.column(odced.res)
odced.res$Comparison <- "Organoid_DM_Chip_EMDM"

DEG.merged <- rbind(rbind(odted.res, cedted.res), odced.res) 

org_counts <- getVarianceStabilizedData(dds_test)
write.table(org_counts, file = "Results/figure_data/model_system_normalized_counts.csv")

colnames(org_counts)

DEG.merged$Organoid_DM_avg_expr    <- rowMeans(org_counts[DEG.merged$ENSEMBL ,grep("O_DM_DM", colnames(org_counts))])
DEG.merged$Chip_EMDM_avg_expr      <- rowMeans(org_counts[DEG.merged$ENSEMBL ,grep("C_EM_DM", colnames(org_counts))])
DEG.merged$Transwell_EMDM_avg_expr <- rowMeans(org_counts[DEG.merged$ENSEMBL ,grep("T_EM_DM", colnames(org_counts))])
write.table(DEG.merged, "Results/figure_data/model_systems_DEG.tsv", sep = "\t", row.names = FALSE)

###############################################################################
# Organoid
###############################################################################
DEG.merged <- read.table("Results/figure_data/organoid_EE_DD_DEG_results.csv")
DEG.merged <- DEG.merged[(abs(DEG.merged$log2FoldChange) > L2FC.filter) & 
                           (DEG.merged$padj < p.adj.filter), ]

org_counts <- read.table("Results/figure_data/organoid_normalized_counts.csv", row.names = 1, header = TRUE)
org_counts <- org_counts[unique(DEG.merged$ENSEMBL) ,order(colnames(org_counts))]

clusters <- hclust(dist(t(scale(t(org_counts)))))
gene_cluster_idents <- cutree(tree = clusters, k = 2)

DEG.merged[DEG.merged$ENSEMBL %in% (names(gene_cluster_idents[gene_cluster_idents == 1])), "Comparison"] <- "cluster_1"
DEG.merged[DEG.merged$ENSEMBL %in% (names(gene_cluster_idents[gene_cluster_idents == 2])), "Comparison"] <- "cluster_2"
write.table(DEG.merged, file = "Results/figure_data/Organoid_sig_DEG_w_clusters.csv", row.names = F)

#GO: Biological process
DEG.bp <- compareCluster(ENTREZID~Comparison, data=DEG.merged, fun="enrichGO", OrgDb = org.Hs.eg.db, ont = "BP", readable=T) 
DEG.bp@compareClusterResult$minuslog10p.adjust <- -log10(DEG.bp@compareClusterResult$p.adjust)
saveRDS(DEG.bp, "Results/figure_data/organoid_clusters_GO_BP_results.RDS")

###############################################################################
# Transwell
###############################################################################
DEG.merged <- read.table("Results/figure_data/transwell_DEG_results.csv")
DEG.merged <- DEG.merged[(abs(DEG.merged$log2FoldChange) > L2FC.filter) & 
                           (DEG.merged$padj < p.adj.filter), ]

org_counts <- read.table("Results/figure_data/transwell_normalized_counts.csv", row.names = 1, header = TRUE)
org_counts <- org_counts[unique(DEG.merged$ENSEMBL) ,order(colnames(org_counts))]

pheatmap(org_counts, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = FALSE, scale = "row", cutree_rows = 2)

clusters <- hclust(dist(t(scale(t(org_counts)))))
gene_cluster_idents <- cutree(tree = clusters, k = 2)

DEG.merged[DEG.merged$ENSEMBL %in% (names(gene_cluster_idents[gene_cluster_idents == 1])), "Comparison"] <- "cluster_1"
DEG.merged[DEG.merged$ENSEMBL %in% (names(gene_cluster_idents[gene_cluster_idents == 2])), "Comparison"] <- "cluster_2"
write.table(DEG.merged, file = "Results/figure_data/Transwell_sig_DEG_w_clusters.csv", row.names = F)

#GO: Biological process
DEG.bp <- compareCluster(ENTREZID~Comparison, data=DEG.merged, fun="enrichGO", OrgDb = org.Hs.eg.db, ont = "BP", readable=T) 
DEG.bp@compareClusterResult$minuslog10p.adjust <- -log10(DEG.bp@compareClusterResult$p.adjust)
saveRDS(DEG.bp, "Results/figure_data/Transwell_clusters_GO_BP_results.RDS")

###############################################################################
# Model systems
###############################################################################
DEG.merged <- read.table("Results/figure_data/model_systems_DEG.tsv", header = TRUE)
DEG.merged <- DEG.merged[(abs(DEG.merged$log2FoldChange) > L2FC.filter) & 
                           (DEG.merged$padj < p.adj.filter), ]

org_counts <- read.table("Results/figure_data/model_system_normalized_counts.csv", row.names = 1, header = TRUE)
org_counts <- org_counts[unique(DEG.merged$ENSEMBL) ,order(colnames(org_counts))]

clusterHM <- pheatmap(org_counts, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = FALSE, scale = "row", cutree_rows = 3)

clusters <- cutree(clusterHM$tree_row, k=3)[clusterHM$tree_row[["order"]]]
annot_row <- data.frame(row.names = names(clusters),
                        cluster = as.factor(clusters))
clusterHM <- pheatmap(org_counts, cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = FALSE, scale = "row", cutree_rows = 3, annotation_row = annot_row)

clusters <- hclust(dist(t(scale(t(org_counts)))))
gene_cluster_idents <- cutree(tree = clusters, k = 3)

DEG.merged[DEG.merged$ENSEMBL %in% (names(gene_cluster_idents[gene_cluster_idents == 1])), "Comparison"] <- "cluster_1"
DEG.merged[DEG.merged$ENSEMBL %in% (names(gene_cluster_idents[gene_cluster_idents == 2])), "Comparison"] <- "cluster_2"
DEG.merged[DEG.merged$ENSEMBL %in% (names(gene_cluster_idents[gene_cluster_idents == 3])), "Comparison"] <- "cluster_3"
DEG.merged[DEG.merged$ENSEMBL %in% (names(gene_cluster_idents[gene_cluster_idents == 4])), "Comparison"] <- "cluster_4"
DEG.merged[DEG.merged$ENSEMBL %in% (names(gene_cluster_idents[gene_cluster_idents == 5])), "Comparison"] <- "cluster_5"

write.table(DEG.merged, file = "Results/figure_data/Model_systems_sig_DEG_w_3_clusters.csv", row.names = F)

#GO: Biological process
DEG.bp <- compareCluster(ENTREZID~Comparison, data=DEG.merged, fun="enrichGO", OrgDb = org.Hs.eg.db, ont = "BP", readable=T) 
DEG.bp@compareClusterResult$minuslog10p.adjust <- -log10(DEG.bp@compareClusterResult$p.adjust)
emapplot(pairwise_termsim(DEG.bp))
saveRDS(DEG.bp, "Results/figure_data/Model_system_clusters_GO_BP_results_3_clusters.RDS")

###############################################################################
# distance matrix counts
geneSymbol <- read.table("Data/features.tsv", sep = "\t")[,1:2]
colnames(geneSymbol) <- c("ensembl", "symbol")
rownames(geneSymbol) <- geneSymbol$ensembl
###############################################################################
ms_counts <- as.matrix(read.csv(file = "Data/COT_MC_Final_Iris-Moerkens.expression.genelevel.v75.htseq.txt.table", sep="\t", row.names = "probe"))

metadata <- read.csv(file = "Data/MetaData_COTMCFinal_RNAseq_csv.csv",sep=";",stringsAsFactors = TRUE)
# make a single coded name for the samples containing useful info
metadata$codename <- paste0(str_sub(as.character(metadata$Model_system),1,1),  # first character of the model system
                            "_",
                            as.character(metadata$Medium_Bottom_Top),
                            "_",
                            str_sub(as.character(metadata$Cell_line), 6, 7))
metadata$codename
rownames(metadata) <- metadata$codename

metadata$conc.condition <- factor(paste0(str_sub(as.character(metadata$Model_system),1,1),  # first character of the model system
                                         "_",
                                         as.character(metadata$Medium_Bottom_Top)))

metadata <- metadata[,c('Cell_line', 'conc.condition', 'Model_system', 'Medium_Bottom_Top', 'codename', 'Sample')]

ms_counts <- ms_counts[,as.character(metadata$Sample)]
colnames(ms_counts) <- metadata$codename

oslo_counts <- as.matrix(read.csv("Data/RNA_matrix.table", sep="\t", row.names = 1))
oslo_counts <- oslo_counts[,grep("CON", colnames(oslo_counts))]

ct_perc <- read.table("Data/proportion_epithelial_immune.csv", sep = "\t", row.names = 1, header = TRUE)


el_counts <- as.matrix(read.table("Data/pseudobulk_single_cell_2.csv", sep = "\t"))
colnames(el_counts) <- str_replace(colnames(el_counts), "RNA.", "sc")
colnames(el_counts) <- str_replace_all(colnames(el_counts), "\\.", "_")
el_counts <- el_counts[,!(colnames(el_counts) %in% c("scT036POS_Pediatric", "scT110NEG_Pediatric"))]
###############################################################################
# single cell model system pseudo bulk
sc_counts <- as.matrix(read.table("Data/pseudobulk_single_cell.csv", sep = "\t"))[,1:9]
colnames(sc_counts) <- str_replace(colnames(sc_counts), "RNA.Epithelial.cell", "sc")

geneSymbol <- geneSymbol[geneSymbol$symbol %in% rownames(el_counts),]
commonGenes <- intersect(geneSymbol$ensembl, intersect(rownames(oslo_counts), rownames(ms_counts)))
commonSymbols <- geneSymbol[commonGenes, 2]

ms_counts <- ms_counts[commonGenes,]
rownames(ms_counts) <- commonSymbols

oslo_counts <- oslo_counts[commonGenes,]
rownames(oslo_counts) <- commonSymbols

el_counts <- el_counts[commonSymbols,]

sc_counts <- sc_counts[commonSymbols,]
###############################################################################

counts <- cbind(ms_counts, sc_counts, cbind(oslo_counts, el_counts))

sampletype <- as.factor(c(rep("Model_system", 21), rep("scRNA_model", 9), rep("biopsy", 25), rep("scRNA", 29) ))
mergedmeta <- data.frame(sampletype)

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = mergedmeta,
                              design= ~ sampletype)

keep <- rowSums(counts(dds) >= 10) >= 3
table(keep)
dds <- dds[keep,]

dds <- DESeq(dds)

normalized_counts <- getVarianceStabilizedData(dds)
topBiopsyGenes <- rownames(normalized_counts[order(rowSums(normalized_counts[,22:46]), decreasing = TRUE),])[1:5000]
###############################################################################

deMarkersAll <- read.csv("Data/GCA_epithelial_markers_full.csv", sep = "\t", header = TRUE)
epithelialMarkerGenes <- rownames(deMarkersAll[deMarkersAll$p_val_adj < 0.05 & deMarkersAll$avg_log2FC > 2.6871011,])
epithelialMarkerGenes <- epithelialMarkerGenes[epithelialMarkerGenes %in% rownames(normalized_counts)]

figuresubsetcounts <- normalized_counts[epithelialMarkerGenes,]
write.table(figuresubsetcounts, file = "Results/figure_data/distance_matrix_counts.csv")




###############################################################################
ms_counts <- as.matrix(read.csv(file = "Data/COT_MC_Final_Iris-Moerkens.expression.genelevel.v75.htseq.txt.table", sep="\t", row.names = "probe"))

metadata <- read.csv(file = "Data/MetaData_COTMCFinal_RNAseq_csv.csv",sep=";",stringsAsFactors = TRUE)
# make a single coded name for the samples containing useful info
metadata$codename <- paste0(str_sub(as.character(metadata$Model_system),1,1),  # first character of the model system
                            "_",
                            as.character(metadata$Medium_Bottom_Top),
                            "_",
                            str_sub(as.character(metadata$Cell_line), 6, 7))
metadata$codename
rownames(metadata) <- metadata$codename

metadata$conc.condition <- factor(paste0(str_sub(as.character(metadata$Model_system),1,1),  # first character of the model system
                                         "_",
                                         as.character(metadata$Medium_Bottom_Top)))

metadata <- metadata[,c('Cell_line', 'conc.condition', 'Model_system', 'Medium_Bottom_Top', 'codename', 'Sample')]
metadata <- metadata[!(metadata$conc.condition == 'C_EM_EM'),]
ms_counts <- ms_counts[,as.character(metadata$Sample)]

colnames(ms_counts) <- metadata$codename

dds <- DESeqDataSetFromMatrix(countData = ms_counts,
                              colData = metadata,
                              design= ~ conc.condition)

# keep <- rowSums(counts(dds) >= 10) >= 3
# table(keep)
# dds <- dds[keep,]

dds <- DESeq(dds)

normalized_counts <- getVarianceStabilizedData(dds)
write.table(normalized_counts, file = "Results/figure_data/supp_fig1_counts.csv")
