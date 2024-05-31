
###############################################################################
# Generate the transwell figure for Renee's manuscript
# Marijn Berg m.berg@umcg.nl
###############################################################################
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggplot2)
library(stringr)
library(ComplexHeatmap)
library(Hmisc)
library(enrichplot)
library(msigdbr)
library(cowplot)
library(RColorBrewer)
library(gridExtra)
library(magick)
library(ggplotify)
library(ggrepel)
library(scales)
library(circlize)
library(ggforce)
library(shadowtext)
###############################################################################
geneSymbol <- read.table("Data/features.tsv", sep = "\t")[,1:2]
colnames(geneSymbol) <- c("ensembl", "symbol")

# I checked if the removed genes are needed below, they're not
geneSymbol <- geneSymbol[!duplicated(geneSymbol$symbol),]
rownames(geneSymbol) <- geneSymbol$symbol

###############################################################################
#####Gene selection Transwells
DNA_replication <- c('TSPYL2', 'STN1', 'CCDC88A', 'HMGA1', 'CDC7', 'CDC6', 'E2F7', 'POLD2', 'EREG', 'DDX21', 'RRM2B', 'MCM5')
Cell_division <- c('SMOC2','LGR5','RGMB', 'PLK3', 'EIF4EBP1', 'HES1', 'PLK2', 'WEE1', 'MDC1', 'ASNS', 'PDGFB', 'GJC2')
Neurogenesis <- c('TUBB2B','NCAM1', 'SEMA7A', 'ID1', 'LRP8', 'LPAR3', 'CRABP2', 'EPHB2', 'TSPO', 'NRP1', 'DBN1', 'SEMA3F', 'ETV5', 'SERPINE2', 'RGMA', 'MYB')
Nutrient_metabolism <- c('ANPEP', 'LCT', 'ENPEP', 'PEPD', 'ALDOB')

Nutrient_transport <- c("APOA1","APOB", "APOA2", "APOH", "APOC3", 'SLC2A2','SLC5A1',"SLC5A9", "SLC5A11", 'SLC15A1','SLC7A9','SLC36A1',"SLC26A2","SLC26A3","SLC13A2","SLC34A3", "SLC23A1", "SLC23A3", 'SLC11A2', "SLC30A2")
Xenobiotic_metabolism <- c('CYP1A1', "CYP3A4","CYP2C8", "CYP2C9",  "CYP2C19",  "CYP2C58P",'CYP2D6', "CYP2U1", "CYP27A1",'CYP2W1',  "CYP4F2",'ABCB1','ABCC2', 'MAOB', "UGT2A3",  "UGT2B15", "UGT2B17", "UGT2B4",  "UGT2B10")
Hormone_response <- c('CHGA', 'GCG', 'SCT')

DNA_replication <- DNA_replication[order(DNA_replication)]
Cell_division <- Cell_division[order(Cell_division)]
Neurogenesis <- Neurogenesis[order(Neurogenesis)]
Nutrient_metabolism <- Nutrient_metabolism[order(Nutrient_metabolism)]


Nutrient_transport <- Nutrient_transport[order(Nutrient_transport)]
Xenobiotic_metabolism <- Xenobiotic_metabolism[order(Xenobiotic_metabolism)]
Hormone_response <- Hormone_response[order(Hormone_response)]


all_genes <- c(DNA_replication, Cell_division, Neurogenesis, Nutrient_metabolism, 
               Nutrient_transport, Xenobiotic_metabolism, Hormone_response)

gene_sets <- c(rep('DNA\nreplication', length(DNA_replication)),
               rep('Cell\ndivision', length(Cell_division)),
               rep('Neurogenesis', length(Neurogenesis)),
               rep('Nutrient\nmetabolism', length(Nutrient_metabolism)),
               rep('Nutrient\ntransport', length(Nutrient_transport)),
               rep('Xenobiotic\nmetabolism', length(Xenobiotic_metabolism)),
               rep('Hormone\nresponse', length(Hormone_response))
)

gene_table <- data.frame(all_genes, gene_sets)
gene_table <- gene_table[gene_table$all_genes %in% geneSymbol$symbol,]

gene_table$ensemble <- geneSymbol[gene_table$all_genes, 1]
###############################################################################
# FIGURE SETTINGS
###############################################################################
plotseed <- 61    # emapplot uses some randomization so if you want the same figure
bigfontsize <- 10
medfontsize <- 8
smlfontsize <- 5
col_fun = colorRamp2(c(-2, 0, 2), c("purple", "white", "darkgreen"))
col_fun(seq(-3, 3))

###############################################################################
counts  <- read.table("Results/figure_data/transwell_normalized_counts.csv")
deRes   <- read.table("Results/figure_data/Transwell_sig_DEG_w_clusters.csv", header = TRUE)
pathway <- readRDS("Results/figure_data/Transwell_clusters_GO_BP_results.RDS")
###############################################################################
plotRes <- pairwise_termsim(pathway)


set.seed(plotseed)
pathplot <- emapplot(plotRes,
                     cex_label_category=0.8,
                     cex_category=1, # !RELATIVE! size of the dots
                     legend_n=2, showCategory = 20,
                     group_legend=FALSE) +
  theme(legend.position= c(0.1, 0.2)) # + 

pathplot
###############################################################################
clusthmcounts <- t(scale(t(counts[unique(deRes$ENSEMBL),order(colnames(counts))])))

clusters <- hclust(dist(clusthmcounts))
gene_cluster_idents <- cutree(tree = clusters, k = 2)
clusterset <- paste0("Cluster ", gene_cluster_idents[clusters$order])
clusthmcounts <- clusthmcounts[clusters$order,]

annot_group <- data.frame(row.names = colnames(clusthmcounts),
                          sampletype = factor(as.factor(c(rep("DM+D+P", 3), rep("EM-DM+D+P", 3), rep("EM", 3))), 
                                              ordered = TRUE, levels=c("EM", "EM-DM+D+P", "DM+D+P")))
clusterHM <- as.ggplot(
  Heatmap(clusthmcounts, name="expression", cluster_rows = FALSE, cluster_columns = FALSE, #col = col_fun,
          show_column_names = FALSE,
          show_row_names = FALSE, 
          column_split = annot_group$sampletype, 
          column_title_gp = gpar(fontsize=6),
          row_title_gp = gpar(fill = c("#ff6e67", "#2596be")),
          row_split = clusterset,
          heatmap_legend_param = list(title_position = "lefttop-rot", grid_width = unit(2, "mm")))
)
clusterHM
###############################################################################
hmcounts <- scale(t(counts[gene_table$ensemble, order(colnames(counts))]))
colnames(hmcounts) <- gene_table$all_genes

annot_gene <- data.frame(row.names = gene_table$all_genes,
                         process = factor(gene_table$gene_sets, ordered = TRUE, 
                                          levels=c("DNA\nreplication", "Cell\ndivision",  "Neurogenesis", "Nutrient\nmetabolism",
                                                   "Nutrient\ntransport", "Xenobiotic\nmetabolism", 'Hormone\nresponse')
                         ))

annot_group <- data.frame(row.names = rownames(hmcounts),
                          sampletype = factor(as.factor(c(rep("DM+D+P", 3), rep("EM-DM+D+P", 3), rep("EM", 3))), 
                          ordered = TRUE, levels=c("EM", "EM-DM+D+P", "DM+D+P")))

markerHM <- as.ggplot(
  Heatmap(hmcounts, name="expression", cluster_rows = FALSE, cluster_columns = FALSE, 
          column_names_gp = gpar(fontsize = smlfontsize), column_names_rot = 90,
          show_row_names = FALSE, 
          column_split = annot_gene$process, 
          row_split = annot_group$sampletype, 
          column_title_gp = gpar(fontsize=medfontsize),
          row_title_gp = gpar(fontsize=medfontsize),
          row_title_rot = 0, 
          heatmap_legend_param = list(title_position = "lefttop-rot", grid_width = unit(2, "mm")))
)
markerHM

###############################################################################
toprow <- plot_grid(clusterHM, pathplot, nrow = 1, labels = c("A", "B"), rel_widths = c(1,3))
botrow <- plot_grid(markerHM, labels = c("C"))

organoidfig <- plot_grid(toprow, botrow, nrow = 2, rel_heights = c(3,1) )

png("Figures/Figure_transwell_test3.png", width = 20, height = 20, units = 'cm', res = 300)
grid.arrange(organoidfig) # Make plot
dev.off()

pdf("Figures/Transwell_figure.pdf", width = 20/2.54, height = 20/2.54)
grid.arrange(organoidfig) # Make plot
dev.off()

