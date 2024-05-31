
###############################################################################
# Generate the model system figure for Renee's manuscript
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
####Gene selection Model systems
####Gene selection Model systems
Mineral_metabolism <- c('MT1M', 'MT1F', 'MT1X', 'MT1H', "MT1L")
Embryonic_development <- c("HOXA9",    "HOXB9",    "HOXB8",    "HOXB4",    "HOXB2",    "HOXB6",    "HOXA1", 'SHH', 'IHH', 'SMO', 'HES1', 'SUFU', 'VANGL2', 'FOXA1', 'SOX4', 'SOX17', 'WNT5A', 'NOTCH1')
ECM_organization <- c('LAMB1','LAMB2', 'COL2A1',"COL4A1", "COL4A2",'COL4A5', 'TGFB1', 'TGFB2', 'TIMP1', 'MMP19' , 'ACTA2')
Neurogenesis <- c('TUBB2B','STMN1', 'NCAM1', 'CRABP2', 'RGMA', 'NRP1', 'SEMA3F', 'EPHA6', 'EFNA1')
Cell_division <- c('SMOC2', 'MKI67', 'TOP2A', 'CENPF', 'UBE2C', 'NEK6', 'TTK', 'MAD2L1', 'PTTG1', 'BUB1B', 'KNTC1', 'CCNB1')
Nutrient_metabolism <- c('ANPEP', 'SI', 'LCT', 'MGAM', 'ENPEP', 'DPP4', 'PEPD', 'ALPI', 'LIPA', 'TMPRSS15', 'ALDOB')
Nutrient_transport <- c("APOA4",  "APOA2",   "APOA1",   "APOB",    "APOC2", 'FABP2', "SLC2A2","SLC2A5","SLC2A7", 'SLC5A1',"SLC5A4","SLC5A12",'SLC15A1', 'SLC7A9',"SLC6A19", "SLC13A1",'SLC26A3', "SLC28A2","SLC10A5", "SLC46A3", 'SLC39A4')
Xenobiotic_metabolism <- c("CYP1A1","CYP3A4", 'CYP3A5', "CYP3A7", 'CYP2C9', "CYP2C18",'CYP2D6','CYP2J2', "CYP7A1",  "CYP4F2",  "CYP2U1",  "CYP4F12",  "CYP17A1",'ABCB1', 'ABCC2','ABCG2', 'MAOB', "UGT2B4",  "UGT2B7",  "UGT2B17", "UGT1A1")

all_genes <- c(Embryonic_development, 
               ECM_organization, 
               Neurogenesis, 
               Cell_division,
               Mineral_metabolism,
               Nutrient_metabolism,
               Nutrient_transport,
               Xenobiotic_metabolism)

gene_sets <- c(rep('Embryonic\ndevelopment', length(Embryonic_development)),
               rep('ECM\norganization', length(ECM_organization)),
               rep('Neurogenesis', length(Neurogenesis)),
               rep('Cell\ndivision', length(Cell_division)),
               
               rep('Mineral\nmetabolism', length(Mineral_metabolism)),
               rep('Nutrient\nmetabolism', length(Nutrient_metabolism)),
               rep('Nutrient\ntransport', length(Nutrient_transport)),
               rep('Xenobiotic\nmetabolism', length(Xenobiotic_metabolism))
)

gene_table <- data.frame(all_genes, gene_sets)
gene_table <- gene_table[gene_table$all_genes %in% geneSymbol$symbol,]

gene_table$ensemble <- geneSymbol[gene_table$all_genes, 1]
table(duplicated(gene_table$all_genes))
###############################################################################
# FIGURE SETTINGS
###############################################################################
plotseed <-  19   # emapplot uses some randomization so if you want the same figure
bigfontsize <- 10
medfontsize <- 8
smlfontsize <- 5
scale_col_fun = colorRamp2(c(-2, 0, 2), c("yellow", "blue", 'darkblue'))
# scale_col_fun = colorRamp2(c(0, 2.5, 5), c("yellow", "blue", 'darkblue'))

###############################################################################
counts  <- read.table("Results/figure_data/model_system_normalized_counts.csv")
deRes   <- read.table("Results/figure_data/Model_systems_sig_DEG_w_3_clusters.csv", header = TRUE)
pathway <- readRDS("Results/figure_data/Model_system_clusters_GO_BP_results_3_clusters.RDS")
###############################################################################
plotRes <- pairwise_termsim(pathway)
# 104, 132, 160, 158
plotseed <-  132
set.seed(plotseed)
pathplot <- emapplot(plotRes,
                     cex_label_category=0.8,
                     cex_category=1, # !RELATIVE! size of the dots
                     legend_n=2, showCategory = 20,
                     group_legend=FALSE) +
  theme(legend.position= c(0.1, 0.175)) #+ 
pathplot

# for(i in 100:1100){
#   set.seed(i)
#   png(filename = paste0("Results/testEMAPs/modelsystems/test_", i, ".png"), width = 500, height = 500)
#   pathplot <- emapplot(plotRes,
#                        cex_label_category=0.8,
#                        cex_category=1, # !RELATIVE! size of the dots
#                        legend_n=2, showCategory = 20,
#                        group_legend=FALSE) +
#     theme(legend.position= c(0.1, 0.175)) #+ 
#   
#   print(pathplot)
#   dev.off()
# }

###############################################################################
clusthmcounts <- t(scale(t(counts[unique(deRes$ENSEMBL),order(colnames(counts))])))

clusters <- hclust(dist(clusthmcounts))
gene_cluster_idents <- cutree(tree = clusters, k = 3)
clusterset <- factor(paste0("Cluster ", gene_cluster_idents[clusters$order]), ordered = TRUE, levels = c("Cluster 3", "Cluster 1", "Cluster 2"))
clusthmcounts <- clusthmcounts[clusters$order,]

annot_group <- data.frame(row.names = colnames(clusthmcounts),
                          sampletype = factor(as.factor(c(rep("Chip", 3), rep("Organoid", 3), rep("Transwell", 3))),
                                              ordered = TRUE, levels=c("Organoid", "Transwell", "Chip")))
clusterHM <- as.ggplot(
  Heatmap(clusthmcounts, name="expression", cluster_rows = FALSE, cluster_columns = FALSE, #col = col_fun,
          show_column_names = FALSE,
          show_row_names = FALSE, 
          column_split = annot_group$sampletype, 
          column_title_gp = gpar(fontsize=medfontsize),
          row_title_gp = gpar(fill = c("#2596be", "#ff6e67", "#7cae00")),
          row_split = clusterset, cluster_row_slices = FALSE,
          heatmap_legend_param = list(title_position = "lefttop-rot", grid_width = unit(2, "mm")))
)
clusterHM
###############################################################################
hmcounts <- scale(t(counts[unique(gene_table$ensemble), order(colnames(counts))]))
# hmcounts <- t(counts[gene_table$ensemble, order(colnames(counts))])

colnames(hmcounts) <- gene_table$all_genes

annot_gene <- data.frame(row.names = gene_table$all_genes,
                         process = factor(gene_table$gene_sets, ordered = TRUE, 
                                          levels=c("Embryonic\ndevelopment", 
                                                   "ECM\norganization", "Neurogenesis", "Cell\ndivision","Mineral\nmetabolism",  
                                                   "Nutrient\ntransport",  "Nutrient\nmetabolism",  "Xenobiotic\nmetabolism")
                         ))

annot_group <- data.frame(row.names = rownames(hmcounts),
                          sampletype = factor(as.factor(c(rep("Chip", 3), rep("Organoid", 3), rep("Transwell", 3))),
                                              ordered = TRUE, levels=c("Organoid", "Transwell", "Chip")))

markerHM <- as.ggplot(
  Heatmap(hmcounts, name="expression", cluster_rows = FALSE, cluster_columns = FALSE, 
          column_names_gp = gpar(fontsize = smlfontsize), 
          show_row_names = FALSE, 
          column_split = annot_gene$process, 
          row_split = annot_group$sampletype, 
          column_title_gp = gpar(fontsize=medfontsize),
          row_title_gp = gpar(fontsize=medfontsize),
          row_title_rot = 0, 
          column_names_rot = 90,
          heatmap_legend_param = list(title_position = "lefttop-rot", grid_width = unit(2, "mm")))
)
markerHM
###############################################################################
distanceCounts <- read.table("Results/figure_data/distance_matrix_counts.csv")
sample_distances <- as.matrix(dist(t(distanceCounts)))[1:21, 31:84]

sample_distances <- sample_distances[c(grep('T_EM_DM', rownames(sample_distances)),
                                       grep('C_EM_DM', rownames(sample_distances)),
                                       grep('O_DM_DM', rownames(sample_distances))),]
sample_distances <- scale(sample_distances)
sample_distances <- sample_distances[order(rownames(sample_distances)),]


rownames(sample_distances)
modeltype <- factor(as.factor(c(rep("Chip", 3), rep("Organoid", 3), rep("Transwell", 3))),
                    ordered = TRUE, levels=c("Organoid", "Transwell", "Chip"))

sampletype <- factor(as.factor(c(rep('Epithelial layer duodenal biopsy(adult)', 25),
                rep('adult', 4),
                rep('first\ntrimester', 11),
                rep('second\ntrimester', 5),
                rep('pediatric', 9))), 
                ordered = TRUE, levels = c('first\ntrimester', 'second\ntrimester', 'pediatric', 'adult', 'Epithelial layer duodenal biopsy(adult)'))

distHM <- as.ggplot(
  Heatmap(sample_distances, name="distance", cluster_rows = FALSE, cluster_columns = FALSE, col = scale_col_fun,
          show_column_names = FALSE,
          show_row_names = FALSE, 
          row_split = modeltype, 
          row_title_rot = 0, 
          column_title_gp = gpar(fontsize=medfontsize),
          row_title_gp = gpar(fontsize=medfontsize),
          column_split = sampletype,
          heatmap_legend_param = list(title_position = "lefttop-rot", 
                                      grid_width = unit(2, "mm"), 
                                      legend_height = unit(1, "cm"),
                                      at = c(-2, 2), labels = c("close", "distant")))
)
distHM

# 104, 132, 160, 158
plotseed <-  132
set.seed(plotseed)
pathplot <- emapplot(plotRes,
                     cex_label_category=0.8,
                     cex_category=1, # !RELATIVE! size of the dots
                     legend_n=2, showCategory = 20,
                     group_legend=FALSE) +
  theme(legend.position= c(0.2, 0.12)) #+ 
pathplot

###############################################################################
toprow <- plot_grid(clusterHM, pathplot, nrow = 1, labels = c("B", "C"), rel_widths = c(1,3))
midrow <- plot_grid(markerHM, labels = c("D"))
botrow <- plot_grid(distHM, labels = c("E"))

organoidfig <- plot_grid(toprow, midrow, botrow, nrow = 3, rel_heights = c(3,1, 0.5) )

png("Figures/model_system_figure_8.png", width = 20, height = 23, units = 'cm', res = 300)
grid.arrange(organoidfig) # Make plot
dev.off()

pdf("Figures/model_system_figure.pdf", width = 20/2.54, height = 23/2.54)
grid.arrange(organoidfig) # Make plot
dev.off()
