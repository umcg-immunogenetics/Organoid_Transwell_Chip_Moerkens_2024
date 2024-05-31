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

bigfontsize <- 10
medfontsize <- 8
smlfontsize <- 5


geneSymbol <- read.table("Data/features.tsv", sep = "\t")[,1:2]
colnames(geneSymbol) <- c("ensembl", "symbol")

# I checked if the removed genes are needed below, they're not
geneSymbol <- geneSymbol[!duplicated(geneSymbol$symbol),]
rownames(geneSymbol) <- geneSymbol$symbol



#####NEW lists
#Epithelial subtypes, Supplementary Figure
#Epithelial subtypes, Supplementary Figure
Intestinal_epithelial_cell <- c('EPCAM', 'CDX2', 'CDH1', 'TJP1')
Neuronal_cell <- c('TUBB2B', 'CRABP1', 'NCAM1')
Mesenchymal_cell <- c('VIM', 'ACTA2')
Stem_cell <- c('SMOC2', 'LGR5')
Transit_amplifying_cell <- c('MKI67', 'TOP2A', 'CENPF')
Enterocyte  <- c('FABP2', 'RBP2', 'APOA4', 'ALDOB', 'ANPEP', 'SI')
Paneth_cell <- c('LYZ', 'PRSS2', 'PLA2G2A')
Goblet_cell <- c('MUC2', 'ZG16', 'CLCA1')
Enteroendocrine_cell <- c('CHGA', 'NEUROD1')

all_genes <- c(Intestinal_epithelial_cell, 
               Neuronal_cell,
               Mesenchymal_cell,
               Stem_cell, 
               Transit_amplifying_cell,
               
               Enterocyte,
               Paneth_cell, 
               Goblet_cell, 
               Enteroendocrine_cell)

gene_sets <- c(rep('Intestinal\nepithelial cell', length(Intestinal_epithelial_cell)),
               rep('Neuronal cell', length(Neuronal_cell)),
               rep('Mesenchymal cell', length(Mesenchymal_cell)),
               rep('Stem cell', length(Stem_cell)),
               rep('Transit\namplifying cell', length(Transit_amplifying_cell)),
               rep('Enterocyte', length(Enterocyte)),
               
               rep('Paneth cell', length(Paneth_cell)),
               rep('Goblet cell', length(Goblet_cell)),
               rep('Enteroendocrine\ncell', length(Enteroendocrine_cell))
               )
               

gene_table <- data.frame(all_genes, gene_sets)
gene_table <- gene_table[gene_table$all_genes %in% geneSymbol$symbol,]

gene_table$ensemble <- geneSymbol[gene_table$all_genes, 1]
table(duplicated(gene_table$all_genes))
###############################################################################

counts <- read.table("Results/figure_data/supp_fig1_counts.csv")
gene_table = gene_table[gene_table$ensemble %in% rownames(counts),]


# hmcounts <- scale(t(counts[unique(gene_table$ensemble), order(colnames(counts))]))
hmcounts <- counts[gene_table$ensemble, order(colnames(counts))]
hmcounts <- t(scale(t(counts[unique(gene_table$ensemble), order(colnames(counts))])))

rownames(hmcounts) <- gene_table$all_genes

annot_gene <- data.frame(row.names = gene_table$all_genes,
                         process = factor(gene_table$gene_sets, ordered = TRUE, 
                                          levels=c("Intestinal\nepithelial cell", 
                                                   "Neuronal cell", 
                                                   "Mesenchymal cell", 
                                                   "Stem cell", 
                                                   "Transit\namplifying cell",
                                                   "Enterocyte",
                                                   "Paneth cell", 
                                                   "Goblet cell", 
                                                   "Enteroendocrine\ncell")
                         ))

annot_group <- data.frame(row.names = colnames(hmcounts),
                          sampletype = factor(as.factor(c(rep("Chip\nEM-DM+D+P", 3), 
                                                          rep("Organoid\nDM+D+P", 3), 
                                                          rep("Organoid\nEM", 3), 
                                                          rep("Transwell\nDM+D+P", 3),
                                                          rep("Transwell\nEM-DM+D+P", 3),
                                                          rep("Transwell\nEM", 3))),
                                              ordered = TRUE, levels=c("Organoid\nEM", "Organoid\nDM+D+P",
                                                                       "Transwell\nEM", "Transwell\nEM-DM+D+P", "Transwell\nDM+D+P",
                                                                       "Chip\nEM-DM+D+P")))
scale_col_fun = colorRamp2(c(0, 5, 20), c("yellow", "orange", 'red'))
markerHM <- as.ggplot(
  Heatmap(hmcounts, name="expression", cluster_rows = FALSE, cluster_columns = FALSE, 
          column_names_gp = gpar(fontsize = medfontsize),
          col = scale_col_fun,
          row_names_gp = gpar(fontsize = medfontsize),
          show_row_names = TRUE, 
          row_split = annot_gene$process, 
          column_split = annot_group$sampletype, 
          column_title_gp = gpar(fontsize=medfontsize),
          row_title_gp = gpar(fontsize=medfontsize),
          row_title_rot = 0, 
          column_names_rot = 90,
          column_title_rot = 90,
          heatmap_legend_param = list(title_position = "lefttop-rot", grid_width = unit(2, "mm")))
)
markerHM

markerHM <- as.ggplot(
  Heatmap(hmcounts, name="expression", cluster_rows = FALSE, cluster_columns = FALSE, 
          column_names_gp = gpar(fontsize = medfontsize),
#          col = scale_col_fun,
          row_names_gp = gpar(fontsize = medfontsize),
          show_row_names = TRUE, 
          row_split = annot_gene$process, 
          column_split = annot_group$sampletype, 
          column_title_gp = gpar(fontsize=medfontsize),
          row_title_gp = gpar(fontsize=medfontsize),
          row_title_rot = 0, 
          column_names_rot = 90,
          column_title_rot = 90,
          heatmap_legend_param = list(title_position = "lefttop-rot", grid_width = unit(2, "mm")))
)
markerHM
