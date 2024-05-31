
###############################################################################
# Generate the organoid figure for Renee's manuscript
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

plotseed <-  8   # emapplot uses some randomization so if you want the same figure
bigfontsize <- 10
medfontsize <- 8
smlfontsize <- 5
scale_col_fun = colorRamp2(c(-2, 0, 2), c("yellow", "blue", 'darkblue'))
# scale_col_fun = colorRamp2(c(0, 2.5, 5), c("yellow", "blue", 'darkblue'))

###############################################################################

distanceCounts <- read.table("Results/figure_data/distance_matrix_counts.csv")
sample_distances <- as.matrix(dist(t(distanceCounts)))[1:21, 31:84]
sample_distances <- sample_distances[grep("C_EM_EM", rownames(sample_distances), invert = TRUE),]

sample_distances <- scale(sample_distances)
sample_distances <- sample_distances[order(rownames(sample_distances)),]


rownames(sample_distances)
modeltype <- factor(as.factor(c(rep("Chip EM-DM+D+P", 3), 
                                rep("Organoid DM+D+P", 3),
                                rep("Organoid EM", 3),
                                rep("Transwell DM+D+P", 3),
                                rep("Transwell EM-DM+D+P", 3),
                                rep("Transwell EM", 3))),
                    ordered = TRUE, levels=c("Organoid EM", "Organoid DM+D+P",
                                             "Transwell EM", "Transwell EM-DM+D+P", "Transwell DM+D+P", 
                                             "Chip EM-DM+D+P"))

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


###############################################################################

distanceCounts <- read.table("Results/figure_data/distance_matrix_counts.csv")
sample_distances <- as.matrix(dist(t(distanceCounts)))[1:21, 31:84]
sample_distances <- sample_distances[grep("C_", rownames(sample_distances), invert = TRUE),]

sample_distances <- scale(sample_distances)
sample_distances <- sample_distances[order(rownames(sample_distances)),]


rownames(sample_distances)
modeltype <- factor(as.factor(c(rep("Organoid DM+D+P", 3),
                                rep("Organoid EM", 3),
                                rep("Transwell DM+D+P", 3),
                                rep("Transwell EM-DM+D+P", 3),
                                rep("Transwell EM", 3))),
                    ordered = TRUE, levels=c("Organoid EM", "Organoid DM+D+P",
                                             "Transwell EM", "Transwell EM-DM+D+P", "Transwell DM+D+P"))

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
