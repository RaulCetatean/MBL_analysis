library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(dendextend)
library(dplyr)

setwd("/Users/Raul/Desktop/PhD Docs/Publication_works/MBL/")
amr_data <-read.csv("AMR_MBL_poster.csv", header=TRUE, encoding = "UTF-8",check.names=FALSE)

amr <- amr_data[,-1]
rownames(amr) <- amr_data[,1]
amr_mat <- as.matrix(amr)
set.seed(42)

# Separate data in groups
beta_lactam <- amr[,1:9]
porins <- amr[,10:11]
virulence <- amr[,12:17]
plasmids <- amr[,18:30]
iron <- amr[,31:42]

# Create heatmap for each group

#Beta-lactams
col_1 = colorRamp2(c(0, 1), c("white", "#db4260"))
ht1 <- Heatmap(beta_lactam, column_names_side = "top",
               rect_gp = gpar(col = "black", lwd = 2),
               cluster_rows = FALSE,
               show_column_dend = FALSE,
               cluster_columns = FALSE,
               show_row_names = FALSE,
               col=col_1,
               column_title ="Beta-lactams",
               #column_title_rot = 45, 
               column_title_gp = gpar(fontsize = 18, fontface = "bold"),
               column_names_gp = gpar(fontsize = 16),
               height = unit(200, "mm"))

# Porins
col_2 =colorRamp2(c(0,1),c("white","#1E88E5"))
ht2 <- Heatmap(porins, column_names_side = "top",
               rect_gp = gpar(col = "black", lwd = 2),
               cluster_rows = FALSE,
               show_column_dend = FALSE,
               cluster_columns = FALSE,
               show_row_names = FALSE,
               col=col_2,
               column_title ="Porins",
               #column_title_rot = 45,
               column_title_gp = gpar(fontsize = 18, fontface = "bold"),
               column_names_gp = gpar(fontsize = 16),
               height = unit(200, "mm"))

col_3 <- colorRamp2(c(0,1),c("white","#58cfab"))
ht3 <- Heatmap(virulence, column_names_side = "top",
               rect_gp = gpar(col = "black", lwd = 2),
               cluster_rows = FALSE,
               show_column_dend = FALSE,
               cluster_columns = FALSE,
               show_row_names = FALSE,
               col=col_3,
               column_title ="Virulence",
               #column_title_rot = 45,
               column_title_gp = gpar(fontsize = 18, fontface = "bold"),
               column_names_gp = gpar(fontsize = 16),
               height = unit(200, "mm"))

col_4 =colorRamp2(c(0,1),c("white","#FFC107"))
ht4 <- Heatmap(plasmids, column_names_side = "top",
               rect_gp = gpar(col = "black", lwd = 2),
               cluster_rows = FALSE,
               show_column_dend = FALSE,
               cluster_columns = FALSE,
               show_row_names = FALSE,
               col=col_4,
               column_title ="Plasmids",
               #column_title_rot = 45,
               column_title_gp = gpar(fontsize = 18, fontface = "bold"),
               column_names_gp = gpar(fontsize = 16),
               height = unit(200, "mm"))

col_5 =colorRamp2(c(0,1),c("white","#66615C"))
ht5 <- Heatmap(iron, column_names_side = "top",
               rect_gp = gpar(col = "black", lwd = 2),
               cluster_rows = FALSE,
               show_column_dend = FALSE,
               cluster_columns = FALSE,
               show_row_names = FALSE,
               col=col_5,
               column_title ="Iron Uptake/Transport",
               #column_title_rot = 45,
               column_title_gp = gpar(fontsize = 18, fontface = "bold"),
               column_names_gp = gpar(fontsize = 16),
               height = unit(200, "mm"))

ht_list = ht1 + ht2 + ht3 + ht4 + ht5

ht <- draw(ht_list, ht_gap = unit(5, "mm"), show_heatmap_legend = FALSE)
# I usually save it with height 1125, trying width 950 Rplot_8
# for rplot_9 ill try with width 900 not bad
# for Rplot10 width 850px