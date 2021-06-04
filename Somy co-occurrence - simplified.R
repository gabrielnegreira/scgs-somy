#OBS: to get the MixPop clusters, run the script "HU3 in MixPop" first!

#libraries
library(ggplot2)
library(tidyr)
library(grid)
library(gridExtra)
library(corrplot)
library(ggbiplot)
library(ggfortify)
library(cluster)
library(RColorBrewer)
library(plotly)
library(igraph)
library(circlize)
library(ComplexHeatmap)

cells_rare <- filter(all_cells, karyo_occurrence == 1 & baseline_ploidy == 2)
cells_common <- filter(equalmix, karyo_occurrence > 1 & baseline_ploidy == 2)

data <- cells_rare
data <- cells_common[sample(c(1:nrow(cells_common), size = nrow(cells_rare))),] #common karyotypes


#parameters
clust_method <- "ward.D2"
heatmap <- TRUE
pearson_cutoff <- 0.4 #only links with pearson correlation score higher than the cutoff or lower than -cutoff will be displayed in the choord diagram. 
pvalue_cutoff <- 0.05 #the pvalue cutoff for both the choord diagram and the correlation matrix plot. 
input_data <- "strains" #choose "equalmix" to use the equalmix or "strains" to use the BGS data of the 205 strains of L.donovani
round <- FALSE #if TRUE will round somy values. 
update_main_somies <- FALSE #keep it TRUE.
pca_clust_method <- "clara" #choose clara or kmeans to make the clustering in the PCA plot. 
pca_k <- 3 #number of groups in the PCA plot.

#making the equal mix of cells:

#getting BPK282 and BPK081 data
BPK282 <- BPK282_integer_somies
BPK081 <- BPK081_integer_somies

#selecting the columns that are common to all data frames
cnames <- intersect(colnames(BPK282), colnames(MixPop_A))

BPK282 <- BPK282[,which(colnames(BPK282) %in% cnames)]
BPK081 <- BPK081[,which(colnames(BPK081) %in% cnames)]
MixPop_A <- MixPop_A[,which(colnames(MixPop_A) %in% cnames)]
MixPop_B <- MixPop_B[,which(colnames(MixPop_B) %in% cnames)]
MixPop_C <- MixPop_C[,which(colnames(MixPop_C) %in% cnames)]
MixPop_D <- MixPop_D[,which(colnames(MixPop_D) %in% cnames)]

#choosing the number of cells to randomly pick from each data frame
num_cells <- min(c(nrow(BPK282), 
                   nrow(BPK081), 
                   nrow(MixPop_A), 
                   nrow(MixPop_B), 
                   nrow(MixPop_C), 
                   nrow(MixPop_D)))

#adding information 
BPK282$sample <- "BPK282"
BPK081$sample <- "BPK081"
MixPop_A$sample <- "Cluster A"
MixPop_B$sample <- "Cluster B"
MixPop_C$sample <- "Cluster C"
MixPop_D$sample <- "Cluster D"

#randomly selecting cells from each data frame
equalmix <- rbind(BPK282[sample(c(1:nrow(BPK282)), size = num_cells),],
                  BPK081[sample(c(1:nrow(BPK081)), size = num_cells),], 
                  MixPop_A[sample(c(1:nrow(MixPop_A)), size = num_cells),], 
                  MixPop_B[sample(c(1:nrow(MixPop_B)), size = num_cells),], 
                  MixPop_C[sample(c(1:nrow(MixPop_C)), size = num_cells),], 
                  MixPop_D[sample(c(1:nrow(MixPop_D)), size = num_cells),])

#getting the BGS data of the strains

df_204_strains <- read_delim("inputs/eLife_paper_donovani_strains.csv", 
                             ";", escape_double = FALSE, locale = locale(decimal_mark = ",", 
                                                                         grouping_mark = "."), trim_ws = TRUE)

#inputs
if(input_data == "strains"){
  input <- t(df_204_strains[,c(2:206)]) #change
}else{
  if(input_data == "equalmix"){
    input <- equalmix
  }else{
    input <- input_data
  }
}

#code
if(input_data != "strains"){
  cells <- filter(input, mean_somy < 2.5 & mean_somy > 2) #change
  somies <- cells[,2:37] #change
  somies <- somies[,1:36] #change
}else{
  somies <- input
}


somies <- apply(somies, 2, unlist(as.numeric))

if(round == TRUE){
  somies <- round(somies)
}

chromo_var_percent <- c()
for(i in c(1:ncol(somies))){
  chromo_var_percent <- c(chromo_var_percent, length(somies[which(somies[,i] == 2), i])/nrow(somies))
}

colnames(somies) <- c(1:ncol(somies))
total <- nrow(somies)

if(update_main_somies!= FALSE){
    main_somies <- apply(round(somies), 2, Mode)
    main_somies <- as.vector(unlist(main_somies))
}


pearson_matrix <- cor(somies, method = "pearson")
matrix_to_plot <- pearson_matrix

#plotting

#heatmap
if(heatmap == TRUE){
  to_plot <- somies
  to_plot <- to_plot[,1:36]
  to_plot <- numeric.matrix(data.frame(to_plot))
  colnames(to_plot) <- c(1:ncol(to_plot))
  min_scale <- 0
  max_scale <- round(max(to_plot))
  max_scale <- 8
  
  heatmap.2(to_plot, key.xlab = "Somy",
            main = "All samples merged",
            #Rowv = F, 
            #Colv = F,
            trace="none", 
            keysize = 1,
            #key = F,
            density = "density", 
            linecol = 1,
            col= color_palette(200),
            #col = heat_palette(200),
            adjCol = c(0.5,0.5),
            breaks = seq(min_scale,max_scale, length = 201), 
            dendrogram = "both",
            cexCol=2, 
            labRow = FALSE, 
            key.title = "Color Key", 
            hclust=function(x) hclust(x,method = "ward.D2")
  )
  
}

#correlation plot

matrix_to_plot[which(is.na(matrix_to_plot))] <- 0
sig <- cor.mtest(somies, method = "pearson", conf.level = 1-pvalue_cutoff)

to_plot <- matrix_to_plot
diag(to_plot) <- NA

corrplot(to_plot, 
         p.mat = sig$p,
         sig.level = 0.05,
         insig = "blank",
         method = "circle", 
         #addrect = 3,
         na.label = "o",
         #type = "upper", 
         order = "hclust",
         hclust.method = "complete",
         tl.col = "black")

 
#circular relationship correlogram using pearson

ad_matrix <- matrix_to_plot
diag(ad_matrix) <- 0
ad_matrix[which(sig$p > pvalue_cutoff)] <- 0.0008
ad_matrix[which(ad_matrix < pearson_cutoff)] <- 0.0008
#ad_matrix[which(ad_matrix > -pearson_cutoff & ad_matrix < pearson_cutoff)] <- 0.0008



grid_col <- integer_col[main_somies+1]
names(grid_col) <- c(1:36)

intervals <- c(round(min(ad_matrix[which(ad_matrix > 0.0008)]), 1),
               round(mean(ad_matrix[which(ad_matrix > 0.0008)]), 1),
               round(max(ad_matrix[which(ad_matrix > 0.0008)]),1))

col_fun <- colorRamp2(intervals, c("#fcf2f2","#fc0000", "#630000"), transparency = 0.2)

group <- as.character(main_somies)
names(group) <- c(1:36)
group <- sort(group)

gap <- rep(1.6, times = length(group))
gap[c(which(diff(as.numeric(group))!=0), length(group))] <- 20

  
circos.par(start.degree = 90, gap.after = gap)
chordDiagram(ad_matrix, 
             #group = group,
             symmetric = TRUE,
             annotationTrack = c("grid"),
             preAllocateTracks = list(track.height = 0.1),
             scale = FALSE,
             order = names(group), 
             grid.col = grid_col, 
             link.visible = ad_matrix > pearson_cutoff,
             col = col_fun
             )

circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)

legend_links <- Legend(at = round(intervals, 3), col_fun = col_fun, title_position = "topleft", title = "Pearson\nCorrelation")
legend_tracks <- Legend(at = unique(main_somies), type = "grid", legend_gp = gpar(fill = unique(grid_col)), title_position = "topleft", title = "Main\nSomy")
  
  
legend_to_plot <- packLegend(legend_links, legend_tracks, row_gap = unit(10, "mm"))

draw(legend_to_plot, x = unit(0.85, "npc"),  just = c("right", "center"))

#rm(ad_matrix)
rm(legend_to_plot)
rm(legend_links)
rm(legend_tracks)

circos.clear()

#pca
data_to_pca <- t(somies[,-31])
pca <- prcomp(data_to_pca)

if(pca_clust_method == "clara"){
  groups <- clara(data_to_pca, k = pca_k)$clustering
}else{
  groups <- kmeans(data_to_pca, centers = pca_k)$cluster
}

#group_eq <- c(2,1) #this is used to know which group in BGS correspond to the group in scgs pca

pca_out <- as.data.frame(pca$x)
pca_out$lab <- rownames(pca_out)
pca_out$group <- groups
#pca_out <- mutate(pca_out, group = group_eq[group])
#pca_out <- mutate(pca_out, lab = ifelse(group == 1, NA, lab))
pca_importance <- summary(pca)$importance[2,]

pca_out$main_somy <- main_somies[as.numeric(pca_out$lab)]

plot_title <- ifelse(input_data == "strains", paste(nrow(input), "isolates (BGS)"), paste(nrow(input), "cells (SCGS)"))
plot_colors <- c("orange", "blue", "red", "light blue", "purple", "yellow")
plot_colors <- plot_colors[pca_out$group]
names(plot_colors) <- pca_out$group


p1 <- ggplot(pca_out, aes(x = PC1, y = PC2, label = lab, shape = factor(main_somy), color = factor(group)))+
  geom_text(color = "black", hjust = -0.3)+
  geom_point(alpha = 0.5, size = 3)+
  scale_x_continuous(name = paste("PC1 (", pca_importance[[1]]*100, "% explained var.)", sep = ""))+
  scale_y_continuous(name = paste("PC2 (", pca_importance[[2]]*100, "% explained var.)", sep = ""))+
  scale_color_manual(name = "Group", values = plot_colors)+
  scale_shape(name = "Main Somy")+
  labs(title = plot_title)+ #change
  theme(axis.title = element_text(size = 14, face = 2), 
        legend.title = element_text(size = 14, face = 2), 
        plot.title = element_text(size = 16, face = 2, hjust = 0.5),
        legend.text = element_text(size = 12))

p1
