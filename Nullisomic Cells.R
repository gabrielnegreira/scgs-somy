#inputs
BPK282 <- BPK282_raw_somies[,1:64]
BPK081 <- BPK081_raw_somies[,1:64]
MixPop <- MixPop_raw_somies[,1:64]

BPK282_integers <- BPK282_integer_somies
BPK081_integers <- BPK081_integer_somies
MixPop_integers <- MixPop_integer_somies

BPK282$sample <- "BPK282/0 cl4"
BPK081$sample <- "BPK081/0 cl8"
MixPop$sample <- "Sup-Mosaic"


BPK282[match(BPK282_integers$cell, BPK282$cell_id), 2:37] <- BPK282_integers[,2:37]
BPK081[match(BPK081_integers$cell, BPK081$cell_id), 2:37] <- BPK081_integers[,2:37]
MixPop[match(MixPop_integers$cell, MixPop$cell_id), 2:37] <- MixPop_integers[,2:37]

BPK282$is_integer <- FALSE
BPK081$is_integer <- FALSE
MixPop$is_integer <- FALSE

BPK282$is_integer[match(BPK282_integers$cell, BPK282$cell_id)] <- TRUE
BPK081$is_integer[match(BPK081_integers$cell, BPK081$cell_id)] <- TRUE
MixPop$is_integer[match(MixPop_integers$cell, MixPop$cell_id)] <- TRUE


all_cells <- rbind(BPK282, BPK081, MixPop)

# rm(BPK282)
# rm(BPK081)
# rm(MixPop)

#custom functions
has_nullisomy <- function(x, threshold = 0.2){
  if(length(which(x <= threshold)) > 0){
    x <- TRUE
  }else{
    x <- FALSE
  }
  return(x)
}


#code
all_cells$nullisomic <- apply(all_cells[,2:37], 1, has_nullisomy)
to_plot <- filter(all_cells, nullisomic)
to_plot <- to_plot[,c(2:37, which(colnames(all_cells) == "cell"), which(colnames(all_cells) == "sample"), which(colnames(all_cells) == "nullisomic"), which(colnames(all_cells) == "is_integer"))]
to_plot <- gather(to_plot, key = "chromosome", value = "somy", -cell, -sample, -nullisomic, -is_integer)


min_scale <- 0
max_scale <- 8

cell_levels <- unique(arrange(to_plot, sample)$cell)
sample_levels <- c(unique(BPK282$sample), unique(BPK081$sample), unique(MixPop$sample))
to_plot <- mutate(to_plot, label_color = ifelse(somy == 0, "red", "grey"))



  p1 <- to_plot %>% 
    filter(is_integer) %>%
      ggplot(aes(x = factor(cell, levels = cell_levels), y = factor(chromosome, levels = c(36:1)), fill = somy, label = round(somy, 1)))+
        facet_grid(col = vars(factor(sample, levels = sample_levels)), scales = "free", space = "free")+
        geom_tile(size = 1)+
        geom_text(color = "grey")+
        guides(col = F, fill = F)+
        ggtitle("Cells maintained in analysis (integer somies)")+
        scale_y_discrete(name = "Chromosome")+
        scale_x_discrete(name = "Cell id")+
        scale_fill_gradientn(name = "Somy", colors = heat_col, limits = c(min_scale, max_scale))+
        scale_color_discrete(name = "Sample")+
        theme(axis.text.x = element_text(angle = 90), axis.title = element_text(size = 14, face = 2), legend.title = element_text(size = 14, face = 2))
      

  p2 <- to_plot %>% 
          filter(!is_integer) %>%
          ggplot(aes(x = factor(cell, levels = cell_levels), y = factor(chromosome, levels = c(36:1)), fill = somy, label = round(somy, 1)))+
            facet_grid(col = vars(factor(sample, levels = sample_levels)), scales = "free", space = "free")+   
            geom_tile(size = 1)+
            geom_text(color = "grey")+
            guides(col = F) + 
            ggtitle("Cells removed from analysis (raw somies)")+
            scale_y_discrete(name = "Chromosome")+
            scale_x_discrete(name = "Cell id")+
            scale_fill_gradientn(name = "Somy", colors = heat_col, limits = c(min_scale, max_scale))+
            scale_color_discrete(name = "Sample")+
            theme(axis.text.x = element_text(angle = 90), axis.title = element_text(size = 14, face = 2), legend.title = element_text(size = 14, face = 2))



ratio <- c(length(unique(to_plot$cell[which(to_plot$is_integer)])),length(unique(to_plot$cell[which(!to_plot$is_integer)])))
ratio <- ratio[1]/ratio[2]

if(ratio > 1){
  ratio <- 1/ratio
  ratio <- ratio - 0.07
  plot_grid(p1, p2, ncol = 2, rel_widths = c(ratio, 1-ratio))
}else{
  ratio <- ratio - 0.07
  plot_grid(p1, p2, ncol = 2, rel_widths = c(ratio, 1-ratio))
}

# rm(BPK282_integers)
# rm(BPK081_integers)
# rm(MixPop_integers)
   
   
  
