#inputs
sample_name <- "BPK081"
reads_raw <- "BPK081 - reads_raw.csv"
cell_metrics <- "BPK081 - per_cell_summary_metrics.csv" #this is used to know how many cells are present in the sample

#parameters
plot_cells <- F #Will make a plot with several information for each cell. Takes about 1 second per cell, so it is quite time-consuming. 
chromo_var_threshold <- 1.7 #I used 2 for BPK282 and 1.7 for the others (BPK282 had lower depth than other samples and thus, was more noisy)
scale_factor_min <- 1.8 #determines the minimum scale_factor value
scale_factor_max <- 5 #determines the maximum scale_factor value
scale_factor_filter_min <- 1.8 #will remove cells resolved with scale factor lower than this number
scale_factor_filter_max <- 5 #will remove cells resolved with scale factor higher than this number
polyploids_var_threshold <- 1.7 #polyploid cells will only be considered true polyploid if chromo_var is lower than this number
polyploids_distance_to_integers_limit <- 0.4 #polyploid cells will only be considered true polyploid if raw somies intermediate values are at distances lower than this to their closest integer. 
chromo_var_method <- "max_5" #will determine how intrachromosomal variation will be calculated. "average" or "mean" will use the average between the variation of all chromosomes in a cell. "max_N" will use the average between the N chromosomes with the highest variation.
chromo_var_segments <- 3 #to calculate the intrachromosomal variation, the code divide each chromosome into segments. This parameter determines how many segments each chromosome will be divided into. 
karyo_top_number <- "all" #number of karyotypes to display in the karyo_top_plot. Use "all" to display all karyotypes or use a numeric value to limit the number of karyotypes to display. 
color_option <- 1 #two options for coloring the heatmaps. 1 uses a "rainbow"-like pallet, while 2 uses a "brewer"-like pallet. In 1, color differences are more pronounced but 2 works better if there are somies higher than 6. 
gc_correction <- FALSE #if TRUE will correct depth by GC content. PS: Not implemented yet. Keep it FALSE. 
depth_threshold <- c(0, 2000) #cells with depth lower than i1 or higher than i2 are removed (use it to try to remove doublets)
#############################################START#######################################################################
#NTS: Add comments to each step. Solve issues with limited number of color in color palettes vectors. 
#code 
library(readr)
library(gridExtra)
#cleaning environment
raw_somies <- NULL
integer_somies <- NULL
karyotypes <- NULL
number_of_cells <- NULL
y_max <- NULL
rnames <- NULL


#getting inputs
reads_raw <- read_csv(paste(getwd(),"/inputs/", reads_raw, sep = ""))
cell_metrics <- read_csv(paste(getwd(),"/inputs/", cell_metrics, sep = ""))
bins_gc <- read_csv("inputs/bins_gc_content.csv")
number_of_cells <- nrow(cell_metrics)
reads_raw <- reads_raw[1:number_of_cells,]

#creating working directory
original_wd <- getwd()
current_wd <- getwd()
dir.create("outputs")
sample_folder_path <- paste(current_wd, "/outputs/", sample_name, sep = "")
dir.create(sample_folder_path)



#gc correction: OBS not working yet. 
if(gc_correction == TRUE){
  bins_gc$index <- c(1:nrow(bins_gc))-1
  bins_gc$norm_depth_cell <- colMeans(scgs_normalize_reads(reads_raw)[,-1])
  bins_gc$norm_depth_chromo <- colMeans(scgs_normalize_reads(reads_raw, "chromosome")[,-1])
  bins_gc_calc <- filter(bins_gc,
                         norm_depth_cell != 0
                         & norm_depth_cell > quantile(bins_gc$norm_depth_cell, c(0.05, 0.5, 0.99))[[1]] 
                         & norm_depth_cell <= quantile(bins_gc$norm_depth_cell, c(0.05, 0.5, 0.99))[[3]]
                         & gc >= quantile(bins_gc$gc, c(0.05, 0.5, 0.95))[[1]] 
                         & gc <= quantile(bins_gc$gc, c(0.05, 0.5, 0.95))[[3]]
                         & gc != 0
                         )

  gc_loess <- loess(norm_depth_chromo ~ gc, data = bins_gc_calc, span = 0.3)
  bins_gc_calc$correction_factor <- predict(gc_loess, bins_gc_calc$gc)
  bins_gc_calc$corrected_depth <- bins_gc_calc$norm_depth_cell * (1/bins_gc_calc$correction_factor)

  gc_function <- lm(bins_gc_calc$correction_factor ~ bins_gc_calc$gc)

  j <- order(bins_gc_calc$gc)

  plot(bins_gc_calc$gc, bins_gc_calc$norm_depth_cell)
  lines(bins_gc_calc$gc[j], gc_loess$fitted[j], col = "red", lwd = 3)

  p1 <- ggplot(bins_gc_calc, aes(x = index, y = norm_depth_cell))+
    geom_point(size = 1, alpha = 0.5, aes(color = chrom))+
    #geom_boxplot(aes(group = chrom))+
    #geom_smooth(method = "loess")+
    geom_smooth(method = "lm", color = "red")+
    scale_color_manual(values = rep(c("orange", "black"), times = 100))+
    scale_y_continuous(limits = c(0, 3))

  p2 <- ggplot(bins_gc_calc, aes(x = index, y = corrected_depth))+
    geom_point(size = 1, alpha = 0.5, aes(color = chrom))+
   #geom_smooth(method = "loess")+
    geom_smooth(method = "lm", color = "red")+
    scale_color_manual(values = rep(c("orange", "black"), times = 100))+
    scale_y_continuous(limits = c(0, 3))

  plot_grid(p1, p2, ncol = 1)

  ggplot(bins_gc_calc, aes(x = gc, y = norm_depth_cell))+
   geom_point()+
   geom_smooth(method = "loess")

  ggplot(bins_gc_calc, aes(x = gc, y = 1/correction_factor))+
   geom_point()+
    geom_smooth()
}

#raw somy calculation
raw_somies <- scgs_calc_somy(reads_raw, 
                            scale_factor_min = scale_factor_min, 
                            scale_factor_max = scale_factor_max, 
                            method = "mean", 
                            chromo_var_method = chromo_var_method, 
                            chromo_var_segments = chromo_var_segments,
                            plot_path = sample_folder_path)

raw_values <- round(raw_somies[,2:37])
raw_somies$baseline_ploidy <- apply(raw_values, 1, Mode)
rm(raw_values)
raw_somies <- cbind(raw_somies, cell_metrics, sample_name)


integer_somies <- raw_somies

#filtering cells to establish integer_somies
integer_somies <- integer_somies %>%
  filter(scale_factor >= scale_factor_filter_min & scale_factor <= scale_factor_filter_max) %>% #filter by scale_factor
  filter(chromo_var <= chromo_var_threshold) %>% #filter by chromo_var (ICV score)
  filter(effective_depth_of_coverage >= depth_threshold[1] & effective_depth_of_coverage <= depth_threshold[2]) #filter by depth


#determining the most common ploidy
ploidies <- as.numeric(names(sort(table(integer_somies$baseline_ploidy), decreasing = TRUE)))
main_ploidy <- ploidies[1]

#filtering the polyploid cells which doesn't meet the criteria. 
integer_somies <- integer_somies %>%
                    mutate(to_remove = ifelse(baseline_ploidy == 2, FALSE,
                           ifelse(chromo_var >= polyploids_var_threshold, TRUE,
                                  ifelse(scaled_furthest_distance_to_integers >= polyploids_distance_to_integers_limit, TRUE, FALSE))))

integer_somies <- filter(integer_somies, !to_remove)

#formating the integer_somies data frame
integer_somies$to_remove <- NULL
integer_somies$diff <- NA
integer_somies$karyotype <- NA
integer_somies$karyo_occurrence <- NA
#converting raw somies into integers
GMM_path <- paste(sample_folder_path, "GMMs", sep = "/")
dir.create(GMM_path)
for(i in ploidies){
  if(nrow(filter(integer_somies, baseline_ploidy == i)) == 0){
    next
  }
  integers <- filter(integer_somies, baseline_ploidy == i)
  if(i == main_ploidy){
    print(paste("establishing integer somies of", nrow(integers), "cells with baseline ploidy", i, "using GMMs"))
  path = paste(GMM_path, paste("baseline_ploidy", i, "- GMM"), sep = "/")
  dir.create(path)
  integers <- scgs_solve_somy(integers, integer_method = "GMM", sample_name = sample_name, plot_path = path)
  }else{
    print(paste("establishing integer somies of", nrow(integers), "cells with baseline ploidy", i, "by rounding to closest integer"))
    path = paste(GMM_path, paste("baseline_ploidy", i, "- rounding"), sep = "/")
    dir.create(path)
    integers <- scgs_solve_somy(integers, integer_method = "round", sample_name = sample_name, plot_path = path)
  }
  integer_somies[match(integers$cell, integer_somies$cell),] <- integers
}

raw_somies <- mutate(raw_somies, removed = ifelse(cell %in% integer_somies$cell, FALSE, TRUE))

assign(paste(sample_name,"_raw_somies", sep = ""), raw_somies)
assign(paste(sample_name,"_integer_somies", sep = ""), integer_somies)

#creating the karyotypes list

karyotypes <- integer_somies
karyotypes <- karyo_uniques(karyotypes[,2:37], add_node_id = FALSE, decluster = FALSE)

f <- function(x){ #this function will estimate the baseline ploidy of a karyotype
  x <- x[[1]]
  x <- as.character(x)
  x <- strsplit(x, split = "_")
  x <- unlist(x)
  x <- Mode(x)
  x <- as.numeric(x)
  return(x)
}

karyotypes$baseline_ploidy <- apply(karyotypes, 1, f)
rm(f)

assign(paste(sample_name,"_karyotypes", sep = ""), karyotypes)
assign(paste(sample_name, "_metrics", sep = ""), cell_metrics)


write.csv(reads_raw, paste(sample_folder_path,"/", sample_name, "_reads_raw.csv", sep = ""))
write.csv(raw_somies, paste(sample_folder_path,"/", sample_name, "_raw_somies.csv", sep = ""))
write.csv(integer_somies, paste(sample_folder_path,"/", sample_name, "_integer_somies.csv", sep = ""))
write.csv(cell_metrics, paste(sample_folder_path,"/", sample_name, "_metrics.csv", sep = ""))
write.csv(karyotypes, paste(sample_folder_path,"/", sample_name, "_karyotypes.csv", sep = ""))

print(paste(nrow(raw_somies) - nrow(integer_somies), "cells removed"))
print(paste(length(unique(integer_somies$karyotype)), "karyotypes"))


#plots

#1
print("plotting chromo_var")

to_plot <- raw_somies
pdf(paste(sample_folder_path, "/", sample_name,"_chromo_var_plot",".pdf", sep = ""), width = 22, height = 8.5)
p1 <- ggplot(to_plot, aes(x = chromo_var))+
  geom_density()+
  geom_vline(xintercept = chromo_var_threshold, color = "red")+
  scale_x_continuous(limits = c(1,5), breaks = c(1:5))+
  labs(title = paste(sample_name,
                     ": intrachromosomal variation calculated with method ",
                     chromo_var_method,
                     " and segmentation = ",
                     chromo_var_segments,
                     "; ",
                     length(raw_somies$chromo_var[which(raw_somies$chromo_var > 5)]), 
                     " values out of scale not shown.", sep = ""))

to_plot <- filter(raw_somies, effective_depth_of_coverage >= depth_threshold[1] & effective_depth_of_coverage <= depth_threshold[2])
p2 <- ggplot(to_plot, aes(x = chromo_var))+
  geom_density()+
  geom_vline(xintercept = chromo_var_threshold, color = "red")+
  scale_x_continuous(limits = c(1,5), breaks = c(1:5))+
  labs(title = paste("same as before cells with depth lower than ", depth_threshold[1], " or higher than ", depth_threshold[2]," were removed. ",
                     length(raw_somies$chromo_var[which(raw_somies$chromo_var > 5)]), 
                     " values out of scale not shown.", sep = ""))
plot_grid(p1,p2, ncol = 2)  
  
dev.off()



print("plotting baseline somy vs depth")
to_plot <- raw_somies
baseline_depth_plot <-ggplot(to_plot, aes(x = factor(baseline_ploidy), y = effective_depth_of_coverage, color = chromo_var))+
  geom_boxplot()+
  geom_jitter(width = 0.1, height = 0, size = 0.7)+
  scale_x_discrete(name = "Baseline somy", breaks = c(1:1000))+
  scale_color_viridis_c(option = "A")+
  labs(title = paste(sample_name, "- All cells - scale_factor set to be between", scale_factor_min, "and", scale_factor_max))

assign(paste(sample_name,"_baseline_depth_plot", sep = ""), baseline_depth_plot)
pdf(paste(sample_folder_path, "/", sample_name,"_baseline_depth_plot",".pdf", sep = ""), width = 11.69, height = 8.27)
print(baseline_depth_plot)
dev.off()


#2
print("plotting baseline somy vs depth after removing cells wiht high chromo_var")

to_plot <- filter(raw_somies, chromo_var < chromo_var_threshold)
baseline_depth_plot_filtered <-ggplot(to_plot, aes(x = factor(baseline_ploidy), y = effective_depth_of_coverage, color = chromo_var))+
  geom_boxplot()+
  geom_jitter(width = 0.1, height = 0, size = 0.7)+
  scale_x_discrete(name = "Baseline somy")+
  scale_color_viridis_c(option = "A")+
  labs(title = paste(sample_name, "- All cells - scale_factor set to be between", scale_factor_min,"and", scale_factor_max, "cells with chromo_var higher than", chromo_var_threshold, "removed"))

assign(paste(sample_name,"_baseline_depth_plot_filtered", sep = ""), baseline_depth_plot_filtered)
pdf(paste(sample_folder_path, "/", sample_name,"_baseline_depth_plot_filtered",".pdf", sep = ""), width = 11.69, height = 8.27)
print(baseline_depth_plot_filtered)
dev.off()

#3
print("plotting karyotype frequencies")
plot_colors <- c("#4287f5", "#666666", "#f2d027", "orange", "red")
to_plot <- karyotypes

plot_colors <- plot_colors[sort(unique(to_plot$baseline_ploidy))]
karyo_freq_plot <- ggplot(to_plot, aes(x = position, y = frequency, fill = factor(baseline_ploidy)))+
  geom_col()+
  geom_line(color = "red", data = to_plot, aes(y = cummulative/(1/max(to_plot$frequency))))+
  labs(title = paste(sample_name, "- number of identified karyotypes:", nrow(karyotypes)))+
  scale_fill_manual(name = "baseline ploidy", values = plot_colors)+
  scale_y_continuous(name = "Karyotype Frequency", breaks = seq(0, round(max(to_plot$frequency), 1), length = ((round(max(to_plot$frequency), 1) * 100)/5)+1), labels = scales::percent, sec.axis = sec_axis(~. *(1/max(to_plot$frequency)), name = "Cummulative Frequency", labels = scales::percent, breaks = seq(0,1, length = 11)))+
  scale_x_continuous(name = "Karyotype", breaks = fix_seq(1, 10, max(to_plot$position)/10))+
  theme_linedraw()+
  theme(axis.title = element_text(size = 14, face = "bold"), axis.text = element_text(size = 7), panel.grid = element_line(color = "light grey"))

assign(paste(sample_name,"_karyo_freq_plot", sep = ""), karyo_freq_plot)
pdf(paste(sample_folder_path, "/", sample_name,"_karyo_freq_plot",".pdf", sep = ""), width = 11.69, height = 8.27)
print(karyo_freq_plot)
dev.off()

rm(plot_colors)

#4
print("ploting heatmap with unfiltered raw somies")
to_plot <- raw_somies
to_plot <- to_plot[,-1]
to_plot <- to_plot[,1:36]
to_plot <- numeric.matrix(to_plot)
to_plot <- t(to_plot)
min_scale <- 0
max_scale <- round(max(to_plot))
max_scale <- 8


if(color_option == 1){
  integer_col <- c("#001221","#002342", "#014175","#00c3ff", "#33ff00", "#fffa00", "#ffa600", "#D73027", "#A50026", "#541b1b", "#4d0600")
  heat_col <- c("#001221", "#002342", "#002342", "#014175", "#035ba3", "#00c3ff", "#00ffee", "#33ff00", "#ccff00", "#fffa00","#ffa600", "#D73027", "#A50026", "#541b1b", "#4d0600")
}else{
  integer_col <- c("#001221", "#4575B4", "#a4cbe0", "#E0F3F8", "#FFFFBF", "#FDAE61", "#F46D43", "#D73027", "#A50026", "#541b1b")
  heat_col <- c("#313695", "#4575B4", "#74ADD1", "#ABD9E9", "#E0F3F8", "#FFFFBF", "#FEE090", "#FDAE61", "#F46D43", "#D73027", "#A50026", "#541b1b")
}

color_palette <- colorRampPalette(heat_col)

par(bg = "white")
#pdf(paste(sample_folder_path, "/", sample_name,"_heatmap_unfiltered.pdf", sep = ""), width = 11.69, height = 8.27)
png(paste(sample_folder_path, "/", sample_name,"_heatmap_unfiltered.png", sep = ""), width = 1527, height = 1080)
heatmap.2(to_plot, 
          key.xlab = "Somy",
          main = paste(sample_name, " - scale_factor lower than", scale_factor_filter_min, "or higher than", scale_factor_filter_max, "and chromo_var higher than", chromo_var_threshold,"removed"), 
          density = "density", 
          denscol = "black",
          keysize = 1,
          key.title = "Color Key", 
          xlab = paste(ncol(to_plot), "cells"),
          Rowv = F, 
          #Colv = F,
          trace="none", 
          col= color_palette(200),
          #adjCol = c(0.5,0.5),
          breaks = seq(min_scale,max_scale, length = 201), 
          dendrogram = "both",
          cexRow=2, 
          labCol = FALSE, 
          hclust = function(x) hclust(x, method = "ward.D2")
)

dev.off()

#5
print("ploting heatmap with filtered raw somies")
to_plot <- raw_somies[match(integer_somies$cell, raw_somies$cell),]
to_plot <- to_plot[,-1]
to_plot <- to_plot[,1:36]
to_plot <- numeric.matrix(to_plot)
to_plot <- t(to_plot)
min_scale <- 0
max_scale <- round(max(to_plot))
max_scale <- 8

par(bg = "white")
#pdf(paste(sample_folder_path, "/", sample_name,"_heatmap_filtered.pdf", sep = ""), width = 11.69, height = 8.27)
png(paste(sample_folder_path, "/", sample_name,"_heatmap_filtered.png", sep = ""), width = 1527, height = 1080)
par(cex.main = 0.6)
heatmap.2(to_plot, 
          key.xlab = "Somy",
          main = paste(sample_name, " - scale_factor lower than", scale_factor_filter_min, "or higher than", scale_factor_filter_max, ", chromo_var higher than", chromo_var_threshold, "depth of coverage lower than", depth_threshold[1],"or higher than", depth_threshold[2], "removed"), 
          density = "density", 
          denscol = "black",
          keysize = 1,
          key.title = "Color Key", 
          xlab = paste(ncol(to_plot), "cells"),
          Rowv = F, 
          #Colv = F,
          trace="none", 
          col= color_palette(200),
          #adjCol = c(0.5,0.5),
          breaks = seq(min_scale,max_scale, length = 201), 
          dendrogram = "both",
          cexRow=2, 
          labCol = FALSE, 
          hclust = function(x) hclust(x, method = "ward.D2")
)

dev.off()

#6
print("ploting heatmap with integer somies")
to_plot <- integer_somies
to_plot <- to_plot[,-1]
to_plot <- to_plot[,1:36]
to_plot <- numeric.matrix(to_plot)
to_plot <- t(to_plot)
min_scale <- 0
max_scale <- round(max(to_plot))
max_scale <- 8

par(bg = "white")
#pdf(paste(sample_folder_path, "/", sample_name,"_heatmap_integers.pdf", sep = ""), width = 11.69, height = 8.27)
png(paste(sample_folder_path, "/", sample_name,"_heatmap_integers.png", sep = ""), width = 1527, height = 1080)
heatmap.2(to_plot, 
          key.xlab = "Somy",
          density = "density", 
          denscol = "black",
          keysize = 1,
          key.title = "Color Key", 
          main = paste(sample_name, " - integer somies"),
          xlab = paste(ncol(to_plot), "cells"),
          Rowv = F, 
          #Colv = F,
          trace="none", 
          col= color_palette(200),
          #adjCol = c(0.5,0.5),
          breaks = seq(min_scale,max_scale, length = 201), 
          dendrogram = "both",
          cexRow=2, 
          labCol = FALSE, 
          hclust = function(x) hclust(x, method = "ward.D2")
)

dev.off()


#7
to_plot <- karyo_gather(integer_somies[,2:37], decluster = FALSE)
number <- karyo_top_number

if(tolower(karyo_top_number) == "all"){
  number <- nrow(karyotypes)
}

to_plot <- to_plot[1:(number*36),]
karyo_label_cells <- c()
karyo_label_percent <- c()
for(i in c(1:nrow(karyotypes))){
  karyo_label_cells <- c(karyo_label_cells, paste(i,  paste("(",karyotypes[i,2]," cells)", sep = "")))
  
}

legend_size <- ifelse(number > 20, 10, 16)
legend_face <- ifelse(number <= 20, "plain", "bold")

somy_col <- c(integer_col, rep(integer_col[length(integer_col)], times = 100))
karyo_top_plot <- ggplot(to_plot, aes(x = factor(karyo_id), y = as.numeric(Chromosome), label = diff_from_the_first, fill = as.factor(Somy), color = 0)) +
  geom_tile(width = 0.98) +
  geom_text(y = 0.5, size = 4) +
  #geom_text(data = to_plot, aes(label = Marker), color = "grey")+
  scale_fill_manual(values = somy_col[c(min(to_plot$Somy, na.rm = TRUE):max(to_plot$Somy, na.rm = TRUE))+1],
                    breaks = c(min(to_plot$Somy, na.rm = TRUE):max(to_plot$Somy, na.rm = TRUE))) +
  scale_y_continuous(trans = "reverse", breaks = c(1:36)) +
  scale_x_discrete(label = karyo_label_cells) +
  guides(col = FALSE) +
  labs(title = paste(sample_name, " - Ploidy Profile of the", number,"Most Abundant Karyotypes"), fill = "Somy", x = "Karyotype", y = "Chromosome")+
  theme_bw()+
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5, margin = margin(0,0,20,0)), 
        legend.title = element_text(size = 16, face = "bold"), 
        legend.text = element_text(size = 16), 
        axis.text.x = element_text(size = legend_size, face = legend_face, angle = 90, vjust = 0.5) , 
        axis.title.x = element_text(size = 16, face = "bold", margin = margin(20,0,0,0)), 
        panel.grid = element_blank(), 
        axis.title.y = element_text(size=16, face = "bold", margin = margin(0,20,0,0)))

assign(paste(sample_name,"_karyo_top_plot", sep = ""), karyo_top_plot)
pdf(paste(sample_folder_path, "/", sample_name,"_karyo_top_plot",".pdf", sep = ""), width = 11.69, height = 8.27)
print(karyo_top_plot)
dev.off()


#8
to_plot <- karyo_calc_prop(integer_somies[,2:37], TRUE, FALSE)
percent_scale <- seq(0, 1, length = 11)
chromo_prop_plot <- ggplot(to_plot, aes(x = reorder(Chromosome, stability_position), label = Somy, group = Order, y = Proportion, fill = Somy))+
  geom_col()+
  labs(title = sample_name)+
  scale_fill_manual(values = bar_plot_colors[c(min(to_plot$Somy, na.rm = TRUE):max(to_plot$Somy, na.rm = TRUE))+1])+
  scale_y_continuous(breaks = percent_scale, labels = scales::percent, name = "Proportion of Cells")+
  scale_x_discrete(name = "Chromosome")+
  theme_minimal()

assign(paste(sample_name,"_chromo_prop_plot", sep = ""), chromo_prop_plot)
pdf(paste(sample_folder_path, "/", sample_name,"_chromo_prop_plot",".pdf", sep = ""), width = 11.69, height = 8.27)
print(chromo_prop_plot)
dev.off()



#9

reads <- reads_raw
before <- ncol(reads)-1
reads2 <- scgs_remove_empty_bins(reads)
reads2 <- scgs_remove_dev_bins(reads2)
after <- ncol(reads2)-1
bins_removed <- before - after

to_plot<- scgs_gather_bin_depth(reads)
y_max <- max(to_plot$reads_per_bin)

p1 <- ggplot(to_plot, aes(x = bin_num, y = reads_per_bin, color = factor(chromosome)))+
  geom_point(size = 0.1)+
  guides(col = FALSE)+
  scale_color_manual(values = rep(manhattan_plot_colors[c(1,2)], times = 36))+
  labs(title = paste(sample_name, "- Before empty and outlier bins removal"))+
  scale_y_continuous(limits = c(0,y_max))+
  scale_x_continuous(name = "Chromosome", breaks = unique(to_plot$chromo_center), labels = unique(to_plot$chromosome))+
  theme(axis.title = element_text(size = 26, face = "bold"), axis.text = element_text(size = 22))

to_plot <- scgs_gather_bin_depth(reads2)


p2 <- ggplot(to_plot, aes(x = bin_num, y = reads_per_bin, color = factor(chromosome)))+
  geom_point(size = 0.1)+
  guides(col = FALSE)+
  scale_color_manual(values = rep(manhattan_plot_colors[c(1,2)], times = 36))+
  labs(title = paste(sample_name, "- After empty and outlier bins removal:", bins_removed, "bins removed of a total of", before))+
  scale_y_continuous(limits = c(0,y_max))+
  scale_x_continuous(name = "Chromosome", breaks = unique(to_plot$chromo_center), labels = unique(to_plot$chromosome))+
  theme(axis.title = element_text(size = 26, face = "bold"), axis.text = element_text(size = 22))

p3 <- plot_grid(p1, p2, nrow = 2, align = "v", rel_heights = c(0.5, 0.5))

assign(paste(sample_name,"_bin_removal", sep = ""), p3)
png(paste(sample_folder_path, "/", sample_name,"_bin_removal",".png", sep = ""), width = 1524, height = 1080)
print(p3)
dev.off()

#10

to_plot <- raw_somies
cells_with_integer <- integer_somies$cell
to_plot <- to_plot[which(to_plot$cell %in% cells_with_integer),]
relevant_cols <- c("baseline_ploidy", "cell", "barcode")
to_plot <- cbind(to_plot[,2:37], to_plot[,which(colnames(to_plot) %in% relevant_cols)])
to_plot <- filter(to_plot, baseline_ploidy != main_ploidy)

hclust_order <- to_plot$barcode[hclust(dist(to_plot[,2:37]), method = "complete")$order]
to_plot <- gather(to_plot, key = "chromosome", value = "somy", -baseline_ploidy, -cell, -barcode)

min_scale <- 0
max_scale <- ifelse(round(max(to_plot$somy)) > 10, round(max(to_plot$somy)), 8)


colors <- heat_col
colors <- c(colors, rep(colors[length(colors)], times = 100))

  if(max_scale == 8){
    colors <- colors[c(0:(max_scale+6))+1]
  }else{
    colors <- colors[c(0:(max_scale+8))+1]
}



pdf(paste(sample_folder_path, "/", sample_name,"_Potential polyploid cells.pdf", sep = ""), width = length(unique(to_plot$cell)), height = 8.27)

ggplot(to_plot, aes(x = factor(barcode, levels = hclust_order), y = factor(chromosome, levels = c(36:1)), fill = somy, label = round(somy, 1)))+
  facet_grid(col = vars(baseline_ploidy), scales = "free", space = "free")+
  geom_tile(size = 1)+
  geom_text(color = "grey")+
  #guides(col = F, fill = F)+
  ggtitle("Potential Polyploid Cells")+
  scale_y_discrete(name = "Chromosome")+
  scale_x_discrete(name = "Cell Barcode")+
  scale_fill_gradientn(name = "Somy", colors = colors, limits = c(min_scale, max_scale))+
  scale_color_discrete(name = "Sample")+
  theme(axis.text.x = element_text(angle = 90), axis.title = element_text(size = 10, face = 2), legend.title = element_text(size = 14, face = 2))


dev.off()



#11
if(plot_cells == TRUE){
folder_path <- paste(sample_folder_path, "/per_cell_boxplot", sep = "")
unlink(folder_path, recursive = TRUE)
dir.create(folder_path)

reads <- reads2
reads <- scgs_normalize_reads(reads)
data <- scgs_gather_bin_depth(reads)

for(i in c(1:number_of_cells)-1){
cell_num <- i
print(paste("plotting cell", cell_num))

to_plot <- filter(data, data$cell == cell_num)
barcode <- cell_metrics$barcode[which(cell_metrics$cell_id == cell_num)]
scale_factor <- raw_somies$scale_factor[which(raw_somies$cell == cell_num)]
scale_factor <- round(scale_factor, 3)
chromo_var <- round(raw_somies$chromo_var[which(raw_somies$cell == cell_num)], 3)
baseline_ploidy <- raw_somies$baseline_ploidy[which(raw_somies$cell == cell_num)]
depth <- round(cell_metrics$effective_depth_of_coverage[which(cell_metrics$cell_id == cell_num)], 3)

cell_position <- which(raw_somies$cell == cell_num)
raw <- raw_somies[cell_position, c(2:37)]
raw <- unlist(raw)
raw <- as.numeric(raw)


cell_position <- which(integer_somies$cell == cell_num)
integers <- integer_somies[cell_position, c(2:37)]
integers <- unlist(integers)
integers <- as.numeric(integers)
karyotype_id <- karyotypes$position[which(karyotypes$karyotype == integer_somies$karyotype[which(integer_somies$cell == cell_num)])]
norm_avrg <- raw/scale_factor
int_pos <- integers/scale_factor

to_plot <- to_plot %>%
  mutate(raw_somy = raw[chromosome]) %>%
  mutate(integer_somy = integers[chromosome]) %>%
  mutate(size  = ifelse(chromosome <= 2, as.numeric(chromosome)+1,
                        ifelse(chromosome <= 10, sqrt(chromo_center)-2,
                        ifelse(chromosome <= 26, 10, sqrt(chromo_center))))) %>%
  mutate(average_depth = norm_avrg[chromosome])%>%
  mutate(integer_position = int_pos[chromosome])

  
p1 <- ggplot(to_plot, aes(x = bin_num, y = reads_per_bin, group = factor(chromosome), color = factor(chromosome)))+
  geom_boxplot(color = "black")+
  geom_line(color = "red", size = 1, aes(y = average_depth))+
  geom_text(data = to_plot, color = "black", y = max(to_plot$reads_per_bin)+0.5, aes(x = chromo_center, size = size, label = round(raw_somy, 2)))+
  labs(title = paste(sample_name, "- Cell ", cell_num, " - ", barcode, "; karyotype_id: ",karyotype_id, "; scale_factor: ",scale_factor, "; chromo_var: ", chromo_var, "; baseline_ploidy: ", baseline_ploidy, "; cell_depth: ", depth ,sep = ""))+
  guides(col = FALSE)+
  guides(size = FALSE)+
  scale_color_manual(values = rep(manhattan_plot_colors[4:5], times = 36))+
  scale_y_continuous(name = "normalized Reads/20kb")+
  scale_x_continuous(name = "Chromosome", label = unique(to_plot$chromosome), breaks = unique(to_plot$chromo_center))+
  coord_cartesian(clip = "off")+
  expand_limits(y = max(to_plot$reads_per_bin)+0.5)+
  theme_linedraw()+
  theme(axis.title = element_text(size = "16", face = "bold"), axis.text = element_text(size = 12), panel.grid.major.x = element_blank(), panel.grid.major = element_line(color = "grey"), panel.grid.minor = element_line(color = "grey"))



p2 <- ggplot(to_plot, aes(x = bin_num, y = reads_per_bin, group = factor(chromosome), color = factor(chromosome), label = integer_somy))+
  geom_point(size = 1, alpha = 0.75)+
  geom_line(color = "black", size = 1, aes(y = integer_position))+
  geom_text(color = "black", y = max(to_plot$reads_per_bin)+0.5, aes(x = chromo_center))+
  guides(col = FALSE)+
  scale_color_manual(values = rep(manhattan_plot_colors[c(1,2)], times = 36))+
  scale_y_continuous(name = "normalized Reads/20kb")+
  scale_x_continuous(name = "Chromosome", label = unique(to_plot$chromosome), breaks = unique(to_plot$chromo_center))+
  coord_cartesian(clip = "off")+
  expand_limits(y = max(to_plot$reads_per_bin)+0.5)+
  theme_linedraw()+
  theme(axis.title = element_text(size = "10", face = "bold"), axis.text = element_text(size = 12), panel.grid.major.x = element_blank(), panel.grid.major = element_line(color = "grey"), panel.grid.minor = element_blank())

p3 <- plot_grid(p1, p2, nrow = 2, align = "v", rel_heights = c(0.75, 0.25))


pdf(paste(folder_path, "/", sample_name, " - Karyotype ", karyotype_id," - Cell ", cell_num, " - ", barcode, ".pdf", sep = ""), width = 20, height = 14.15)
print(p3)
dev.off()

}
}




