#Author: Gabriel Negreira
#Description: This is a group of functions created to analyze the outputs from cellranger-dna
# The main input is the reads_raw.csv file from the cellranger, which is a table containing the number of reads mapped to each 20kb bin in each cell, as well in each 10X node. 


#####################################################################LIBRARIES######################################################################################################
library(tidyr)
library(dplyr)
library(stringr)
library(tibble)
library(ggplot2)
library(gplots)
library(RColorBrewer)
library(miscTools)
library(mixtools)
library(forcats)
library(cowplot)
library(Hmisc)
#####################################################################FUNCTIONS######################################################################################################

#setting color pallets
manhattan_plot_colors <- c("#000000", "#d28751", "#a4a4a4", "#eec15b", "#739bcc", "#88a95f", "#2a3347")
bar_plot_palette <- colorRampPalette(manhattan_plot_colors)

heat_col <- c("#001221", "#002342", "#01417a", "#0098f7", "#00ff33", "#fffa00", "#ff0000")
heat_palette <- colorRampPalette(heat_col)

heat_col <- c("#002342", "#01417a",  "#0098f7", "#00ff33", "#fffa00")
heat_palette2 <- colorRampPalette(heat_col)

heat_col <- c("#0f57ff", "#b8dae0","#ffffff", "#ffcfbf", "#ff4d4d")
heat_palette3 <- colorRampPalette(heat_col)


heat_col <- c("#001221", "#002342", "#01417a", "#0098f7", "#00ff33", "#fffa00", "#ff0000")
heat_palette4 <- colorRampPalette(rev(brewer.pal(11,"RdBu")))

heat_col <- c("#001221", "#002342", "#01417a", "#0098f7", "#03fcf8", "#00ff33", "#fffa00", "#fc8c03", "#ff0000","#eb3c25", "#b2271b", "#761e11", "#3e1408")
heat_palette5 <- colorRampPalette(heat_col)

heat_col <- c("#f5f5f5", "#1c1c1c")
bw_palette <- colorRampPalette(heat_col)

heat_col <- c("#01417a", "#2f6fab", "#92afcf", "#f8d6bd", "#f5bf92", "#f2a66c", "#f0884c", "#eb3c25", "#b2271b", "#761e11", "#3e1408", "#000000")
loupe_palette <- colorRampPalette(heat_col)

heat_col <- rev(c("#541b1b", "#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695"))
brewer_palette <- colorRampPalette(heat_col)

heat_col <- rev(c("#541b1b", "#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFFFBF", "#74ADD1", "#313695"))
somy_col <- c("black", heat_col, rep("#313695", times = 100))
bar_plot_colors <- somy_col
bar_plot_colors[3] <- "#a4a4a4"

################################################MAIN FUNCTIONS#################################################################
#_______________________________________________________________________________________________________________________________________
#1st FUNCTION: scgs_calc_somy will calculate the somy values based on the reads/20kb bin.
#if bin_filter is set to TRUE, it will apply several bin removal steps and will normalize reads by the average reads/bin of the cell. It can be set to FALSE for debugging. 
#if scale is set to TRUE, it will scale the values by a scale factor to stablish cells somies.  It can be set to FALSE for debugging. 
#scale_factor_min and scale_factor_max will determine the range of scale factors that the algorithm will try to find a scale factor for each cell. This is only used if scale is set to TRUE
#if round is set to TRUE, it will round somy values to the closest integers. 
#if method is "average" or "mean" it will use the average reads/bin of each chromosome to calculate somy. If it is set to "median", it will use the median instead.
#Use as input (x) the reads_raw file. This file should be a data frame with bins as columns (from column 2 to beyond) and cells as rows. Column 1 should contain cells id. Bins column names should contain the chromosome from which the bins are.   

scgs_calc_somy <- function(x, scale_factor_min = 1.8, scale_factor_max = 5, bin_filter = TRUE, scale = TRUE, adjust_distributions = TRUE, round = FALSE, method = "average", chromo_var_method = "max_5", chromo_var_segments = 3, plot_path){

  if(missing(plot_path)){
    plot_path <- getwd()
  }  
  
  if(bin_filter == TRUE){ #useful if dealing directly with reads_raw
    x <- x %>% 
      scgs_remove_empty_bins() %>%
      scgs_remove_dev_bins() %>%
      scgs_normalize_reads()
  }
  
  chromo_var <- scgs_calc_bin_var(x, method = chromo_var_method, number_of_segments = chromo_var_segments) #here it will calculate the intra-chromosomal variation
  
  #here it will look for the position of the bins of each chromosome
  bins_position <- scgs_find_bins(x, FALSE)
  first_chromo_bin <- min(bins_position$chromo_first_bin)
  last_chromo_bin <- max(bins_position$chromo_last_bin)
  
  #determining how many chromosomes the data has
  last_chromo <- max(bins_position$chromosome)
  
  
  
  chromo_average_depth <- c()
  chromosome <- c(1:last_chromo)
  
  print(paste("calculating the read depth of each chromosome using method ", method, "...", sep = ""))
  for(i in chromosome){#for each chromosome...
    first_bin <- bins_position$chromo_first_bin[i] #determine which is the first bin of that chromosome
    last_bin <- bins_position$chromo_last_bin[i] #determine which is the last bin of that chromosome
    
    if(method == "average" | method == "mean"){
      chromo_average_depth <- cbind(chromo_average_depth, rowMeans(x[,first_bin:last_bin], na.rm = TRUE)) #saves the average of reads per bin between first and last bin of that chromosome
      next
    }else{
      if(method == "median"){
        chromo_average_depth <- cbind(chromo_average_depth, rowMedians(x[,first_bin:last_bin], na.rm = TRUE)) #saves the median of reads per bin between first and last bin of that chromosome
        
        next
      }else{
        stop("argument method should be 'average', 'mean' or 'median'")
      }
    }
    
  }
  
  cell <- x[,-c(first_chromo_bin:last_chromo_bin)] #eliminate all bin collumns. It will be replaced by the calculated somies.
  
  #calculating the baseline somy of the entire dataset
  baseline_somy <- which.max(density(unlist(chromo_average_depth))$y) #will first determine which is the position of the maximum peak in the density of somy values
  baseline_somy <- density(unlist(chromo_average_depth))$x[baseline_somy] #will find which is the x value of the maximum peak
  baseline_somy <- 2/baseline_somy #will assume that the most frequent somy value should be 2. This is not used for somy estimation, just for adjusting somy distributions of each chromosome
  print(paste("baseline value is", baseline_somy))
  
  if(adjust_distributions == TRUE){
    print("adjusting the somy distributions to fit closest integer...")
    
    df <- data.frame(chromo_average_depth*baseline_somy)
    colnames(df) <- c(1:ncol(df))
    before <- gather(df, key = "chromosome", value = "somy")
    
    for(i in c(1:ncol(chromo_average_depth))){
      #for each chromosome, values will be multiplied by baseline somy, adjusted to closest integer and divided back by baseline somy
      chromo_average_depth[,i] <- scgs_adjust_somy_distribution(chromo_average_depth[,i]*baseline_somy)/baseline_somy
    }
    
    df <- data.frame(chromo_average_depth*baseline_somy)
    colnames(df) <- c(1:ncol(df))
    after <- gather(df, key = "chromosome", value = "somy")
    
    #will make a plot comparing before and after adjusting peaks.
    heat_col <- c("#001221", "#002342", "#01417a", "#0098f7", "#03fcf8", "#00ff33", "#fffa00", "#fc8c03", "#ff0000","#eb3c25", "#b2271b", "#761e11", "#3e1408")
    heat_palette5 <- colorRampPalette(heat_col)
    
    p1 <- ggplot(before, aes(x = somy, color = factor(chromosome, levels = sort(as.numeric(unique(chromosome))))))+
      geom_density()+
      scale_x_continuous(breaks = seq(0, round(max(before$somy))), name = "normalized average depth x baseline factor")+
      scale_color_manual(name = "Chromosome", values = heat_palette5(length(unique(before$chromosome))))+
      labs(title = paste("Before adjusting peaks - baseline factor:", round(baseline_somy, 4)))+
      theme(legend.text = element_text(size = 10), legend.key.size = unit(2, "mm"))
    p2 <- ggplot(after, aes(x = somy, color = factor(chromosome, levels = sort(as.numeric(unique(chromosome))))))+
      geom_density()+
      scale_x_continuous(breaks = seq(0, round(max(after$somy))), name = "normalized average depth x baseline factor")+
      scale_color_manual(name = "Chromosome", values = heat_palette5(length(unique(after$chromosome))))+
      labs(title = "After adjusting peaks")+
      theme(legend.text = element_text(size = 10), legend.key.size = unit(2, "mm"))
    
    plot <- plot_grid(p1, p2, nrow = 2, align = "v", rel_heights = c(0.50, 0.50))
    
    pdf(paste(plot_path, "/Before and After adjusting peaks.pdf", sep = ""), width = 11.69, height = 8.27)
    print(plot)
    dev.off()
    print("Saving plot: Before and After adjusting peaks.pdf")
  }
  
  unscaled_distance_to_integers <- apply(chromo_average_depth[,1:last_chromo]*baseline_somy, 1, integer_dist)
  unscaled_furthest_distance_to_integers <- apply(chromo_average_depth[,1:last_chromo]*baseline_somy, 1, integer_furthest)
  
  
  if(scale == FALSE){ #will return the average normalized reads/bin of each chromosome in each cell without scaling them.
    colnames(chromo_average_depth) <- c(1:ncol(chromo_average_depth))
    x <- cbind(cell, chromo_average_depth, chromo_var)
    cnames <- colnames(x)
    x <- data.frame(x)
    colnames(x) <- cnames
    print("Done")
    return(x)
  }
  
  #if scale is not FALSE, then it will apply the scaling algorithm
  chromo_average_depth <- scgs_scale_cells(chromo_average_depth, scale_factor_min = scale_factor_min, scale_factor_max = scale_factor_max)
  
  scaled_distance_to_integers <- apply(chromo_average_depth[,1:last_chromo], 1, integer_dist)
  scaled_furthest_distance_to_integers <- apply(chromo_average_depth[,1:last_chromo], 1, integer_furthest)
  #formating the data
  chromo_average_depth <- data.frame(chromo_average_depth)
  colnames(chromo_average_depth)[1:last_chromo] <- c(1:last_chromo)
  mean_somy <- rowMeans(chromo_average_depth[,1:last_chromo])
  print("calculating Gini factors of each cell...")
  Gini <- apply(chromo_average_depth[,1:last_chromo],1, gini, na.rm = TRUE)
  x <- cbind(cell, chromo_average_depth, mean_somy, chromo_var, scaled_distance_to_integers, scaled_furthest_distance_to_integers,unscaled_distance_to_integers, unscaled_furthest_distance_to_integers, Gini)
  
  print("Done")
  return(x)
}
#_________________________________________________________________________________________________________________________________
#2nd Function: scgs_solve_somy
#This function will apply gaussian mixture models (GMMs) to convert continuous somy values into integer somy values. 
#it only works well if distribution peaks were adjusted to integers. 
#There is also the option to change integer_method argument to "round". In this case, integers will simply be the result of round() function. (not recommended)
#Use as input the output of scgs_calc_somy() function

scgs_solve_somy <- function(x, integer_method = "GMM", save_plot = TRUE, sample_name, plot_path){
  
  plot_path <- ifelse(missing(plot_path), getwd(), plot_path)
  
  mode <- "heuristic"
  
  #if sample_name is not give, it will assume the name of x
  if(missing(sample_name)){
    sample_name <- deparse(substitute(x))
  }
  
  #defining the columns that represent the chromosomes...
  first_chromo <- which(colnames(x) == "1")
  last_chromo <- NULL
  
  if("scale_factor" %in% colnames(x)){
    last_chromo <- which(colnames(x) == "scale_factor")-1
  }else{
    last_chromo <- which(colnames(x) == "mean_somy")-1
  }
  
  if(is.null(last_chromo)){
    stop("Could not define which column represent the last chromosome")
  } 
  
  
  data <- x[,c(first_chromo:last_chromo)] #subseting the somy values
  if(integer_method == "round"){
    print("rounding somy values to the closest integer...")
    somies <- round(data)
    
    if(save_plot == TRUE){
      chromosomes <- c(1:ncol(data))
      for(i in chromosomes){
        to_plot <- data.frame(somy = data[,i])
        plot <- ggplot(to_plot, aes(x = somy))+
          geom_histogram(aes(y = ..density..), bins = 100)+
          geom_vline(xintercept = c(0:max(round(to_plot$somy)))+0.5, color = "red", linetype = "dashed")+
          scale_x_continuous(limits = c(0, max(round(to_plot$somy))+1), breaks = c(0:max(round(to_plot$somy))+1))
          labs(title = paste(sample_name, "rounding to the closest integer"))+
          theme(title = element_text(size = 7))
        
        #saving the plot
        plot_name <- paste("Chromosome", i)
        print(paste("Saving plot:", plot_name))
        pdf(paste(plot_path, "/", plot_name, ".pdf", sep = ""), width = 11.69, height = 8.27)
        print(plot)
        dev.off()
      }
    }
  }
  if(integer_method == "GMM"){
    # #creating a folder to store the data of the sample
    
    print("Applying Gaussian Mixture Models to solve somies...")
    somies <- c()
    for(i in c(1:ncol(data))){ #for each chromosome...
      #will run the scgs_solve_by_GMM function to solve the somy of each chromosome
      values <- try(scgs_solve_by_GMM(data[,i], save_plot = save_plot, plot_path = plot_path, mode = mode, plot_name = paste("Chromosome", i), sample_name = sample_name))
      #sometimes the mixtools function used in to create the GMM returns an error, so it will repeat in that case.
      attempts <- 0
      while(grepl("Error", values) & attempts < 20){
        print("retrying...")
        attempts <- attempts + 1
        values <- try(scgs_solve_by_GMM(data[,i], save_plot = save_plot, plot_path = plot_path, mode = mode, plot_name = paste("Chromosome", i), sample_name = sample_name))
      }
      if(grepl("Error", values)){
        setwd(current_wd)
        stop("Could not solve distributions with GMM")
      }
      somies <- cbind(somies, values)
    }
    setwd(current_wd)
  }
  
  diff <- rowSums(abs((data/rowMeans(data))-(somies/rowMeans(somies)))) #here it will calculate how much the data changed
  
  karyotype <- c()
  for(i in c(1:nrow(somies))){
    karyotype <- c(karyotype, paste(somies[i,], collapse = "_"))
  }
  
  karyo_table <- data.frame(table(karyotype), stringsAsFactors = FALSE)
  karyo_table$frequency <- karyo_table$Freq/sum(karyo_table$Freq)
  karyo_table <- arrange(karyo_table, 1-frequency)
  karyo_table$position <- c(1:nrow(karyo_table))
  
  karyo_frequency <- rep("none", times = length(karyotype))
  karyo_occurrence <- rep("none", times = length(karyotype))
  karyo_position <- rep("none", times = length(karyotype))
  
  for(i in c(1:nrow(karyo_table))){
    karyo_frequency[c(which(karyotype %in% karyo_table$karyotype[i]))] <- karyo_table$frequency[i]
    karyo_occurrence[c(which(karyotype %in% karyo_table$karyotype[i]))] <- karyo_table$Freq[i]
    karyo_position[c(which(karyotype %in% karyo_table$karyotype[i]))] <- karyo_table$position[i]
  }
  
  x[,first_chromo:last_chromo] <- somies
  
  karyo_frequency <- as.numeric(karyo_frequency)
  karyo_occurrence <- as.numeric(karyo_occurrence)
  
  if("karyotype" %in% colnames(x)){
    x[,which(colnames(x) == "karyotype")] <- karyotype
  }
  
  if("karyo_occurrence" %in% colnames(x)){
    x[,which(colnames(x) == "karyo_occurrence")] <- karyo_occurrence
  }
  
  
  return(x)
}

########################################################Sub-functions of the main functions#################################################################
#These functions are often not usefull for separate use. They are called inside the main functions above. 

#scgs_find_bins
#this function will determine the intervals where the bins of each chromosome are. 
#Works for any data set where the bin names (column names) mention the chromosome number.
#It will return a dataframe with the chromosome number, the first bin and the last bin of that chromosome

scgs_find_bins <- function(x, verbose = TRUE){
  #First it will determine in which column is the first bin. It does so by looking for the first bin name with "01_":
  bin_start <- min(grep("01_", colnames(x)))
  
  #now the code will determine how many chromosomes are present. 
  #since it already know which is the name of the first bin of chromosome 01, it will use the characters that flank "01" to determine the other chromosomes:
  
  flank <- strsplit(colnames(x[bin_start]), split = "01") #first it determines which are the characters that flank the chromosome number
  flank <- unlist(flank)
  flank[2] <- unlist(strsplit(flank[2], split = ""))[1] #here it takes the first character that flanks the chromosome number on the right
  flank <- unlist(flank)
  
  #once the flanking characters were determined, it will use this information to get the number of the last chromosome
  last_chromo <- max(grep(flank[1], colnames(x))) #it will assume that the last bin that contains the flanking characters is the last bin of the last chromosome
  last_chromo <- colnames(x)[last_chromo]
  last_chromo <- sub(flank[1], "", last_chromo) #will replace the flanking characters on the left by nothing
  last_chromo <- unlist(strsplit(last_chromo, split = flank[2]))[1] #will return the number before the character flanking the chromosome number on the right
  last_chromo <- as.numeric(last_chromo)
  
  #Here, it will look for the first bin of each chromosome
  number_of_bins <- c()
  chromosome <- c()
  chromo_first_bin <- c()
  chromo_last_bin <- c()
  
  if(verbose == TRUE){
    print(paste("this data contains", last_chromo, "chromosomes"))
  }
  
  for(i in c(1:last_chromo)){ #for each chromosome, it will first define a pattern to determine which bins correspond to that chromosome
    if(i < 10){
      pattern <- paste(c(0, i, flank[2]), collapse = "")
    }else{
      pattern <- paste(c(i, flank[2]), collapse = "")
    }
    
    first_bin <- min(grep(pattern, colnames(x))) #determines the column index of the first bin of that chromosome
    last_bin <- max(grep(pattern, colnames(x))) #determines the column index of the last bin of that chromosome
    number_of_bins <- c(number_of_bins, last_bin-first_bin)
    
    
    if(verbose == TRUE){
      print(paste("Chromosome", i, "is located between column", first_bin, "and", last_bin, "showing", last_bin - first_bin, "bins" ,sep = " "))
    }
    
    chromosome <- c(chromosome, i)
    chromo_first_bin <- c(chromo_first_bin, first_bin)
    chromo_last_bin <- c(chromo_last_bin, last_bin)
    
  }
  
  df <- data.frame(chromosome, chromo_first_bin, chromo_last_bin, number_of_bins, stringsAsFactors = FALSE)
  return(df)
}



#scgs_normalize_reads
#this function will simply divide the number of reads of each bin by the average reads/bin of the cell
scgs_normalize_reads <- function(x, norm = "cell", method = "mean"){
  bins_position <- scgs_find_bins(x, FALSE)
  first_bin <- min(bins_position$chromo_first_bin)
  last_bin <- max(bins_position$chromo_last_bin)
  
  if(norm == "cell"){
    if(method == "median"){
      print("normalizing reads by cells' median...")
      average <- rowMedians(x[,first_bin:last_bin])
    }else{
      print("normalizing reads by cells' average...")
      average <- rowMeans(x[,first_bin:last_bin])
    }
    
    x[first_bin:last_bin] <- x[first_bin:last_bin]/average
  }else{
    if(norm == "chromosome"){
      norm_depth <- c(1:nrow(x))-1
      for(i in c(1:nrow(bins_position))){
        first_chromo_bin <- bins_position$chromo_first_bin[i]
        last_chromo_bin <- bins_position$chromo_last_bin[i]
        if(method == "median"){
          print("normalizing reads by chromosomes median...")
          ave_depth <- rowMedians(x[,first_chromo_bin:last_chromo_bin])
        }else{
          print("normalizing reads by chromosome average...")
          ave_depth <- rowMeans(x[,first_chromo_bin:last_chromo_bin])
        }
        
        norm_values <- x[,first_chromo_bin:last_chromo_bin]/ave_depth
        norm_depth <- cbind(norm_depth, norm_values)
      }
      
      norm_depth <- norm_depth[,-1]
      x[,first_bin:last_bin] <- norm_depth
    }else{
    stop("norm should be \"cell\" or \"chromosome\"")
    }
  }
  return(x)
}

#scgs_calc_bin_dev
#this function will calculate the deviation from the chromosomal median of each bin and return the values in a log2 scale. 
#there is the option to skip log2 transformation
#if gather is set to true, it will return a gathered data frame ready for ggplot
#set_zero_to replaces values = 0 to any other value
#This is useful for downstream removal of highly deviating bins as well as to determine coverage eveness

scgs_calc_bin_dev <- function(reads, set_zero_to = NA, return_log2 = TRUE, gather = FALSE, verbose = TRUE){
  
  print("calculating the deviation of each bin from the chromosomal median...")
  
  #First it will determine the position of the first and the last bin of each chromosome
  bins_position <- scgs_find_bins(reads, verbose = FALSE)
  last_chromo <- max(bins_position$chromosome) #determines which is the last chromosome
  
  #getting the cell ids:
  if("cell" %in% colnames(reads)){
    cell_ids <- reads$cell
    cell_ids <- unlist(cell_ids) 
  }else{
    if("X1" %in% colnames(reads)){
      cell_ids <- reads$X1
      cell_ids <- unlist(cell_ids)
    }else{
      cell_ids <- c(1:nrow(reads))-1
      print(paste("assuming that cell_ids go from", cell_ids[1], "to", cell_ids[tail(cell_ids, n = 1)]))
    }
  }
  
  median_read_count <- c()
  chromo_center <- c()
  variation <- cell_ids
  
  
  for(i in c(1:last_chromo)){ #for each chromosome, it will first determine the first and last bin
    first_bin <- bins_position$chromo_first_bin[i]
    last_bin <- bins_position$chromo_last_bin[i]
    
    chromo_center <- c(chromo_center, median(c(first_bin, last_bin)))
    
    if(verbose == TRUE){
      print(paste("Calculating median between bins", first_bin, "and", last_bin))
    }
    
    #it will then comput the average read count per bin of the chromosome in each cell
    median_read_count <- cbind(median_read_count, rowMedians(reads[,first_bin:last_bin]))
    
    #here it will comput the deviation from the average of each bin
    for(j in c(first_bin:last_bin)){
      var <- reads[,j]/median_read_count[,i]
      variation <- cbind(variation, var)
    }
  }
  df <- variation[,-1]
  
  first_bin <- bins_position$chromo_first_bin[1]
  last_bin <- bins_position$chromo_last_bin[nrow(bins_position)]
  
  colnames(df) <- colnames(reads)[first_bin:last_bin]
  df <- data.frame(df, stringsAsFactors = FALSE)
  df[df == 0] <- set_zero_to
  
  if(return_log2 == TRUE){
    df <- log(df, 2)
  }
  
  df <- cbind(cell_ids, df)
  colnames(df)[1] <- "cell"
  
  if(gather == FALSE){
    return(df)
  }else{
    df <- gather(df, key = "bin", value = "log_variation", -cell)
    df$chromosome <- sub("_.*", "", df$bin)
    df$chromosome <- sub(flank[1], "", df$chromosome)
    df$chromosome <- as.numeric(df$chromosome)
    df$chromo_center <- chromo_center[df$chromosome]
    bin_num <- c(1:length(unique(df$bin)))
    bin_num <- rep(bin_num-1, each = length(unique(df$cell)))
    df$bin_num <- bin_num
    return(df)
  }
}

#scgs_remove_empty_bins
#this function will simply remove any bin where almost no read is maping to it. 

scgs_remove_empty_bins <- function(x){
  bins_position <- scgs_find_bins(x, FALSE) #will define which are the bins for each chromosome
  first_bin <- min(bins_position$chromo_first_bin) #define the first bin of the first chromosome
  last_bin <- max(bins_position$chromo_last_bin) #define the last bin of the last chromosome
  
  print("Removing empty bins...")
  cell <- x
  x <- x[,first_bin:last_bin] #subset only the part of the dataframe that contains the  bins
  before <- ncol(x)
  x <- x[,colSums(x) > nrow(x)/2] #here, bins that have 0 or very few reads mapping to them will be removed
  x <- x[,colSums(x) != nrow(x)*128] #this will remove bins composed by CNV values = 128 in case x is the cnv_data. 
  after <- ncol(x)
  print(paste(before-after, "bins removed"))
  cell <- cell[,-c(first_bin:last_bin)]
  x <- cbind(cell, x)
  
  return(x)
}

#scgs_remove_dev_bins
#this function will remove bins that deviate from the chromosomal median reads/bin.
#it is adiviced to run this function only after removing unmapable bins with scgs_remove_empty_bins()

scgs_remove_dev_bins <- function(reads, set_zero_to = 0.00000001){
  
  #calculating the deviation from chromosomal median of each bin
  reads_dev <- scgs_calc_bin_dev(reads, set_zero_to = set_zero_to, gather = FALSE, verbose = FALSE)
  
  print("Removing bins that deviate much from chromosomal median...")
  #determining the position of the first and last bin
  bins_position <- scgs_find_bins(reads_dev, FALSE) #gets the column index of each bin
  first_bin <- min(bins_position$chromo_first_bin) #define the column index of the first bin of the first chromosome
  last_bin <- max(bins_position$chromo_last_bin) #define the column index of the last bin of the last chromosome
  
  #subseting data to bins
  dev_values <- reads_dev[,first_bin:last_bin] #subset only the columns containing the bins
  
  #calculating the average deviation from the median of each bin
  dev_values <- t(dev_values) 
  dev_values <- cbind(dev_values, rowMedians(dev_values, na.rm = TRUE))
  dev_values <- data.frame(dev_values)
  last_col <- ncol(dev_values)
  colnames(dev_values)[last_col] <- "median"
  
  #determining outliers
  outliers <- boxplot.stats(dev_values$median)$out
  
  #removing outlier bins
  before <- nrow(dev_values)
  
  bin <- rownames(dev_values)
  dev_values <- cbind(bin, dev_values)#this is done to preserve row names after filtering
  dev_values <- filter(dev_values, !is.na(dev_values$median))
  dev_values <- filter(dev_values, !dev_values$median %in% outliers)
  rownames(dev_values) <- dev_values$bin
  dev_values$bin <- NULL
  
  after <- nrow(dev_values)
  
  print(paste(before - after, "bins removed"))
  
  #returning bins that passed the outliers removal
  dev_values$median <- NULL #remove the "average" column
  dev_values <- t(dev_values) 
  dev_values <- data.frame(dev_values)
  
  to_remove <- which(!colnames(reads) %in% colnames(dev_values)) #determine the bins that were removed
  
  other_data <- reads[, -c(first_bin:last_bin)] #will save any information that is not a bin
  other_data <- data.frame(other_data)
  if(ncol(other_data) == 1){
    colnames(other_data) <- "cell"
  }
  
  reads[,to_remove] <- NULL #remove bins that are not present in the dev_values data (the filtered bins)
  
  reads <- cbind(other_data, reads)
  
  rownames(reads) <- c(1:nrow(reads))
  
  if(!"cell" %in% colnames(reads)){
    cell <- c(1:nrow(reads))-1
    reads <- cbind(cell, reads)
    warning(paste("No column with name 'cell' was present. A column with cell_ids was introduced based on row numbers"))
  }
  
  return(reads)
}


#scgs_find_scale_factor
#this function will find a factor of which the average normalized reads/bin values of each chromosome should be multiplied by it in orther to stablish the somy values
#it takes a numeric vector containing the normalized reads/bin values of each chromomosome and find a factor that makes the numbers closest to integers
#in most cases, this function is not used individually, but as a subfunction of scgs_calc_somy
scgs_find_scale_factor <- function(x, range = c(2,4), plot = FALSE, plot_name = "scale_factor"){
  length <- (range[2]-range[1]) * 500
  #x must be a vector containing the normalized average reads/bin of each chromosome
  x <- unlist(x) 
  values_to_try <- seq(range[1], range[length(range)], length = length) #creates a vector containing 'length' values between 'range'
  
  values <- c()
  integer_distances <- c()
  
  for(i in c(1:length(values_to_try))){ #for each value in the created vector...
    values <- c(values, values_to_try[i]) #saves the value that is being tested...
    test <- x*values_to_try[i] #multiplies the normalized reads by the value...
    distance_to_integer <- integer_dist(test) #calculate how far from integers the normalized reads get when multiplied by the value being tested...
    integer_distances <- c(integer_distances, sum(distance_to_integer)) #save the sum of the distances to integers of all chromosomes when multiplied
  }
  

  scale_factor <- which(integer_distances == min(integer_distances)) #will select the index of the values that result in the minimum distances to integers in the created vector...
  scale_factor <- min(values[scale_factor]) #it will select the minimum value that resuts in the minimum distance to integers as a scale factor
  
  if(plot == TRUE){
    path <- paste(getwd(), "/scale_factors", sep = "")
    dir.create(path = path, showWarnings = FALSE)
    to_plot <- data.frame(scale_factor = values, mean_distance_to_integers = integer_distances)
     plot <- ggplot(to_plot, aes(x = scale_factor, y = mean_distance_to_integers))+
              geom_point(alpha = 0.3)+
              labs(title = plot_name)+
              geom_vline(xintercept = scale_factor, color = "red", linetype = "dashed")+
              scale_x_continuous(limits = range, breaks = c(range[1]:range[2]))
    
    pdf(paste(path, "/", plot_name, ".pdf", sep = ""),  width = 5.83, height = 4.13)
    print(plot)
    dev.off()
  }
  
  
  return(scale_factor)
}

#scgs_scale_cells
#this function scales the normalized read depth of each cell by a factor in order to stablish the somy values
#The factor is determined as a factor that minimizes intermediate somy values. see the function scgs_find_scale_factor
#in most cases, this function is not used individually, but as a subfunction of scgs_calc_somy
scgs_scale_cells <- function(x, scale_factor_min = 2, scale_factor_max = 3){ 
  print("scalling cells to stablish somy values...")
  scale_factor <- c()
  for(i in c(1:nrow(x))){
    #print(paste("finding value for cell", i))
    scale <- scgs_find_scale_factor(x[i,], range = c(scale_factor_min,scale_factor_max)) # here it will find the best scale factor to fit the data in each cell
    scale_factor <- c(scale_factor, scale) #saving the scale factor of each cell
    x[i,] <- x[i,]*scale #multiply the normalized read depth values by the scale factor
  }
    x <- cbind(x, scale_factor)
  return(x)
}


#scgs_calc_bin_var
#thins function will calculate the intrachromosomal variation score of each cell and return a vector with the caculated values
#method can be "max", which will assign the maximum intrachromosomal variation to the score of each cell, or "average" which will take the average intrachromosomal variation of all chromosomes for each cell. Method "max_N" will calculate the average between the N chromosomes highest variations (replace N by a number).
scgs_calc_bin_var <- function(x, method = "max_5", number_of_segments = 3){
  print(paste("calculating intrachromosomal variation of each cell by dividing chromosomes in", number_of_segments, "segments using method:",method, "..."))
  #here it will look for the position of the bins of each chromosome
  bins_position <- scgs_find_bins(x, FALSE)
  first_chromo_bin <- min(bins_position$chromo_first_bin)
  last_chromo_bin <- max(bins_position$chromo_last_bin)
  
  #determine how many chromosomes the data has
  last_chromo <- max(bins_position$chromosome)
  #create a vector with the chomosomes
  chromosomes <- c(1:last_chromo)
  #pre-allocate vectors 
  bin_var <- numeric(length = length(chromosomes))
  cell_values <- numeric(length = length(chromosomes))
  #for each chromosome...
  for(i in chromosomes){
    #determine the index of the first and the last bin corresponding to that chromosome...
    chromo_first_bin <- bins_position$chromo_first_bin[i]
    chromo_last_bin <- bins_position$chromo_last_bin[i]
    #determine how many bins that chromosome has...
    number_of_bins <- length(c(chromo_first_bin:chromo_last_bin))
    #determine how many bins correspond to the size o a intrachromosomal segment...
    window_size <- number_of_bins/number_of_segments
    #prevent that rounding the window_size number leads to more windows thant the actual size of the chromosome
    if((window_size-round(window_size)) < 0){
      window_size <- round(window_size)-1
    }else{
      window_size <- round(window_size)
    }
    first_bin <- chromo_first_bin
    #pre-allocate the chromo_values vector
    chromo_values <- c()
    #create a vector that will be used to move the window across the chromosomoe
    window_move <- c(0, rep(1, times = number_of_segments - 1))
    #for each window...
    for(j in window_move){
      #define the first bin of that window
      first_bin <- first_bin + (window_size*j)
      #define the last bin of that window
      last_bin <- first_bin + window_size
      #in case the last bin index is further than the end of the chromosome... 
      if(last_bin > chromo_last_bin){
        #use the last bin of the chromosome instead
        last_bin <- chromo_last_bin
        first_bin <- last_bin-window_size
      }
      #get the depth values of each bin inside that window
      values <- x[, c(first_bin:last_bin)]
      #calculate the average between the bins inside the window. 
      values <- rowMeans(values)
      #append the mean value to the chromo_values vector
      chromo_values <- cbind(chromo_values, values)
    } #end of windows loop
    
    #if there are 0 values in the chromo_values, convert them to the lowest non-0 value divided by 10. 
    if(0 %in% chromo_values){
      chromo_values[chromo_values == 0] <- min(chromo_values[chromo_values != 0])
    }
    #check which is the minimum value obtained in the chromo_values vector
    min <- apply(chromo_values, 1, min, na.rm = TRUE)
    #if the minimum is 0, convert it to minimum value found other than 0 (to allow downstream division)
    #check which is the maximum value obtained in the chromo_values vector
    max <- apply(chromo_values, 1, max, na.rm = TRUE)
    #calculate the ratio between the window with the maximum and the window with minimum mean depth. 
    chromo_values <- max/min
    #append the chromo_var value of chromosome i to the cell_values vector.
    cell_values <- cbind(cell_values, chromo_values)
  }#end of chromosome loop
  
  #resume the cell_values vector to a single value for each cell by using the method...
  #if after 'max' there is a '_' and a number, it will get the mean between the N chromosomes with the highest variation
  if(grepl("max_", method)){
    method <- unlist(strsplit(method, split = "_"))
    method <- method[length(method)]
    method <- as.numeric(method)
    method <- round(method)
    if(is.na(method)){
      stop("when using method max, ensure there is a number after \'_\'")
    }
    #create a function to calculate the mean between the N highest values 
    top_average <- function(x, na.rm = FALSE, top = 5){
      if(na.rm){
        if(NA %in% x){
          x <- x[-which(is.na(x))]
        }
      }
      x <- sort(x, decreasing = TRUE)[c(1:top)]
      x <- mean(x)
      
      return(x)
    }
    bin_var <- apply(cell_values, 1, top_average, na.rm = TRUE, top = method)
  }else{
    if(method == "max"){
      #the value of the chromosome with the maximum variation will be assigned to the cell. 
      bin_var <- apply(cell_values, 1, max, na.rm = TRUE)
    }else{
      #the variation of all chromosomes in a cell will be averaged and used as a summary to the overall intrachromosomal variation of each cell.
      bin_var <- rowMeans(cell_values, na.rm = TRUE)
    }
  }

  #NA values will be replaced by 0s. 
  bin_var[is.na(bin_var)] <- 0
  #return the values 
  return(bin_var)
}

#scgs_adjust_somy_distribution
#this function will calculate a density distribution of somy values for a given chromosome, and adjust that distribution so peaks suround the closest integer

scgs_adjust_somy_distribution <- function(x){
  
  density <- density(x)
  y_max <- which(density$y == max(density$y))
  peak1 <- density$x[y_max]
  
  adjust_factor <- round(peak1)/peak1 #will calculate how far the peak is from its closest integer
  x <- x*adjust_factor #will multiply x to make the peak of the distribution to became an integer

  return(x)
}

#scgs_solve_by_GMM 
#this function will take a vector containing somy values and will apply a gaussian mixture model to cluster the values to integers.
#This is a more robust way of clustering somy values to integers than rounding them. 
#This function is mainly used as a subfunction of scgs_solve_somy

scgs_solve_by_GMM <- function(x, return = "vector", mode = "heuristic", save_plot = FALSE, plot_name, plot_path, sample_name){
  plot_name <- ifelse(missing(plot_name), "non-defined", plot_name)
  plot_path <- ifelse(missing(plot_path), getwd(), plot_path)
  sample_name <- ifelse(missing(sample_name), "non-defined", sample_name)

  somies <- unique(round(x)) #will determine which are the integer somy values for the vector
  #somies <- sort(somies)
  
  if(0 %in% somies){
    somies <- somies[which(somies != 0)] #will remove zeros. 
  }
  
  cells_in_somies <- c() #here it will count how many cells display a somy value close to the integer somies. This is used to decide if a gaussian will be representative or not
  for(i in c(1:length(somies))){
    min <- somies[i]-0.2
    max <- somies[i]+0.2
    cells_in_somies <- c(cells_in_somies, length(which(x >= min & x <= max)))
  }
  
  ncells <- length(x)
  reliable_gaussian <- ifelse(cells_in_somies > round(0.005 * ncells), TRUE, FALSE)
  #print(reliable_gaussian)
  print(somies)
  states <- length(somies) #will determine the number of somy states (k)
  
  #this will create a vector to define means.constrain argument so the means of the generated gaussians should follow the ratios between integer somies.
  #previously I was defining the mean.constrain as fixed in the integer somies. The use of `ratios` allow a flexibility in defining the means surounding the integer somies.
  ratios <- c() 
  sd.constr <- c()
  for(i in c(1:length(somies))){
    ratios <- c(ratios, somies[i]/somies[(1)])
    sd.constr <- c(sd.constr, ifelse(reliable_gaussian[i] == TRUE, NA, 0.1))
  }
  
  
  ratios <- paste(ratios, "a", sep = "")#so what it will "tell" the EM algorithm is "mean 2 should be x times higher than mean 1, mean 3 should be y times higher than mean 1" and so on)
  
  sds <- 0
  while(length(which(sds < 0.01)) != 0){ #sometimes the EM algorithm creates a gaussian with a very small stdev, ecompassing only few observations. If this happens, the EM will be repeated
    #this will create a Gaussian mixture model where the starting means are equal to the possible integer somies in x, but it will allow means to change as long as they keep the ratios (defined in `ratios`) between them.
    MModel <- normalmixEM(x, k = states, mu = somies, sd.constr = sd.constr, mean.constr = ratios, verb = FALSE, maxit = 10000) #will create a GMM for the data with fixed means
    
    while(length(MModel$all.loglik) < 5){ #sometimes the EM algorithm does only 1 iteration and this creates a lot of errors. Here, if few iterations were perfomed, it will repeat the algorithm
      MModel <- normalmixEM(x, k = states, mu = somies, sd.constr = sd.constr, mean.constr = ratios, verb = FALSE, maxit = 10000)
    }
    
    #sometimes the EM algorithm finds solutions that respect the `ratios` but shift means too far  from the integer somies. If this happens, it will fix the means precisely in the integer somies. 
    #In this case, the means of the gaussians are not flexible anymore. 
    if(round(MModel$mu) != somies){
      print("OBS: mean.constr had to be fixed at integer somies.")
      MModel <- normalmixEM(x, k = states, mu = somies, sd.constr = sd.constr, mean.constr = somies, verb = FALSE, maxit = 10000)
      
      while(length(MModel$all.loglik) < 5){ 
        MModel <- normalmixEM(x, k = states, mu = somies, sd.constr = sd.constr, mean.constr = somies, verb = FALSE, maxit = 10000)
      }
    }
    sds <- MModel$sigma
  }
  nclusters <- length(MModel$mu) #determine the number of clusters in the model. Should be the same as `states`
  means <- MModel$mu #determine the integer somy values. Should be close to `somies`
  sds <- MModel$sigma #determine the standard deviation of each gaussian. Useful for plots.
  proportions <- MModel$lambda #determine the proportion of each gaussian. Useful for plots
  print(sds)
  
  if(save_plot == TRUE){
    #making the plot
    x <- data.frame(somy = x)
    plot <- ggplot(x, aes(x = somy))+
      geom_histogram(aes(y = ..density..), bins = 100)+
      labs(title = paste(sample_name, plot_name, paste("means =", paste(round(means,3), collapse = "; ")),paste("proportions =", paste(round(proportions,3), collapse = "; ")), paste("sds =", paste(round(sds,3), collapse = "; ")), sep = " - " ))+
      theme(title = element_text(size = 7))
    for(i in c(1:nclusters)){#for each cluster, it will add the Gaussian to the plot
      plot <- plot + stat_function(geom = "line", fun = fun_prop, args = list(mean = means[i], sd = sds[i], proportion = proportions [i]), color = i)
    }
    #saving the plot
    print(paste("Saving plot:", plot_name))
    pdf(paste(plot_path, "/", plot_name, ".pdf", sep = ""), width = 11.69, height = 8.27)
    print(plot)
    dev.off()
  }
  
  output <- data.frame(somy = MModel$x)#save the raw input values
  output2 <- MModel$posterior#save the posterior values, i.e., the probability of each value to belong to each cluster
  
  colnames(output2) <- round(MModel$mu) #this function uses the names of the columns to find the name of each cluster
  
  for(i in c(1:length(colnames(output2)))){
    #since colnames are strings, to order them properly, numbers lower than 10 must have a 0 in front of them.
    if(as.numeric(colnames(output2)[i]) > 9){ 
      next
    }else{
      number <- colnames(output2)[i]
      number <- paste(0, number, sep = "")
      colnames(output2)[i] <-number
    }
  }
  output2 <- output2[,sort(colnames(output2))] #arranges the columns
  
  flank1 <- min(as.numeric(colnames(output2)))-1
  if(flank1 < 10){
    flank1 <- paste(0, flank1, sep = "")
  }
  
  flank2 <- max(as.numeric(colnames(output2)))+1
  if(flank2 < 10){
    flank2 <- paste(0, flank2, sep = "")
  }
  
  output2 <- cbind(0,output2, 0) #this is done to flank the data with probabilities 0, so when the heuristic algorithm tries to decide between the extremes of the data, it will always choose the extreme, not data beyond the extreme.
  colnames(output2)[1] <- flank1
  colnames(output2)[ncol(output2)] <- flank2
  output <- cbind(output, output2)
  
  #here it will create a data.frame with several informations about the calculations. This is more used to debug. That's what returned when return == "complete"
  output <- output %>%
    mutate(by_rounding = round(somy)) %>%
    mutate(by_clustering = as.numeric(colnames(output[c(1:nclusters)+2])[max.col(output[,c(1:nclusters)+2])])) %>% #+2 is used because the first 2 data are not probabilities
    mutate(agree = ifelse(by_rounding == by_clustering, TRUE, FALSE)) %>%
    mutate(distance_to_round = somy-by_rounding) %>%
    mutate(distance_to_cluster = somy-by_clustering)
  
  #here it will return the the values clustered to the group with the higher posterior, regardless of the proximity between the value and the group.
  if(mode == "clustering"){
    if(return == "vector"){
      return(as.numeric(output$by_clustering))
    }
  }
  
  #if mode == "heurisitic", it will assign each value to one of the integer numbers that surround that value, based on the posterior value of the cluster representing each integer.
  
  if(mode == "heuristic"){
    by_heuristic <- c()
    
    for(i in c(1:nrow(output))){ #for each cell...
      somy <- output[i,1]
      
      if(somy <= 0.1){#if raw somy is too low, it will be considered nullisomic
        by_heuristic <- c(by_heuristic, 0)
        next
        
      }
      integers_col <- which(as.numeric(colnames(output)) == round(somy)) #in the generated 'output' data.frame, it will define the column where the cluster that correspond to the rounded somy value is  
      if((somy - round(somy) < 0)){ #if raw somy is lower than the closest integer, the comparison will be done between the rounded integer of the somy value and the integer before it 
        integers_col <- c(integers_col-1, integers_col)
      }else{
        if((somy - round(somy)) == 0){ #if there is no difference between the raw somy and the rounded somy, it will not compare the values with any other column
          integers_col <- c(integers_col, integers_col)
        }else{ #if raw somy is higher than the closest integer, the comparison will be done between the rounded integer of the somy value and the integer after it
          integers_col <- c(integers_col, integers_col+1)
        }
        
      }
      probs <- output[i, integers_col] #selecting the probabilities of the somy value beloging to group to each group
      position <- min(which(probs == max(probs))) #chooses which is the highest probability
      position <- integers_col[position] #choose the position of the column that represents the cluster where the somy value has the highest probability to belong
      by_heuristic <- c(by_heuristic, as.numeric(colnames(output)[position])) #save the integer value represented by the chosen cluster
    }
    
    output$by_heuristic <- by_heuristic
    
    if(return == "vector"){ #this will return only the calculated integer values
      return(as.numeric(output$by_heuristic))
    }
    if(return == "complete"){ #complete this will return the complete dataframe. Usefull for debugging and to understand better the calculations.
      return(output)
    }
  }else{
    stop(paste("mode should be equal to 'heuristic' or 'clustering'"))
  }
}


########################################Functions to work with karyotypes####################################################
#karyo_uniques:

#the function will give as an output a dataframe with all different karyotypes ordered by the number of times each karyotype appears in the population together with other useful info.
#the function round the somy values of each chromosome and concatenates all somies of each cell in a single string from chromosome 1 to 36. Each karyotype will be representended by a unique numeric code formed by 36 numbers representing the somy values of each chromosome in that karyotype from chr1 to 36
#OBS: if not using a table generated by 10X, give as input a table containing cells as collumns and somies as rows, and set arguments to FALSE.  

karyo_uniques <- function(x, round_method = 1, decluster = TRUE, add_node_id = TRUE){
  
  backup <- x
  
  if (decluster == TRUE){
    x <- karyo_decluster(x, TRUE, TRUE, TRUE)
  }else{
    rownames(x) <- x$barcode
  }
  x$barcode <- NULL
  x$node_id <- NULL
  x$barcodes <- NULL
  x$num_cells <- NULL
  x$num_noisy <- NULL
  x$average_somy <- NULL
  x <- t(x)
  
  x <- data.frame(x, stringsAsFactors = FALSE)
  
  x <- numeric.matrix(x)
  
  if (round_method == 1){
    print("karyo_uniques is running with round method 1")
    x <- round(x)
  }
  
  if (round_method == 2){
    print("karyo_uniques is running with round method 2")
    x <- round2(x)
  }
  
  v <- c()
  occurrences <- c()
  diff_from_the_first <- c()
  
  for (i in c(1:ncol(x))){
    c <- as.character(paste(x[,i], collapse = "_"))
    v <- c(c, v)
  }
  
  v <- table(v)
  v <- data.frame(v, stringsAsFactors = FALSE)
  colnames(v) <- c("karyotype", "number_of_cells")
  v <- df_order(v, 2)
  rownames(v) <- c(1:nrow(v))
  v$frequency <- v$number_of_cells/ncol(x)
  
  main_karyo <- unlist(strsplit(as.character(v[1,1]), "_"))
  
  for(karyotype in c(1:nrow(v))){
    karyo_to_compare <- unlist(strsplit(as.character(v[karyotype,1]), "_"))
    count <- 0
    for (i in c(1:length(main_karyo))){
      if (karyo_to_compare[i] != main_karyo[i]){
        count <- count+1
      }
    }
    diff_from_the_first <- c(diff_from_the_first, count)
  }
  
  v$diff_from_the_first <- diff_from_the_first
  
  
  b <- 0
  cummulative <- c()
  
  for(i in c(1:nrow(v))){
    b <- v[i,4]+b
    cummulative <- c(cummulative, b)
  }
  
  v$cummulative <- cummulative
  
  if(add_node_id == TRUE){
    v <- karyo_add_node_id(v, backup, round_method)
  }
  
  return(v)
}


#karyo_gather
#this function will create a dataframe with all relevant info for all karyotypes in a format ready for ggplot. It applies the gather function from the dyplr library to the output of the karyo_somies function.

karyo_gather <- function(x, round_method = 1, decluster = TRUE){
  
  print("karyo_gather is running")
  
  x2 <- x
  x$node_id <- NULL
  x$barcode <- NULL
  x$barcodes <- NULL
  x$num_cells <- NULL
  x$num_noisy <- NULL
  x$average_somy <-NULL
  
  somies <- karyo_somies(x2, round_method, decluster)
  
  
  if(decluster == TRUE){
    karyotypes_info <- karyo_uniques(x2, round_method)
  }else{
    karyotypes_info <- karyo_uniques(x2, round_method, decluster = FALSE, add_node_id = FALSE) 
  }
  
  
  
  #determining main somies and proportion of cells with main somies:
  
  main_somy <- c()
  prop_with_main_somy <- c()
  
  if(decluster == TRUE){
    print("declustering cell nodes")
    x <- karyo_decluster(x2, TRUE, TRUE, TRUE)
  }
  
  #return(x)
  
  x <- numeric.matrix(x)
  x <- round(x)
  ncell <- nrow(x)
  
  print("Getting proportion of cells displaying main somies...")
  
  #return(x)
  
  for(i in c(1:ncol(x))){
    main <- Mode(x[,i])
    prop <- sum(x[,i] == main)
    print(paste("Chromo", i , "has",  prop, "cells with main somy."))
    prop <- prop/ncell
    main_somy <- c(main_somy, main)
    prop_with_main_somy <- c(prop_with_main_somy, prop)
  }
  
  chromo_number <- nrow(somies)
  
  somies$Chromosome <- c(1:chromo_number)
  somies$main_somy <- main_somy
  somies$prop_main_somy <- prop_with_main_somy 
  somies <- arrange(somies, prop_with_main_somy)
  somies$variability_position <- c(1:nrow(somies))
  
  #return(somies)  
  
  dataframe <- gather(somies, key = karyo_id, value = Somy, -Chromosome, -main_somy, -prop_main_somy, -variability_position)
  dataframe$karyo_id <- as.numeric(dataframe$karyo_id)
  dataframe <- arrange(dataframe, karyo_id)
  
  
  Chromosome <- c()
  Occurrency <- c()
  Frequency <- c()
  diff_from_the_first <- c()
  number_of_nodes <- c()
  
  is_somy_different_from_first <- c()
  Marker <- c()
  
  print("Adding Karyotypes info to the dataframe...")
  
  for (i in c(1:nrow(karyotypes_info))){
    for(j in c(1:chromo_number)){
      Occurrency <- c(Occurrency ,karyotypes_info[i,2])
      Frequency <- c(Frequency ,karyotypes_info[i,4])
      diff_from_the_first <- c(diff_from_the_first ,karyotypes_info[i,5])
      number_of_nodes <- c(number_of_nodes, karyotypes_info[i,6])
    }
  }
  
  #dataframe <- arrange(dataframe, karyo_id)
  for(karyotype in c(0:(nrow(karyotypes_info)-1))){
    for (chromosome in as.numeric(levels(as.factor(dataframe$Chromosome)))){
      is_somy_different_from_first <- c(is_somy_different_from_first,dataframe[chromosome,6] != dataframe[chromosome+(36*karyotype),6])
    }
  }
  
  
  print("Finalizing...")
  
  if(decluster == TRUE){
    dataframe$karyo_id <- as.numeric(dataframe$karyo_id)
    dataframe$number_of_nodes <- number_of_nodes
  }
  
  
  dataframe$Somy <- as.numeric(dataframe$Somy)
  dataframe$number_of_cells <- Occurrency
  dataframe$Frequency <- Frequency
  dataframe$diff_from_the_first <- diff_from_the_first
  dataframe$is_somy_different_from_first <- is_somy_different_from_first
  
  number <- which(colnames(dataframe) == "is_somy_different_from_first")
  
  for(i in c(1:nrow(dataframe))){
    if(dataframe[i,number] == TRUE){
      Marker <- c(Marker, "*")
    }else{
      Marker <- c(Marker, NA)
    }
  }
  dataframe$Marker <- Marker
  
  dataframe <- arrange(dataframe, karyo_id, Chromosome)
  print("Done!")
  return(dataframe)
}

#karyo_calc_prop will calculate the proportions of somies for each chromosome and give them the stability_position. If gather is set to true, the output will be ready for ggplot.
karyo_calc_prop <- function(x, gather = FALSE, decluster = TRUE){
  
  if(decluster == TRUE){
    x <- karyo_decluster(x, TRUE, TRUE, TRUE)
  }
  
  x$barcode <- NULL
  x$node_id <- NULL
  x$barcodes <- NULL
  x$num_cells <- NULL
  x$num_noisy <- NULL
  
  
  karyo_info <- karyo_uniques(x, decluster = FALSE, add_node_id = FALSE)
  
  
  x <- numeric.matrix(x)
  x <- round(x)
  x <- t(x)
  
  #return(x)
  
  min_somy <- min(x)
  max_somy <- max(x)
  somies <- c(min_somy:max_somy)
  num_of_cells <- ncol(x)
  num_of_chromo <- nrow(x)
  dataframe <- tibble(c(1:num_of_chromo))
  
  
  
  print(paste("min_somy is", min_somy))
  print(paste("max_somy is", max_somy))
  print(paste("num_of_cells is", num_of_cells))
  print(paste("num_of_chromo", num_of_chromo))
  
  
  
  for (somy in somies){
    vec <- c()
    for (chromosome in c(1:nrow(x))){
      vec <- c(vec, sum(x[chromosome,] == somy)/num_of_cells)
    }
    dataframe <- cbind(dataframe, vec)
  }
  colnames(dataframe) <- c("Chromosome", somies)
  #return(dataframe)
  
  main_somy <- c()
  for (i in c(1:nrow(dataframe))){
    main_somy <- c(main_somy, which(dataframe[i,2:6] == max(dataframe[i,2:6])))
  }
  
  dataframe$main_somy <- main_somy
  
  prop_main_somy <- c()
  for (chromosome in c(1:nrow(dataframe))){
    main <- dataframe[chromosome,ncol(dataframe)]
    prop <- dataframe[chromosome,main+1]
    prop_main_somy <- c(prop_main_somy, prop)
  }
  
  dataframe <- cbind(dataframe, prop_main_somy)
  dataframe <- arrange(dataframe, prop_main_somy)
  dataframe$stability_position <- rev(c(1:nrow(dataframe)))
  dataframe <- arrange(dataframe, Chromosome)
  dataframe$Chromosome <- factor(dataframe$Chromosome)
  dataframe$Chromosome <- factor(dataframe$Chromosome, levels = rev(levels(factor(dataframe$Chromosome))))
  
  if(gather == TRUE){
    dataframe <- gather(dataframe, key = "Somy", value = "Proportion", -Chromosome, -main_somy, -prop_main_somy, -stability_position)
    dataframe <- dataframe %>%
      arrange(stability_position, Proportion) %>%
      mutate(Order = factor(paste(stability_position, Somy), levels = paste(stability_position, Somy)))
  }
  
  
  return(dataframe)
}



####################################OTHER FUNCTIONS################################
#scgs_remove_nonCNV_bins
#This function will remove remove bins in the reads/bin file that were not used by cellranger to estimate CNVs.
#These bins display a CNV value = 128, so any bin composed by 128-only values in the cnv_data will be removed in the reads data.
#It needs the cnv_data file to know which bins to remove

scgs_remove_nonCNV_bins <- function(reads, cnv_data){   
  #firt part will compare both datas to check if they have the same number of columns (meaning the same bins)
  if(ncol(cnv_data) != ncol(reads)){
    stop("cnv_data and reads have different number of collumns.")
  }
  
  #then it will check if column names are equal in both data
  if(paste(colnames(cnv_data), collapse = "") != paste(colnames(reads), collapse = "")){
    stop("cnv_data and reads have different collumn names.")
  }
  
  #now it will determine which bins were used to calculate cnv. Non-used bins will have only CNV values = 128
  
  print("removing bins that were not used for CNV estimation by Cell Ranger DNA...")
  
  print(paste("reads had", ncol(reads), "columns"))
  bin_remove_list <- c()
  
  for(i in c(1:ncol(cnv_data))){
    
    if(length(which(cnv_data[,i] != 128)) == 0) {#here it is checking which bins have only 128 values
      bin_remove_list <- c(bin_remove_list, i) #here it will save the bins (columns) that will be removed (because they were not used for CNV calculation)
    }
    
  }
  
  reads <- reads[,-bin_remove_list]
  print(paste("now it has", ncol(reads), "columns"))
  print(paste(length(bin_remove_list), "bins were removed"))
  return(reads)
}


#scgs_gather_bin_depth:
#This function returns a gathered data frame with the number of reads per bin. This is usefull for ggplot. 

scgs_gather_bin_depth <- function(reads){
  
  #first it will get which are the bins of each chromosome
  bins_position <- scgs_find_bins(reads, FALSE)
  last_chromo <- max(bins_position$chromosome) #determines which is the last chromosome
  
  #getting the cell ids:
  if("cell" %in% colnames(reads)){
    cell_ids <- reads$cell
    cell_ids <- unlist(cell_ids) 
  }else{
    if("X1" %in% colnames(reads)){
      cell_ids <- reads$X1
      cell_ids <- unlist(cell_ids)
    }else{
      cell_ids <- c(1:nrow(reads))-1
      print(paste("assuming that cell_ids go from", cell_ids[1], "to", cell_ids[length(cell_ids)]))
    }
  }
  read_count <- cell_ids
  chromo_center <- c()
  #now, for each chromosome, it will bind the read count of the bins of that chromosome. This is done to remove any other type of information
  for(i in c(1:last_chromo)){
    first_bin <- bins_position$chromo_first_bin[i]
    last_bin <- bins_position$chromo_last_bin[i]
    read_count <- cbind(read_count, reads[,first_bin:last_bin])
    chromo_center <- c(chromo_center, median(c(first_bin, last_bin)))
  }
  
  
  #gathering the data
  reads <- read_count
  flank <- c("Ld")
  
  colnames(reads)[1] <- "cell" #since the read_count started with a vector containing the cell numbers, now it just need to chance the name of the first column to "cell
  reads <- gather(reads, key = "bin", value = "reads_per_bin", -cell)
  bin_num <- c(1:length(unique(reads$bin)))
  bin_num <- rep(bin_num-1, each = length(unique(reads$cell)))
  reads$chromosome <- sub("_.*", "", reads$bin)
  reads$chromosome <- sub(flank[1], "", reads$chromosome)
  reads$chromosome <- as.numeric(reads$chromosome)
  reads$chromo_center <- chromo_center[reads$chromosome]
  reads$bin_num <- bin_num
  return(reads)
}



################################SMALL FUNCTIONS####################################
#function to quickly save a daframe to current directory:
save_file <- function(x){
  file_name <-deparse(substitute(x))
  write.csv(x, paste(getwd(),"/", file_name, ".csv", sep = ""))
}

#small function that converts dataframes into numeric matrix
numeric.matrix <- function(x){
  x <- apply(as.matrix.noquote(x),2,as.numeric)
  return(x)
}

#Function to order a dataframe based on a Column adding a "Position" column  
df_order <- function(x, colnumber = 3, reverse = TRUE){
  
  if (reverse == TRUE){
    x <- x[rev(order(x[,colnumber])),]
  }else{
    x <- x[order(x[,colnumber]),]
  }
  
  x$position <- c(1:nrow(x))
  
  
  return(x)
}

#karyo_somies
#Function to creat a dataframe with the rounded somy values of all karyotypes

karyo_somies <- function(x, round_method = 1, decluster = TRUE){
  
  print("karyo_somies is running")
  
  y <- karyo_uniques(x, round_method, decluster, add_node_id =  FALSE)
  
  top <- nrow(y)
  
  df <- vector()
  
  for (i in c(1:top)){
    a <- as.vector(y[i,1])
    a <- strsplit(a, "_")
    df <- c(df, a)
  }
  df <- data.frame(df)
  names(df) <- c(1:ncol(df))
  
  return(df)
}

#small function usefull for plots of gaussians
fun_prop <-function(x, mean, sd, proportion){
  proportion * dnorm(x = x, mean = mean, sd = sd)
}

#this function will return the most frequent value in a vector (the mode)
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


#this function will create a table with the karyotypes 
scgs_karyo_table <- function(x){
  karyotype <- unique(x$karyotypes)
  karyo_occurrence <- unique(x$karyotypes)
  
}

integer_dist <- function(x){
  x <- sum(abs(x-round(x)))/(length(x)/2)
  
  return(x)
}

integer_furthest <- function(x){
  x <- abs(x - round(x))
  x <- max(x)
  return(x)
}


gini <- function (x, weights = rep(1, length = length(x)), na.rm = FALSE){
  x <- unlist(x)
  
  if (na.rm) {
    k <- is.na(x)
    if (any(k)) 
      x <- x[!k]
  }
  
  ox <- order(x)
  x <- x[ox]
  weights <- weights[ox]/sum(weights)
  p <- cumsum(weights)
  nu <- cumsum(weights * x)
  n <- length(nu)
  nu <- nu/nu[n]
  sum(nu[-1] * p[-n]) - sum(nu[-n] * p[-1])
}

fix_seq <-function(start, space, times){
  
  vec <- c(start)
  for(i in c(1:times)){
    vec <- c(vec, (vec[i] + space))
  }
  return(vec)
}


################################BETA Functions#######################################

#this is a function still in Beta. It try to calculate somies based on a Poisson mixture model of the read depth of each bin
scgs_calc_somy_MM <- function(x, cell = 1153){
  
  
  bins_position <- scgs_find_bins(x, FALSE)
  first_bin <- bins_position$chromo_first_bin[1]
  last_chromo <- nrow(bins_position)
  last_bin <- bins_position$chromo_last_bin[last_chromo]
  
  cell_id <- x[,-c(first_bin:last_bin)]
  data <- x[,first_bin:last_bin]
  
  
  data <- t(data)
  data <- data.frame(data)
  colnames(data) <- cell_id
  
  cell_model <- stepFlexmix(data[,cell] ~ 1, k = 2:6, nrep = 3, model = FLXMCmvpois())
  cell_model <- getModel(cell_model)
  
  number_of_clusters <- length(prior(cell_model))
  proportions <- prior(cell_model)
  lambdas <- parameters(cell_model)
  
  vec <- c()
  
  for(i in c(1:number_of_clusters)){
    vec <- c(vec, lambdas[[i]])
  }
  
  lambdas <- vec
  
  min_cluster <- which(lambdas == min(lambdas))
  
  normalized_lambdas <- lambdas/lambdas[min_cluster]
  
  
  
  cluster_summary <- data.frame(cluster = c(1:number_of_clusters), proportion = proportions, lambdas = lambdas, normalized_lambdas = normalized_lambdas)
  
  x <- clusters(cell_model)
  
  for(i in c(1:number_of_clusters)){
    x[x == i] <- cluster_summary$normalized_lambdas[i]
  }
  
  somies <- c()
  for(i in c(1:last_chromo)){
    first_bin <- bins_position$chromo_first_bin[i]
    last_bin <- bins_position$chromo_last_bin[i]
    somies <- c(somies, Mode(x[first_bin:last_bin]))
  }
  
  somies <- somies * scgs_find_scale_factor(somies, range = c(2:2.5))
  return(somies)
}




