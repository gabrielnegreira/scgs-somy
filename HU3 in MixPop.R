#loading tables:
library(readr)
library(readxl)

BGS_BPK475 <- read_delim("~/OneDrive - ITG/ITM/PhD/Results/WGS/2020-06-30 - WGS for somy estimation of 10X samples and barcoded BPK282/Pieter Analysis/BPK475 - 104123-001-021.mapq30.removedups.proper_paired.somy.csv", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)

BGS_BPK498 <- read_delim("~/OneDrive - ITG/ITM/PhD/Results/WGS/2020-06-30 - WGS for somy estimation of 10X samples and barcoded BPK282/Pieter Analysis/BPK498 - HG7FFDSXY_104123-001-022_GGCAAGTT-CTCAGAAG.mapq30.removedups.proper_paired.somy.csv", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)

BGS_HU3 <- read_delim("~/OneDrive - ITG/ITM/PhD/Results/WGS/2020-06-30 - WGS for somy estimation of 10X samples and barcoded BPK282/Pieter Analysis/HU3 - HG7FFDSXY_104123-001-019_TGGAAGCA-GAGAGTAC.mapq30.removedups.proper_paired.somy.csv", 
                      "\t", escape_double = FALSE, trim_ws = TRUE)

BGS_BPK506 <- read_delim("~/OneDrive - ITG/ITM/PhD/Results/WGS/2020-06-30 - WGS for somy estimation of 10X samples and barcoded BPK282/Pieter Analysis/BPK506 - HG7FFDSXY_104123-001-020_CAATAGCC-TCTACGCA.mapq30.removedups.proper_paired.somy.csv", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)

BGS_BPK081 <- read_delim("~/OneDrive - ITG/ITM/PhD/Results/WGS/2020-06-30 - WGS for somy estimation of 10X samples and barcoded BPK282/Pieter Analysis/BPK081 - HG7FFDSXY_104123-001-018_CACAGACT-CTGTACCA.mapq30.removedups.proper_paired.somy.csv", 
                         "\t", escape_double = FALSE, trim_ws = TRUE)

#getting the barcodes of the cells marked as HU3 based on SNPs
barcodes_HU3 <- barcodes_HU3 <- read_excel("inputs/barcodes_HU3.xlsx")

#gertting the barcodes of the doublets detected by the SNP method
doublets <- c("AATCCAGTCTGGTGGC-1", "ACATACGGTACCACTA-1", "ACTGCTCCAGGTCTGC-1", "AGCAGCCCACGTCGGT-1", "AGTGAGGGTGGCCACT-1", "ATCCGAAGTGCGAACA-1", "CACAAACTCATTCGCC-1", "CACCACTAGTCACGCC-1", "CAGAATCTCGTTATCT-1", "CAGATCACAGTAGAAT-1", "CCCTCCTCAGTTTGCA-1", "CGTGTCTCATCGCACG-1", "CTCGTACGTTCGGTGC-1", "CTCGTACCAACTCACA-1", "CTGGTCTCATCCCTTG-1", "GACGTTAAGATCTGAA-1", "GGCGACTAGCAGGTCA-1", "GGCTGGTAGTGATGAT-1", "TAGCCGGAGGCTCATT-1", "TGCACCTCAGGCCTTG-1", "TGGACGCGTGTGTATC-1")

#getting the BGS of all 204 strains from the eLife paper
df_204_strains <- read_excel("inputs/eLife_paper_donovani_strains.xlsx")
df_204_strains[,c(2:ncol(df_204_strains))] <- df_204_strains[,c(2:ncol(df_204_strains))]*2 #the values are normalized to 1. So here it multiplies everything by 2 to get the somy values.

BGS_BPK475 <- BGS_BPK475[-37,]
BGS_BPK498 <- BGS_BPK498[-37,]     
BGS_HU3 <- BGS_HU3[-37,]
BGS_BPK506 <- BGS_BPK506[-37,]
BGS_BPK081 <- BGS_BPK081[-37,]

#inputs
HU3_barcodes <- barcodes_HU3$barcode
MixPop <- MixPop_integer_somies
MixPop_data <- MixPop_metrics

#adding barcodes to integer somies
integer_cells_barcodes <- MixPop_data$barcode[which(MixPop_data$cell_id %in% MixPop$cell)]
MixPop$barcode <- integer_cells_barcodes

#checking which cells are HU3 and doublets
MixPop <- MixPop %>%
  mutate(is_HU3 = ifelse(barcode %in% HU3_barcodes, TRUE, FALSE)) %>%
  mutate(is_doublet = ifelse(barcode %in% doublets, TRUE, FALSE))

MixPop_backup <- MixPop
MixPop <- filter(MixPop, baseline_ploidy == 2)


#clustering the cells and dividing in 4 groups
distances <- dist(MixPop[,2:37])
clusters <- hclust(distances, method = "ward.D2")
MixPop$cluster <- cutree(clusters, k = 4)
group_vec <- c("A", "D", "B", "C")
MixPop <- mutate(MixPop, group = group_vec[cluster])



MixPop_HU3 <- filter(MixPop, barcode %in% HU3_barcodes)

MixPop_A <- filter(MixPop, group == "A")
MixPop_B <- filter(MixPop, group == "B")
MixPop_C <- filter(MixPop, group == "C")
MixPop_D <- filter(MixPop, group == "D")

#ploting
ggplot(MixPop, aes(x = group, y = 1, color = is_HU3))+
  geom_jitter(height = 0.5, size = 1.5, alpha = 0.75)+
  scale_color_manual(name = "SNP Profile", labels = c("BPK", "HU3"), values = c("#383838", "orange")) +
  scale_x_discrete(name = "Cluster")+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.x = element_text(size = 14, face = 2))


BGSVsSCGS_HU3 <- data.frame(chromosome = c(1:36), BGS = BGS_HU3$somy_bin20000, SCGS = colMeans(MixPop_C[,2:37]), sample = "HU3")
BGSVsSCGS_BPK475 <- data.frame(chromosome = c(1:36), BGS = BGS_BPK475$somy_bin20000, SCGS = colMeans(MixPop_A[,2:37]), sample = "BPK475")
BGSVsSCGS_BPK498 <- data.frame(chromosome = c(1:36), BGS = BGS_BPK498$somy_bin20000, SCGS = colMeans(MixPop_B[,2:37]), sample = "BPK498")
BGSVsSCGS_BPK506 <- data.frame(chromosome = c(1:36), BGS = BGS_BPK506$somy_bin20000, SCGS = colMeans(MixPop_D[,2:37]), sample = "BPK506")

to_plot <- rbind(BGSVsSCGS_HU3, BGSVsSCGS_BPK475, BGSVsSCGS_BPK498, BGSVsSCGS_BPK506)
to_plot2 <- cbind(BGSVsSCGS_HU3, BGSVsSCGS_BPK475, BGSVsSCGS_BPK498, BGSVsSCGS_BPK506)
to_plot <- gather(to_plot, key = "method", value = "somy", -sample, - chromosome)
to_plot <- to_plot %>%
  mutate(sample_label = ifelse(sample == "BPK475", "BPK475/Cluster A", ifelse(sample == "BPK498", "BPK498/Cluster B", ifelse(sample == "HU3", "HU3/Cluster C", "BPK506/Cluster D"))))

ggplot(to_plot, aes(x = factor(chromosome), y = somy, fill = method, facets = sample))+
  geom_col(position = "dodge", color = "black")+
  scale_x_discrete(name = "Chromosome") +
  scale_fill_manual(name = "Data", values = c("black", "#68afc4"), labels = str_wrap(c("Strain Average Somies (BGS)", "Cluster Cumulative Somies (SCGS)"), 20))+
  facet_grid(rows = vars(factor(sample_label, levels = c("BPK475/Cluster A", "BPK498/Cluster B", "HU3/Cluster C", "BPK506/Cluster D"))))+
  theme(axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15), 
        legend.title = element_text(size = 15, face = 2), 
        legend.text = element_text(size = 15, margin = margin(t = 10)), 
        legend.key = element_rect(size = 3, color = "white"), 
        legend.key.height = unit(20, "pt"), 
        legend.spacing.y = unit(15, "pt"), 
        strip.text = element_text(size = 13, face = 2))
  
df_compare <- data.frame(BPK081_old = df_204_strains$BPK081A1, 
                         BPK081_new = BGS_BPK081$somy_bin20000, 
                         BPK081_SCGS = colMeans(BPK081_integer_somies[,2:37]), 
                         BPK475_old = df_204_strains$BPK475A1, 
                         BPK475_new = BGS_BPK475$somy_bin20000, 
                         cluster_A = BGSVsSCGS_BPK475$SCGS, 
                         BPK498_old = df_204_strains$BPK498A1, 
                         BPK498_new = BGS_BPK498$somy_bin20000, 
                         cluster_B = BGSVsSCGS_BPK498$SCGS, 
                         HU3_old = df_204_strains$LdonLV9, 
                         HU3_new = BGS_HU3$somy_bin20000, 
                         cluster_C = BGSVsSCGS_HU3$SCGS, 
                         BPK506_old = df_204_strains$BPK506A1, 
                         BPK506_new = BGS_BPK506$somy_bin20000, 
                         cluster_D = BGSVsSCGS_BPK506$SCGS)

to_plot <- numeric.matrix(df_compare)

par(bg = "white")
heatmap.2(to_plot, key.xlab = "Somy",
          Rowv = F, 
          Colv = F,
          trace="none", 
          density = "density",
          colsep = c(3,6,9,12),
          linecol = 1,
          col= brewer_palette(200),
          #col = heat_palette(200),
          breaks = seq(min_scale,max_scale, length = 201), 
          dendrogram = "both",
          margins = c(7, 3),
          cexCol=1,
          key.title = "Color Key", 
          hclust=function(x) hclust(x,method = "ward.D2")
)


#HU3 SCGS Vs BGS
to_plot <- data.frame(chromosome = c(1:36), HU3_BGS = BGS_HU3$somy_bin20000, HU3_SCGS = colMeans(HU3_integer_somies[,2:37]))
to_plot <- gather(to_plot, key = "sample", value = "somy", -chromosome)

ggplot(to_plot, aes(x = factor(chromosome), y = somy, fill = sample))+
  geom_col(position = "dodge")+
  scale_fill_manual(values = c("black", "#68afc4"))
    

MixPop_A$sample <- "Cluster A"
MixPop_B$sample <- "Cluster B"
MixPop_C$sample <- "Cluster C"
MixPop_D$sample <- "Cluster D"
BPK081_integer_somies$sample <- "BPK081"
BPK282_integer_somies$sample <- "BPK282"

equalmix <- rbind(MixPop_A[sample(c(1:nrow(MixPop_A)),  size = 250),which(colnames(MixPop_A) %in% colnames(BPK282_integer_somies))],
                  MixPop_B[sample(c(1:nrow(MixPop_B)),  size = 250),which(colnames(MixPop_A) %in% colnames(BPK282_integer_somies))],
                  MixPop_C[sample(c(1:nrow(MixPop_C)),  size = 250),which(colnames(MixPop_A) %in% colnames(BPK282_integer_somies))],
                  MixPop_D[sample(c(1:nrow(MixPop_D)),  size = 250),which(colnames(MixPop_A) %in% colnames(BPK282_integer_somies))],
                  BPK081_integer_somies[sample(c(1:nrow(BPK081_integer_somies)),  size = 250),],
                  BPK282_integer_somies[sample(c(1:nrow(BPK282_integer_somies)),  size = 250),]
                  )


MixPop_backup <- MixPop_backup %>%
  mutate(color = (ifelse(is_doublet == TRUE, "purple", 
                         ifelse(is_HU3 == TRUE, "orange", "black"))))

to_cluster <- MixPop_backup
to_cluster <- hclust(dist(to_cluster[,2:37], method = "euclidean"), method = "ward.D2")
colors <- MixPop_backup$color[to_cluster$order]

to_plot <- MixPop_backup[to_cluster$order,]
to_plot <- to_plot[,-1]
to_plot <- to_plot[,1:36]
to_plot <- numeric.matrix(to_plot)
to_plot <- t(to_plot)
min_scale <- 0
max_scale <- round(max(to_plot))
max_scale <- 8

heatmap.2(to_plot,
          main = paste(ncol(to_plot), "Cells from 4 strains"),
          key.xlab = "Somy",
          density = "density", 
          denscol = "black",
          keysize = 1,
          key.title = "Color Key", 
          Rowv = F, 
          #Colv = as.dendrogram(to_cluster),
          trace="none", 
          #key = F,
          linecol = 1,
          col= color_palette(200),
          ColSideColors = colors,
          #col = heat_palette(200),
          #adjCol = c(0.5,0.5),
          breaks = seq(min_scale,max_scale, length = 201), 
          dendrogram = "both",
          cexRow=1.8, 
          labCol = FALSE,
          #labRow = FALSE,
          hclust = function(x) hclust(x, method = "ward.D2")
)



