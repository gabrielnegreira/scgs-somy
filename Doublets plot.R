#inputs
doublets <- c("AATCCAGTCTGGTGGC-1", "ACATACGGTACCACTA-1", "ACTGCTCCAGGTCTGC-1", "AGCAGCCCACGTCGGT-1", "AGTGAGGGTGGCCACT-1", "ATCCGAAGTGCGAACA-1", "CACAAACTCATTCGCC-1", "CACCACTAGTCACGCC-1", "CAGAATCTCGTTATCT-1", "CAGATCACAGTAGAAT-1", "CCCTCCTCAGTTTGCA-1", "CGTGTCTCATCGCACG-1", "CTCGTACGTTCGGTGC-1", "CTCGTACCAACTCACA-1", "CTGGTCTCATCCCTTG-1", "GACGTTAAGATCTGAA-1", "GGCGACTAGCAGGTCA-1", "GGCTGGTAGTGATGAT-1", "TAGCCGGAGGCTCATT-1", "TGCACCTCAGGCCTTG-1", "TGGACGCGTGTGTATC-1")
population_raw <- MixPop_raw_somies
population_integers <- MixPop_integer_somies
columns <- c("cell", as.character(c(1:36)), "barcode")

#getting doublets
to_plot <- filter(population_raw, barcode %in% doublets)
to_plot <- to_plot[,which(colnames(to_plot) %in% columns)]
to_plot <- to_plot %>%
  mutate(removed = if_else(cell %in% population_integers$cell, FALSE, TRUE))

cells_with_integers <- filter(population_integers, cell %in% to_plot$cell)

to_plot[which(to_plot$cell %in% cells_with_integers$cell), c(2:37)] <- cells_with_integers[which(cells_with_integers$cell %in% to_plot$cell), c(2:37)]
to_plot$karyo_occurrence <- NA

to_plot$karyo_occurrence[which(to_plot$cell %in% cells_with_integers$cell)] <- cells_with_integers$karyo_occurrence[which(cells_with_integers$cell %in% to_plot$cell)]

to_plot <- gather(to_plot, key = "chromosome", value = "somy", -cell, -barcode, -removed, -karyo_occurrence)

to_plot <- mutate(to_plot, label = ifelse(is.na(karyo_occurrence),
                                          paste("cell: ", cell, sep = ""),
                                          ifelse(karyo_occurrence > 2,
                                                 paste("karyotype found in \n", karyo_occurrence-1, " other cells", sep = ""),
                                                 ifelse(karyo_occurrence > 1,
                                                        paste("karyotype found in \n", karyo_occurrence-1, " other cell)", sep = ""),
                                                        paste("karyotype found only\n in this cell", sep = "")))))

to_plot$barcode <- reorder(to_plot$barcode, -to_plot$karyo_occurrence)
to_plot$barcode <- reorder(to_plot$barcode, to_plot$removed)


#to_plot <- mutate(to_plot, group = ifelse(!removed, "Non-removed (integer somies)", "Removed"))

to_plot <- mutate(to_plot, group = ifelse(removed, "removed", ifelse(karyo_occurrence > 2, "found in other cells", "new karyotypes")))


min_scale <- 0
max_scale <- 8

barcodes <- to_plot$barcode
barcodes <- unique(barcodes)

labels <- numeric(length(barcodes))

for(i in c(1:length(barcodes))){
  labels[i] <- unique(as.character(to_plot$label[which(to_plot$barcode == barcodes[i])]))
}

names(labels) <- barcodes 

#labels_to_plot <- unique(labels_to_plot)

ggplot(to_plot, aes(y = factor(chromosome, levels = c(36:1)), x = barcode, fill = somy, label = round(somy, 2) ))+
  geom_tile()+
  geom_text()+
  facet_grid(col = vars(group), scales = "free", space = "free")+
  scale_fill_gradientn(name = "Copy\nNumber", colors = heat_col, limits = c(0,8.25))+
  scale_x_discrete(labels = labels)+
  scale_y_discrete(name = "Chromosome")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), axis.title.x = element_blank(), legend.title = element_text(face = "bold"))



