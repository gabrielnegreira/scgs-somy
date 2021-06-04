#inputs
a <- BPK282_raw_somies
b <- BPK081_raw_somies
c <- MixPop_raw_somies

a$sample <- "BPK282"
b$sample <- "BPK081"
c$sample <- "Super-mosaic"

#parameters
method = "glm"
quartiles <- c(0.01, 0.99)

#merging data
common_collumns <- Reduce(intersect, list(colnames(a), colnames(b), colnames(c)))

a <- a[,which(colnames(a) %in% common_collumns)]
b <- b[,which(colnames(b) %in% common_collumns)]
c <- c[,which(colnames(c) %in% common_collumns)]

to_plot <- rbind(a, b, c)

p1 <- ggplot(to_plot, aes(y = effective_depth_of_coverage, x = baseline_ploidy, group = baseline_ploidy))+
  geom_jitter(width = 0.1, height = 0, alpha = 0.4, size = 0.5)+
  geom_violin(alpha = 0.4)+
  #geom_hline(yintercept = 0.45, color = "red")+
  scale_x_continuous(breaks = c(1:200))+
  scale_y_continuous(limits = c(0,1.5), breaks = c(1:200))+
  facet_grid(cols = vars(factor(sample, levels = c("BPK282", "BPK081", "Super-mosaic"))))
p2 <- ggplot(to_plot, aes(y = chromo_var, x = baseline_ploidy, group = baseline_ploidy))+
  geom_jitter(width = 0.1, height = 0, alpha = 0.4, size = 0.5)+
  geom_violin(alpha = 0.4)+
  #geom_hline(yintercept = 1.7, color = "red")+
  scale_x_continuous(breaks = c(1:200))+
  scale_y_continuous(limits = c(0,3), breaks = c(1:200))+
  facet_grid(cols = vars(factor(sample, levels = c("BPK282", "BPK081", "Super-mosaic"))))

plot_grid(p1,p2, nrow = 2)

p3 <- ggplot(to_plot, aes(y = chromo_var, x = effective_depth_of_coverage))+
          #geom_point(alpha = 0.4, aes(color = factor(sample, levels = c("BPK282", "BPK081", "Super-mosaic"))))+
          geom_point(alpha = 0.4, size = 0.5, aes(color = removed))+
          geom_smooth(method = method, formula = y ~ log(x))+
          #geom_vline(xintercept = 0.1, color = "red")+
          scale_color_discrete(name = "Sample")+
          scale_x_continuous(name = "Effective Depth of Coverage (x)", 
                             #breaks = seq(0, 200, length = 401), 
                             limits = quantile(to_plot$effective_depth_of_coverage, quartiles)
                             )+
          scale_y_continuous(name = "ICV score",
                             limits = quantile(to_plot$chromo_var, quartiles),  
                             #breaks = round(seq(0, 3, length = 16), 2)
                             )+
          scale_color_manual(values = c("black", "red"))+
          facet_grid(cols = vars(factor(sample, levels = c("BPK282", "BPK081", "Super-mosaic"))))

p4 <- ggplot(to_plot, aes(y = scaled_distance_to_integers, x = effective_depth_of_coverage))+
  #geom_point(alpha = 0.4, aes(color = factor(sample, levels = c("BPK282", "BPK081", "Super-mosaic"))))+
  geom_point(alpha = 0.4, size = 0.5, aes(color = removed))+
  geom_smooth(method = method, formula = y ~ log(x))+
  #geom_vline(xintercept = 0.1, color = "red")+
  scale_x_continuous(name = "Effective Depth of Coverage (x)", 
                    # breaks = seq(0, 200, length = 401), 
                      limits = quantile(to_plot$effective_depth_of_coverage, quartiles)
                     )+
  scale_y_continuous(name = "Mean distance to integers",  
                     #limits = c(0,0.5),
                     #breaks = round(seq(0, 0.5, length = 11), 2)
                     )+
  scale_color_manual(values = c("black", "red"))+
  facet_grid(cols = vars(factor(sample, levels = c("BPK282", "BPK081", "Super-mosaic"))))

p5 <- ggplot(to_plot, aes(x = factor(baseline_ploidy), y = chromo_var))+
  geom_jitter(width = 0.3, heigth = 0, alpha = 0.5, size = 0.5, aes(color = removed))+
  geom_violin(alpha = 0.3)+
  scale_y_continuous(name = "ICV-Score", limits = c(1,3))+
  scale_x_discrete(name = "Baseline Ploidy")+
  scale_color_manual(values = c("black", "red"))+
  facet_grid(cols = vars(factor(sample, levels = c("BPK282", "BPK081", "Super-mosaic"))))

plot_grid(p3,p4, p5, nrow = 3)
