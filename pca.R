suppressMessages(library(tidyverse))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))

################################
######### Load in Data #########
################################

load("../results/umapdf.RData")

# perform PCA
gexpr_df <- as.matrix(gexpr_df)
gexpr_pca <- prcomp(gexpr_df)

isoform_df <- as.matrix(isoform_df)
isoform_pca <- prcomp(isoform_df)

ciri_df <- as.matrix(ciri_df)
ciri_pca <- prcomp(ciri_df)

circ_df <- as.matrix(circ_df)
circ_pca <- prcomp(circ_df)

lung_df <- as.matrix(lung_df)
lung_pca <- prcomp(lung_df)


######################################
### get eigenvalues for scree plot ###
######################################

gexpr_eigenvalues <- gexpr_pca$sdev^2
gexpr_eigenvalues <- tibble(PC = factor(1:length(gexpr_eigenvalues)), 
                              variance = gexpr_eigenvalues) %>% 
  mutate(pct = variance/sum(variance)*100) %>% 
  mutate(pct_cum = cumsum(pct))
gexpr_eigenvalues <- as_tibble(gexpr_eigenvalues)

isoform_eigenvalues <- isoform_pca$sdev^2
isoform_eigenvalues <- tibble(PC = factor(1:length(isoform_eigenvalues)), 
                         variance = isoform_eigenvalues) %>% 
  mutate(pct = variance/sum(variance)*100) %>% 
  mutate(pct_cum = cumsum(pct))
isoform_eigenvalues <- as_tibble(isoform_eigenvalues)

ciri_eigenvalues <- ciri_pca$sdev^2
ciri_eigenvalues <- tibble(PC = factor(1:length(ciri_eigenvalues)), 
                              variance = ciri_eigenvalues) %>% 
  mutate(pct = variance/sum(variance)*100) %>% 
  mutate(pct_cum = cumsum(pct))
ciri_eigenvalues <- as_tibble(ciri_eigenvalues)

circ_eigenvalues <- circ_pca$sdev^2
circ_eigenvalues <- tibble(PC = factor(1:length(circ_eigenvalues)), 
                              variance = circ_eigenvalues) %>% 
  mutate(pct = variance/sum(variance)*100) %>% 
  mutate(pct_cum = cumsum(pct))
circ_eigenvalues <- as_tibble(circ_eigenvalues)

lung_eigenvalues <- lung_pca$sdev^2
lung_eigenvalues <- tibble(PC = factor(1:length(lung_eigenvalues)), 
                           variance = lung_eigenvalues) %>% 
  mutate(pct = variance/sum(variance)*100) %>% 
  mutate(pct_cum = cumsum(pct))
lung_eigenvalues <- as_tibble(lung_eigenvalues)

##################################
### Scree plot of first 20 PCs ###
##################################

png("../results/scree_plot.png", width=400, height=200, units='mm', res = 600, pointsize=80)
p1 <- gexpr_eigenvalues[1:20,] %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained", title = "Genes") +
  scale_y_continuous(limits = c(0, 75), expand=c(0,0)) +
  theme_classic()
p2 <- isoform_eigenvalues[1:20,] %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained", title = "Isoforms") +
  scale_y_continuous(limits = c(0, 65), expand=c(0,0)) +
  theme_classic()
p3 <- ciri_eigenvalues[1:18,] %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained", title = "CIRI2 circRNAs") +
  scale_y_continuous(limits = c(0, 105), expand=c(0,0)) +
  theme_classic()
p4 <- circ_eigenvalues[1:20,] %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained", title = "Circexplorer2 circRNAs") +
  scale_y_continuous(limits = c(0,105), expand=c(0,0)) +
  theme_classic()
ggarrange(p1, p2, p3, p4,
          ncol = 2, nrow = 2,
          common.legend = TRUE,
          legend = "right")
dev.off()


png("../results/scree_plot_lung.png", width=200, height=100, units='mm', res = 600, pointsize=80)
lung_eigenvalues[1:20,] %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained", title = "Lung Tumours") +
  scale_y_continuous(limits = c(0, 75), expand=c(0,0)) +
  theme_classic()
dev.off()


#####################
### get PC scores ###
#####################

gexpr_scores <- as_tibble(gexpr_pca$x)
gexpr_scores$cell_line <- cell_line_gexpr
gexpr_scores$dataset <- c(rep("gCSI", 48), rep("CCLE", 48), rep("GDSC", 48))

isoform_scores <- as_tibble(isoform_pca$x)
isoform_scores$cell_line <- cell_line_isoforms
isoform_scores$dataset <- c(rep("gCSI", 48), rep("CCLE", 48), rep("GDSC", 48))

ciri_scores <- as_tibble(ciri_pca$x)
ciri_scores$cell_line <- cell_line_ciri
ciri_scores$dataset <- c(rep("gCSI", 48), rep("CCLE", 48), rep("GDSC", 48))

circ_scores <- as_tibble(circ_pca$x)
circ_scores$cell_line <- cell_line_circ
circ_scores$dataset <- c(rep("gCSI", 48), rep("CCLE", 48), rep("GDSC", 48))

lung_scores <- as_tibble(lung_pca$x)
lung_scores$tumourID <- c(rep(paste0("tumor",1:51), 4))
lung_scores$Selection <- c(rep("poly(A)", 51), rep("RiboZero", 51), rep("poly(A)", 51), rep("RiboZero", 51))
lung_scores$Method <- c(rep("CIRI2", 102), rep("CIRCexplorer2", 102))

################
### plot PCs ###
################

png("../results/pca_gene_expression.png", width=300, height=200, units='mm', res = 600, pointsize=80)
p1 <- ggplot(data = gexpr_scores, aes(x = PC1, y = PC2, group = cell_line)) + 
  geom_point(aes(color = dataset, shape = dataset), size = 3) + 
  scale_color_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), labels=c("CCLE", "gCSI", "GDSC2"), values = c("#392C57", "#4CC5AB", "#3670A0")) + 
  guides(shape = guide_legend(ncol = 1)) + #make legend one column instead of 3
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "Gene Expression PC1&PC2") +
  scale_shape_discrete(name = "Dataset", labels = c("CCLE", "gCSI", "GDSC2"))
p2 <- ggplot(data = gexpr_scores, aes(x = PC2, y = PC3, group = cell_line)) + 
  geom_point(aes(color = dataset, shape = dataset), size = 3) + 
  scale_color_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), labels=c("CCLE", "gCSI", "GDSC2"), values = c("#392C57", "#4CC5AB", "#3670A0")) + 
  guides(shape = guide_legend(ncol = 1)) + #make legend one column instead of 3
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "Gene Expression PC2&PC3") +
  scale_shape_discrete(name = "Dataset", labels = c("CCLE", "gCSI", "GDSC2"))
p3 <- ggplot(data = gexpr_scores, aes(x = PC1, y = PC2, group = cell_line)) + 
  geom_line(show.legend = F) + 
  geom_point(aes(color = dataset, shape = dataset), size = 3) + 
  scale_color_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), labels=c("CCLE", "gCSI", "GDSC2"), values = c("#392C57", "#4CC5AB", "#3670A0")) + 
  guides(shape = guide_legend(ncol = 1)) + #make legend one column instead of 3
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "Gene Expression PC1&PC2") +
  scale_shape_discrete(name = "Dataset", labels = c("CCLE", "gCSI", "GDSC2"))
p4 <- ggplot(data = gexpr_scores, aes(x = PC2, y = PC3, group = cell_line)) + 
  geom_line(show.legend = F) + 
  geom_point(aes(color = dataset, shape = dataset), size = 3) + 
  scale_color_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), labels=c("CCLE", "gCSI", "GDSC2"), values = c("#392C57", "#4CC5AB", "#3670A0")) + 
  guides(shape = guide_legend(ncol = 1)) + #make legend one column instead of 3
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "Gene Expression PC2&PC3") +
  scale_shape_discrete(name = "Dataset", labels = c("CCLE", "gCSI", "GDSC2"))
ggarrange(p1, p2, p3, p4,
          ncol = 2, nrow = 2,
          common.legend = TRUE,
          legend = "bottom")
dev.off()


png("../results/pca_isoforms.png", width=300, height=200, units='mm', res = 600, pointsize=80)
p1 <- ggplot(data = isoform_scores, aes(x = PC1, y = PC2, group = cell_line)) + 
  geom_point(aes(color = dataset, shape = dataset), size = 3) + 
  scale_color_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), labels=c("CCLE", "gCSI", "GDSC2"), values = c("#392C57", "#4CC5AB", "#3670A0")) + 
  guides(shape = guide_legend(ncol = 1)) + #make legend one column instead of 3
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "Isoforms PC1&PC2") +
  scale_shape_discrete(name = "Dataset", labels = c("CCLE", "gCSI", "GDSC2"))
p2 <- ggplot(data = isoform_scores, aes(x = PC2, y = PC3, group = cell_line)) + 
  geom_point(aes(color = dataset, shape = dataset), size = 3) + 
  scale_color_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), labels=c("CCLE", "gCSI", "GDSC2"), values = c("#392C57", "#4CC5AB", "#3670A0")) + 
  guides(shape = guide_legend(ncol = 1)) + #make legend one column instead of 3
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "Isoforms PC2&PC3") +
  scale_shape_discrete(name = "Dataset", labels = c("CCLE", "gCSI", "GDSC2"))
p3 <- ggplot(data = isoform_scores, aes(x = PC1, y = PC2, group = cell_line)) + 
  geom_line(show.legend = F) + 
  geom_point(aes(color = dataset, shape = dataset), size = 3) + 
  scale_color_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), labels=c("CCLE", "gCSI", "GDSC2"), values = c("#392C57", "#4CC5AB", "#3670A0")) + 
  guides(shape = guide_legend(ncol = 1)) + #make legend one column instead of 3
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "Isoforms PC1&PC2") +
  scale_shape_discrete(name = "Dataset", labels = c("CCLE", "gCSI", "GDSC2"))
p4 <- ggplot(data = isoform_scores, aes(x = PC2, y = PC3, group = cell_line)) + 
  geom_line(show.legend = F) + 
  geom_point(aes(color = dataset, shape = dataset), size = 3) + 
  scale_color_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), labels=c("CCLE", "gCSI", "GDSC2"), values = c("#392C57", "#4CC5AB", "#3670A0")) + 
  guides(shape = guide_legend(ncol = 1)) + #make legend one column instead of 3
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "Isoforms PC2&PC3") +
  scale_shape_discrete(name = "Dataset", labels = c("CCLE", "gCSI", "GDSC2"))
ggarrange(p1, p2, p3, p4,
          ncol = 2, nrow = 2,
          common.legend = TRUE,
          legend = "bottom")
dev.off()


png("../results/pca_ciri.png", width=300, height=200, units='mm', res = 600, pointsize=80)
p1 <- ggplot(data = ciri_scores, aes(x = PC1, y = PC2, group = cell_line)) + 
  geom_point(aes(color = dataset, shape = dataset), size = 3) + 
  scale_color_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), labels=c("CCLE", "gCSI", "GDSC2"), values = c("#392C57", "#4CC5AB", "#3670A0")) + 
  guides(shape = guide_legend(ncol = 1)) + #make legend one column instead of 3
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRI2 PC1&PC2") +
  scale_shape_discrete(name = "Dataset", labels = c("CCLE", "gCSI", "GDSC2"))
p2 <- ggplot(data = ciri_scores, aes(x = PC2, y = PC3, group = cell_line)) + 
  geom_point(aes(color = dataset, shape = dataset), size = 3) + 
  scale_color_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), labels=c("CCLE", "gCSI", "GDSC2"), values = c("#392C57", "#4CC5AB", "#3670A0")) + 
  guides(shape = guide_legend(ncol = 1)) + #make legend one column instead of 3
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRI2 PC2&PC3") +
  scale_shape_discrete(name = "Dataset", labels = c("CCLE", "gCSI", "GDSC2"))
p3 <- ggplot(data = ciri_scores, aes(x = PC1, y = PC2, group = cell_line)) + 
  geom_line(show.legend = F) + 
  geom_point(aes(color = dataset, shape = dataset), size = 3) + 
  scale_color_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), labels=c("CCLE", "gCSI", "GDSC2"), values = c("#392C57", "#4CC5AB", "#3670A0")) + 
  guides(shape = guide_legend(ncol = 1)) + #make legend one column instead of 3
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRI2 PC1&PC2") +
  scale_shape_discrete(name = "Dataset", labels = c("CCLE", "gCSI", "GDSC2"))
p4 <- ggplot(data = ciri_scores, aes(x = PC2, y = PC3, group = cell_line)) + 
  geom_line(show.legend = F) + 
  geom_point(aes(color = dataset, shape = dataset), size = 3) + 
  scale_color_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), labels=c("CCLE", "gCSI", "GDSC2"), values = c("#392C57", "#4CC5AB", "#3670A0")) + 
  guides(shape = guide_legend(ncol = 1)) + #make legend one column instead of 3
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRI2 PC2&PC3") +
  scale_shape_discrete(name = "Dataset", labels = c("CCLE", "gCSI", "GDSC2"))
ggarrange(p1, p2, p3, p4,
          ncol = 2, nrow = 2,
          common.legend = TRUE,
          legend = "bottom")
dev.off()

png("../results/pca_circ.png", width=300, height=200, units='mm', res = 600, pointsize=80)
p1 <- ggplot(data = circ_scores, aes(x = PC1, y = PC2, group = cell_line)) + 
  geom_point(aes(color = dataset, shape = dataset), size = 3) + 
  scale_color_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), labels=c("CCLE", "gCSI", "GDSC2"), values = c("#392C57", "#4CC5AB", "#3670A0")) + 
  guides(shape = guide_legend(ncol = 1)) + #make legend one column instead of 3
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRCexplorer2 PC1&PC2") +
  scale_shape_discrete(name = "Dataset", labels = c("CCLE", "gCSI", "GDSC2"))
p2 <- ggplot(data = circ_scores, aes(x = PC2, y = PC3, group = cell_line)) + 
  geom_point(aes(color = dataset, shape = dataset), size = 3) + 
  scale_color_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), labels=c("CCLE", "gCSI", "GDSC2"), values = c("#392C57", "#4CC5AB", "#3670A0")) + 
  guides(shape = guide_legend(ncol = 1)) + #make legend one column instead of 3
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRCexplorer2 PC2&PC3") +
  scale_shape_discrete(name = "Dataset", labels = c("CCLE", "gCSI", "GDSC2"))
p3 <- ggplot(data = circ_scores, aes(x = PC1, y = PC2, group = cell_line)) + 
  geom_line(show.legend = F) + 
  geom_point(aes(color = dataset, shape = dataset), size = 3) + 
  scale_color_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), labels=c("CCLE", "gCSI", "GDSC2"), values = c("#392C57", "#4CC5AB", "#3670A0")) + 
  guides(shape = guide_legend(ncol = 1)) + #make legend one column instead of 3
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRCexplorer2 PC1&PC2") +
  scale_shape_discrete(name = "Dataset", labels = c("CCLE", "gCSI", "GDSC2"))
p4 <- ggplot(data = circ_scores, aes(x = PC2, y = PC3, group = cell_line)) + 
  geom_line(show.legend = F) + 
  geom_point(aes(color = dataset, shape = dataset), size = 3) + 
  scale_color_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), labels=c("CCLE", "gCSI", "GDSC2"), values = c("#392C57", "#4CC5AB", "#3670A0")) + 
  guides(shape = guide_legend(ncol = 1)) + #make legend one column instead of 3
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRCexplorer2 PC2&PC3") +
  scale_shape_discrete(name = "Dataset", labels = c("CCLE", "gCSI", "GDSC2"))
ggarrange(p1, p2, p3, p4,
          ncol = 2, nrow = 2,
          common.legend = TRUE,
          legend = "bottom")
dev.off()

png("../results/pca_lung.png", width=150, height=125, units='mm', res = 600, pointsize=80)
ggplot(data = lung_scores, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = Selection, shape = Method), size = 3) + 
  guides(shape = guide_legend(ncol = 1), color = guide_legend(override.aes=list(shape=15, size = 8))) + #make legend one column instead of 3, change shape to square
  scale_color_manual(guide = guide_legend(reverse = FALSE), labels=c("poly(A)", "RiboZero"), values = c("#4CC5AB", "#392C57")) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "Lung PC1&PC2") 
dev.off()


# tumour2 is an outlier (at least with polyA), redo with tumour2 removed
lung_scores <- lung_scores[-which(lung_scores$tumourID == "tumor2"),]

png("../results/pca_lung_outlier.png", width=150, height=125, units='mm', res = 600, pointsize=80)
ggplot(data = lung_scores, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = Selection, shape = Method), size = 3) + 
  guides(shape = guide_legend(ncol = 1), color = guide_legend(override.aes=list(shape=15, size = 8))) + #make legend one column instead of 3, change shape to square
  scale_color_manual(guide = guide_legend(reverse = FALSE), labels=c("poly(A)", "RiboZero"), values = c("#4CC5AB", "#392C57")) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "Lung PC1&PC2") 
dev.off()



###################################
### Identify circRNA cell lines ###
###################################


gdsc_ciri <- ciri_scores[which(ciri_scores$PC1 < -5),]$cell_line
ciri_scores_keep <- unique(ciri_scores[-which(ciri_scores$cell_line %in% gdsc_ciri),]$cell_line)
gdsc_circ <- circ_scores[which(circ_scores$PC1 > 5),]$cell_line
circ_scores_keep <- unique(circ_scores[-which(circ_scores$cell_line %in% gdsc_circ),]$cell_line)

#read in current cell annotations
cell_all <- read.csv(file = "../data/rnaseq_meta/cell_annotation_all.csv", na.strings=c("", " ", "NA"))

# Keep cell lines of interest common cell lines
ciri_annot <- cell_all[which(cell_all$X %in% ciri_scores_keep),]
ciri_annot <- ciri_annot[,c("X", "Cellosaurus.Disease.Type", "CCLE.tissueid")]

# Keep cell lines of GDSC2 cluster
gdsc_group <- unique(circ_scores[which(circ_scores$cell_line %in% gdsc_circ),]$cell_line)
gdsc_annot <- cell_all[which(cell_all$X %in% gdsc_group),]
gdsc_annot <- gdsc_annot[,c("X", "Cellosaurus.Disease.Type", "CCLE.tissueid")]

png("../results/pca_circRNAs_group.png", width=300, height=200, units='mm', res = 600, pointsize=80)
p1 <- ggplot(data = ciri_scores, aes(x = PC1, y = PC2, group = cell_line)) + 
  geom_point(aes(shape = check, color = dataset), size = 3) + 
  #scale_color_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), values = c("#392C57", "#4CC5AB")) + 
  guides(shape = guide_legend(ncol = 1)) + #make legend one column instead of 3
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRI2 PC1&PC2") #+
  #scale_shape_discrete(name = "Dataset", labels = c("CCLE", "gCSI", "GDSC2"))
p3 <- ggplot(data = ciri_scores_group, aes(x = PC1, y = PC2, group = cell_line)) + 
  geom_line(show.legend = F) + 
  geom_point(aes(color = dataset, shape = dataset), size = 3) + 
  scale_color_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), labels=c("CCLE", "gCSI", "GDSC2"), values = c("#392C57", "#4CC5AB", "#3670A0")) + 
  guides(shape = guide_legend(ncol = 1)) + #make legend one column instead of 3
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRI2 PC1&PC2") +
  scale_shape_discrete(name = "Dataset", labels = c("CCLE", "gCSI", "GDSC2"))
p4 <- ggplot(data = ciri_scores, aes(x = PC2, y = PC3, group = cell_line)) + 
  geom_line(show.legend = F) + 
  geom_point(aes(color = dataset, shape = dataset), size = 3) + 
  scale_color_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), labels=c("CCLE", "gCSI", "GDSC2"), values = c("#392C57", "#4CC5AB", "#3670A0")) + 
  guides(shape = guide_legend(ncol = 1)) + #make legend one column instead of 3
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRI2 PC2&PC3") +
  scale_shape_discrete(name = "Dataset", labels = c("CCLE", "gCSI", "GDSC2"))
ggarrange(p1, p2, p3, p4,
          ncol = 2, nrow = 2,
          common.legend = TRUE,
          legend = "bottom")
dev.off()



################
### Lung PCA ###
################

lung_polyA <- rbind(lung_df[1:51,], lung_df[103:153,])
lung_polyA <- prcomp(lung_polyA)

lung_ribo <- rbind(lung_df[52:102,], lung_df[154:204,])
lung_ribo <- prcomp(lung_ribo)


polyA_lung_eigenvalues <- lung_polyA$sdev^2
polyA_lung_eigenvalues <- tibble(PC = factor(1:length(polyA_lung_eigenvalues)), 
                           variance = polyA_lung_eigenvalues) %>% 
  mutate(pct = variance/sum(variance)*100) %>% 
  mutate(pct_cum = cumsum(pct))
polyA_lung_eigenvalues <- as_tibble(polyA_lung_eigenvalues)

ribo_lung_eigenvalues <- lung_ribo$sdev^2
ribo_lung_eigenvalues <- tibble(PC = factor(1:length(ribo_lung_eigenvalues)), 
                              variance = ribo_lung_eigenvalues) %>% 
  mutate(pct = variance/sum(variance)*100) %>% 
  mutate(pct_cum = cumsum(pct))
ribo_lung_eigenvalues <- as_tibble(ribo_lung_eigenvalues)



png("../results/scree_plot_lung2.png", width=200, height=100, units='mm', res = 600, pointsize=80)
p1 <- polyA_lung_eigenvalues[1:20,] %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained", title = "polyA Lung Tumours") +
  scale_y_continuous(limits = c(0, 58), expand=c(0,0)) +
  theme_classic()
p2 <- ribo_lung_eigenvalues[1:20,] %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained", title = "Ribo Lung Tumours") +
  scale_y_continuous(limits = c(0, 58), expand=c(0,0)) +
  theme_classic()
ggarrange(p1, p2, 
          ncol = 2, nrow = 1,
          common.legend = TRUE,
          legend = "bottom")
dev.off()




lung_polyA <- as_tibble(lung_polyA$x)
lung_ribo <- as_tibble(lung_ribo$x)

lung_polyA$tumourID <- c(rep(paste0("tumor",1:51), 2))
lung_polyA$method <- c(rep("CIRI2", 51), rep("CIRCexplorer2", 51))
lung_ribo$tumourID <- c(rep(paste0("tumor",1:51), 2))
lung_ribo$method <- c(rep("CIRI2", 51), rep("CIRCexplorer2", 51))



annotation_lung <- read.csv("../data/survival_data/IMPACT13list_09Jan2019_with_survival_info_55_ADC_original.csv")
#keep only IMPACT tumor samples (51 samples)
annotation_lung <- annotation_lung[-c(1,46,47,49),]
annotation_lung$helab_id <- c(paste0("tumor",1:51))

lung_polyA$sex <- c(rep(annotation_lung$sex, 2))
lung_polyA$stage <- c(rep(annotation_lung$pathStage, 2))
lung_polyA$stage <- gsub("A", "", lung_polyA$stage)
lung_polyA$stage <- gsub("B", "", lung_polyA$stage)

lung_ribo$sex <- c(rep(annotation_lung$sex, 2))
lung_ribo$stage <- c(rep(annotation_lung$pathStage, 2))
lung_ribo$stage <- gsub("A", "", lung_ribo$stage)
lung_ribo$stage <- gsub("B", "", lung_ribo$stage)

# remove outlier
lung_polyA <- lung_polyA[-which(lung_polyA$tumourID == "tumor2"),]
lung_ribo <- lung_ribo[-which(lung_ribo$tumourID == "tumor2"),]

png("../results/pca_lung2_outlier.png", width=450, height=200, units='mm', res = 600, pointsize=80)
p1 <- ggplot(data = lung_polyA, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = method), size = 3) + 
  guides(color = guide_legend(override.aes=list(shape=15, size = 8))) + #make legend one column instead of 3, change shape to square
  scale_color_manual(guide = guide_legend(reverse = FALSE), values = c("#4CC5AB", "#392C57")) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "polyA Lung PC1&PC2") 
p2 <- ggplot(data = lung_polyA, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = sex), size = 3) + 
  guides(color = guide_legend(override.aes=list(shape=15, size = 8))) + #make legend one column instead of 3, change shape to square
  scale_color_manual(guide = guide_legend(reverse = FALSE), values = c("#4CC5AB", "#392C57")) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "polyA Lung PC1&PC2") 
p3 <- ggplot(data = lung_polyA, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = stage), size = 3) + 
  guides(color = guide_legend(override.aes=list(shape=15, size = 8))) + #make legend one column instead of 3, change shape to square
  scale_color_manual(guide = guide_legend(reverse = FALSE), values = c("#4CC5AB", "#392C57", "darkgreen", "darkblue")) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "polyA Lung PC1&PC2") 
p4 <- ggplot(data = lung_ribo, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = method), size = 3) + 
  guides(color = guide_legend(override.aes=list(shape=15, size = 8))) + #make legend one column instead of 3, change shape to square
  scale_color_manual(guide = guide_legend(reverse = FALSE), values = c("#4CC5AB", "#392C57")) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "Ribo Lung PC1&PC2") 
p5 <- ggplot(data = lung_ribo, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = sex), size = 3) + 
  guides(color = guide_legend(override.aes=list(shape=15, size = 8))) + #make legend one column instead of 3, change shape to square
  scale_color_manual(guide = guide_legend(reverse = FALSE), values = c("#4CC5AB", "#392C57")) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "Ribo Lung PC1&PC2") 
p6 <- ggplot(data = lung_ribo, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = stage), size = 3) + 
  guides(color = guide_legend(override.aes=list(shape=15, size = 8))) + #make legend one column instead of 3, change shape to square
  scale_color_manual(guide = guide_legend(reverse = FALSE), values = c("#4CC5AB", "#392C57", "darkgreen", "darkblue")) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "Ribo Lung PC1&PC2") 

ggarrange(p1, p2, p3, p4, p5, p6,
          ncol = 3, nrow = 2)
dev.off()

##########################################
### Liquid vs Solid Tumours Cell Lines ###
##########################################

# Classify cell lines as liquid or solid
load("../results/inter.RData")

#read in current cell annotations
cell_all <- read.csv(file = "../data/rnaseq_meta/cell_annotation_all.csv", na.strings=c("", " ", "NA"))

# Keep 48 common cell lines
cell_lines <- cell_all[which(cell_all$X %in% intersected_rnacells),]
cell_lines <- cell_lines[,c("X", "Cellosaurus.Disease.Type", "CCLE.tissueid")]

# label liquid and solid tumours
cell_lines$tumour.type <- ifelse(cell_lines$CCLE.tissueid == "haematopoietic_and_lymphoid_tissue",
                                 "liquid", "solid")
cell_lines$tumour.type[is.na(cell_lines$tumour.type)] <- "solid"


# add labels to PCA scores
gexpr_scores$tumour.type <- cell_lines[match(gexpr_scores$cell_line, cell_lines$X),]$tumour.type
isoform_scores$tumour.type <- cell_lines[match(isoform_scores$cell_line, cell_lines$X),]$tumour.type
ciri_scores$tumour.type <- cell_lines[match(ciri_scores$cell_line, cell_lines$X),]$tumour.type
circ_scores$tumour.type <- cell_lines[match(circ_scores$cell_line, cell_lines$X),]$tumour.type

# add labels to original data frames
gexpr_df$cell_line <- cell_line_gexpr
isoform_df$cell_line <- cell_line_isoforms
ciri_df$cell_line <- cell_line_ciri
circ_df$cell_line <- cell_line_circ

gexpr_df$tumour.type <- cell_lines[match(gexpr_df$cell_line, cell_lines$X),]$tumour.type
isoform_df$tumour.type <- cell_lines[match(isoform_df$cell_line, cell_lines$X),]$tumour.type
ciri_df$tumour.type <- cell_lines[match(ciri_df$cell_line, cell_lines$X),]$tumour.type
circ_df$tumour.type <- cell_lines[match(circ_df$cell_line, cell_lines$X),]$tumour.type

gexpr_df <- gexpr_df[which(gexpr_df$tumour.type == "solid"),]
isoform_df <- isoform_df[which(isoform_df$tumour.type == "solid"),]
ciri_df <- ciri_df[which(ciri_df$tumour.type == "solid"),]
circ_df <- circ_df[which(circ_df$tumour.type == "solid"),]

# save new data frames
save(gexpr_df, isoform_df, ciri_df, circ_df, file = "../results/solid_tumours.RData")

## plot PCs ##
png("../results/pca_gene_expression_tumour_type.png", width=300, height=200, units='mm', res = 600, pointsize=80)
p1 <- ggplot(data = gexpr_scores, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = tumour.type), size = 3) + 
  guides(color = guide_legend(ncol = 1)) + #make legend one column instead of 3
  scale_color_manual(guide = guide_legend(reverse = FALSE), name = "Tumour Type", labels = c("Liquid", "Solid"), values = c("#392C57", "#4CC5AB")) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "Gene Expression PC1&PC2") 
p2 <- ggplot(data = gexpr_scores, aes(x = PC2, y = PC3)) + 
  geom_point(aes(color = tumour.type), size = 3) + 
  guides(color = guide_legend(ncol = 1)) + #make legend one column instead of 3
  scale_color_manual(guide = guide_legend(reverse = FALSE), name = "Tumour Type", labels = c("Liquid", "Solid"), values = c("#392C57", "#4CC5AB")) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "Gene Expression PC2&PC3") 
p3 <- ggplot(data = gexpr_scores, aes(x = PC1, y = PC2, group = cell_line)) + 
  geom_line(show.legend = F) + 
  geom_point(aes(color = tumour.type), size = 3) + 
  guides(color = guide_legend(ncol = 1)) + #make legend one column instead of 3
  scale_color_manual(guide = guide_legend(reverse = FALSE), name = "Tumour Type", labels = c("Liquid", "Solid"), values = c("#392C57", "#4CC5AB")) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "Gene Expression PC1&PC2") 
p4 <- ggplot(data = gexpr_scores, aes(x = PC2, y = PC3, group = cell_line)) + 
  geom_line(show.legend = F) + 
  geom_point(aes(color = tumour.type), size = 3) + 
  guides(color = guide_legend(ncol = 1)) + #make legend one column instead of 3
  scale_color_manual(guide = guide_legend(reverse = FALSE), name = "Tumour Type", labels = c("Liquid", "Solid"), values = c("#392C57", "#4CC5AB")) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "Gene Expression PC2&PC3") 
ggarrange(p1, p2, p3, p4,
          ncol = 2, nrow = 2,
          common.legend = TRUE,
          legend = "bottom")
dev.off()

png("../results/pca_isoforms_tumour_type.png", width=300, height=200, units='mm', res = 600, pointsize=80)
p1 <- ggplot(data = isoform_scores, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = tumour.type), size = 3) + 
  guides(color = guide_legend(ncol = 1)) + #make legend one column instead of 3
  scale_color_manual(guide = guide_legend(reverse = FALSE), name = "Tumour Type", labels = c("Liquid", "Solid"), values = c("#392C57", "#4CC5AB")) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "Isoforms PC1&PC2") 
p2 <- ggplot(data = isoform_scores, aes(x = PC2, y = PC3)) + 
  geom_point(aes(color = tumour.type), size = 3) + 
  guides(color = guide_legend(ncol = 1)) + #make legend one column instead of 3
  scale_color_manual(guide = guide_legend(reverse = FALSE), name = "Tumour Type", labels = c("Liquid", "Solid"), values = c("#392C57", "#4CC5AB")) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "Isoforms PC2&PC3") 
p3 <- ggplot(data = isoform_scores, aes(x = PC1, y = PC2, group = cell_line)) + 
  geom_line(show.legend = F) + 
  geom_point(aes(color = tumour.type), size = 3) + 
  guides(color = guide_legend(ncol = 1)) + #make legend one column instead of 3
  scale_color_manual(guide = guide_legend(reverse = FALSE), name = "Tumour Type", labels = c("Liquid", "Solid"), values = c("#392C57", "#4CC5AB")) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "Isoforms PC1&PC2") 
p4 <- ggplot(data = isoform_scores, aes(x = PC2, y = PC3, group = cell_line)) + 
  geom_line(show.legend = F) + 
  geom_point(aes(color = tumour.type), size = 3) + 
  guides(color = guide_legend(ncol = 1)) + #make legend one column instead of 3
  scale_color_manual(guide = guide_legend(reverse = FALSE), name = "Tumour Type", labels = c("Liquid", "Solid"), values = c("#392C57", "#4CC5AB")) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "Isoforms PC2&PC3") 
ggarrange(p1, p2, p3, p4,
          ncol = 2, nrow = 2,
          common.legend = TRUE,
          legend = "bottom")
dev.off()


png("../results/pca_ciri_tumour_type.png", width=300, height=200, units='mm', res = 600, pointsize=80)
p1 <- ggplot(data = ciri_scores, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = tumour.type), size = 3) + 
  guides(color = guide_legend(ncol = 1)) + #make legend one column instead of 3
  scale_color_manual(guide = guide_legend(reverse = FALSE), name = "Tumour Type", labels = c("Liquid", "Solid"), values = c("#392C57", "#4CC5AB")) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRI2 circRNA PC1&PC2") 
p2 <- ggplot(data = ciri_scores, aes(x = PC2, y = PC3)) + 
  geom_point(aes(color = tumour.type), size = 3) + 
  guides(color = guide_legend(ncol = 1)) + #make legend one column instead of 3
  scale_color_manual(guide = guide_legend(reverse = FALSE), name = "Tumour Type", labels = c("Liquid", "Solid"), values = c("#392C57", "#4CC5AB")) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRI2 circRNA PC2&PC3") 
p3 <- ggplot(data = ciri_scores, aes(x = PC1, y = PC2, group = cell_line)) + 
  geom_line(show.legend = F) + 
  geom_point(aes(color = tumour.type), size = 3) + 
  guides(color = guide_legend(ncol = 1)) + #make legend one column instead of 3
  scale_color_manual(guide = guide_legend(reverse = FALSE), name = "Tumour Type", labels = c("Liquid", "Solid"), values = c("#392C57", "#4CC5AB")) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRI2 circRNA PC1&PC2") 
p4 <- ggplot(data = ciri_scores, aes(x = PC2, y = PC3, group = cell_line)) + 
  geom_line(show.legend = F) + 
  geom_point(aes(color = tumour.type), size = 3) + 
  guides(color = guide_legend(ncol = 1)) + #make legend one column instead of 3
  scale_color_manual(guide = guide_legend(reverse = FALSE), name = "Tumour Type", labels = c("Liquid", "Solid"), values = c("#392C57", "#4CC5AB")) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRI2 circRNA PC2&PC3") 
ggarrange(p1, p2, p3, p4,
          ncol = 2, nrow = 2,
          common.legend = TRUE,
          legend = "bottom")
dev.off()

png("../results/pca_circ_tumour_type.png", width=300, height=200, units='mm', res = 600, pointsize=80)
p1 <- ggplot(data = circ_scores, aes(x = PC1, y = PC2)) + 
  geom_point(aes(color = tumour.type), size = 3) + 
  guides(color = guide_legend(ncol = 1)) + #make legend one column instead of 3
  scale_color_manual(guide = guide_legend(reverse = FALSE), name = "Tumour Type", labels = c("Liquid", "Solid"), values = c("#392C57", "#4CC5AB")) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRCexplorer2 circRNA PC1&PC2") 
p2 <- ggplot(data = circ_scores, aes(x = PC2, y = PC3)) + 
  geom_point(aes(color = tumour.type), size = 3) + 
  guides(color = guide_legend(ncol = 1)) + #make legend one column instead of 3
  scale_color_manual(guide = guide_legend(reverse = FALSE), name = "Tumour Type", labels = c("Liquid", "Solid"), values = c("#392C57", "#4CC5AB")) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRCexplorer2 circRNA PC2&PC3") 
p3 <- ggplot(data = circ_scores, aes(x = PC1, y = PC2, group = cell_line)) + 
  geom_line(show.legend = F) + 
  geom_point(aes(color = tumour.type), size = 3) + 
  guides(color = guide_legend(ncol = 1)) + #make legend one column instead of 3
  scale_color_manual(guide = guide_legend(reverse = FALSE), name = "Tumour Type", labels = c("Liquid", "Solid"), values = c("#392C57", "#4CC5AB")) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRCexplorer2 circRNA PC1&PC2") 
p4 <- ggplot(data = circ_scores, aes(x = PC2, y = PC3, group = cell_line)) + 
  geom_line(show.legend = F) + 
  geom_point(aes(color = tumour.type), size = 3) + 
  guides(color = guide_legend(ncol = 1)) + #make legend one column instead of 3
  scale_color_manual(guide = guide_legend(reverse = FALSE), name = "Tumour Type", labels = c("Liquid", "Solid"), values = c("#392C57", "#4CC5AB")) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRCexplorer2 circRNA PC2&PC3") 
ggarrange(p1, p2, p3, p4,
          ncol = 2, nrow = 2,
          common.legend = TRUE,
          legend = "bottom")
dev.off()
######################################
### Liquid vs Solid Tumours Counts ###
######################################

# get total circRNA counts per cell line per data set
gexpr_df$count <- rowSums(gexpr_df)
isoform_df$count <- rowSums(isoform_df)
ciri_df$count <- rowSums(ciri_df)
circ_df$count <- rowSums(circ_df)

# add cell line labels to original data frames
gexpr_df$cell_line <- cell_line_gexpr
isoform_df$cell_line <- cell_line_isoforms
ciri_df$cell_line <- cell_line_ciri
circ_df$cell_line <- cell_line_circ

# add tumour type label
gexpr_df$tumour.type <- cell_lines[match(gexpr_df$cell_line, cell_lines$X),]$tumour.type
isoform_df$tumour.type <- cell_lines[match(isoform_df$cell_line, cell_lines$X),]$tumour.type
ciri_df$tumour.type <- cell_lines[match(ciri_df$cell_line, cell_lines$X),]$tumour.type
circ_df$tumour.type <- cell_lines[match(circ_df$cell_line, cell_lines$X),]$tumour.type

# add dataset label
gexpr_df$dataset <- c(rep("gCSI", 48), rep("CCLE", 48), rep("GDSC", 48))
isoform_df$dataset <- c(rep("gCSI", 48), rep("CCLE", 48), rep("GDSC", 48))
ciri_df$dataset <- c(rep("gCSI", 48), rep("CCLE", 48), rep("GDSC", 48))
circ_df$dataset <- c(rep("gCSI", 48), rep("CCLE", 48), rep("GDSC", 48))


# get mean
# code from http://www.sthda.com/english/wiki/ggplot2-error-bars-quick-start-guide-r-software-and-data-visualization
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

genes <- data_summary(gexpr_df, varname="count", 
                    groupnames=c("dataset", "tumour.type"))
isoforms <- data_summary(isoform_df, varname="count", 
                      groupnames=c("dataset", "tumour.type"))
ciri <- data_summary(ciri_df, varname="count", 
                      groupnames=c("dataset", "tumour.type"))
circ <- data_summary(circ_df, varname="count", 
                      groupnames=c("dataset", "tumour.type"))

# graph comparisons
p1 <- ggplot(data = genes, aes(x = tumour.type, y = count, fill = dataset)) + 
  geom_bar(stat = "identity", size = 0.5, width = 0.5, position = position_dodge(), color = "black") + 
  scale_fill_manual(values=c("#7F557D", "#D84727", "#6883BA")) + 
  geom_errorbar(aes(ymin=count-sd, ymax=count+sd), width=.2,
                position=position_dodge(0.5)) +
  theme_classic() + 
  scale_y_continuous(limits = c(0, 200000), expand=c(0,0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),    # keep panel borders (removed from theme_classic())
        text = element_text(size = 17), 
        legend.key.size = unit(0.5, 'cm'),
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(x = "Tumour Type", y = "No. of circRNAs", title = "Gene Expression")
p2 <- ggplot(data = isoforms, aes(x = tumour.type, y = count, fill = dataset)) + 
  geom_bar(stat = "identity", size = 0.5, width = 0.5, position = position_dodge(), color = "black") + 
  scale_fill_manual(values=c("#7F557D", "#D84727", "#6883BA")) + 
  geom_errorbar(aes(ymin=count-sd, ymax=count+sd), width=.2,
                position=position_dodge(0.5)) +
  theme_classic() + 
  scale_y_continuous(limits = c(0, 4.5e+05), expand=c(0,0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),    # keep panel borders (removed from theme_classic())
        text = element_text(size = 17), 
        legend.key.size = unit(0.5, 'cm'),
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(x = "Tumour Type", y = "No. of circRNAs", title = "Isoform Expression")
p3 <- ggplot(data = ciri, aes(x = tumour.type, y = count, fill = dataset)) + 
  geom_bar(stat = "identity", size = 0.5, width = 0.5, position = position_dodge(), color = "black") + 
  scale_fill_manual(values=c("#7F557D", "#D84727", "#6883BA")) + 
  geom_errorbar(aes(ymin=count-sd, ymax=count+sd), width=.2,
                position=position_dodge(0.5)) +
  theme_classic() + 
  scale_y_continuous(limits = c(-5, 60), expand=c(0,0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),    # keep panel borders (removed from theme_classic())
        text = element_text(size = 17), 
        legend.key.size = unit(0.5, 'cm'),
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(x = "Tumour Type", y = "No. of circRNAs", title = "CIRI2 circRNA Expression")
p4 <- ggplot(data = circ, aes(x = tumour.type, y = count, fill = dataset)) + 
  geom_bar(stat = "identity", size = 0.5, width = 0.5, position = position_dodge(), color = "black") + 
  scale_fill_manual(values=c("#7F557D", "#D84727", "#6883BA")) + 
  geom_errorbar(aes(ymin=count-sd, ymax=count+sd), width=.2,
                position=position_dodge(0.5)) +
  theme_classic() + 
  scale_y_continuous(limits = c(0, 60), expand=c(0,0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),    # keep panel borders (removed from theme_classic())
        text = element_text(size = 17), 
        legend.key.size = unit(0.5, 'cm'),
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(x = "Tumour Type", y = "No. of circRNAs", title = "CIRCexplorer2 circRNA Expression")
png("../results/compare_solid_liquid_counts.png", width=300, height=200, units='mm', res = 600, pointsize=80)
ggarrange(p1, p2, p3, p4,
          ncol = 2, nrow = 2,
          common.legend = TRUE,
          legend = "right")
dev.off()


# tables of distributions #

gene_distribution <- data.frame(table(gexpr_df$tumour.type), c(sum(gexpr_df[which(gexpr_df$tumour.type == "liquid"),]$count), sum(gexpr_df[which(gexpr_df$tumour.type == "solid"),]$count)))
colnames(gene_distribution) <- c("tumour.type", "sample.count", "gene.count")
gene_distribution$avg.count <- gene_distribution$gene.count / gene_distribution$sample.count

isoform_distribution <- data.frame(table(isoform_df$tumour.type), c(sum(isoform_df[which(isoform_df$tumour.type == "liquid"),]$count), sum(isoform_df[which(isoform_df$tumour.type == "solid"),]$count)))
colnames(isoform_distribution) <- c("tumour.type", "sample.count", "isoform.count")
isoform_distribution$avg.count <- isoform_distribution$isoform.count / isoform_distribution$sample.count

ciri_distribution <- data.frame(table(ciri_df$tumour.type), c(sum(ciri_df[which(ciri_df$tumour.type == "liquid"),]$count), sum(ciri_df[which(ciri_df$tumour.type == "solid"),]$count)))
colnames(ciri_distribution) <- c("tumour.type", "sample.count", "circRNA.count")
ciri_distribution$avg.count <- ciri_distribution$circRNA.count / ciri_distribution$sample.count

circ_distribution <- data.frame(table(circ_df$tumour.type), c(sum(circ_df[which(circ_df$tumour.type == "liquid"),]$count), sum(circ_df[which(circ_df$tumour.type == "solid"),]$count)))
colnames(circ_distribution) <- c("tumour.type", "sample.count", "circRNA.count")
circ_distribution$avg.count <- circ_distribution$circRNA.count / circ_distribution$sample.count

gene_distribution
isoform_distribution
ciri_distribution
circ_distribution

#############################
### circRNA Batch Effects ###
#############################

load("../results/Kallisto_isoforms_counts.RData")

gCSI <- data.frame(gCSI_rna@colData@listData)
CCLE <- data.frame(CCLE_rna@colData@listData)
GDSC <- data.frame(GDSC_rna@colData@listData)


# match ENA and EGA experiments with cell lines
ENA_EXPERIMENT <- GDSC[match(GDSC$cellid, ciri_scores[1:48,]$cell_line),]$Comment.ENA_EXPERIMENT
EGA_EXPERIMENT <- GDSC[match(GDSC$cellid, ciri_scores[1:48,]$cell_line),]$Comment.EGA_EXPERIMENT

ciri_scores$ENA_EXPERIMENT <- rep(ENA_EXPERIMENT, 3)
ciri_scores$EGA_EXPERIMENT <- rep(EGA_EXPERIMENT, 3)

ENA_EXPERIMENT <- GDSC[match(GDSC$cellid, circ_scores[1:48,]$cell_line),]$Comment.ENA_EXPERIMENT
EGA_EXPERIMENT <- GDSC[match(GDSC$cellid, circ_scores[1:48,]$cell_line),]$Comment.EGA_EXPERIMENT

circ_scores$ENA_EXPERIMENT <- rep(ENA_EXPERIMENT, 3)
circ_scores$EGA_EXPERIMENT <- rep(EGA_EXPERIMENT, 3)


png("../results/pca_ciri_EGA_EXPERIMENT.png", width=300, height=200, units='mm', res = 600, pointsize=80)
p1 <- ggplot(data = ciri_scores, aes(x = PC1, y = PC2, group = cell_line)) + 
  geom_point(aes(color = EGA_EXPERIMENT), size = 3) + 
  #scale_color_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), labels=c("CCLE", "gCSI", "GDSC2"), values = c("#392C57", "#4CC5AB", "#3670A0")) + 
  guides(shape = guide_legend(ncol = 1)) + #make legend one column instead of 3
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRI2 PC1&PC2") 
p2 <- ggplot(data = ciri_scores, aes(x = PC2, y = PC3, group = cell_line)) + 
  geom_point(aes(color = EGA_EXPERIMENT), size = 3) + 
  #scale_color_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), labels=c("CCLE", "gCSI", "GDSC2"), values = c("#392C57", "#4CC5AB", "#3670A0")) + 
  guides(shape = guide_legend(ncol = 1)) + #make legend one column instead of 3
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRI2 PC2&PC3") 
p3 <- ggplot(data = ciri_scores, aes(x = PC1, y = PC2, group = cell_line)) + 
  geom_line(show.legend = F) + 
  geom_point(aes(color = EGA_EXPERIMENT), size = 3) + 
  #scale_color_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), labels=c("CCLE", "gCSI", "GDSC2"), values = c("#392C57", "#4CC5AB", "#3670A0")) + 
  guides(shape = guide_legend(ncol = 1)) + #make legend one column instead of 3
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRI2 PC1&PC2") 
p4 <- ggplot(data = ciri_scores, aes(x = PC2, y = PC3, group = cell_line)) + 
  geom_line(show.legend = F) + 
  geom_point(aes(color = EGA_EXPERIMENT), size = 3) + 
  #scale_color_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), labels=c("CCLE", "gCSI", "GDSC2"), values = c("#392C57", "#4CC5AB", "#3670A0")) + 
  guides(shape = guide_legend(ncol = 1)) + #make legend one column instead of 3
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRI2 PC2&PC3") 
ggarrange(p1, p2, p3, p4,
          ncol = 2, nrow = 2,
          common.legend = TRUE,
          legend = "bottom")
dev.off()

png("../results/pca_circ_EGA_EXPERIMENT.png", width=300, height=200, units='mm', res = 600, pointsize=80)
p1 <- ggplot(data = circ_scores, aes(x = PC1, y = PC2, group = cell_line)) + 
  geom_point(aes(color = EGA_EXPERIMENT), size = 3) + 
  #scale_color_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), labels=c("CCLE", "gCSI", "GDSC2"), values = c("#392C57", "#4CC5AB", "#3670A0")) + 
  guides(shape = guide_legend(ncol = 1)) + #make legend one column instead of 3
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRCexplorer2 PC1&PC2") 
p2 <- ggplot(data = circ_scores, aes(x = PC2, y = PC3, group = cell_line)) + 
  geom_point(aes(color = EGA_EXPERIMENT), size = 3) + 
  #scale_color_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), labels=c("CCLE", "gCSI", "GDSC2"), values = c("#392C57", "#4CC5AB", "#3670A0")) + 
  guides(shape = guide_legend(ncol = 1)) + #make legend one column instead of 3
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRCexplorer2 PC2&PC3") 
p3 <- ggplot(data = circ_scores, aes(x = PC1, y = PC2, group = cell_line)) + 
  geom_line(show.legend = F) + 
  geom_point(aes(color = EGA_EXPERIMENT), size = 3) + 
  #scale_color_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), labels=c("CCLE", "gCSI", "GDSC2"), values = c("#392C57", "#4CC5AB", "#3670A0")) + 
  guides(shape = guide_legend(ncol = 1)) + #make legend one column instead of 3
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRCexplorer2 PC1&PC2") 
p4 <- ggplot(data = circ_scores, aes(x = PC2, y = PC3, group = cell_line)) + 
  geom_line(show.legend = F) + 
  geom_point(aes(color = EGA_EXPERIMENT), size = 3) + 
  #scale_color_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), labels=c("CCLE", "gCSI", "GDSC2"), values = c("#392C57", "#4CC5AB", "#3670A0")) + 
  guides(shape = guide_legend(ncol = 1)) + #make legend one column instead of 3
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRCexplorer2 PC2&PC3") 
ggarrange(p1, p2, p3, p4,
          ncol = 2, nrow = 2,
          common.legend = TRUE,
          legend = "bottom")
dev.off()