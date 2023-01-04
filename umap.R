suppressMessages(library(plyr))
suppressMessages(library(tidyverse))
suppressMessages(library(umap))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(viridis))


###########################################
######### Load in Expression Data #########
###########################################

######### Load and Isoform and Gene Expression Data #########

# load in expression data
load("../results/gene_expression.RData")
load("../results/isoform_expression.RData")


# save labels for cell lines
cell_line_gexpr <- c(colnames(expr_gcsi_p), colnames(expr_ccle_p), colnames(expr_gdsc_p))
cell_line_isoforms <- c(colnames(expr_gcsi_i), colnames(expr_ccle_i), colnames(expr_gdsc_i))


# format gene expression data
expr_gcsi_p <- as.data.frame(t(expr_gcsi_p))
expr_ccle_p <- as.data.frame(t(expr_ccle_p))
expr_gdsc_p <- as.data.frame(t(expr_gdsc_p))

rownames(expr_gcsi_p) <- NULL
rownames(expr_ccle_p) <- NULL
rownames(expr_gdsc_p) <- NULL


# format isoform expression data
expr_gcsi_i <- as.data.frame(t(expr_gcsi_i))
expr_ccle_i <- as.data.frame(t(expr_ccle_i))
expr_gdsc_i <- as.data.frame(t(expr_gdsc_i))

rownames(expr_gcsi_i) <- NULL
rownames(expr_ccle_i) <- NULL
rownames(expr_gdsc_i) <- NULL


# merge all expression values from each data set 
gexpr_df <- rbind.fill(expr_gcsi_p, expr_ccle_p, expr_gdsc_p)
isoform_df <- rbind.fill(expr_gcsi_i, expr_ccle_i, expr_gdsc_i)


######### Load and Filter circRNA Expression Data #########

load("../results/circRNAcounts.RData")

### filter circRNA transcripts with low detection rates
ciri_gcsi_filtered <- ciri_gcsi[-which(rownames(ciri_gcsi) %in% names(which(rowSums(ciri_gcsi == 0) > 36))),]
ciri_ccle_filtered <- ciri_ccle[-which(rownames(ciri_ccle) %in% names(which(rowSums(ciri_gcsi == 0) > 36))),]
ciri_gdsc_filtered <- ciri_gdsc[-which(rownames(ciri_gdsc) %in% names(which(rowSums(ciri_gdsc == 0) > 36))),]
circ_gcsi_filtered <- circ_gcsi[-which(rownames(circ_gcsi) %in% names(which(rowSums(circ_gcsi == 0) > 36))),]
circ_ccle_filtered <- circ_ccle[-which(rownames(circ_ccle) %in% names(which(rowSums(circ_ccle == 0) > 36))),]
circ_gdsc_filtered <- circ_gdsc[-which(rownames(circ_gdsc) %in% names(which(rowSums(circ_gdsc == 0) > 36))),]

# save labels for cell lines
cell_line_ciri <- c(colnames(ciri_gcsi_filtered), colnames(ciri_ccle_filtered), colnames(ciri_gdsc_filtered))
cell_line_circ <- c(colnames(circ_gcsi_filtered), colnames(circ_ccle_filtered), colnames(circ_gdsc_filtered))

## get common circRNA transcripts
transcripts <- data.frame(c(rownames(ciri_gcsi_filtered), rownames(ciri_ccle_filtered), rownames(ciri_gdsc_filtered),
                            rownames(circ_gcsi_filtered), rownames(circ_ccle_filtered), rownames(circ_gdsc_filtered)))
colnames(transcripts) <- "circRNAID"
transcript_counts <- transcripts %>% count(circRNAID)

# remove circRNA transcripts that are only in one method_dataset object
transcript_counts <- transcript_counts[-which(transcript_counts$n == 1),] 
common_transcripts <- transcript_counts$circRNAID

## format circRNA data
ciri_gcsi_filtered <- as.data.frame(t(ciri_gcsi_filtered[which(rownames(ciri_gcsi_filtered) %in% common_transcripts),]))
ciri_ccle_filtered <- as.data.frame(t(ciri_ccle_filtered[which(rownames(ciri_ccle_filtered) %in% common_transcripts),]))
ciri_gdsc_filtered <- as.data.frame(t(ciri_gdsc_filtered[which(rownames(ciri_gdsc_filtered) %in% common_transcripts),]))
circ_gcsi_filtered <- as.data.frame(t(circ_gcsi_filtered[which(rownames(circ_gcsi_filtered) %in% common_transcripts),]))
circ_ccle_filtered <- as.data.frame(t(circ_ccle_filtered[which(rownames(circ_ccle_filtered) %in% common_transcripts),]))
circ_gdsc_filtered <- as.data.frame(t(circ_gdsc_filtered[which(rownames(circ_gdsc_filtered) %in% common_transcripts),]))

rownames(ciri_gcsi_filtered) <- NULL
rownames(ciri_ccle_filtered) <- NULL
rownames(ciri_gdsc_filtered) <- NULL
rownames(circ_gcsi_filtered) <- NULL
rownames(circ_ccle_filtered) <- NULL
rownames(circ_gdsc_filtered) <- NULL


ciri_df <- rbind.fill(ciri_gcsi_filtered, ciri_ccle_filtered, ciri_gdsc_filtered)
ciri_df[is.na(ciri_df)] <- 0
ciri_df[] <- lapply(ciri_df, as.double)


circ_df <- rbind.fill(circ_gcsi_filtered, circ_ccle_filtered, circ_gdsc_filtered)
circ_df[is.na(circ_df)] <- 0
circ_df[] <- lapply(circ_df, as.double)



######### Load in lung circRNA expression data #########
load("../results/lungcounts.RData")

rownames(lung_ribozero_ciri) <- NULL
rownames(lung_polyA_ciri) <- NULL
rownames(lung_ribozero_circ) <- NULL
rownames(lung_polyA_circ) <- NULL

lung_ribozero_ciri <- as.data.frame(lung_ribozero_ciri)
lung_polyA_ciri <- as.data.frame(lung_polyA_ciri)
lung_ribozero_circ <- as.data.frame(lung_ribozero_circ)
lung_polyA_circ <- as.data.frame(lung_polyA_circ)

lung_df <- rbind.fill(lung_ribozero_ciri, lung_polyA_ciri, lung_ribozero_circ, lung_polyA_circ)
lung_df[is.na(lung_df)] <- 0
lung_df[] <- lapply(lung_df, as.double)


save(gexpr_df, isoform_df, ciri_df, circ_df, lung_df,
     cell_line_gexpr, cell_line_isoforms, cell_line_ciri, cell_line_circ,
     file="../results/umapdf.RData")

###########################################
######### Create UMAP Projections #########
###########################################

# load in solid tumours
load("../results/solid_tumours.RData")
cell_line_gexpr <- gexpr_df$cell_line
cell_line_isoforms <- isoform_df$cell_line
cell_line_ciri <- ciri_df$cell_line
cell_line_circ <- circ_df$cell_line

gexpr_df <- gexpr_df[,-which(names(gexpr_df) %in% c("cell_line", "tumour.type"))]
isoform_df <- isoform_df[,-which(names(isoform_df) %in% c("cell_line", "tumour.type"))]
ciri_df <- ciri_df[,-which(names(ciri_df) %in% c("cell_line", "tumour.type"))]
circ_df <- circ_df[,-which(names(circ_df) %in% c("cell_line", "tumour.type"))]

# create gene expression umap projection
gexpr_umap <- umap(gexpr_df)

gexpr_umap <- as.data.frame(gexpr_umap$layout)
colnames(gexpr_umap) <- c("UMAP1", "UMAP2")
gexpr_umap$cell_line <- cell_line_gexpr
gexpr_umap$dataset <- c(rep("gCSI", 48), rep("CCLE", 48), rep("GDSC", 48)) # change 48 to 39 if doing solid tumours


# create isoforms umap projection
isoform_umap <- umap(isoform_df)

isoform_umap <- as.data.frame(isoform_umap$layout)
colnames(isoform_umap) <- c("UMAP1", "UMAP2")
isoform_umap$cell_line <- cell_line_isoforms
isoform_umap$dataset <- c(rep("gCSI", 48), rep("CCLE", 48), rep("GDSC", 48))



# create circRNA umap projection
ciri_umap <- umap(ciri_df)

ciri_umap <- as.data.frame(ciri_umap$layout)
colnames(ciri_umap) <- c("UMAP1", "UMAP2")
ciri_umap$cell_line <- cell_line_ciri
ciri_umap$dataset <- c(rep("gCSI", 48), rep("CCLE", 48), rep("GDSC", 48))



circ_umap <- umap(circ_df)

circ_umap <- as.data.frame(circ_umap$layout)
colnames(circ_umap) <- c("UMAP1", "UMAP2")
circ_umap$cell_line <- cell_line_circ
circ_umap$dataset <- c(rep("gCSI", 48), rep("CCLE", 48), rep("GDSC", 48))


# create lung umap projection
lung_umap <- umap(lung_df)

lung_umap <- as.data.frame(lung_umap$layout)
colnames(lung_umap) <- c("UMAP1", "UMAP2")
lung_umap$tumourID <- c(rep(paste0("tumor",1:51), 4))
lung_umap$Selection <- c(rep("poly(A)", 51), rep("RiboZero", 51), rep("poly(A)", 51), rep("RiboZero", 51))
lung_umap$Method <- c(rep("CIRI2", 102), rep("CIRCexplorer2", 102))


##############################################
######### Visualize UMAP Projections #########
##############################################

######### Lines to Connect Cell Line Data Points #########
p1 <- ggplot(data = gexpr_umap, aes(x = UMAP1, y = UMAP2, group = cell_line)) + 
  geom_line(show.legend = F) + 
  geom_point(aes(color = dataset, shape = dataset), size = 3) + 
  scale_color_manual(guide = guide_legend(reverse = FALSE, title = "Dataset"), 
                     labels=c("CCLE", "gCSI", "GDSC2"), 
                     values = c("#392C57", "#4CC5AB", "#3670A0")) + 
  guides(shape = guide_legend(ncol = 1)) + #make legend one column instead of 3
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(0.7, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "Gene Expression") +
  scale_shape_discrete(name = "Dataset", labels = c("CCLE", "gCSI", "GDSC2"))

p2 <- ggplot(data = isoform_umap, aes(x = UMAP1, y = UMAP2, group = cell_line)) + 
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
  labs(title = "Isoform Expression") 


p3 <- ggplot(data = ciri_umap, aes(x = UMAP1, y = UMAP2, group = cell_line)) + 
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
  labs(title = "CIRI2 circRNA Expression") +
  scale_shape_discrete(name = "Dataset", labels = c("CCLE", "gCSI", "GDSC2"))



p4 <- ggplot(data = circ_umap, aes(x = UMAP1, y = UMAP2, group = cell_line)) + 
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
  labs(title = "CIRCexplorer2 circRNA Expression") +
  scale_shape_discrete(name = "Dataset", labels = c("CCLE", "gCSI", "GDSC2"))


png("../results/figure1.png", width=300, height=225, units='mm', res = 600, pointsize=80)
ggarrange(p1, p2, p3, p4,
          ncol = 2, nrow = 2,
          common.legend = TRUE,
          legend = "right")
dev.off()



######### Facets to Group by Cell Lines #########

png("../results/supplementary_figure3a.png", width=300, height=225, units='mm', res = 600, pointsize=80)
ggplot() + 
  geom_point(data = gexpr_umap[,1:2], aes(x = UMAP1, y = UMAP2), color = "grey", size = 3) +
  geom_point(data = gexpr_umap, aes(x = UMAP1, y = UMAP2, color = cell_line, shape = dataset), size = 4) +
  facet_wrap(~cell_line, ncol = 6, nrow = 8) +
  guides(color = "none") +
  scale_color_viridis(begin = 0.2, end = 0.8, discrete = TRUE, option = "mako") +
  guides(shape = guide_legend(ncol = 1)) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(1, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 17), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "Gene Expression") +
  scale_shape_discrete(name = "Dataset", labels = c("CCLE", "gCSI", "GDSC2"))
dev.off()

png("../results/supplementary_figure3b.png", width=300, height=225, units='mm', res = 600, pointsize=80)
ggplot() + 
  geom_point(data = isoform_umap[,1:2], aes(x = UMAP1, y = UMAP2), color = "grey", size = 3) +
  geom_point(data = isoform_umap, aes(x = UMAP1, y = UMAP2, color = cell_line, shape = dataset), size = 4) +
  facet_wrap(~cell_line, ncol = 6, nrow = 8) +
  guides(color = "none") +
  scale_color_viridis(begin = 0.2, end = 0.8, discrete = TRUE, option = "mako") +
  guides(shape = guide_legend(ncol = 1)) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),    
        text = element_text(size = 15), 
        legend.key.size = unit(1, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 17), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "Isoform Expression") +
  scale_shape_discrete(name = "Dataset", labels = c("CCLE", "gCSI", "GDSC2"))
dev.off()

png("../results/supplementary_figure3c.png", width=300, height=225, units='mm', res = 600, pointsize=80)
ggplot() + 
  geom_point(data = ciri_umap[,1:2], aes(x = UMAP1, y = UMAP2), color = "grey", size = 3) +
  geom_point(data = ciri_umap, aes(x = UMAP1, y = UMAP2, color = cell_line, shape = dataset), size = 4) +
  facet_wrap(~cell_line, ncol = 6, nrow = 8) +
  guides(color = "none") +
  scale_color_viridis(begin = 0.2, end = 0.8, discrete = TRUE, option = "mako") +
  guides(shape = guide_legend(ncol = 1)) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(1, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 17), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRI2 circRNA Expression") +
  scale_shape_discrete(name = "Dataset", labels = c("CCLE", "gCSI", "GDSC2"))
dev.off()

png("../results/supplementary_figure3d.png", width=300, height=225, units='mm', res = 600, pointsize=80)
ggplot() + 
  geom_point(data = circ_umap[,1:2], aes(x = UMAP1, y = UMAP2), color = "grey", size = 3) +
  geom_point(data = circ_umap, aes(x = UMAP1, y = UMAP2, color = cell_line, shape = dataset), size = 4) +
  facet_wrap(~cell_line, ncol = 6, nrow = 8) +
  guides(color = "none") +
  scale_color_viridis(begin = 0.2, end = 0.8, discrete = TRUE, option = "mako") +
  guides(shape = guide_legend(ncol = 1)) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(1, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 17), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRCexplorer2 circRNA Expression") +
  scale_shape_discrete(name = "Dataset", labels = c("CCLE", "gCSI", "GDSC2"))
dev.off()



######### Lung UMAP #########

png("../results/supplementary_figure7.png", width=150, height=125, units='mm', res = 600, pointsize=80)
ggplot(data = lung_umap, aes(x = UMAP1, y = UMAP2)) + 
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
  labs(title = "Lung Samples circRNA Expression") 
dev.off()



