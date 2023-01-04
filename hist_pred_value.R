suppressMessages(library(PharmacoGx))
suppressMessages(library(survcomp))
suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))

##### All Genes #####

load("../results/inter.RData")

#read in PSets isoforms.counts
load("../results/Kallisto_isoforms_counts.RData")
#read in PSets
gCSI <- readRDS("/data/PSets/gCSI.rds")



drugs <- c("Crizotinib","Erlotinib","Irinotecan","Lapatinib","Paclitaxel","PD-0325901")
transcripts <- read.csv("../results/transcript_stability.csv")
transcripts <- transcripts$transcript_id


##compute mean stability for all transcripts (genome_wide)

#gcsi/ccle
gcsi_ccle_stability <- as.data.frame(sapply(transcripts, function(i) cor(as.numeric(gCSI_rna@assays@data$expr[i,]),as.numeric(CCLE_rna@assays@data$expr[i,]), method="spearman")))

#gcsi/gdsc
gcsi_gdsc_stability <- as.data.frame(sapply(transcripts, function(i) cor(as.numeric(gCSI_rna@assays@data$expr[i,]),as.numeric(GDSC_rna@assays@data$expr[i,]), method="spearman")))

#gdsc/ccle
gdsc_ccle_stability <- as.data.frame(sapply(transcripts, function(i) cor(as.numeric(GDSC_rna@assays@data$expr[i,]),as.numeric(CCLE_rna@assays@data$expr[i,]), method="spearman")))

mean_stability <- as.data.frame(rowMeans(cbind(gcsi_ccle_stability, gcsi_gdsc_stability, gdsc_ccle_stability)))
colnames(mean_stability) <- "mean_stability"

load("/data/biomarker/meta_pearson_transcript.RData")


#Add gene ids to df's
genes <- gCSI@molecularProfiles$Kallisto_0.46.1.isoforms.counts@elementMetadata@listData$gene_name[match(rownames(meta_pearson), gCSI@molecularProfiles$Kallisto_0.46.1.isoforms.counts@elementMetadata@listData$transcript_id)]

meta_pearson$gene_id <- genes
mean_stability$gene_id <- genes 


#create final df for spearman of mean stability + pearson for all transcripts per gene
spearman_df <- data.frame(matrix(data=NA, nrow = length(unique(genes)), ncol=length(drugs)))
colnames(spearman_df) <- drugs
rownames(spearman_df) <- unique(genes)


spearman_df_randomize <- data.frame(matrix(data=NA, nrow = length(unique(genes)), ncol=length(drugs)))
colnames(spearman_df_randomize) <- drugs
rownames(spearman_df_randomize) <- unique(genes)


for (g in rownames(spearman_df_randomize)){
  for(d in drugs){
    
    if (!length(which(mean_stability$gene_id == g)) < 5){
      s <- cor(sample(mean_stability$mean_stability[which(mean_stability$gene_id == g)]), meta_pearson[which(meta_pearson$gene_id == g),d], method = "spearman")
      spearman_df_randomize[g,d] <- s
    }
  }
}

for (g in rownames(spearman_df)){
  for(d in drugs){
    
    if (!length(which(mean_stability$gene_id == g)) < 5){
      s <- cor(mean_stability$mean_stability[which(mean_stability$gene_id == g)], meta_pearson[which(meta_pearson$gene_id == g),d], method = "spearman")
      spearman_df[g,d] <- s
    }
  }
}

spearman_df <- spearman_df[which(!is.na(spearman_df$Crizotinib)),]
spearman_df_randomize <- spearman_df_randomize[which(!is.na(spearman_df_randomize$Crizotinib)),]


# plot histograms
spearman_df$gene <- rownames(spearman_df)
spearman_df <- melt(spearman_df)
spearman_df$label <- c(rep("Non-randomized SI", dim(spearman_df)[1]))

spearman_df_randomize$gene <- rownames(spearman_df_randomize)
spearman_df_randomize <- melt(spearman_df_randomize)
spearman_df_randomize$label <- c(rep("Randomized SI", dim(spearman_df_randomize)[1]))

hist_pred_val <- rbind(spearman_df, spearman_df_randomize)

png("../results/supplementary_figure11.png", width = 9, height = 9, res = 600, units = "in")
ggplot() + 
  geom_histogram(data = hist_pred_val[which(hist_pred_val$label == "Non-randomized SI"),], aes(x = value, fill = "#71A2B6"), color = "black") + 
  geom_histogram(data = hist_pred_val[which(hist_pred_val$label == "Randomized SI"),], aes(x = value, fill = "#CA5A7D"), color = "black", alpha = 0.45) +
  facet_wrap(~variable, ncol = 3, nrow = 2) +
  scale_y_continuous(limits = c(0, 200), expand=c(0,0)) +
  scale_fill_manual(values=c("#71A2B6", "#CA5A7D"), labels = c("Non-randomized SI", "Randomized SI")) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(1, 'cm'),
        legend.position = "bottom",
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5, size = 17), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  guides(fill = guide_legend(ncol = 2)) +
  labs(x = "Spearman Correlation Coefficient", y = "Count")
dev.off()

