suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))
suppressMessages(library(matrixStats))

###########################################################
########## Get circRNA CIRI transcript stability ##########
###########################################################

load("../results/circRNAcounts.RData")

#function to compute spearman correlation for each dataset pair (gCSI/CCLE; gCSI/GDSC; GDSC/CCLE)                                 
compute_CIRIspearman  <- function(x){
  i <- x
  gcsi_ccle_spearman <- suppressWarnings(cor(x=as.numeric(ciri_gcsi[i,]), y=as.numeric(ciri_ccle[i,]), method="spearman")) #gCSI vs CCLE
  gcsi_gdsc_spearman <- suppressWarnings(cor(x=as.numeric(ciri_gcsi[i,]), y=as.numeric(ciri_gdsc[i,]), method="spearman")) #gCSI vs GDSC
  gdsc_ccle_spearman <- suppressWarnings(cor(x=as.numeric(ciri_gdsc[i,]), y=as.numeric(ciri_ccle[i,]), method="spearman")) #GDSC vs CCLE
  
  combined <- c(gcsi_ccle_spearman, gcsi_gdsc_spearman, gdsc_ccle_spearman)
  return(combined)
}                                         

##compute spearman
transcripts <- rownames(ciri_gcsi)                                     

spearman_compute <- sapply(transcripts, compute_CIRIspearman)
circRNA_CIRI_stability <- as.data.frame(t(as.data.frame(spearman_compute)))
colnames(circRNA_CIRI_stability) <- c("gCSI_CCLE", "gCSI_GDSC", "GDSC_CCLE")  



###########################################################
########## Get circRNA CIRC transcript stability ##########
###########################################################

#function to compute spearman correlation for each dataset pair (gCSI/CCLE; gCSI/GDSC; GDSC/CCLE)                                 
compute_CIRCspearman  <- function(x){
  i <- x
  gcsi_ccle_spearman <- suppressWarnings(cor(x=as.numeric(circ_gcsi[i,]), y=as.numeric(circ_ccle[i,]), method="spearman")) #gCSI vs CCLE
  gcsi_gdsc_spearman <- suppressWarnings(cor(x=as.numeric(circ_gcsi[i,]), y=as.numeric(circ_gdsc[i,]), method="spearman")) #gCSI vs GDSC
  gdsc_ccle_spearman <- suppressWarnings(cor(x=as.numeric(circ_gdsc[i,]), y=as.numeric(circ_ccle[i,]), method="spearman")) #GDSC vs CCLE
  
  combined <- c(gcsi_ccle_spearman, gcsi_gdsc_spearman, gdsc_ccle_spearman)
  return(combined)
}                                         

##compute spearman
transcripts <- rownames(circ_gcsi)                                     

spearman_compute <- sapply(transcripts, compute_CIRCspearman)
circRNA_CIRC_stability <- as.data.frame(t(as.data.frame(spearman_compute)))
colnames(circRNA_CIRC_stability) <- c("gCSI_CCLE", "gCSI_GDSC", "GDSC_CCLE")


###########################
#### Isoform stability ####
###########################

load("../results/isoform_expression.RData")

#function to compute spearman correlation for each dataset pair (gCSI/CCLE; gCSI/GDSC; GDSC/CCLE)                                 
compute_spearman  <- function(x){
  i <- x
  gcsi_ccle_spearman <- suppressWarnings(cor(x=as.numeric(expr_gcsi_i[i,]), y=as.numeric(expr_ccle_i[i,]), method="spearman"))
  gcsi_gdsc_spearman <- suppressWarnings(cor(x=as.numeric(expr_gcsi_i[i,]), y=as.numeric(expr_gdsc_i[i,]), method="spearman"))
  gdsc_ccle_spearman <- suppressWarnings(cor(x=as.numeric(expr_gdsc_i[i,]), y=as.numeric(expr_ccle_i[i,]), method="spearman"))
  
  combined <- c(gcsi_ccle_spearman, gcsi_gdsc_spearman, gdsc_ccle_spearman)
  return(combined)
}                                         

##compute pearson & spearman
transcripts <- rownames(expr_gcsi_i)                                     

spearman_compute <- sapply(transcripts, compute_spearman)
transcript_stability <- as.data.frame(t(as.data.frame(spearman_compute)))
colnames(transcript_stability) <- c("gcsi_ccle_spearman", "gcsi_gdsc_spearman", "gdsc_ccle_spearman")  


##compute median expression for each transcript across biological replicates for each dataset
gcsi_matrix <- as.matrix(expr_gcsi_i)                                 
ccle_matrix <- as.matrix(expr_ccle_i)  
gdsc_matrix <- as.matrix(expr_gdsc_i)

gcsi_median <- as.numeric(rowMedians(gcsi_matrix))                                 
ccle_median <- as.numeric(rowMedians(ccle_matrix))                                  
gdsc_median <- as.numeric(rowMedians(gdsc_matrix))                                      


##compile everything into one data frame                                 
transcript_stability$transcript_id <- transcripts
rownames(transcript_stability) <- transcript_stability$transcript_id
transcript_stability$gcsi_median <- gcsi_median
transcript_stability$ccle_median <- ccle_median
transcript_stability$gdsc_median <- gdsc_median


#remove transcripts that have median expression of 0 across all datasets
transcript_stability <- transcript_stability[-which(transcript_stability$gcsi_median == 0 & transcript_stability$ccle_median == 0 & transcript_stability$gdsc_median == 0),]  

save(transcript_stability, file = "../results/transcript_stability.RData")



########################
#### Gene Stability ####
########################

load("../results/gene_expression.RData")

#function to compute spearman correlation for each dataset pair (gCSI/CCLE; gCSI/GDSC; GDSC/CCLE)                                 
compute_spearman  <- function(x){
  i <- x
  gcsi_ccle_spearman <- suppressWarnings(cor(x=as.numeric(expr_gcsi_p[i,]), y=as.numeric(expr_ccle_p[i,]), method="spearman"))
  gcsi_gdsc_spearman <- suppressWarnings(cor(x=as.numeric(expr_gcsi_p[i,]), y=as.numeric(expr_gdsc_p[i,]), method="spearman"))
  gdsc_ccle_spearman <- suppressWarnings(cor(x=as.numeric(expr_gdsc_p[i,]), y=as.numeric(expr_ccle_p[i,]), method="spearman"))
  
  combined <- c(gcsi_ccle_spearman, gcsi_gdsc_spearman, gdsc_ccle_spearman)
  return(combined)
}                                         

##compute pearson & spearman
genes <- rownames(expr_gcsi_p)                                     


spearman_compute <- sapply(genes, compute_spearman)
gene_stability <- as.data.frame(t(as.data.frame(spearman_compute)))
colnames(gene_stability) <- c("gcsi_ccle_spearman", "gcsi_gdsc_spearman", "gdsc_ccle_spearman")  


##compute median expression for each gene across biological replicates for each dataset
gcsi_matrix <- as.matrix(expr_gcsi_p)                                 
ccle_matrix <- as.matrix(expr_ccle_p)  
gdsc_matrix <- as.matrix(expr_gdsc_p)

gcsi_median <- as.numeric(rowMedians(gcsi_matrix))                                 
ccle_median <- as.numeric(rowMedians(ccle_matrix))                                  
gdsc_median <- as.numeric(rowMedians(gdsc_matrix))                                      


##compile everything into one data frame                                 
gene_stability$gene_id <- genes
rownames(gene_stability) <- gene_stability$gene_id
gene_stability$gcsi_median <- gcsi_median
gene_stability$ccle_median <- ccle_median
gene_stability$gdsc_median <- gdsc_median


#######################
### Plot Stability ####
#######################

# format circRNA stability data
colnames(circRNA_CIRI_stability) <- c("gCSI/CCLE", "gCSI/GDSC", "GDSC/CCLE")
circRNA_CIRI_stability <- melt(circRNA_CIRI_stability)
colnames(circRNA_CIRI_stability) <- c("dataset", "stability")


colnames(circRNA_CIRC_stability) <- c("gCSI/CCLE", "gCSI/GDSC", "GDSC/CCLE")
circRNA_CIRC_stability <- melt(circRNA_CIRC_stability)
colnames(circRNA_CIRC_stability) <- c("dataset", "stability")


# format mRNA transcript stability data
transcript_stability <- transcript_stability[,c("gcsi_ccle_spearman", "gcsi_gdsc_spearman", "gdsc_ccle_spearman")]
colnames(transcript_stability) <- c("gCSI/CCLE", "gCSI/GDSC", "GDSC/CCLE")
transcript_stability <- melt(transcript_stability)
colnames(transcript_stability) <- c("dataset", "stability")

# format gene transcript stability data
gene_stability <- gene_stability[,c("gcsi_ccle_spearman", "gcsi_gdsc_spearman", "gdsc_ccle_spearman")]
colnames(gene_stability) <- c("gCSI/CCLE", "gCSI/GDSC", "GDSC/CCLE")
gene_stability <- melt(gene_stability)
colnames(gene_stability) <- c("dataset", "stability")
gene_stability[is.na(gene_stability)] <- 0

# add stable / unstable labels 
label_stable <- function(df) {
  label <- c()
  for (i in 1:length(rownames(df))) {
    if (df[i,2] > 0.8) {
      label <- c(label, "stable")
    }
    else {
      label <- c(label, "unstable")
    }
  }
  df$label <- label
  return(df)
}

circRNA_CIRI_stability <- label_stable(circRNA_CIRI_stability)
circRNA_CIRC_stability <- label_stable(circRNA_CIRC_stability)
transcript_stability <- label_stable(transcript_stability)
gene_stability <- label_stable(gene_stability)


# plot histograms of SI distribution 
p1 <- ggplot(gene_stability, aes(x = stability, fill = label, color = label)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 50) +
  geom_vline(xintercept = 0.8, color = "#982649", size = 1, linetype = "dotted") +
  facet_grid(dataset ~ .) + scale_x_continuous(breaks = seq(-0.2, 1.0 , 0.2)) +
  scale_color_manual(values=c("black", "black")) +
  scale_fill_manual(values=c("#982649", "white")) +
  theme_classic() +
  theme(text = element_text(size = 15), 
        legend.key.size = unit(0.4, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 17), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15),
        legend.title=element_blank()) +
  labs(title = "Gene Expression", x = "Stability Index \n", y = "SI Count")

p2 <- ggplot(transcript_stability, aes(x = stability, fill = label, color = label)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 50) +
  geom_vline(xintercept = 0.8, color = "#982649", size = 1, linetype = "dotted") +
  facet_grid(dataset ~ .) + scale_x_continuous(breaks = seq(-0.2, 1.0 , 0.2)) +
  scale_color_manual(values=c("black", "black")) +
  scale_fill_manual(values=c("#982649", "white")) +
  theme_classic() +
  theme(text = element_text(size = 15), 
        plot.title = element_text(hjust = 0.5, size = 17), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15),
        legend.title=element_blank()) +
  labs(title = "Isoform Expression", x = "Stability Index \n", y = "SI Count")

p3 <- ggplot(circRNA_CIRI_stability, aes(x = stability, fill = label, color = label)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 50) +
  geom_vline(xintercept = 0.8, color = "#982649", size = 1, linetype = "dotted") +
  facet_grid(dataset ~ .) + scale_x_continuous(breaks = seq(-0.2, 1.0 , 0.2)) +
  ylim(0, 200) +
  scale_color_manual(values=c("black", "black")) +
  scale_fill_manual(values=c("#982649", "white")) +
  theme_classic() +
  theme(text = element_text(size = 15), 
        plot.title = element_text(hjust = 0.5, size = 17), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15),
        legend.title=element_blank()) +
  labs(title = "CIRI2 circRNA Expression", x = "Stability Index", y = "SI Count")

p4 <- ggplot(circRNA_CIRC_stability, aes(x = stability, fill = label, color = label)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 50) +
  geom_vline(xintercept = 0.8, color = "#982649", size = 1, linetype = "dotted") +
  facet_grid(dataset ~ .) + scale_x_continuous(breaks = seq(-0.2, 1.0 , 0.2)) +
  ylim(0, 200) +
  scale_color_manual(values=c("black", "black")) +
  scale_fill_manual(values=c("#982649", "white")) +
  theme_classic() +
  theme(text = element_text(size = 15), 
        plot.title = element_text(hjust = 0.5, size = 17), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15),
        legend.title=element_blank()) +
  labs(title = "CIRCexplorer2 circRNA Expression", x = "Stability Index", y = "SI Count")


png("../results/supplementary_figure4.png", width=300, height=275, units='mm', res = 600, pointsize=80)
ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2,
          common.legend=TRUE, legend = "right")
dev.off()


