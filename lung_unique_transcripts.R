suppressMessages(library(reshape2))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))


#===========================================================================
##### circRNA IDs as Identifiers #####

# load in circRNA counts from lung data
load("../results/lungcounts.RData")


sum_transcripts <- function(polyA_df, ribozero_df) {
  
  # format data frames and compute counts for each circRNA transcript
  polyA_df <- data.frame(t(polyA_df))
  ribozero_df <- data.frame(t(ribozero_df))
  
  polyA_df$count <- rowSums(polyA_df)
  ribozero_df$count <- rowSums(ribozero_df)
  
  # get common circRNA transcript IDs
  transcripts <- intersect(rownames(polyA_df), rownames(ribozero_df))
  df <- data.frame(transcript_id = transcripts)
  
  polyA <- c()
  ribozero <- c()
  
  # get transcript counts
  for (i in 1:length(df$transcript_id)) {
    transcript <- df$transcript_id[i]
    
    count_polyA <- polyA_df[which(rownames(polyA_df) == transcript),]$count
    polyA <- c(polyA, count_polyA)
    
    count_ribozero <- ribozero_df[which(rownames(ribozero_df) == transcript),]$count
    ribozero <- c(ribozero, count_ribozero)
  }
  
  df$polyA <- polyA
  df$ribozero <- ribozero
  
  return(df)
}

CIRI_counts <- sum_transcripts(lung_polyA_ciri, lung_ribozero_ciri)
CIRC_counts <- sum_transcripts(lung_polyA_circ, lung_ribozero_circ)



### format data.frames for histogram ###
CIRI_counts <- melt(CIRI_counts)
CIRC_counts <- melt(CIRC_counts)

colnames(CIRI_counts) <- c("transcript_id", "Method", "Count")
colnames(CIRC_counts) <- c("transcript_id", "Method", "Count")

# log2 + 1 normalize circRNA counts
CIRI_counts$Count <- log2(CIRI_counts$Count + 1)
CIRC_counts$Count <- log2(CIRC_counts$Count + 1)

# plot bar plots
p1 <- ggplot(data = CIRI_counts, aes(x = reorder(transcript_id, -Count), y = Count, fill = Method, alpha = Method)) + 
  geom_bar(stat = "identity", position = "identity", width = 1) +
  scale_y_continuous(limits = c(0, 9), expand=c(0,0)) +
  scale_fill_manual(values=c("#4CC5AB", "#392C57")) +
  scale_alpha_manual(values=c(0.8, 0.6)) +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),    
        text = element_text(size = 15), 
        legend.key.size = unit(1, 'cm'),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 17)) +
  labs(title = "CIRI2 Transcript Counts", 
       y = "log2 Normalized \nCounts", x = "circRNA Transcript")

p2 <- ggplot(data = CIRC_counts, aes(x = reorder(transcript_id, -Count), y = Count, fill = Method, alpha = Method)) + 
  geom_bar(stat = "identity", position = "identity", width = 1) +
  scale_y_continuous(limits = c(0, 9), expand=c(0,0)) +
  scale_fill_manual(values=c("#4CC5AB", "#392C57")) +
  scale_alpha_manual(values=c(0.8, 0.6)) +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),    
        text = element_text(size = 15), 
        legend.key.size = unit(1, 'cm'),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 17)) +
  labs(title = "CIRCexplorer2 Transcript Counts", 
       y = "log2 Normalized \nCounts", x = "circRNA Transcript")


png("../results/supplementary_figure8b.png", width=275, height=200, units='mm', res = 600, pointsize=80)
ggarrange(p1, p2, ncol = 1, common.legend = TRUE, legend = "bottom")
dev.off()


# plot density plots
p1 <- ggplot(data = CIRI_counts) +
  geom_density(aes(x = Count, fill = Method), color = "black", alpha = 0.5) +
  scale_fill_manual(values=c("#4CC5AB", "#392C57")) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),    
        text = element_text(size = 15), 
        legend.key.size = unit(1, 'cm'),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 17), 
        axis.text.x = element_text(size=12, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRI2 Transcript Counts", 
       y = "Density", x = "log2 Normalized Counts") +
  guides(fill = guide_legend(title = "Method", ncol = 2))

p2 <- ggplot(data = CIRC_counts) +
  geom_density(aes(x = Count, fill = Method), color = "black", alpha = 0.5) +
  scale_fill_manual(values=c("#4CC5AB", "#392C57")) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),    
        text = element_text(size = 15), 
        legend.key.size = unit(1, 'cm'),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 17), 
        axis.text.x = element_text(size=12, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRCexplorer2 Transcript Counts", 
       y = "Density", x = "log2 Normalized Counts") +
  guides(fill = guide_legend(title = "Method", ncol = 2))

png("../results/supplementary_figure8c.png", width=275, height=200, units='mm', res = 600, pointsize=80)
ggarrange(p1, p2, ncol = 1, common.legend = TRUE, legend = "bottom")
dev.off()



# proportion of unique transcripts

# get list of circRNA transcripts detected
polyA_ciri_transcripts <- colnames(lung_polyA_ciri)
polyA_circ_transcripts <- colnames(lung_polyA_circ)
ribozero_ciri_transcripts <- colnames(lung_ribozero_ciri)
ribozero_circ_transcripts <- colnames(lung_ribozero_circ)

CIRI_intersect <- intersect(polyA_ciri_transcripts, ribozero_ciri_transcripts)
CIRC_intersect <- intersect(polyA_circ_transcripts, ribozero_circ_transcripts)

# record number of unique transcripts and number of transcripts detected by both polyA and ribozero
CIRI_unique <- data.frame(category = c("poly(A)", "RiboZero", "Both"),
                          count = c(length(polyA_ciri_transcripts) - length(CIRI_intersect),
                                    length(ribozero_ciri_transcripts) - length(CIRI_intersect),
                                    length(CIRI_intersect)))
CIRC_unique <- data.frame(category = c("poly(A)", "RiboZero", "Both"),
                          count = c(length(polyA_circ_transcripts) - length(CIRC_intersect),
                                    length(ribozero_circ_transcripts) - length(CIRC_intersect),
                                    length(CIRC_intersect)))

# compute percentages
CIRI_unique$fraction = CIRI_unique$count / sum(CIRI_unique$count)
CIRC_unique$fraction = CIRC_unique$count / sum(CIRC_unique$count)

# compute the cumulative percentages (top of each rectangle)
CIRI_unique$ymax <- cumsum(CIRI_unique$fraction)
CIRI_unique$ymin <- c(0, head(CIRI_unique$ymax, n = -1))

CIRC_unique$ymax <- cumsum(CIRC_unique$fraction)
CIRC_unique$ymin <- c(0, head(CIRC_unique$ymax, n = -1))

p1 <- ggplot(CIRI_unique, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect(color = "black") +
  coord_polar(theta="y") + 
  guides(fill = guide_legend(title = "Proportion of \nDetected Transcripts", ncol = 1)) + #make legend one column instead of 3
  scale_fill_manual(labels=c("Both poly(A) and RiboZero", "poly(A) Only", "RiboZero Only"), values = c("#3670A0", "#4CC5AB", "#392C57")) + 
  xlim(c(2, 4)) + # makes donut chart
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),    
        text = element_text(size = 15), 
        legend.key.size = unit(1, 'cm'),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, size = 17), 
        axis.text.x = element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  labs(title = "CIRI2 Unique Transcript Detection")

p2 <- ggplot(CIRC_unique, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=category)) +
  geom_rect(color = "black") +
  coord_polar(theta="y") + 
  scale_fill_manual(labels=c("Both poly(A) and RiboZero", "poly(A) Only", "RiboZero Only"), values = c("#3670A0", "#4CC5AB", "#392C57")) + 
  xlim(c(2, 4)) + # makes donut chart
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),    
        text = element_text(size = 15), 
        legend.key.size = unit(1, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 17), 
        axis.text.x = element_blank(), 
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) +
  labs(title = "CIRCexplorer2 Unique Transcript Detection")

png("../results/supplementary_figure8a.png", width=250, height=150, units='mm', res = 600, pointsize=80)
ggarrange(p1, p2, ncol = 2, common.legend = TRUE, legend = "bottom")
dev.off()