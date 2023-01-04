### load packages ###

suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(viridis))
suppressMessages(library(ggpubr))

options(stringsAsFactors = F)


#########################################################################
##### 1. Summarize circRNA detection (summarizeCIRI, summarizeCIRC) #####
#########################################################################

### CIRI2 (summarizeCIRI function) ###

summarizeCIRI <- function(dir_path){
  
  circRNA_ids <- c() #holds all circRNA_ids for dataset
  ciri_files <- list.files(path=dir_path,  pattern = "\\.tsv$", full.names = T) 
  circ_counts_df <- data.frame(matrix(ncol=2, nrow = length(ciri_files)))
  colnames(circ_counts_df) <- c("sample", "count")
  for (f in 1:length(ciri_files)) {
    sample <- read.table(file = ciri_files[f], 
                         sep = '\t', 
                         skip = 1,
                         header = FALSE)
    
    colnames(sample) <- c("circRNA_ID","chr","circRNA_start","circRNA_end","junction_reads", "SM_MS_SMS", "non_junction_reads", "junction_reads_ratio", "circRNA_type", "gene_id", "strand","junction_reads_ID")
    #sample_name <- gsub("\\..*","", ciri_files[f])
    sample_name <- gsub(".*/","", ciri_files[f])
    sample_name <- gsub("*.tsv", "", sample_name)
    
    sample <- sample[which(sample$circRNA_type=="exon"),] #keep only exonic circRNA's
    sample <- sample[which(sample$junction_reads >= 2),] #keep only circRNA's with >= 2 junction reads
    
    #convert coordinate from 1-based to 0-based to match CIRCexplorer2 coordinates
    sample$circRNA_start <- sample$circRNA_start - 1
    sample$circRNA_ID <- paste0(sample$chr, ":", sample$circRNA_start, "|", sample$circRNA_end)
    
    circRNA_ids <- c(sample$circRNA_ID, circRNA_ids)
    count <- sum(sample$junction_reads) #total number of junction reads per sample
    circ_counts_df[f,] <- c(sample_name,count)
  }
  xx <- list(circ_counts_df, circRNA_ids)
  names(xx)[1] <- "circ_counts_df"
  names(xx)[2] <- "circRNA_ids"
  return(xx)
} 


#run summarizeCIRI function
lung_ribozero_ciri <- summarizeCIRI(dir_path = "/data/processed_lung/CIRI2/RiboZero/result")
lung_polyA_ciri <- summarizeCIRI(dir_path = "/data/processed_lung/CIRI2/polyA/result")

save(lung_ribozero_ciri, 
     lung_polyA_ciri, 
     file="../results/summarizeCIRI_lung.RData")

# read in lung metadata and match RiboZero with polyA samples
lung_ribozero_meta <- read.csv("/data/rnaseq_meta/lung_ribozero.csv")
lung_polyA_meta <- read.csv("/data/rnaseq_meta/lung_polyA.tsv", sep = "\t")

tt <- lung_polyA_meta$sample[match(lung_ribozero_meta$TB_id, lung_polyA_meta$TB_id)]
lung_ribozero_meta$polyA <- tt

rownames(lung_ribozero_ciri$circ_counts_df) <- lung_ribozero_ciri$circ_counts_df$sample
rownames(lung_polyA_ciri$circ_counts_df) <- lung_polyA_ciri$circ_counts_df$sample

lung_ribozero_ciri$circ_counts_df <- lung_ribozero_ciri$circ_counts_df[lung_ribozero_meta$helab_id,]
lung_polyA_ciri$circ_counts_df <- lung_polyA_ciri$circ_counts_df[lung_ribozero_meta$polyA,]

combined_lung_ciri_df <- lung_ribozero_ciri$circ_counts_df
combined_lung_ciri_df$bhk_id <- paste0("tumor",1:51)
combined_lung_ciri_df$count_polyA <- lung_polyA_ciri$circ_counts_df$count



### CIRCexplorer2 (summarizeCIRC function) ###

summarizeCIRCexplorer <- function(dir_path){
  
  circRNA_ids <- c() #holds all circRNA_ids for dataset
  circ_files <- list.files(path=dir_path,  pattern = ".txt", full.names = TRUE, recursive = TRUE)
  circ_counts_df <- data.frame(matrix(ncol=2, nrow = length(circ_files)))
  colnames(circ_counts_df) <- c("sample", "count")
  for (f in 1:length(circ_files)) {
    sample <- read.table(file = circ_files[f], 
                         sep = '\t',
                         header = FALSE)
    
    colnames(sample) <- c("chrom","start","end", "name","score","strand","thickStart","thickEnd","itemRgb","exonCount","exonSizes", "exonOffsets", "readNumber", "circType", "geneName", "isoformName", "index", "flankIntron")
    
    sample_name <- sub(".*/ *(.*?) *.txt.*", "\\1", circ_files[f])
    sample <- sample[which(sample$circType=="circRNA"),] #keep only exonic circRNA's
    sample <- sample[which(sample$readNumber >= 2),] #keep only circRNA's with junction reads >= 2
    
    sample$circRNA_ID <- paste0(sample$chrom,":",sample$start,"|",sample$end) #create circRNA_id using chr, start & end coordinates
    circRNA_ids <- c(sample$circRNA_ID, circRNA_ids)
    count <- sum(sample$readNumber) #total # of junction reads per sample
    circ_counts_df[f,] <- c(sample_name,count)
  }
  xx <- list(circ_counts_df, circRNA_ids)
  names(xx)[1] <- "circ_counts_df"
  names(xx)[2] <- "circRNA_ids"
  return(xx)
} 

lung_ribozero_circ <- summarizeCIRCexplorer(dir_path = "/data/processed_lung/CIRCexplorer2/RiboZero/annotate")
lung_polyA_circ <- summarizeCIRCexplorer(dir_path = "/data/processed_lung/CIRCexplorer2/polyA/annotate")


save(lung_ribozero_circ, 
     lung_polyA_circ, 
     file="../results/summarizeCIRC_lung.RData")


# match RiboZero with polyA samples
rownames(lung_ribozero_circ$circ_counts_df) <- lung_ribozero_circ$circ_counts_df$sample
rownames(lung_polyA_circ$circ_counts_df) <- lung_polyA_circ$circ_counts_df$sample

lung_ribozero_circ$circ_counts_df <- lung_ribozero_circ$circ_counts_df[lung_ribozero_meta$helab_id,]
lung_polyA_circ$circ_counts_df <- lung_polyA_circ$circ_counts_df[lung_ribozero_meta$polyA,]

combined_lung_circ_df <- lung_ribozero_circ$circ_counts_df
combined_lung_circ_df$bhk_id <- paste0("tumor",1:51)
combined_lung_circ_df$count_polyA <- lung_polyA_circ$circ_counts_df$count


# format data.frame objects
combined_lung_ciri_df$count <- as.numeric(combined_lung_ciri_df$count)
combined_lung_ciri_df$count_polyA <- as.numeric(combined_lung_ciri_df$count_polyA)
rownames(combined_lung_ciri_df) <- combined_lung_ciri_df$bhk_id
ciri_combined <- data.frame("rRNA-depleted"=as.numeric(combined_lung_ciri_df$count), 
                            "poly-A"=as.numeric(combined_lung_ciri_df$count_polyA))
rownames(ciri_combined) <- rownames(combined_lung_ciri_df)


combined_lung_circ_df$count <- as.numeric(combined_lung_circ_df$count)
combined_lung_circ_df$count_polyA <- as.numeric(combined_lung_circ_df$count_polyA)
rownames(combined_lung_circ_df) <- combined_lung_circ_df$bhk_id
circ_combined <- data.frame("rRNA-depleted"=as.numeric(combined_lung_circ_df$count), 
                            "poly-A"=as.numeric(combined_lung_circ_df$count_polyA))

rownames(circ_combined) <- rownames(combined_lung_circ_df)


####################
##### 2. Plots #####
####################

tumourID <- c(rownames(ciri_combined), rownames(circ_combined))
method <- c(rep("CIRI2", 51), rep("CIRCexplorer2", 51))
lung_combined <- rbind(ciri_combined, circ_combined)
lung_combined$tumourID <- tumourID
lung_combined$method <- method
lung_combined <- melt(lung_combined)
lung_combined$value <- log2(lung_combined$value + 1)

png("../results/figure2.png", width=125, height=200, units='mm', res = 600, pointsize=80)
ggplot(lung_combined, 
       aes(x = variable, y = tumourID, fill = value)) + 
  geom_tile() + 
  facet_grid(~factor(method, levels=c("CIRI2", "CIRCexplorer2"))) +
  theme_void() +
  theme(text = element_text(size = 17), 
        axis.text.x = element_text(size=15, vjust = 0.5, angle = 45),
        axis.title.x = element_text(size=17),
        axis.title.y = element_text(size=17, angle = 90),
        legend.title=element_text(size=12),
        legend.text=element_text(size=12)) +
  scale_fill_viridis(limits=c(0, 17), 
                     option="mako", 
                     direction = -1) +
  scale_x_discrete(labels = c("rRNA-\ndepleted", "poly(A)")) +
  guides(fill = guide_colourbar(barwidth = 0.2, 
                                barheight = 3, 
                                title = "log2 Normalized \ncircRNA Counts")) +
  labs(x = "Method", y = "Tumour Sample")
dev.off()


#################################
##### 3. Unique Transcripts #####
#################################

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



### format data.frames for plots ###
CIRI_counts <- melt(CIRI_counts)
CIRC_counts <- melt(CIRC_counts)

colnames(CIRI_counts) <- c("transcript_id", "Method", "Count")
colnames(CIRC_counts) <- c("transcript_id", "Method", "Count")

# log2 + 1 normalize circRNA counts
CIRI_counts$Count <- log2(CIRI_counts$Count + 1)
CIRC_counts$Count <- log2(CIRC_counts$Count + 1)


#### proportion of unique transcripts ####

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



# plot bar plots
p1 <- ggplot(data = CIRI_counts, aes(x = reorder(transcript_id, -Count), y = Count, fill = Method, alpha = Method)) + 
  geom_bar(stat = "identity", position = "identity", width = 1, color = "black") +
  scale_y_continuous(limits = c(0, 9), expand=c(0,0)) +
  scale_fill_manual(values=c("#4CC5AB", "#392C57"), labels=c("poly(A)", "RiboZero")) +
  scale_alpha_manual(values=c(0.8, 0.6), guide = "none") +
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
  geom_bar(stat = "identity", position = "identity", width = 1, color = "black") +
  scale_y_continuous(limits = c(0, 9), expand=c(0,0)) +
  scale_fill_manual(values=c("#4CC5AB", "#392C57"), labels=c("poly(A)", "RiboZero")) +
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

