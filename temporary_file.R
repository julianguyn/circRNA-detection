suppressMessages(library(edgeR))
suppressMessages(library(reshape2))
suppressMessages(library(dplyr))
suppressMessages(library(plyr))
suppressMessages(library(ggplot2))
suppressMessages(library(DESeq2))

options(stringsAsFactors = F)

#match cell.id to unique.cellid from PharmacoGx
matchToIDTable <- function(ids,tbl, column, returnColumn="unique.cellid") {
  sapply(ids, function(x) {
    myx <- grep(paste0("((///)|^)",Hmisc::escapeRegex(x),"((///)|$)"), tbl[,column])
    if(length(myx) > 1){
      stop("Something went wrong in curating ids, we have multiple matches")
    }
    if(length(myx) == 0){return(NA_character_)}
    return(tbl[myx, returnColumn])
  })
}

##############################################################
##### 1. Determine biological replicates across datasets #####
##############################################################

load("../results/inter.RData")

#read in current cell annotations
cell_all <- read.csv(file = "../data/rnaseq_meta/cell_annotation_all.csv", na.strings=c("", " ", "NA"))

#read in gdsc cell annotations
gdsc <- read.csv(file = "../data/rnaseq_meta/gdsc_rnaseq_meta.txt", sep = "\t")
gdsc <- gdsc[which(!gdsc$Comment.SUBMITTED_FILE_NAME. == "15552_5.cram"),]
gdsc$cellid <- matchToIDTable(ids=gdsc$Source.Name, tbl = cell_all, column = "GDSC_rnaseq.cellid", returnColumn = "unique.cellid")
gdsc$files <- gsub(".cram","",gdsc$Comment.SUBMITTED_FILE_NAME.)
rownames(gdsc) <- gdsc$files


circRNA_counts_CIRI <- function(dir_path){
  
  sample_names <- c() #holds all sample names
  circRNA_counts_df <- data.frame(matrix(ncol=0, nrow=0))
  pairs_df <- data.frame(matrix(ncol=0, nrow=0))
  
  ciri_files <- list.files(path=dir_path,  pattern = "\\.tsv$", full.names = T) 
  for (f in 1:length(ciri_files)) {
    sample <- read.table(file = ciri_files[f], 
                         sep = '\t', 
                         skip = 1,
                         header = FALSE)
    # store sample names
    #sample_name <- gsub("\\..*","", ciri_files[f])
    sample_name <- gsub(".*/","", ciri_files[f])
    sample_name <- gsub("*.tsv", "", sample_name)
    sample_names <- c(sample_names, sample_name)
    
    # keep only exonic circRNAs and circRNAs with >=2 junction reads
    colnames(sample) <- c("circRNA_ID","chr","circRNA_start","circRNA_end","junction_reads", "SM_MS_SMS", "non_junction_reads", "junction_reads_ratio", "circRNA_type", "gene_id", "strand","junction_reads_ID")
    sample <- sample[which(sample$circRNA_type=="exon"),] 
    sample <- sample[which(sample$junction_reads >= 2),] 
    
    #convert coordinate from 1-based to 0-based to match CIRCexplorer2 coordinates
    sample$circRNA_start <- sample$circRNA_start - 1
    sample$circRNA_ID <- paste0(sample$chr, ":", sample$circRNA_start, "|", sample$circRNA_end)
    
    #keep only # junction reads (counts), gene ID, and circRNA ID
    sample <- sample[,c("circRNA_ID", "junction_reads", "gene_id")]
    
    #keep gene - circRNA pairings
    pairs <- data.frame("circRNA_id" = sample$circRNA_ID, "gene_id" = sample$gene_id)
    
    # format counts table
    sample <- t(sample)
    colnames(sample) <- sample[1,]
    sample <- as.data.frame(sample)
    sample <- sample[-c(1,3),]
    
    circRNA_counts_df <- rbind.fill(circRNA_counts_df, sample)
    
    if (f > 1) {
      pairs <- pairs[-which(pairs$circRNA_id %in% pairs_df$circRNA_id),] 
    }
    pairs_df <- rbind(pairs_df, pairs)
    
  }
  circRNA_counts_df <- as.data.frame(circRNA_counts_df)
  circRNA_counts_df[] <- sapply(circRNA_counts_df, as.integer)
  rownames(circRNA_counts_df) <- sample_names 
  
  circRNA_counts_df <- as.data.frame(t(circRNA_counts_df))
  circRNA_counts_df$gene_id <- pairs_df$gene_id
  return(circRNA_counts_df)
} 



circRNA_counts_CIRCexplorer <- function(dir_path){
  
  sample_names <- c() #holds all sample names
  circRNA_counts_df <- data.frame(matrix(ncol=0, nrow = 0))
  pairs_df <- data.frame(matrix(ncol=0, nrow=0))
  
  circ_files <- list.files(path=dir_path,  pattern = ".txt", full.names = TRUE, recursive = TRUE)
  for (f in 1:length(circ_files)) {
    sample <- read.table(file = circ_files[f], 
                         sep = '\t',
                         header = FALSE)
    # store sample names
    sample_name <- sub(".*/ *(.*?) *.txt.*", "\\1", circ_files[f])
    sample_names <- c(sample_names, sample_name)
    
    # keep only exonic circRNAs and circRNAs with >=2 junction reads
    colnames(sample) <- c("chrom","start","end", "name","score","strand","thickStart","thickEnd","itemRgb","exonCount","exonSizes", "exonOffsets", "readNumber", "circType", "geneName", "isoformName", "index", "flankIntron")
    sample <- sample[which(sample$circType=="circRNA"),] 
    sample <- sample[which(sample$readNumber >= 2),] 
    
    sample$circRNA_ID <- paste0(sample$chrom,":",sample$start,"|",sample$end) #create circRNA_id using chr, start & end coordinates
    
    # keep only #reads (counts), gene IDs and circRNA IDs
    sample <- sample[,c("readNumber", "circRNA_ID", "geneName")]
    
    #keep gene - circRNA pairings
    pairs <- data.frame("circRNA_id" = sample$circRNA_ID, "gene_id" = sample$geneName)
    
    # format counts table
    sample <- t(sample)
    colnames(sample) <- sample[2,]
    sample <- as.data.frame(sample)
    sample <- sample[-c(2,3),]
    
    
    circRNA_counts_df <- rbind.fill(circRNA_counts_df, sample)
    
    if (f > 1) {
      pairs <- pairs[-which(pairs$circRNA_id %in% pairs_df$circRNA_id),] 
    }
    pairs_df <- rbind(pairs_df, pairs)
    
  }
  circRNA_counts_df <- as.data.frame(circRNA_counts_df)
  circRNA_counts_df[] <- sapply(circRNA_counts_df, as.integer)
  rownames(circRNA_counts_df) <- sample_names 
  
  circRNA_counts_df <- as.data.frame(t(circRNA_counts_df))
  circRNA_counts_df$gene_id <- pairs_df$gene_id
  return(circRNA_counts_df)
} 

ciri_gdsc_polyA <- circRNA_counts_CIRI(dir_path = "../data/processed_cellline/CIRI2/GDSC/result")
ciri_gdsc_ribozero <- circRNA_counts_CIRI(dir_path = "../data/processed_cellline/CIRI2/GDSC_Ribo/result") #GDSC Ribo-Zero samples (HCT-116, HeLa)

circ_gdsc_polyA <- circRNA_counts_CIRCexplorer(dir_path = "../data/processed_cellline/CIRCexplorer2/GDSC/annotate")
circ_gdsc_ribozero <- circRNA_counts_CIRCexplorer(dir_path = "../data/processed_cellline/CIRCexplorer2/GDSC_Ribo/annotate") #GDSC Ribo-Zero samples (HCT-116, HeLa)


#clean up sample names (remove 'gdsc' strings from samples)
colnames(ciri_gdsc_polyA)[1:48] <- gsub("gdsc*", "\\1", colnames(ciri_gdsc_polyA))[1:48]
colnames(ciri_gdsc_ribozero) <- c("HCT 116","HeLa", "gene_id")
colnames(circ_gdsc_polyA)[1:48] <- gsub("gdsc*", "\\1", colnames(circ_gdsc_polyA))[1:48]
colnames(circ_gdsc_polyA)[1:48] <- gsub("*_circularRNA_known", "\\1", colnames(circ_gdsc_polyA))[1:48]
colnames(circ_gdsc_ribozero) <- c("HCT 116","HeLa", "gene_id")

#get unique.cellids for samples + order CCLE & GDSC & GDSC-ribo unique.cellid by gCSI unique.cellid
colnames(ciri_gdsc_polyA)[1:48] <- gdsc$cellid[match(colnames(ciri_gdsc_polyA)[1:48], rownames(gdsc))]
ciri_gdsc_polyA <- ciri_gdsc_polyA[,c("HCT 116","HeLa", "gene_id")]

colnames(circ_gdsc_polyA)[1:48] <- gdsc$cellid[match(colnames(circ_gdsc_polyA)[1:48], rownames(gdsc))]
circ_gdsc_polyA <- circ_gdsc_polyA[,c("HCT 116","HeLa", "gene_id")]


## progress checkpoint
save(ciri_gdsc_polyA, ciri_gdsc_ribozero, circ_gdsc_polyA, circ_gdsc_ribozero,
     file = "../temporary_file.RData")



load("../temporary_file.RData")

# load rnaseq.counts
load("../data/2019-01-03_circData_STAR.rda")
load("../results/gene_expression.RData")

##un-normalize counts 
#gdsc_counts <- (2^(expr_gdsc_p)) - 1

# use linear library size to calculate normalization factors
linear_library <- DGEList(counts=expr_gdsc_p)
linear_normfactor <- calcNormFactors(linear_library, method="TMM")

# simplify object to just be the samples table
linear_library <- linear_library$samples
linear_normfactor <- linear_normfactor$samples

# get circRNA DGEList object
ciri_gdsc_polyA <- ciri_gdsc_polyA[,1:2]
ciri_gdsc_polyA[is.na(ciri_gdsc_polyA)] <- 0
#ciri_gdsc_polyA <- (ciri_gdsc_polyA + 1)^2
ciri_gdsc_polyA <- DGEList(counts = ciri_gdsc_polyA)

ciri_gdsc_ribozero <- ciri_gdsc_ribozero[,1:2]
ciri_gdsc_ribozero[is.na(ciri_gdsc_ribozero)] <- 0
#ciri_gdsc_ribozero <- (ciri_gdsc_ribozero + 1)^2
ciri_gdsc_ribozero <- DGEList(counts = ciri_gdsc_ribozero)

circ_gdsc_polyA <- circ_gdsc_polyA[,1:2]
circ_gdsc_polyA[is.na(circ_gdsc_polyA)] <- 0
#circ_gdsc_polyA <- (circ_gdsc_polyA + 1)^2
circ_gdsc_polyA <- DGEList(counts = circ_gdsc_polyA)

circ_gdsc_ribozero <- circ_gdsc_ribozero[,1:2]
circ_gdsc_ribozero[is.na(circ_gdsc_ribozero)] <- 0
#circ_gdsc_ribozero <- (circ_gdsc_ribozero + 1)^2
circ_gdsc_ribozero <- DGEList(counts = circ_gdsc_ribozero)


# set library size and normalization factor from linear counts
ciri_gdsc_polyA_normfactor <- linear_normfactor[match(colnames(ciri_gdsc_polyA),rownames(linear_normfactor)), ]
ciri_gdsc_ribozero_normfactor <- linear_normfactor[match(colnames(ciri_gdsc_ribozero),rownames(linear_normfactor)), ]
circ_gdsc_polyA_normfactor <- linear_normfactor[match(colnames(circ_gdsc_polyA),rownames(linear_normfactor)), ]
circ_gdsc_ribozero_normfactor <- linear_normfactor[match(colnames(circ_gdsc_ribozero),rownames(linear_normfactor)), ]


ciri_gdsc_polyA$samples$lib.size <- ciri_gdsc_polyA_normfactor$lib.size
ciri_gdsc_polyA$samples$norm.factors <- ciri_gdsc_polyA_normfactor$norm.factors

ciri_gdsc_ribozero$samples$lib.size <- ciri_gdsc_ribozero_normfactor$lib.size
ciri_gdsc_ribozero$samples$norm.factors <- ciri_gdsc_ribozero_normfactor$norm.factors

circ_gdsc_polyA$samples$lib.size <- circ_gdsc_polyA_normfactor$lib.size
circ_gdsc_polyA$samples$norm.factors <- circ_gdsc_polyA_normfactor$norm.factors

circ_gdsc_ribozero$samples$lib.size <- circ_gdsc_ribozero_normfactor$lib.size
circ_gdsc_ribozero$samples$norm.factors <- circ_gdsc_ribozero_normfactor$norm.factors

# normalize circRNA counts
ciri_gdsc_polyA_normalized <- data.frame(cpm(ciri_gdsc_polyA, normalized.lib.sizes = TRUE))
ciri_gdsc_ribozero_normalized <- data.frame(cpm(ciri_gdsc_ribozero, normalized.lib.sizes = TRUE))
circ_gdsc_polyA_normalized <- data.frame(cpm(circ_gdsc_polyA, normalized.lib.sizes = TRUE))
circ_gdsc_ribozero_normalized <- data.frame(cpm(circ_gdsc_ribozero, normalized.lib.sizes = TRUE))

# make data.frame of total circRNA counts
gdsc_combined <- data.frame(matrix(ncol=4, nrow = 0))
colnames(gdsc_combined) <- c("Replicate", "Count", "Method", "Tool")

gdsc_combined <- rbind(gdsc_combined, data.frame(Replicate = "HCT 116", Count = sum(ciri_gdsc_polyA_normalized$HCT.116), Method = "poly(A)", Tool = "CIRI2"))
gdsc_combined <- rbind(gdsc_combined, data.frame(Replicate = "HeLa", Count = sum(ciri_gdsc_polyA_normalized$HeLa), Method = "poly(A)", Tool = "CIRI2"))

gdsc_combined <- rbind(gdsc_combined, data.frame(Replicate = "HCT 116", Count = sum(ciri_gdsc_ribozero_normalized$HCT.116), Method = "RiboZero", Tool = "CIRI2"))
gdsc_combined <- rbind(gdsc_combined, data.frame(Replicate = "HeLa", Count = sum(ciri_gdsc_ribozero_normalized$HeLa), Method = "RiboZero", Tool = "CIRI2"))

gdsc_combined <- rbind(gdsc_combined, data.frame(Replicate = "HCT 116", Count = sum(circ_gdsc_polyA_normalized$HCT.116), Method = "poly(A)", Tool = "CIRCexplorer2"))
gdsc_combined <- rbind(gdsc_combined, data.frame(Replicate = "HeLa", Count = sum(circ_gdsc_polyA_normalized$HeLa), Method = "poly(A)", Tool = "CIRCexplorer2"))

gdsc_combined <- rbind(gdsc_combined, data.frame(Replicate = "HCT 116", Count = sum(circ_gdsc_ribozero_normalized$HCT.116), Method = "RiboZero", Tool = "CIRCexplorer2"))
gdsc_combined <- rbind(gdsc_combined, data.frame(Replicate = "HeLa", Count = sum(circ_gdsc_ribozero_normalized$HeLa), Method = "RiboZero", Tool = "CIRCexplorer2"))


# plot counts
png("../results/gdsc_genes_normalized.png", width=145, height=100, units='mm', res = 600, pointsize=80)
ggplot(data = gdsc_combined, aes(x = Replicate, y = Count, fill = Method)) + 
  geom_bar(stat = "identity", size = 0.5, width = 0.5, position = position_dodge(), color = "black") + 
  facet_grid(~factor(Tool, levels=c("CIRI2", "CIRCexplorer2"))) +
  scale_fill_manual(values=c("#8B7B96", "#71A2B6"), labels=c("poly(A)", "RiboZero")) + 
  theme_classic() + 
  #scale_y_continuous(limits = c(0, 520000), expand=c(0,0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),    # keep panel borders (removed from theme_classic())
        text = element_text(size = 17), 
        legend.key.size = unit(0.5, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 30), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(x = "Cell line", y = "No. of circRNAs")
dev.off()






###### change gene names
x <- gsub("_", "", ciri_gdsc_polyA$gene_id)
x <- gsub(",", "", x)
x <- gsub("\\.[0-9]*$","",x)

y <- gsub("\\.[0-9]*$","",rownames(linear_library$counts))
#

# create DGElist object for circRNA counts




linear_normfactor$genes$Length <- annotation[match(rownames(linear_normfactor$counts), annotation$Gene.stable.ID), ]$exon_length


#####################################
############### rpkm ################
#####################################

load("../temporary_file.RData")

# get gene lengths
annotations <- read.delim(file = "../data/processed_cellline/h38_exons.txt", header = TRUE)
annotations$exon_length <- with(annotations, Exon.region.end..bp. - Exon.region.start..bp.)
annotations <- annotations[,c("Gene.stable.ID", "Gene.name", "exon_length")]

# take average exon lengths per gene
annotations_ciri <- aggregate(exon_length ~ Gene.stable.ID, annotations, mean)
annotations_circ <- aggregate(exon_length ~ Gene.name, annotations, mean)


# get stable gene id from ciri results
ciri_gdsc_polyA$gene_id <- gsub("\\..*", "", ciri_gdsc_polyA$gene_id)
ciri_gdsc_ribozero$gene_id <- gsub("\\..*", "", ciri_gdsc_ribozero$gene_id)


# save gene lengths
genes <- c(annotations_ciri$exon_length[match(ciri_gdsc_polyA$gene_id, annotations_ciri$Gene.stable.ID)],
           annotations_circ$exon_length[match(circ_gdsc_polyA$gene_id, annotations_circ$Gene.name)],
           annotations_ciri$exon_length[match(ciri_gdsc_ribozero$gene_id, annotations_ciri$Gene.stable.ID)],
           annotations_circ$exon_length[match(circ_gdsc_ribozero$gene_id, annotations_circ$Gene.name)])


# merge circRNA counts into one data.frame
rownames(ciri_gdsc_polyA) <- NULL
rownames(circ_gdsc_polyA) <- NULL
rownames(ciri_gdsc_ribozero) <- NULL
rownames(circ_gdsc_ribozero) <- NULL

gdsc_combined <- rbind(ciri_gdsc_polyA[,1:2], circ_gdsc_polyA[,1:2], ciri_gdsc_ribozero[,1:2], circ_gdsc_ribozero[,1:2])
gdsc_combined[is.na(gdsc_combined)] <- 0


###############
##### tpm #####
###############

# divide each column of counts matrix by estimate of gene length
x <- gdsc_combined / genes
x[is.na(x)] <- 0

# get TPM 
gdsc_combined_tpm <- data.frame(t(t(x) * 1e6 / colSums(x)))


# annotate values and save sum per group
gdsc_combined_tpm$label <- c(rep("ciri_gdsc_polyA", 4020), rep("circ_gdsc_polyA", 4902), rep("ciri_gdsc_ribozero", 10335), rep("circ_gdsc_ribozero", 10580))
HCT_116_counts <- aggregate(HCT.116 ~ label, gdsc_combined_tpm, sum)
HeLa_counts <- aggregate(HeLa ~ label, gdsc_combined_tpm, sum)

# merge results into one data.frame
HCT_116_counts$HeLa <- HeLa_counts$HeLa
gdsc_combined <- HCT_116_counts
colnames(gdsc_combined) <- c("label", "HCT 116", "HeLa")

# add labels for plot
gdsc_combined$Tool <- c(rep("CIRCexplorer2", 2), rep("CIRI2", 2))
gdsc_combined$Method <- c("poly.A", "RiboZero", "poly.A", "RiboZero")

gdsc_combined <- melt(gdsc_combined)
gdsc_combined <- gdsc_combined[,-1]

colnames(gdsc_combined) <- c("Tool", "Method", "Replicate", "Count")


# fold difference calculations
# (new value / original value) - 1
((gdsc_combined[which(gdsc_combined$Tool == "CIRI2" & gdsc_combined$Replicate == "HCT 116" & gdsc_combined$Method == "RiboZero"),]$Count) /
    (gdsc_combined[which(gdsc_combined$Tool == "CIRI2" & gdsc_combined$Replicate == "HCT 116" & gdsc_combined$Method == "poly.A"),]$Count)) -1
((gdsc_combined[which(gdsc_combined$Tool == "CIRI2" & gdsc_combined$Replicate == "HeLa" & gdsc_combined$Method == "RiboZero"),]$Count) /
    (gdsc_combined[which(gdsc_combined$Tool == "CIRI2" & gdsc_combined$Replicate == "HeLa" & gdsc_combined$Method == "poly.A"),]$Count)) -1
((gdsc_combined[which(gdsc_combined$Tool == "CIRCexplorer2" & gdsc_combined$Replicate == "HCT 116" & gdsc_combined$Method == "RiboZero"),]$Count) /
    (gdsc_combined[which(gdsc_combined$Tool == "CIRCexplorer2" & gdsc_combined$Replicate == "HCT 116" & gdsc_combined$Method == "poly.A"),]$Count)) -1
((gdsc_combined[which(gdsc_combined$Tool == "CIRCexplorer2" & gdsc_combined$Replicate == "HeLa" & gdsc_combined$Method == "RiboZero"),]$Count) /
    (gdsc_combined[which(gdsc_combined$Tool == "CIRCexplorer2" & gdsc_combined$Replicate == "HeLa" & gdsc_combined$Method == "poly.A"),]$Count)) -1



# annotations for geom_signif
gdsc_annotations <- data.frame(Tool = c("CIRI2", "CIRCexplorer2"), Method = c("poly.A", "RiboZero"),
                               y = c(564152, 465248, 510918, 483577), start = c(0.85, 0.85, 1.85, 1.85),
                               end = c(1.15, 1.15, 2.15, 2.15),labels = c("", "", "", ""))
gdsc_labels <- data.frame(label = c("**", "*", "*", "*"),
                          cyl = c("CIRI2", "CIRCexplorer2"), Tool = c("CIRI2", "CIRCexplorer2"), Method = c("poly.A", "RiboZero"))

# plot circRNA counts
png("../results/aa_supplementary_figure5_tpm.png", width=145, height=100, units='mm', res = 600, pointsize=80)
ggplot(data = gdsc_combined, aes(x = Replicate, y = Count, fill = Method)) + 
  geom_bar(stat = "identity", size = 0.5, width = 0.5, position = position_dodge(), color = "black") + 
  geom_signif(data = gdsc_annotations,
              aes(y_position = y, 
                  xmin = start,
                  xmax = end, 
                  annotations = labels),
              tip_length = c(0.05, 0.01),
              manual = TRUE) +
  geom_text(data = gdsc_labels, 
            mapping = aes(x = c(0.9, 0.95, 1.95, 1.95),
                          y = c(544152, 445248, 490918, 463577),
                          label = label),
            hjust = -0.1, vjust = -1) +
  facet_grid(~factor(Tool, levels=c("CIRI2", "CIRCexplorer2"))) +
  scale_fill_manual(values=c("#8B7B96", "#71A2B6"), labels=c("poly(A)", "RiboZero")) + 
  theme_classic() + 
  scale_y_continuous(limits = c(0, 650000), expand=c(0,0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),    # keep panel borders (removed from theme_classic())
        text = element_text(size = 15), 
        legend.key.size = unit(0.5, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 30), 
        axis.text.x = element_text(size=12, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(x = "Cell line", y = "No. of circRNAs")
dev.off()


######################
##### tpm log2+1 #####
######################

gdsc_combined$Count <- log2(gdsc_combined$Count + 1)

# fold difference calculations
# (new value / original value) - 1
((gdsc_combined[which(gdsc_combined$Tool == "CIRI2" & gdsc_combined$Replicate == "HCT 116" & gdsc_combined$Method == "RiboZero"),]$Count) /
   (gdsc_combined[which(gdsc_combined$Tool == "CIRI2" & gdsc_combined$Replicate == "HCT 116" & gdsc_combined$Method == "poly.A"),]$Count)) -1
((gdsc_combined[which(gdsc_combined$Tool == "CIRI2" & gdsc_combined$Replicate == "HeLa" & gdsc_combined$Method == "RiboZero"),]$Count) /
    (gdsc_combined[which(gdsc_combined$Tool == "CIRI2" & gdsc_combined$Replicate == "HeLa" & gdsc_combined$Method == "poly.A"),]$Count)) -1
((gdsc_combined[which(gdsc_combined$Tool == "CIRCexplorer2" & gdsc_combined$Replicate == "HCT 116" & gdsc_combined$Method == "RiboZero"),]$Count) /
    (gdsc_combined[which(gdsc_combined$Tool == "CIRCexplorer2" & gdsc_combined$Replicate == "HCT 116" & gdsc_combined$Method == "poly.A"),]$Count)) -1
((gdsc_combined[which(gdsc_combined$Tool == "CIRCexplorer2" & gdsc_combined$Replicate == "HeLa" & gdsc_combined$Method == "RiboZero"),]$Count) /
    (gdsc_combined[which(gdsc_combined$Tool == "CIRCexplorer2" & gdsc_combined$Replicate == "HeLa" & gdsc_combined$Method == "poly.A"),]$Count)) -1


# plot circRNA counts
png("../results/aa_supplementary_figure5_tpm_log2.png", width=145, height=100, units='mm', res = 600, pointsize=80)
ggplot(data = gdsc_combined, aes(x = Replicate, y = Count, fill = Method)) + 
  geom_bar(stat = "identity", size = 0.5, width = 0.5, position = position_dodge(), color = "black") + 
  facet_grid(~factor(Tool, levels=c("CIRI2", "CIRCexplorer2"))) +
  scale_fill_manual(values=c("#8B7B96", "#71A2B6"), labels=c("poly(A)", "RiboZero")) + 
  theme_classic() + 
  scale_y_continuous(limits = c(0, 22), expand=c(0,0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),    # keep panel borders (removed from theme_classic())
        text = element_text(size = 15), 
        legend.key.size = unit(0.5, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 30), 
        axis.text.x = element_text(size=12, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(x = "Cell line", y = "No. of circRNAs")
dev.off()





################
##### rpkm #####
################


# get circRNA DGEList objects
gdsc_combined_DGE <- DGEList(counts = gdsc_combined)

# add gene lengths
gdsc_combined_DGE$genes$Length <- genes

# compute RPKM values
gdsc_combined_rpkm <- data.frame(rpkm(gdsc_combined_DGE))


# annotate values and save sum per group
gdsc_combined_rpkm$label <- c(rep("ciri_gdsc_polyA", 4020), rep("circ_gdsc_polyA", 4902), rep("ciri_gdsc_ribozero", 10335), rep("circ_gdsc_ribozero", 10580))
HCT_116_counts <- aggregate(HCT.116 ~ label, gdsc_combined_rpkm, sum)
HeLa_counts <- aggregate(HeLa ~ label, gdsc_combined_rpkm, sum)

# merge results into one data.frame
HCT_116_counts$HeLa <- HeLa_counts$HeLa
gdsc_combined <- HCT_116_counts
colnames(gdsc_combined) <- c("label", "HCT 116", "HeLa")


# add labels for plot
gdsc_combined$Tool <- c(rep("CIRCexplorer2", 2), rep("CIRI2", 2))
gdsc_combined$Method <- c("poly.A", "RiboZero", "poly.A", "RiboZero")

gdsc_combined <- melt(gdsc_combined)
gdsc_combined <- gdsc_combined[,-1]

colnames(gdsc_combined) <- c("Tool", "Method", "Replicate", "Count")


# fold difference calculations
# (new value / original value) - 1
((gdsc_combined[which(gdsc_combined$Tool == "CIRI2" & gdsc_combined$Replicate == "HCT 116" & gdsc_combined$Method == "RiboZero"),]$Count) /
    (gdsc_combined[which(gdsc_combined$Tool == "CIRI2" & gdsc_combined$Replicate == "HCT 116" & gdsc_combined$Method == "poly.A"),]$Count)) -1
((gdsc_combined[which(gdsc_combined$Tool == "CIRI2" & gdsc_combined$Replicate == "HeLa" & gdsc_combined$Method == "RiboZero"),]$Count) /
    (gdsc_combined[which(gdsc_combined$Tool == "CIRI2" & gdsc_combined$Replicate == "HeLa" & gdsc_combined$Method == "poly.A"),]$Count)) -1
((gdsc_combined[which(gdsc_combined$Tool == "CIRCexplorer2" & gdsc_combined$Replicate == "HCT 116" & gdsc_combined$Method == "RiboZero"),]$Count) /
    (gdsc_combined[which(gdsc_combined$Tool == "CIRCexplorer2" & gdsc_combined$Replicate == "HCT 116" & gdsc_combined$Method == "poly.A"),]$Count)) -1
((gdsc_combined[which(gdsc_combined$Tool == "CIRCexplorer2" & gdsc_combined$Replicate == "HeLa" & gdsc_combined$Method == "RiboZero"),]$Count) /
    (gdsc_combined[which(gdsc_combined$Tool == "CIRCexplorer2" & gdsc_combined$Replicate == "HeLa" & gdsc_combined$Method == "poly.A"),]$Count)) -1



# annotations for geom_signif
gdsc_annotations <- data.frame(Tool = c("CIRI2", "CIRCexplorer2"), Method = c("poly.A", "RiboZero"),
                               y = c(2475118, 2051525, 2118361, 2008508), start = c(0.85, 0.85, 1.85, 1.85),
                               end = c(1.15, 1.15, 2.15, 2.15),labels = c("", "", "", ""))
gdsc_labels <- data.frame(label = c("**", "*", "*", "*"),
                          cyl = c("CIRI2", "CIRCexplorer2"), Tool = c("CIRI2", "CIRCexplorer2"), Method = c("poly.A", "RiboZero"))

# plot circRNA counts
png("../results/aa_supplementary_figure5_rpkm.png", width=145, height=100, units='mm', res = 600, pointsize=80)
ggplot(data = gdsc_combined, aes(x = Replicate, y = Count, fill = Method)) + 
  geom_bar(stat = "identity", size = 0.5, width = 0.5, position = position_dodge(), color = "black") + 
  geom_signif(data = gdsc_annotations,
              aes(y_position = y, 
                  xmin = start,
                  xmax = end, 
                  annotations = labels),
              tip_length = c(0.05, 0.01),
              manual = TRUE) +
  geom_text(data = gdsc_labels, 
            mapping = aes(x = c(0.9, 0.95, 1.95, 1.95),
                          y = c(2365118, 1941525, 2008361, 1898508),
                          label = label),
            hjust = -0.1, vjust = -1) +
  facet_grid(~factor(Tool, levels=c("CIRI2", "CIRCexplorer2"))) +
  scale_fill_manual(values=c("#8B7B96", "#71A2B6"), labels=c("poly(A)", "RiboZero")) + 
  theme_classic() + 
  scale_y_continuous(limits = c(0, 2700000), expand=c(0,0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),    # keep panel borders (removed from theme_classic())
        text = element_text(size = 15), 
        legend.key.size = unit(0.5, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 30), 
        axis.text.x = element_text(size=12, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(x = "Cell line", y = "No. of circRNAs")
dev.off()


################################
### Unique gene distribution ###
################################

# ciri_geneHCT_polyA
ciri_geneHCT_polyA <- ciri_gdsc_polyA[, c("HCT 116", "gene_id")]
colnames(ciri_geneHCT_polyA) <- c("count", "gene")

ciri_geneHCT_polyA[is.na(ciri_geneHCT_polyA)] <- 0
ciri_geneHCT_polyA <- ciri_geneHCT_polyA[-which(ciri_geneHCT_polyA$count == 0),]

ggplot(data = ciri_geneHCT_polyA, aes(x = reorder(gene, -count), y = count)) + 
  geom_bar(stat = "identity", color = "black", "fill" = "grey") +
  scale_y_continuous(limits = c(0, 67), expand=c(0,0)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ggtitle("CIRI2 HCT 116 poly(A) Gene Frequency")


# circ_geneHCT_polyA
circ_geneHCT_polyA <- circ_gdsc_polyA[, c("HCT 116", "gene_id")]
colnames(circ_geneHCT_polyA) <- c("count", "gene")

circ_geneHCT_polyA[is.na(circ_geneHCT_polyA)] <- 0
circ_geneHCT_polyA <- circ_geneHCT_polyA[-which(circ_geneHCT_polyA$count == 0),]

ggplot(data = circ_geneHCT_polyA, aes(x = reorder(gene, -count), y = count)) + 
  geom_bar(stat = "identity", color = "black", "fill" = "grey") +
  scale_y_continuous(limits = c(0, 75), expand=c(0,0)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ggtitle("CIRCexplorer2 HCT 116 poly(A) Gene Frequency")


# ciri_geneHCT_ribozero
ciri_geneHCT_ribozero <- ciri_gdsc_ribozero[, c("HCT 116", "gene_id")]
colnames(ciri_geneHCT_ribozero) <- c("count", "gene")

ciri_geneHCT_ribozero[is.na(ciri_geneHCT_ribozero)] <- 0
ciri_geneHCT_ribozero <- ciri_geneHCT_ribozero[-which(ciri_geneHCT_ribozero$count == 0),]

ggplot(data = ciri_geneHCT_ribozero, aes(x = reorder(gene, -count), y = count)) + 
  geom_bar(stat = "identity", color = "black", "fill" = "grey") +
  scale_y_continuous(limits = c(0, 250), expand=c(0,0)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ggtitle("CIRI2 HCT 116 RiboZero Gene Frequency")


# circ_geneHCT_ribozero
circ_geneHCT_ribozero <- circ_gdsc_ribozero[, c("HCT 116", "gene_id")]
colnames(circ_geneHCT_ribozero) <- c("count", "gene")

circ_geneHCT_ribozero[is.na(circ_geneHCT_ribozero)] <- 0
circ_geneHCT_ribozero <- circ_geneHCT_ribozero[-which(circ_geneHCT_ribozero$count == 0),]

ggplot(data = circ_geneHCT_ribozero, aes(x = reorder(gene, -count), y = count)) + 
  geom_bar(stat = "identity", color = "black", "fill" = "grey") +
  scale_y_continuous(limits = c(0, 230), expand=c(0,0)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ggtitle("CIRCexplorer2 HCT 116 RiboZero Gene Frequency")






# ciri_geneHeLa_polyA
ciri_geneHeLa_polyA <- ciri_gdsc_polyA[, c("HeLa", "gene_id")]
colnames(ciri_geneHeLa_polyA) <- c("count", "gene")

ciri_geneHeLa_polyA[is.na(ciri_geneHeLa_polyA)] <- 0
ciri_geneHeLa_polyA <- ciri_geneHeLa_polyA[-which(ciri_geneHeLa_polyA$count == 0),]

ggplot(data = ciri_geneHeLa_polyA, aes(x = reorder(gene, -count), y = count)) + 
  geom_bar(stat = "identity", color = "black", "fill" = "grey") +
  scale_y_continuous(limits = c(0, 67), expand=c(0,0)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ggtitle("CIRI2 HeLa poly(A) Gene Frequency")


# circ_geneHeLa_polyA
circ_geneHeLa_polyA <- circ_gdsc_polyA[, c("HeLa", "gene_id")]
colnames(circ_geneHeLa_polyA) <- c("count", "gene")

circ_geneHeLa_polyA[is.na(circ_geneHeLa_polyA)] <- 0
circ_geneHeLa_polyA <- circ_geneHeLa_polyA[-which(circ_geneHeLa_polyA$count == 0),]

ggplot(data = circ_geneHeLa_polyA, aes(x = reorder(gene, -count), y = count)) + 
  geom_bar(stat = "identity", color = "black", "fill" = "grey") +
  scale_y_continuous(limits = c(0, 75), expand=c(0,0)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ggtitle("CIRCexplorer2 HeLa poly(A) Gene Frequency")


# ciri_geneHeLa_ribozero
ciri_geneHeLa_ribozero <- ciri_gdsc_ribozero[, c("HeLa", "gene_id")]
colnames(ciri_geneHeLa_ribozero) <- c("count", "gene")

ciri_geneHeLa_ribozero[is.na(ciri_geneHeLa_ribozero)] <- 0
ciri_geneHeLa_ribozero <- ciri_geneHeLa_ribozero[-which(ciri_geneHeLa_ribozero$count == 0),]

ggplot(data = ciri_geneHeLa_ribozero, aes(x = reorder(gene, -count), y = count)) + 
  geom_bar(stat = "identity", color = "black", "fill" = "grey") +
  scale_y_continuous(limits = c(0, 1800), expand=c(0,0)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ggtitle("CIRI2 HeLa RiboZero Gene Frequency")


# circ_geneHeLa_ribozero
circ_geneHeLa_ribozero <- circ_gdsc_ribozero[, c("HeLa", "gene_id")]
colnames(circ_geneHeLa_ribozero) <- c("count", "gene")

circ_geneHeLa_ribozero[is.na(circ_geneHeLa_ribozero)] <- 0
circ_geneHeLa_ribozero <- circ_geneHeLa_ribozero[-which(circ_geneHeLa_ribozero$count == 0),]

ggplot(data = circ_geneHeLa_ribozero, aes(x = reorder(gene, -count), y = count)) + 
  geom_bar(stat = "identity", color = "black", "fill" = "grey") +
  scale_y_continuous(limits = c(0, 1800), expand=c(0,0)) +
  theme_classic() +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ggtitle("CIRCexplorer2 HeLa RiboZero Gene Frequency")
