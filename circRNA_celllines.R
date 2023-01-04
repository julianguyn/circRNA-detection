### load packages ###
suppressMessages(library(ggplot2))
suppressMessages(library(reshape2))
suppressMessages(library(ggpubr))
suppressMessages(library(ggsignif))

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
cell_all <- read.csv(file = "/data/rnaseq_meta/cell_annotation_all.csv", na.strings=c("", " ", "NA"))

#read in gdsc cell annotations
gdsc <- read.csv(file = "/data/rnaseq_meta/gdsc_rnaseq_meta.txt", sep = "\t")
gdsc <- gdsc[which(!gdsc$Comment.SUBMITTED_FILE_NAME. == "15552_5.cram"),]
gdsc$cellid <- matchToIDTable(ids=gdsc$Source.Name, tbl = cell_all, column = "GDSC_rnaseq.cellid", returnColumn = "unique.cellid")
gdsc$files <- gsub(".cram","",gdsc$Comment.SUBMITTED_FILE_NAME.)
rownames(gdsc) <- gdsc$files


#########################################################################
##### 2. Summarize circRNA detection (summarizeCIRI, summarizeCIRC) #####
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
gdsc_ciri_summ <- summarizeCIRI(dir_path = "/data/processed_cellline/CIRI2/GDSC/result")
hansen_ciri_summ <- summarizeCIRI(dir_path = "/data/processed_cellline/CIRI2/hansen/result") #both Ribo-zero + RNAse-R
hansen_ciri_matched_summ <- summarizeCIRI(dir_path = "/data/processed_cellline/CIRI2/hansen_match/result") #poly-A equivalent to Hansen data
gdsc_RiboMinus_ciri_summ <- summarizeCIRI(dir_path = "/data/processed_cellline/CIRI2/GDSC_Ribo/result") #GDSC Ribo-Zero samples (HCT-116, HeLa)


#extract count df only
gdsc_ciri_counts <- gdsc_ciri_summ$circ_counts_df
hansen_ciri_counts <- hansen_ciri_summ$circ_counts_df
hansen_ciri_matched_counts <- hansen_ciri_matched_summ$circ_counts_df
gdsc_RiboMinus_ciri_counts <- gdsc_RiboMinus_ciri_summ$circ_counts_df

#clean up sample names (remove 'gdsc' strings from samples)
gdsc_ciri_counts$sample <- gsub("gdsc", "", gdsc_ciri_counts$sample)
rownames(gdsc_ciri_counts) <- gdsc_ciri_counts$sample

gdsc_RiboMinus_ciri_counts$sample <- c("HCT 116","HeLa")
rownames(gdsc_RiboMinus_ciri_counts) <- gdsc_RiboMinus_ciri_counts$sample


#get unique.cellids for samples and keep only cell lines with matched RiboMinus profiles
gdsc_ciri_counts$cellid <- gdsc$cellid[match(rownames(gdsc_ciri_counts), rownames(gdsc))]
gdsc_ciri_counts <- gdsc_ciri_counts[which(gdsc_ciri_counts$cellid %in% c("HCT 116","HeLa")),]


#create combined dataframe of CIRI2 counts for each replicate for poly(A) and RiboMinus GDSC
ciri_combined <- data.frame("poly-A"=as.numeric(gdsc_ciri_counts$count),
                            "RiboMinus"=as.numeric(gdsc_RiboMinus_ciri_counts$count))
rownames(ciri_combined) <- gdsc_ciri_counts$cellid


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


gdsc_circ_summ <- summarizeCIRCexplorer(dir_path = "/data/processed_cellline/CIRCexplorer2/GDSC/annotate")
hansen_circ_summ <- summarizeCIRCexplorer(dir_path = "/data/processed_cellline/CIRCexplorer2/hansen/annotate") #both Ribo-zero + RNAse-R
hansen_circ_matched_summ <- summarizeCIRCexplorer(dir_path = "/data/processed_cellline/CIRCexplorer2/hansen_match/annotate") #poly-A equivalent to Hansen data
gdsc_RiboMinus_circ_summ <- summarizeCIRCexplorer(dir_path = "/data/processed_cellline/CIRCexplorer2/GDSC_Ribo/annotate") #GDSC Ribo-Zero samples (HCT-116, HeLa)


#extract count df only
gdsc_circ_counts <- gdsc_circ_summ$circ_counts_df
hansen_circ_counts <- hansen_circ_summ$circ_counts_df
hansen_circ_matched_counts <- hansen_circ_matched_summ$circ_counts_df
gdsc_RiboMinus_circ_counts <- gdsc_RiboMinus_circ_summ$circ_counts_df

#clean up sample names (remove 'gdsc' strings from samples)
gdsc_circ_counts$sample <- sub(".*gdsc *(.*?) *_circular.*", "\\1", gdsc_circ_counts$sample)
rownames(gdsc_circ_counts) <- gdsc_circ_counts$sample

gdsc_RiboMinus_circ_counts$sample <- c("HCT 116","HeLa")
rownames(gdsc_RiboMinus_circ_counts) <- gdsc_RiboMinus_circ_counts$sample



#get unique.cellids for samples + order CCLE & GDSC & GDSC-ribo unique.cellid by gCSI unique.cellid
gdsc_circ_counts$cellid <- gdsc$cellid[match(rownames(gdsc_circ_counts), rownames(gdsc))]
gdsc_circ_counts <- gdsc_circ_counts[which(gdsc_circ_counts$cellid %in% c("HCT 116","HeLa")),]

#create combined dataframe of CIRCexplorer counts for each replicate across gCSI, CCLE, GDSC
circ_combined <- data.frame("poly-A"=as.numeric(gdsc_circ_counts$count),
                            "RiboMinus"=as.numeric(gdsc_RiboMinus_circ_counts$count))
rownames(circ_combined) <- gdsc_circ_counts$cellid


######################################################################################
##### 3. Validate circRNA's with matched RNAse-R sample (filterCIRI, filterCIRC) #####
######################################################################################

### CIRI2 (filerCIRI function) ###

filterCIRI <- function(nonRNAseR_dir, RNAseR_dir, suff){
  #"*neg.tsv$" <- for ribominus
  #"\\.tsv$"<- for poly-A selected data
  ciri_files <- list.files(path=nonRNAseR_dir, pattern = suff, full.names = T)
  circ_counts_df <- data.frame(matrix(ncol=2, nrow = length(ciri_files)))
  colnames(circ_counts_df) <- c("sample", "count")
  for (f in 1:length(ciri_files)) {
    sample <- read.table(file = ciri_files[f], 
                         sep = '\t', 
                         skip = 1,
                         header = FALSE)
    
    #sample_name <- gsub("\\..*","", ciri_files[f])
    sample_name <- gsub(".*/","", ciri_files[f])
    sample_name <- gsub("*.tsv", "", sample_name)
    match_name <- gsub("ccle|gcsi", "",ciri_files[f])
    #match_name <- gsub("\\..*","", match_name)
    match_name <- gsub(".*/","", match_name)
    match_name <- gsub("neg","", match_name)
    match_name <- gsub("*.tsv", "", match_name)
    
    hansen_match <- read.table(file = paste0(RNAseR_dir, match_name,"pos.tsv"), 
                               sep = '\t', 
                               skip = 1,
                               header = FALSE)
    
    colnames(sample) <- c("circRNA_ID","chr","circRNA_start","circRNA_end","junction_reads", "SM_MS_SMS", "non_junction_reads", "junction_reads_ratio", "circRNA_type", "gene_id", "strand","junction_reads_ID")
    colnames(hansen_match) <- c("circRNA_ID","chr","circRNA_start","circRNA_end","junction_reads", "SM_MS_SMS", "non_junction_reads", "junction_reads_ratio", "circRNA_type", "gene_id", "strand","junction_reads_ID")
    
    sample_match <- sample[which(sample$circRNA_type=="exon"),] #keep exonic circRNA's only
    sample_match <- sample_match[which(sample_match$junction_reads >= 2),] #keep only circRNA's with junction reads >= 2
    
    hansen_match <- hansen_match[which(hansen_match$circRNA_type=="exon"),] #keep exonic circRNA's only
    hansen_match <- hansen_match[which(hansen_match$junction_reads >= 2),]  #keep only circRNA's with junction reads >= 2 
    
    sample_match <- sample_match[which(sample_match$circRNA_ID %in% hansen_match$circRNA_ID),]
    hansen_match <- hansen_match[which(hansen_match$circRNA_ID %in% sample_match$circRNA_ID),]
    
    filtered <- sample_match[hansen_match$junction_reads/sample_match$junction_reads >= 1.5,] #keep only circRNA's with at least 1.5-fold enrichment in matching RNAse-R sample
    
    count <- sum(filtered$junction_reads) #total # of junction reads per sample
    circ_counts_df[f,] <- c(sample_name,count)
  }
  return(circ_counts_df)
} 

#validated circRNA's from 22Rv1 (CCLE), LNCaP(CCLE), PC3(CCLE/gCSI) using matched Hansen RNAse-R enriched samples
validated_polyA_ciri <- filterCIRI(nonRNAseR_dir = "/data/processed_cellline/CIRI2/hansen_match/result", 
                                   RNAseR_dir = "/data/processed_cellline/CIRI2/hansen/result/",
                                   suff = "\\.tsv$")

#validated circRNA's from 22Rv1, LNCaP, and PC3 Hansen Ribo-Zero samples using matched Hansen RNAse-R enriched samples
validated_RiboMinus_ciri <- filterCIRI(nonRNAseR_dir = "/data/processed_cellline/CIRI2/hansen/result", 
                                      RNAseR_dir = "/data/processed_cellline/CIRI2/hansen/result/",
                                      suff = "*neg.tsv$")

### CIRCexplorer2 (filerCIRC function) ###

filterCIRC <- function(nonRNAseR_dir, RNAseR_dir, suff){
  #"*neg*" <- for ribominus
  #"\\.txt$"<- for poly-A selected data
  circ_files <- list.files(path=nonRNAseR_dir,  pattern = suff, full.names = T, recursive = TRUE)
  circ_counts_df <- data.frame(matrix(ncol=2, nrow = length(circ_files)))
  colnames(circ_counts_df) <- c("sample", "count")
  for (f in 1:length(circ_files)) {
    sample <- read.table(file = circ_files[f], 
                         sep = '\t',
                         header = FALSE)
    
    sample_name <- gsub("_circularRNA_known.txt","", circ_files[f])
    sample_name <- gsub(".*/","", sample_name)
    match_name <- gsub("neg", "",sample_name)
    match_name <- gsub("ccle|gcsi", "",match_name)
    
    hansen_match <- read.table(file = paste0(RNAseR_dir, "/" ,match_name, "pos", "/", match_name,"pos", "_circularRNA_known.txt"),
                               sep = '\t',
                               header = FALSE)
    
    colnames(sample) <- c("chrom","start","end", "name","score","strand","thickStart","thickEnd","itemRgb","exonCount","exonSizes", "exonOffsets", "readNumber", "circType", "geneName", "isoformName", "index", "flankIntron")
    colnames(hansen_match) <- c("chrom","start","end", "name","score","strand","thickStart","thickEnd","itemRgb","exonCount","exonSizes", "exonOffsets", "readNumber", "circType", "geneName", "isoformName", "index", "flankIntron")
    
    sample_match <- sample[which(sample$circType=="circRNA"),]
    sample_match <- sample_match[which(sample_match$readNumber >= 2),] #keep only circRNA's with junction reads >= 2
    
    hansen_match <- hansen_match[which(hansen_match$circType=="circRNA"),]
    hansen_match <- hansen_match[which(hansen_match$readNumber >= 2),]
    
    sample_match$pos <- paste0(sample_match$chrom,sample_match$start,sample_match$end)
    hansen_match$pos <- paste0(hansen_match$chrom,hansen_match$start,hansen_match$end)
    
    sample_match <- sample_match[which(sample_match$pos %in% hansen_match$pos),]
    hansen_match <- hansen_match[which(hansen_match$pos %in% sample_match$pos),]
    hansen_match  <- hansen_match[match(sample_match$pos, hansen_match$pos),]
    
    filtered <- sample_match[hansen_match$readNumber/sample_match$readNumber >= 1.5,] #keep only circRNA's with at least 1.5-fold enrichment in matching RNAse-R sample
    
    count <- sum(filtered$readNumber)
    circ_counts_df[f,] <- c(sample_name,count)
  }
  return(circ_counts_df)
} 

#RNAse-R validated circRNA's from 22Rv1 (CCLE), LNCaP(CCLE), PC3(CCLE/gCSI) using matched Hansen RNAse-R enriched samples
validated_polyA_circ <- filterCIRC(nonRNAseR_dir = "/data/processed_cellline/CIRCexplorer2/hansen_match/annotate", 
                                   RNAseR_dir = "/data/processed_cellline/CIRCexplorer2/hansen/annotate",
                                   suff = "\\.txt$")


#RNAse-R validated circRNA's from 22Rv1, LNCaP, and PC3 Hansen Ribo-Zero samples using matched Hansen RNAse-R enriched samples
validated_RiboMinus_circ <- filterCIRC(nonRNAseR_dir = "/data/processed_cellline/CIRCexplorer2/hansen/annotate", 
                                      RNAseR_dir = "/data/processed_cellline/CIRCexplorer2/hansen/annotate",
                                      suff = "*neg*")

save(ciri_combined, circ_combined, 
     hansen_ciri_counts, hansen_ciri_matched_counts, validated_RiboMinus_ciri, validated_polyA_ciri,
     hansen_circ_counts, hansen_circ_matched_counts, validated_RiboMinus_circ, validated_polyA_circ,
     file = "../results/cell_lines_counts.RData")

####################
##### 4. Plots #####
####################

## Comparision of GDSC HCT 116 and HeLa poly(A)-selected and rRNA-depleted replicates ##
ciri_combined$replicate <- rownames(ciri_combined)
circ_combined$replicate <- rownames(circ_combined)

ciri_combined <- melt(ciri_combined)
circ_combined <- melt(circ_combined)

ciri_combined$tool <- c(rep("CIRI2", 4))
circ_combined$tool <- c(rep("CIRCexplorer2", 4))

gdsc_combined <- rbind(ciri_combined, circ_combined)
colnames(gdsc_combined) <- c("Replicate", "Method", "Count", "Tool")
gdsc_combined$Count <- log2(gdsc_combined$Count + 1)

# plot circRNA counts
png("../results/supplementary_figure5.png", width=145, height=100, units='mm', res = 600, pointsize=80)
ggplot(data = gdsc_combined, aes(x = Replicate, y = Count, fill = Method)) + 
  geom_bar(stat = "identity", size = 0.5, width = 0.5, position = position_dodge(), color = "black") +
  facet_grid(~factor(Tool, levels=c("CIRI2", "CIRCexplorer2"))) +
  scale_fill_manual(values=c("#8B7B96", "#71A2B6"), labels=c("poly(A)", "RiboMinus")) + 
  theme_classic() + 
  scale_y_continuous(limits = c(0, 18), expand=c(0,0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),    # keep panel borders (removed from theme_classic())
        text = element_text(size = 17), 
        legend.key.size = unit(0.5, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 30), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(x = "Cell line", y = "log2 Normailized \ncircRNA Counts")
dev.off()




## Comparision of CCLE poly(A)-selected and rRNA-depleted replicates before and after RNAse R treatment##

# compile data.frame of CIRI2 countsn  
hansen_ciri_counts <- hansen_ciri_counts[which(hansen_ciri_counts$sample %in% c("22Rv1neg", "LNCaPneg", "PC3neg")),]
hansen_ciri_matched_counts <- hansen_ciri_matched_counts[which(hansen_ciri_matched_counts$sample %in% c("22Rv1ccle", "LNCaPccle", "PC3ccle")),]
validated_polyA_ciri <- validated_polyA_ciri[which(validated_polyA_ciri$sample %in% c("22Rv1ccle", "LNCaPccle", "PC3ccle")),]

hansen_ciri <- rbind(hansen_ciri_counts, hansen_ciri_matched_counts, validated_RiboMinus_ciri, validated_polyA_ciri)
hansen_ciri$count <- as.numeric(hansen_ciri$count)

hansen_ciri$sample <- c(rep(c("22Rv1", "LNCaP", "PC3"), 4))
hansen_ciri$method <- c(rep(c(rep("RiboMinus", 3), rep("poly-A", 3)), 2))
hansen_ciri$validation <- c(rep("RNase R-", 6), rep("RNase R+", 6))


# compile data.frame of CIRCexplorer2 counts
hansen_circ_counts$sample <- gsub("_circularRNA_known", "", hansen_circ_counts$sample)
hansen_circ_matched_counts$sample <- gsub("_circularRNA_known", "", hansen_circ_matched_counts$sample)
hansen_circ_counts <- hansen_circ_counts[which(hansen_circ_counts$sample %in% c("22Rv1neg", "LNCaPneg", "PC3neg")),]
hansen_circ_matched_counts <- hansen_circ_matched_counts[which(hansen_circ_matched_counts$sample %in% c("22Rv1ccle", "LNCaPccle", "PC3ccle")),]
validated_polyA_circ <- validated_polyA_circ[which(validated_polyA_circ$sample %in% c("22Rv1ccle", "LNCaPccle", "PC3ccle")),]

hansen_circ <- rbind(hansen_circ_counts, hansen_circ_matched_counts, validated_RiboMinus_circ, validated_polyA_circ)
hansen_circ$count <- as.numeric(hansen_circ$count)

hansen_circ$sample <- c(rep(c("22Rv1", "LNCaP", "PC3"), 4))
hansen_circ$method <- c(rep(c(rep("RiboMinus", 3), rep("poly-A", 3)), 2))
hansen_circ$validation <- c(rep("RNase R-", 6), rep("RNase R+", 6))


# fold difference calculations
# (new value / original value) - 1
((hansen_ciri[which(hansen_ciri$validation == "RNase R-" & hansen_ciri$sample == "22Rv1" & hansen_ciri$method == "RiboMinus"),]$count) /
    (hansen_ciri[which(hansen_ciri$validation == "RNase R-" & hansen_ciri$sample == "22Rv1" & hansen_ciri$method == "poly-A"),]$count)) -1
((hansen_ciri[which(hansen_ciri$validation == "RNase R-" & hansen_ciri$sample == "LNCaP" & hansen_ciri$method == "RiboMinus"),]$count) /
    (hansen_ciri[which(hansen_ciri$validation == "RNase R-" & hansen_ciri$sample == "LNCaP" & hansen_ciri$method == "poly-A"),]$count)) -1
((hansen_ciri[which(hansen_ciri$validation == "RNase R-" & hansen_ciri$sample == "PC3" & hansen_ciri$method == "RiboMinus"),]$count) /
    (hansen_ciri[which(hansen_ciri$validation == "RNase R-" & hansen_ciri$sample == "PC3" & hansen_ciri$method == "poly-A"),]$count)) -1

((hansen_ciri[which(hansen_ciri$validation == "RNase R+" & hansen_ciri$sample == "22Rv1" & hansen_ciri$method == "RiboMinus"),]$count) /
    (hansen_ciri[which(hansen_ciri$validation == "RNase R+" & hansen_ciri$sample == "22Rv1" & hansen_ciri$method == "poly-A"),]$count)) -1
((hansen_ciri[which(hansen_ciri$validation == "RNase R+" & hansen_ciri$sample == "LNCaP" & hansen_ciri$method == "RiboMinus"),]$count) /
    (hansen_ciri[which(hansen_ciri$validation == "RNase R+" & hansen_ciri$sample == "LNCaP" & hansen_ciri$method == "poly-A"),]$count)) -1
((hansen_ciri[which(hansen_ciri$validation == "RNase R+" & hansen_ciri$sample == "PC3" & hansen_ciri$method == "RiboMinus"),]$count) /
    (hansen_ciri[which(hansen_ciri$validation == "RNase R+" & hansen_ciri$sample == "PC3" & hansen_ciri$method == "poly-A"),]$count)) -1


((hansen_circ[which(hansen_circ$validation == "RNase R-" & hansen_circ$sample == "22Rv1" & hansen_circ$method == "RiboMinus"),]$count) /
    (hansen_circ[which(hansen_circ$validation == "RNase R-" & hansen_circ$sample == "22Rv1" & hansen_circ$method == "poly-A"),]$count)) -1
((hansen_circ[which(hansen_circ$validation == "RNase R-" & hansen_circ$sample == "LNCaP" & hansen_circ$method == "RiboMinus"),]$count) /
    (hansen_circ[which(hansen_circ$validation == "RNase R-" & hansen_circ$sample == "LNCaP" & hansen_circ$method == "poly-A"),]$count)) -1
((hansen_circ[which(hansen_circ$validation == "RNase R-" & hansen_circ$sample == "PC3" & hansen_circ$method == "RiboMinus"),]$count) /
    (hansen_circ[which(hansen_circ$validation == "RNase R-" & hansen_circ$sample == "PC3" & hansen_circ$method == "poly-A"),]$count)) -1

((hansen_circ[which(hansen_circ$validation == "RNase R+" & hansen_circ$sample == "22Rv1" & hansen_circ$method == "RiboMinus"),]$count) /
    (hansen_circ[which(hansen_circ$validation == "RNase R+" & hansen_circ$sample == "22Rv1" & hansen_circ$method == "poly-A"),]$count)) -1
((hansen_circ[which(hansen_circ$validation == "RNase R+" & hansen_circ$sample == "LNCaP" & hansen_circ$method == "RiboMinus"),]$count) /
    (hansen_circ[which(hansen_circ$validation == "RNase R+" & hansen_circ$sample == "LNCaP" & hansen_circ$method == "poly-A"),]$count)) -1
((hansen_circ[which(hansen_circ$validation == "RNase R+" & hansen_circ$sample == "PC3" & hansen_circ$method == "RiboMinus"),]$count) /
    (hansen_circ[which(hansen_circ$validation == "RNase R+" & hansen_circ$sample == "PC3" & hansen_circ$method == "poly-A"),]$count)) -1


# annotations for geom_signif
hansen_ciri_annotations <- data.frame(validation = c("RNase R-", "RNase R+"), method = c("poly-A", "RiboMinus"),
                                      y = c(5297, 4291, 11079, 9700, 10913, 9542), start = c(0.85, 0.85, 1.85, 1.85, 2.85, 2.85),
                                      end = c(1.15, 1.15, 2.15, 2.15, 3.15, 3.15), labels = c("", "", "", "", "", ""))
hansen_ciri_labels <- data.frame(label = c("", "****", "*", "***", "*", "****"),
                                 cyl = c("RNase R-", "RNase R+"), validation = c("RNase R-", "RNase R+"), method = c("poly-A", "RiboMinus"))

hansen_circ_annotations <- data.frame(validation = c("RNase R-", "RNase R+"), method = c("poly-A", "RiboMinus"),
                                      y = c(9534, 5838, 13158, 10423, 14522, 11011), start = c(0.85, 0.85, 1.85, 1.85, 2.85, 2.85),
                                      end = c(1.15, 1.15, 2.15, 2.15, 3.15, 3.15), labels = c("", "", "", "", "", ""))
hansen_circ_labels <- data.frame(label = c("", "**", "", "***", "", "**"),
                                 cyl = c("RNase R-", "RNase R+"), validation = c("RNase R-", "RNase R+"), method = c("poly-A", "RiboMinus"))


# plot circRNA counts
p1 <- ggplot(hansen_ciri, aes(x=sample, y=count, fill=method)) + 
  geom_bar(stat="identity", size = 0.5, width = 0.7, position = position_dodge(), color = "black") + 
  geom_signif(data = hansen_ciri_annotations,
              aes(y_position = y, 
                  xmin = start,
                  xmax = end, 
                  annotations = labels),
              tip_length = c(0.05, 0.01),
              manual = TRUE) +
  geom_text(data = hansen_ciri_labels,
            mapping = aes(x = c(0.9, 0.75, 1.9, 1.8, 2.9, 2.75),
                          y = c(5297, 4291, 11079, 9700, 10913, 9542),
                          label = label),
            hjust = -0.1, vjust = -1) +
  facet_grid(~validation) + 
  scale_y_continuous(limits = c(0, 15500), expand=c(0,0)) +
  scale_fill_manual("Method",values=c("#8B7B96", "#71A2B6"), limits=c("poly-A", "RiboMinus"), labels=c("poly(A)", "RiboMinus")) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        text = element_text(size = 17), 
        legend.key.size = unit(0.4, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5, angle = 45), 
        axis.text.y = element_text(size=15), 
        legend.text = element_text(size=15)) +
  labs(x = "Cell line", y = "No. of circRNAs", title = "CIRI2 circRNA Counts")

p2 <- ggplot(hansen_circ, aes(x=sample, y=count, fill=method)) + 
  geom_bar(stat="identity", size = 0.5, width = 0.7, position = position_dodge(), color = "black") + 
  geom_signif(data = hansen_circ_annotations,
              aes(y_position = y, 
                  xmin = start,
                  xmax = end, 
                  annotations = labels),
              tip_length = c(0.05, 0.01),
              manual = TRUE) +
  geom_text(data = hansen_circ_labels,
            mapping = aes(x = c(0.9, 0.85, 1.9, 1.8, 2.9, 2.85),
                          y = c(9534, 5838, 13158, 10423, 14522, 11011),
                          label = label),
            hjust = -0.1, vjust = -1) +
  facet_grid(~validation) + 
  scale_y_continuous(limits = c(0, 15500), expand=c(0,0)) +
  scale_fill_manual("Method",values=c("#8B7B96", "#71A2B6"), limits=c("poly-A", "RiboMinus"), labels=c("poly(A)", "RiboMinus")) + 
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        text = element_text(size = 17), 
        legend.key.size = unit(0.4, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 18), 
        axis.text.x = element_text(size=15, vjust = 0.5, angle = 45), 
        axis.text.y = element_text(size=15), 
        legend.text = element_text(size=15)) +
  labs(x = "Cell line", y = "No. of circRNAs", title = "CIRCexplorer2 circRNA Counts")

png("../results/supplementary_figure6.png", width=150, height=200, units='mm', res = 600, pointsize=80)
ggarrange(p1, p2, ncol = 1, nrow = 2)
dev.off()
