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

###########################################################
##### Determine biological replicates across datasets #####
###########################################################

load("../results/inter.RData")

#read in current cell annotations
cell_all <- read.csv(file = "../data/rnaseq_meta/cell_annotation_all.csv", na.strings=c("", " ", "NA"))

#read in gdsc cell annotations
gdsc <- read.csv(file = "../data/rnaseq_meta/gdsc_rnaseq_meta.txt", sep = "\t")
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
gdsc_ciri_summ <- summarizeCIRI(dir_path = "../data/processed_cellline/CIRI2/GDSC/result")
gdsc_ribozero_ciri_summ <- summarizeCIRI(dir_path = "../data/processed_cellline/CIRI2/GDSC_Ribo/result") #GDSC Ribo-Zero samples (HCT-116, HeLa)


#extract count df only
gdsc_ciri_counts <- gdsc_ciri_summ$circ_counts_df
gdsc_ribozero_ciri_counts <- gdsc_ribozero_ciri_summ$circ_counts_df

#clean up sample names (remove 'gdsc' strings from samples)
gdsc_ciri_counts$sample <- gsub("gdsc", "", gdsc_ciri_counts$sample)
rownames(gdsc_ciri_counts) <- gdsc_ciri_counts$sample

gdsc_ribozero_ciri_counts$sample <- c("HCT 116","HeLa")
rownames(gdsc_ribozero_ciri_counts) <- gdsc_ribozero_ciri_counts$sample


#get unique.cellids for samples and keep only cell lines with matched ribozero profiles
gdsc_ciri_counts$cellid <- gdsc$cellid[match(rownames(gdsc_ciri_counts), rownames(gdsc))]
gdsc_ciri_counts <- gdsc_ciri_counts[which(gdsc_ciri_counts$cellid %in% c("HCT 116","HeLa")),]


#create combined dataframe of CIRI2 counts for each replicate for poly(A) and ribozero GDSC
ciri_combined <- data.frame("poly-A"=as.numeric(gdsc_ciri_counts$count),
                            "RiboZero"=as.numeric(gdsc_ribozero_ciri_counts$count))
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


gdsc_circ_summ <- summarizeCIRCexplorer(dir_path = "../data/processed_cellline/CIRCexplorer2/GDSC/annotate")
gdsc_ribozero_circ_summ <- summarizeCIRCexplorer(dir_path = "../data/processed_cellline/CIRCexplorer2/GDSC_Ribo/annotate") #GDSC Ribo-Zero samples (HCT-116, HeLa)


#extract count df only
gdsc_circ_counts <- gdsc_circ_summ$circ_counts_df
gdsc_ribozero_circ_counts <- gdsc_ribozero_circ_summ$circ_counts_df

#clean up sample names (remove 'gdsc' strings from samples)
gdsc_circ_counts$sample <- sub(".*gdsc *(.*?) *_circular.*", "\\1", gdsc_circ_counts$sample)
rownames(gdsc_circ_counts) <- gdsc_circ_counts$sample

gdsc_ribozero_circ_counts$sample <- c("HCT 116","HeLa")
rownames(gdsc_ribozero_circ_counts) <- gdsc_ribozero_circ_counts$sample



#get unique.cellids for samples + order CCLE & GDSC & GDSC-ribo unique.cellid by gCSI unique.cellid
gdsc_circ_counts$cellid <- gdsc$cellid[match(rownames(gdsc_circ_counts), rownames(gdsc))]
gdsc_circ_counts <- gdsc_circ_counts[which(gdsc_circ_counts$cellid %in% c("HCT 116","HeLa")),]

#create combined dataframe of CIRCexplorer counts for each replicate across gCSI, CCLE, GDSC
circ_combined <- data.frame("poly-A"=as.numeric(gdsc_circ_counts$count),
                            "RiboZero"=as.numeric(gdsc_ribozero_circ_counts$count))
rownames(circ_combined) <- gdsc_circ_counts$cellid


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
gdsc_combined$Count <- c(log2(gdsc_combined$Count + 1))
colnames(gdsc_combined) <- c("Replicate", "Method", "Count", "Tool")


png("../results/aa_supplementary_figure5_log2.png", width=145, height=100, units='mm', res = 600, pointsize=80)
ggplot(data = gdsc_combined, aes(x = Replicate, y = Count, fill = Method)) + 
  geom_bar(stat = "identity", size = 0.5, width = 0.5, position = position_dodge(), color = "black") + 
  facet_grid(~factor(Tool, levels=c("CIRI2", "CIRCexplorer2"))) +
  scale_fill_manual(values=c("#8B7B96", "#71A2B6"), labels=c("poly(A)", "RiboZero")) + 
  theme_classic() + 
  scale_y_continuous(limits = c(0, 18), expand=c(0,0)) +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),  
        text = element_text(size = 15), 
        legend.key.size = unit(0.5, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 30), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(x = "Cell line", y = "log2 Normalized \n circRNA Counts")
dev.off()


