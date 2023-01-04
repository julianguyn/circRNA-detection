#load packages
suppressMessages(library(survival))
suppressMessages(library(survminer))
suppressMessages(library(survcomp))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ggpubr))

options(stringsAsFactors = F)

######################################################################
#### Create circRNA log2 normalized count matrix for each circRNA ####
######################################################################

#read in circRNA data
load("../results/summarizeCIRC_lung.RData")
load("../results/summarizeCIRI_lung.RData")


###CIRI2

#read in unique circRNA ids for matrix rows
unique_ids <- unique(lung_ribozero_ciri$circRNA_ids)

summarizeCIRIMatrix <- function(dir_path, unique_circ_ids){
  
  ciri_files <- list.files(path=dir_path,  pattern = "\\.tsv$", full.names = T)
  circid_matrix <- data.frame(matrix(ncol=length(ciri_files), nrow = length(unique_circ_ids)))
  rownames(circid_matrix) <- unique_circ_ids
  
  for (f in 1:length(ciri_files)) {
    sample <- read.table(file = ciri_files[f], 
                         sep = '\t', 
                         skip = 1,
                         header = FALSE)
    
    colnames(sample) <- c("circRNA_ID","chr","circRNA_start","circRNA_end","junction_reads", "SM_MS_SMS", "non_junction_reads", "junction_reads_ratio", "circRNA_type", "gene_id", "strand","junction_reads_ID")
    
    #convert coordinate from 1-based to 0-based to match CIRCexplorer2 coordinates
    sample$circRNA_start <- sample$circRNA_start - 1
    sample$circRNA_ID <- paste0(sample$chr, ":", sample$circRNA_start, "|", sample$circRNA_end)
    
    sample_name <- gsub("\\..*","", ciri_files[f])
    sample_name <- gsub(".*/","", sample_name)
    sample <- sample[which(sample$circRNA_type=="exon"),]
    sample <- sample[which(sample$junction_reads >= 2),]
    
    #group junction reads by circRNA_ID
    circ_reads <- sample %>% 
      group_by(circRNA_ID) %>% 
      summarise(junction_reads = sum(junction_reads))
    
    circ_reads <- as.data.frame(circ_reads)
    rownames(circ_reads) <- circ_reads$circRNA_ID
    
    #add sample to circRNA_ID count matrix
    circid_matrix[rownames(circ_reads), f] <- circ_reads$junction_reads
    names(circid_matrix)[f] <- sample_name
    
  }
  
  #circid_matrix <- circid_matrix[rowSums(is.na(circid_matrix)) != ncol(circid_matrix), ]#remove rows that have NO exp (NA) for any sample
  circid_matrix[is.na(circid_matrix)] <- 0
  return(circid_matrix)
  
}

#run function on CIRI2 data
lung_ribo_ciri_matrix <- summarizeCIRIMatrix(dir_path = "/data/processed_lung/CIRI2/RiboZero/result", unique_circ_ids = unique_ids)

#log2 normalize CIRI2 count matrix
lung_ribo_ciri_matrix <- log2(lung_ribo_ciri_matrix  + 1)



###CIRCexplorer2

#read in unique circRNA ids for matrix rows
unique_ids <- unique(lung_ribozero_circ$circRNA_ids)

summarizeCIRCMatrix <- function(dir_path, unique_circ_ids, RNAseR_dir=FALSE, suff=FALSE){
  
  circ_files <- list.files(path=dir_path,  pattern = ".txt", full.names = T, recursive = TRUE)
  circid_matrix <- data.frame(matrix(ncol=length(circ_files), nrow = length(unique_circ_ids)))
  rownames(circid_matrix) <- unique_circ_ids
  
  for (f in 1:length(circ_files)) {
    sample <- read.table(file = circ_files[f], 
                         sep = '\t', 
                         header = FALSE)
    
    colnames(sample) <- c("chrom","start","end", "name","score","strand","thickStart","thickEnd","itemRgb","exonCount","exonSizes", "exonOffsets", "readNumber", "circType", "geneName", "isoformName", "index", "flankIntron")
    sample_name <- sub(".*/ *(.*?) *.txt.*", "\\1", circ_files[f])
    sample_name <- gsub("_circularRNA_known", "", sample_name)
    sample <- sample[which(sample$circType=="circRNA"),]
    sample <- sample[which(sample$readNumber >= 2),]
    
    #create circRNA position ID, since not given, only in CIRI2 is it given automatically
    sample$circRNA_ID <- paste0(sample$chrom,":",sample$start,"|",sample$end)
    
    #group junction reads by circRNA_ID
    circ_reads <- sample %>% 
      group_by(circRNA_ID) %>% 
      summarise(readNumber = sum(readNumber))
    
    
    circ_reads <- as.data.frame(circ_reads)
    rownames(circ_reads) <- circ_reads$circRNA_ID
    
    #add sample to circRNA_ID count matrix
    circid_matrix[rownames(circ_reads), f] <- circ_reads$readNumber
    names(circid_matrix)[f] <- sample_name
    
    
  }
  
  #circid_matrix <- circid_matrix[rowSums(is.na(circid_matrix)) != ncol(circid_matrix), ]#remove rows that have NO exp (NA) for any sample
  circid_matrix[is.na(circid_matrix)] <- 0
  return(circid_matrix)
}

#run function on CIRCexplorer2 data
lung_ribo_circ_matrix <- summarizeCIRCMatrix(dir_path = "/data/processed_lung/CIRCexplorer2/RiboZero/annotate", unique_circ_ids = unique_ids)

#log2 normalize data
lung_ribo_circ_matrix <- log2(lung_ribo_circ_matrix + 1)


#######################################################################
##### perform survival analysis (Overall) using Concordance Index #####
#######################################################################

#keep only samples where no more than 26 samples have 0 circRNA count for a given circRNA coordinate
lung_ribo_ciri_matrix <- lung_ribo_ciri_matrix[rowSums(lung_ribo_ciri_matrix == 0) <= 26, ]
lung_ribo_circ_matrix <- lung_ribo_circ_matrix[rowSums(lung_ribo_circ_matrix == 0) <= 26, ]

#transpose circRNA expression for integration into survival (CIRI2, CIRCexplorer2)
lung_ribo_ciri_matrix_t <- as.data.frame(t(lung_ribo_ciri_matrix))
lung_ribo_circ_matrix_t <- as.data.frame(t(lung_ribo_circ_matrix))

#read in survival data + circRNA data
survival <- read.csv("/data/survival_data/IMPACT13list_09Jan2019_with_survival_info_55_ADC.csv")
survival$status[which(survival$status == "A")] = 0
survival$status[which(survival$status == "D")] = 1
survival$status <- as.numeric(survival$status)

#survival$rec[which(survival$rec == "N")] = 0
#survival$rec[which(survival$rec == "Y")] = 1
#survival$rec <- as.numeric(survival$rec)

#keep only IMPACT tumor samples (51 samples)
survival <- survival[-c(1,46,47,49),]

#create overall survival df
survival_overall <- survival[,c("surv_time", "status")]
rownames(survival_overall) <- rownames(lung_ribo_ciri_matrix_t)

#CIRI2
ciri_survival <- data.frame(matrix(ncol=5, nrow = ncol(lung_ribo_ciri_matrix_t)))
rownames(ciri_survival) <- colnames(lung_ribo_ciri_matrix_t)
colnames(ciri_survival) <- c("c.index", "se","lower","upper","p.value")

for (i in 1:ncol(lung_ribo_ciri_matrix_t)){
  
  con <- concordance.index(x=lung_ribo_ciri_matrix_t[,i], surv.time=as.numeric(as.character(survival_overall$surv_time)), 
                           surv.event=as.numeric(as.character(survival_overall$status)), 
                           na.rm=TRUE, method="noether")
  
  ciri_survival[i,"c.index"] <- con$c.index
  ciri_survival[i,"se"] <- con$se
  ciri_survival[i,"lower"] <- con$lower
  ciri_survival[i,"upper"] <- con$upper
  ciri_survival[i,"p.value"] <- con$p.value
  
}


#CIRCexplorer2
circ_survival <- data.frame(matrix(ncol=5, nrow = ncol(lung_ribo_circ_matrix_t)))
rownames(circ_survival) <- colnames(lung_ribo_circ_matrix_t)
colnames(circ_survival) <- c("c.index", "se","lower","upper","p.value")

for (i in 1:ncol(lung_ribo_circ_matrix_t)){
  
  con <- concordance.index(x=lung_ribo_circ_matrix_t[,i], surv.time=as.numeric(as.character(survival_overall$surv_time)), 
                           surv.event=as.numeric(as.character(survival_overall$status)), 
                           na.rm=TRUE, method="noether")
  
  circ_survival[i,"c.index"] <- con$c.index
  circ_survival[i,"se"] <- con$se
  circ_survival[i,"lower"] <- con$lower
  circ_survival[i,"upper"] <- con$upper
  circ_survival[i,"p.value"] <- con$p.value
  
}

#multiple testing correction (FDR). Looking for less than 5% or 10%
#False discovery rate (FDR) is designed to control the proportion of false positives (that's the proportion of circRNA's that you call as candidate biomarkers that you expect really aren't candidate biomarkers)
# FDR: If you call 500 circRNA's as candidate biomarkers with an FDR cutoff of 0.1 (10%), you expect 50 of them to be false positives.

ciFDR_ciri <- p.adjust(ciri_survival[,5], method="fdr")
ciFDR_circ <- p.adjust(circ_survival[,5], method="fdr")

ciri_survival$FDR <- ciFDR_ciri
circ_survival$FDR <- ciFDR_circ

#keep circRNA's with (CI > 0.6; p-value < 0.05) OR CI(CI < 0.35 - opposite direction)
ciri_survival_keep_for <- ciri_survival[which(ciri_survival$c.index > 0.65 & ciri_survival$p.value < 0.05 & ciri_survival$FDR < 0.05),]
circ_survival_keep_for <- circ_survival[which(circ_survival$c.index > 0.65 & circ_survival$p.value < 0.05 & circ_survival$FDR < 0.05),]

ciri_survival_keep_rev <- ciri_survival[which(ciri_survival$c.index < 0.35 & ciri_survival$p.value < 0.05 & ciri_survival$FDR < 0.05),]
circ_survival_keep_rev <- circ_survival[which(circ_survival$c.index < 0.35 & circ_survival$p.value < 0.05 & circ_survival$FDR < 0.05),]

#########################################
####### Overall Survival Plot (KM) ######
#########################################

splots <- list()

##CIRI2 

#top circRNA's (CIRI2) with CI < 0.35; p-value < 0.05; FDR < 0.05 that were found with both CIRI2 & CIRCexplorer2. circRNA's with IDS mean that they are experimentally validated from circBase

top_ciri <- ciri_survival_keep_rev[order(ciri_survival_keep_rev$c.index, decreasing = T),]
top_ciri$gene.symbol <- NA
top_ciri$circRNA.ID <- NA

#get circBase ID

circ_base <- read.csv("/data/survival_data/hsa_hg19_circRNA.txt", sep = "\t")
circ_base$ID <- paste0(circ_base$X..chrom,":",circ_base$start,"|",circ_base$end)
x1 <- circ_base[grep("chr20:40161688\\|40179999", circ_base$ID),c("circRNA.ID","gene.symbol")] #only 1 circRNA matched between CIRI2 & CIRCexplorer2

top_ciri$gene.symbol <- x1$gene.symbol
top_ciri$circRNA.ID <- x1$circRNA.ID

lung_ribo_ciri_matrix_t_p <- cbind(lung_ribo_ciri_matrix_t, survival_overall)

colnames(lung_ribo_ciri_matrix_t_p)[grep("chr20:40161688\\|40179999", colnames(lung_ribo_ciri_matrix_t_p))] <- top_ciri$circRNA.ID[1]

res.cut <- surv_cutpoint(lung_ribo_ciri_matrix_t_p, time = "surv_time", event = "status",
                         variables = "hsa_circ_0001159")

res.cat <- surv_categorize(res.cut)
head(res.cat)

fit1 <- survfit(Surv(surv_time, status) ~hsa_circ_0001159, data = res.cat)

splots[[1]] <- ggsurvplot(fit1, palette = c("#71A2B6", "#CA5A7D"), 
                          legend.labs = c("hsa_circ_0001159 = high", "hsa_circ_0001159 = low"),
                          font.main = c(18, "plain", "black"),
                          ggtheme = theme_classic(),
                          risk.table = TRUE, conf.int = FALSE, 
                          title = "CIRI2", xlab="Time (years)")
splots[[1]]$plot <- splots[[1]]$plot + theme(plot.title = element_text(hjust = .5),
                                             legend.key.size = unit(1, 'cm'))



##CIRCexplorer2

#top circRNA's (CIRCexplorer2). circRNA's with IDS mean that they are experimentally validated from circBase

top_circ <- circ_survival_keep_rev[order(circ_survival_keep_rev$c.index, decreasing = T),]
top_circ$gene.symbol <- NA
top_circ$circRNA.ID <- NA

#get circBase ID

x1 <- circ_base[grep("chr20:40161688\\|40179999", circ_base$ID),c("circRNA.ID","gene.symbol")]

top_circ$gene.symbol <- x1$gene.symbol
top_circ$circRNA.ID <- x1$circRNA.ID

lung_ribo_circ_matrix_t_p <- cbind(lung_ribo_circ_matrix_t, survival_overall)

colnames(lung_ribo_circ_matrix_t_p)[grep("chr20:40161688\\|40179999", colnames(lung_ribo_circ_matrix_t_p))] <- top_circ$circRNA.ID[2]

res.cut <- surv_cutpoint(lung_ribo_circ_matrix_t_p, time = "surv_time", event = "status",
                         variables = "hsa_circ_0001159")

res.cat <- surv_categorize(res.cut)
head(res.cat)

fit1 <- survfit(Surv(surv_time, status) ~hsa_circ_0001159, data = res.cat)

png("../results/figure5.png", width=400, height=150, units='mm', res = 600, pointsize=80)
splots[[2]] <- ggsurvplot(fit1, palette = c("#71A2B6", "#CA5A7D"), 
                          legend.labs = c("hsa_circ_0001159 = high", "hsa_circ_0001159 = low"),
                          font.main = c(18, "plain", "black"),
                          ggtheme = theme_classic(),
                          risk.table = TRUE, conf.int = FALSE, 
                          title = "CIRCexplorer2", xlab="Time (years)")
splots[[2]]$plot <- splots[[2]]$plot + theme(plot.title = element_text(hjust = .5),
                                             legend.key.size = unit(1, 'cm'))
arrange_ggsurvplots(splots, print = TRUE,
                    ncol = 2, nrow = 1, risk.table.height = 0.25)
dev.off()


###########################
####### Volcano Plot ######
###########################

# Add gene symbol to survival data frames
ciri_survival$gene <- circ_base[match(rownames(ciri_survival), circ_base$ID),]$gene.symbol
circ_survival$gene <- circ_base[match(rownames(circ_survival), circ_base$ID),]$gene.symbol

# Colour by FDR (<20%, <10%, <5%)
ciri_survival$color <- ">20%"
ciri_survival$color <- ifelse(ciri_survival$FDR < 0.2, "<20%", ciri_survival$color)
ciri_survival$color <- ifelse(ciri_survival$FDR < 0.1, "<10%", ciri_survival$color)
ciri_survival$color <- ifelse(ciri_survival$FDR < 0.05, "<5%", ciri_survival$color)

circ_survival$color <- ">20%"
circ_survival$color <- ifelse(circ_survival$FDR < 0.2, "<20%", circ_survival$color)
circ_survival$color <- ifelse(circ_survival$FDR < 0.1, "<10%", circ_survival$color)
circ_survival$color <- ifelse(circ_survival$FDR < 0.05, "<5%", circ_survival$color)

p1 <- ggplot() +
  geom_point(data = ciri_survival, aes(x = c.index, y = -log10(p.value), color = color), shape = 19, size = 3) +
  scale_x_continuous(limits=c(0.2, 0.8)) +
  scale_color_manual(guide = guide_legend(reverse = FALSE), labels=c("<5%", ">20%"), values = c("#4CC5AB", "#767676")) +
  geom_text(data = ciri_survival, aes(x = c.index, y = -log10(p.value), label = ifelse(color == "<5%", gene, "")), hjust=0.5, vjust=1.5) +
  guides(color = guide_legend(override.aes=list(shape=15, size = 8), title = "FDR")) + 
  geom_vline(xintercept = 0.35, linetype = "dashed", size = 0.5, color = "#767676") +
  geom_vline(xintercept = 0.65, linetype = "dashed", size = 0.5, color = "#767676") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(1, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 17), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRI2", x = "Concordance Index", y = expression(-log[10]*P)) 
p2 <- ggplot() +
  geom_point(data = circ_survival, aes(x = c.index, y = -log10(p.value), color = color), shape = 19, size = 3) +
  scale_x_continuous(limits=c(0.2, 0.8)) +
  scale_color_manual(guide = guide_legend(reverse = FALSE), labels=c("<5%", ">20%"), values = c("#4CC5AB", "#767676")) +
  geom_text(data = circ_survival, aes(x = c.index, y = -log10(p.value), label = ifelse(color == "<5%", gene, "")), hjust=0.5, vjust=1.5) +
  guides(color = guide_legend(override.aes=list(shape=15, size = 8), title = "FDR")) + 
  geom_vline(xintercept = 0.35, linetype = "dashed", size = 0.5, color = "#767676") +
  geom_vline(xintercept = 0.65, linetype = "dashed", size = 0.5, color = "#767676") +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.5),   
        text = element_text(size = 15), 
        legend.key.size = unit(1, 'cm'),
        plot.title = element_text(hjust = 0.5, size = 17), 
        axis.text.x = element_text(size=15, vjust = 0.5), 
        axis.text.y = element_text(size=15)) +
  labs(title = "CIRCexplorer2", x = "Concordance Index", y = expression(-log[10]*P)) 

png("../results/supplementary_figure12.png", width=250, height=125, units='mm', res = 600, pointsize=80)
ggarrange(p1, p2,
          ncol = 2, 
          common.legend = TRUE,
          legend = "right")
dev.off()

