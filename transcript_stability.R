load("../results/transcript_stability.RData")
load("../results/Kallisto_isoforms_counts.RData")

###########################
## Transcript GC content ##   
###########################

## Code to obtain GC%, results saved under "/data/transcript_stability/transcript_gc.RData"
    #suppressMessages(library(AnnotationHub))                                 
    #require('ensembldb')
    #edb <- AnnotationHub()[["AH83216"]] #ensembl v101 human hg38
    #transcripts <- transcripts(edb) #GC content

load("/data/transcript_stability/transcript_gc.RData")
gc <- transcripts$gc_content[match(transcript_stability$transcript_id, transcripts$tx_id_version)]
gc_mean <- mean(gc[which(!is.na(gc))])
gc[which(is.na(gc))] <- gc_mean
transcript_stability$gc <- gc



#########################################
## Transcript Length + Number of Exons ##
#########################################

## Code to obtain length and no. exons, results saved under "/data/transcript_stability/txdb.RData"
    #library(GenomicFeatures)
    #library(AnnotationDbi)
    #library(rtracklayer)
    #txdb <- makeTxDbFromGFF("~/Desktop/gencode.v33.annotation.gff3") (ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/gencode.v33.annotation.gff3.gz)
    #lengths <- transcriptLengths(txdb = txdb, with.cds_len = FALSE, with.utr5_len = FALSE)

load("/data/transcript_stability/txdb.RData")                                 
transcript_stability$length <- lengths$tx_len[match(transcript_stability$transcript_id, lengths$tx_name)]
transcript_stability$n_exon <- lengths$nexon[match(transcript_stability$transcript_id, lengths$tx_name)]

mean_transcript_length <- mean(transcript_stability$length[which(!is.na(transcript_stability$length))])
transcript_stability$length[which(is.na(transcript_stability$length))] <- mean_transcript_length

mean_n_exons <- mean(transcript_stability$n_exon[which(!is.na(transcript_stability$n_exon))])
transcript_stability$n_exon[which(is.na(transcript_stability$n_exon))] <- mean_n_exons 


#################
## Mappability ##
#################

# Mappability obtained using UMAP, data obtained from the Hoffman Lab
# All code and source files can be found under "/data/transcript_stability/mappability"
mappability <- read.csv("/data/transcript_stability/mappability/gencode_33_umap.bed", sep = "\t", header = F)
mappability$ID <- paste0(mappability$V1,":",mappability$V2,"-", mappability$V3)                                 
colnames(mappability)[4] <- "single_read_mappability"
colnames(mappability)[5] <- "multi_read_mappability"

meta <- as.data.frame(gCSI_rna@elementMetadata)
meta$start <- meta$start - 1 #change to 0-based to match mappability BED
xx <- paste0(meta$seqnames, ":", meta$start, "-", meta$end)
meta$ID <- xx  


#add gene-ids to transcript_stability
transcript_stability$gene_ID <-  meta$ID[match(transcript_stability$transcript_id, meta$transcript_id)]

#add mappability to transcript_stability
transcript_stability$single_read_mappability <- mappability$single_read_mappability[match(transcript_stability$gene_ID, mappability$ID)]
transcript_stability$multi_read_mappability <- mappability$multi_read_mappability[match(transcript_stability$gene_ID, mappability$ID)]

mean_single_mapp <- mean(transcript_stability$single_read_mappability[which(!is.na(transcript_stability$single_read_mappability))])     

transcript_stability$single_read_mappability[which(is.na(transcript_stability$single_read_mappability))] <- mean_single_mapp

mean_multi_mapp <- mean(transcript_stability$multi_read_mappability[which(!is.na(transcript_stability$multi_read_mappability))])

transcript_stability$multi_read_mappability[which(is.na(transcript_stability$multi_read_mappability))] <- mean_multi_mapp                                 


########################
## 3'/5' Mapping Bias ##
########################

# Bias computed from formula derived from Picard's CollectRnaSeqMetrics (broadinstitute.github.io/picard/)
# All code and source files can be found under "/data/transcript_stability/coverage_bias"
load("/data/transcript_stability/coverage_bias/coverage_bias.RData")
transcript_stability$three_prime_bias <- bias$three_prime_bias[match(transcript_stability$transcript_id,bias$transcript_id)]
transcript_stability$five_prime_bias <- bias$five_prime_bias[match(transcript_stability$transcript_id,bias$transcript_id)]




## save data.frame of stability and features                                                                  
write.csv(transcript_stability, file="../results/transcript_stability.csv")
