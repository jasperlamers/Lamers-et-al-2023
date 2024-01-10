setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## DESEQ2 ANALYSIS FIG 1
library(tximport)
library(SummarizedExperiment)
library(DESeq2)
library(tibble)
library(rtracklayer)
library("gplots")
library("RColorBrewer")
library("AnnotationDbi")
library("org.At.tair.db")
library(ggplot2)
library(readr)

annotations <- read.csv("Data/Annotation_JL.csv", stringsAsFactors = F)

#Read samplelist
setup <- read.csv("Data/Fig1_DESeq2/Samplelist.csv", header = TRUE)

dir_salmon<- "Data/Fig1_DESeq2/Salmon"
setup$filename <- Map(paste, dir_salmon, setup$filename, sep = "")

#SELECT SAMPLES
Genotype <- "Col-0"
Timepoint <- "6hrs"
Experiment <- "96"
Folder_name <- paste(Genotype,Experiment,Timepoint,sep="_")
dir.create(Folder_name, showWarnings = FALSE)

setup <- setup[setup$Experiment == Experiment,]
setup <- setup[setup$Timepoint == Timepoint,]
setup <- setup[setup$Genotype == Genotype,]

#Import all geneIDs for transcriptIDs and remove all n/a
if (!exists("gff")){
  gff_file = file.path("Data/Fig1_DESeq2/Arabidopsis_thaliana.TAIR10.42.gff3")
  gff = import(gff_file)
}

#Full name of file and checks whether it exists
files <- file.path(setup$filename)
names(files) <- paste(setup$Genotype, setup$Treatment, sep="_")
all(file.exists(files))

#Write GeneIDs at TRanscripts levels (sum all splice variants)
tx2gene <- tibble(txid = gff$transcript_id, geneid = as.character(gff$Parent))# ) %>% na.omit()
tx2gene2 <- na.omit(tx2gene)

#Make summary files and import tx (from Salmon)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene2)#, txout=TRUE)
txi_counts = txi[2]
write.csv(txi_counts, file=paste(Folder_name,"/","counts_sum_before_DESeq2.csv", sep=""))

#DESeq2
ddsTxi <- DESeqDataSetFromTximport(txi, setup, design = ~ Treatment)
ddsTxi$Treatment <- relevel(ddsTxi$Treatment, "Control")
dds <- DESeq(ddsTxi)
rownames(dds) <- gsub("gene:","",rownames(dds))

#Dispention graph
jpeg(paste(Folder_name,"/",Folder_name,"_Dispention.jpeg", sep=""),quality=100,width = 500, height = 500, pointsize=10)
plotDispEsts(dds)
dev.off()

#Write output
write.csv(assay(dds), paste(Folder_name,"/","dds.csv", sep=""))

#Write log_output
#rld <- rlog(dds)
rld <- vst(dds)
write.csv(assay(rld), paste(Folder_name,"/","dds_log.csv", sep=""))

#Dividing counts by size_factor
dds_SF <- estimateSizeFactors(dds)
dds_normalized <- counts(dds_SF, normalized=TRUE)
write.csv(dds_normalized, paste(Folder_name,"/","dds_normalized.csv", sep=""))

#Make heatmap
sampleDists <- dist( t( assay(rld) ) )
sampleDistMatrix <- as.matrix( sampleDists )
colours = colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

jpeg(paste(Folder_name,"/",Folder_name,"_Heatmap.jpeg", sep=""),quality=100,width = 1000, height = 1000, pointsize=10)
heatmap.2(sampleDistMatrix, col=colours, key=TRUE, symkey=FALSE, density.info="none",cexRow=1.6,cexCol=1.6,margins=c(15,20),trace="none",srtCol=45)
graphics.off()

pcaData <- plotPCA(rld, intgroup = c("Treatment","Genotype"),ntop=30000, returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

max_no <- ceiling(max(rbind(abs(pcaData$PC1),abs(pcaData$PC2))) * 1.1)

plot <- ggplot(pcaData, aes(PC1, PC2, color=Treatment, shape=Genotype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values=c("#cc79a7", "#e69f00", "#56b4e9", "#009e73", "#f0e442", "#0072b2", "#d55e00")) +
  coord_fixed() + ylim(c(-max_no,max_no)) + xlim(c(-max_no,max_no))

ggsave(paste(Folder_name,"/",Folder_name,"_PCA.jpeg", sep=""), plot = plot, device = NULL, path = NULL,
       scale = 1, width = 10, height = 10, dpi = 500, limitsize = TRUE)

# Pair-wise DEGs on total normalized dds
for(TR in resultsNames(dds)[-1]) {
  #Create folder for results
  TR_subfolder <- paste(Folder_name,TR,sep="/" )
  dir.create(TR_subfolder, showWarnings = FALSE)
  
  #Annotate results file with short name AT and ENTREZID
  res <- results(dds, name=TR)
  res_temp <- merge(as.data.frame(res), annotations, by.x="row.names", by.y="GeneID", all.x=T)
  write.csv(res_temp, paste(TR_subfolder,"/",TR,"_DEGs_DESeq2.csv", sep=""))
  
  #Write summary Results
  sink(paste(TR_subfolder,"/",TR,"_DESeq2_result_summary.txt", sep=""), type=c("output", "message"), append = FALSE)
  summary(results(dds, name=TR, alpha=0.05))
  resSig <- subset(res, padj < 0.05)
  paste("DEG = ", nrow(resSig), " Genes (p<0.05)", sep="")
  sink(NULL)
  
  #Volcano plot
  res_LFC <- lfcShrink(dds, coef=TR, type="normal", res=res)
  
  jpeg(paste(TR_subfolder,"/",TR,"_MA_normal.jpeg", sep=""),quality=100)
  plotMA(res_LFC, ylim=c(-8,8), alpha=0.05, main=TR)
  dev.off()
  
  jpeg(paste(TR_subfolder,"/",TR,"_MA_uncorrected.jpeg", sep=""),quality=100)
  plotMA(res, ylim=c(-8,8), alpha=0.05, main=TR)
  dev.off()
  
  res_LFC <- merge(as.data.frame(res_LFC), annotations, by.x="row.names", by.y="GeneID", all.x=T)
  write.csv(res_LFC, paste(TR_subfolder,"/",TR,"_resLFC.csv", sep=""))
  
  #Histogram
  jpeg(paste(TR_subfolder,"/",TR,"_HIST.jpeg", sep=""),quality=100,width = 500, height = 500, pointsize=10)
  hist(res$pvalue[res$baseMean > 1], breaks=20, col="grey50", border="white")
  dev.off()
}

##VENN DIAGRAMS FIG 1A
library("ggplot2")
library(venn)
Timepoint_OI_total <- c("6Hours","24Hours") #"6hrs" or "24hrs"
#Regulation_OI <- "DOWN" #"ALL", "UP" or "DOWN"
p_value <- 0.05
LFC <- 0.0

for(Timepoint_OI in Timepoint_OI_total){
folder <- paste("Data/Fig1a/", Timepoint_OI, "/", sep="")
DEGs_NaCl <- read.csv(paste(folder,"Treatment_NaCl..130mM._vs_Control_resLFC.csv",sep=""), header = TRUE)
DEGs_NaNO3 <- read.csv(paste(folder,"Treatment_NaNO3..130mM._vs_Control_resLFC.csv",sep=""), header = TRUE)
DEGs_Sorbitol <- read.csv(paste(folder,"Treatment_Sorbitol..260mM._vs_Control_resLFC.csv",sep=""), header = TRUE)
DEGs_KCl <- read.csv(paste(folder,"Treatment_KCl..130mM._vs_Control_resLFC.csv",sep=""), header = TRUE)
DEGs_ABA <- read.csv(paste(folder,"Treatment_ABA..25uM._vs_Control_resLFC.csv",sep=""), header = TRUE)

for(Regulation_OI in c("UP", "DOWN")){
#Upregulation
if (Regulation_OI == "UP"){
Genes_NaCl <- DEGs_NaCl[(DEGs_NaCl$padj < p_value & !is.na(DEGs_NaCl$padj) & DEGs_NaCl$log2FoldChange > LFC),"X"]
Genes_NaNO3 <- DEGs_NaNO3[(DEGs_NaNO3$padj < p_value & !is.na(DEGs_NaNO3$padj) & DEGs_NaNO3$log2FoldChange > LFC),"X"]
Genes_Sorbitol <- DEGs_Sorbitol[(DEGs_Sorbitol$padj < p_value & !is.na(DEGs_Sorbitol$padj) & DEGs_Sorbitol$log2FoldChange > LFC),"X"]
Genes_KCl <- DEGs_KCl[(DEGs_KCl$padj < p_value & !is.na(DEGs_KCl$padj) & DEGs_KCl$log2FoldChange > LFC),"X"]
}

#Downregulation
if (Regulation_OI == "DOWN"){
Genes_NaCl <- DEGs_NaCl[(DEGs_NaCl$padj < p_value & !is.na(DEGs_NaCl$padj) & DEGs_NaCl$log2FoldChange < (-1*LFC)),"X"]
Genes_NaNO3 <- DEGs_NaNO3[(DEGs_NaNO3$padj < p_value & !is.na(DEGs_NaNO3$padj) & DEGs_NaNO3$log2FoldChange < (-1*LFC)),"X"]
Genes_Sorbitol <- DEGs_Sorbitol[(DEGs_Sorbitol$padj < p_value & !is.na(DEGs_Sorbitol$padj) & DEGs_Sorbitol$log2FoldChange < (-1*LFC)),"X"]
Genes_KCl <- DEGs_KCl[(DEGs_KCl$padj < p_value & !is.na(DEGs_KCl$padj) & DEGs_KCl$log2FoldChange < (-1*LFC)),"X"]
}

NaCl <-rep(1:1, length(Genes_NaCl), by=1)
NaNO3 <-rep(1:1, length(Genes_NaNO3), by=1)
Sorbitol <-rep(1:1, length(Genes_Sorbitol), by=1)
KCl <-rep(1:1, length(Genes_KCl), by=1)

DF_NaCl <- data.frame(row.names = Genes_NaCl, "NaCl" = NaCl)
DF_NaNO3 <- data.frame(row.names = Genes_NaNO3, "NaNO3" = NaNO3)
DF_Sorbitol <- data.frame(row.names = Genes_Sorbitol, "Sorbitol" = Sorbitol)
DF_KCl <- data.frame(row.names = Genes_KCl, "KCl" = KCl)

DF_Total <- merge(DF_NaCl, DF_NaNO3, by=0, all=TRUE)
rownames(DF_Total) <- DF_Total$Row.names
DF_Total$Row.names <- NULL
DF_Total <- merge(DF_Total, DF_Sorbitol, by=0, all=TRUE)
rownames(DF_Total) <- DF_Total$Row.names
DF_Total$Row.names <- NULL
DF_Total <- merge(DF_Total, DF_KCl, by=0, all=TRUE)
rownames(DF_Total) <- DF_Total$Row.names
DF_Total$Row.names <- NULL
DF_Total[is.na(DF_Total)] <- 0


#Uses the thresholds of Heatmap (see script below Fig. 1b)
x <- as.data.frame(DF_Total)
rownames(x) <- gsub("gene:","",rownames(x))
if(Timepoint_OI == "6Hours"){HeatmapDF <- read.csv("Data/Fig1a/Clusters/9_ward.D_euclidean_Clusters.csv")}
if(Timepoint_OI == "24Hours"){HeatmapDF <- read.csv("Data/Fig1a/Clusters/6_ward.D_euclidean_Clusters.csv")}
x <- x[rownames(x) %in% HeatmapDF$GeneID,]

file_name = paste(Timepoint_OI,Regulation_OI,"Root.csv",sep="_")
write.csv(x, file_name)
file_name = paste(Timepoint_OI,Regulation_OI,"Root.pdf",sep="_")
ggsave(file_name, venn(x, ilab=TRUE, sncs=1.2, ilcs=1.2, size=c(2000,2000), opacity=0.65, box=FALSE, zcolor=c("#56b4e9","#0072b2","#d55e00","#f0e442")), width=5, height=5, bg = "transparent")
}}

##HEATMAP FIG 1B

#Read all packages required for script
library(ggfortify)
library(pheatmap)
Annotation <- read.csv("Data/Annotation_JL.csv", stringsAsFactors = F, row.names = 1)

#Read all packages required for script
directory <- "Data/Fig1a/6Hours"
filelist_all = list.files(path=directory,pattern="*.csv$",full.names=T)
filelist_all <- filelist_all[!grepl("NaCl..100mM", filelist_all, fixed=TRUE)]

ABA <- read.csv(filelist_all[grepl("ABA", filelist_all, fixed=TRUE)], stringsAsFactors = F)
ABA_LFC <- ABA[c("X","log2FoldChange","padj")]

filelist <- filelist_all[!grepl("ABA", filelist_all, fixed=TRUE)]

#Select DEGs of interest
genes_of_interest <- c()
for(file in filelist){
  temp <- read.csv(file)
  temp <- temp[complete.cases(temp$log2FoldChange), ]
  temp <- temp[complete.cases(temp$padj), ]
  temp <- temp[temp$padj < 0.05,]
  temp <- as.vector(temp$X[abs(temp$log2FoldChange) > 0.5])
  genes_of_interest <- unique(c(genes_of_interest, temp))
}

#PICK ALL GENES THAT ARE 0.5<ABS(LFC)
Heatmap.DF <- data.frame("GeneID"=genes_of_interest)
for(file in filelist){
  temp <- read.csv(file)
  temp <- temp[temp$X %in% genes_of_interest,c('X','log2FoldChange')]
  string <- sub(".*/", "", file,perl = T)
  string <- strsplit(string,"_vs_")[[1]][1]
  string <- sub("Treatment_", "",string)
  string <- sub("\\..", " (",string)
  string <- sub("\\.", ")",string)
  colnames(temp) <- c("GeneID",string)
  Heatmap.DF <- merge(Heatmap.DF, temp, by="GeneID", all.x=TRUE)}
rownames(Heatmap.DF) <- Heatmap.DF$GeneID
Heatmap.DF$GeneID <- NULL

Heatmap.DF_Stat <- data.frame("GeneID"=genes_of_interest)
for(file in filelist_all){
  temp <- read.csv(file)
  temp <- temp[temp$X %in% genes_of_interest,c('X','padj')]
  string <- sub(".*/", "", file,perl = T)
  string <- strsplit(string,"_vs_")[[1]][1]
  string <- sub("Treatment_", "",string)
  string <- sub("\\..", " (",string)
  string <- sub("\\.", ")",string)
  colnames(temp) <- c("GeneID",paste(string,"_padj",sep=""))
  Heatmap.DF_Stat <- merge(Heatmap.DF_Stat, temp, by="GeneID", all.x=TRUE)}

Heatmap.DF_Stat <- Heatmap.DF_Stat[,c(1,4,5,3,6,2)]

#SETTINGS AND COLORS OF HEATMAP
fraction_reduced_black = 0
number_of_colors <- 200
max_cutoff <- 4
min_cutoff <- -4
number_clusters <- 9
cluster_m <- "ward.D"
distance_m <- "euclidean"

colors <- c("#56b4e9", "#FFFFFF", "#d55e00", "#FFFFFF")
ramp_neg <- colorRampPalette(colors[1:2])(number_of_colors/2 * (1+fraction_reduced_black))[1:(number_of_colors/2)]
ramp_pos <- colorRampPalette(colors[3:4])(number_of_colors/2 * (1+fraction_reduced_black))[1:(number_of_colors/2)]
cols <- c(ramp_neg, rev(ramp_pos))

#CREATE FIRST HEATMAP
Heatmap.DF <- Heatmap.DF[,c("Sorbitol (260mM)","NaNO3 (130mM)","NaCl (130mM)","KCl (130mM)")]
heatmap <- pheatmap(Heatmap.DF, border_color = NA, cex = 1, treeheight_row=0, cutree_rows = number_clusters,
                    cellwidth = 16, cellheight = 0.1, show_rownames=F , labels_col = colnames(Heatmap.DF) , cluster_cols = TRUE,
                    cluster_rows = TRUE, clustering_distance_rows = distance_m, color = cols, clustering_method=cluster_m)

#GET CLUSTERS FROM HEATMAP AND ASSIGN ANNOTATION COLUMN (TO COLOR CLUSTERS)
ordered <- rownames(Heatmap.DF[heatmap$tree_row[["order"]],])
clusters <- data.frame(sort(cutree(heatmap$tree_row, k=number_clusters)))
clusterOrder <- unique(clusters[ordered,])
colnames(clusters) <- c("clusters")

#Set clusternames as rownames
Ordered_row <- table(clusters)[clusterOrder]
Ordered_row_median <- Ordered_row/2
Ordered_row_median <- cumsum(Ordered_row) - Ordered_row + Ordered_row_median
Ordered_row_median <- round(Ordered_row_median,0)
row_lab <- as.data.frame(ordered)
row_lab["X"] <- ""
for(index in names(Ordered_row_median)){
  row_lab[as.integer(Ordered_row_median[index]),"X"] <- index
}
rownames(row_lab) <- row_lab$ordered
row_lab <- row_lab[rownames(Heatmap.DF),"X"]

clusters <- merge(clusters, ABA_LFC, by.x="row.names", by.y="X", all.x=T)
rownames(clusters) <- clusters$Row.names
clusters[,c("Row.names","clusters","X","padj")] <- NULL
colnames(clusters) <- "ABA (25µM)"
clusters[clusters$`ABA (25µM)` > max_cutoff,"ABA (25µM)"] <- 4.0
clusters[clusters$`ABA (25µM)` < min_cutoff,"ABA (25µM)"] <- -4.0
clusters["ABA (25µM)"] <- round(clusters["ABA (25µM)"], digits = 1)
clusters["ABA (25µM)"] <- as.character(clusters$`ABA (25µM)`)
cluster_unique <- as.character(sort(as.numeric(unique(clusters[ordered,"ABA (25µM)"]))))

colors <- c()
for(item in cluster_unique){
  item <- as.numeric(item)
  if(item <= min_cutoff){colors <- c(colors, "#56b4e9")
  } else if(item >= max_cutoff){colors <- c(colors, "#d55e00")
  } else {colors <- c(colors, cols[(item - min_cutoff)/((max_cutoff-min_cutoff)/number_of_colors)])}}

names(colors) <- cluster_unique
annotation_color <- list(`ABA (25µM)` = colors)

#SET LEGEND OPTIONS
labels_legend <- as.character((seq(min_cutoff,max_cutoff,by=max_cutoff/2)))
labels_legend[1] <- paste("≤ ",labels_legend[1],sep="")
labels_legend[length(labels_legend)] <- paste("≥ ",labels_legend[length(labels_legend)],sep="")

heatmap <- pheatmap(Heatmap.DF, border_color = NA, treeheight_row=0, treeheight_col=0, cutree_rows = number_clusters,
                    cellwidth = 25, fontsize=16, cellheight = 0.1, show_rownames=T , labels_col = colnames(Heatmap.DF), cluster_cols = TRUE,
                    cluster_rows = TRUE, annotation_row = clusters, clustering_distance_rows = distance_m, color = cols, annotation_legend=F, annotation_colors=annotation_color,
                    breaks=seq(min_cutoff,max_cutoff,by=(max_cutoff-min_cutoff)/number_of_colors), legend_breaks = seq(min_cutoff,max_cutoff,by=max_cutoff/2) ,legend_labels=labels_legend, clustering_method=cluster_m, labels_row=row_lab,
                    gaps_row = head(Ordered_row_median, -1))

ggsave(filename = paste(number_clusters, cluster_m,distance_m,"_RNAseq_heatmap.pdf", sep="_"), plot = heatmap, width = 5, height = 17, dpi=400, limitsize = FALSE, bg = "transparent")
dev.off()

ABA_LFC <- ABA_LFC[,c("X","log2FoldChange")]
colnames(ABA_LFC)[2] <- "ABA (25uM)"
Heatmap.DF <- merge(Heatmap.DF,ABA_LFC,by.x="row.names",by.y="X",all.x=T)
Heatmap.DF <- Heatmap.DF[,c(1,4,3,5,2,6)]

#WRITE CSV WITH CLUSTER INFORMATION
ordered <- data.frame("GeneID"=Heatmap.DF[heatmap$tree_row[["order"]],"Row.names"])
ordered["Order"] <- seq(1:nrow(ordered))

Clusters <- as.data.frame(cutree(as.hclust(heatmap$tree_row), number_clusters))
colnames(Clusters) <- "Cluster"
Clusters <- merge(Clusters, Annotation, by.x=0, by.y="GeneID", all.x=T)
colnames(Clusters)[1] <- "GeneID"
Clusters <- merge(Clusters,ordered, by="GeneID")
Clusters <- Clusters[,c("GeneID","GeneName","Cluster","Order","Na_REG")]

Clusters <- merge(Clusters, Heatmap.DF, by.x="GeneID",by.y="Row.names")
Clusters <- merge(Clusters, Heatmap.DF_Stat, by="GeneID")

SummaryDF <- data.frame()
for(i in clusterOrder){
  ClusterLength <- nrow(Clusters[Clusters$Cluster == i, ])
  SummaryDF <- rbind(SummaryDF, c(i, ClusterLength,SSCount,SS_UP,SS_DOWN))}
colnames(SummaryDF) <- c("Cluster","ClusterSize")
SummaryDF$Name <- c("I","II","III","IV","V","VI","VII","VIII","IX")[1:number_clusters]

SummaryDF$Cluster_order <- row.names(SummaryDF)
Clusters <- merge(Clusters, SummaryDF[,c("Cluster","Cluster_order","Name")], by="Cluster")
write.csv(Clusters,paste(number_clusters, cluster_m,distance_m,"Clusters.csv", sep="_"), row.names = F)

##VENN DIAGRAMS FIG 1C

# Author: Jasper Lamers
# Affiliation: Wageningen University & Research
###### adapted script from:
###### "https://knozue.github.io/post/2018/09/26/over-representation-analysis-1-goseq-with-arabidopsis-go-term.html#fn2" ######

# Script to perform GO enrichment analysis with the GOseq package. With this script it is possible to:
#      - extract genes from enriched GO categories: extract_degsINgo <- function(pwf_out, goseq_out)
#      - get GO annotation of DEGs specific to a treatment of interest: get_treatmentSpecific_GO2DEG <- function(pwf_out, goseq_out)
#      - select specific enriched GO terms and extract the relevant DEGs: select_GOannotations <- function(terms_ofInterest, go2term2deg, deg_data)

# Output is redirected to the folder specified in "all_"
##########################################################################################################################################

library(GO.db)
library(tidyverse)
library("Rgraphviz")
library(org.At.tair.db)
library(goseq)
library(topGO)
library(ggplot2)
library(pheatmap)
###############################################################################################################################

# PREPARE GENE-TO-GO ANNOTATION
# Download TAIR GO annotation file and store it in Atgoslim.TAIR object
sink_dir <- "Sink_update/"
dir.create(sink_dir, showWarnings = FALSE)

empty_DF <- data.frame()

Clusters <- read.csv("Data/Fig1b_output/9_ward.D_euclidean_Clusters.csv")
tairGO <- readr::read_tsv("Data/ATH_GO_GOSLIM.txt", skip = 4, col_names = FALSE)
gene2GO <- unique(tapply(tairGO$X6, tairGO$X1, c)) # construct gene-to-GO annotation array

cdna <- Biostrings::readDNAStringSet("Data/Araport11_genes.201606.cdna.fasta")

bias_cdna <- nchar(cdna)  # get length of cdna sequences
names(bias_cdna) <- mapply(function(cdna_name){
  full_name <- strsplit(cdna_name, " | ")[[1]][1]
  cropped_name <- strsplit(full_name, "\\.")[[1]][1]
}, names(cdna))  # assign the geneID to the length data

median_bias <- mapply(function(unique_name){
  median_length <- median(bias_cdna[which(names(bias_cdna) == unique_name)])
}, unique(names(bias_cdna)))

# NOTE: the same 27586 are tested for all 24hrs samples and the same 25260 for all 6hrs samples
lendata <- median_bias
length(lendata)

for(cluster_no in seq(unique(Clusters$Cluster_order))){
  print(cluster_no)
  
  Clusters_DEG <- Clusters[Clusters$Cluster_order == cluster_no, "GeneID"]
  median_biasDF <- as.data.frame(median_bias)
  median_biasDF["DEG"] <- as.integer(names(median_bias) %in% Clusters_DEG)
  
  ###############################################################################################################################
  # make space
  #rm(cdna, allgenes6hrs_csv, allgenes24hrs_csv, degs6hrs_csv, degs24hrs_csv, allgenes_files, degs6hrs_files, degs24hrs_files)
  ###############################################################################################################################
  # GOSEQ ANLALYSIS
  deg_vector <- median_biasDF$DEG
  names(deg_vector) <- rownames(median_biasDF)
  pwf_cluster <- nullp(DEgenes = deg_vector, bias.data = lendata, plot.fit = FALSE)
  
  y <- rownames(pwf_cluster[pwf_cluster$DEgenes == 1,])
  goseq_out <- goseq(pwf = pwf_cluster, "Arabidopsis", "geneSymbol", gene2cat = gene2GO, use_genes_without_cat = TRUE)
  goseq_out <- goseq_out[!is.na(goseq_out$ontology),]
  goseq_out <- goseq_out[goseq_out$ontology == "BP",]
  
  goseq_out[, 8] <- p.adjust(goseq_out[, 2], method = "BH")
  write.csv(goseq_out, file = paste(sink_dir, cluster_no,"_GOterms_all.csv", sep = ""))
  
  colnames(goseq_out)[8] <- "overep_padjust"
  write.csv(goseq_out, file = paste(sink_dir, cluster_no,"_GOterms_all.csv", sep = ""))
  GOenriched <- goseq_out[which(goseq_out[, 8] < 0.05), ]
  
  if(nrow(GOenriched) != 0){
    GOenriched["Total"] <- nrow(pwf_cluster)
    GOenriched["ClusterSize"] <- length(y)
    
    write.csv(GOenriched[order(GOenriched$ontology),], file = paste(sink_dir, cluster_no,"_GOterms.csv", sep = ""))
    degs2GO <- gene2GO[row.names(pwf_cluster[which(pwf_cluster[, 1] == 1),])]
    GO2deg <- inverseList(degs2GO)[GOenriched[,1]]
    if(length(names(GO2deg)) != 0){for (i in 1:length(names(GO2deg))){
      temp <- as.data.frame(unique(GO2deg[[names(GO2deg)[i]]]))
      colnames(temp) <- "GeneID"
      temp["Term"] <- GOenriched[which(GOenriched[,1] == names(GO2deg)[i]), 1]
      temp["Response"] <- GOenriched[which(GOenriched[,1] == names(GO2deg)[i]), 6]
      temp["Ontology"] <- GOenriched[which(GOenriched[,1] == names(GO2deg)[i]), 7]
      if(i == 1){GO2term2deg <- temp}else{GO2term2deg <- rbind(GO2term2deg, temp)}}
      
      GO2term2deg <- merge(GO2term2deg, Clusters[,c("GeneName","GeneID")], by="GeneID", all.x=T)
      write.csv(GO2term2deg[order(GO2term2deg$Ontology),], file = paste(sink_dir, cluster_no,"_DEGSinGO.csv", sep = ""))
      
      GOenriched$cluster <- cluster_no
      empty_DF <- rbind(empty_DF, GOenriched[order(GOenriched$ontology),])
    }}}

write.csv(empty_DF, file = paste(sink_dir,"All_GOterms.csv", sep = ""))

Summary <- read.csv("Data/Fig1b_output/9_ward.D_euclidean_Clusters.csv")
empty_DF <- merge(empty_DF,Summary[,c("X","Name")],by.x="cluster",by.y="X",all.x=T)
empty_DF <- empty_DF[order(empty_DF$cluster),]
empty_DF$new_order <- seq(nrow(empty_DF),1)

empty_DF["GeneRatio"] <- empty_DF$numDEInCat/empty_DF$numInCat
empty_DF$term <- gsub("double-strand break repair via homologous recombination","double-strand break repair via HR",empty_DF$term)
empty_DF$term <- gsub("negative regulation of DNA-templated transcription, initiation","negative regulation of transcription",empty_DF$term)

GOterms <- unique(empty_DF[empty_DF$numDEInCat >= 10,"term"])
empty_DF <- empty_DF[empty_DF$term %in% GOterms,]
empty_DF$ClusterRatio <- empty_DF$numDEInCat/empty_DF$ClusterSize

DP <- ggplot(empty_DF, aes(x = cluster, y = reorder(term,new_order))) +
  geom_point(aes(size = ClusterRatio, color = overep_padjust)) +
  theme_bw(base_size = 7) +
  scale_colour_gradient(name="p-value",low="#56b4e9", high ="#000000",limits = c(0,0.05)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5),panel.grid.minor = element_blank(), text = element_text(size = 7)) + scale_size(name="Cluster ratio", range=c(0,4), breaks = seq(0.05, 0.15, 0.05), limits = c(0, 0.15)) +
  ylab(NULL) + xlab("Cluster") + guides(color=guide_colourbar(title.vjust=3))+ scale_x_continuous(breaks = round(seq(1,9),1), labels=c("I","II","III","IV","V","VI","VII","VIII","IX"), limits=c(.7,9.3))

ggsave(filename = "All_DotPlot_ClusterR.pdf", plot = DP, dpi=400, width = 3.75, height = 3.8 , limitsize = FALSE)  