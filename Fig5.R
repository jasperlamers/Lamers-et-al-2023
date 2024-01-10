setwd("")
dirname(rstudioapi::getActiveDocumentContext()$path)
## DESEQ2 ANALYSIS FIG 5
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
setup <- read.csv("Data/0_Input/Samplelist.csv", header = TRUE)

dir_salmon<- "Data/0_Input/Salmon/"
setup$File <- Map(paste, dir_salmon, setup$File, sep = "")

#SELECT SAMPLES
setup <- setup[setup$timepoint == 3,]
setup <- setup[setup$Treatment == "0.5MS" | setup$Treatment == "NaCl",]
setup <- setup[setup$Genotype %in% c("Col0","moca1","snrk2.2.snrk2.3"),]

#Import all geneIDs for transcriptIDs and remove all n/a
if (!exists("gff")){
  gff_file = file.path("Data/Fig1_DESeq2/Arabidopsis_thaliana.TAIR10.42.gff3")
  gff = import(gff_file)
}

#Full name of file and checks whether it exists
files <- file.path(setup$File)
names(files) <- paste(setup$Genotype, setup$Treatment, sep="_")
all(file.exists(files))

#Write GeneIDs at TRanscripts levels (sum all splice variants)
tx2gene <- tibble(txid = gff$transcript_id, geneid = as.character(gff$Parent))# ) %>% na.omit()
tx2gene2 <- na.omit(tx2gene)

#DESeq2
Folder_name <- "Output.I"
dir.create(Folder_name, showWarnings = FALSE)

#Make summary files and import tx (from Salmon)
txi <- tximport(files, type = "salmon", tx2gene = tx2gene2)#, txout=TRUE)
txi_counts = txi[2]
write.csv(txi_counts, file=paste(Folder_name,"/","counts_sum_before_DESeq2.csv", sep=""))

ddsTxi <- DESeqDataSetFromTximport(txi, setup, design = ~ Genotype + Treatment + Genotype:Treatment)#, ignoreRank=T)
ddsTxi$Treatment <- relevel(ddsTxi$Treatment, "0.5MS")
ddsTxi$Genotype <- relevel(ddsTxi$Genotype, "Col0")
dds <- DESeq(ddsTxi)
rownames(dds) <- gsub("gene:","",rownames(dds))

# topGene <- "AT4G34410"
# data <- plotCounts(dds, gene=topGene, intgroup=c("Treatment"), returnData=TRUE)
# data$comb <- rownames(data)
# data[endsWith(data$comb,".1"),"comb"] <- gsub('.{2}$', '', data[endsWith(data$comb,".1"),"comb"])
# data[endsWith(data$comb,".2"),"comb"] <- gsub('.{2}$', '', data[endsWith(data$comb,".2"),"comb"])
# ggplot(data, aes(x=comb, y=count, fill=Treatment)) + ggtitle(topGene) + scale_y_log10() + geom_dotplot(binaxis="y", stackdir="center") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

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

write.csv(pcaData,paste(Folder_name,"/",Folder_name,"_PCA.csv",sep=""))
plot <- ggplot(pcaData, aes(PC1, PC2, color=Treatment, shape=Genotype)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  scale_color_manual(values=c("#cc79a7", "#e69f00", "#56b4e9", "#009e73", "#f0e442", "#0072b2", "#d55e00")) +
  coord_fixed() + ylim(c(-max_no,max_no)) + xlim(c(-max_no,max_no)) + scale_shape_manual(values=c(5,16,17,18,8,7,15))

ggsave(paste(Folder_name,"/",Folder_name,"_PCA.pdf", sep=""), plot = plot, device = NULL, path = NULL,
       scale = 1, width = 10, height = 10, dpi = 500, limitsize = TRUE)

#tests = c(resultsNames(dds)[-1],list(c("Genotype_moca1_vs_Col0","Genotypemoca1.TreatmentNaCl")))

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
  res_LFC <- lfcShrink(dds, coef=TR, type="apeglm", res=res)
  
  jpeg(paste(TR_subfolder,"/",TR,"_MA.jpeg", sep=""),quality=100)
  plotMA(res_LFC, ylim=c(-8,8), alpha=0.05, main=TR)
  dev.off()
  
  res_LFC <- merge(as.data.frame(res_LFC), annotations, by.x="row.names", by.y="GeneID", all.x=T)
  write.csv(res_LFC, paste(TR_subfolder,"/",TR,"_resLFC.csv", sep=""))
  
  #Histogram
  jpeg(paste(TR_subfolder,"/",TR,"_HIST.jpg", sep=""),quality=100,width = 500, height = 500, pointsize=10)
  hist(res$pvalue[res$baseMean > 1], breaks=20, col="grey50", border="white")
  dev.off()
}

#Figure 5 - Selecting overlap snrk2.2/2.3 vs NaCl
library(venn)  
library(ggplot2)

Annotation <- read.csv("Data/Annotation_JL.csv")
Clusters <- read.csv("Data/Fig1b_output/9_ward.D_euclidean_Clusters.csv")

LFC <- 0
pvalue <- 0.05

snrk <- read.csv("Data/Fig5/Genotypesnrk2.2.snrk2.3.TreatmentNaCl_resLFC.csv", stringsAsFactors = F)
TR <- read.csv("Data/Fig5/Treatment_NaCl_vs_0.5MS_resLFC.csv", stringsAsFactors = F)
TR_UP <- TR[(TR$padj < pvalue & !is.na(TR$padj) & TR$log2FoldChange > LFC),c("Row.names","log2FoldChange")]
TR_DOWN <- TR[(TR$padj < pvalue & !is.na(TR$padj) & TR$log2FoldChange < (-1*LFC)),c("Row.names","log2FoldChange")]
TR <- merge(TR_UP,TR_DOWN,by="Row.names",all=T)
rownames(TR) <- TR$Row.names
TR$Row.names <- NULL
TR[is.na(TR)] <- 0
TR[TR != 0] <- 1
colnames(TR) <- c("TR_UP","TR_DOWN")

DF_Total <- TR
for(REG in c("UP","DOWN")){
  temp <- snrk
  temp <- temp[(temp$padj < pvalue & !is.na(temp$padj)),]

  if(REG == "UP"){temp <- temp[(temp$log2FoldChange > LFC),]}
  if(REG == "DOWN"){temp <- temp[(temp$log2FoldChange < (-1*LFC)),]}
  
  temp2 <- rep(1:1, nrow(temp), by=1)
  temp2 <- data.frame(temp2)
  colnames(temp2) <- paste("snrk",REG,sep="")
  rownames(temp2) <- temp$Row.names
  
  DF_Total <- merge(DF_Total, temp2, by=0, all=TRUE)
  rownames(DF_Total) <- DF_Total$Row.names
  DF_Total$Row.names <- NULL
}
  
DF_Total[is.na(DF_Total)] <- 0

ggsave("snrk_3h_venn.png", venn(DF_Total, ilabels=T, zcolor=c("#56b4e9","#0072b2","#d55e00","#f0e442","#e69f00"), opacity=0.65,box=FALSE), width=6.2, height=6.2, unit= "cm", bg = "transparent")
DF_Total <- merge(Annotation,DF_Total,by.y=0,by.x="GeneID",all.y=T)
DF_Total <- merge(DF_Total, snrk,by.x="GeneID",by.y="Row.names",all.x=T,all.y=F)
DF_Total <- merge(DF_Total,Clusters[,c("GeneID","Cluster_order","Name")],by="GeneID",all.x=T)
write.csv(DF_Total,"snrk_3h_venn.csv")

snrkDF <- DF_Total
Enhancedsnrk <- snrkDF[(snrkDF$TR_UP == 1 & snrkDF$snrkUP == 1) | (snrkDF$TR_DOWN == 1 & snrkDF$snrkDOWN == 1),"GeneID"]
Suppressedsnrk <- snrkDF[(snrkDF$TR_DOWN == 1 & snrkDF$snrkUP == 1) | (snrkDF$TR_UP == 1 & snrkDF$snrkDOWN == 1),"GeneID"]

snrk <- read.csv("Data/Fig5/Genotypesnrk2.2.snrk2.3.TreatmentNaCl_resLFC.csv", stringsAsFactors = F)
snrk <- snrk[snrk$Row.names %in% c(Enhancedsnrk,Suppressedsnrk), c("Row.names","baseMean","log2FoldChange")]
colnames(snrk) <- c("GeneID","baseMean","snrk")

NaCl <- read.csv("Data/Fig5/Treatment_NaCl_vs_0.5MS_resLFC.csv", stringsAsFactors = F)
NaCl <- NaCl[NaCl$X %in% snrk$GeneID,c("Row.names","log2FoldChange")]
colnames(NaCl) <- c("GeneID","NaCl")

merged <- merge(snrk,NaCl,by="GeneID")
merged$snrk_TR <- merged$NaCl + merged$snrk
colnames(merged)[3:5] <- c("snrk2223_IF","Col-0_TR","snrk2223_TR")

merged <- merge(merged,Annotation,by="GeneID")
merged <- merge(merged,Clusters[,c("GeneID","Cluster_order","Name")],by="GeneID",all.x=T)

merged["REG"] <- ""
merged[merged$GeneID %in% Enhancedsnrk,"REG"] <- "Enhanced"
merged[merged$GeneID %in% Suppressedsnrk,"REG"] <- "Suppressed"

merged[is.na(merged$Cluster_order),"Cluster_order"] <- 10  
write.csv(merged,"Regulation_SnRKvsNaCl_IF_3H.csv")

##VENN DIAGRAMS FIG 5B
#Read all packages required for script
library(ggfortify)
library(pheatmap)

#Get clusters
Annotation <- read.csv("Data/Annotation_JL.csv", stringsAsFactors = F, row.names = 1)
Clusters <- read.csv("Data/Fig1b_output/9_ward.D_euclidean_Clusters.csv")

#Get enhanced and suppressed genes (LFC > 0 and p<0.05)
DF <- merged
DF <- merge(DF, Clusters[,c("GeneID","Order")],by="GeneID",all.x=T,all.y=F)
DF <- DF[order(DF$Order),]

#SETTINGS AND COLORS OF HEATMAP
fraction_reduced_black = 0
number_of_colors <- 200
max_cutoff <- 6
min_cutoff <- -6
cluster_m <- "ward.D"
distance_m <- "euclidean"

colors <- c("#56b4e9", "#FFFFFF", "#d55e00", "#FFFFFF")
ramp_neg <- colorRampPalette(colors[1:2])(number_of_colors/2 * (1+fraction_reduced_black))[1:(number_of_colors/2)]
ramp_pos <- colorRampPalette(colors[3:4])(number_of_colors/2 * (1+fraction_reduced_black))[1:(number_of_colors/2)]
cols <- c(ramp_neg, rev(ramp_pos))

Heatmap.DF <- DF[,c("GeneID","Col-0_TR","snrk2223_IF","snrk2223_TR","Cluster_order")]
row.names(Heatmap.DF) <- Heatmap.DF$GeneID
Heatmap.DF$GeneID <- NULL

heatmap <- pheatmap(Heatmap.DF[,c("Col-0_TR","snrk2223_TR"),], border_color = NA, cex = 1, treeheight_row=0, cutree_rows = 1,
                    cellwidth = 16, cellheight = 0.1, show_rownames=F , labels_col = colnames(Heatmap.DF) , cluster_cols = TRUE,
                    cluster_rows = TRUE, clustering_distance_rows = distance_m, color = cols, clustering_method=cluster_m)

#GET CLUSTERS FROM HEATMAP AND ASSIGN ANNOTATION COLUMN (TO COLOR CLUSTERS)
ordered <- rownames(Heatmap.DF[heatmap$tree_row[["order"]],])

Heatmap.DF <- Heatmap.DF[order(ordered),]
Heatmap.DF <- Heatmap.DF[order(Heatmap.DF$Cluster_order),]

write.csv(Heatmap.DF,"Heatmap.csv")

#Set clusternames as rownames
labels_legend <- as.character((seq(min_cutoff,max_cutoff,by=max_cutoff/2)))
Ordered_row <- table(Heatmap.DF$Cluster_order)
breaks <- cumsum(Ordered_row) #- Ordered_row + Ordered_row_median

Ordered_row_median <- Ordered_row/2
Ordered_row_median <- cumsum(Ordered_row) - Ordered_row + Ordered_row_median
Ordered_row_median <- round(Ordered_row_median,0)
row_lab <- as.data.frame(rownames(Heatmap.DF))
row_lab["X"] <- ""
for(index in names(Ordered_row_median)){
  row_lab[as.integer(Ordered_row_median[index]),"X"] <- index
}
rownames(row_lab) <- row_lab$`rownames(Heatmap.DF)`
row_lab <- row_lab[rownames(Heatmap.DF),"X"]

AnnotationDF <- DF[,c("GeneID","Cluster_order")]
row.names(AnnotationDF) <- AnnotationDF$GeneID
AnnotationDF$Cluster_order <- as.factor(AnnotationDF$Cluster_order)
AnnotationDF$GeneID <- NULL

heatmap <- pheatmap(Heatmap.DF[,c("Col-0_TR","snrk2223_IF","snrk2223_TR"),], border_color = NA, treeheight_row=0, treeheight_col=0,
                    cellwidth = 25, fontsize=16, cellheight = 0.15, show_rownames=T , labels_col = colnames(Heatmap.DF), cluster_cols = TRUE,
                    cluster_rows = F, annotation_row = AnnotationDF, clustering_distance_rows = distance_m, color = cols, annotation_legend=F,
                    breaks=seq(min_cutoff,max_cutoff,by=(max_cutoff-min_cutoff)/number_of_colors), legend_breaks = seq(min_cutoff,max_cutoff,by=max_cutoff/2) ,legend_labels=labels_legend, clustering_method=cluster_m, labels_row=row_lab,
                    gaps_row = head(breaks, -1))

ggsave(filename = paste(cluster_m,distance_m,"_RNAseq_heatmap.pdf", sep="_"), plot = heatmap, width = 5, height = 9, dpi=400, limitsize = FALSE, bg = "transparent")
