#!/usr/local/bin/Rscript
# R 3.6.2
#---------------------------------------------------------------------------------
# Rat cell culture (soft, stiff 100pa is stiff, 10kpa is soft) RNASeq, soft faster
# synapse density
# Link to publication
# TO ADD ONCE AVAILABLE
#
# Script available from:
# https://github.com/CTR-BFX/Kreysing_Franze
#
#
# Analysis Performed by Xiaohui Zhao
# Centre for Trophoblast Reseach, University of Cambridge, UK
# Copyright Xiaohui Zhao (xz289@cam.ac.uk)
#
#---------------------------------------------------------------------------------
# License Information
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------


message("+-------------------------------------------------------------------------------+")
message("+                            Calling or install libraries                       +")
message("+-------------------------------------------------------------------------------+")

suppressPackageStartupMessages({
  library('DESeq2')
  library('ggplot2')
  library('ggfortify')
  #library('cluster')
  library('RColorBrewer')
  library("cowplot")
  library("pheatmap")
  library("ggrepel")
  library("reshape2")
  library("biomaRt")
  library("matrixStats")
  library("plyr")
  library("BiocParallel")
  library("dplyr")
  library("ggalt")
  library("limma")
  library("apeglm")
  library("gdata") ## read.xls pacakge.
  library("ComplexHeatmap")
  library("dplyr")
  library("methods")
  library("utils")
  library("Seurat")
  library("Matrix")
  library("useful")
  library("SingleCellExperiment")
  library("bigmemory")
  library("mltools")
  library("rhdf5")
  library("recommenderlab")
  library("edgeR")
  library("GetoptLong")
  library("UpSetR")
  library("circlize")
  library("VennDiagram")
  library("eulerr")
  ## GO package
  library("clusterProfiler")
  library("DOSE")
  library("GSEABase")
  library("AnnotationHub")
  library("org.Hs.eg.db")
  library("gage")
  library("gageData")
  library("enrichplot")
  library("ggraph")
  library("ggforce")
  ## secreted and tissueEnrich
  library("TissueEnrich")
  library("UniProt.ws")
  library("readxl")
  library("corrplot") 
  #library("clusterSim") 
  #library("sva")
  library("cqn")
})


register(MulticoreParam(2))


message("+-------------------------------------------------------------------------------+")
message("+            Set up some constants e.g. base directories                        +")
message("+-------------------------------------------------------------------------------+")

Project         <- "CTR_kf284_0004"
significance    <- 0.05
l2fc            <- 1 
elementTextSize <- 10

Base.dir <- "/storage/CTR-Projects/CTR_kf284/CTR_kf284_0004"
setwd(Base.dir)

Count.dir  <- "./count_files"
Out.dir <-"./Analysis/Rm4R_Analysis"


message("+-------------------------------------------------------------------------------+")
message("+ Retrieve ensEMBL annotations                                                  +")
message("+-------------------------------------------------------------------------------+")

## BiocManager::install('grimbough/biomaRt')
listMarts()
ensembl    <- useMart("ensembl")
datasets   <- listDatasets(ensembl)
datasets[grep("Rnor", datasets$version),]
ensembl    <- useDataset("rnorvegicus_gene_ensembl", mart=ensembl)

ensEMBL2id <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'description', 'chromosome_name'), 
                    mart = ensembl,
                    useCache = FALSE)          
head(ensEMBL2id)

ensEMBL2id$description <- gsub("..Source.*", "", ensEMBL2id$description)

message("+-------------------------------------------------------------------------------+")
message("+            Set up the sample table                                            +")
message("+-------------------------------------------------------------------------------+")

sampleFiles      <- grep('*count.txt',list.files(Count.dir),value=TRUE)

sampleNames      <- gsub("_count.txt", "", sampleFiles)

sampleLabels     <- gsub("_[0-9]R", "", sampleNames)

sampleReps       <- unlist(lapply(sampleNames, function(x) strsplit(x, split="_")[[1]][2]))
sampleReps       <- gsub("R", "", sampleReps)

sampleOri        <- sampleLabels
sampleOri        <- gsub("[0-9]_", "_", sampleOri)      

sampleGroup      <- unlist(lapply(sampleNames, function(x) strsplit(x, split="_")[[1]][3]))
sampleGroup      <- ifelse(sampleGroup=="10", "fast", "slow")

sampleCon        <- unlist(lapply(sampleOri, function(x) strsplit(x, split="_")[[1]][1]))

sampleTable      <- data.frame(sampleName=sampleNames, 
                          fileName=sampleFiles,  
                          condition=sampleCon,
                          lane=sampleReps,
                          origin=sampleOri,
                          batch=sampleGroup,
                          label=sampleLabels
)

print(sampleTable)
nrow(sampleTable)
str(sampleTable)

sampleTable.rmR4 <- sampleTable[-grep("4R", sampleTable$sampleName), ] 


message("+-------------------------------------------------------------------------------+")
message("+ Produce ddHTSeq object without collapse the replicates, rm R4                 +")
message("+-------------------------------------------------------------------------------+")



ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable.rmR4, 
                                       directory=Count.dir, 
                                       design=~ lane + label )

dds <- DESeq(ddsHTSeq, parallel=TRUE)


vsd <- vst(dds,     blind=F)
colData(vsd)
assay(vsd) <- limma::removeBatchEffect(assay(vsd), vsd$lane)

message("+-------------------------------------------------------------------------------+")
message("+--------                   Custom PCA                                 ---------+")
message("+-------------------------------------------------------------------------------+")

TOPNUM <- 500
pcaData    <- plotPCA(vsd, ntop=TOPNUM, intgroup=c("label"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf(paste0(Out.dir, "/", Project, "_rm4R_Fig.PCA.T", TOPNUM, ".pdf"),width=10,height=10)
par(bg=NA)
plt1 <- ggplot(pcaData, aes(PC1, PC2, color=label)) +
  geom_mark_ellipse(aes(fill = NULL, color=paste0(label), group=paste0(label), 
                        label=label), alpha=0.1, label.fontsize = 6, label.buffer = unit(15, 'mm')) +
  geom_point(size=2, alpha=0.75) +
  geom_text_repel(aes(label=name), show.legend = FALSE, size=2)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  labs(title = paste0("PCA (Top ", TOPNUM, " Most Variable Genes, all samples)"),
       subtitle = "2 Sequencing Runs for H, 4 lanes each (8 points per sample)",
       caption = "CTR_kf284_0004: Eva Kreysing" ) +
  theme(text = element_text(size=12), 
        plot.subtitle=element_text(size=8,  face="italic", color="darkgrey"),
        legend.position='top', aspect.ratio=1) 

plt2 <- ggplot(pcaData[1:16,], aes(PC1, PC2, color=label)) +
  geom_mark_ellipse(aes(fill = NULL, color=paste0(label), group=paste0(label), 
                        label=label), alpha=0.1, label.fontsize = 6, label.buffer = unit(15, 'mm')) +
  geom_point(size=2, alpha=0.75) +
  geom_text_repel(aes(label=name), show.legend = FALSE, size=2)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  labs(title = paste0("PCA (Top ", TOPNUM, " Most Variable Genes, H samples)"),
       subtitle = "2 Sequencing Runs for H, 4 lanes each (8 points per sample)",
       caption = "CTR_kf284_0004: Eva Kreysing" ) +
  theme(text = element_text(size=12), 
        plot.subtitle=element_text(size=8,  face="italic", color="darkgrey"),
        legend.position='top', aspect.ratio=1) 
plt3 <- ggplot(pcaData[1:8,], aes(PC1, PC2, color=label)) +
  geom_mark_ellipse(aes(fill = NULL, color=paste0(label), group=paste0(label), 
                        label=label), alpha=0.1, label.fontsize = 6, label.buffer = unit(15, 'mm')) +
  geom_point(size=2, alpha=0.75) +
  geom_text_repel(aes(label=name), show.legend = FALSE, size=2)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  labs(title = paste0("PCA (Top ", TOPNUM, " Most Variable Genes, H1 samples)"),
       subtitle = "2 Sequencing Runs for H, 4 lanes each (8 points per sample)",
       caption = "CTR_kf284_0004: Eva Kreysing" ) +
  theme(text = element_text(size=12), 
        plot.subtitle=element_text(size=8,  face="italic", color="darkgrey"),
        legend.position='top', aspect.ratio=1) 
plt4 <- ggplot(pcaData[9:24,], aes(PC1, PC2, color=label)) +
  geom_mark_ellipse(aes(fill = NULL, color=paste0(label), group=paste0(label), 
                        label=label), alpha=0.1, label.fontsize = 6, label.buffer = unit(15, 'mm')) +
  geom_point(size=2, alpha=0.75) +
  geom_text_repel(aes(label=name), show.legend = FALSE, size=2)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  labs(title = paste0("PCA (Top ", TOPNUM, " Most Variable Genes, H2+WT samples)"),
       subtitle = "2 Sequencing Runs for H, 4 lanes each (8 points per sample)",
       caption = "CTR_kf284_0004: Eva Kreysing" ) +
  theme(text = element_text(size=12), 
        plot.subtitle=element_text(size=8,  face="italic", color="darkgrey"),
        legend.position='top', aspect.ratio=1) 
plot_grid(plt1,plt2,plt3,plt4, nrow=2,ncol=2)
dev.off()


message("+---May,2021 add individual genes plot for WT_10 vs WT_100 DGEs and also PCA for these five genes------ +")

genes.WT    <- c("Pcsk9", "Adam3a", "RGD1560775", "Optc", "Ttr", "Thbs1", "Gfap")
elementTextSize   <- 8
makeGeneCountPlot_ord <- function(DDS, ensEMBL2id, CONDITION, gene2plot) {
  #
  # Plot the normalised read counts for a specified gene
  #
  
  genename2plot <- ensEMBL2id[ensEMBL2id$ensembl_gene_id == gene2plot, ]$external_gene_name
  t2            <- plotCounts(DDS, gene=gene2plot, intgroup=c(CONDITION), normalized=TRUE, returnData=TRUE)
  colnames(t2)  <- c("count", "condition")
  t2$count      <- log2(t2$count+1)
  
  t2$condition  <- factor(t2$condition,levels = c("WT_10", "H1_10", "H2_10", "WT_100", "H1_100", "H2_100"))
  
  
  plt.grp <- ggplot(t2, aes(x=condition, y=count, fill=condition)) + 
    geom_boxplot(width = 0.5, color=c("purple","lightgreen", "cyan",  "purple4", "green", "blue"), outlier.shape=NA, alpha=0.5) + 
    geom_point(aes(fill=condition),position=position_jitterdodge(),size=0.8) +
    scale_fill_manual(name="Group", values = c("purple","lightgreen", "cyan",  "purple4", "green", "blue")) +
    theme(text = element_text(size=elementTextSize), legend.position="none") +
    ggtitle(genename2plot) + 
    xlab("") + ylab("log2(Normalised count)") +
    theme_classic() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size=elementTextSize))
  
  
  t2$samples   <- rownames(t2)
  colnames(t2) <- c("count", "group", "samples")
  t2           <- t2[order(t2$group),]
  t2$samples2 <- factor(t2$samples, as.character(t2$samples))
  
  
  plt.ind <- ggplot(t2, aes(x=samples2, y=count, fill=group, group=group)) + 
    geom_bar(stat="identity", alpha=0.5) +
    scale_fill_manual(name="Group", values = c("purple","lightgreen", "cyan",  "purple4", "green", "blue")) +
    xlab("") + ylab("log2(Normalised count)") +
    ggtitle(genename2plot) +
    theme_classic() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  
  print(paste("Created plot for", gene2plot), sep=" ")
  
  return(list(plt.grp, plt.ind))
  
}
makeGeneCountPlot <- function(DDS, ensEMBL2id, CONDITION, gene2plot) {
  #
  # Plot the normalised read counts for a specified gene
  #
  
  genename2plot <- ensEMBL2id[ensEMBL2id$ensembl_gene_id == gene2plot, ]$external_gene_name
  t2            <- plotCounts(DDS, gene=gene2plot, intgroup=c(CONDITION), normalized=TRUE, returnData=TRUE)
  colnames(t2) <- c("count", "condition")
  t2$count <- log2(t2$count+1)
  
  
  plt.grp <- ggplot(t2, aes(x=condition, y=count, fill=condition)) + 
    geom_boxplot(width = 0.5, color=c("purple","purple4", "lightgreen", "green", "cyan", "blue"), outlier.shape=NA, alpha=0.5) + 
    geom_point(aes(fill=condition),position=position_jitterdodge(),size=0.8) +
    scale_fill_manual(name="Group", values = c("purple","purple4", "lightgreen", "green", "cyan", "blue")) +
    theme(text = element_text(size=elementTextSize), legend.position="none") +
    ggtitle(genename2plot) + 
    xlab("") + ylab("log2(Normalised count)") +
    theme_classic() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size=elementTextSize))
  
  
  t2$samples   <- rownames(t2)
  colnames(t2) <- c("count", "group", "samples")
  t2           <- t2[order(t2$group),]
  t2$samples2 <- factor(t2$samples, as.character(t2$samples))
  
  
  plt.ind <- ggplot(t2, aes(x=samples2, y=count, fill=group, group=group)) + 
    geom_bar(stat="identity", alpha=0.5) +
    scale_fill_manual(name="Group", values = c("purple","purple4", "lightgreen", "green", "cyan", "blue")) +
    xlab("") + ylab("log2(Normalised count)") +
    ggtitle(genename2plot) +
    theme_classic() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  
  print(paste("Created plot for", gene2plot), sep=" ")
  
  return(list(plt.grp, plt.ind))
  
}
for(i in 1:length(genes.WT)){
  genes.WT.id     <- ensEMBL2id[ensEMBL2id$external_gene_name == genes.WT[i], ]$ensembl_gene_id
  plt.count.ord   <- makeGeneCountPlot_ord(dds, ensEMBL2id, "label", genes.WT.id)
  plt.count       <- makeGeneCountPlot(dds, ensEMBL2id, "label", genes.WT.id)
  genename2plot   <- ensEMBL2id[ensEMBL2id$ensembl_gene_id == genes.WT.id, ]$external_gene_name
  pdf(paste0(Out.dir, "/", Project, "-", genename2plot, "_Count_plot.pdf"), width=6, height=4 )
  plot_grid(plt.count[[1]], plt.count[[2]], nrow= 2)
  plot_grid(plt.count.ord[[1]], plt.count.ord[[2]], nrow= 2)
  dev.off()
}

## PCA for these genes.
customPCA_selGenes     <- function(sampleTBL, RLD, selGenes, ensEMBL2id, TOPNUM) {
  
  RLDs   <- assay(RLD)
  select <- ensEMBL2id[ensEMBL2id$external_gene_name%in%selGenes, ]$ensembl_gene_id
  pca    <- prcomp(t(RLDs[select, ]))
  
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  
  scores    <- data.frame(sampleName=sampleTBL$sampleName, pca$x, lane=sampleTBL$lane, 
                          condition=sampleTBL$label)
  
  plt.pca <- ggplot(scores, aes(x = PC1, y = PC2, colour=condition, fill=condition, label=sampleName) ) +
    geom_point(size = 3, alpha=0.75) + 
    geom_text_repel(aes(label=sampleName), show.legend = FALSE, size=2) +
    xlab(pc1lab) + ylab(pc2lab) + 
    ggtitle(paste0("PCA Top ", TOPNUM, " MV")) +             
    scale_shape_manual(name="Group", values = c(17, 16)) + 
    theme(text = element_text(size=elementTextSize)) + 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  plt.pca.nl <- ggplot(scores, aes(x = PC1, y = PC2, colour=condition, fill=condition, label=sampleName) ) +
    geom_point(size = 3, alpha=0.75 ) + 
    geom_text_repel(aes(label=sampleName), show.legend = FALSE, size=2) +
    geom_encircle(alpha = 0.1, show.legend = FALSE, aes(colour=condition, fill=condition, group = condition)) +
    xlab(pc1lab) + ylab(pc2lab) + 
    scale_shape_manual(name="Group", values = c(17, 16)) +
    ggtitle(paste0("PCA Top ", TOPNUM, " MV")) +
    theme(text = element_text(size=elementTextSize))+ 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  
  return(list(plt.pca, plt.pca.nl) )
  
}
pdf(paste0(Out.dir, "/", Project, "-selWT_10vs100_7genes_PCA_plot.pdf"), width=8, height=4)
selG.PCA <- customPCA_selGenes(sampleTable.rmR4, vsd, genes.WT, ensEMBL2id, 7)
plot_grid(selG.PCA[[1]], selG.PCA[[2]], ncol = 2)
dev.off()
pdf(paste0(Out.dir, "/", Project, "-selWT_10vs100_5genes_PCA_plot.pdf"), width=8, height=4)
selG.PCA <- customPCA_selGenes(sampleTable.rmR4, vsd, genes.WT[1:5], ensEMBL2id, 5)
plot_grid(selG.PCA[[1]], selG.PCA[[2]], ncol = 2)
dev.off()

## Heatmap for these genes
normS.mat            <- as.data.frame(assay(vsd))
normS.mat$ensembl_gene_id <- rownames(normS.mat)
normS.matE           <- merge(normS.mat, ensEMBL2id, by = "ensembl_gene_id")
normS.matE           <- normS.matE[,c("external_gene_name", "WT_1R_10", "WT_2R_10", "WT_3R_10", 
                                      "H1_1R_10","H1_2R_10", "H1_3R_10", "H2_1R_10","H2_2R_10", "H2_3R_10",
                                      "WT_1R_100", "WT_2R_100", "WT_3R_100",  "H1_1R_100","H1_2R_100", "H1_3R_100",
                                      "H2_1R_100","H2_2R_100", "H2_3R_100")]
normS.matE1          <- normS.matE[normS.matE[,1]%in%genes.WT,]
rownames(normS.matE1)<- normS.matE1[,1]
normS.hmatE1         <- normS.matE1[,-1]

## heatmap [2] "lung development" [3] "respiratory tube development" [6] "respiratory system development"  
breaksList   = seq(6, 12, by = 1)
ScaleCols    <- colorRampPalette(colors = c("purple4","white","darkgreen"))(length(breaksList))
pdf(paste0(Out.dir, "/", Project,"-WT_10vs100_Rm4_sel7Genes_heatmap.pdf"))
pht_listS1 = Heatmap(as.matrix(normS.hmatE1), col = ScaleCols, 
                    name = "WT_10vsWT_100", show_row_names=T,
                    show_column_names = T, width = unit(8, "cm"),
                    heatmap_legend_param = list(title = "NormExp"),
                    cluster_rows = T,show_row_dend = T,
                    column_title="WT_10vsWT_100",
                    cluster_columns = T,
                    #column_km=3,
                    #column_order=levels(as.factor(colnames(normS.hmatE1))),
                    row_title_rot = 0,
                    row_names_gp = gpar( fontsize = 10))
pht_listS2 = Heatmap(as.matrix(normS.hmatE1[-c(1,7),]), col = ScaleCols, 
                     name = "WT_10vsWT_100", show_row_names=T,
                     show_column_names = T, width = unit(8, "cm"),
                     heatmap_legend_param = list(title = "NormExp"),
                     cluster_rows = T,show_row_dend = T,
                     column_title="WT_10vsWT_100",
                     cluster_columns = T,
                     #column_km=4,
                     #column_order=levels(as.factor(colnames(normS.hmatE1[-c(1,7),]))),
                     row_title_rot = 0,
                     row_names_gp = gpar( fontsize = 10))
print(pht_listS1)

print(pht_listS2)
dev.off()
 
 
  
  
  

message("+----- -----------------------------------------------------------------------------------+")
message("+----- Paired comparison (9) and use ComBa_Seq rm replicates effect          -------------+")
message("+----- -----------------------------------------------------------------------------------+")

message("+----- basic defined the paired and general PCA plot function                -------------+")

customPCA_pair     <- function(sampleTBL, RLD, TOPNUM) {
  
  RLDs   <- assay(RLD)
  rv     <- rowVars(RLDs)
  select <- order(rv, decreasing = TRUE)[seq_len(min(TOPNUM, length(rv)))]
  pca    <- prcomp(t(RLDs[select, ]))
  
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  
  scores    <- data.frame(sampleName=sampleTBL$sampleName, pca$x, lane=sampleTBL$lane, 
                          condition=sampleTBL$label)
  
  plt.pca <- ggplot(scores, aes(x = PC1, y = PC2, colour=condition, fill=condition, label=sampleName) ) +
    geom_point(size = 3, alpha=0.75) + 
    geom_text_repel(aes(label=sampleName), show.legend = FALSE, size=2) +
    xlab(pc1lab) + ylab(pc2lab) + 
    ggtitle(paste0("PCA Top ", TOPNUM, " MV")) +             
    scale_shape_manual(name="Group", values = c(17, 16)) + 
    theme(text = element_text(size=elementTextSize)) + 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  plt.pca.nl <- ggplot(scores, aes(x = PC1, y = PC2, colour=condition, fill=condition, label=sampleName) ) +
    geom_point(size = 3, alpha=0.75 ) + 
    geom_text_repel(aes(label=sampleName), show.legend = FALSE, size=2) +
    geom_encircle(alpha = 0.1, show.legend = FALSE, aes(colour=condition, fill=condition, group = condition)) +
    xlab(pc1lab) + ylab(pc2lab) + 
    scale_shape_manual(name="Group", values = c(17, 16)) +
    ggtitle(paste0("PCA Top ", TOPNUM, " MV")) +
    theme(text = element_text(size=elementTextSize))+ 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  
  return(list(plt.pca, plt.pca.nl) )
  
}

ensemble_new       <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'start_position',
                                         'end_position', 'percentage_gene_gc_content'), 
                            mart = ensembl,
                            useCache = FALSE)

control            <- c("H1_100", "H2_100", "H1_10",  "H1_100",  "WT_10", "WT_10", "WT_100", "WT_100", "WT_100")
treatment          <- c("H1_10", "H2_10", "H2_10", "H2_100", "H1_10", "H2_10", "H1_100", "H2_100", "WT_10")
models             <- paste0(treatment, "vs", control)
pairColum          <- list(1:6, 7:12, c(1,3,5,7,9,11), c(2,4,6,8,10,12),
                           c(1,3,5,13,15,17), c(7,9,11,13,15,17),
                           c(2,4,6,14,16,18), c(8,10,12,14,16,18),
                           c(13:18))
sumTab             <- NULL
TOPNUM             <- 500
rmRep              <- "4R"          
 
## 
corTab             <- list()
FDRL.23            <- list()


message("+----- Paired comparison DEGs using edgeR (glm, LR test & PCA analysis       -------------+")

for(i in 1:length(control)){
  print(i)
  sampleTable.Hs   <- subset(sampleTable.rmR4, label==control[i]|label==treatment[i])
  batch            <- sampleTable.Hs$lane
  group            <- factor(sampleTable.Hs$label, levels = c(control[i], treatment[i]))
  count_matrix.Hs  <- assay(ddsHTSeq)[,pairColum[[i]]]
  adj_counts.Hs    <- ComBat_seq(count_matrix.Hs, batch=batch, group=group)
  
  sampleTable.Hs$label <- factor(sampleTable.Hs$label, levels = c(control[i], treatment[i]))
  ddsHTSeq.Hs      <- DESeqDataSetFromMatrix(countData = adj_counts.Hs,
                                             colData = sampleTable.Hs,
                                             design = ~ lane+label)
  ddsHTSeq.Hs      <- estimateSizeFactors(ddsHTSeq.Hs)
  keep             <- rowSums(counts(ddsHTSeq.Hs, normalized=TRUE) >= 5 ) >= 4
  ddsHTSeq.Hs      <- ddsHTSeq.Hs[keep,] 
  ddsHTSeq.Hs      <- DESeq(ddsHTSeq.Hs, parallel=TRUE)
  vsd.Hs           <- vst(ddsHTSeq.Hs, blind=FALSE)
  
  assay(vsd.Hs)    <- limma::removeBatchEffect(assay(vsd.Hs), vsd.Hs$lane)
  
  pdf(paste0(Out.dir, "/", Project, "-PCAplot_", models[i], "_Top", TOPNUM, "_rm", rmRep,".pdf"), width= 8, height = 4)
  pca.plot         <- customPCA_pair(sampleTable.Hs, vsd.Hs, TOPNUM)
  plot_grid(pca.plot[[1]], pca.plot[[2]])
  dev.off()
  
  edgeHTSeq.Hs     <- DGEList(counts = adj_counts.Hs,group  = group, samples = sampleTable.Hs) 
  
  y                <- edgeHTSeq.Hs
  keep             <- filterByExpr(y,  min.count = 5, min.total.count = 15, large.n = 3)
  #keep             <- rowSums(cpm(adj_counts.Hs) >= 5) >=3
  
  table(keep)
  y                <- y[keep, , keep.lib.sizes=FALSE]
  
  AveLogCPM        <- aveLogCPM(y)
  
  ## TMM normalisation
  y                <- calcNormFactors(y)
  y$sample
  cpm              <- cpm(y, log=TRUE, lib.size=libsize, prior.count=2)
  ## add cqn to Correct for GC content and gene length bias
  count.set                 <- as.data.frame(y$counts)
  count.set$ensembl_gene_id <- rownames(count.set)
  count.set.mer             <- merge(count.set, ensemble_new, by = "ensembl_gene_id")
  count.set.mer$length      <- count.set.mer$end_position-count.set.mer$start_position + 1
  
  uCovar                    <- count.set.mer[,c("length", "percentage_gene_gc_content")]
  rownames(uCovar)          <- count.set.mer$ensembl_gene_id
  
  sizeFactors.subset        <- y$samples[,2]
  names(sizeFactors.subset) <- rownames(y$samples)
  
  
  cqn.y                     <- cqn(y$counts, lengths = uCovar$length,x = uCovar$percentage_gene_gc_content, 
                                   sizeFactors = sizeFactors.subset,verbose = TRUE)
  
  design                    <- model.matrix(~batch+batch:group)
  logFC                     <- predFC(y,design,prior.count=1,dispersion=0.05)
  
  ## check correlation 
  corTab[[i]]               <- cor(logFC[,4:6]); 
  
  design                    <- model.matrix(~batch+group)
  #y$offset                  <- cqn.y$glm.offset
  
  norm.counts               <- cqn.y$y + cqn.y$offset
  y                         <- estimateDisp(y, design, robust=TRUE)
  fit                       <- glmQLFit(y, design, robust=TRUE)
  ## QL F-test
  qlf.23           <- glmLRT(fit, coef=2:3) ## batch 2 & 3 sigtest
  
  
  FDR.23           <- p.adjust(qlf.23$table$PValue, method="BH")
  FDRL.23[[i]]     <- sum(FDR.23 < 0.05)
  
  
  qlf              <- glmLRT(fit)
  topTags(qlf)
  top              <- rownames(topTags(qlf))
  cpm(y)[top,]
  
  FDR              <- p.adjust(qlf$table$PValue, method="BH")
  qlf.dat          <- cbind(qlf$table, FDR)
  qlf.dat$ensembl_gene_id <- rownames(qlf.dat)
  qlf.dat.ann      <- merge(qlf.dat, ensEMBL2id, by="ensembl_gene_id")
  qlf.dat.ann      <- qlf.dat.ann[order(qlf.dat.ann$FDR),]
  qlf.ann.sig      <- subset(qlf.dat.ann, abs(logFC) >= 1 & FDR < 0.05)
  numsig           <- dim(qlf.ann.sig)[1]
  
  sumTab           <- c(sumTab, numsig)
  print(numsig)
  
  #if(numsig!=0){
    #write.csv(qlf.ann.sig, file = paste0(Out.dir, "/", Project, "-edgeR_comseq_batch_group_", models[i], "_rm", rmRep, "_l2fc1_p0.05_N", numsig, ".csv"),
    #          row.names=F, quote=T)}
  
  ## save RData for recalling
  save(qlf.dat.ann, norm.counts, file = paste0(Out.dir, "/", Project, "-edgeR_comseq_batch_group_cqnNorm_", models[i], "_rm", rmRep,"_DESeq_Object.RData"))
  
  print(dim(qlf.dat.ann))

  }
sumTabP <- NULL; subTdat <- list(); ssTdat <- list()
for(i in 1:9){
  load(file = paste0(Out.dir, "/", Project, "-edgeR_comseq_batch_group_cqnNorm_", models[i], "_rm", rmRep,"_DESeq_Object.RData"))
  subTdat[[i]] <- subset(qlf.dat.ann, FDR <=0.05 & abs(logFC) >=1)
  write.csv(subTdat[[i]], file = paste0(Out.dir, "/", Project, "-", models[i], "_rm", rmRep, "DGEs_FDR_l2fc_N", dim(subTdat[[i]])[1], "_new.csv"), row.names=F)
  ss <- sum(qlf.dat.ann$FDR <=0.05)
  ssTdat[[i]]  <- subset(qlf.dat.ann, FDR <=0.05)
  write.csv(ssTdat[[i]], file = paste0(Out.dir, "/", Project, "-", models[i], "_rm", rmRep, "DGEs_FDR_N", dim(ssTdat[[i]])[1], "_new.csv"), row.names=F)
  
  sumTabP <- c(sumTabP, ss)
  sumTabP
}


message("+---------------volcano plot------------------------------------------------------+")
functionPlotDEVolcano <- function(results, sig_cut, logfc_cut, title,  xrange, yrange, topN, xlabel, ylabel) {
  
  results       <- as.data.frame(results)
  results$genes <- results$external_gene_name
  results <- subset(results, !is.na(FDR))
  results <- results[order(-results$logFC),]
  
  volc.plt <- ggplot(data=results, aes(x=logFC, y=-log10(FDR), label=genes)) +
    geom_vline(xintercept = logfc_cut,     colour="black", linetype = "dashed", alpha=0.5) +
    geom_vline(xintercept = -(logfc_cut),  colour="black", linetype = "dashed", alpha=0.5) +
    geom_hline(yintercept = -log10(sig_cut), colour="black", linetype = "dashed", alpha=0.5) +
    
    geom_point(data=subset(results, abs(logFC) < logfc_cut | FDR > sig_cut), alpha=0.75, size=0.3, colour="grey") +
    geom_point(data=subset(results, FDR<=sig_cut & logFC >= logfc_cut),      alpha=0.75, size=0.8, colour="red") +
    geom_point(data=subset(results, FDR<=sig_cut & logFC <= -(logfc_cut)),   alpha=0.75, size=0.8, colour="blue") +
    geom_text_repel( data= subset(results, logFC > logfc_cut & FDR<= sig_cut & external_gene_name!="")[1:topN,],
                     show.legend = FALSE, nudge_x=0.1, nudge_y=0.1, segment.size = 0.25, size=3 ) +
    geom_text_repel( data= tail(subset(results, logFC < (-logfc_cut) & FDR<= sig_cut & external_gene_name!=""), n=topN),
                     show.legend = FALSE, nudge_x=0.1, nudge_y=0.1, segment.size = 0.25, size=3 ) +
    xlab(xlabel) + ylab(ylabel) +
    scale_x_continuous(limits=c(xrange[1],xrange[2]), breaks=seq(xrange[1],xrange[2],xrange[3])) +
    scale_y_continuous(limits=c(yrange[1],yrange[2]), breaks=seq(yrange[1],yrange[2],yrange[3])) +
    theme(aspect.ratio=1) +
    ggtitle(title) +
    theme_update(plot.title = element_text(size=16, face="bold", hjust=0.5),
                 axis.title.x = element_text(size=12, face= "bold"),
                 axis.text.x = element_text(size=12, face="bold"),
                 axis.title.y.left = element_text(size=12, face= "bold"),
                 axis.text.y = element_text(size=12, face="bold")) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = "white", colour = NA)) 
  
  
  return(volc.plt)
  
}

sig_cut   <- 0.05
logfc_cut <- 1
topN      <- 50

xlabels   <- c("log2FC (H1_10/H1_100)", "log2FC (H2_10/H2_100)", "log2FC (H2_10/H1_10)","log2FC (H2_100/H1_100)",
               "log2FC (H1_10/WT_10)", "log2FC (H2_10/WT_10)", "log2FC (H1_100/WT_100)","log2FC (H2_100/WT_100)",
               "log2FC (WT_10/WT_100)")
ylabel    <- bquote("-log"[10]~"(adj.p.value)")

xranges   <- list(c(-8,8,2), c(-8,8,2), c(-8,8,2), c(-8,8,2), c(-10,10,2), c(-8,8,2), c(-6,6,2), c(-10,10,1), c(-10,10,1))
yranges   <- list(c(0,10,2), c(0,10,2), c(0,10,2), c(0,10,2), c(0,60,10), c(0,60,10), c(0,40,10), c(0,40,10), c(0,30,10))

for( i in 1:9){
  load(paste0(Out.dir, "/", Project, "-edgeR_comseq_batch_group_cqnNorm_", models[i], "_rm", rmRep, "_DESeq_Object.RData"))
  title      <- models[i]
  xlabel     <- xlabels[i]
  results    <- qlf.dat.ann
  upnum      <- sum(results$logFC>=1 & results$FDR <0.05)
  downnum    <- sum(results$logFC<=-1 & results$FDR <0.05)
  #print(upnum); print(downnum)
  gupnum     <- sum(results$logFC>=0.6 & results$FDR <0.05)
  gdownnum   <- sum(results$logFC<=-0.6 & results$FDR <0.05)
  print(gupnum+gdownnum)
  
  volc.plt   <- functionPlotDEVolcano(results, sig_cut, logfc_cut, title,  xranges[[i]], yranges[[i]], topN, xlabel, ylabel)
  #volc.plt   <- functionPlotDEVolcano(results, sig_cut, logfc_cut, title,  xranges, yranges, topN, xlabel, ylabel)
  pdf(paste0(Out.dir, "/", Project, "-edgeR_", models[i],"_rm", rmRep,"_volcano_plot_new.pdf"))
  print(volc.plt)
  dev.off()
  
}
message("+---------------volcano plot with CV2 correction-----------------------------------------------------+")


functionApplyCV2Filtering  <- function(results.ann, cpm, samT, cont, trt){
  
  t.results                         <- results.ann
  t.counts                          <- as.data.frame(cpm)
  t.counts$ensembl_gene_id          <- rownames(t.counts)
  t.table                           <- merge(t.results, t.counts, by="ensembl_gene_id")
  t.table.undup                     <- t.table[!duplicated(t.table$external_gene_name),]
  rownames(t.table.undup)           <- t.table.undup$external_gene_name
  sel.dim                           <- dim(t.table.undup)[2]
  t.table.undup.filt                <- t.table.undup[, c(1,2,3,6,10:sel.dim)]
  #colnames(t.table.undup.filt)[6:13]<- colData(dds)$condition
  t.table.undup.filt                <- t.table.undup.filt[,-1]
  
  t.table.undup.filt.m              <- melt(t.table.undup.filt, id=c("logFC","logCPM","FDR"))
  colnames(t.table.undup.filt.m)[4] <- "sampleName"
  sample_groups                     <- subset(samT, label==cont|label==trt)
  t.table.undup.filt.m.ann          <- merge(t.table.undup.filt.m, sample_groups, by="sampleName")
  t.table.undup.filt.m.ann          <- t.table.undup.filt.m.ann[, c(1:5,11)]
  t.table.undup.filt.m.ann$external_gene_name <- rep(as.character(rownames(t.table.undup)),dim(sample_groups)[1])
  t.table.undup.filt.m.ann.agg      <- do.call(data.frame, aggregate(t.table.undup.filt.m.ann$value,
                                                                     by = list(t.table.undup.filt.m.ann$external_gene_name,
                                                                               t.table.undup.filt.m.ann$label), 
                                                                     FUN=function(x) c(mean=mean(x), sd=sd(x))) )
  
  colnames(t.table.undup.filt.m.ann.agg) <- c("external_gene_name", "label", "Mean", "SD")
  t.table.undup.filt.m.ann.agg.mrg       <- merge(t.table.undup.filt.m.ann, t.table.undup.filt.m.ann.agg,
                                                  by=c("external_gene_name", "label"))
  
  test.data               <- t.table.undup.filt.m.ann.agg.mrg 
  test.data.u             <- test.data
  test.data.u             <- test.data.u[!duplicated(test.data.u),]
  test.data.u$CV2         <- test.data.u$SD/test.data.u$Mean
  test.data.u             <- test.data.u[-which(test.data.u$external_gene_name==""),]
  return(test.data.u)
}

functionMakeCV2FilteredTable <- function(xx.filt, Project, GROUP1, GROUP2, padj_cut, l2fc_cut){
  
  xx.filt$GROUP1_CV2  <- with(xx.filt, ifelse(label==GROUP1, CV2, NA))
  xx.filt$GROUP2_CV2  <- with(xx.filt, ifelse(label==GROUP2, CV2, NA))
  xx.filt$GROUP1_Mean <- with(xx.filt, ifelse(label==GROUP1, Mean, NA))
  xx.filt$GROUP2_Mean <- with(xx.filt, ifelse(label==GROUP2, Mean, NA))
  xx.filt             <- xx.filt[,c("external_gene_name", "label", "logFC", "logCPM", "FDR", "GROUP1_CV2", "GROUP2_CV2", "GROUP1_Mean", "GROUP2_Mean" )]
  colnames(xx.filt)   <- c("external_gene_name", "label", "logFC", "logCPM", "FDR", paste0(GROUP1,"_CV2"), paste0(GROUP2,"_CV2"), paste0(GROUP1,"_Mean"), paste0(GROUP2,"_Mean") )
  
  xx.1                    <- subset(xx.filt, label==GROUP1)
  xx.2                    <- subset(xx.filt, label==GROUP2)
  xxx                     <- merge(xx.1, xx.2, by=c("external_gene_name", "logFC","logCPM","FDR"))
  xxx                     <- xxx[, !colnames(xxx) %in% c(paste0(GROUP1,"_CV2.y"), paste0(GROUP1,"_Mean.y"), 
                                                         paste0(GROUP2,"_CV2.x"), paste0(GROUP2,"_Mean.x"),
                                                         "label.x", "label.y")]
  xxx.new                 <- unique(xxx)
  print(head(xxx.new))
  
  
  gzf <- gzfile(paste0(Out.dir, "/", Project, "-DESeq2_DEGs_CV2Filtered_padj_", padj_cut, '_l2fc_', l2fc_cut, "_CV2.", GROUP1, "_vs_", GROUP2, '_rmR4.ann.csv.gz'), "w" )
  write.csv(xxx.new[order(abs(xxx.new$logFC),decreasing=TRUE),], gzf)
  close(gzf)
  
  
  return(xxx.new)
}

res.cv2        <- list();
res.cv2FT      <- list();

for(i in 5:8){
   load(paste0(Out.dir, "/", Project, "-edgeR_comseq_batch_group_", models[i], "_rm", rmRep,"_DESeq_Object.RData"))
   results.ann <- qlf.dat.ann
   cpm         <- cpm
   samT        <- sampleTable.rmR4
   cont        <- control[i]
   trt         <- treatment[i]
   res.cv2[[i]]   <- functionApplyCV2Filtering(results.ann, cpm, samT, cont, trt)
   res.cv2FT[[i]] <- functionMakeCV2FilteredTable(res.cv2[[i]], Project, treatment[i], control[i], 0.05, 1)
}

head(res.cv2FT[[5]]); head(res.cv2FT[[6]]); 
head(res.cv2FT[[7]]); head(res.cv2FT[[8]]); 

message("+-------------------------------------------------------------------------------+")
message("+                      Make Volcano and correlation Plots                       +")
message("+-------------------------------------------------------------------------------+")

l2fc_cut         <- 1
padj_cut         <- 0.05
CV2_cut          <- 0.75
topN             <- 10
topN.grey        <- 3
# reads_cut        <- 10
elementTextSize  <- 10
ctrl.col1        <- 7
trt.col1         <- 5
ctrl.colm        <- 8
trt.colm         <- 6


CV2_volcano_plt_function <- function(res.ann.tab, CV2_cut, l2fc_cut, padj_cut, topN, topN.grey, elementTextSize,
                                     ctrl.col1, trt.col1, ctrl.colm, trt.colm, ctrl, trt){
  res.ann.tab                 <- as.data.frame(res.ann.tab)
  res.table.ann.ord.u         <- res.ann.tab[order(res.ann.tab$logFC,decreasing=T),]
  
  res.table.ann.ord.u         <- subset(res.table.ann.ord.u, logFC >= l2fc_cut & FDR <= padj_cut & 
                                        res.table.ann.ord.u[,ctrl.col1]  < CV2_cut   &   
                                        res.table.ann.ord.u[,trt.col1] < CV2_cut)
  
  res.table.ann.ord.d         <- res.ann.tab[order(res.ann.tab$logFC,decreasing=F),]
  res.table.ann.ord.d         <- subset(res.table.ann.ord.d, logFC <= -1*(l2fc_cut) & FDR <= padj_cut & 
                                          res.table.ann.ord.u[,ctrl.col1]  < CV2_cut   &   
                                          res.table.ann.ord.u[,trt.col1] < CV2_cut
                                        )
  
  res.table.ann.ord.u.grey    <- res.ann.tab[order(res.ann.tab$logFC,decreasing=T),]
  res.table.ann.ord.u.grey    <- subset(res.table.ann.ord.u.grey, logFC >= l2fc_cut & FDR <= padj_cut & 
                                          (res.table.ann.ord.u[,ctrl.col1] >= CV2_cut   | res.table.ann.ord.u[,trt.col1]  >= CV2_cut))
  
  res.table.ann.ord.d.grey    <- res.ann.tab[order(res.ann.tab$logFC,decreasing=F),]
  res.table.ann.ord.d.grey    <- subset(res.table.ann.ord.d.grey, logFC <= -1*(l2fc_cut) & FDR <= padj_cut & 
                                          (res.table.ann.ord.u[,ctrl.col1] >= CV2_cut   | res.table.ann.ord.u[,trt.col1]  >= CV2_cut))
  
  plt.volcano.CV2             <- ggplot(data=res.ann.tab, aes(x=logFC, y=-log(FDR), label=external_gene_name)) +
    geom_vline(xintercept = l2fc_cut,       colour="black", linetype = "dashed", alpha=0.5) +
    geom_vline(xintercept = -(l2fc_cut),    colour="black", linetype = "dashed", alpha=0.5) +
    geom_hline(yintercept = -log(padj_cut), colour="black", linetype = "dashed", alpha=0.5) +
    geom_point(data=subset(res.ann.tab, abs(logFC) < l2fc_cut | FDR > padj_cut),
               alpha=0.75, size=1.5, colour="grey85") +
    geom_point(data=subset(res.ann.tab, FDR<=padj_cut & abs(logFC) >= l2fc_cut & 
                             ((res.table.ann.ord.u[,ctrl.col1] >= CV2_cut   | res.table.ann.ord.u[,trt.col1]  >= CV2_cut)) ),
               alpha=0.75, size=1.5, colour="grey40") +
    geom_point(data=subset(res.ann.tab, FDR<=padj_cut & abs(logFC) >= l2fc_cut & 
                             res.table.ann.ord.u[,ctrl.col1]  < CV2_cut   &   
                             res.table.ann.ord.u[,trt.col1] < CV2_cut) ,
               alpha=0.75, size=1.5, colour="firebrick") +
    geom_label_repel(data=res.table.ann.ord.u[c(1:topN),], fill='white', colour='firebrick', point.padding = unit(0.1, "lines"),  
                     size=2,  segment.color = 'firebrick', nudge_x = 0, nudge_y=0) +
    geom_label_repel(data=res.table.ann.ord.d[c(1:topN),], fill='white', colour='firebrick', point.padding = unit(0.1, "lines"),  
                     size=2,  segment.color = 'firebrick', nudge_x = 0, nudge_y=0) +
    
    geom_label_repel(data=res.table.ann.ord.u.grey[c(1:topN.grey),], fill='white', colour='grey40', 
                     point.padding = unit(0.1, "lines"), size=2,  segment.color = 'grey40', nudge_x = 0, nudge_y=0) +
    geom_label_repel(data=res.table.ann.ord.d.grey[c(1:topN.grey),], fill='white', colour='grey40', 
                     point.padding = unit(0.1, "lines"),  size=2,  segment.color = 'grey40', nudge_x = 0, nudge_y=0) +
    
    xlab(bquote("log"[2]~"(Fold Change)")) + 
    ylab(bquote("-log"[10]~"(adj.p.value)")) + 
    ggtitle(paste0(ctrl, "vs", trt, " CV2_", CV2_cut," filtering)")) +
    theme_bw() +
    theme(aspect.ratio=1,
          text=element_text(size=elementTextSize), plot.title=element_text(size=elementTextSize),
          axis.text.x = element_text(size=elementTextSize), axis.text.y = element_text(size=elementTextSize))
  
   return(plt.volcano.CV2)
}

plt.CV2.m5 <- CV2_volcano_plt_function(res.cv2FT[[5]], CV2_cut, l2fc_cut, padj_cut, topN, topN.grey, elementTextSize,
                                       7, 5, 8, 6, control[5], treatment[5])

plt.CV2.m6 <- CV2_volcano_plt_function(res.cv2FT[[6]], CV2_cut, l2fc_cut, padj_cut, topN, topN.grey, elementTextSize,
                                       7, 5, 8, 6, control[6], treatment[6])

plt.CV2.m7 <- CV2_volcano_plt_function(res.cv2FT[[7]], CV2_cut, l2fc_cut, padj_cut, topN, topN.grey, elementTextSize,
                                       7, 5, 8, 6, control[7], treatment[7])

plt.CV2.m8 <- CV2_volcano_plt_function(res.cv2FT[[8]], CV2_cut, l2fc_cut, padj_cut, topN, topN.grey, elementTextSize,
                                       7, 5, 8, 6, control[8], treatment[8])

pdf(paste0(Out.dir, "/", Project, "-CV2_volcano_corr_l2fc_",l2fc_cut,"_CV2_",CV2_cut,"_plot_new.pdf"), width=10, height=10)

plot_grid(plt.CV2.m5, plt.CV2.m6, plt.CV2.m7,plt.CV2.m8, nrow=2)

dev.off()

message("+-------------------------------------------------------------------------------+")
message("+               DEGs summary between H1 vs WT, H2 vs WT                         +")
message("+-------------------------------------------------------------------------------+")
res.ann.list <- list()
resdf.list   <- list()

for(i in 1:4){
  load(paste0(Out.dir, "/", Project, "-edgeR_comseq_batch_group_cqnNorm_", models[i+4], "_rm", rmRep,"_DESeq_Object.RData"))
  res.ann.list[[i]] <- qlf.dat.ann
  resdf.list[[i]]   <-  subset(res.ann.list[[i]], FDR <= significance & abs(logFC)>= l2fc)
}

resdf_10           <- resdf.list[[1]]$external_gene_name[resdf.list[[1]]$external_gene_name%in%resdf.list[[2]]$external_gene_name==T]
resdf_10_logFC1    <- resdf.list[[1]][resdf.list[[1]]$external_gene_name%in%resdf_10==T,c(2,7)]
resdf_10_logFC2    <- resdf.list[[2]][resdf.list[[2]]$external_gene_name%in%resdf_10==T,c(2,7)]
resdf_10.mer       <- merge(resdf_10_logFC1, resdf_10_logFC2, by = "external_gene_name")
colnames(resdf_10.mer) <- c("gene_name", "H1.10.L2FC", "H2.10.L2FC")
resdf_10.mer$class <- ifelse(resdf_10.mer[,2]>0,"up","down")
resdf_10.mer$class <- factor(resdf_10.mer$class, levels = c("up", "down"))



pdf(paste0(Out.dir, "/", Project, "-VennDiagram_H1_H2_WT_10_new.pdf"))
venn.plot  <- draw.pairwise.venn(315,179, 117, c("H1vsWT.10", "H2vsWT.10"), 
                                 col=c("red", "blue"),
                                 fill=c("red", "blue"),
                                 alpha=0.95,
                                 cex=2.5);
grid.draw(venn.plot)
dev.off()

pdf(paste0(Out.dir, "/", Project, "-CommonDEGs_H1_H2_WT_10_new.pdf"))
plt.common <- ggplot(resdf_10.mer, aes(H1.10.L2FC, H2.10.L2FC, colour=class, label=gene_name)) + 
              geom_point() +
              geom_vline(xintercept = l2fc_cut,       colour="black", linetype = "dashed", alpha=0.5) +
              geom_vline(xintercept = -(l2fc_cut),    colour="black", linetype = "dashed", alpha=0.5) +
              geom_hline(yintercept = l2fc_cut,       colour="black", linetype = "dashed", alpha=0.5) +
              geom_hline(yintercept = -(l2fc_cut),    colour="black", linetype = "dashed", alpha=0.5) +
              geom_label_repel(data=resdf_10.mer, fill='white', colour='firebrick', point.padding = unit(0.1, "lines"),  
                               size=2,  segment.color = 'firebrick', nudge_x = 0, nudge_y=0) +
              xlab(bquote("log"[2]~"(Fold Change) H1_10 vs WT_10")) + 
              ylab(bquote("log"[2]~"(Fold Change) H2_10 vs WT_10")) + 
              theme_bw() +
              theme(aspect.ratio=1,
                    text=element_text(size=elementTextSize), plot.title=element_text(size=elementTextSize),
                    axis.text.x = element_text(size=elementTextSize), axis.text.y = element_text(size=elementTextSize))
print(plt.common)
dev.off()

## 100
resdf_100           <- resdf.list[[3]]$external_gene_name[resdf.list[[3]]$external_gene_name%in%resdf.list[[4]]$external_gene_name==T]
resdf_100_logFC1    <- resdf.list[[3]][resdf.list[[3]]$external_gene_name%in%resdf_100==T,c(2,7)]
resdf_100_logFC2    <- resdf.list[[4]][resdf.list[[4]]$external_gene_name%in%resdf_100==T,c(2,7)]
resdf_100.mer       <- merge(resdf_100_logFC1, resdf_100_logFC2, by = "external_gene_name")
colnames(resdf_100.mer) <- c("gene_name", "H1.100.L2FC", "H2.100.L2FC")
resdf_100.mer$class <- ifelse(resdf_100.mer[,2]>0,"up","down")
resdf_100.mer$class <- factor(resdf_100.mer$class, levels = c("up", "down"))

pdf(paste0(Out.dir, "/", Project, "-VennDiagram_H1_H2_WT_100_new.pdf"))
venn.plot  <- draw.pairwise.venn(102,74, 41, c("H1vsWT.100", "H2vsWT.100"), 
                                 col=c("red", "blue"),
                                 fill=c("red", "blue"),
                                 alpha=0.95,
                                 cex=2.5);
grid.draw(venn.plot)
dev.off()

pdf(paste0(Out.dir, "/", Project, "-CommonDEGs_H1_H2_WT_100_new.pdf"))
plt.common <- ggplot(resdf_100.mer, aes(H1.100.L2FC, H2.100.L2FC, colour=class, label=gene_name)) + 
  geom_point() +
  geom_vline(xintercept = l2fc_cut,       colour="black", linetype = "dashed", alpha=0.5) +
  geom_vline(xintercept = -(l2fc_cut),    colour="black", linetype = "dashed", alpha=0.5) +
  geom_hline(yintercept = l2fc_cut,       colour="black", linetype = "dashed", alpha=0.5) +
  geom_hline(yintercept = -(l2fc_cut),    colour="black", linetype = "dashed", alpha=0.5) +
  geom_label_repel(data=resdf_100.mer, fill='white', colour='firebrick', point.padding = unit(0.1, "lines"),  
                   size=2,  segment.color = 'firebrick', nudge_x = 0, nudge_y=0) +
  xlab(bquote("log"[2]~"(Fold Change) H1_100 vs WT_100")) + 
  ylab(bquote("log"[2]~"(Fold Change) H2_100 vs WT_100")) + 
  theme_bw() +
  theme(aspect.ratio=1,
        text=element_text(size=elementTextSize), plot.title=element_text(size=elementTextSize),
        axis.text.x = element_text(size=elementTextSize), axis.text.y = element_text(size=elementTextSize))
print(plt.common)
dev.off()


message("+-------------------------------------------------------------------------------+")
message("+               Gene Ontology Analysis with ClusterProfiler                     +")
message("+-------------------------------------------------------------------------------+")


ensEMBL2id_GO <- getBM(attributes=c('ensembl_gene_id',  'entrezgene_id'), 
                    mart = ensembl,
                    useCache = FALSE)          
head(ensEMBL2id_GO)

#ensEMBL2id_GO$description <- gsub("..Source.*", "", ensEMBL2id_GO$description)

res.ann.list <- list()
resdf.list   <- list()
significance <- 0.05
l2fc_go      <- 1

for(i in 1:4){
  load(paste0(Out.dir, "/", Project, "-edgeR_comseq_batch_group_cqnNorm_", models[i+4], "_rm", rmRep,"_DESeq_Object.RData"))
  res.ann.list[[i]] <- qlf.dat.ann
  resdf.list[[i]]   <-  subset(res.ann.list[[i]], FDR <= significance & abs(logFC)>= l2fc_go)
  resdf.list[[i]]   <-  merge(resdf.list[[i]], ensEMBL2id_GO, by = "ensembl_gene_id")
  resdf.list[[i]]   <- resdf.list[[i]][,c("entrezgene_id", "ensembl_gene_id", "external_gene_name", "logFC")]
  colnames(resdf.list[[i]]) <- c("ENTREZID", "ENSEMBL", "SYMBOL", "L2FC") 
  resdf.list[[i]] <- resdf.list[[i]][order(-resdf.list[[i]]$L2FC),]
  print(dim(resdf.list[[i]]))
  print(length(unique(resdf.list[[i]]$SYMBOL))) ;  
  print(length(unique(resdf.list[[i]]$ENTREZID))) 
}

message("+---------------Run GO analysis: BP, CC, MF + KEGG -----------------------------+")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Rn.eg.db")

library("org.Rn.eg.db")
download_KEGG("rno", keggType = "KEGG", keyType = "kegg")

GO.BP   <- list()
GO.MF   <- list()
GO.CC   <- list()
GO.KEGG <- list()


for(a in 1:4){
  GO.BP[[a]] <- enrichGO(gene          = resdf.list[[a]]$SYMBOL,
                         OrgDb         = org.Rn.eg.db,
                         keyType       = 'SYMBOL',
                         ont           = "BP",
                         pAdjustMethod = "BH",
                         minGSSize     = 10,
                         pvalueCutoff  = 0.05 )
  GO.CC[[a]] <- enrichGO(gene          = resdf.list[[a]]$SYMBOL,
                         OrgDb         = org.Rn.eg.db,
                         keyType       = 'SYMBOL',
                         ont           = "CC",
                         pAdjustMethod = "BH",
                         minGSSize     = 10,
                         pvalueCutoff  = 0.05 )
  GO.MF[[a]] <- enrichGO(gene          = resdf.list[[a]]$SYMBOL,
                         OrgDb         = org.Rn.eg.db,
                         keyType       = 'SYMBOL',
                         ont           = "MF",
                         pAdjustMethod = "BH",
                         minGSSize     = 10,
                         pvalueCutoff  = 0.05 )
  GO.KEGG[[a]] <-  enrichKEGG(resdf.list[[a]]$ENTREZID, organism = "rno")  
  GO.KEGG[[a]] <- setReadable(GO.KEGG[[a]], OrgDb = org.Rn.eg.db, keyType = "ENTREZID")
  
}

save(GO.BP, GO.CC, GO.MF, GO.KEGG,file=paste0(Out.dir, "/", Project, "-KEGG_GO_p0.05_l2fc1_Analysis_new.RData"))


message("+---------------Run GO analysis simply for BP -----------------------------+")
GO.BP2       <- list()

for(a in 1:4){
  GO.BP2[[a]] <- simplify(GO.BP[[a]], cutoff=0.7, by="p.adjust", select_fun=min)
}

write.csv(as.data.frame(GO.BP2[[1]]), file = paste0(Out.dir,"/", Project, "-H1_10vsWT_10_GO_BP_N", dim(GO.BP2[[1]])[1], "_new.csv"), row.names=F)
write.csv(as.data.frame(GO.BP2[[2]]), file = paste0(Out.dir,"/", Project, "-H2_10vsWT_10_GO_BP_N", dim(GO.BP2[[2]])[1], "_new.csv"), row.names=F)
write.csv(as.data.frame(GO.BP2[[3]]), file = paste0(Out.dir,"/", Project, "-H1_100vsWT_100_GO_BP_N", dim(GO.BP2[[3]])[1], "_new.csv"), row.names=F)
write.csv(as.data.frame(GO.BP2[[4]]), file = paste0(Out.dir,"/", Project, "-H2_100vsWT_100_GO_BP_N", dim(GO.BP2[[4]])[1], "_new.csv"), row.names=F)

write.csv(as.data.frame(GO.CC[[1]]), file = paste0(Out.dir,"/", Project, "-H1_10vsWT_10_GO_CC_N", dim(GO.CC[[1]])[1], "_new.csv"), row.names=F)
write.csv(as.data.frame(GO.CC[[2]]), file = paste0(Out.dir,"/", Project, "-H2_10vsWT_10_GO_CC_N", dim(GO.CC[[2]])[1], "_new.csv"), row.names=F)
write.csv(as.data.frame(GO.CC[[3]]), file = paste0(Out.dir,"/", Project, "-H1_100vsWT_100_GO_CC_N", dim(GO.CC[[3]])[1], "_new.csv"), row.names=F)
write.csv(as.data.frame(GO.CC[[4]]), file = paste0(Out.dir,"/", Project, "-H2_100vsWT_100_GO_CC_N", dim(GO.CC[[4]])[1], "_new.csv"), row.names=F)

write.csv(as.data.frame(GO.MF[[2]]), file = paste0(Out.dir,"/", Project, "-H2_10vsWT_10_GO_MF_N", dim(GO.MF[[2]])[1], "_new.csv"), row.names=F)
write.csv(as.data.frame(GO.MF[[3]]), file = paste0(Out.dir,"/", Project, "-H1_100vsWT_100_GO_MF_N", dim(GO.MF[[3]])[1], "_new.csv"), row.names=F)
write.csv(as.data.frame(GO.MF[[4]]), file = paste0(Out.dir,"/", Project, "-H2_100vsWT_100_GO_MF_N", dim(GO.MF[[4]])[1], "_new.csv"), row.names=F)

write.csv(as.data.frame(GO.KEGG[[1]]), file = paste0(Out.dir,"/", Project, "-H1_10vsWT_10_KEGG_N",dim(GO.KEGG[[1]])[1], "_new.csv"), row.names=F)
write.csv(as.data.frame(GO.KEGG[[2]]), file = paste0(Out.dir,"/", Project, "-H2_10vsWT_10_KEGG_N",dim(GO.KEGG[[2]])[1], "_new.csv"), row.names=F)
write.csv(as.data.frame(GO.KEGG[[3]]), file = paste0(Out.dir,"/", Project, "-H1_100vsWT_100_KEGG_N",dim(GO.KEGG[[3]])[1],"_new.csv"), row.names=F)
write.csv(as.data.frame(GO.KEGG[[4]]), file = paste0(Out.dir,"/", Project, "-H2_100vsWT_100_KEGG_N",dim(GO.KEGG[[4]])[1],"_new.csv"), row.names=F)




message("+---------------Run GO analysis up/down regulated compare -----------------------------+")

GO.Comp.BP   <- list()
GO.Comp.CC   <- list()
GO.Comp.MF   <- list()
KEGG.Comp    <- list()
for(a in c(1:4)){
  print(a)
  resdfC         <- resdf.list[[a]]
  resdfC         <- data.frame(Entrez=resdfC$ENTREZID, FC=resdfC$L2FC)
  resdfC$group   <- "upregulated"
  resdfC$group[resdfC$FC < 0] <- "downregulated"
  print(table(resdfC$group))
  KEGG.Comp[[a]]    <- compareCluster(Entrez~group, data=resdfC, fun="enrichKEGG", organism="rno")
  GO.Comp.BP[[a]]   <- compareCluster(Entrez~group, data=resdfC, fun="enrichGO", ont='BP', OrgDb='org.Rn.eg.db')
  GO.Comp.BP[[a]]   <- simplify(GO.Comp.BP[[a]], cutoff=0.7, by="p.adjust", select_fun=min)
  GO.Comp.CC[[a]]   <- compareCluster(Entrez~group, data=resdfC, fun="enrichGO", ont='CC', OrgDb='org.Rn.eg.db')
 # if(a!=2){
  GO.Comp.MF[[a]]   <- compareCluster(Entrez~group, data=resdfC, fun="enrichGO", ont='MF', OrgDb='org.Rn.eg.db')
  #}
}

save(GO.Comp.BP,GO.Comp.CC,GO.Comp.MF, KEGG.Comp, file = paste0(Out.dir,"/", Project, "-GeneOntology_Compare_updown_BP_CC_MF_KEGG.RData"))
write.csv(as.data.frame(GO.Comp.BP[[1]]), file = paste0(Out.dir,"/", Project, "-H1_10vsWT_10_GO_BP_updown_Compare_N", dim(GO.Comp.BP[[1]])[1], "_new.csv"), row.names=F)
write.csv(as.data.frame(GO.Comp.BP[[2]]), file = paste0(Out.dir,"/", Project, "-H2_10vsWT_10_GO_BP_updown_Compare_N", dim(GO.Comp.BP[[2]])[1], "_new.csv"), row.names=F)
write.csv(as.data.frame(GO.Comp.BP[[3]]), file = paste0(Out.dir,"/", Project, "-H1_100vsWT_100_GO_BP_updown_Compare_N", dim(GO.Comp.BP[[3]])[1], "_new.csv"), row.names=F)
write.csv(as.data.frame(GO.Comp.BP[[4]]), file = paste0(Out.dir,"/", Project, "-H2_100vsWT_100_GO_BP_updown_Compare_N", dim(GO.Comp.BP[[4]])[1], "_new.csv"), row.names=F)

write.csv(as.data.frame(GO.Comp.CC[[1]]), file = paste0(Out.dir,"/", Project, "-H1_10vsWT_10_GO_CC_updown_Compare_N", dim(GO.Comp.CC[[1]])[1], "_new.csv"), row.names=F)
write.csv(as.data.frame(GO.Comp.CC[[2]]), file = paste0(Out.dir,"/", Project, "-H2_10vsWT_10_GO_CC_updown_Compare_N", dim(GO.Comp.CC[[2]])[1], "_new.csv"), row.names=F)
write.csv(as.data.frame(GO.Comp.CC[[3]]), file = paste0(Out.dir,"/", Project, "-H1_100vsWT_100_GO_CC_updown_Compare_N", dim(GO.Comp.CC[[3]])[1], "_new.csv"), row.names=F)
write.csv(as.data.frame(GO.Comp.CC[[4]]), file = paste0(Out.dir,"/", Project, "-H2_100vsWT_100_GO_CC_updown_Compare_N", dim(GO.Comp.CC[[4]])[1], "_new.csv"), row.names=F)

write.csv(as.data.frame(GO.Comp.MF[[1]]), file = paste0(Out.dir,"/", Project, "-H1_10vsWT_10_GO_MF_updown_Compare_N", dim(GO.Comp.MF[[1]])[1], "_new.csv"), row.names=F)
write.csv(as.data.frame(GO.Comp.MF[[2]]), file = paste0(Out.dir,"/", Project, "-H2_10vsWT_10_GO_MF_updown_Compare_N", dim(GO.Comp.MF[[2]])[1], "_new.csv"), row.names=F)
write.csv(as.data.frame(GO.Comp.MF[[3]]), file = paste0(Out.dir,"/", Project, "-H1_100vsWT_100_GO_MF_updown_Compare_N", dim(GO.Comp.MF[[3]])[1], "_new.csv"), row.names=F)
write.csv(as.data.frame(GO.Comp.MF[[4]]), file = paste0(Out.dir,"/", Project, "-H2_100vsWT_100_GO_MF_updown_Compare_N", dim(GO.Comp.MF[[4]])[1], "_new.csv"), row.names=F)

write.csv(as.data.frame(KEGG.Comp[[1]]), file = paste0(Out.dir,"/", Project, "-H1_10vsWT_10_KEGG_updown_Compare_N", dim(KEGG.Comp[[1]])[1], "_new.csv"), row.names=F)
write.csv(as.data.frame(KEGG.Comp[[2]]), file = paste0(Out.dir,"/", Project, "-H2_10vsWT_10_KEGG_updown_Compare_N", dim(KEGG.Comp[[2]])[1], "_new.csv"), row.names=F)
write.csv(as.data.frame(KEGG.Comp[[3]]), file = paste0(Out.dir,"/", Project, "-H1_100vsWT_100_KEGG_updown_Compare_N", dim(KEGG.Comp[[3]])[1], "_new.csv"), row.names=F)
write.csv(as.data.frame(KEGG.Comp[[4]]), file = paste0(Out.dir,"/", Project, "-H2_100vsWT_100_KEGG_updown_Compare_N", dim(KEGG.Comp[[4]])[1], "_new.csv"), row.names=F)


message("+---------------Generate bar/dot plot for selected BP, MF, CC -----------------------------+")

for(a in 1:4){
  i <- a+4
  print(a)
  pdf(paste0(Out.dir, "/", Project, "-", models[i], "_UpDown_Compare_GOKEGG_dotplot_new.pdf"), width=15, height=10)
  dotplot(KEGG.Comp[[a]], showCategory= 20, includeAll=F, title = "KEGG")
  dev.off()
  
  pdf(paste0(Out.dir, "/", Project,  "-", models[i], "_UpDown_Compare_GOBP_dotplot_new.pdf"), width=15, height=10)
  dotplot(GO.Comp.BP[[a]],  showCategory = 20, includeAll=F, title="Biological Process")
  dev.off()
  
  #if(a!=1){
  pdf(paste0(Out.dir, "/", Project, "-", models[i], "_UpDown_Compare_GOMF_dotplot_new.pdf"), width=15, height=10)
  dotplot(GO.Comp.MF[[a]],  showCategory = 20, includeAll=T, title ="Molecular Function")
  dev.off()
  #}
  pdf(paste0(Out.dir, "/", Project, "-", models[i], "_UpDown_Compare_GOCC_dotplot_new.pdf"), width=15, height=10)
  dotplot(GO.Comp.CC[[a]],  showCategory = 20, includeAll=F, title ="Cellular Component")
  dev.off()
  
  
  ##
  pdf(paste0(Out.dir, "/", Project, "-", models[i], "_GOKEGG_dotplot_new.pdf"), width=8, height=10)
  dotplot(GO.KEGG[[a]], showCategory= 20, title = "KEGG")
  dev.off()
  
  pdf(paste0(Out.dir, "/", Project,  "-", models[i], "_GOBP_dotplot_new.pdf"), width=8, height=10)
  dotplot(GO.BP[[a]],  showCategory = 20, title="Biological Process")
  dev.off()
  
  #if(a!=1){
  pdf(paste0(Out.dir, "/", Project, "-", models[i], "_GOMF_dotplot_new.pdf"), width=8, height=10)
  dotplot(GO.MF[[a]],  showCategory = 20, title ="Molecular Function")
  dev.off() #}
  
  pdf(paste0(Out.dir, "/", Project, "-", models[i], "_GOCC_dotplot_new.pdf"), width=8, height=10)
  dotplot(GO.CC[[a]],  showCategory = 20, title ="Cellular Component")
  dev.off()
  
  
}

message("+----- GeneOntology Key words to search the common genes and pathways ------------+")

keywords <- c("channel",
              "membrane",
              "resting potential",
              "synapsogenesis",
              "synaptic",
              "neuro",
              "maturation",
              "mechanosensitivity", "mechanical stimuli",
              "ECM", "cytoskeleton",
              "adhesion",
              "development",
              "growth",
              "neuron", "brain",
              "microtubules",
              "actin",
              "integrin",
              "ncam", 
              "fate", "stem cell", "precursor cell",
              "CNS","FGF", "BMP", "Notch", "WNT",
              "glucoseaminoglykans",
              "inflammation", "wound ",
              "epithelia",
              "BMP")
load(paste0(Out.dir, "/", Project, "-KEGG_GO_p0.05_l2fc1_Analysis_new.RData"))
selBP <- list(); selKEGG <- list(); selCC <- list(); selMF <- list()
for(i in 1:4){
  selBP.sub      <- lapply(keywords, function(x) GO.BP[[i]][grep(x, GO.BP[[i]]$Description),])
  testlen        <- unlist(lapply(selBP.sub, function(x) dim(x)[1]))
  testlen.rmind  <- which(testlen==0)
  selBP.new      <- selBP.sub[-testlen.rmind]
  selBP[[i]]     <- selBP.new
  
  selKEGG.sub      <- lapply(keywords, function(x) GO.KEGG[[i]][grep(x, GO.KEGG[[i]]$Description),])
  testlen          <- unlist(lapply(selKEGG.sub, function(x) dim(x)[1]))
  testlen.rmind    <- which(testlen==0)
  selKEGG.new      <- selKEGG.sub[-testlen.rmind]
  selKEGG[[i]]       <- selKEGG.new
  
  selCC.sub      <- lapply(keywords, function(x) GO.CC[[i]][grep(x, GO.CC[[i]]$Description),])
  testlen        <- unlist(lapply(selCC.sub, function(x) dim(x)[1]))
  testlen.rmind  <- which(testlen==0)
  selCC.new      <- selCC.sub[-testlen.rmind]
  selCC[[i]]     <- selCC.new
  
  selMF.sub      <- lapply(keywords, function(x) GO.MF[[i]][grep(x, GO.MF[[i]]$Description),])
  testlen        <- unlist(lapply(selMF.sub, function(x) dim(x)[1]))
  testlen.rmind  <- which(testlen==0)
  selMF.new      <- selMF.sub[-testlen.rmind]
  selMF[[i]]     <- selMF.new
  }

selBP.keys1       <- unique(purrr::invoke(rbind, selBP[[1]])) ## 192
selBP.keys2       <- unique(purrr::invoke(rbind, selBP[[2]])) ## 138
dim(selBP.keys1);dim(selBP.keys2)
sum(selBP.keys1$ID%in%selBP.keys2$ID==T)
## common 104
selBP.keys3       <- unique(purrr::invoke(rbind, selBP[[3]])) ## 99
selBP.keys4       <- unique(purrr::invoke(rbind, selBP[[4]])) ## 61
dim(selBP.keys3);dim(selBP.keys4)
sum(selBP.keys3$ID%in%selBP.keys4$ID==T)
s1234 <- sum(selBP.keys1$ID[selBP.keys1$ID%in%selBP.keys2$ID==T]%in%selBP.keys3$ID[selBP.keys3$ID%in%selBP.keys4$ID==T]==T)
## common 34
## all common 26

selKEGG.keys1       <- unique(purrr::invoke(rbind, selKEGG[[1]])) ## 3, ECM, Notch,Focal
selKEGG.keys2       <- unique(purrr::invoke(rbind, selKEGG[[2]])) ## 3, ECM, Notch,Focal
## common 3
selKEGG.keys3       <- unique(purrr::invoke(rbind, selKEGG[[3]])) ## 1, ECM
## all common ECM


selCC.keys1       <- unique(purrr::invoke(rbind, selCC[[1]])) ## 4, basement membrane,basolateral plasma membrane,synaptic cleft,apical plasma membrane
selCC.keys2       <- unique(purrr::invoke(rbind, selCC[[2]])) ## 6, membrane raft, membrane microdomain,membrane region,basolateral plasma membrane,synaptic cleft,apical plasma membrane
dim(selCC.keys1);dim(selCC.keys2)
sum(selCC.keys1$ID%in%selCC.keys2$ID==T)
## common 3

selCC.keys3       <- unique(purrr::invoke(rbind, selCC[[3]])) ## 3, external side of plasma membrane,basement membrane,basolateral plasma membrane
selCC.keys4       <- unique(purrr::invoke(rbind, selCC[[4]])) ## 1, external side of plasma membrane
dim(selCC.keys3);dim(selCC.keys4)
sum(selCC.keys3$ID%in%selCC.keys4$ID==T)
## common 1
s1234 <- sum(selCC.keys1$ID[selCC.keys1$ID%in%selCC.keys2$ID==T]%in%selCC.keys3$ID[selCC.keys3$ID%in%selCC.keys4$ID==T]==T)
## all common None.

selMF.keys2       <- unique(purrr::invoke(rbind, selMF[[2]])) ## 1,oxidoreductase activity, acting on the CH-NH2 group of donors
selMF.keys3       <- unique(purrr::invoke(rbind, selMF[[3]])) ## 1,amide transmembrane transporter activity
selMF.keys4       <- unique(purrr::invoke(rbind, selMF[[4]])) ## 1,eoxidoreductase activity, acting on the CH-NH2 group of donors
## all common 0

## calling the pathways which selected by Eva's interests
int.Pathways     <- read.csv(paste0(Out.dir, "/descriptions_reduced.csv"), header=F)
colnames(int.Pathways) <- c("Description", "Rank")
int.BP <- list(); int.KEGG <- list(); int.CC <- list(); int.MF <- list()
for(i in 1:4){
 int.BP[[i]]   <- GO.BP[[i]][GO.BP[[i]]$Description%in%int.Pathways[,1],] ## 141,108,87,35, 1&2 common 85, 3&4 common 24, all common 15
 int.KEGG[[i]] <- GO.KEGG[[i]][GO.KEGG[[i]]$Description%in%int.Pathways[,1],] ## 4,5,1,0. 1&2 common 3,(ECM,Notch and Focal) all common ECM
 int.CC[[i]]   <- GO.CC[[i]][GO.CC[[i]]$Description%in%int.Pathways[,1],] ## 8,10,6,3, 1&2 common 6, 3&4 common 3, all common 2 (extracellular matrix & collagen trimer)
 int.MF[[i]]   <- GO.MF[[i]][GO.MF[[i]]$Description%in%int.Pathways[,1],] ## 0,0,0,0
}

int.BP.dim     <- lapply(int.BP, dim); int.CC.dim       <- lapply(int.CC, dim); 
int.MF.dim     <- lapply(int.MF, dim); int.KEGG.dim     <- lapply(int.KEGG, dim); 

sel.BP_H1_WT_10  <- unique(rbind(selBP.keys1, int.BP[[1]]))
sel.BP_H2_WT_10  <- unique(rbind(selBP.keys2, int.BP[[2]]))
sel.BP_H1_WT_100 <- unique(rbind(selBP.keys3, int.BP[[3]]))
sel.BP_H2_WT_100 <- unique(rbind(selBP.keys4, int.BP[[4]]))

sel.CC_H1_WT_10  <- unique(rbind(selCC.keys1, int.CC[[1]]))
sel.CC_H2_WT_10  <- unique(rbind(selCC.keys2, int.CC[[2]]))
sel.CC_H1_WT_100 <- unique(rbind(selCC.keys3, int.CC[[3]]))
sel.CC_H2_WT_100 <- unique(rbind(selCC.keys4, int.CC[[4]]))

sel.KEGG_H1_WT_10  <- unique(rbind(selKEGG.keys1, int.KEGG[[1]]))
sel.KEGG_H2_WT_10  <- unique(rbind(selKEGG.keys2, int.KEGG[[2]]))
sel.KEGG_H1_WT_100 <- unique(rbind(selKEGG.keys3, int.KEGG[[3]]))

sel.MF_H2_WT_10    <- unique(selMF.keys2)
sel.MF_H1_WT_100   <- unique(selMF.keys3)
sel.MF_H2_WT_100   <- unique(selMF.keys4)

sel.BP_H1_WT_10.mer  <- merge(sel.BP_H1_WT_10, int.Pathways, by = "Description", all.x=T) ## 1-72;2-26;3-43;NA-142; total 283
sel.BP_H2_WT_10.mer  <- merge(sel.BP_H2_WT_10, int.Pathways, by = "Description", all.x=T) ## 1-54;2-23;3-31;NA-99; total 207
sel.BP_H1_WT_100.mer <- merge(sel.BP_H1_WT_100, int.Pathways, by = "Description", all.x=T) ## 1-43;2-14;3-30;NA-71; total 158
sel.BP_H2_WT_100.mer <- merge(sel.BP_H2_WT_100, int.Pathways, by = "Description", all.x=T) ## 1-10;2-8;3-17;NA-54; total 89

sel.CC_H1_WT_10.mer  <- merge(sel.CC_H1_WT_10, int.Pathways, by = "Description", all.x=T) ## 1-6;2-2;total 8
sel.CC_H2_WT_10.mer  <- merge(sel.CC_H2_WT_10, int.Pathways, by = "Description", all.x=T) ## 1-6;2-4;total 10
sel.CC_H1_WT_100.mer <- merge(sel.CC_H1_WT_100, int.Pathways, by = "Description", all.x=T) ## 1-5;2-1;total 6
sel.CC_H2_WT_100.mer <- merge(sel.CC_H2_WT_100, int.Pathways, by = "Description", all.x=T) ## 1-2;2-1;total 3

sel.KEGG_H1_WT_10.mer  <- merge(sel.KEGG_H1_WT_10, int.Pathways, by = "Description", all.x=T) ## 1-3;2-1; total 4
sel.KEGG_H2_WT_10.mer  <- merge(sel.KEGG_H2_WT_10, int.Pathways, by = "Description", all.x=T) ## 1-2;2-2;3-1; total 5
sel.KEGG_H1_WT_100.mer <- merge(sel.KEGG_H1_WT_100, int.Pathways, by = "Description", all.x=T) ## 1-1; total 1

sel.MF_H2_WT_10.mer  <- merge(sel.MF_H2_WT_10, int.Pathways, by = "Description", all.x=T) ## NA 10; total 10
sel.MF_H1_WT_100.mer <- merge(sel.MF_H1_WT_100, int.Pathways, by = "Description", all.x=T) ## NA 1; total 1
sel.MF_H2_WT_100.mer <- merge(sel.MF_H2_WT_100, int.Pathways, by = "Description", all.x=T) ## NA 1; total 1


testNoRank <- unique(c(sel.BP_H1_WT_10.mer[is.na(sel.BP_H1_WT_10.mer$Rank),]$Description,
                sel.BP_H2_WT_10.mer[is.na(sel.BP_H2_WT_10.mer$Rank),]$Description,
                sel.BP_H1_WT_100.mer[is.na(sel.BP_H1_WT_100.mer$Rank),]$Description,
                sel.BP_H2_WT_100.mer[is.na(sel.BP_H2_WT_100.mer$Rank),]$Description))
write.csv(testNoRank, file = paste0(Out.dir, "/SelPathway_NoRank_list_for_Eva_new.csv"), row.names=F)

## Generate the pathways with genes list for Eva
sel.BP_H1_WT_10.merN  <- sel.BP_H1_WT_10.mer[order(sel.BP_H1_WT_10.mer$Rank),c("Description", "ID", "geneID", "Rank")]
sel.BP_H2_WT_10.merN  <- sel.BP_H2_WT_10.mer[order(sel.BP_H2_WT_10.mer$Rank),c("Description", "ID", "geneID", "Rank")]
sel.BP_H1_WT_100.merN <- sel.BP_H1_WT_100.mer[order(sel.BP_H1_WT_100.mer$Rank),c("Description", "ID", "geneID", "Rank")]
sel.BP_H2_WT_100.merN <- sel.BP_H2_WT_100.mer[order(sel.BP_H2_WT_100.mer$Rank),c("Description", "ID", "geneID", "Rank")]

sel.CC_H1_WT_10.merN  <- sel.CC_H1_WT_10.mer[order(sel.CC_H1_WT_10.mer$Rank),c("Description", "ID", "geneID", "Rank")]
sel.CC_H2_WT_10.merN  <- sel.CC_H2_WT_10.mer[order(sel.CC_H2_WT_10.mer$Rank),c("Description", "ID", "geneID", "Rank")]
sel.CC_H1_WT_100.merN <- sel.CC_H1_WT_100.mer[order(sel.CC_H1_WT_100.mer$Rank),c("Description", "ID", "geneID", "Rank")]
sel.CC_H2_WT_100.merN <- sel.CC_H2_WT_100.mer[order(sel.CC_H2_WT_100.mer$Rank),c("Description", "ID", "geneID", "Rank")]

sel.MF_H2_WT_10.merN  <- sel.MF_H2_WT_10.mer[order(sel.MF_H2_WT_10.mer$Rank),c("Description", "ID", "geneID", "Rank")]
sel.MF_H1_WT_100.merN <- sel.MF_H1_WT_100.mer[order(sel.MF_H1_WT_100.mer$Rank),c("Description", "ID", "geneID", "Rank")]
sel.MF_H2_WT_100.merN <- sel.MF_H2_WT_100.mer[order(sel.MF_H2_WT_100.mer$Rank),c("Description", "ID", "geneID", "Rank")]

sel.KEGG_H1_WT_10.merN  <- sel.KEGG_H1_WT_10.mer[order(sel.KEGG_H1_WT_10.mer$Rank),c("Description", "ID", "geneID", "Rank")]
sel.KEGG_H1_WT_100.merN <- sel.KEGG_H1_WT_100.mer[order(sel.KEGG_H1_WT_100.mer$Rank),c("Description", "ID", "geneID", "Rank")]
sel.KEGG_H2_WT_10.merN  <- sel.KEGG_H2_WT_10.mer[order(sel.KEGG_H2_WT_10.mer$Rank),c("Description", "ID", "geneID", "Rank")]

sel.H1_WT_10.merN    <- rbind(sel.BP_H1_WT_10.merN, sel.CC_H1_WT_10.merN,  sel.KEGG_H1_WT_10.merN)
sel.H2_WT_10.merN    <- rbind(sel.BP_H2_WT_10.merN, sel.CC_H2_WT_10.merN, sel.MF_H2_WT_10.merN, sel.KEGG_H2_WT_10.merN)
sel.H1_WT_100.merN   <- rbind(sel.BP_H1_WT_100.merN, sel.CC_H1_WT_100.merN, sel.MF_H1_WT_100.merN, sel.KEGG_H1_WT_100.merN)
sel.H2_WT_100.merN   <- rbind(sel.BP_H2_WT_100.merN, sel.CC_H2_WT_100.merN, sel.MF_H2_WT_100.merN)

sel.H1_WT_10.merN$GO <- rep(c("BP", "CC", "KEGG"), c(dim(sel.BP_H1_WT_10.merN)[1], dim(sel.CC_H1_WT_10.merN)[1],dim(sel.KEGG_H1_WT_10.merN)[1]))
sel.H2_WT_10.merN$GO  <- rep(c("BP", "CC", "MF", "KEGG"), c(dim(sel.BP_H2_WT_10.merN)[1], dim(sel.CC_H2_WT_10.merN)[1], dim(sel.MF_H2_WT_10.merN)[1], dim(sel.KEGG_H2_WT_10.merN)[1]))
sel.H1_WT_100.merN$GO <- rep(c("BP", "CC", "MF", "KEGG"), c(dim(sel.BP_H1_WT_100.merN)[1], dim(sel.CC_H1_WT_100.merN)[1], dim(sel.MF_H1_WT_100.merN)[1], dim(sel.KEGG_H1_WT_100.merN)[1]))
sel.H2_WT_100.merN$GO <- rep(c("BP", "CC", "MF"), c(dim(sel.BP_H2_WT_100.merN)[1], dim(sel.CC_H2_WT_100.merN)[1], dim(sel.MF_H2_WT_100.merN)[1]))

write.csv(sel.H1_WT_10.merN, file = paste0(Out.dir, "/SelPathway_H1_10vsWT_10_BP283_CC8_KEGG4_N295_list_for_Eva_new.csv"), row.names=F)
write.csv(sel.H2_WT_10.merN, file = paste0(Out.dir, "/SelPathway_H2_10vsWT_10_BP207_CC10_MF10_KEGG5_N232_list_for_Eva_new.csv"), row.names=F)
write.csv(sel.H1_WT_100.merN, file = paste0(Out.dir, "/SelPathway_H1_100vsWT_100_BP158_CC6_MF1_KEGG1_N165_list_for_Eva_new.csv"), row.names=F)
write.csv(sel.H2_WT_100.merN, file = paste0(Out.dir, "/SelPathway_H2_100vsWT_100_BP89_CC3_MF1_N93_list_for_Eva_new.csv"), row.names=F)

save(sel.BP_H1_WT_10.mer, sel.BP_H2_WT_10.mer,sel.BP_H1_WT_100.mer, sel.BP_H2_WT_100.mer,
     sel.CC_H1_WT_10.mer, sel.CC_H2_WT_10.mer,sel.CC_H1_WT_100.mer, sel.CC_H2_WT_100.mer,
     sel.MF_H2_WT_10.mer, sel.MF_H1_WT_100.mer, sel.MF_H2_WT_100.mer,
     sel.KEGG_H1_WT_10.mer, sel.KEGG_H1_WT_100.mer,sel.KEGG_H2_WT_100.mer,
     file= paste0(Out.dir, "/SelPathways_Compare4_BP_CC_MF_KEGG_new.RData"))


message("+---- Add information for the direction of the pathways up/down regulated based on genes direction------------------+")

selPathways.list <- list(sel.H1_WT_10.merN, sel.H2_WT_10.merN, sel.H1_WT_100.merN, sel.H2_WT_100.merN)
Direct_Genes_fn  <- function(data, gene.col, refd ){
  GS            <- data[,gene.col]
  GSall         <- lapply(GS, function(x) strsplit(x, split="/")[[1]])
  GSall.res     <- lapply(GSall, function(x) refd[refd$external_gene_name%in%x==T,2])
  GSall.res.U   <- lapply(GSall.res, function(x) ifelse(x>0,1,0))
  UpGenes       <- list()
  for(ind in 1:length(GSall)){
    UpGenes[[ind]] <- paste(GSall[[ind]][which(GSall.res.U[[ind]]==1)],collapse="/")
    UpGenes
  }
  DwGenes       <- list()
  for(ind in 1:length(GSall)){
    DwGenes[[ind]] <- paste(GSall[[ind]][which(GSall.res.U[[ind]]==0)],collapse="/")
    DwGenes
  }
  UpGenesRatio  <- unlist(lapply(GSall.res.U, function(x) sum(x)/length(x)))
  UpGenes       <- unlist(UpGenes)
  DwGenes       <- unlist(DwGenes)
  PathwayDirec  <- ifelse(UpGenesRatio>0.5, "Up", "Down")
  PathwayDirec  <- ifelse(UpGenesRatio==0.5, "Unknow", PathwayDirec)
  dataN         <- cbind(data, UpGenes, DwGenes, UpGenesRatio, PathwayDirec) 
  dataN
}
selPathways.Nlist <- list()
for(i in 1:4){
  selPathways.Nlist[[i]] <- Direct_Genes_fn(selPathways.list[[i]], 3, res.ann.list[[i]])
}

save(selPathways.Nlist, file = paste0(Out.dir, "/SelPathways_Compare4_BP_CC_MF_KEGG_Up_Dw_PDir_new.RData"))

for(i in 1:4){
  write.csv(selPathways.Nlist[[i]], file = paste0(Out.dir, "/SelPathways_", models[i+4], "_BP_CC_MF_KEGG_Up_Dw_PDir_N", dim(selPathways.Nlist[[i]])[1], "_new.csv"), row.names=F)
}


message("+---May,2021 add individual genes plot for WT_10 vs WT_100 DGEs and also PCA for these five genes------ +")

## Generate normcounts matrix after filter and QC control using edgeR, why the pairwise did not filter it out???
edgeHTSeq.Hs     <- DGEList(counts = rawmat, group  = group, samples = sampleTable.rmR4) 

y                <- edgeHTSeq.Hs

AveLogCPM        <- aveLogCPM(y)

## TMM normalisation
y                <- calcNormFactors(y)
y$sample
cpm              <- cpm(y, log=TRUE, lib.size=libsize, prior.count=2)
## add cqn to Correct for GC content and gene length bias
count.set                 <- as.data.frame(y$counts)
count.set$ensembl_gene_id <- rownames(count.set)
count.set.mer             <- merge(count.set, ensemble_new, by = "ensembl_gene_id")
count.set.mer$length      <- count.set.mer$end_position-count.set.mer$start_position + 1

uCovar                    <- count.set.mer[,c("length", "percentage_gene_gc_content")]
rownames(uCovar)          <- count.set.mer$ensembl_gene_id

sizeFactors.subset        <- y$samples[,2]
names(sizeFactors.subset) <- rownames(y$samples)


cqn.y                     <- cqn(y$counts, lengths = uCovar$length,x = uCovar$percentage_gene_gc_content, 
                                 sizeFactors = sizeFactors.subset,verbose = TRUE)


norm.counts               <- cqn.y$y + cqn.y$offset
write.csv(norm.counts, file = paste0(Out.dir, "/", Project, "-edgeR_batch_combat_cqn_normalisedCounts_rmR4.csv"))

## individual gene counts plot 
genes.WT.new    <- c("Ripk3", "Fcgr2b", "Optc", "Ttr", "Thbs1", "Gfap")
elementTextSize   <- 8

## check dds.counts and edgeR norm.counts.

for(i in 1:length(genes.WT.new)){
  genes.WT.id     <- ensEMBL2id[ensEMBL2id$external_gene_name == genes.WT.new[i], ]$ensembl_gene_id
  print(i);
  print(genes.WT.id)
  plt.count.ord   <- makeGeneCountPlot_ord(dds, ensEMBL2id, "label", genes.WT.id)
  plt.count       <- makeGeneCountPlot(dds, ensEMBL2id, "label", genes.WT.id)
  genename2plot   <- ensEMBL2id[ensEMBL2id$ensembl_gene_id == genes.WT.id, ]$external_gene_name
  pdf(paste0(Out.dir, "/", Project, "-", genename2plot, "_Count_plot_new.pdf"), width=6, height=4 )
  plot_grid(plt.count[[1]], plt.count[[2]], nrow= 2)
  plot_grid(plt.count.ord[[1]], plt.count.ord[[2]], nrow= 2)
  dev.off()
}


## Check the sigGenes in WT_10 vs WT_100 in the other four comparisons H1/2 vs WT
subT <- list()
for(i in 1:8){
  load(paste0(Out.dir, "/", Project, "-edgeR_comseq_batch_group_cqnNorm_", models[i], "_rm", rmRep,"_DESeq_Object.RData"))
  subT[[i]]   <- qlf.dat.ann[qlf.dat.ann$external_gene_name%in%genes.WT.new==T,]
  subT
}


## PCA for these genes.

pdf(paste0(Out.dir, "/", Project, "-selWT_10vs100_6genes_PCA_plot_new.pdf"), width=8, height=4)
selG6.PCA <- customPCA_selGenes(sampleTable.rmR4, vsd, genes.WT.new, ensEMBL2id, 6)
plot_grid(selG6.PCA[[1]], selG6.PCA[[2]], ncol = 2)
dev.off()
pdf(paste0(Out.dir, "/", Project, "-selWT_10vs100_4genes_PCA_plot_new.pdf"), width=8, height=4)
selG4.PCA <- customPCA_selGenes(sampleTable.rmR4, vsd, genes.WT.new[1:4], ensEMBL2id, 4)
plot_grid(selG4.PCA[[1]], selG4.PCA[[2]], ncol = 2)
dev.off()

## Heatmap for these genes
normS.mat                 <- as.data.frame(assay(vsd))
normS.mat$ensembl_gene_id <- rownames(normS.mat)
normS.matE           <- merge(normS.mat, ensEMBL2id, by = "ensembl_gene_id")
normS.matE           <- normS.matE[,c("external_gene_name", "WT_1R_10", "WT_2R_10", "WT_3R_10", 
                                      "H1_1R_10","H1_2R_10", "H1_3R_10", "H2_1R_10","H2_2R_10", "H2_3R_10",
                                      "WT_1R_100", "WT_2R_100", "WT_3R_100",  "H1_1R_100","H1_2R_100", "H1_3R_100",
                                      "H2_1R_100","H2_2R_100", "H2_3R_100")]
normS.matE1          <- normS.matE[normS.matE[,1]%in%genes.WT.new,]
rownames(normS.matE1)<- normS.matE1[,1]
normS.hmatE1         <- normS.matE1[,-1]

breaksList   = seq(6, 12, by = 1)
ScaleCols    <- colorRampPalette(colors = c("purple4","white","darkgreen"))(length(breaksList))
pdf(paste0(Out.dir, "/", Project,"-WT_10vs100_Rm4_sel6Genes_heatmap_new.pdf"), width=6, height=5)
pht_listS1 = Heatmap(as.matrix(normS.hmatE1), col = ScaleCols, 
                     name = "WT_10vsWT_100", show_row_names=T,
                     show_column_names = T, width = unit(8, "cm"),
                     heatmap_legend_param = list(title = "NormExp"),
                     cluster_rows = T,show_row_dend = T,
                     column_title="WT_10vsWT_100",
                     cluster_columns = T,
                     #column_km=3,
                     #column_order=levels(as.factor(colnames(normS.hmatE1))),
                     row_title_rot = 0,
                     row_names_gp = gpar( fontsize = 10))
pht_listS2 = Heatmap(as.matrix(normS.hmatE1[-c(1,5),]), col = ScaleCols, 
                     name = "WT_10vsWT_100", show_row_names=T,
                     show_column_names = T, width = unit(8, "cm"),
                     heatmap_legend_param = list(title = "NormExp"),
                     cluster_rows = T,show_row_dend = T,
                     column_title="WT_10vsWT_100",
                     cluster_columns = T,
                     #column_km=4,
                     #column_order=levels(as.factor(colnames(normS.hmatE1[-c(1,7),]))),
                     row_title_rot = 0,
                     row_names_gp = gpar( fontsize = 10))
print(pht_listS1)

print(pht_listS2)
dev.off()



message("+-----           Common GO pathways between WT_10vsH_10, WT_100vsH_100           ------+")

load(paste0(Out.dir, "/", Project, "-KEGG_GO_p0.05_l2fc1_Analysis_new.RData"))
GO_common_BP_10      <- GO.BP[[1]]$Description[GO.BP[[1]]$Description%in%GO.BP[[2]]$Description==T]
res_common_10        <- resdf.list[[1]]$SYMBOL[resdf.list[[1]]$SYMBOL%in%resdf.list[[2]]$SYMBOL]
GOcommon_10_gene     <- lapply(GO.BP[[1]]$geneID[GO.BP[[1]]$Description%in%GO.BP[[2]]$Description==T],
                               function(x) strsplit(x, split="/")[[1]])
GOres_common_10_gene <- lapply(GOcommon_10_gene, function(x) x[x%in%res_common_10])
GOres_common_10.len  <- unlist(lapply(GOres_common_10_gene, length))
GOresCBP_ratio       <- GOres_common_10.len/GO.BP[[1]]$Count[GO.BP[[1]]$Description%in%GO.BP[[2]]$Description==T]
GO_common_BP_10.sel  <- GO_common_BP_10[GOresCBP_ratio>0.6&GOres_common_10.len>=8] 

GO_common_BP.hmat    <- GO.BP[[1]][GO.BP[[1]]$Description%in%GO_common_BP_10.sel]

## complexheatmap for the genes involved in the above pathways.
sel10_genes          <- lapply(GO_common_BP.hmat$geneID, function(x) strsplit(x, split="/")[[1]])

## group 1: relating to neuron & cell development
sel10_genes1         <- unique(c(unlist(sel10_genes[[3]]), unlist(sel10_genes[[4]]), unlist(sel10_genes[[6]]), unlist(sel10_genes[[7]])))
normC.mat            <- as.data.frame(assay(vsd))
normC.mat$ensembl_gene_id <- rownames(normC.mat)
normC.matE           <- merge(normC.mat, ensEMBL2id, by = "ensembl_gene_id")
normC.matE           <- normC.matE[,c("external_gene_name", "WT_1R_10", "WT_2R_10", "WT_3R_10", 
                                      "H1_1R_10","H1_2R_10", "H1_3R_10", "H2_1R_10","H2_2R_10", "H2_3R_10")]
normC.matE1          <- normC.matE[normC.matE[,1]%in%sel10_genes1,1:10]
rownames(normC.matE1)<- normC.matE1[,1]
normC.hmatE1         <- normC.matE1[,-1]
normC.hmatE1$WT      <- rowMeans(normC.hmatE1[,1:3]); 
normC.hmatE1$H1      <- rowMeans(normC.hmatE1[,4:6]);
normC.hmatE1$H2      <- rowMeans(normC.hmatE1[,7:9])
## group 2 : cell/leukocyte activation
sel10_genes2         <- unique(c(unlist(sel10_genes[[10]]), unlist(sel10_genes[[13]])))
normC.matE2          <- normC.matE[normC.matE[,1]%in%sel10_genes2,1:10]
rownames(normC.matE2)<- normC.matE2[,1]
normC.hmatE2         <- normC.matE2[,-1]
normC.hmatE2$WT      <- rowMeans(normC.hmatE2[,1:3]); 
normC.hmatE2$H1      <- rowMeans(normC.hmatE2[,4:6]);
normC.hmatE2$H2      <- rowMeans(normC.hmatE2[,7:9])

## "epithelial cell proliferation" [1]
sel10_genes3         <- unique(c(unlist(sel10_genes[[1]])))
normC.matE3          <- normC.matE[normC.matE[,1]%in%sel10_genes3,1:10]
rownames(normC.matE3)<- normC.matE3[,1]
normC.hmatE3         <- normC.matE3[,-1]
normC.hmatE3$WT      <- rowMeans(normC.hmatE3[,1:3]); 
normC.hmatE3$H1      <- rowMeans(normC.hmatE3[,4:6]);
normC.hmatE3$H2      <- rowMeans(normC.hmatE3[,7:9])

## [2] "cellular response to interleukin-1"   [5] "response to retinoic acid" [8] "regeneration" 
sel10_genes4           <- list(sel10_genes[[2]], sel10_genes[[5]],sel10_genes[[8]])

normC.matE4            <- lapply(sel10_genes4, function(x) normC.matE[normC.matE[,1]%in%x,1:10])

for(i in 1:length(normC.matE4)){
  rownames(normC.matE4[[i]]) <- normC.matE4[[i]]$external_gene_name
  normC.matE4[[i]]           <- normC.matE4[[i]][,-1]
  normC.matE4[[i]]$WT        <- rowMeans(normC.matE4[[i]][,1:3]); 
  normC.matE4[[i]]$H1        <- rowMeans(normC.matE4[[i]][,4:6]);
  normC.matE4[[i]]$H2        <- rowMeans(normC.matE4[[i]][,7:9])
}

## [9] "skin development" [11] "small molecule catabolic process"   [12] "positive regulation of growth"

sel10_genes5           <- list(sel10_genes[[9]], sel10_genes[[11]],sel10_genes[[12]])

normC.matE5            <- lapply(sel10_genes5, function(x) normC.matE[normC.matE[,1]%in%x,1:10])

for(i in 1:length(normC.matE5)){
  rownames(normC.matE5[[i]]) <- normC.matE5[[i]]$external_gene_name
  normC.matE5[[i]]           <- normC.matE5[[i]][,-1]
  normC.matE5[[i]]$WT        <- rowMeans(normC.matE5[[i]][,1:3]); 
  normC.matE5[[i]]$H1        <- rowMeans(normC.matE5[[i]][,4:6]);
  normC.matE5[[i]]$H2        <- rowMeans(normC.matE5[[i]][,7:9])
}


message("+-----     Heatmap: Common GO pathways between WT_10vsH_10, WT_100vsH_100           ------+")

breaksList   = seq(6, 12, by = 1)
ScaleCols    <- colorRampPalette(colors = c("purple4","white","darkgreen"))(length(breaksList))

pdf(paste0(Out.dir, "/", Project,"-H_WT_10_Rm4_BP_commonG_sel_nureon_heatmap.pdf"))
pht_list1 = Heatmap(as.matrix(normC.hmatE1[,10:12]), col = ScaleCols, 
                    name = "H_10vsWT_10_C1", show_row_names=T,
                    show_column_names = T, width = unit(4, "cm"),
                    heatmap_legend_param = list(title = "NormExp"),
                    cluster_rows = T,show_row_dend = T,
                    column_title="Neuron",
                    cluster_columns = T,
                    column_km=2,
                    column_order=levels(as.factor(colnames(normC.hmatE1[,10:12]))),
                    row_title_rot = 0,
                    row_names_gp = gpar( fontsize = 10))
print(pht_list1)
dev.off()

pdf(paste0(Out.dir, "/", Project,"-H_WT_10_Rm4_BP_commonG_sel_activation_heatmap.pdf"))

pht_list2 = Heatmap(as.matrix(normC.hmatE2[,10:12]), col = ScaleCols, 
                    name = "H_10vsWT_10_C1", show_row_names=T,
                    show_column_names = T, width = unit(4, "cm"),
                    heatmap_legend_param = list(title = "NormExp"),
                    cluster_rows = T,show_row_dend = T,
                    column_title="Cell/leukocyte Activation",
                    cluster_columns = T,
                    column_km=2,
                    column_order=levels(as.factor(colnames(normC.hmatE1[,10:12]))),
                    row_title_rot = 0,
                    row_names_gp = gpar( fontsize = 10))
print(pht_list2)
dev.off()


pdf(paste0(Out.dir, "/", Project,"-H_WT_10_Rm4_BP_commonG_sel_epithelia_heatmap.pdf"))

pht_list3 = Heatmap(as.matrix(normC.hmatE3[,10:12]), col = ScaleCols, 
                    name = "H_10vsWT_10_C1", show_row_names=T,
                    show_column_names = T, width = unit(4, "cm"),
                    heatmap_legend_param = list(title = "NormExp"),
                    cluster_rows = T,show_row_dend = T,
                    column_title="epithelial cell proliferation",
                    cluster_columns = T,
                    column_km=2,
                    column_order=levels(as.factor(colnames(normC.hmatE3[,10:12]))),
                    row_title_rot = 0,
                    row_names_gp = gpar( fontsize = 10))
print(pht_list3)
dev.off()

pdf(paste0(Out.dir, "/", Project,"-H1_WT_10_Rm4_BP_commonG_sel_response_regeneration_heatmap_mean.pdf"))
pht_list4.1 = Heatmap(as.matrix(normC.matE4[[1]][,10:12]), col = ScaleCols, 
                       name = "H_WT_10", show_row_names=T,
                       show_column_names = T, width = unit(4, "cm"),
                       heatmap_legend_param = list(title = "NormExp"),
                       cluster_rows = T,show_row_dend = T,
                       column_title="Response",
                       row_title="cellular response to \ninterleukin-1",
                       cluster_columns = T,
                       column_km=2,
                       column_order=levels(as.factor(colnames(normC.matE4[[1]][,10:12]))),
                       row_title_rot = 90,
                       row_names_gp = gpar( fontsize = 10))

pht_list4.2 = Heatmap(as.matrix(normC.matE4[[2]][,10:12]), col = ScaleCols, 
                      name = "H_WT_10", show_row_names=T,
                      show_column_names = T, width = unit(4, "cm"),
                      heatmap_legend_param = list(title = "NormExp"),
                      cluster_rows = T,show_row_dend = T,
                      column_title="Response",
                      row_title="response to \nretinoic acid",
                      cluster_columns = T,
                      column_km=2,
                      column_order=levels(as.factor(colnames(normC.matE4[[2]][,10:12]))),
                      row_title_rot = 90,
                      row_names_gp = gpar( fontsize = 10))
pht_list4.3 = Heatmap(as.matrix(normC.matE4[[3]][,10:12]), col = ScaleCols, 
                       name = "H_WT_10", show_row_names=T,
                       show_column_names = T, width = unit(4, "cm"),
                       heatmap_legend_param = list(title = "NormExp"),
                       cluster_rows = T,show_row_dend = T,
                       column_title="Response",
                       row_title="regeneration",
                       cluster_columns = T,
                       column_km=2,
                       column_order=levels(as.factor(colnames(normC.matE4[[3]][,10:12]))),
                       row_title_rot = 90,
                       row_names_gp = gpar( fontsize = 10))

ht_list4 = pht_list4.1 %v% pht_list4.2 %v% pht_list4.3
draw(ht_list4)

dev.off()


pdf(paste0(Out.dir, "/", Project,"-H1_WT_10_Rm4_BP_commonG_sel_development_growth_heatmap_mean.pdf"))
pht_list5.1 = Heatmap(as.matrix(normC.matE5[[1]][,10:12]), col = ScaleCols, 
                      name = "H_WT_10", show_row_names=T,
                      show_column_names = T, width = unit(4, "cm"),
                      heatmap_legend_param = list(title = "NormExp"),
                      cluster_rows = T,show_row_dend = T,
                      column_title="Development",
                      row_title="skin development",
                      cluster_columns = T,
                      column_km=2,
                      column_order=levels(as.factor(colnames(normC.matE5[[1]][,10:12]))),
                      row_title_rot = 90,
                      row_names_gp = gpar( fontsize = 10))

pht_list5.2 = Heatmap(as.matrix(normC.matE5[[2]][,10:12]), col = ScaleCols, 
                      name = "H_WT_10", show_row_names=T,
                      show_column_names = T, width = unit(4, "cm"),
                      heatmap_legend_param = list(title = "NormExp"),
                      cluster_rows = T,show_row_dend = T,
                      column_title="Development",
                      row_title="small molecule \ncatabolic process",
                      cluster_columns = T,
                      column_km=2,
                      column_order=levels(as.factor(colnames(normC.matE5[[2]][,10:12]))),
                      row_title_rot = 90,
                      row_names_gp = gpar( fontsize = 10))
pht_list5.3 = Heatmap(as.matrix(normC.matE5[[3]][,10:12]), col = ScaleCols, 
                      name = "H_WT_10", show_row_names=T,
                      show_column_names = T, width = unit(4, "cm"),
                      heatmap_legend_param = list(title = "NormExp"),
                      cluster_rows = T,show_row_dend = T,
                      column_title="Development",
                      row_title="positive regulation \nof growth",
                      cluster_columns = T,
                      column_km=2,
                      column_order=levels(as.factor(colnames(normC.matE5[[3]][,10:12]))),
                      row_title_rot = 90,
                      row_names_gp = gpar( fontsize = 10))

ht_list5 = pht_list5.1 %v% pht_list5.2 %v% pht_list5.3
draw(ht_list5)

dev.off()


message("+------------- Barplot to show up/down regulated for selected pathways for H vs WT BP     -------+")

H1_WT_10.BP   <- GO.BP[[1]][GO.BP[[1]]$Description%in%GO_common_BP_10.sel,]
H1_WT_10.BPG  <- lapply(H1_WT_10.BP$geneID, function(x) unlist(strsplit(x, split="/")))

BP_list_up <- list()

for (i in 1:nrow(H1_WT_10.BP)){
  df_tmp <- resdf.list[[1]][resdf.list[[1]]$SYMBOL %in% H1_WT_10.BPG[[i]],]
  tmp_up <- length(subset(df_tmp[,4], df_tmp[,4] > 0))
  BP_list_up[[i]] <-tmp_up
}

H1_WT_10.BP$UP     <- as.numeric(unlist(BP_list_up))
H1_WT_10.BP$DOWN   <- H1_WT_10.BP$Count - H1_WT_10.BP$UP
H1_WT_10.BP$zscore <- (H1_WT_10.BP$UP-H1_WT_10.BP$DOWN)/sqrt(H1_WT_10.BP$Count)
H1_WT_10.BP$DOWN   <- -H1_WT_10.BP$DOWN


H1_WT_10.BP_molten <- melt(H1_WT_10.BP[,c(2,6:7,9:11)],
                           id.vars=c("Description","p.adjust","qvalue","Count") )

brks <- seq(-5, 25, 5)
lbls = as.character(c(seq(5, 0, -5), seq(5, 25, 5)))
zscore <- round(H1_WT_10.BP$zscore, digits = 3)

pdf(paste0(Out.dir, "/", Project, "-BP_H1_WT_10_CommonselPathways_GeneCount_Barplot_N13.pdf"), width=12, height=6)
barplotSum1 <- ggplot(H1_WT_10.BP_molten, aes(x = reorder(Description, -qvalue), y = value, fill = variable)) +   # Fill column
  geom_bar(stat = "identity", width = 0.8) +   # draw the bars
  scale_y_continuous(breaks = brks,   # Breaks
                     labels = lbls,
                     limits=c(-5,25)) + # Labels
  
  coord_flip() +  # Flip axes
  #labs(title="Biological Process") +
  xlab("") +
  ylab("") +
  guides(fill=FALSE) +
  theme(plot.title = element_text(hjust = .2),
        axis.ticks = element_blank()) +   # Centre plot title
  scale_fill_manual(values = c("red", "blue"))+  # Color palette
  theme_bw() +
  theme_update(axis.title.x = element_text(size=14, face= "bold"),
               axis.text.x = element_text(size=14, face="bold"),
               axis.title.y = element_text(size=14, face= "bold"),
               axis.text.y = element_text(size=14, face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),panel.background = element_rect(fill = "white", colour = NA)) 


barplotSum1

dev.off()

## H2_10 vs WT_10
H2_WT_10.BP   <- GO.BP[[2]][GO.BP[[2]]$Description%in%GO_common_BP_10.sel,]
H2_WT_10.BPG  <- lapply(H2_WT_10.BP$geneID, function(x) unlist(strsplit(x, split="/")))

BP_list_up <- list()

for (i in 1:nrow(H2_WT_10.BP)){
  df_tmp <- resdf.list[[2]][resdf.list[[2]]$SYMBOL %in% H2_WT_10.BPG[[i]],]
  tmp_up <- length(subset(df_tmp[,4], df_tmp[,4] > 0))
  BP_list_up[[i]] <-tmp_up
}

H2_WT_10.BP$UP     <- as.numeric(unlist(BP_list_up))
H2_WT_10.BP$DOWN   <- H2_WT_10.BP$Count - H2_WT_10.BP$UP
H2_WT_10.BP$zscore <- (H2_WT_10.BP$UP-H2_WT_10.BP$DOWN)/sqrt(H2_WT_10.BP$Count)
H2_WT_10.BP$DOWN   <- -H2_WT_10.BP$DOWN


H2_WT_10.BP_molten <- melt(H2_WT_10.BP[,c(2,6:7,9:11)],
                           id.vars=c("Description","p.adjust","qvalue","Count") )


pdf(paste0(Out.dir, "/", Project, "-BP_H2_WT_10_CommonselPathways_GeneCount_Barplot_N13.pdf"), width=12, height=6)
barplotSum2 <- ggplot(H2_WT_10.BP_molten, aes(x = reorder(Description, -qvalue), y = value, fill = variable)) +   # Fill column
  geom_bar(stat = "identity", width = 0.8) +   # draw the bars
  scale_y_continuous(breaks = brks,   # Breaks
                     labels = lbls,
                     limits=c(-5,25)) + # Labels
  
  coord_flip() +  # Flip axes
  #labs(title="Biological Process") +
  xlab("") +
  ylab("") +
  guides(fill=FALSE) +
  theme(plot.title = element_text(hjust = .2),
        axis.ticks = element_blank()) +   # Centre plot title
  scale_fill_manual(values = c("red", "blue"))+  # Color palette
  theme_bw() +
  theme_update(axis.title.x = element_text(size=14, face= "bold"),
               axis.text.x = element_text(size=14, face="bold"),
               axis.title.y = element_text(size=14, face= "bold"),
               axis.text.y = element_text(size=14, face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),panel.background = element_rect(fill = "white", colour = NA)) 


barplotSum2

dev.off()


message("+----------     Heatmap for common GO KEGG, MF and CC with Common Genes    ----------------------+")

GO_common_CC_10      <- GO.CC[[1]]$Description[GO.CC[[1]]$Description%in%GO.CC[[2]]$Description==T]
GOcommonC_10_gene    <- lapply(GO.CC[[1]]$geneID[GO.CC[[1]]$Description%in%GO.CC[[2]]$Description==T],
                                function(x) strsplit(x, split="/")[[1]])
GOres_common_10_gene <- lapply(GOcommonC_10_gene, function(x) x[x%in%res_common_10])
GOres_common_10.len  <- unlist(lapply(GOres_common_10_gene, length))

normC.matC           <- lapply(GOres_common_10_gene, function(x) normC.matE[normC.matE[,1]%in%x,1:10])

for(i in 1:length(normC.matC)){
  rownames(normC.matC[[i]]) <- normC.matC[[i]]$external_gene_name
  normC.matC[[i]]           <- normC.matC[[i]][,-1]
  normC.matC[[i]]$WT        <- rowMeans(normC.matC[[i]][,1:3]); 
  normC.matC[[i]]$H1        <- rowMeans(normC.matC[[i]][,4:6]);
  normC.matC[[i]]$H2        <- rowMeans(normC.matC[[i]][,7:9])
}

pht_listCC <- list();
pdf(paste0(Out.dir, "/", Project,"-H1_WT_10_Rm4_CC_commonG_sel_N7_heatmap_mean.pdf"))
for(i in 1:7){
  pht_listCC[[i]] = Heatmap(as.matrix(normC.matC[[i]][,10:12]), col = ScaleCols, 
                         name = "H1_H2_WT_10", show_row_names=T,
                         show_column_names = T, width = unit(4, "cm"),
                         heatmap_legend_param = list(title = "NormExp"),
                         cluster_rows = T,show_row_dend = T,
                         column_title="Cellular Component",
                         row_title=GO_common_CC_10[i],
                         cluster_columns = T,
                         column_km=2,
                         column_order=levels(as.factor(colnames(normC.matC[[i]][,10:12]))),
                         row_title_rot = 0,
                         row_names_gp = gpar( fontsize = 10))
}


ht_list = pht_listCC[[1]] %v% pht_listCC[[2]] %v% pht_listCC[[3]] %v% pht_listCC[[4]] %v% pht_listCC[[5]] %v% pht_listCC[[6]] %v% pht_listCC[[7]]
draw(ht_list)

dev.off()

### KEGG
GO_common_KEGG_10      <- GO.KEGG[[1]]$Description[GO.KEGG[[1]]$Description%in%GO.KEGG[[2]]$Description==T]
GOcommonK_10_gene      <- lapply(GO.KEGG[[1]]$geneID[GO.KEGG[[1]]$Description%in%GO.KEGG[[2]]$Description==T],
                                function(x) strsplit(x, split="/")[[1]])

GOres_common_10_gene <- lapply(GOcommonK_10_gene, function(x) x[x%in%res_common_10])
GOres_common_10.len  <- unlist(lapply(GOres_common_10_gene, length))

normC.matK           <- lapply(GOres_common_10_gene, function(x) normC.matE[normC.matE[,1]%in%x,1:10])

for(i in 1:length(normC.matK)){
  rownames(normC.matK[[i]]) <- normC.matK[[i]]$external_gene_name
  normC.matK[[i]]           <- normC.matK[[i]][,-1]
  normC.matK[[i]]$WT        <- rowMeans(normC.matK[[i]][,1:3]); 
  normC.matK[[i]]$H1        <- rowMeans(normC.matK[[i]][,4:6]);
  normC.matK[[i]]$H2        <- rowMeans(normC.matK[[i]][,7:9])
}

pht_listKK <- list();
pdf(paste0(Out.dir, "/", Project,"-H1_WT_10_Rm4_KEGG_commonG_sel_N5_heatmap_mean.pdf"))
for(i in 3:7){
  pht_listKK[[i]] = Heatmap(as.matrix(normC.matK[[i]][,10:12]), col = ScaleCols, 
                            name = "H1_H2_WT_10", show_row_names=T,
                            show_column_names = T, width = unit(4, "cm"),
                            heatmap_legend_param = list(title = "NormExp"),
                            cluster_rows = T,show_row_dend = T,
                            column_title="KEGG",
                            row_title=GO_common_KEGG_10[i],
                            cluster_columns = T,
                            column_km=2,
                            column_order=levels(as.factor(colnames(normC.matK[[i]][,10:12]))),
                            row_title_rot = 0,
                            row_names_gp = gpar( fontsize = 10))
}


ht_list = pht_listKK[[3]] %v% pht_listKK[[4]] %v% pht_listKK[[5]] %v% pht_listKK[[6]] %v% pht_listKK[[7]]
draw(ht_list)

dev.off()

message("+--------------------       H1/2_100 vs WT_100 Common Pathways             ---------------------+")

GO_common_BP_100      <- GO.BP[[3]]$Description[GO.BP[[3]]$Description%in%GO.BP[[4]]$Description==T]
res_common_100        <- resdf.list[[3]]$SYMBOL[resdf.list[[3]]$SYMBOL%in%resdf.list[[4]]$SYMBOL]
GOcommon_100_gene     <- lapply(GO.BP[[3]]$geneID[GO.BP[[3]]$Description%in%GO.BP[[4]]$Description==T],
                               function(x) strsplit(x, split="/")[[1]])
GOres_common_100_gene <- lapply(GOcommon_100_gene, function(x) x[x%in%res_common_100])
GOres_common_100.len  <- unlist(lapply(GOres_common_100_gene, length))
GOresCBP_ratio        <- GOres_common_100.len/GO.BP[[3]]$Count[GO.BP[[3]]$Description%in%GO.BP[[4]]$Description==T]

GO_common_BP_100.sel  <- GO_common_BP_100[GOresCBP_ratio>0.6&GOres_common_100.len>=5] 

message("+------------- Barplot to show up/down regulated for selected pathways for H vs WT BP     -------+")

H1_WT_100.BP   <- GO.BP[[3]][GO.BP[[3]]$Description%in%GO_common_BP_100.sel,]
H1_WT_100.BPG  <- lapply(H1_WT_100.BP$geneID, function(x) unlist(strsplit(x, split="/")))

BP_list_up <- list()

for (i in 1:nrow(H1_WT_100.BP)){
  df_tmp <- resdf.list[[3]][resdf.list[[3]]$SYMBOL %in% H1_WT_100.BPG[[i]],]
  tmp_up <- length(subset(df_tmp[,4], df_tmp[,4] > 0))
  BP_list_up[[i]] <-tmp_up
}

H1_WT_100.BP$UP     <- as.numeric(unlist(BP_list_up))
H1_WT_100.BP$DOWN   <- H1_WT_100.BP$Count - H1_WT_100.BP$UP
H1_WT_100.BP$zscore <- (H1_WT_100.BP$UP-H1_WT_100.BP$DOWN)/sqrt(H1_WT_100.BP$Count)
H1_WT_100.BP$DOWN   <- -H1_WT_100.BP$DOWN


H1_WT_100.BP_molten <- melt(H1_WT_100.BP[,c(2,6:7,9:11)],
                           id.vars=c("Description","p.adjust","qvalue","Count") )

brks <- seq(-5, 10, 5)
lbls = as.character(c(seq(5, 0, -5), seq(5, 10, 5)))
zscore <- round(H1_WT_100.BP$zscore, digits = 3)

pdf(paste0(Out.dir, "/", Project, "-BP_H1_WT_100_CommonselPathways_GeneCount_Barplot_N27.pdf"), width=12, height=16)
barplotSum1 <- ggplot(H1_WT_100.BP_molten, aes(x = reorder(Description, -qvalue), y = value, fill = variable)) +   # Fill column
  geom_bar(stat = "identity", width = 0.8) +   # draw the bars
  scale_y_continuous(breaks = brks,   # Breaks
                     labels = lbls,
                     limits=c(-5,25)) + # Labels
  
  coord_flip() +  # Flip axes
  #labs(title="Biological Process") +
  xlab("") +
  ylab("") +
  guides(fill=FALSE) +
  theme(plot.title = element_text(hjust = .2),
        axis.ticks = element_blank()) +   # Centre plot title
  scale_fill_manual(values = c("red", "blue"))+  # Color palette
  theme_bw() +
  theme_update(axis.title.x = element_text(size=14, face= "bold"),
               axis.text.x = element_text(size=14, face="bold"),
               axis.title.y = element_text(size=14, face= "bold"),
               axis.text.y = element_text(size=14, face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),panel.background = element_rect(fill = "white", colour = NA)) 


barplotSum1

dev.off()

## H2_100 vs WT_100
H2_WT_100.BP   <- GO.BP[[4]][GO.BP[[4]]$Description%in%GO_common_BP_100.sel,]
H2_WT_100.BPG  <- lapply(H2_WT_100.BP$geneID, function(x) unlist(strsplit(x, split="/")))

BP_list_up <- list()

for (i in 1:nrow(H2_WT_100.BP)){
  df_tmp <- resdf.list[[4]][resdf.list[[4]]$SYMBOL %in% H2_WT_100.BPG[[i]],]
  tmp_up <- length(subset(df_tmp[,4], df_tmp[,4] > 0))
  BP_list_up[[i]] <-tmp_up
}

H2_WT_100.BP$UP     <- as.numeric(unlist(BP_list_up))
H2_WT_100.BP$DOWN   <- H2_WT_100.BP$Count - H2_WT_100.BP$UP
H2_WT_100.BP$zscore <- (H2_WT_100.BP$UP-H2_WT_100.BP$DOWN)/sqrt(H2_WT_100.BP$Count)
H2_WT_100.BP$DOWN   <- -H2_WT_100.BP$DOWN


H2_WT_100.BP_molten <- melt(H2_WT_100.BP[,c(2,6:7,9:11)],
                           id.vars=c("Description","p.adjust","qvalue","Count") )


pdf(paste0(Out.dir, "/", Project, "-BP_H2_WT_100_CommonselPathways_GeneCount_Barplot_N27.pdf"), width=12, height=16)
barplotSum2 <- ggplot(H2_WT_100.BP_molten, aes(x = reorder(Description, -qvalue), y = value, fill = variable)) +   # Fill column
  geom_bar(stat = "identity", width = 0.8) +   # draw the bars
  scale_y_continuous(breaks = brks,   # Breaks
                     labels = lbls,
                     limits=c(-5,25)) + # Labels
  
  coord_flip() +  # Flip axes
  #labs(title="Biological Process") +
  xlab("") +
  ylab("") +
  guides(fill=FALSE) +
  theme(plot.title = element_text(hjust = .2),
        axis.ticks = element_blank()) +   # Centre plot title
  scale_fill_manual(values = c("red", "blue"))+  # Color palette
  theme_bw() +
  theme_update(axis.title.x = element_text(size=14, face= "bold"),
               axis.text.x = element_text(size=14, face="bold"),
               axis.title.y = element_text(size=14, face= "bold"),
               axis.text.y = element_text(size=14, face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),panel.background = element_rect(fill = "white", colour = NA)) 


barplotSum2

dev.off()

message("+----------     Heatmap for common GO KEGG, MF and CC with Common Genes    ----------------------+")

GO_common_CC_100      <- GO.CC[[3]]$Description[GO.CC[[3]]$Description%in%GO.CC[[4]]$Description==T]
GOcommonC_100_gene    <- lapply(GO.CC[[3]]$geneID[GO.CC[[3]]$Description%in%GO.CC[[4]]$Description==T],
                               function(x) strsplit(x, split="/")[[1]])
GOres_common_100_gene <- lapply(GOcommonC_100_gene, function(x) x[x%in%res_common_100])
GOres_common_100.len  <- unlist(lapply(GOres_common_100_gene, length))

normC.matC           <- lapply(GOres_common_100_gene, function(x) normC.matE[normC.matE[,1]%in%x,1:10])

for(i in 1:length(normC.matC)){
  rownames(normC.matC[[i]]) <- normC.matC[[i]]$external_gene_name
  normC.matC[[i]]           <- normC.matC[[i]][,-1]
  normC.matC[[i]]$WT        <- rowMeans(normC.matC[[i]][,1:3]); 
  normC.matC[[i]]$H1        <- rowMeans(normC.matC[[i]][,4:6]);
  normC.matC[[i]]$H2        <- rowMeans(normC.matC[[i]][,7:9])
}

pht_listCC <- list();
pdf(paste0(Out.dir, "/", Project,"-H1_WT_100_Rm4_CC_commonG_sel_N3_heatmap_mean.pdf"))
for(i in 1:3){
  pht_listCC[[i]] = Heatmap(as.matrix(normC.matC[[i]][,10:12]), col = ScaleCols, 
                            name = "H1_H2_WT_100", show_row_names=T,
                            show_column_names = T, width = unit(4, "cm"),
                            heatmap_legend_param = list(title = "NormExp"),
                            cluster_rows = T,show_row_dend = T,
                            column_title="Cellular Component",
                            row_title=GO_common_CC_100[i],
                            cluster_columns = T,
                            column_km=2,
                            column_order=levels(as.factor(colnames(normC.matC[[i]][,10:12]))),
                            row_title_rot = 0,
                            row_names_gp = gpar( fontsize = 10))
}


ht_list = pht_listCC[[1]] %v% pht_listCC[[2]] %v% pht_listCC[[3]] 
draw(ht_list)

dev.off()

### KEGG
GO_common_KEGG_100      <- GO.KEGG[[3]]$Description[GO.KEGG[[3]]$Description%in%GO.KEGG[[4]]$Description==T]
GOcommonK_100_gene      <- lapply(GO.KEGG[[3]]$geneID[GO.KEGG[[3]]$Description%in%GO.KEGG[[4]]$Description==T],
                                 function(x) strsplit(x, split="/")[[1]])

GOres_common_100_gene <- lapply(GOcommonK_100_gene, function(x) x[x%in%res_common_100])
GOres_common_100.len  <- unlist(lapply(GOres_common_100_gene, length))

normC.matK           <- lapply(GOres_common_100_gene, function(x) normC.matE[normC.matE[,1]%in%x,1:10])

for(i in 1:length(normC.matK)){
  rownames(normC.matK[[i]]) <- normC.matK[[i]]$external_gene_name
  normC.matK[[i]]           <- normC.matK[[i]][,-1]
  normC.matK[[i]]$WT        <- rowMeans(normC.matK[[i]][,1:3]); 
  normC.matK[[i]]$H1        <- rowMeans(normC.matK[[i]][,4:6]);
  normC.matK[[i]]$H2        <- rowMeans(normC.matK[[i]][,7:9])
}

pht_listKK <- list();
GO_common_KEGG_100_names <- c("PPAR signaling \npathway", "AGE-RAGE signaling \npathway in \ndiabetic complications")
pdf(paste0(Out.dir, "/", Project,"-H1_WT_100_Rm4_KEGG_commonG_sel_N2_heatmap_mean.pdf"))
for(i in 1:2){
  pht_listKK[[i]] = Heatmap(as.matrix(normC.matK[[i]][,10:12]), col = ScaleCols, 
                            name = "H1_H2_WT_100", show_row_names=T,
                            show_column_names = T, width = unit(4, "cm"),
                            heatmap_legend_param = list(title = "NormExp"),
                            cluster_rows = T,show_row_dend = T,
                            column_title="KEGG",
                            row_title=GO_common_KEGG_100_names[i],
                            cluster_columns = T,
                            column_km=2,
                            column_order=levels(as.factor(colnames(normC.matK[[i]][,10:12]))),
                            row_title_rot = 90,
                            row_names_gp = gpar( fontsize = 8))
}


ht_list = pht_listKK[[1]] %v% pht_listKK[[2]] 
draw(ht_list)

dev.off()


message("+-----------------          Transcriptor Factor enrichment Analysis               -------------------------+")
## https://github.com/Dowell-Lab/TFEA
## https://www.biorxiv.org/content/10.1101/2020.01.25.919738v3.full
## Transcription factor enrichment analysis (TFEA): Quantifying the activity of hundreds of transcription factors from a single experiment
## http://gtrd.biouml.org/#! download the TF data

TFs  <- read.csv("Rat_TFlist_N1064.csv", header=T)[,1:9]



message("+-----------------         Finish 26/10/2023                          -------------------------+")


























message("+----- Try to make the network for BP selected pathways -------------+")

GO.BP.new   <- list()
GO.MF.new   <- list()
GO.CC.new   <- list()
GO.KEGG <- list()


for(a in 1:4){
  GO.BP.new[[a]] <- enrichGO(gene          = resdf.list[[a]]$ENTREZID,
                             OrgDb         = org.Rn.eg.db,
                             keyType       = 'ENTREZID',
                             ont           = "BP",
                             pAdjustMethod = "BH",
                             minGSSize     = 10,
                             pvalueCutoff  = 0.05 )
  GO.BP.new[[a]] <- setReadable(GO.BP.new[[a]], OrgDb = org.Rn.eg.db, keyType = "ENTREZID")
  
  GO.CC.new[[a]] <- enrichGO(gene          = resdf.list[[a]]$ENTREZID,
                             OrgDb         = org.Rn.eg.db,
                             keyType       = 'ENTREZID',
                             ont           = "CC",
                             pAdjustMethod = "BH",
                             minGSSize     = 10,
                             pvalueCutoff  = 0.05 )
  GO.CC.new[[a]] <- setReadable(GO.CC.new[[a]], OrgDb = org.Rn.eg.db, keyType = "ENTREZID")
  GO.MF.new[[a]] <- enrichGO(gene          = resdf.list[[a]]$ENTREZID,
                             OrgDb         = org.Rn.eg.db,
                             keyType       = 'ENTREZID',
                             ont           = "MF",
                             pAdjustMethod = "BH",
                             minGSSize     = 10,
                             pvalueCutoff  = 0.05 )
  GO.MF.new[[a]] <- setReadable(GO.MF.new[[a]], OrgDb = org.Rn.eg.db, keyType = "ENTREZID")
  GO.KEGG[[a]] <-  enrichKEGG(resdf.list[[a]]$ENTREZID, organism = "rno")  
  GO.KEGG[[a]] <- setReadable(GO.KEGG[[a]], OrgDb = org.Rn.eg.db, keyType = "ENTREZID")
  
}

save(GO.BP.new, GO.CC.new, GO.MF.new, GO.KEGG,file=paste0(Out.dir, "/", Project, "-KEGG_GO_p0.05_l2fc1_Analysis_EntrezID.RData"))

resfold.List     <- lapply(resdf.list, function(x) x$L2FC)
for(i in 1:4){
  names(resfold.List[[1]]) <- resdf.list[[1]]$ENTREZID
}
pdf(paste0(Out.dir, "/", Project,"-H1_WT_10_Rm4_BP_commonG_sel_cnetplot.pdf"), width=15, height = 10)
pH110.BP            <- cnetplot(GO.BP.new[[1]], foldChange=resfold.List[[1]],
                                showCategory=GO_common_BP_10.sel[-c(1:2,6,11)])
print(pH110.BP)
dev.off()

pdf(paste0(Out.dir, "/", Project,"-H1_WT_10_Rm4_CC_commonG_sel_cnetplot.pdf"), width=15, height = 10)
pH110.CC            <- cnetplot(GO.CC.new[[1]], foldChange=resfold.List[[1]],
                                showCategory=GO_common_CC_10)
print(pH110.CC)
dev.off()

pdf(paste0(Out.dir, "/", Project,"-H1_WT_10_Rm4_KEGG_commonG_sel_cnetplot.pdf"), width=15, height = 10)
pH110.KEGG            <- cnetplot(GO.KEGG[[1]], foldChange=resfold.List[[1]],
                                  showCategory=GO_common_KEGG_10)
print(pH110.KEGG)
dev.off()

message("+-----------WT_100,H1,H2       -------------------+")
GO_common_BP_100      <- GO.BP[[3]]$Description[GO.BP[[3]]$Description%in%GO.BP[[4]]$Description==T]
res_common_100        <- resdf.list[[3]]$SYMBOL[resdf.list[[3]]$SYMBOL%in%resdf.list[[4]]$SYMBOL]
GOcommon_100_gene     <- lapply(GO.BP[[3]]$geneID[GO.BP[[3]]$Description%in%GO.BP[[4]]$Description==T],
                                function(x) strsplit(x, split="/")[[1]]) ## 147
GOres_common_100_gene <- lapply(GOcommon_100_gene, function(x) x[x%in%res_common_100])
GOres_common_100.len  <- unlist(lapply(GOres_common_100_gene, length))
GOresCBP_ratio        <- GOres_common_100.len/GO.BP[[3]]$Count[GO.BP[[3]]$Description%in%GO.BP[[4]]$Description==T]
GO_common_BP_100.sel  <- GO_common_BP_100[GOresCBP_ratio>0.6&GOres_common_100.len>=5] 
GO_common_BP.hmat    <- GO.BP[[3]][GO.BP[[3]]$Description%in%GO_common_BP_100.sel]

sel100_genes          <- lapply(GO_common_BP.hmat$geneID, function(x) strsplit(x, split="/")[[1]])

## merge development, metabolic process, cell proliferation, 
sel100_genes1         <- unique(c(unlist(sel100_genes[[7]]), unlist(sel100_genes[[8]]), unlist(sel100_genes[[12]]),
                                  unlist(sel100_genes[[17]]), unlist(sel100_genes[[19]]),unlist(sel100_genes[[20]]),
                                  unlist(sel100_genes[[21]]), unlist(sel100_genes[[22]])))
sel100_genes2         <- unique(c(unlist(sel100_genes[[2]]), unlist(sel100_genes[[3]]), unlist(sel100_genes[[4]]),
                                  unlist(sel100_genes[[5]]), unlist(sel100_genes[[15]]),unlist(sel100_genes[[16]]),
                                  unlist(sel100_genes[[18]])))
sel100_genes3         <- unique(c(unlist(sel100_genes[[9]]), unlist(sel100_genes[[10]]), unlist(sel100_genes[[11]])))

sel100_genes4         <- unlist(sel100_genes[[6]])
sel100_genes5         <- unlist(sel100_genes[[13]])
sel100_genes6         <- unlist(sel100_genes[[14]])

sel100_geneBP_List    <- list(sel100_genes1,sel100_genes2,sel100_genes3,sel100_genes4,sel100_genes5, sel100_genes6)
pht.100BP <- list()
rtitles.BP <- c("Development", "Metabolic Process", "Cell Proliferation", "response to estradiol",
                "skeletal system \nmorphogenesis", "actomyosin structure \norganization")
for(i in 1:6){
  mean.hmatE1           <- normC.matE[normC.matE[,1]%in%sel100_geneBP_List[[i]],]
  rownames(mean.hmatE1) <- mean.hmatE1[,1]
  mean.hmatE2       <- mean.hmatE1[,-1]
  mean.hmatE2$WT    <- rowMeans(mean.hmatE2[,1:3])
  mean.hmatE2$H1    <- rowMeans(mean.hmatE2[,4:6])
  mean.hmatE2$H2    <- rowMeans(mean.hmatE2[,7:9])
  mean.hmatE2.plt   <- mean.hmatE2[,10:12]
  
  pht.100BP[[i]] = Heatmap(as.matrix(mean.hmatE2.plt), col = ScaleCols, 
                           name = "H1_H2_WT_100", show_row_names=T,
                           show_column_names = T, width = unit(4, "cm"),
                           heatmap_legend_param = list(title = "NormExp"),
                           cluster_rows = T,show_row_dend = T,
                           column_title="Biological Process",
                           row_title=rtitles.BP[i],
                           cluster_columns = T,
                           column_km=2,
                           column_order=levels(as.factor(colnames(mean.hmatE2.plt))),
                           row_title_rot = 0,
                           row_names_gp = gpar( fontsize = 10))
  pht.100BP
}

pdf(paste0(Out.dir, "/", Project,"-H1_WT_100_Rm4_BP_commonG_development_metabolic_heatmap_mean.pdf"))
htlist = pht.100BP[[1]] %v% pht.100BP[[2]] 
draw(htlist)
dev.off()



pdf(paste0(Out.dir, "/", Project,"-H1_WT_100_Rm4_BP_commonG_cellProf_estradiol_morph_orga_heatmap_mean.pdf"))
htlist = pht.100BP[[3]] %v% pht.100BP[[4]] %v% pht.100BP[[5]] %v% pht.100BP[[6]]
draw(htlist)
dev.off()


message("+--- CC heatmap----------+")

GO_common_CC_100      <- GO.CC[[3]]$Description[GO.CC[[3]]$Description%in%GO.CC[[4]]$Description==T]
GOcommonC_100_gene    <- lapply(GO.CC[[3]]$geneID[GO.CC[[3]]$Description%in%GO.CC[[4]]$Description==T],
                                function(x) strsplit(x, split="/")[[1]])
sel100_genes1         <- unlist(GOcommonC_100_gene[[1]])
sel100_genes2         <- unlist(GOcommonC_100_gene[[2]])
sel100_genes3         <- unlist(GOcommonC_100_gene[[3]])

sel100_geneCC_List    <- list(sel100_genes1,sel100_genes2,sel100_genes3)
pht.100CC <- list()
rtitles.CC <- c("extracellular matrix", "external side of plasma membrane", 
                "collagen trimer")
for(i in 1:3){
  mean.hmatE1           <- normC.matE[normC.matE[,1]%in%sel100_geneCC_List[[i]],]
  rownames(mean.hmatE1) <- mean.hmatE1[,1]
  mean.hmatE2       <- mean.hmatE1[,-1]
  mean.hmatE2$WT    <- rowMeans(mean.hmatE2[,1:3])
  mean.hmatE2$H1    <- rowMeans(mean.hmatE2[,4:6])
  mean.hmatE2$H2    <- rowMeans(mean.hmatE2[,7:9])
  mean.hmatE2.plt   <- mean.hmatE2[,10:12]
  
  pht.100CC[[i]] = Heatmap(as.matrix(mean.hmatE2.plt), col = ScaleCols, 
                           name = "H1_H2_WT_100", show_row_names=T,
                           show_column_names = T, width = unit(4, "cm"),
                           heatmap_legend_param = list(title = "NormExp"),
                           cluster_rows = T,show_row_dend = T,
                           column_title="Cellular Component",
                           row_title=rtitles.CC[i],
                           cluster_columns = T,
                           column_km=2,
                           column_order=levels(as.factor(colnames(mean.hmatE2.plt))),
                           row_title_rot = 0,
                           row_names_gp = gpar( fontsize = 10))
  pht.100BP
}

pdf(paste0(Out.dir, "/", Project,"-H1_WT_100_Rm4_CC_commonG_3GO_heatmap_mean.pdf"))
htlist = pht.100CC[[1]] %v% pht.100CC[[2]] %v% pht.100CC[[3]]
draw(htlist)
dev.off()



















message("+---Check the Rank 1 first, common genes-common Pathways H1_10vsWT_10, H2_10vsWT_10, ----------------+")

set.seed(123)
library(circlize)
df = generateRandomBed(30)[ , 1:3]
df[, 3] = df[, 2]
df$gene = paste0("gene_", 1:nrow(df))
df$pathway = paste0("pathway_", sample(10, nrow(df), replace = TRUE))
head(df)

cytoband = read.cytoband()$df
cytoband = rbind(cytoband,
                 data.frame(V1 = "pathway", V2 = 1,  V3 = 2e8, V4 = "", V5 = "")
)
tail(cytoband)
foo = round(seq(1, 2e8, length = 11))
pathway_mid = structure(as.integer((foo[1:10] + foo[2:11])/2), 
                        names = paste0("pathway_", 1:10))

df$pathway_chr = "pathway"
df$pathway_start = pathway_mid[ df$pathway ]
df$pathway_end = pathway_mid[ df$pathway ]
head(df)
circos.initializeWithIdeogram(cytoband)
circos.genomicLink(df[, 1:3], df[, 6:8])

##
circos.par(gap.after = c(rep(1, 23), 5, 5))
circos.genomicInitialize(cytoband, plotType = NULL)

# the labels track
label_df = rbind(df[, 1:4], setNames(df[, c(6:8, 5)], colnames(df)[1:4]))
label_df = unique(label_df)
circos.genomicLabels(label_df, labels.column = 4, side = "outside", cex = 0.6)

# the chromosome names track
circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              niceFacing = TRUE, adj = c(0.5, 0), cex = 0.8)
}, track.height = strheight("fj", cex = 0.8)*1.2, bg.border = NA, cell.padding = c(0, 0, 0, 0))

# ideogram track
circos.genomicIdeogram(cytoband)

# add points to the cell "pathway"
tb = table(df$pathway)
pathway_col = structure(1:10, names = paste0("pathway_", 1:10))
set.current.cell(sector.index = "pathway", track.index = get.current.track.index())
circos.points(x = pathway_mid[names(tb)], y = CELL_META$ycenter, pch = 16,
              cex = tb/5, col = pathway_col[names(tb)])

# genomic links
circos.genomicLink(df[, 1:3], df[, 6:8], col = pathway_col[df[, 5]])
circos.clear()
message("+----- Group some pathways or use Enrich pathways to select the interesting ones and generate heatmap and network or barplot -------+")
install.packages("magick")
library("magick")
library(simplifyEnrichment)
go_id1 = sel.BP_H1_WT_10.mer$ID
go_id2 = sel.BP_H2_WT_10.mer$ID
go_id3 = sel.BP_H1_WT_100.mer$ID
go_id4 = sel.BP_H2_WT_100.mer$ID
mat1 = GO_similarity(go_id1,ont="BP", db = 'org.Rn.eg.db', measure = "Rel")
mat2 = GO_similarity(go_id2,ont="BP", db = 'org.Rn.eg.db', measure = "Rel") 
mat3 = GO_similarity(go_id3,ont="BP", db = 'org.Rn.eg.db', measure = "Rel")
mat4 = GO_similarity(go_id3,ont="BP", db = 'org.Rn.eg.db', measure = "Rel")

pdf(paste0(Out.dir, "/SimplifyGO_BP_H12_WT_10_100.pdf"))
simplifyGO(mat1)
simplifyGO(mat2)
simplifyGO(mat3)
simplifyGO(mat4)
dev.off()










message("+---Check the Rank 1 first, common genes-common Pathways H1_10vsWT_10, H2_10vsWT_10, ----------------+")











message("+-----   March, 2021 check DEGs in WT_10 vs WT_100 with H1, H2 DEGs between 10 vs 100  ---------------+")
SgeneList  <- c("Ddx3y","Ca3","Fgf21","Cyp2c12","Cyp2c13", "Kdm5c","Igfbp1","Kdm5d")
## Ddx3y, kdm5d, Ca3,Fgf21,Cyp2C13-male
## Cyp2C12, Igfbp1, Kdm5c-Female
vsd.assay  <- assay(vsd.1)
Sgindex    <- ensEMBL2id[ensEMBL2id$external_gene_name%in%SgeneList,]
vsd.selXY  <- vsd.assay[rownames(vsd.assay)%in%Sgindex[,1],]

hk.1            <- plotCounts(dds, intgroup=c("label"), normalized=TRUE, returnData=TRUE,gene="ENSRNOG00000002501")
colnames(hk.1)  <- c(SgeneList[1], "label")
hk.2            <- plotCounts(dds, intgroup=c("label"), normalized=TRUE, returnData=TRUE,gene="ENSRNOG00000010079")
colnames(hk.2)  <- c(SgeneList[2], "label")
hk.3            <- plotCounts(dds, intgroup=c("label"), normalized=TRUE, returnData=TRUE,gene="ENSRNOG00000020990")
colnames(hk.3)  <- c(SgeneList[3], "label")
hk.4            <- plotCounts(dds, intgroup=c("label"), normalized=TRUE, returnData=TRUE,gene="ENSRNOG00000047945")
colnames(hk.4)  <- c(SgeneList[4], "label")
hk.5            <- plotCounts(dds, intgroup=c("label"), normalized=TRUE, returnData=TRUE,gene="ENSRNOG00000049464")
colnames(hk.5)  <- c(SgeneList[5], "label")
hk.6            <- plotCounts(dds, intgroup=c("label"), normalized=TRUE, returnData=TRUE,gene="ENSRNOG00000057706")
colnames(hk.6)  <- c(SgeneList[6], "label")
hk.7            <- plotCounts(dds, intgroup=c("label"), normalized=TRUE, returnData=TRUE,gene="ENSRNOG00000058780")
colnames(hk.7)  <- c(SgeneList[7], "label")
hk.8            <- plotCounts(dds, intgroup=c("label"), normalized=TRUE, returnData=TRUE,gene="ENSRNOG00000060496")
colnames(hk.8)  <- c(SgeneList[8], "label")


xtest           <- merge(hk.1, hk.2, by=c("row.names", "label") )
colnames(xtest)[3:4] <- c(SgeneList[1:2])
xtest$Fgf21     <- hk.3$Fgf21 
xtest$Cyp2c12   <- hk.4$Cyp2c12
xtest$Cp2c13    <- hk.5$Cp2c13
xtest$Kdm5c     <- hk.6$Kdm5c
xtest$Igfbp1    <- hk.7$Igfbp1
xtest$Kdm5d     <- hk.8$Kdm5d

xtest           <- xtest[c(13,15,17,1,3,5,7,9,11,14,16,18,2,4,6,8,10,12),]

xtest.m           <- melt(xtest)
colnames(xtest.m) <- c("Sample", "group", "Gene", "normCount")
xtest.m$Sample    <- factor(xtest.m$Sample,levels = c(paste0("WT_", c(1:3),"R_10"),
                                                      paste0("H1_", c(1:3),"R_10"),
                                                      paste0("H2_", c(1:3),"R_10"),
                                                      paste0("WT_", c(1:3),"R_100"),
                                                      paste0("H1_", c(1:3),"R_100"),
                                                      paste0("H2_", c(1:3),"R_100")))
xtest.m$ngroup    <- paste0(xtest.m$group,"_",xtest.m$Gene)

pdf(paste0(Out.dir,"/", Project, "-GenderGeneExpression.pdf"),width=30,height=12.5, onefile=T)
par(bg=NA)
ggplot(data=xtest.m, aes(x=Sample, y=log2(normCount+1), group=ngroup, shape=group, colour=Gene, label=Gene)) +
  geom_line() +
  geom_text_repel(aes(label=Gene), show.legend =FALSE, size=2)+
  geom_point(size=3) +
  theme_bw()
dev.off()

message("+--- Housekeeping genes checking -----------------+")
HgeneList  <- c("Piezo1","Ywhaz","Tbp","Actb","Hprt1", "Pgk1","Gusb","B2m", "Ppia")
    
vsd.assay  <- assay(vsd)
Hgindex    <- ensEMBL2id[ensEMBL2id$external_gene_name%in%HgeneList,]
vsd.selXY  <- vsd.assay[rownames(vsd.assay)%in%Hgindex[,1],]

hk.1            <- plotCounts(dds, intgroup=c("label"), normalized=TRUE, returnData=TRUE,gene="ENSRNOG00000056786")
colnames(hk.1)  <- c(HgeneList[1], "label")
hk.2            <- plotCounts(dds, intgroup=c("label"), normalized=TRUE, returnData=TRUE,gene="ENSRNOG00000008195")
colnames(hk.2)  <- c(HgeneList[2], "label")
hk.3            <- plotCounts(dds, intgroup=c("label"), normalized=TRUE, returnData=TRUE,gene="ENSRNOG00000001489")
hk.4            <- plotCounts(dds, intgroup=c("label"), normalized=TRUE, returnData=TRUE,gene="ENSRNOG00000034254")
hk.5            <- plotCounts(dds, intgroup=c("label"), normalized=TRUE, returnData=TRUE,gene="ENSRNOG00000031367")
hk.6            <- plotCounts(dds, intgroup=c("label"), normalized=TRUE, returnData=TRUE,gene="ENSRNOG00000058249")
hk.7            <- plotCounts(dds, intgroup=c("label"), normalized=TRUE, returnData=TRUE,gene="ENSRNOG00000000913")
hk.8            <- plotCounts(dds, intgroup=c("label"), normalized=TRUE, returnData=TRUE,gene="ENSRNOG00000017123")
hk.9            <- plotCounts(dds, intgroup=c("label"), normalized=TRUE, returnData=TRUE,gene="ENSRNOG00000027864")


xtest           <- merge(hk.1, hk.2, by=c("row.names", "label") )
colnames(xtest)[3:4] <- c(SgeneList[1:2])
xtest$Tbp       <- hk.3$count 
xtest$Actb      <- hk.4$count 
xtest$Hprt1     <- hk.5$count 
xtest$Pgk1      <- hk.6$count 
xtest$Gusb      <- hk.7$count 
xtest$B2m       <- hk.8$count 
xtest$Ppia      <- hk.9$count 

xtest           <- xtest[c(13,15,17,1,3,5,7,9,11,14,16,18,2,4,6,8,10,12),]

xtest.m           <- melt(xtest)
colnames(xtest.m) <- c("Sample", "group", "Gene", "normCount")
xtest.m$Sample    <- factor(xtest.m$Sample,levels = c(paste0("WT_", c(1:3),"R_10"),
                                                      paste0("H1_", c(1:3),"R_10"),
                                                      paste0("H2_", c(1:3),"R_10"),
                                                      paste0("WT_", c(1:3),"R_100"),
                                                      paste0("H1_", c(1:3),"R_100"),
                                                      paste0("H2_", c(1:3),"R_100")))
xtest.m$ngroup    <- paste0(xtest.m$group,"_",xtest.m$Gene)

pdf(paste0(Out.dir,"/", Project, "-HouseKeepingGeneExpression_rmR4.pdf"),width=30,height=12.5, onefile=T)
par(bg=NA)
ggplot(data=xtest.m, aes(x=Sample, y=log2(normCount+1), group=ngroup, shape=group, colour=Gene, label=Gene)) +
  geom_line() +
  geom_text_repel(aes(label=Gene), show.legend =FALSE, size=2)+
  geom_point(size=3) +
  theme_bw()
dev.off()


message("+WT special compare to cas9 KO------+")

DEGs.tlist   <- list()
for(i in c(1:9)){
  load(paste0(Out.dir, "/", Project, "-edgeR_comseq_batch_group_", models[i], "_rm", rmRep, "_DESeq_Object.RData"))
  DEGs.tlist[[i]]  <- subset(qlf.dat.ann, FDR < 0.05)
  DEGs.tlist
}
DEGs.tlist.len  <- unlist(lapply(DEGs.tlist, function(x) dim(x)[1]))
DEGs.tlist.comb <- rbind(DEGs.tlist[[1]], DEGs.tlist[[2]], DEGs.tlist[[3]],
                         DEGs.tlist[[4]], DEGs.tlist[[5]], DEGs.tlist[[6]],
                         DEGs.tlist[[7]], DEGs.tlist[[8]], DEGs.tlist[[9]]) ## 1718, unique 1046. 
DEGs.tlist.comb$model <- rep(models, DEGs.tlist.len)
## Fos and Npas4 duplicated, which common in H1_10vsH1_100 and H2_10vsH2_100

wgene.uni       <- DEGs.tlist[[9]]$external_gene_name
overD.list      <- list()
for( i in c(1:8)){
  load(paste0(Out.dir, "/", Project, "-edgeR_comseq_batch_group_", models[i], "_rm", rmRep, "_DESeq_Object.RData"))
  overD.list[[i]] <- qlf.dat.ann[qlf.dat.ann$external_gene_name%in%wgene.uni,]
  overD.list
  
}


