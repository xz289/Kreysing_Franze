#!/usr/local/bin/Rscript
# R 4.2.3
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
  #library("ggalt")
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
  library("openxlsx")
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
  library("cqn")
  library("BiocFileCache")
  ## save figures as svg format
  library("svglite")
})


register(MulticoreParam(2))


message("+-------------------------------------------------------------------------------+")
message("+            Set up some constants e.g. base directories                        +")
message("+-------------------------------------------------------------------------------+")

Project         <- "CTR_kf284_0004"
significance    <- 0.05
l2fc            <- 1 
elementTextSize <- 10

Base.dir <- "/Users/xz289/Documents/CTR_kf284_0004"
setwd(Base.dir)

Count.dir  <- "./count_files"
Out.dir <-"./Rm4R_Analysis_New_2023/DESeq2"


message("+-------------------------------------------------------------------------------+")
message("+ Load ensEMBL annotations with gene length and gc percentage                   +")
message("+-------------------------------------------------------------------------------+")

load("ensembl_gc_summary.RData")

## ensEMBL2id and uCovar
uCovar$ensembl_gene_id <- rownames(uCovar)
ensEMBL2id_mer <- merge(ensEMBL2id, uCovar, by = "ensembl_gene_id")

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

sampleSurface      <- unlist(lapply(sampleOri, function(x) strsplit(x, split="_")[[1]][2]))
sampleSurface      <- ifelse(sampleSurface=="10", "fast", "slow")

sampleCon        <- unlist(lapply(sampleOri, function(x) strsplit(x, split="_")[[1]][1]))
sampleGroup      <- unlist(lapply(sampleNames, function(x) strsplit(x, split="_")[[1]][1]))

sampleTable      <- data.frame(sampleName=sampleNames, 
                               fileName=sampleFiles,  
                               condition=factor(sampleCon, levels=c("WT", "H")),
                               replicate=sampleReps,
                               origin=factor(sampleOri, levels=c("WT_10", "WT_100",
                                                                 "H_10", "H_100")),
                               surface=factor(sampleSurface, levels=c("slow", "fast")),
                               treatment=factor(sampleGroup, levels=c("WT", "H1", "H2")),
                               group=factor(sampleLabels, levels=c("WT_10", "WT_100", 
                                                                   "H1_10", "H1_100",
                                                                   "H2_10", "H2_100"))
)

print(sampleTable)
nrow(sampleTable)
str(sampleTable)

sampleTable.rmR4 <- sampleTable[-grep("4R", sampleTable$sampleName), ] 
write.csv(sampleTable, file = paste0(Out.dir, "/SampleTable_Jan_2024.csv"))


message("+-------------------------------------------------------------------------------+")
message("+     Normalised raw counts using DESeq2 analysis pipeline                      +")
message("+-------------------------------------------------------------------------------+")


ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable=sampleTable.rmR4, 
                                       directory=Count.dir, 
                                       design=~0+group )
dds <- DESeq(ddsHTSeq, parallel=TRUE)
vsd <- vst(dds,     blind=F)
colData(vsd)
normC <- counts(dds, normalized=TRUE)
normC <- as.data.frame(normC)
normC$ensembl_gene_id <- rownames(normC)
normC.merRE <- merge(normC, ensEMBL2id, by = "ensembl_gene_id")
normC.merRE.log2 <- apply(normC.merRE[,2:19], 2, function(x) log2(x+1))
rownames(normC.merRE.log2) <- normC.merRE$ensembl_gene_id
normC.merRE.log2 <- as.data.frame(normC.merRE.log2)
normC.merRE.log2 <- normC.merRE.log2[,c("WT_1R_10", "WT_2R_10", "WT_3R_10","WT_1R_100", "WT_2R_100", "WT_3R_100",
                                        "H1_1R_10", "H1_2R_10", "H1_3R_10","H1_1R_100", "H1_2R_100", "H1_3R_100",
                                        "H2_1R_10", "H2_2R_10", "H2_3R_10","H2_1R_100", "H2_2R_100", "H2_3R_100")]
normC.merRE.log2F <- cbind(normC.merRE.log2, normC.merRE[,20:22])
write.csv(normC.merRE.log2F, file = paste0(Out.dir, "/Normalised_Counts_dds.csv"))

message("+-------------------------------------------------------------------------------+")
message("+-------Normalised individual genes box-dot plot                       ---------+")
message("+-------------------------------------------------------------------------------+")

makeGeneCountPlot_ord_new <- function(DDS, ensEMBL2id, CONDITION, gene2plot) {
  #
  # Plot the normalised read counts for a specified gene
  #
  
  genename2plot <- ensEMBL2id[ensEMBL2id$ensembl_gene_id == gene2plot, ]$external_gene_name
  t2            <- plotCounts(DDS, gene=gene2plot, intgroup=c(CONDITION), normalized=TRUE, returnData=TRUE)
  colnames(t2)  <- c("count", "condition")
  t2$count      <- log2(t2$count+1)
  
  t2$condition  <- factor(t2$condition,levels = c("WT_10", "WT_100", "H1_10", "H1_100", "H2_10", "H2_100"))
  
  
  plt.grp <- ggplot(t2, aes(x=condition, y=count, fill=condition)) + 
    geom_boxplot(width = 0.5, color=c("tomato","tomato", "steelblue",  "steelblue", "darkorange", "darkorange"), outlier.shape=NA, alpha=0.5) + 
    geom_point(aes(fill=condition),position=position_jitterdodge(),size=0.8) +
    scale_fill_manual(name="Group", values = c("tomato","tomato", "steelblue",  "steelblue", "darkorange", "darkorange")) +
    theme(text = element_text(size=elementTextSize), legend.position="none") +
    ggtitle(genename2plot) + 
    xlab("") + ylab("log2(Normalised count)") +
    theme_classic() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size=elementTextSize, face="bold"),
          axis.text.y = element_text(size=elementTextSize, face="bold"))

  return(plt.grp)

}

## Ripk3, Fcgr2b, Adam3a, and Pcsk9 are having low raw counts across the samples. 
genenames2plot   <- c("Optc", "Ttr", "Thbs1", "Gfap", "Ripk3", "Fcgr2b", "Adam3a", "Pcsk9")
for(i in 1:length(genenames2plot)){
  genes.WT.id     <- ensEMBL2id[ensEMBL2id$external_gene_name == genenames2plot[i], ]$ensembl_gene_id
  plt.count.ord   <- makeGeneCountPlot_ord_new(dds, ensEMBL2id, "group", genes.WT.id)
  genename2plot   <- ensEMBL2id[ensEMBL2id$ensembl_gene_id == genes.WT.id, ]$external_gene_name
  pdf(paste0(Out.dir, "/", Project, "-", genename2plot, "_Count_plot_DESeq2_new.pdf"), width=6, height=3 )
  print(plt.count.ord)
  dev.off()
  ggsave(file = paste0(Out.dir, "/", Project, "-", genename2plot, "_Count_plot_DESeq2_new.svg"), plot = plt.count.ord, width=6, height=3 )
}



message("+-------------------------------------------------------------------------------+")
message("+--------               Customised PCA                                 ---------+")
message("+-------------------------------------------------------------------------------+")

TOPNUM <- 500
pcaData    <- plotPCA(vsd, ntop=TOPNUM, intgroup=c("group"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf(paste0(Out.dir, "/", Project, "_rm4R_Fig.PCA.T", TOPNUM, "_DESeq2_Jan_2023.pdf"), width = 5, height = 6)
par(bg=NA)
plt1 <- ggplot(pcaData, aes(PC1, PC2, color=group)) +
  geom_mark_ellipse(aes(fill = NULL, color=group, group=group, 
                        label=group), alpha=0.1, label.fontsize = 6, label.buffer = unit(15, 'mm')) +
  geom_point(size=2, alpha=0.75) +
  geom_text_repel(aes(label=name), show.legend = FALSE, size=2)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  labs(title = paste0("PCA (Top ", TOPNUM, " Most Variable Genes, all samples)")) +
  theme(text = element_text(size=12), 
        plot.subtitle=element_text(size=8,  face="italic", color="darkgrey"),
        legend.position='top', aspect.ratio=1) 


print(plt1)
dev.off()

ggsave(file = paste0(Out.dir, "/", Project, "-PCAplot_deseq2_new_1.svg"), plot = plt1, width=5, height=6)

message("+-------------------------------------------------------------------------------+")
message("+-------pairwise DESeq2 analysis for selected possible pairs           ---------+")
message("+-------------------------------------------------------------------------------+")
## useful function for customised paired pca
customPCA_pair     <- function(sampleTBL, RLD, TOPNUM) {
  
  RLDs   <- assay(RLD)
  rv     <- rowVars(RLDs, useNames = T)
  select <- order(rv, decreasing = TRUE)[seq_len(min(TOPNUM, length(rv)))]
  pca    <- prcomp(t(RLDs[select, ]))
  
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  
  scores    <- data.frame(sampleName=sampleTBL$sampleName, pca$x, 
                          condition=sampleTBL$group)
  
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
    geom_mark_ellipse(alpha = 0.1, show.legend = FALSE, aes(colour=condition, fill=condition, group = condition)) +
    xlab(pc1lab) + ylab(pc2lab) + 
    scale_shape_manual(name="Group", values = c(17, 16)) +
    ggtitle(paste0("PCA Top ", TOPNUM, " MV")) +
    theme(text = element_text(size=elementTextSize))+ 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  
  return(list(plt.pca, plt.pca.nl) )
  
}

## define paired
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

for(i in 1:length(control)){
  print(i)
  sampleTable.Hs   <- sampleTable.rmR4[sampleTable.rmR4$group%in%c(control[i], treatment[i]),]
  group            <- factor(sampleTable.Hs$group, levels = c(control[i], treatment[i]))
  count_matrix.Hs  <- assay(ddsHTSeq)[,pairColum[[i]]]
  sampleTable.Hs$group <- factor(sampleTable.Hs$group, levels = c(control[i], treatment[i]))
  ddsHTSeq.Hs      <- DESeqDataSetFromMatrix(countData = count_matrix.Hs,
                                             colData = sampleTable.Hs,
                                             design = ~group)
  ddsHTSeq.Hs      <- estimateSizeFactors(ddsHTSeq.Hs)
  ddsHTSeq.Hs      <- DESeq(ddsHTSeq.Hs, parallel=TRUE)
  vsd.Hs           <- vst(ddsHTSeq.Hs, blind=FALSE)
  
  pdf(paste0(Out.dir, "/", Project, "-PCAplot_", models[i], "_DESeq2_Top", TOPNUM, "_rm", rmRep,"_Jan_2024.pdf"), width=6, height=5)
  pca.plot         <- customPCA_pair(sampleTable.Hs, vsd.Hs, TOPNUM)
  print(pca.plot[[1]])
  dev.off()
  
  ggsave(file = paste0(Out.dir, "/", Project, "-PCAplot_", models[i], "_DESeq2_new.svg"), plot = pca.plot[[1]], width=6, height=5)
  
  ## results DEGs
  res <- results(ddsHTSeq.Hs, name = resultsNames(ddsHTSeq.Hs)[2])
  res_dat <- as.data.frame(res)
  res_dat$ensembl_gene_id <- rownames(res_dat)
  res_dat_merE <- merge(res_dat, ensEMBL2id, by = "ensembl_gene_id")
  res_dat_merE <- res_dat_merE[!is.na(res_dat_merE$padj),]
  openxlsx::write.xlsx(res_dat_merE, file = paste0(Out.dir, "/", Project, "-PCAplot_", models[i], "_DESeq2_results_rm", rmRep,"_Jan_2024.xlsx"))
  gc()
  
}

message("+-------------------------------------------------------------------------------+")
message("+-----volcano plot generation, based on either pvalue or padj for yaxis---------+")
message("+-------------------------------------------------------------------------------+")

Resfiles <- paste0(Out.dir, "/", Project, "-PCAplot_", models, "_DESeq2_results_rm", rmRep,"_Jan_2024.xlsx")

ResfileList <- lapply(Resfiles, function(x) openxlsx::read.xlsx(x, sheet = 1))

unlist(lapply(ResfileList, function(x) sum(x$padj <= 0.05)))
## [1]  0  0  0  0 59 13 31  4  0

## models 5,6,7,8 will use padj as yaxis, the rest of the model will apply pvalue.

functionPlotDEVolcano <- function(results, sig_cut, logfc_cut, title, topN, xlabel, ylabel) {
  
  results       <- as.data.frame(results)
  results$genes <- results$external_gene_name
  padjcount     <- sum(results$padj<=0.05)

  if(padjcount!=0){
    pvalnew <- results$padj
    }else{
      pvalnew <- results$pvalue
    }
  results$pval <- pvalnew
  
  xrange_sum <- round(max(abs(results$log2FoldChange))) + 1
  xrange <- c(-xrange_sum, xrange_sum, xrange_sum/4)
  yrange_sum <- round(max(-log10(results$pval))) + 1
  yrange <- c(0, yrange_sum,  5)
  
  volc.plt <- ggplot(data=results, aes(x=log2FoldChange, y=-log10(pval), label=genes)) +
    geom_vline(xintercept = logfc_cut,     colour="black", linetype = "dashed", alpha=0.5) +
    geom_vline(xintercept = -(logfc_cut),  colour="black", linetype = "dashed", alpha=0.5) +
    geom_hline(yintercept = -log10(sig_cut), colour="black", linetype = "dashed", alpha=0.5) +
    
    geom_point(data=subset(results, abs(log2FoldChange) < logfc_cut | pval > sig_cut), alpha=0.75, size=0.3, colour="grey") +
    geom_point(data=subset(results, pval<=sig_cut & log2FoldChange >= logfc_cut),      alpha=0.75, size=0.8, colour="red") +
    geom_point(data=subset(results, pval<=sig_cut & log2FoldChange <= -(logfc_cut)),   alpha=0.75, size=0.8, colour="blue") +
    geom_text_repel( data= subset(results, log2FoldChange > logfc_cut & pval<= sig_cut & external_gene_name!="")[1:topN,],
                     show.legend = FALSE, nudge_x=0.1, nudge_y=0.1, segment.size = 0.25, size=3 ) +
    geom_text_repel( data= tail(subset(results, log2FoldChange< (-logfc_cut) & pval<= sig_cut & external_gene_name!=""), n=topN),
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

xlabels   <- c("log2FC (H1_10/H1_100)", "log2FC (H2_10/H2_100)", "log2FC (H2_10/H1_10)","log2FC (H2_100/H1_100)",
               "log2FC (H1_10/WT_10)", "log2FC (H2_10/WT_10)", "log2FC (H1_100/WT_100)","log2FC (H2_100/WT_100)",
               "log2FC (WT_10/WT_100)")
ylabels    <- c("-log10(pval)","-log10(pval)","-log10(pval)","-log10(pval)",
                "-log10(padj)", "-log10(padj)","-log10(padj)","-log10(padj)",
                "-log10(pval)")
for(i in 1:length(models)){
  print(i)
  volplt <- functionPlotDEVolcano(ResfileList[[i]], sig_cut=0.05, logfc_cut=0.6, title=models[i], topN=20, 
                                  xlabel=xlabels[i], ylabel=ylabels[i])
  pdf(paste0(Out.dir, "/", Project, "-Volcanoplot_", models[i], "_DESeq2_results_rm", rmRep,"_Jan_2024.pdf"))
  print(volplt)
  dev.off()
  ggsave(file = paste0(Out.dir, "/", Project, "-Volcanoplot_", models[i], "_DESeq2_new.svg"), plot = volplt, width=6, height=6)
}

## Finish 03/01/2024 ##