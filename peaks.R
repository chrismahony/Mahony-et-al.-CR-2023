library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(EnsDb.Mmusculus.v79)
#library(EnsDb.Mmusculus.v75)
library(BSgenome.Mmusculus.UCSC.mm10)
library(ggplot2)
library(cowplot)
library(patchwork) #latest version is required!
library(TFBSTools)
library(JASPAR2020)
library(gsfisher)
library(harmony)
library(clustree)
library(plyr)
library(dplyr)
library(stringr)


setwd("/rds/projects/2018/monteirr-lab/2019-12-19.SingleCellRNAseq/analysis_from_loupe_CM")
load("/rds/projects/2018/monteirr-lab/2019-12-19.SingleCellRNAseq/analysis_from_loupe_CM/moncle_ATACseq.RData")

aggr$CMclusters <- aggr@active.ident
Idents(aggr) <- "CMclusters"
aggr$cluster_sample <- paste(Idents(aggr), aggr$sample, sep = "_")

Idents(aggr) <- "cluster_sample"
peaks_ery <- FindMarkers(aggr, ident.1 = "ery_MUT", ident.2 = "ery_WT", test.use = 'LR', latent.vars = 'nCount_peaks', logfc.threshold= -Inf, min.pct = -Inf, only.pos = FALSE)
peaks_eryprog <- FindMarkers(aggr, ident.1 = "eryProg_MUT", ident.2 = "eryProg_WT", test.use = 'LR', latent.vars = 'nCount_peaks', logfc.threshold= -Inf, min.pct = -Inf, only.pos = FALSE)
peaks_dHSPCs <- FindMarkers(aggr, ident.1 = "dHSPCs_MUT", ident.2 = "dHSPCs_WT", test.use = 'LR', latent.vars = 'nCount_peaks', logfc.threshold= -Inf, min.pct = -Inf, only.pos = FALSE)
peaks_DN2 <- FindMarkers(aggr, ident.1 = "DN2_MUT", ident.2 = "DN2_WT", test.use = 'LR', latent.vars = 'nCount_peaks', logfc.threshold= -Inf, min.pct = -Inf, only.pos = FALSE)

peaks_ery <- peaks_ery %>%  arrange(desc(avg_log2FC))
output <- peaks_ery 
output <- str_split_fixed(rownames(output), "-", 3)
output <- as.data.frame(output)
output$avg_log2FC <- peaks_ery$avg_log2FC
output$avg_log2FC = '^'(2, output$avg_log2FC)
output <- output[!output$avg_log2FC == "1", ]
write.table(output,"./peaks_eryBDG.txt",sep="\t",row.names=FALSE, col.names=FALSE, quote = F)
peaks_ery <- output
output = subset(output, select = -c(avg_log2FC))
write.table(output,"./peaks_eryBED.txt",sep="\t",row.names=FALSE, col.names=FALSE, quote = F)
rm(output)

peaks_eryprog  <- peaks_eryprog  %>%  arrange(desc(avg_log2FC))
output <- peaks_eryprog 
output <- str_split_fixed(rownames(output), "-", 3)
output <- as.data.frame(output)
output$avg_log2FC <- peaks_eryprog$avg_log2FC
output$avg_log2FC = '^'(2, output$avg_log2FC)
output <- output[!output$avg_log2FC == "1", ]
write.table(output,"./peaks_eryprogBDG.txt",sep="\t",row.names=FALSE, col.names=FALSE, quote = F)
peaks_eryprog <- output
output = subset(output, select = -c(avg_log2FC))
write.table(output,"./peaks_eryprogBED.txt",sep="\t",row.names=FALSE, col.names=FALSE, quote = F)
rm(output)

peaks_dHSPCs <- peaks_dHSPCs %>%  arrange(desc(avg_log2FC))
output <- peaks_dHSPCs 
output <- str_split_fixed(rownames(output), "-", 3)
output <- as.data.frame(output)
output$avg_log2FC <- peaks_dHSPCs$avg_log2FC
output$avg_log2FC = '^'(2, output$avg_log2FC)
output <- output[!output$avg_log2FC == "1", ]
write.table(output,"./peaks_dHSPCsBDG.txt",sep="\t",row.names=FALSE, col.names=FALSE, quote = F)
peaks_dHSPCs <- output
output = subset(output, select = -c(avg_log2FC))
write.table(output,"./peaks_dHSPCsBED.txt",sep="\t",row.names=FALSE, col.names=FALSE, quote = F)
rm(output)

peaks_DN2 <- peaks_DN2 %>%  arrange(desc(avg_log2FC))
output <- peaks_DN2 
output <- str_split_fixed(rownames(output), "-", 3)
output <- as.data.frame(output)
output$avg_log2FC <- peaks_DN2$avg_log2FC
output$avg_log2FC = '^'(2, output$avg_log2FC)
output <- output[!output$avg_log2FC == "1", ]
write.table(output,"./peaks_DN2BDG.txt",sep="\t",row.names=FALSE, col.names=FALSE, quote = F)
peaks_DN2 <- output
output = subset(output, select = -c(avg_log2FC))
write.table(output,"./peaks_DN2BED.txt",sep="\t",row.names=FALSE, col.names=FALSE, quote = F)
rm(output)

save.image("moncle_ATACseq_NEW.RData")


