##########################################################################################
##########################################################################################
##
## Intron retention -- FeatureCounts - HPV with introns - HUMAN
##
##########################################################################################
##########################################################################################

genes <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/HPV_counts_introns.mx",sep="\t", skip=1, header = T, row.names=1) 
colnames(genes) <- sub("X.home.mmacchie.human.encode.infection_TEtrans.HPV.bams_no_virus_UCSC.", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.sorted_DS.bam", "", colnames(genes))
colnames(genes) <- sub("...bams_no_virus_UCSC.", "", colnames(genes))

hg38.viruses <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/virus.genes.hg38.txt",sep="\t",header=FALSE)

# separate counts columns from annotation columns
genes <- genes[!rownames(genes) %in% hg38.viruses$V1,]
gene.anno <- genes[, c(1:5)]
genes <- genes[, c(6:dim(genes)[2])]

library(edgeR)
y <- DGEList(counts=genes,group=c(1,1,2,2))
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

DEgenes <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) >1 ,]
DEintrons <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 1 ,]
DEintrons_up <- DEintrons[DEintrons$logFC < 0,]
DEintrons_down <- DEintrons[DEintrons$logFC > 0,]

DEgenes2 <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes2 <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) > 2,]
DEintrons2 <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons2 <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 2 ,] 

dim(DEgenes)
dim(DEintrons)
dim(DEintrons_up)
dim(DEintrons_down)

write.table(DEintrons,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_HPV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEintrons2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_HPV_FDR0.5FC2.txt",sep="\t",quote=FALSE)
write.table(DEgenes,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_HPV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEgenes2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_HPV_FDR0.5FC2.txt",sep="\t",quote=FALSE)



gcs <- subset(genes,!grepl("_intron$",rownames(genes)))
gcs <- colSums(gcs)/colSums(genes)
ics <- subset(genes,grepl("_intron$",rownames(genes)))
ics <- colSums(ics)/colSums(genes)
data.frame("HPV",gcs,ics)



##########################################################################################
##########################################################################################
##
## Intron retention -- FeatureCounts - HCV with introns - HUMAN
##
##########################################################################################
##########################################################################################

genes <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/HCV_counts_introns.mx",sep="\t", skip=1, header = T, row.names=1) 
colnames(genes) <- sub("X.home.mmacchie.human.encode.infection_TEtrans.HCV.bams_no_virus_UCSC.", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.sorted_DS.bam", "", colnames(genes))
colnames(genes) <- sub("...bams_no_virus_UCSC.", "", colnames(genes))

hg38.viruses <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/virus.genes.hg38.txt",sep="\t",header=FALSE)

# separate counts columns from annotation columns
genes <- genes[!rownames(genes) %in% hg38.viruses$V1,]
gene.anno <- genes[, c(1:5)]
genes <- genes[, c(6:dim(genes)[2])]

library(edgeR)
y <- DGEList(counts=genes,group=c(1,1,1,1,1,2,2,2,2,2))
keep <- rowSums(cpm(y)>1) >= 5
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

DEgenes <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) >1 ,]
DEintrons <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 1 ,] 
DEintrons_up <- DEintrons[DEintrons$logFC < 0,]
DEintrons_down <- DEintrons[DEintrons$logFC > 0,]

DEgenes2 <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes2 <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) > 2,]
DEintrons2 <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons2 <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 2 ,] 

dim(DEgenes)
dim(DEintrons)
dim(DEintrons_up)
dim(DEintrons_down)


write.table(DEintrons,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_HCV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEintrons2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_HCV_FDR0.5FC2.txt",sep="\t",quote=FALSE)
write.table(DEgenes,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_HCV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEgenes2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_HCV_FDR0.5FC2.txt",sep="\t",quote=FALSE)



gcs <- subset(genes,!grepl("_intron$",rownames(genes)))
gcs <- colSums(gcs)/colSums(genes)
ics <- subset(genes,grepl("_intron$",rownames(genes)))
ics <- colSums(ics)/colSums(genes)
data.frame("HCV",gcs,ics)



##########################################################################################
##########################################################################################
##
## Intron retention -- FeatureCounts - HIVrest with introns - HUMAN
##
##########################################################################################
##########################################################################################

genes <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/HIV_counts_introns.mx",sep="\t", skip=1, header = T, row.names=1) 
colnames(genes) <- sub("X.home.mmacchie.human.encode.infection_TEtrans.HIV.bams_no_virus_UCSC.", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.sorted_DS.bam", "", colnames(genes))
colnames(genes) <- sub("...bams_no_virus_UCSC.", "", colnames(genes))

hg38.viruses <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/virus.genes.hg38.txt",sep="\t",header=FALSE)

# separate counts columns from annotation columns
genes <- genes[!rownames(genes) %in% hg38.viruses$V1,]
gene.anno <- genes[, c(1:5)]
genes <- genes[, c(6:dim(genes)[2])]

genes <- genes[,c(4,8,12,10,2,6)]
library(edgeR)
y <- DGEList(counts=genes,group=c(1,1,1,2,2,2))
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

DEgenes <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) >1 ,]
DEintrons <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 1 ,] 
DEintrons_up <- DEintrons[DEintrons$logFC < 0,]
DEintrons_down <- DEintrons[DEintrons$logFC > 0,]

DEgenes2 <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes2 <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) > 2,]
DEintrons2 <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons2 <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 2 ,] 

dim(DEgenes)
dim(DEintrons)
dim(DEintrons_up)
dim(DEintrons_down)


write.table(DEintrons,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_HIVrest_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEintrons2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_HIVrest_FDR0.5FC2.txt",sep="\t",quote=FALSE)
write.table(DEgenes,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_HIVrest_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEgenes2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_HIVrest_FDR0.5FC2.txt",sep="\t",quote=FALSE)


gcs <- subset(genes,!grepl("_intron$",rownames(genes)))
gcs <- colSums(gcs)/colSums(genes)
ics <- subset(genes,grepl("_intron$",rownames(genes)))
ics <- colSums(ics)/colSums(genes)
data.frame("HIVrest",gcs,ics)


##########################################################################################
##########################################################################################
##
## Intron retention -- FeatureCounts - HIVact with introns - HUMAN
##
##########################################################################################
##########################################################################################

genes <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/HIV_counts_introns.mx",sep="\t", skip=1, header = T, row.names=1) 
colnames(genes) <- sub("X.home.mmacchie.human.encode.infection_TEtrans.HIV.bams_no_virus_UCSC.", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.sorted_DS.bam", "", colnames(genes))
colnames(genes) <- sub("...bams_no_virus_UCSC.", "", colnames(genes))

hg38.viruses <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/virus.genes.hg38.txt",sep="\t",header=FALSE)

# separate counts columns from annotation columns
genes <- genes[!rownames(genes) %in% hg38.viruses$V1,]
gene.anno <- genes[, c(1:5)]
genes <- genes[, c(6:dim(genes)[2])]

genes <- genes[,c(3,7,11,1,5,9)]
library(edgeR)
y <- DGEList(counts=genes,group=c(1,1,1,2,2,2))
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

DEgenes <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) >1 ,]
DEintrons <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 1 ,] 
DEintrons_up <- DEintrons[DEintrons$logFC < 0,]
DEintrons_down <- DEintrons[DEintrons$logFC > 0,]

DEgenes2 <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes2 <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) > 2,]
DEintrons2 <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons2 <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 2 ,] 

dim(DEgenes)
dim(DEintrons)
dim(DEintrons_up)
dim(DEintrons_down)


write.table(DEintrons,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_HIVact_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEintrons2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_HIVact_FDR0.5FC2.txt",sep="\t",quote=FALSE)
write.table(DEgenes,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_HIVact_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEgenes2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_HIVact_FDR0.5FC2.txt",sep="\t",quote=FALSE)



gcs <- subset(genes,!grepl("_intron$",rownames(genes)))
gcs <- colSums(gcs)/colSums(genes)
ics <- subset(genes,grepl("_intron$",rownames(genes)))
ics <- colSums(ics)/colSums(genes)
data.frame("HIVact",gcs,ics)


##########################################################################################
##########################################################################################
##
## Intron retention -- FeatureCounts - SeV with introns - HUMAN
##
##########################################################################################
##########################################################################################

genes <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/SeV_counts_introns.mx",sep="\t", skip=1, header = T, row.names=1) 
colnames(genes) <- sub("X.home.mmacchie.human.encode.infection_TEtrans.sendai.bams_no_virus_UCSC.", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.sorted_DS.bam", "", colnames(genes))
colnames(genes) <- sub("...bams_no_virus_UCSC.", "", colnames(genes))

hg38.viruses <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/virus.genes.hg38.txt",sep="\t",header=FALSE)

# separate counts columns from annotation columns
genes <- genes[!rownames(genes) %in% hg38.viruses$V1,]
gene.anno <- genes[, c(1:5)]
genes <- genes[, c(6:dim(genes)[2])]
genes <- genes[,c(3,4,1,2)]

library(edgeR)
y <- DGEList(counts=genes,group=c(1,1,2,2))
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

DEgenes <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) >1 ,]
DEintrons <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 1 ,] 
DEintrons_up <- DEintrons[DEintrons$logFC < 0,]
DEintrons_down <- DEintrons[DEintrons$logFC > 0,]

DEgenes2 <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes2 <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) > 2,]
DEintrons2 <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons2 <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 2 ,] 

dim(DEgenes)
dim(DEintrons)
dim(DEintrons_up)
dim(DEintrons_down)


write.table(DEintrons,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_SeV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEintrons2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_SeV_FDR0.5FC2.txt",sep="\t",quote=FALSE)
write.table(DEgenes,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_SeV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEgenes2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_SeV_FDR0.5FC2.txt",sep="\t",quote=FALSE)


gcs <- subset(genes,!grepl("_intron$",rownames(genes)))
gcs <- colSums(gcs)/colSums(genes)
ics <- subset(genes,grepl("_intron$",rownames(genes)))
ics <- colSums(ics)/colSums(genes)
data.frame("SeV",gcs,ics)


##########################################################################################
##########################################################################################
##
## Intron retention -- FeatureCounts - KSHV in HUVEC cells - introns - HUMAN
##
##########################################################################################
##########################################################################################

genes <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/KSHV_counts_introns.mx",sep="\t", skip=1, header = T, row.names=1) 
colnames(genes) <- sub("X.home.mmacchie.human.encode.infection_TEtrans.KSHV.bams_no_virus_UCSC.", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.sorted_DS.bam", "", colnames(genes))
colnames(genes) <- sub("...bams_no_virus_UCSC.", "", colnames(genes))

hg38.viruses <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/virus.genes.hg38.txt",sep="\t",header=FALSE)

# separate counts columns from annotation columns
genes <- genes[!rownames(genes) %in% hg38.viruses$V1,]
gene.anno <- genes[, c(1:5)]
genes <- genes[, c(6:dim(genes)[2])]

genes <- genes[,c(1:4)]
library(edgeR)
y <- DGEList(counts=genes,group=c(1,1,2,2))
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

DEgenes <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) >1 ,]
DEintrons <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 1 ,] 
DEintrons_up <- DEintrons[DEintrons$logFC < 0,]
DEintrons_down <- DEintrons[DEintrons$logFC > 0,]

DEgenes2 <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes2 <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) > 2,]
DEintrons2 <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons2 <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 2 ,] 

dim(DEgenes)
dim(DEintrons)
dim(DEintrons_up)
dim(DEintrons_down)


write.table(DEintrons,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_KSHV_HUVEC_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEintrons2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_KSHV_HUVEC_FDR0.5FC2.txt",sep="\t",quote=FALSE)
write.table(DEgenes,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_KSHV_HUVEC_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEgenes2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_KSHV_HUVEC_FDR0.5FC2.txt",sep="\t",quote=FALSE)


gcs <- subset(genes,!grepl("_intron$",rownames(genes)))
gcs <- colSums(gcs)/colSums(genes)
ics <- subset(genes,grepl("_intron$",rownames(genes)))
ics <- colSums(ics)/colSums(genes)
data.frame("KSHVhuvec",gcs,ics)


##########################################################################################
##########################################################################################
##
## Intron retention -- FeatureCounts - KSHV in MC116 cells - introns - HUMAN
##
##########################################################################################
##########################################################################################

genes <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/KSHV_counts_introns.mx",sep="\t", skip=1, header = T, row.names=1) 
colnames(genes) <- sub("X.home.mmacchie.human.encode.infection_TEtrans.KSHV.bams_no_virus_UCSC.", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.sorted_DS.bam", "", colnames(genes))
colnames(genes) <- sub("...bams_no_virus_UCSC.", "", colnames(genes))

hg38.viruses <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/virus.genes.hg38.txt",sep="\t",header=FALSE)

# separate counts columns from annotation columns
genes <- genes[!rownames(genes) %in% hg38.viruses$V1,]
gene.anno <- genes[, c(1:5)]
genes <- genes[, c(6:dim(genes)[2])]

genes <- genes[,c(5:10)]
library(edgeR)
y <- DGEList(counts=genes,group=c(1,1,1,2,2,2))
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

DEgenes <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) >1 ,]
DEintrons <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 1 ,] 
DEintrons_up <- DEintrons[DEintrons$logFC < 0,]
DEintrons_down <- DEintrons[DEintrons$logFC > 0,]

DEgenes2 <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes2 <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) > 2,]
DEintrons2 <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons2 <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 2 ,] 

dim(DEgenes)
dim(DEintrons)
dim(DEintrons_up)
dim(DEintrons_down)


write.table(DEintrons,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_KSHV_MC116_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEintrons2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_KSHV_MC116_FDR0.5FC2.txt",sep="\t",quote=FALSE)
write.table(DEgenes,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_KSHV_MC116_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEgenes2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_KSHV_MC116_FDR0.5FC2.txt",sep="\t",quote=FALSE)

gcs <- subset(genes,!grepl("_intron$",rownames(genes)))
gcs <- colSums(gcs)/colSums(genes)
ics <- subset(genes,grepl("_intron$",rownames(genes)))
ics <- colSums(ics)/colSums(genes)
data.frame("KSHVmc116",gcs,ics)




##########################################################################################
##########################################################################################
##
## Intron retention -- FeatureCounts - ZIKV with introns - HUMAN
##
##########################################################################################
##########################################################################################

genes <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/ZIKV_counts_introns.mx",sep="\t", skip=1, header = T, row.names=1) 
colnames(genes) <- sub("X.home.mmacchie.human.encode.infection_TEtrans.ZIKV.bams_no_virus_UCSC.", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.sorted_DS.bam", "", colnames(genes))
colnames(genes) <- sub("...bams_no_virus_UCSC.", "", colnames(genes))

hg38.viruses <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/virus.genes.hg38.txt",sep="\t",header=FALSE)

# separate counts columns from annotation columns
genes <- genes[!rownames(genes) %in% hg38.viruses$V1,]
gene.anno <- genes[, c(1:5)]
genes <- genes[, c(6:dim(genes)[2])]


library(edgeR)
y <- DGEList(counts=genes,group=c(1,1,1,2,2,2))
keep <- rowSums(cpm(y)>1) >= 3
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

DEgenes <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) >1 ,]
DEintrons <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 1 ,] 
DEintrons_up <- DEintrons[DEintrons$logFC < 0,]
DEintrons_down <- DEintrons[DEintrons$logFC > 0,]

DEgenes2 <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes2 <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) > 2,]
DEintrons2 <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons2 <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 2 ,] 

dim(DEgenes)
dim(DEintrons)
dim(DEintrons_up)
dim(DEintrons_down)


write.table(DEintrons,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_ZIKV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEintrons2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_ZIKV_FDR0.5FC2.txt",sep="\t",quote=FALSE)
write.table(DEgenes,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_ZIKV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEgenes2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_ZIKV_FDR0.5FC2.txt",sep="\t",quote=FALSE)


gcs <- subset(genes,!grepl("_intron$",rownames(genes)))
gcs <- colSums(gcs)/colSums(genes)
ics <- subset(genes,grepl("_intron$",rownames(genes)))
ics <- colSums(ics)/colSums(genes)
data.frame("ZIKV",gcs,ics)


##########################################################################################
##########################################################################################
##
## Intron retention -- FeatureCounts - EBOV with introns - HUMAN
##
##########################################################################################
##########################################################################################

genes <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/EBOV_counts_introns.mx",sep="\t", skip=1, header = T, row.names=1) 
colnames(genes) <- sub("X.home.mmacchie.human.encode.infection_TEtrans.EBOV.bams_no_virus_UCSC.", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.sorted_DS.bam", "", colnames(genes))
colnames(genes) <- sub("...bams_no_virus_UCSC.", "", colnames(genes))

hg38.viruses <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/virus.genes.hg38.txt",sep="\t",header=FALSE)

# separate counts columns from annotation columns
genes <- genes[!rownames(genes) %in% hg38.viruses$V1,]
gene.anno <- genes[, c(1:5)]
genes <- genes[, c(6:dim(genes)[2])]


library(edgeR)
y <- DGEList(counts=genes,group=c(1,1,1,2,2,2))
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

DEgenes <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) >1 ,]
DEintrons <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 1 ,] 
DEintrons_up <- DEintrons[DEintrons$logFC < 0,]
DEintrons_down <- DEintrons[DEintrons$logFC > 0,]

DEgenes2 <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes2 <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) > 2,]
DEintrons2 <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons2 <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 2 ,] 

dim(DEgenes)
dim(DEintrons)
dim(DEintrons_up)
dim(DEintrons_down)


write.table(DEintrons,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_EBOV_ARPE19_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEintrons2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_EBOV_ARPE19_FDR0.5FC2.txt",sep="\t",quote=FALSE)
write.table(DEgenes,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_EBOV_ARPE19_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEgenes2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_EBOV_ARPE19_FDR0.5FC2.txt",sep="\t",quote=FALSE)

gcs <- subset(genes,!grepl("_intron$",rownames(genes)))
gcs <- colSums(gcs)/colSums(genes)
ics <- subset(genes,grepl("_intron$",rownames(genes)))
ics <- colSums(ics)/colSums(genes)
data.frame("EBOV",gcs,ics)


##########################################################################################
##########################################################################################
##
## Intron retention -- FeatureCounts - EBOV in MDM with introns - HUMAN
##
##########################################################################################
##########################################################################################

genes <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/EBOV_MDM_counts_introns.mx",sep="\t", skip=1, header = T, row.names=1) 
colnames(genes) <- sub("X.home.mmacchie.human.encode.infection_TEtrans.EBOV_MDM.bams_no_virus_UCSC.", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.sorted_DS.bam", "", colnames(genes))
colnames(genes) <- sub("...bams_no_virus_UCSC.", "", colnames(genes))

hg38.viruses <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/virus.genes.hg38.txt",sep="\t",header=FALSE)

# separate counts columns from annotation columns
genes <- genes[!rownames(genes) %in% hg38.viruses$V1,]
gene.anno <- genes[, c(1:5)]
genes <- genes[, c(6:dim(genes)[2])]


library(edgeR)
y <- DGEList(counts=genes,group=c(1,1,1,2,2,2))
keep <- rowSums(cpm(y)>1) >= 3
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

DEgenes <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) >1 ,]
DEintrons <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 1 ,] 
DEintrons_up <- DEintrons[DEintrons$logFC < 0,]
DEintrons_down <- DEintrons[DEintrons$logFC > 0,]

DEgenes2 <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes2 <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) > 2,]
DEintrons2 <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons2 <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 2 ,] 

dim(DEgenes)
dim(DEintrons)
dim(DEintrons_up)
dim(DEintrons_down)


write.table(DEintrons,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_EBOV_MDM_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEintrons2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_EBOV_MDM_FDR0.5FC2.txt",sep="\t",quote=FALSE)
write.table(DEgenes,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_EBOV_MDM_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEgenes2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_EBOV_MDM_FDR0.5FC2.txt",sep="\t",quote=FALSE)

gcs <- subset(genes,!grepl("_intron$",rownames(genes)))
gcs <- colSums(gcs)/colSums(genes)
ics <- subset(genes,grepl("_intron$",rownames(genes)))
ics <- colSums(ics)/colSums(genes)
data.frame("EBOV_MDM",gcs,ics)


##########################################################################################
##########################################################################################
##
## Intron retention -- FeatureCounts - EBOV in CD4Tcells with introns - HUMAN
##
##########################################################################################
##########################################################################################

genes <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/EBOV_CD4Tcells_counts_introns.mx",sep="\t", skip=1, header = T, row.names=1) 
colnames(genes) <- sub("X.home.mmacchie.human.encode.infection_TEtrans.EBOV_CD4Tcells.bams_no_virus_UCSC.", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.sorted_DS.bam", "", colnames(genes))
colnames(genes) <- sub("...bams_no_virus_UCSC.", "", colnames(genes))

hg38.viruses <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/virus.genes.hg38.txt",sep="\t",header=FALSE)

# separate counts columns from annotation columns
genes <- genes[!rownames(genes) %in% hg38.viruses$V1,]
gene.anno <- genes[, c(1:5)]
genes <- genes[, c(6:dim(genes)[2])]


library(edgeR)
y <- DGEList(counts=genes,group=c(1,1,1,1,2,2,2,2))
keep <- rowSums(cpm(y)>1) >= 4
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

DEgenes <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) >1 ,]
DEintrons <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 1 ,] 
DEintrons_up <- DEintrons[DEintrons$logFC < 0,]
DEintrons_down <- DEintrons[DEintrons$logFC > 0,]

DEgenes2 <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes2 <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) > 2,]
DEintrons2 <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons2 <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 2 ,] 

dim(DEgenes)
dim(DEintrons)
dim(DEintrons_up)
dim(DEintrons_down)


write.table(DEintrons,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_EBOV_CD4Tcell_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEintrons2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_EBOV_CD4Tcell_FDR0.5FC2.txt",sep="\t",quote=FALSE)
write.table(DEgenes,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_EBOV_CD4Tcell_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEgenes2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_EBOV_CD4Tcell_FDR0.5FC2.txt",sep="\t",quote=FALSE)


gcs <- subset(genes,!grepl("_intron$",rownames(genes)))
gcs <- colSums(gcs)/colSums(genes)
ics <- subset(genes,grepl("_intron$",rownames(genes)))
ics <- colSums(ics)/colSums(genes)
data.frame("EBOVtcell",gcs,ics)




##########################################################################################
##########################################################################################
##
## Intron retention -- FeatureCounts - HRSV with introns - HUMAN
##
##########################################################################################
##########################################################################################

genes <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/HRSV_counts_introns.mx",sep="\t", skip=1, header = T, row.names=1) 
colnames(genes) <- sub("X.home.mmacchie.human.encode.infection_TEtrans.HRSV.bams_no_virus_UCSC.", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.sorted_DS.bam", "", colnames(genes))
colnames(genes) <- sub("...bams_no_virus_UCSC.", "", colnames(genes))

hg38.viruses <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/virus.genes.hg38.txt",sep="\t",header=FALSE)

# separate counts columns from annotation columns
genes <- genes[!rownames(genes) %in% hg38.viruses$V1,]
gene.anno <- genes[, c(1:5)]
genes <- genes[, c(6:dim(genes)[2])]
head(genes)

library(edgeR)
y <- DGEList(counts=genes,group=c(1,1,1,2,2))
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

DEgenes <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) >1 ,]
DEintrons <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 1 ,] 
DEintrons_up <- DEintrons[DEintrons$logFC < 0,]
DEintrons_down <- DEintrons[DEintrons$logFC > 0,]

DEgenes2 <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes2 <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) > 2,]
DEintrons2 <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons2 <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 2 ,] 

dim(DEgenes)
dim(DEintrons)
dim(DEintrons_up)
dim(DEintrons_down)


write.table(DEintrons,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_HRSV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEintrons2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_HRSV_FDR0.5FC2.txt",sep="\t",quote=FALSE)
write.table(DEgenes,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_HRSV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEgenes2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_HRSV_FDR0.5FC2.txt",sep="\t",quote=FALSE)

gcs <- subset(genes,!grepl("_intron$",rownames(genes)))
gcs <- colSums(gcs)/colSums(genes)
ics <- subset(genes,grepl("_intron$",rownames(genes)))
ics <- colSums(ics)/colSums(genes)
data.frame("HRSV",gcs,ics)



##########################################################################################
##########################################################################################
##
## Intron retention -- FeatureCounts - HSV-1 with introns - HUMAN
##
##########################################################################################
##########################################################################################

genes <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/HSV-1_counts_introns.mx",sep="\t", skip=1, header = T, row.names=1) 
colnames(genes) <- sub("X.home.mmacchie.human.encode.infection_TEtrans.HSV.1.bams_no_virus_UCSC.", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.sorted_DS.bam", "", colnames(genes))
colnames(genes) <- sub("...bams_no_virus_UCSC.", "", colnames(genes))

hg38.viruses <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/virus.genes.hg38.txt",sep="\t",header=FALSE)

# separate counts columns from annotation columns
genes <- genes[!rownames(genes) %in% hg38.viruses$V1,]
gene.anno <- genes[, c(1:5)]
genes <- genes[, c(6:dim(genes)[2])]
head(genes)

library(edgeR)
y <- DGEList(counts=genes,group=c(1,1,2,2))
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

DEgenes <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) >1 ,]
DEintrons <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 1 ,] 
DEintrons_up <- DEintrons[DEintrons$logFC < 0,]
DEintrons_down <- DEintrons[DEintrons$logFC > 0,]

DEgenes2 <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes2 <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) > 2,]
DEintrons2 <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons2 <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 2 ,] 

dim(DEgenes)
dim(DEintrons)
dim(DEintrons_up)
dim(DEintrons_down)


write.table(DEintrons,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_HSV-1_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEintrons2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_HSV-1_FDR0.5FC2.txt",sep="\t",quote=FALSE)
write.table(DEgenes,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_HSV-1_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEgenes2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_HSV-1_FDR0.5FC2.txt",sep="\t",quote=FALSE)

gcs <- subset(genes,!grepl("_intron$",rownames(genes)))
gcs <- colSums(gcs)/colSums(genes)
ics <- subset(genes,grepl("_intron$",rownames(genes)))
ics <- colSums(ics)/colSums(genes)
data.frame("HSV-1",gcs,ics)

##########################################################################################
##########################################################################################
##
## Intron retention -- FeatureCounts - ORFV with introns - HUMAN
##
##########################################################################################
##########################################################################################

genes <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/ORFV_counts_introns.mx",sep="\t", skip=1, header = T, row.names=1) 
colnames(genes) <- sub("X.home.mmacchie.human.encode.infection_TEtrans.ORFV.bams_no_virus_UCSC.", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.sorted_DS.bam", "", colnames(genes))
colnames(genes) <- sub("...bams_no_virus_UCSC.", "", colnames(genes))

hg38.viruses <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/virus.genes.hg38.txt",sep="\t",header=FALSE)

# separate counts columns from annotation columns
genes <- genes[!rownames(genes) %in% hg38.viruses$V1,]
gene.anno <- genes[, c(1:5)]
genes <- genes[, c(6:dim(genes)[2])]
head(genes)


library(edgeR)
y <- DGEList(counts=genes,group=c(1,1,1,2,2,2))
keep <- rowSums(cpm(y)>1) >= 3
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

DEgenes <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) >1 ,]
DEintrons <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 1 ,] 
DEintrons_up <- DEintrons[DEintrons$logFC < 0,]
DEintrons_down <- DEintrons[DEintrons$logFC > 0,]

DEgenes2 <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes2 <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) > 2,]
DEintrons2 <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons2 <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 2 ,] 

dim(DEgenes)
dim(DEintrons)
dim(DEintrons_up)
dim(DEintrons_down)


write.table(DEintrons,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_ORFV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEintrons2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_ORFV_FDR0.5FC2.txt",sep="\t",quote=FALSE)
write.table(DEgenes,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_ORFV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEgenes2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_ORFV_FDR0.5FC2.txt",sep="\t",quote=FALSE)

gcs <- subset(genes,!grepl("_intron$",rownames(genes)))
gcs <- colSums(gcs)/colSums(genes)
ics <- subset(genes,grepl("_intron$",rownames(genes)))
ics <- colSums(ics)/colSums(genes)
data.frame("ORFV",gcs,ics)

##########################################################################################
##########################################################################################
##
## Intron retention -- FeatureCounts - VZV with introns - HUMAN
##
##########################################################################################
##########################################################################################

genes <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/VZV_counts_introns.mx",sep="\t", skip=1, header = T, row.names=1) 
colnames(genes) <- sub("X.home.mmacchie.human.encode.infection_TEtrans.VZV.bams_no_virus_UCSC.", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.sorted_DS.bam", "", colnames(genes))
colnames(genes) <- sub("...bams_no_virus_UCSC.", "", colnames(genes))

hg38.viruses <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/virus.genes.hg38.txt",sep="\t",header=FALSE)

# separate counts columns from annotation columns
genes <- genes[!rownames(genes) %in% hg38.viruses$V1,]
gene.anno <- genes[, c(1:5)]
genes <- genes[, c(6:dim(genes)[2])]
head(genes)


library(edgeR)
y <- DGEList(counts=genes,group=c(1,1,2,2))
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

DEgenes <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) >1 ,]
DEintrons <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 1 ,] 
DEintrons_up <- DEintrons[DEintrons$logFC < 0,]
DEintrons_down <- DEintrons[DEintrons$logFC > 0,]

DEgenes2 <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes2 <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) > 2,]
DEintrons2 <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons2 <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 2 ,] 

dim(DEgenes)
dim(DEintrons)
dim(DEintrons_up)
dim(DEintrons_down)


write.table(DEintrons,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_VZV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEintrons2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_VZV_FDR0.5FC2.txt",sep="\t",quote=FALSE)
write.table(DEgenes,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_VZV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEgenes2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_VZV_FDR0.5FC2.txt",sep="\t",quote=FALSE)

gcs <- subset(genes,!grepl("_intron$",rownames(genes)))
gcs <- colSums(gcs)/colSums(genes)
ics <- subset(genes,grepl("_intron$",rownames(genes)))
ics <- colSums(ics)/colSums(genes)
data.frame("VZV",gcs,ics)





##########################################################################################
##########################################################################################
##
## Intron retention -- FeatureCounts - DFV with introns - HUMAN
##
##########################################################################################
##########################################################################################

genes <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/DFV_counts_introns.mx",sep="\t", skip=1, header = T, row.names=1) 
colnames(genes) <- sub("X.home.mmacchie.human.encode.infection_TEtrans.DFV.bams_no_virus_UCSC.", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.sorted_DS.bam", "", colnames(genes))
colnames(genes) <- sub("...bams_no_virus_UCSC.", "", colnames(genes))

hg38.viruses <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/virus.genes.hg38.txt",sep="\t",header=FALSE)

# separate counts columns from annotation columns
genes <- genes[!rownames(genes) %in% hg38.viruses$V1,]
gene.anno <- genes[, c(1:5)]
genes <- genes[, c(6:dim(genes)[2])]
head(genes)


library(edgeR)
y <- DGEList(counts=genes,group=c(1,1,1,2,2,2))
keep <- rowSums(cpm(y)>1) >= 3
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

DEgenes <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) >1 ,]
DEintrons <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 1 ,] 
DEintrons_up <- DEintrons[DEintrons$logFC < 0,]
DEintrons_down <- DEintrons[DEintrons$logFC > 0,]

DEgenes2 <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes2 <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) > 2,]
DEintrons2 <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons2 <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 2 ,] 

dim(DEgenes)
dim(DEintrons)
dim(DEintrons_up)
dim(DEintrons_down)


write.table(DEintrons,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_DFV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEintrons2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_DFV_FDR0.5FC2.txt",sep="\t",quote=FALSE)
write.table(DEgenes,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_DFV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEgenes2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_DFV_FDR0.5FC2.txt",sep="\t",quote=FALSE)


gcs <- subset(genes,!grepl("_intron$",rownames(genes)))
gcs <- colSums(gcs)/colSums(genes)
ics <- subset(genes,grepl("_intron$",rownames(genes)))
ics <- colSums(ics)/colSums(genes)
data.frame("DFV",gcs,ics)



##########################################################################################
##########################################################################################
##
## Intron retention -- FeatureCounts - HCV with introns - HUMAN
##
##########################################################################################
##########################################################################################

genes <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/HCV_counts_introns.mx",sep="\t", skip=1, header = T, row.names=1) 
colnames(genes) <- sub("X.home.mmacchie.human.encode.infection_TEtrans.HCV.bams_no_virus_UCSC.", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.sorted_DS.bam", "", colnames(genes))
colnames(genes) <- sub("...bams_no_virus_UCSC.", "", colnames(genes))

hg38.viruses <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/virus.genes.hg38.txt",sep="\t",header=FALSE)

# separate counts columns from annotation columns
genes <- genes[!rownames(genes) %in% hg38.viruses$V1,]
gene.anno <- genes[, c(1:5)]
genes <- genes[, c(6:dim(genes)[2])]
head(genes)


library(edgeR)
y <- DGEList(counts=genes,group=c(1,1,1,1,1,2,2,2,2,2))
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

DEgenes <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) >1 ,]
DEintrons <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 1 ,] 
DEintrons_up <- DEintrons[DEintrons$logFC < 0,]
DEintrons_down <- DEintrons[DEintrons$logFC > 0,]

DEgenes2 <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes2 <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) > 2,]
DEintrons2 <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons2 <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 2 ,] 

dim(DEgenes)
dim(DEintrons)
dim(DEintrons_up)
dim(DEintrons_down)


write.table(DEintrons,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_HCV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEintrons2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_HCV_FDR0.5FC2.txt",sep="\t",quote=FALSE)
write.table(DEgenes,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_HCV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEgenes2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_HCV_FDR0.5FC2.txt",sep="\t",quote=FALSE)


gcs <- subset(genes,!grepl("_intron$",rownames(genes)))
gcs <- colSums(gcs)/colSums(genes)
ics <- subset(genes,grepl("_intron$",rownames(genes)))
ics <- colSums(ics)/colSums(genes)
data.frame("HCV",gcs,ics)



##########################################################################################
##########################################################################################
##
## Intron retention -- FeatureCounts - IAV with introns - HUMAN
##
##########################################################################################
##########################################################################################

genes <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/IAV_counts_introns.mx",sep="\t", skip=1, header = T, row.names=1) 
colnames(genes) <- sub("X.home.mmacchie.human.encode.infection_TEtrans.IAV.bams_no_virus_UCSC.", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.sorted_DS.bam", "", colnames(genes))
colnames(genes) <- sub("...bams_no_virus_UCSC.", "", colnames(genes))

hg38.viruses <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/virus.genes.hg38.txt",sep="\t",header=FALSE)

# separate counts columns from annotation columns
genes <- genes[!rownames(genes) %in% hg38.viruses$V1,]
gene.anno <- genes[, c(1:5)]
genes <- genes[, c(6:dim(genes)[2])]
head(genes)


library(edgeR)
y <- DGEList(counts=genes,group=c(1,1,2,2))
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

DEgenes <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) >1 ,]
DEintrons <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 1 ,] 
DEintrons_up <- DEintrons[DEintrons$logFC < 0,]
DEintrons_down <- DEintrons[DEintrons$logFC > 0,]

DEgenes2 <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes2 <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) > 2,]
DEintrons2 <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons2 <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 2 ,] 

dim(DEgenes)
dim(DEintrons)
dim(DEintrons_up)
dim(DEintrons_down)


write.table(DEintrons,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_IAV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEintrons2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_IAV_FDR0.5FC2.txt",sep="\t",quote=FALSE)
write.table(DEgenes,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_IAV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEgenes2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_IAV_FDR0.5FC2.txt",sep="\t",quote=FALSE)

gcs <- subset(genes,!grepl("_intron$",rownames(genes)))
gcs <- colSums(gcs)/colSums(genes)
ics <- subset(genes,grepl("_intron$",rownames(genes)))
ics <- colSums(ics)/colSums(genes)
data.frame("IAV",gcs,ics)

##########################################################################################
##########################################################################################
##
## Intron retention -- FeatureCounts - EBV with introns - HUMAN
##
##########################################################################################
##########################################################################################

genes <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/EBV_counts_introns.mx",sep="\t", skip=1, header = T, row.names=1) 
colnames(genes) <- sub("X.home.mmacchie.human.encode.infection_TEtrans.EBV.bams_no_virus_UCSC.", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.sorted_DS.bam", "", colnames(genes))
colnames(genes) <- sub("...bams_no_virus_UCSC.", "", colnames(genes))

hg38.viruses <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/virus.genes.hg38.txt",sep="\t",header=FALSE)

# separate counts columns from annotation columns
genes <- genes[!rownames(genes) %in% hg38.viruses$V1,]
gene.anno <- genes[, c(1:5)]
genes <- genes[, c(6:dim(genes)[2])]

library(edgeR)
y <- DGEList(counts=genes,group=c(1,1,1,1,1,2,2,2,2,2))
keep <- rowSums(cpm(y)>1) >= 5
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

DEgenes <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) >1 ,]
DEintrons <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 1 ,]
DEintrons_up <- DEintrons[DEintrons$logFC < 0,]
DEintrons_down <- DEintrons[DEintrons$logFC > 0,]

DEgenes2 <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes2 <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) > 2,]
DEintrons2 <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons2 <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 2 ,] 

dim(DEgenes)
dim(DEintrons)
dim(DEintrons_up)
dim(DEintrons_down)

write.table(DEintrons,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_EBV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEintrons2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_EBV_FDR0.5FC2.txt",sep="\t",quote=FALSE)
write.table(DEgenes,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_EBV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEgenes2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_EBV_FDR0.5FC2.txt",sep="\t",quote=FALSE)


gcs <- subset(genes,!grepl("_intron$",rownames(genes)))
gcs <- colSums(gcs)/colSums(genes)
ics <- subset(genes,grepl("_intron$",rownames(genes)))
ics <- colSums(ics)/colSums(genes)
data.frame("EBV",gcs,ics)

##########################################################################################
##########################################################################################
##
## Intron retention -- FeatureCounts - MeV with introns - HUMAN
##
##########################################################################################
##########################################################################################

genes <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/MeV_counts_introns.mx",sep="\t", skip=1, header = T, row.names=1) 
colnames(genes) <- sub("X.home.mmacchie.human.encode.infection_TEtrans.MeV.bams_no_virus_UCSC.", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.sorted_DS.bam", "", colnames(genes))
colnames(genes) <- sub("...bams_no_virus_UCSC.", "", colnames(genes))

hg38.viruses <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/virus.genes.hg38.txt",sep="\t",header=FALSE)

# separate counts columns from annotation columns
genes <- genes[!rownames(genes) %in% hg38.viruses$V1,]
gene.anno <- genes[, c(1:5)]
genes <- genes[, c(6:dim(genes)[2])]

library(edgeR)
y <- DGEList(counts=genes,group=c(1,1,1,2,2,2))
keep <- rowSums(cpm(y)>1) >= 3
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

DEgenes <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) >1 ,]
DEintrons <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 1 ,]
DEintrons_up <- DEintrons[DEintrons$logFC < 0,]
DEintrons_down <- DEintrons[DEintrons$logFC > 0,]

DEgenes2 <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes2 <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) > 2,]
DEintrons2 <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons2 <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 2 ,] 

dim(DEgenes)
dim(DEintrons)
dim(DEintrons_up)
dim(DEintrons_down)

write.table(DEintrons,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_MeV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEintrons2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_MeV_FDR0.5FC2.txt",sep="\t",quote=FALSE)
write.table(DEgenes,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_MeV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEgenes2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_MeV_FDR0.5FC2.txt",sep="\t",quote=FALSE)

gcs <- subset(genes,!grepl("_intron$",rownames(genes)))
gcs <- colSums(gcs)/colSums(genes)
ics <- subset(genes,grepl("_intron$",rownames(genes)))
ics <- colSums(ics)/colSums(genes)
data.frame("MeV",gcs,ics)


##########################################################################################
##########################################################################################
##
## Intron retention -- FeatureCounts - SINV with introns - HUMAN
##
##########################################################################################
##########################################################################################

genes <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/SINV_counts_introns.mx",sep="\t", skip=1, header = T, row.names=1) 
colnames(genes) <- sub("X.home.mmacchie.human.encode.infection_TEtrans.SINV.bams_no_virus_UCSC.", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.sorted_DS.bam", "", colnames(genes))
colnames(genes) <- sub("...bams_no_virus_UCSC.", "", colnames(genes))

hg38.viruses <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/virus.genes.hg38.txt",sep="\t",header=FALSE)

# separate counts columns from annotation columns
genes <- genes[!rownames(genes) %in% hg38.viruses$V1,]
gene.anno <- genes[, c(1:5)]
genes <- genes[, c(6:dim(genes)[2])]

library(edgeR)
y <- DGEList(counts=genes,group=c(1,1,1,2,2,2))
keep <- rowSums(cpm(y)>1) >= 3
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

DEgenes <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) >1 ,]
DEintrons <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 1 ,]
DEintrons_up <- DEintrons[DEintrons$logFC < 0,]
DEintrons_down <- DEintrons[DEintrons$logFC > 0,]

DEgenes2 <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes2 <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) > 2,]
DEintrons2 <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons2 <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 2 ,] 

dim(DEgenes)
dim(DEintrons)
dim(DEintrons_up)
dim(DEintrons_down)

write.table(DEintrons,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_SINV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEintrons2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_SINV_FDR0.5FC2.txt",sep="\t",quote=FALSE)
write.table(DEgenes,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_SINV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEgenes2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_SINV_FDR0.5FC2.txt",sep="\t",quote=FALSE)



gcs <- subset(genes,!grepl("_intron$",rownames(genes)))
gcs <- colSums(gcs)/colSums(genes)
ics <- subset(genes,grepl("_intron$",rownames(genes)))
ics <- colSums(ics)/colSums(genes)
data.frame("SINV",gcs,ics)


##########################################################################################
##########################################################################################
##
## Intron retention -- FeatureCounts - HCMV with introns - HUMAN
##
##########################################################################################
##########################################################################################

genes <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/HCMV_counts_introns.mx",sep="\t", skip=1, header = T, row.names=1) 
colnames(genes) <- sub("X.home.mmacchie.human.encode.infection_TEtrans.HCMV.bams_no_virus_UCSC.", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.sorted_DS.bam", "", colnames(genes))
colnames(genes) <- sub("...bams_no_virus_UCSC.", "", colnames(genes))

hg38.viruses <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/virus.genes.hg38.txt",sep="\t",header=FALSE)

# separate counts columns from annotation columns
genes <- genes[!rownames(genes) %in% hg38.viruses$V1,]
gene.anno <- genes[, c(1:5)]
genes <- genes[, c(6:dim(genes)[2])]

library(edgeR)
y <- DGEList(counts=genes,group=c(1,1,1,2,2,2))
keep <- rowSums(cpm(y)>1) >= 3
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

DEgenes <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) >1 ,]
DEintrons <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 1 ,]
DEintrons_up <- DEintrons[DEintrons$logFC < 0,]
DEintrons_down <- DEintrons[DEintrons$logFC > 0,]

DEgenes2 <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes2 <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) > 2,]
DEintrons2 <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons2 <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 2 ,] 

dim(DEgenes)
dim(DEintrons)
dim(DEintrons_up)
dim(DEintrons_down)

write.table(DEintrons,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_HCMV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEintrons2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_HCMV_FDR0.5FC2.txt",sep="\t",quote=FALSE)
write.table(DEgenes,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_HCMV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEgenes2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_HCMV_FDR0.5FC2.txt",sep="\t",quote=FALSE)


gcs <- subset(genes,!grepl("_intron$",rownames(genes)))
gcs <- colSums(gcs)/colSums(genes)
ics <- subset(genes,grepl("_intron$",rownames(genes)))
ics <- colSums(ics)/colSums(genes)
data.frame("HCMV",gcs,ics)






##########################################################################################
##########################################################################################
##
## Intron retention -- FeatureCounts - RESTV with introns - HUMAN
##
##########################################################################################
##########################################################################################

genes <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/RESTV_counts_introns.mx",sep="\t", skip=1, header = T, row.names=1) 
colnames(genes) <- sub("X.home.mmacchie.human.encode.infection_TEtrans.RESTV.bams_no_virus_UCSC.", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.sorted_DS.bam", "", colnames(genes))
colnames(genes) <- sub("...bams_no_virus_UCSC.", "", colnames(genes))

hg38.viruses <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/virus.genes.hg38.txt",sep="\t",header=FALSE)

# separate counts columns from annotation columns
genes <- genes[!rownames(genes) %in% hg38.viruses$V1,]
gene.anno <- genes[, c(1:5)]
genes <- genes[, c(6:dim(genes)[2])]

library(edgeR)
y <- DGEList(counts=genes,group=c(1,1,1,2,2,2))
keep <- rowSums(cpm(y)>1) >= 3
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

DEgenes <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) >1 ,]
DEintrons <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 1 ,]
DEintrons_up <- DEintrons[DEintrons$logFC < 0,]
DEintrons_down <- DEintrons[DEintrons$logFC > 0,]

DEgenes2 <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes2 <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) > 2,]
DEintrons2 <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons2 <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 2 ,] 

dim(DEgenes)
dim(DEintrons)
dim(DEintrons_up)
dim(DEintrons_down)

write.table(DEintrons,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_RESTV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEintrons2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_RESTV_FDR0.5FC2.txt",sep="\t",quote=FALSE)
write.table(DEgenes,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_RESTV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEgenes2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_RESTV_FDR0.5FC2.txt",sep="\t",quote=FALSE)

gcs <- subset(genes,!grepl("_intron$",rownames(genes)))
gcs <- colSums(gcs)/colSums(genes)
ics <- subset(genes,grepl("_intron$",rownames(genes)))
ics <- colSums(ics)/colSums(genes)
data.frame("RESTV",gcs,ics)



##########################################################################################
##########################################################################################
##
## Intron retention -- FeatureCounts - RVFV with introns - HUMAN
##
##########################################################################################
##########################################################################################

genes <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/RVFV_counts_introns.mx",sep="\t", skip=1, header = T, row.names=1) 
colnames(genes) <- sub("X.home.mmacchie.human.encode.infection_TEtrans.RVFV.bams_no_virus_UCSC.", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.sorted_DS.bam", "", colnames(genes))
colnames(genes) <- sub("...bams_no_virus_UCSC.", "", colnames(genes))

hg38.viruses <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/virus.genes.hg38.txt",sep="\t",header=FALSE)

# separate counts columns from annotation columns
genes <- genes[!rownames(genes) %in% hg38.viruses$V1,]
gene.anno <- genes[, c(1:5)]
genes <- genes[, c(6:dim(genes)[2])]

genes <-genes[,-c(1:3)]
library(edgeR)
y <- DGEList(counts=genes,group=c(1,1,1,2,2,2))
keep <- rowSums(cpm(y)>1) >= 3
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

DEgenes <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) >1 ,]
DEintrons <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 1 ,]
DEintrons_up <- DEintrons[DEintrons$logFC < 0,]
DEintrons_down <- DEintrons[DEintrons$logFC > 0,]

DEgenes2 <- subset(tt,!grepl("_intron$",rownames(tt))) 
DEgenes2 <- DEgenes[DEgenes$FDR < 0.05 & abs(DEgenes$logFC) > 2,]
DEintrons2 <- subset(tt,grepl("_intron$",rownames(tt))) 
DEintrons2 <- DEintrons[DEintrons$FDR < 0.05 & abs(DEintrons$logFC) > 2 ,] 

dim(DEgenes)
dim(DEintrons)
dim(DEintrons_up)
dim(DEintrons_down)

write.table(DEintrons,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_RVFV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEintrons2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEintrons_RVFV_FDR0.5FC2.txt",sep="\t",quote=FALSE)
write.table(DEgenes,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_RVFV_FDR0.5FC1.txt",sep="\t",quote=FALSE)
write.table(DEgenes2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/DEgenes_RVFV_FDR0.5FC2.txt",sep="\t",quote=FALSE)

gcs <- subset(genes,!grepl("_intron$",rownames(genes)))
gcs <- colSums(gcs)/colSums(genes)
ics <- subset(genes,grepl("_intron$",rownames(genes)))
ics <- colSums(ics)/colSums(genes)
data.frame("RVFV",gcs,ics)





##########################################################################################
##########################################################################################
##
## Compare DE intron in human
##
##########################################################################################
##########################################################################################

setwd("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/intron_retention/results/")
fullmaster <- read.delim("/Users/mmacchie/Desktop/TE_project/human_encode/annot/master.hg38.gene.to.te.table.txt")
#d <- read.delim("hg38.intron.retention.heatmap.wsummary.txt")
d <- read.delim("hg38.intron.retention.heatmap.genecorrected.wsummary.txt")
rownames(d) <- d$X
d <- d[,-c(1,2)]
d <- d[,1:16]
library(pheatmap)
pheatmap(d)

# keep genes that have IR changes in at least 2 data sets
d <- read.delim("hg38.intron.retention.heatmap.genecorrected.wsummary.txt")
rownames(d) <- d$X
d <- d[,-c(1,2)]
d <- d[d$totDE >=2,]
d <- d[,1:16]
pheatmap(d)


#intupINF <- read.delim("hg38.genes.splicing.upINF.min3datasets.txt",header=FALSE)
#intdownINF <- read.delim("hg38.genes.splicing.downINF.min3datasets.txt",header=FALSE)

intupINF <- read.delim("hg38.genes.splicing.upINF.min3datasets.genecorrected.txt",header=FALSE)
intdownINF <- read.delim("hg38.genes.splicing.downINF.min3datasets.genecorrected.txt",header=FALSE)

library(biomaRt)
listMarts()
ensembl = useMart("ensembl")
ensembl = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
listFilters(ensembl)
listAttributes(ensembl)

master_gene_hg38 <- read.delim("/Users/mmacchie/Desktop/TE_project/human_encode/annot/master.hg38.gene.to.te.table.txt",sep="\t",header=F)
colnames(master_gene_hg38) <- c("gene","gene_length","eff_gene_length","num_introns","cum_intron_length","num_exons","cum_exon_length","gene_strand","gene_biotype","TE_IDs","TE-gene-relations","TE-directions","num_TEs","num_TEs_exons","num_TEs_introns")

library(org.Hs.eg.db)
library(goseq)
library(GO.db)

#test <- intdownINF$V1
test <- intupINF$V1
genes2 <- unique(master_gene_hg38$gene) #58243
genes <- genes2
genes <- as.integer(genes %in% as.character(as.factor(test)))
names(genes) <- genes2 
pwf=nullp(genes,"hg38","ensGene")
GO.wall=goseq(pwf,"hg38","ensGene")
GO.enrich = GO.wall$category[p.adjust(GO.wall$over_represented_pvalue,method = "BH")<0.05]
GO.terms <- as.data.frame(Term(GO.enrich))
GO.terms <- cbind(GO.terms, adj.pvalue=p.adjust(GO.wall$over_represented_pvalue, method="BH")[match(rownames(GO.terms), GO.wall$category)])
GO.terms <- GO.terms[GO.terms$adj.pvalue < 0.05,]
fdrs <- -log(GO.terms$adj.pvalue)
names(fdrs) <- GO.terms$`Term(GO.enrich)`
#pdf(file="GOterm_intron_downinfection_1543genes_3mindata.genecorrected.pdf",width=10,height=12)
pdf(file="GOterm_intron_upinfection_1543genes_3mindata.genecorrected.pdf",width=10,height=12)
par(mar=c(5,30,5,5)) 
barplot(sort(fdrs,decreasing=FALSE),horiz=T, las=1,col="grey", cex.names = 0.7,xlab = "-log(adjusted p-value)",xlim = c(0,30))
dev.off()



##########################################################################################
##########################################################################################
##
## Compare DE gene in human
##
##########################################################################################
##########################################################################################

fullmaster <- read.delim("/Users/mmacchie/Desktop/TE_project/human_encode/annot/master.te.to.gene.table.hg38.exonannot.txt",header=FALSE)

d <- read.delim("shared.DEgenes.FDR0.05.FC1.humanviruses.wsummary.txt")
#test<- d[,1:16]
test <- d[d$totDE >= 5,1:16] # in at least 5 virus data sets (628 genes)
library(pheatmap)
pheatmap(test)

#d <- read.delim("shared.DEgenes.FDR0.05.FC2.humanviruses.wsummary.txt")
#test <- d[d$totDE >= 3,1:16] # in at least 5 virus data sets (628 genes)

back <- unique(fullmaster[!fullmaster$V15=="autonomous",c(10,12,14)])
sum(table(back$V12))
back2 <- table(back$V12)/sum(table(back$V12))

DE_hostgenes <- rownames(test)
fsub <- fullmaster[fullmaster$V14 %in% rownames(test),c(10,12,14)]
fsub2 <- table(fsub$V12)/sum(table(fsub$V12))

x <- table(fsub$V12)[36:45]#/sum(table(fsub$V12))

##-------- bootstrapping --- test-----#
allgenes <- unique(fullmaster$V14)
sample_gene <- sample(1:length(allgenes),628,replace=FALSE)

bootstrap <- data.frame(LTR=integer(),LTRq=integer(),LTR_ERV1=integer(),LTR_ERV1q=integer(),LTR_ERVK=integer(),LTR_ERVL=integer(),LTR_ERVL_MaLR=integer(),LTR_ERVLq=integer(),LTR_Gypsy=integer(),LTR_Gypsyq=integer())
for (i in 1:1000){
	set.seed(i)
	sample_gene <- sample(1:length(allgenes),628,replace=FALSE)
	sample_gene <- allgenes[sample_gene]
	sample_rep <- unique(fullmaster[fullmaster$V14 %in% sample_gene,c(10,12,14)])
	x <- table(sample_rep$V12)[36:45]
	bootstrap <- rbind(bootstrap,x)
}

colnames(bootstrap) <- c("LTR","LTRq","LTR_ERV1","LTR_ERV1q","LTR_ERVK","LTR_ERVL","LTR_ERVL_MaLR","LTR_ERVLq","LTR_Gypsy","LTR_Gypsyq")

##-------end of bootstraping test--------#

# plot the levels of repeats around these genes
barplot(table(fsub$V12),las=2)

table(fsub$V12)[36:45][3]
hist(bootstrap[,3],breaks=2000,col="orange",xlim=c(500,2000))


##-----Fisher's exact enrichment test - ERVs only-----#

Esubs = names(table(fsub$V12)[36:45])
count = 0
for (Esub in table(fsub$V12)[36:45]){
	count = count + 1
	b1 <- Esub # repeat is ERV subtype X, repeat is not autonomous, repeat is around DE gene
	b2 <- sum(table(fsub$V12)) - b1 # repeat is not ERV subtype X, repeat is not autonomous, repeat is around DE gene
	pre_b3 <- dim(fullmaster[fullmaster$V12 == Esubs[count] & !fullmaster$V15 == "autonomous" & !fullmaster$V14 %in% DE_hostgenes,])[1] # repeat is ERV subtype X, repeat is not autonomous, and repeat is not in DE gene
	b3 <- pre_b3 - b1
	b4 <- dim(fullmaster[!fullmaster$V12 == Esubs[count] & !fullmaster$V15 == "autonomous" & !fullmaster$V14 %in% DE_hostgenes,])[1] # repeat is not ERV subtype X, repeat is not autonomous, and repeat is not in DE gene
	cmat <- matrix(c(b1,b3,b2,b4),ncol=2)
	p = fisher.test(cmat)
	#print(c(Esubs[count],p))
	print(c(Esubs[count],p$p.value))
}

##-----end of Fisher's exact enrichment test - ERVs only-----#

##-----Fisher's exact enrichment test - all repeats-----#

Esubs = names(table(fsub$V12)[1:66])
rep_table <- data.frame(matrix(ncol = 6, nrow = 1))
colnames(rep_table) <- c("Repeat_family","p-value","DE_rep","DE_otherrep","nonDE_rep","nonDE_otherrep")

count = 0
for (Esub in table(fsub$V12)[1:66]){
	count = count + 1
	b1 <- Esub # repeat is ERV subtype X, repeat is not autonomous, repeat is around DE gene
	b2 <- sum(table(fsub$V12)) - b1 # repeat is not ERV subtype X, repeat is not autonomous, repeat is around DE gene
	pre_b3 <- dim(fullmaster[fullmaster$V12 == Esubs[count] & !fullmaster$V15 == "autonomous" & !fullmaster$V14 %in% DE_hostgenes,])[1] # repeat is ERV subtype X, repeat is not autonomous, and repeat is not in DE gene
	b3 <- pre_b3 - b1
	b4 <- dim(fullmaster[!fullmaster$V12 == Esubs[count] & !fullmaster$V15 == "autonomous" & !fullmaster$V14 %in% DE_hostgenes,])[1] # repeat is not ERV subtype X, repeat is not autonomous, and repeat is not in DE gene
	cmat <- matrix(c(b1,b3,b2,b4),ncol=2)
	p = fisher.test(cmat,alternative="greater")
	#print(c(Esubs[count],p))
	print(c(Esubs[count],p$p.value,b1,b2,b3,b4))
	rep_table <- rbind(rep_table,c(Esubs[count],p$p.value,b1,b2,b3,b4))
}
print(rep_table)

##-----Fisher's exact enrichment test - all repeats-----#


sessionInfo()
