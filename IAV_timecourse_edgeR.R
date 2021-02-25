##########################################################################################
##########################################################################################
##
## IAV 7-day timecourse - Differential expression analysis
##
##########################################################################################
##########################################################################################

# This code identifies differentially expressed genes, repeats/TEs by comparing 0h (mock) to each infection time point
# Time points that are within 1-2 hours apart are grouped as biological replicates

# load repeat-gene relationship table
fullmaster <- read.delim("~/Desktop/TE_project/mouse_encode/counts/ENCODE-processed-totalRNA/MM_EP_STAR_FC_DFAM_finalresults/master.te.to.gene.table.FULL.exonannot2.txt",skip=1,header=FALSE)

# load gene and repeat counts matrices
genes <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/featurecounts/IAV_timecourse_gene_counts3.mx",sep="\t", skip=1, header = T, row.names=1)
reps <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/featurecounts/IAV_timecourse_repeat_counts3.mx",sep="\t", skip=1, header = T, row.names=1)

# update the sample IDs on the gene counts matrix
colnames(genes) <- sub("X.home.mmacchie.mouse.encode.virus.IAV_timecourse.bams.", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.sorted_DS.bam", "", colnames(genes))
colnames(genes) <- sub("...bams_no_virus_UCSC.", "", colnames(genes))
colnames(genes) <- sub("bams_no_virus_UCSC.","",colnames(genes))

# update the sample IDs on the repeat counts matrix
colnames(reps) <- sub("X.home.mmacchie.mouse.encode.virus.IAV_timecourse.bams.", "", colnames(reps))
colnames(reps) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(reps))
colnames(reps) <- sub("Aligned.sortedByCoord.out.sorted_DS.bam", "", colnames(reps))
colnames(reps) <- sub("...bams_no_virus_UCSC.", "", colnames(reps))
colnames(reps) <- sub("bams_no_virus_UCSC.","",colnames(reps))

# Remove simple repeat elements from the repeat counts matrix
SRs <- fullmaster[fullmaster$V12 =="Simple_repeat",]$V10
reps <- reps[!rownames(reps) %in% SRs,]
reps <- na.omit(reps)
gene.anno <- genes[, c(1:5)]
genes <- genes[, c(6:dim(genes)[2])]
rep.anno <- reps[, c(1:5)]
reps <- reps[, c(6:dim(reps)[2])]

# Combine the repeat and gene counts matrices
#anno <- rbind(gene.anno,rep.anno)
genes <- rbind(genes,reps)
library(edgeR)
rpkms <- rpkm(genes,anno$Length)
cpms <- cpm(genes)
# Reorder the samples in the combined counts matrix
neworder = c("mock","mock_2","inf_3h","inf_4h","inf_11h","inf_12h","inf_26h","inf_28h","inf_49h","inf_50h","inf_74h","inf_75h","inf_98h","inf_99h","inf_122h","inf_123h","inf_148h","inf_149h")
cpms <- cpms[,match(neworder,colnames(cpms))]

# Subset the combined counts matrix for all the downstream infection time point comparisons to 0h
e1 <- genes[,c(9,11,19,20)] #3h
e2 <- genes[,c(1,4,19,20)] #11h
e3 <- genes[,c(7,8,19,20)] #26h
e4 <- genes[,c(10,12,19,20)] #49h
e5 <- genes[,c(13,14,19,20)] #74h
e6 <- genes[,c(17,18,19,20)] #98h
e7 <- genes[,c(2,3,19,20)] #122h
e8 <- genes[,c(5,6,19,20)] #149h

# Comparison of 3-4h to 0h
library(edgeR)
y <- DGEList(counts=e1,group=c(1,1,2,2))
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

tt <- tt[tt$FDR < 0.25 ,]
dgs1 <- subset(tt,grepl('^EN',rownames(tt)))
dreps1 <- subset(tt,grepl('^DF',rownames(tt)))
dim(dgs1)
dim(dreps1)

dgs1_up <- dgs1[dgs1$logFC < 0,]
dgs1_down <- dgs1[dgs1$logFC > 0,]
dreps1_up <- dreps1[dreps1$logFC < 0,]
dreps1_down <- dreps1[dreps1$logFC > 0,]

# Comparison of 11-12h to 0h
y <- DGEList(counts=e2,group=c(1,1,2,2))
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

tt <- tt[tt$FDR < 0.25 ,]
dgs2 <- subset(tt,grepl('^EN',rownames(tt)))
dreps2 <- subset(tt,grepl('^DF',rownames(tt)))
dim(dgs2)
dim(dreps2)

dgs2_up <- dgs2[dgs2$logFC < 0,]
dgs2_down <- dgs2[dgs2$logFC > 0,]
dreps2_up <- dreps2[dreps2$logFC < 0,]
dreps2_down <- dreps2[dreps2$logFC > 0,]

# Comparison of 26-28h to 0h
y <- DGEList(counts=e3,group=c(1,1,2,2))
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

tt <- tt[tt$FDR < 0.25 ,]
dgs3 <- subset(tt,grepl('^EN',rownames(tt)))
dreps3 <- subset(tt,grepl('^DF',rownames(tt)))
dim(dgs3)
dim(dreps3)

dgs3_up <- dgs3[dgs3$logFC < 0,]
dgs3_down <- dgs3[dgs3$logFC > 0,]
dreps3_up <- dreps3[dreps3$logFC < 0,]
dreps3_down <- dreps3[dreps3$logFC > 0,]

# Comparison of 49-50h to 0h
y <- DGEList(counts=e4,group=c(1,1,2,2))
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

tt <- tt[tt$FDR < 0.25 ,]
dgs4 <- subset(tt,grepl('^EN',rownames(tt)))
dreps4 <- subset(tt,grepl('^DF',rownames(tt)))
dim(dgs4)
dim(dreps4)

dgs4_up <- dgs4[dgs4$logFC < 0,]
dgs4_down <- dgs4[dgs4$logFC > 0,]
dreps4_up <- dreps4[dreps4$logFC < 0,]
dreps4_down <- dreps4[dreps4$logFC > 0,]

# Comparison of 74-75h to 0h
y <- DGEList(counts=e5,group=c(1,1,2,2))
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

tt <- tt[tt$FDR < 0.25 ,]
dgs5 <- subset(tt,grepl('^EN',rownames(tt)))
dreps5 <- subset(tt,grepl('^DF',rownames(tt)))
dim(dgs5)
dim(dreps5)

dgs5_up <- dgs5[dgs5$logFC < 0,]
dgs5_down <- dgs5[dgs5$logFC > 0,]
dreps5_up <- dreps5[dreps5$logFC < 0,]
dreps5_down <- dreps5[dreps5$logFC > 0,]

# Comparison of 98-99h to 0h
y <- DGEList(counts=e6,group=c(1,1,2,2))
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

tt <- tt[tt$FDR < 0.25 ,]
dgs6 <- subset(tt,grepl('^EN',rownames(tt)))
dreps6 <- subset(tt,grepl('^DF',rownames(tt)))
dim(dgs6)
dim(dreps6)

dgs6_up <- dgs6[dgs6$logFC < 0,]
dgs6_down <- dgs6[dgs6$logFC > 0,]
dreps6_up <- dreps6[dreps6$logFC < 0,]
dreps6_down <- dreps6[dreps6$logFC > 0,]

# Comparison of 122-123h to 0h
y <- DGEList(counts=e7,group=c(1,1,2,2))
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

tt <- tt[tt$FDR < 0.25 ,]
dgs7 <- subset(tt,grepl('^EN',rownames(tt)))
dreps7 <- subset(tt,grepl('^DF',rownames(tt)))
dim(dgs7)
dim(dreps7)

dgs7_up <- dgs7[dgs7$logFC < 0,]
dgs7_down <- dgs7[dgs7$logFC > 0,]
dreps7_up <- dreps7[dreps7$logFC < 0,]
dreps7_down <- dreps7[dreps7$logFC > 0,]

# Comparison of 149-150h to 0h
y <- DGEList(counts=e8,group=c(1,1,2,2))
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table

tt <- tt[tt$FDR < 0.25 ,]
dgs8 <- subset(tt,grepl('^EN',rownames(tt)))
dreps8 <- subset(tt,grepl('^DF',rownames(tt)))
dim(dgs8)
dim(dreps8)

dgs8_up <- dgs8[dgs8$logFC < 0,]
dgs8_down <- dgs8[dgs8$logFC > 0,]
dreps8_up <- dreps8[dreps8$logFC < 0,]
dreps8_down <- dreps8[dreps8$logFC > 0,]


totdgs <- c(1,dim(dgs1)[1],dim(dgs2)[1],dim(dgs3)[1],dim(dgs4)[1],dim(dgs5)[1],dim(dgs6)[1],dim(dgs7)[1],dim(dgs8)[1])
totdreps <- c(1,dim(dreps1)[1],dim(dreps2)[1],dim(dreps3)[1],dim(dreps4)[1],dim(dreps5)[1],dim(dreps6)[1],dim(dreps7)[1],dim(dreps8)[1])

TEs <- subset(fullmaster,grepl("^L",fullmaster$V12))
TEs <- subset(TEs,!grepl("^Low",TEs$V12))
TEs <- unique(TEs$V10)
dreps1_TEsonly <- dreps1[rownames(dreps1) %in% TEs,]
dreps2_TEsonly <- dreps2[rownames(dreps2) %in% TEs,]
dreps3_TEsonly <- dreps3[rownames(dreps3) %in% TEs,]
dreps4_TEsonly <- dreps4[rownames(dreps4) %in% TEs,]
dreps5_TEsonly <- dreps5[rownames(dreps5) %in% TEs,]
dreps6_TEsonly <- dreps6[rownames(dreps6) %in% TEs,]
dreps7_TEsonly <- dreps7[rownames(dreps7) %in% TEs,]
dreps8_TEsonly <- dreps8[rownames(dreps8) %in% TEs,]
totdreps_TEsonly <- c(1,dim(dreps1_TEsonly)[1],dim(dreps2_TEsonly)[1],dim(dreps3_TEsonly)[1],dim(dreps4_TEsonly)[1],dim(dreps5_TEsonly)[1],dim(dreps6_TEsonly)[1],dim(dreps7_TEsonly)[1],dim(dreps8_TEsonly)[1])
test <- c(rownames(dreps1_TEsonly),rownames(dreps2_TEsonly),rownames(dreps3_TEsonly),rownames(dreps4_TEsonly),rownames(dreps5_TEsonly),rownames(dreps6_TEsonly),rownames(dreps7_TEsonly),rownames(dreps8_TEsonly))
test <- fullmaster[fullmaster$V10 %in% test,c(1,4,5,7,10)]
write.table(test,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/MHC_analysis/all.DE.TEs.IAVtimecourse.FDR0.25.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)


early_TEs <- rbind(dreps1_TEsonly,dreps2_TEsonly,dreps3_TEsonly)
mid_TEs <- rbind(dreps4_TEsonly,dreps5_TEsonly)
late_TEs <- rbind(dreps6_TEsonly,dreps7_TEsonly,dreps8_TEsonly)


early_TEs <- unique(rownames(early_TEs))
early_TEs<- unique(fullmaster[fullmaster$V10 %in% early_TEs,c(1,4,5,12)])
mid_TEs <- unique(rownames(mid_TEs))
mid_TEs<- unique(fullmaster[fullmaster$V10 %in% mid_TEs,c(1,4,5,12)])
late_TEs <- unique(rownames(late_TEs))
late_TEs<- unique(fullmaster[fullmaster$V10 %in% late_TEs,c(1,4,5,12)])

write.table(early_TEs,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/GREAT/TE.coordinates.IAV.timecourse.earlystage.for.GREAT.FDR0.25.UCSC.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(mid_TEs,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/GREAT/TE.coordinates.IAV.timecourse.midstage.for.GREAT.UCSC.FDR0.25.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(late_TEs,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/GREAT/TE.coordinates.IAV.timecourse.latestage.for.GREAT.UCSC.FDR0.25.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)



updgs <- c(1,dim(dgs1_up)[1],dim(dgs2_up)[1],dim(dgs3_up)[1],dim(dgs4_up)[1],dim(dgs5_up)[1],dim(dgs6_up)[1],dim(dgs7_up)[1],dim(dgs8_up)[1])
upreps <- c(1,dim(dreps1_up)[1],dim(dreps2_up)[1],dim(dreps3_up)[1],dim(dreps4_up)[1],dim(dreps5_up)[1],dim(dreps6_up)[1],dim(dreps7_up)[1],dim(dreps8_up)[1])
downdgs <- c(1,dim(dgs1_down)[1],dim(dgs2_down)[1],dim(dgs3_down)[1],dim(dgs4_down)[1],dim(dgs5_down)[1],dim(dgs6_down)[1],dim(dgs7_down)[1],dim(dgs8_down)[1])
downreps <- c(1,dim(dreps1_down)[1],dim(dreps2_down)[1],dim(dreps3_down)[1],dim(dreps4_down)[1],dim(dreps5_down)[1],dim(dreps6_down)[1],dim(dreps7_down)[1],dim(dreps8_down)[1])

write.table(dgs1,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/featurecounts/results_flexthres/DEgenes_IAV3h_FDR0.25.txt",sep="\t",quote=FALSE)
write.table(dreps1,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/featurecounts/results_flexthres/DEreps_IAV3h_FDR0.25.txt",sep="\t",quote=FALSE)
write.table(dgs2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/featurecounts/results_flexthres/DEgenes_IAV11h_FDR0.25.txt",sep="\t",quote=FALSE)
write.table(dreps2,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/featurecounts/results_flexthres/DEreps_IAV11h_FDR0.25.txt",sep="\t",quote=FALSE)
write.table(dgs3,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/featurecounts/results_flexthres/DEgenes_IAV26h_FDR0.25.txt",sep="\t",quote=FALSE)
write.table(dreps3,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/featurecounts/results_flexthres/DEreps_IAV26h_FDR0.25.txt",sep="\t",quote=FALSE)
write.table(dgs4,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/featurecounts/results_flexthres/DEgenes_IAV49h_FDR0.25.txt",sep="\t",quote=FALSE)
write.table(dreps4,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/featurecounts/results_flexthres/DEreps_IAV49h_FDR0.25.txt",sep="\t",quote=FALSE)
write.table(dgs5,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/featurecounts/results_flexthres/DEgenes_IAV74h_FDR0.25.txt",sep="\t",quote=FALSE)
write.table(dreps5,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/featurecounts/results_flexthres/DEreps_IAV74h_FDR0.25.txt",sep="\t",quote=FALSE)
write.table(dgs6,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/featurecounts/results_flexthres/DEgenes_IAV99h_FDR0.25.txt",sep="\t",quote=FALSE)
write.table(dreps6,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/featurecounts/results_flexthres/DEreps_IAV99h_FDR0.25.txt",sep="\t",quote=FALSE)
write.table(dgs7,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/featurecounts/results_flexthres/DEgenes_IAV122h_FDR0.25.txt",sep="\t",quote=FALSE)
write.table(dreps7,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/featurecounts/results_flexthres/DEreps_IAV122h_FDR0.25.txt",sep="\t",quote=FALSE)
write.table(dgs8,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/featurecounts/results_flexthres/DEgenes_IAV149h_FDR0.25.txt",sep="\t",quote=FALSE)
write.table(dreps8,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/featurecounts/results_flexthres/DEreps_IAV149h_FDR0.25.txt",sep="\t",quote=FALSE)

save.image("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/featurecounts/IAV_timecourse_featurecounts_looserthres.Rda")
