##########################################################################################
##########################################################################################
##
## Murine Norovirus (MNV) time course in RAW 264.7 cells
##
##########################################################################################
##########################################################################################

fullmaster <- read.delim("~/Desktop/TE_project/mouse_encode/counts/ENCODE-processed-totalRNA/MM_EP_STAR_FC_DFAM_finalresults/master.te.to.gene.table.FULL.exonannot2.txt",skip=1,header=FALSE)
genes <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/featurecounts/MNV_timecourse_RAW_gene_counts.mx",sep="\t", skip=1, header = T, row.names=1)
reps <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/featurecounts/MNV_timecourse_RAW_repeat_counts.mx",sep="\t", skip=1, header = T, row.names=1)

# Reformat the sample IDs in gene counts matrix
colnames(genes) <- sub("X.home.mmacchie.mouse.encode.virus.MNV_timecourse_RAW.bams.", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(genes))
colnames(genes) <- sub("bams_no_virus_UCSC.","",colnames(genes))

# Reformat the sample IDs in repeat counts matrix
colnames(reps) <- sub("X.home.mmacchie.mouse.encode.virus.MNV_timecourse_RAW.bams.", "", colnames(reps))
colnames(reps) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(reps))
colnames(reps) <- sub("bams_no_virus_UCSC.","",colnames(reps))

# Remove Simple Repeats from the repeat counts matrix
SRs <- fullmaster[fullmaster$V12 =="Simple_repeat",]$V10
reps <- reps[!rownames(reps) %in% SRs,]
reps <- na.omit(reps)
gene.anno <- genes[, c(1:5)]
genes <- genes[, c(6:dim(genes)[2])]
rep.anno <- reps[, c(1:5)]
reps <- reps[, c(6:dim(reps)[2])]

# Combine the repeat and genes counts matrix
genes <- rbind(genes,reps)
library(edgeR)

# Normalize repeat and gene counts matrix by the total number of mapped and quantified reads
cpms <- cpm(genes)

# Separate normalized virus counts from genes and repeat counts
virus <- subset(cpms,grepl("^MNV1",rownames(cpms)))
virus <- colSums(virus)
genes <- subset(genes,!grepl("^MNV1",rownames(genes)))

# reorder the samples in the gene and repeat raw counts matrix ("genes"), the normalized counts matrix ("cpms"), and the virus normalized counts matrix ("virus")
neworder = c("mock_0h_1","mock_0h_2","inf_0h_1","inf_0h_2","inf_0h_3","inf_8h_1","inf_8h_2","inf_8h_3","inf_14h_1","inf_14h_2","inf_14h_3","inf_20h_1","inf_20h_2","inf_20h_3")
cpms <- cpms[,match(neworder,colnames(cpms))]
genes <- genes[,match(neworder,colnames(genes))]
virus <- virus[match(neworder,colnames(virus))]

# create subsetted raw counts matrices for each comparison - each infection time point compared to 0-hour
e1 <- genes[,c(3:5,1:2)] #0
e2 <- genes[,c(6:8,1:2)] #8
e3 <- genes[,c(9:11,1:2)] #14
e4 <- genes[,c(12:14,1:2)] #20


# Compare 0h infected vs 0 hour uninfected
library(edgeR)
y <- DGEList(counts=e1,group=c(1,1,1,2,2))
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table
tt <- tt[tt$FDR < 0.05 ,]
dgs1 <- subset(tt,grepl('^EN',rownames(tt)))
dreps1 <- subset(tt,grepl('^DF',rownames(tt)))
dim(dgs1)
dim(dreps1)
dgs1_up <- dgs1[dgs1$logFC < 0,]
dgs1_down <- dgs1[dgs1$logFC > 0,]
dreps1_up <- dreps1[dreps1$logFC < 0,]
dreps1_down <- dreps1[dreps1$logFC > 0,]

# Compare 8h infected vs 0 hour uninfected
y <- DGEList(counts=e2,group=c(1,1,1,2,2))
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table
tt <- tt[tt$FDR < 0.05 ,]
dgs2 <- subset(tt,grepl('^EN',rownames(tt)))
dreps2 <- subset(tt,grepl('^DF',rownames(tt)))
dim(dgs2)
dim(dreps2)
dgs2_up <- dgs2[dgs2$logFC < 0,]
dgs2_down <- dgs2[dgs2$logFC > 0,]
dreps2_up <- dreps2[dreps2$logFC < 0,]
dreps2_down <- dreps2[dreps2$logFC > 0,]

# Compare 14h infected vs 0 hour uninfected
y <- DGEList(counts=e3,group=c(1,1,1,2,2))
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table
tt <- tt[tt$FDR < 0.05 ,]
dgs3 <- subset(tt,grepl('^EN',rownames(tt)))
dreps3 <- subset(tt,grepl('^DF',rownames(tt)))
dim(dgs3)
dim(dreps3)
dgs3_up <- dgs3[dgs3$logFC < 0,]
dgs3_down <- dgs3[dgs3$logFC > 0,]
dreps3_up <- dreps3[dreps3$logFC < 0,]
dreps3_down <- dreps3[dreps3$logFC > 0,]

# Compare 20h infected vs 0 hour uninfected
library(edgeR)
y <- DGEList(counts=e4,group=c(1,1,1,2,2))
keep <- rowSums(cpm(y)>1) >= 2
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(genes)[1])
tt <- tt$table
tt <- tt[tt$FDR < 0.05 ,]
dgs4 <- subset(tt,grepl('^EN',rownames(tt)))
dreps4 <- subset(tt,grepl('^DF',rownames(tt)))
dim(dgs4)
dim(dreps4)
dgs4_up <- dgs4[dgs4$logFC < 0,]
dgs4_down <- dgs4[dgs4$logFC > 0,]
dreps4_up <- dreps4[dreps4$logFC < 0,]
dreps4_down <- dreps4[dreps4$logFC > 0,]


# Create a list of TE IDs (specifically ERVs and LINEs)
TEs <- subset(fullmaster,grepl("^L",fullmaster$V12))
TEs <- subset(TEs,!grepl("^Low",TEs$V12))
TEs <- unique(TEs$V10)

# Extract TEs from DE repeats from each comparison
dreps1_up_TEsonly <- dreps1_up[rownames(dreps1_up) %in% TEs,]
dreps2_up_TEsonly <- dreps2_up[rownames(dreps2_up) %in% TEs,]
dreps3_up_TEsonly <- dreps3_up[rownames(dreps3_up) %in% TEs,]
dreps4_up_TEsonly <- dreps4_up[rownames(dreps4_up) %in% TEs,]
dreps1_down_TEsonly <- dreps1_down[rownames(dreps1_down) %in% TEs,]
dreps2_down_TEsonly <- dreps2_down[rownames(dreps2_down) %in% TEs,]
dreps3_down_TEsonly <- dreps3_down[rownames(dreps3_down) %in% TEs,]
dreps4_down_TEsonly <- dreps4_down[rownames(dreps4_down) %in% TEs,]

# Write out bed file coordinates of all MNV DE TEs for GREAT
dreps_TEs_all <- rbind(dreps1_up_TEsonly,dreps2_up_TEsonly,dreps3_up_TEsonly,dreps4_up_TEsonly,dreps1_down_TEsonly,dreps2_down_TEsonly,dreps3_down_TEsonly,dreps4_down_TEsonly)
dreps_TEs_all <- unique(rownames(dreps_TEs_all))
dreps_TEs_all <- unique(fullmaster[fullmaster$V10 %in% dreps_TEs_all,c(1,4,5,7,12,10)])
write.table(dreps_TEs_all,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/GREAT/TE.coordinates.MNV.timecourse.for.GREAT.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)



totdgs <- c(dim(dgs1)[1],dim(dgs2)[1],dim(dgs3)[1],dim(dgs4)[1])
totdreps <- c(dim(dreps1)[1],dim(dreps2)[1],dim(dreps3)[1],dim(dreps4)[1])

updgs <- c(dim(dgs1_up)[1],dim(dgs2_up)[1],dim(dgs3_up)[1],dim(dgs4_up)[1])
upreps <- c(dim(dreps1_up)[1],dim(dreps2_up)[1],dim(dreps3_up)[1],dim(dreps4_up)[1])
downdgs <- c(dim(dgs1_down)[1],dim(dgs2_down)[1],dim(dgs3_down)[1],dim(dgs4_down)[1])
downreps <- c(dim(dreps1_down)[1],dim(dreps2_down)[1],dim(dreps3_down)[1],dim(dreps4_down)[1])
upreps_TEsonly <- c(dim(dreps1_up_TEsonly)[1],dim(dreps2_up_TEsonly)[1],dim(dreps3_up_TEsonly)[1],dim(dreps4_up_TEsonly)[1])
downreps_TEsonly <- c(dim(dreps1_down_TEsonly)[1],dim(dreps2_down_TEsonly)[1],dim(dreps3_down_TEsonly)[1],dim(dreps4_down_TEsonly)[1])


time2 <- c(0,8,14,20)
mock <- mean(virus[1],virus[2])
v <- c(mean(virus[3],virus[4],virus[5]),mean(virus[6],virus[7],virus[8]),mean(virus[9],virus[10],virus[11]),mean(virus[12],virus[13],virus[14]))
v <- v/mock

intB <- cpms[rownames(cpms) == "ENSMUSG00000048806",]
intB <- c(mean(intB[3],intB[4],intB[5]),mean(intB[6],intB[7],intB[8]),mean(intB[9],intB[10],intB[11]),mean(intB[12],intB[13],intB[14]))


#--------- normalized total DE genes and repeats
plot(time2,totdgs/sum(totdgs),pch=19,xaxt="n",cex=0.5,col="red",ylab="% DE genes and repeats",xlab="Time (h)",ylim=c(0,0.7))#,ylim=c(0,175))#,log="y")
lines(time2,totdgs/sum(totdgs),col="red")
#lines(time2,upreps/sum(upreps),col="red",lty=2)
#lines(time2,downreps/sum(downreps),col="blue",lty=2)
lines(time2,totdreps/sum(totdreps),col="darkgreen")
lines(time2,v/sum(v),col="black")
lines(time2,intB/sum(intB),col="orange")
axis(side=1, at=time2, labels = TRUE)
abline(v=0,lty=2,col="grey")
abline(v=8,lty=2,col="grey")
abline(v=14,lty=2,col="grey")
abline(v=20,lty=2,col="grey")


#---------  absolute total # DE UP AND DOWN genes and repeats -- dashed line = down, solid line = up; red line= genes; green = repeats
plot(time2,downreps,pch=19,xaxt="n",cex=0.5,col="darkgreen",ylab="number of DE genes and repeats",xlab="Time (h)",ylim=c(0,5200))
lines(time2,downreps,col="darkgreen",lty=2)
lines(time2,upreps,col="darkgreen")
lines(time2,updgs,col="red")
lines(time2,downdgs,col="red",lty=2)
axis(side=1, at=time2, labels = TRUE)
abline(v=0,lty=2,col="grey")
abline(v=8,lty=2,col="grey")
abline(v=14,lty=2,col="grey")
abline(v=20,lty=2,col="grey")


#---------  % DE UP AND DOWN genes and ERVs and LINEs -- dashed line = down, solid line = up; red line= genes; green = repeats
plot(time2,downreps_TEsonly/sum(downreps_TEsonly),pch=19,xaxt="n",cex=0.5,col="darkgreen",ylab="% DE genes and TEs",xlab="Time (h)",ylim=c(0,0.7))
lines(time2,downreps_TEsonly/sum(downreps_TEsonly),col="darkgreen",lty=2)
lines(time2,upreps_TEsonly/sum(upreps_TEsonly),col="darkgreen")
lines(time2,updgs/sum(updgs),col="red")
lines(time2,downdgs/sum(downdgs),col="red",lty=2)
axis(side=1, at=time2, labels = TRUE)
abline(v=0,lty=2,col="grey")
abline(v=8,lty=2,col="grey")
abline(v=14,lty=2,col="grey")
abline(v=20,lty=2,col="grey")
lines(time2,v/sum(v),col="black")
lines(time2,intB/sum(intB),col="orange")

save.image("MNV_timecourse_RAWcells_FDR0.05.Rda")
