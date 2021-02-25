##--------------- IAV time course - Brisbane -------------------------#

load("IAV_timecourse_Brisbane.Rda")
fullmaster <- read.delim("/Users/mmacchie/Desktop/TE_project/human_encode/annot/master.te.to.gene.table.hg38.exonannot.txt",header=FALSE)
genes <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/featurecounts/IAV_Brisbane_gene_counts_virus.mx",sep="\t", skip=1, header = T)
reps <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/featurecounts/IAV_Brisbane_repeat_counts2_virus.mx",sep="\t", skip=1, header = T, row.names=1)

colnames(genes) <- sub("X.home.mmacchie.human.encode.infection_TEtrans.IAV_Brisbane.bams.", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(genes))

colnames(reps) <- sub("X.home.mmacchie.human.encode.infection_TEtrans.IAV_Brisbane.bams.", "", colnames(reps))
colnames(reps) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(reps))


gene.anno <- genes[, c(2:6)]
genes <- genes[, c(1,7:dim(genes)[2])]
virus <- subset(genes,!grepl("^EN",genes$Geneid))
genes <- subset(genes,grepl("^EN",genes$Geneid))

rownames(genes) <- genes$Geneid
genes <- genes[,-1]
virus <- virus[,-1]
virus <- colSums(virus)
genes_virus <- rbind(virus,genes)
cpms <- cpm(genes_virus)
cpms_virus <- cpms[1,]
neworder <- c("inf_1_1h","inf_2_1h","inf_3_1h","inf_1_6h","inf_2_6h","inf_3_6h","inf_1_24h","inf_2_24h","inf_3_24h")
cpms_virus <- cpms_virus[match(neworder,names(cpms_virus))]
neworder <- c("mock_1_1h","mock_2_1h","mock_3_1h","mock_1_6h","mock_2_6h","mock_3_6h","mock_1_24h","mock_2_24h","mock_3_24h")
cpms_virus_mock <- cpms[1,]
cpms_virus_mock <- cpms_virus_mock[match(neworder,names(cpms_virus_mock))]


SRs <- fullmaster[fullmaster$V12 =="Simple_repeat",]$V10
reps <- reps[!rownames(reps) %in% SRs,]
reps <- na.omit(reps)
#gene.anno <- genes[, c(1:5)]
#genes <- genes[, c(6:dim(genes)[2])]
rep.anno <- reps[, c(1:5)]
reps <- reps[, c(6:dim(reps)[2])]


#anno <- rbind(gene.anno,rep.anno)
genes <- rbind(genes,reps)
library(edgeR)
#rpkms <- rpkm(genes,anno$Length)

cpms <- cpm(genes)



neworder = c("inf_1_1h","inf_2_1h","inf_3_1h","mock_1_1h","mock_2_1h","mock_3_1h")
e1 <- genes[,match(neworder,colnames(genes))] #1
neworder = c("inf_1_6h","inf_2_6h","inf_3_6h","mock_1_6h","mock_2_6h","mock_3_6h")
e2 <- genes[,match(neworder,colnames(genes))] #3
neworder = c("inf_1_24h","inf_2_24h","inf_3_24h","mock_1_24h","mock_2_24h","mock_3_24h")
e3 <- genes[,match(neworder,colnames(genes))] #24


library(edgeR)
y <- DGEList(counts=e1,group=c(1,1,1,2,2,2))
keep <- rowSums(cpm(y)>1) >= 3
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




library(edgeR)
y <- DGEList(counts=e2,group=c(1,1,1,2,2,2))
keep <- rowSums(cpm(y)>1) >= 3
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


library(edgeR)
y <- DGEList(counts=e3,group=c(1,1,1,2,2,2))
keep <- rowSums(cpm(y)>1) >= 3
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



TEs <- subset(fullmaster,grepl("^L",fullmaster$V12))
TEs <- subset(TEs,!grepl("^Low",TEs$V12))
TEs <- unique(TEs$V10)

dreps1_up_TEsonly <- dreps1_up[rownames(dreps1_up) %in% TEs,]
dreps2_up_TEsonly <- dreps2_up[rownames(dreps2_up) %in% TEs,]
dreps3_up_TEsonly <- dreps3_up[rownames(dreps3_up) %in% TEs,]


dreps1_down_TEsonly <- dreps1_down[rownames(dreps1_down) %in% TEs,]
dreps2_down_TEsonly <- dreps2_down[rownames(dreps2_down) %in% TEs,]
dreps3_down_TEsonly <- dreps3_down[rownames(dreps3_down) %in% TEs,]

dreps_TEs_all <- rbind(dreps1_up_TEsonly,dreps2_up_TEsonly,dreps3_up_TEsonly,dreps1_down_TEsonly,dreps2_down_TEsonly,dreps3_down_TEsonly)
dreps_TEs_all <- unique(rownames(dreps_TEs_all))
dreps_TEs_all <- unique(fullmaster[fullmaster$V10 %in% dreps_TEs_all,c(1,4,5,7,12,10)])
write.table(dreps_TEs_all,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/GREAT/TE.coordinates.IAVbrisbane.timecourse.for.GREAT.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

fullmaster$V1 <- paste("chr",fullmaster$V1,sep='')
test<- unique(fullmaster[fullmaster$V10 %in% rownames(dreps1_up_TEsonly),c(1,4,5,7,10)])
write.table(test,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/MHC_analysis/IAV/TE.coordinates.IAVbrisbane.timecourse.for.MHCregion.0hUP.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
test<- unique(fullmaster[fullmaster$V10 %in% rownames(dreps1_down_TEsonly),c(1,4,5,7,10)])
write.table(test,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/MHC_analysis/IAV/TE.coordinates.IAVbrisbane.timecourse.for.MHCregion.0hDOWN.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
test<- unique(fullmaster[fullmaster$V10 %in% rownames(dreps2_up_TEsonly),c(1,4,5,7,10)])
write.table(test,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/MHC_analysis/IAV/TE.coordinates.IAVbrisbane.timecourse.for.MHCregion.8hUP.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
test<- unique(fullmaster[fullmaster$V10 %in% rownames(dreps2_down_TEsonly),c(1,4,5,7,10)])
write.table(test,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/MHC_analysis/IAV/TE.coordinates.IAVbrisbane.timecourse.for.MHCregion.8hDOWN.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
test<- unique(fullmaster[fullmaster$V10 %in% rownames(dreps3_up_TEsonly),c(1,4,5,7,10)])
write.table(test,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/MHC_analysis/IAV/TE.coordinates.IAVbrisbane.timecourse.for.MHCregion.14hUP.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
test<- unique(fullmaster[fullmaster$V10 %in% rownames(dreps3_down_TEsonly),c(1,4,5,7,10)])
write.table(test,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/MHC_analysis/IAV/TE.coordinates.IAVbrisbane.timecourse.for.MHCregion.14hDOWN.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)


save.image("IAV_timecourse_Brisbane.Rda")

# write out TEs changing at 0 hours
test <- dreps1[rownames(dreps1) %in% TEs,]
test2 <- fullmaster[fullmaster$V10 %in% rownames(test),c(1,4,5,10,7,12,13,14,15)]
test3 <- test[match(test2$V10,rownames(test)),]
test <- cbind(test2,test3)
write.table(test,"IAV.brisbane.DETEs.0h.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)


# write out TEs changing at 6 hours
test <- dreps2[rownames(dreps2) %in% TEs,]
test2 <- fullmaster[fullmaster$V10 %in% rownames(test),c(1,4,5,10,7,12,13,14,15)]
test3 <- test[match(test2$V10,rownames(test)),]
test <- cbind(test2,test3)
write.table(test,"IAV.brisbane.DETEs.6h.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

# write out TEs changing at 24 hours
test <- dreps3[rownames(dreps3) %in% TEs,]
test2 <- fullmaster[fullmaster$V10 %in% rownames(test),c(1,4,5,10,7,12,13,14,15)]
test3 <- test[match(test2$V10,rownames(test)),]
test <- cbind(test2,test3)
write.table(test,"IAV.brisbane.DETEs.24h.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)




totdgs <- c(dim(dgs1)[1],dim(dgs2)[1],dim(dgs3)[1])
totdreps <- c(dim(dreps1)[1],dim(dreps2)[1],dim(dreps3)[1])

updgs <- c(dim(dgs1_up)[1],dim(dgs2_up)[1],dim(dgs3_up)[1])
upreps <- c(dim(dreps1_up)[1],dim(dreps2_up)[1],dim(dreps3_up)[1])
downdgs <- c(dim(dgs1_down)[1],dim(dgs2_down)[1],dim(dgs3_down)[1])
downreps <- c(dim(dreps1_down)[1],dim(dreps2_down)[1],dim(dreps3_down)[1])
upreps_TEsonly <- c(dim(dreps1_up_TEsonly)[1],dim(dreps2_up_TEsonly)[1],dim(dreps3_up_TEsonly)[1])
downreps_TEsonly <- c(dim(dreps1_down_TEsonly)[1],dim(dreps2_down_TEsonly)[1],dim(dreps3_down_TEsonly)[1])
dreps_TEsonly <- colSums(rbind(upreps_TEsonly,downreps_TEsonly))






time2 <- c(0,6,24)
mock <- c(mean(cpms_virus_mock[1],cpms_virus_mock[2],cpms_virus_mock[3]),mean(cpms_virus_mock[4],cpms_virus_mock[5],cpms_virus_mock[6]),mean(cpms_virus_mock[7],cpms_virus_mock[8],cpms_virus_mock[9]))
v <- c(mean(cpms_virus[1],cpms_virus[2],cpms_virus[3]),mean(cpms_virus[4],cpms_virus[5],cpms_virus[6]),mean(cpms_virus[7],cpms_virus[8],cpms_virus[9]))
v <- v/(mock+1)


intB_e1 <- cpms[rownames(cpms) == "ENSG00000171855",match(colnames(e1),colnames(cpms))]
intB_e2 <- cpms[rownames(cpms) == "ENSG00000171855",match(colnames(e2),colnames(cpms))]
intB_e3 <- cpms[rownames(cpms) == "ENSG00000171855",match(colnames(e3),colnames(cpms))]
intB <- c(mean(intB_e1[1],intB_e1[2],intB_e1[3]),mean(intB_e2[1],intB_e2[2],intB_e2[3]),mean(intB_e3[1],intB_e3[2],intB_e3[3]))
intB <- intB/sum(intB)



#--------- normalized total DE genes and repeats
plot(time2,totdgs/sum(totdgs),pch=19,xaxt="n",cex=0.5,col="red",ylab="% DE genes and repeats",xlab="Time (h)",ylim=c(0,0.9))#,ylim=c(0,175))#,log="y")
lines(time2,totdgs/sum(totdgs),col="red")
#lines(time2,upreps/sum(upreps),col="red",lty=2)
#lines(time2,downreps/sum(downreps),col="blue",lty=2)
#lines(time2,totdreps/sum(totdreps),col="darkgreen")
lines(time2,v/sum(v),col="black")
lines(time2,intB/sum(intB),col="orange")
axis(side=1, at=time2, labels = TRUE)
abline(v=0,lty=2,col="grey")
abline(v=6,lty=2,col="grey")
abline(v=24,lty=2,col="grey")


#add normalized down and up DE genes and reps
lines(time2,downreps/sum(downreps),col="blue")
lines(time2,upreps/sum(upreps),col="red")
lines(time2,downdgs/sum(downdgs),col="blue",lty=2)
lines(time2,updgs/sum(updgs),col="red",lty=2)



#--------- normalized total DE genes and TEs
plot(time2,totdgs/sum(totdgs),pch=19,xaxt="n",cex=0.5,col="red",ylab="% DE genes and TEs",xlab="Time (h)",ylim=c(0,0.9))#,ylim=c(0,175))#,log="y")
lines(time2,totdgs/sum(totdgs),col="red")
#lines(time2,updgs/sum(updgs),col="red")
#lines(time2,downdgs/sum(downdgs),col="red",lty=2)
#lines(time2,upreps/sum(upreps),col="red",lty=2)
#lines(time2,downreps/sum(downreps),col="blue",lty=2)
lines(time2,dreps_TEsonly/sum(dreps_TEsonly),col="darkgreen")
#lines(time2,upreps_TEsonly/sum(upreps_TEsonly),col="darkgreen")
#lines(time2,downreps_TEsonly/sum(downreps_TEsonly),col="darkgreen",lty=2)
lines(time2,v/sum(v),col="black")
lines(time2,intB/sum(intB),col="orange")
axis(side=1, at=time2, labels = TRUE)
abline(v=0,lty=2,col="grey")
abline(v=6,lty=2,col="grey")
abline(v=24,lty=2,col="grey")



#---------  absolute total # DE genes and repeats
plot(time2,totdreps,pch=19,xaxt="n",cex=0.5,col="darkgreen",ylab="number of DE genes and repeats",xlab="Time (h)")
lines(time2,totdreps,col="darkgreen")
lines(time2,totdgs,col="red")
axis(side=1, at=time2, labels = TRUE)
abline(v=0,lty=2,col="grey")
abline(v=6,lty=2,col="grey")

abline(v=24,lty=2,col="grey")

#---------  absolute total # DE UP AND DOWN genes and repeats -- dashed line = up, solid line = down
plot(time2,downreps,pch=19,xaxt="n",cex=0.5,col="darkgreen",ylab="number of DE genes and repeats",xlab="Time (h)",ylim=c(0,1500))
lines(time2,downreps,col="darkgreen",lty=2)
lines(time2,upreps,col="darkgreen")
lines(time2,updgs,col="red")
lines(time2,downdgs,col="red",lty=2)
axis(side=1, at=time2, labels = TRUE)
abline(v=0,lty=2,col="grey")
abline(v=6,lty=2,col="grey")
abline(v=20,lty=2,col="grey")









##--------------- IAV time course - Udorn -------------------------#

load("IAV_timecourse_Udorn.Rda")
fullmaster <- read.delim("/Users/mmacchie/Desktop/TE_project/human_encode/annot/master.te.to.gene.table.hg38.exonannot.txt",header=FALSE)
genes <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/featurecounts/IAV_Udorn_gene_counts_virus.mx",sep="\t", skip=1, header = T)
reps <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/featurecounts/IAV_Udorn_repeat_counts2_virus.mx",sep="\t", skip=1, header = T, row.names=1)

colnames(genes) <- sub("X.home.mmacchie.human.encode.infection_TEtrans.IAV_Udorn.bams.", "", colnames(genes))
colnames(genes) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(genes))

colnames(reps) <- sub("X.home.mmacchie.human.encode.infection_TEtrans.IAV_Udorn.bams.", "", colnames(reps))
colnames(reps) <- sub("Aligned.sortedByCoord.out.bam", "", colnames(reps))


gene.anno <- genes[, c(2:6)]
genes <- genes[, c(1,7:dim(genes)[2])]
virus <- subset(genes,!grepl("^EN",genes$Geneid))
genes <- subset(genes,grepl("^EN",genes$Geneid))

rownames(genes) <- genes$Geneid
genes <- genes[,-1]
virus <- virus[,-1]
virus <- colSums(virus)
genes_virus <- rbind(virus,genes)
cpms <- cpm(genes_virus)
cpms_virus <- cpms[1,]
neworder <- c("inf_1_1h","inf_2_1h","inf_3_1h","inf_1_6h","inf_2_6h","inf_3_6h","inf_1_24h","inf_2_24h","inf_3_24h")
cpms_virus <- cpms_virus[match(neworder,names(cpms_virus))]
neworder <- c("mock_1_1h","mock_2_1h","mock_3_1h","mock_1_6h","mock_2_6h","mock_3_6h","mock_1_24h","mock_2_24h","mock_3_24h")
cpms_virus_mock <- cpms[1,]
cpms_virus_mock <- cpms_virus_mock[match(neworder,names(cpms_virus_mock))]


SRs <- fullmaster[fullmaster$V12 =="Simple_repeat",]$V10
reps <- reps[!rownames(reps) %in% SRs,]
reps <- na.omit(reps)
rep.anno <- reps[, c(1:5)]
reps <- reps[, c(6:dim(reps)[2])]


genes <- rbind(genes,reps)
library(edgeR)


cpms <- cpm(genes)


neworder = c("inf_1_1h","inf_2_1h","inf_3_1h","mock_1_1h","mock_2_1h","mock_3_1h")
e1 <- genes[,match(neworder,colnames(genes))] #1
neworder = c("inf_1_6h","inf_2_6h","inf_3_6h","mock_1_6h","mock_2_6h","mock_3_6h")
e2 <- genes[,match(neworder,colnames(genes))] #3
neworder = c("inf_1_24h","inf_2_24h","inf_3_24h","mock_1_24h","mock_2_24h","mock_3_24h")
e3 <- genes[,match(neworder,colnames(genes))] #24


library(edgeR)
y <- DGEList(counts=e1,group=c(1,1,1,2,2,2))
keep <- rowSums(cpm(y)>1) >= 3
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
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




library(edgeR)
y <- DGEList(counts=e2,group=c(1,1,1,2,2,2))
keep <- rowSums(cpm(y)>1) >= 3
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
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


library(edgeR)
y <- DGEList(counts=e3,group=c(1,1,1,2,2,2))
keep <- rowSums(cpm(y)>1) >= 3
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
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


TEs <- subset(fullmaster,grepl("^L",fullmaster$V12))
TEs <- subset(TEs,!grepl("^Low",TEs$V12))
TEs <- unique(TEs$V10)

dreps1_up_TEsonly <- dreps1_up[rownames(dreps1_up) %in% TEs,]
dreps2_up_TEsonly <- dreps2_up[rownames(dreps2_up) %in% TEs,]
dreps3_up_TEsonly <- dreps3_up[rownames(dreps3_up) %in% TEs,]


dreps1_down_TEsonly <- dreps1_down[rownames(dreps1_down) %in% TEs,]
dreps2_down_TEsonly <- dreps2_down[rownames(dreps2_down) %in% TEs,]
dreps3_down_TEsonly <- dreps3_down[rownames(dreps3_down) %in% TEs,]

dreps_TEs_all <- rbind(dreps1_up_TEsonly,dreps2_up_TEsonly,dreps3_up_TEsonly,dreps1_down_TEsonly,dreps2_down_TEsonly,dreps3_down_TEsonly)
dreps_TEs_all <- unique(rownames(dreps_TEs_all))
dreps_TEs_all <- unique(fullmaster[fullmaster$V10 %in% dreps_TEs_all,c(1,4,5,7,12,10)])
write.table(dreps_TEs_all,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/GREAT/TE.coordinates.IAVUdorn.timecourse.for.GREAT.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

fullmaster$V1 <- paste("chr",fullmaster$V1,sep='')
test<- unique(fullmaster[fullmaster$V10 %in% rownames(dreps1_up_TEsonly),c(1,4,5,7,10)])
write.table(test,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/MHC_analysis/IAV/TE.coordinates.IAVUdorn.timecourse.for.MHCregion.0hUP.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
test<- unique(fullmaster[fullmaster$V10 %in% rownames(dreps1_down_TEsonly),c(1,4,5,7,10)])
write.table(test,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/MHC_analysis/IAV/TE.coordinates.IAVUdorn.timecourse.for.MHCregion.0hDOWN.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
test<- unique(fullmaster[fullmaster$V10 %in% rownames(dreps2_up_TEsonly),c(1,4,5,7,10)])
write.table(test,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/MHC_analysis/IAV/TE.coordinates.IAVUdorn.timecourse.for.MHCregion.8hUP.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
test<- unique(fullmaster[fullmaster$V10 %in% rownames(dreps2_down_TEsonly),c(1,4,5,7,10)])
write.table(test,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/MHC_analysis/IAV/TE.coordinates.IAVUdorn.timecourse.for.MHCregion.8hDOWN.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
test<- unique(fullmaster[fullmaster$V10 %in% rownames(dreps3_up_TEsonly),c(1,4,5,7,10)])
write.table(test,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/MHC_analysis/IAV/TE.coordinates.IAVUdorn.timecourse.for.MHCregion.14hUP.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
test<- unique(fullmaster[fullmaster$V10 %in% rownames(dreps3_down_TEsonly),c(1,4,5,7,10)])
write.table(test,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/MHC_analysis/IAV/TE.coordinates.IAVUdorn.timecourse.for.MHCregion.14hDOWN.bed",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)


save.image("IAV_timecourse_Udorn.Rda")


# write out TEs changing at 0 hours
test <- dreps1[rownames(dreps1) %in% TEs,]
test2 <- fullmaster[fullmaster$V10 %in% rownames(test),c(1,4,5,10,7,12,13,14,15)]
test3 <- test[match(test2$V10,rownames(test)),]
test <- cbind(test2,test3)
write.table(test,"IAV.Udorn.DETEs.0h.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)


# write out TEs changing at 6 hours
test <- dreps2[rownames(dreps2) %in% TEs,]
test2 <- fullmaster[fullmaster$V10 %in% rownames(test),c(1,4,5,10,7,12,13,14,15)]
test3 <- test[match(test2$V10,rownames(test)),]
test <- cbind(test2,test3)
write.table(test,"IAV.Udorn.DETEs.6h.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

# write out TEs changing at 24 hours
test <- dreps3[rownames(dreps3) %in% TEs,]
test2 <- fullmaster[fullmaster$V10 %in% rownames(test),c(1,4,5,10,7,12,13,14,15)]
test3 <- test[match(test2$V10,rownames(test)),]
test <- cbind(test2,test3)
write.table(test,"IAV.Udorn.DETEs.24h.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)




totdgs <- c(dim(dgs1)[1],dim(dgs2)[1],dim(dgs3)[1])
totdreps <- c(dim(dreps1)[1],dim(dreps2)[1],dim(dreps3)[1])

updgs <- c(dim(dgs1_up)[1],dim(dgs2_up)[1],dim(dgs3_up)[1])
upreps <- c(dim(dreps1_up)[1],dim(dreps2_up)[1],dim(dreps3_up)[1])
downdgs <- c(dim(dgs1_down)[1],dim(dgs2_down)[1],dim(dgs3_down)[1])
downreps <- c(dim(dreps1_down)[1],dim(dreps2_down)[1],dim(dreps3_down)[1])
upreps_TEsonly <- c(dim(dreps1_up_TEsonly)[1],dim(dreps2_up_TEsonly)[1],dim(dreps3_up_TEsonly)[1])
downreps_TEsonly <- c(dim(dreps1_down_TEsonly)[1],dim(dreps2_down_TEsonly)[1],dim(dreps3_down_TEsonly)[1])
dreps_TEsonly <- colSums(rbind(upreps_TEsonly,downreps_TEsonly))


time2 <- c(0,6,24)
mock <- c(mean(cpms_virus_mock[1],cpms_virus_mock[2],cpms_virus_mock[3]),mean(cpms_virus_mock[4],cpms_virus_mock[5],cpms_virus_mock[6]),mean(cpms_virus_mock[7],cpms_virus_mock[8],cpms_virus_mock[9]))
v <- c(mean(cpms_virus[1],cpms_virus[2],cpms_virus[3]),mean(cpms_virus[4],cpms_virus[5],cpms_virus[6]),mean(cpms_virus[7],cpms_virus[8],cpms_virus[9]))
v <- v/(mock+1)


intB_e1 <- cpms[rownames(cpms) == "ENSG00000171855",match(colnames(e1),colnames(cpms))]
intB_e2 <- cpms[rownames(cpms) == "ENSG00000171855",match(colnames(e2),colnames(cpms))]
intB_e3 <- cpms[rownames(cpms) == "ENSG00000171855",match(colnames(e3),colnames(cpms))]
intB <- c(mean(intB_e1[1],intB_e1[2],intB_e1[3]),mean(intB_e2[1],intB_e2[2],intB_e2[3]),mean(intB_e3[1],intB_e3[2],intB_e3[3]))
intB <- intB/sum(intB)



#--------- normalized total DE genes and repeats
plot(time2,totdgs/sum(totdgs),pch=19,xaxt="n",cex=0.5,col="red",ylab="% DE genes and repeats",xlab="Time (h)",ylim=c(0,0.9))#,ylim=c(0,175))#,log="y")
lines(time2,totdgs/sum(totdgs),col="red")
#lines(time2,upreps/sum(upreps),col="red",lty=2)
#lines(time2,downreps/sum(downreps),col="blue",lty=2)
lines(time2,totdreps/sum(totdreps),col="darkgreen")
lines(time2,v/sum(v),col="black")
lines(time2,intB/sum(intB),col="orange")
axis(side=1, at=time2, labels = TRUE)
abline(v=0,lty=2,col="grey")
abline(v=6,lty=2,col="grey")
abline(v=24,lty=2,col="grey")


#add normalized down and up DE genes and reps
lines(time2,downreps/sum(downreps),col="blue")
lines(time2,upreps/sum(upreps),col="red")
lines(time2,downdgs/sum(downdgs),col="blue",lty=2)
lines(time2,updgs/sum(updgs),col="red",lty=2)



#--------- normalized total DE genes and TEs
plot(time2,totdgs/sum(totdgs),pch=19,xaxt="n",cex=0.5,col="red",ylab="% DE genes and TEs",xlab="Time (h)",ylim=c(0,0.9))#,ylim=c(0,175))#,log="y")
lines(time2,totdgs/sum(totdgs),col="red")
#lines(time2,updgs/sum(updgs),col="red")
#lines(time2,downdgs/sum(downdgs),col="red",lty=2)
#lines(time2,upreps/sum(upreps),col="red",lty=2)
#lines(time2,downreps/sum(downreps),col="blue",lty=2)
lines(time2,dreps_TEsonly/sum(dreps_TEsonly),col="darkgreen")
#lines(time2,upreps_TEsonly/sum(upreps_TEsonly),col="darkgreen")
#lines(time2,downreps_TEsonly/sum(downreps_TEsonly),col="darkgreen",lty=2)
lines(time2,v/sum(v),col="black")
lines(time2,intB/sum(intB),col="orange")
axis(side=1, at=time2, labels = TRUE)
abline(v=0,lty=2,col="grey")
abline(v=6,lty=2,col="grey")
abline(v=24,lty=2,col="grey")


#---------  absolute total # DE genes and repeats
plot(time2,totdreps,pch=19,xaxt="n",cex=0.5,col="darkgreen",ylab="number of DE genes and repeats",xlab="Time (h)")
lines(time2,totdreps,col="darkgreen")
lines(time2,totdgs,col="red")
axis(side=1, at=time2, labels = TRUE)
abline(v=0,lty=2,col="grey")
abline(v=6,lty=2,col="grey")

abline(v=24,lty=2,col="grey")

#---------  absolute total # DE UP AND DOWN genes and repeats -- dashed line = up, solid line = down
plot(time2,downreps,pch=19,xaxt="n",cex=0.5,col="darkgreen",ylab="number of DE genes and repeats",xlab="Time (h)",ylim=c(0,1500))
lines(time2,downreps,col="darkgreen",lty=2)
lines(time2,upreps,col="darkgreen")
lines(time2,updgs,col="red")
lines(time2,downdgs,col="red",lty=2)
axis(side=1, at=time2, labels = TRUE)
abline(v=0,lty=2,col="grey")
abline(v=6,lty=2,col="grey")
abline(v=20,lty=2,col="grey")


sessionInfo()
