##########################################################################################
##########################################################################################
##
## TEtranscripts - GSE115904 - MOUSE
##
##########################################################################################
##########################################################################################

### Load master gene table
master_gene <- read.delim("~/Desktop/TE_project/mouse_encode/counts/ENCODE-processed-totalRNA/MM_EP_STAR_FC_DFAM_finalresults/master.gene.to.te.table.txt",sep="\t",header=F)
colnames(master_gene) <- c("gene","gene_length","eff_gene_length","num_introns","cum_intron_length","num_exons","cum_exon_length","gene_strand","gene_biotype","TE_IDs","TE-gene-relations","TE-directions","num_TEs","num_TEs_exons","num_TEs_introns")
fullmaster2 <- read.delim("~/Desktop/TE_project/mouse_encode/counts/ENCODE-processed-totalRNA/MM_EP_STAR_FC_DFAM_finalresults/master.te.to.gene.table.FULL.exonannot2.txt",skip=1,header=FALSE)

### Define TPM function
tpm <- function(counts_matrix,master_gene_table){
	master_gene_table <- master_gene_table[match(rownames(counts_matrix),master_gene_table$gene),]
	cum_exon_lens <- master_gene_table$cum_exon_length/1000
	counts_lennorm <- sweep(counts_matrix,MARGIN=1,cum_exon_lens,FUN="/")
	counts_sizenorm <- sweep(counts_lennorm,MARGIN=2,colSums(counts_lennorm)/1000000,"/")
	counts_sizenorm
}

### Load data
#data <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/noSRs/TEtrans_influenza_GSE115904_novDoGsintrons.cntTable",sep="\t",header=TRUE,row.names="gene.TE")
data <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/noSRs/TEtrans_influenza_GSE115904.cntTable",sep="\t",header=TRUE,row.names="gene.TE")
colnames(data) <- sub("bams.","",colnames(data))
colnames(data) <- sub("Aligned.sortedByCoord.out.bam.C","",colnames(data))
colnames(data) <- sub("Aligned.sortedByCoord.out.bam.T","",colnames(data))


###  Virus IDs
mm10.viruses <- read.delim("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/virus.genes.mm10.txt",sep="\t",header=FALSE)


rRNAs <- subset(rownames(data),grepl("rRNA$",rownames(data)))

###  Create repeat table
test <- subset(data,!grepl("^EN",rownames(data)))
test <- test[!rownames(test) %in% mm10.viruses$V1,]
rep_table <- do.call(rbind,strsplit(as.character(rownames(test)),":"))
rep_table <- data.frame(rep_table)
rownames(rep_table) <- rownames(test)
colnames(rep_table) <- c("type","fam","class")



###  Breakdown of reads mapping to gene, repeats (no SRs), and virus
library(reshape)
library(dplyr)
inf_virus <- data[rownames(data) %in% mm10.viruses$V1,]
inf_virus <- colSums(inf_virus)/colSums(data)
inf_genes <- subset(data,grepl("^EN",rownames(data)))
inf_genes <- colSums(inf_genes)/colSums(data)
inf_reps <- subset(data,!grepl("^EN",rownames(data)))
inf_reps <- inf_reps[!rownames(inf_reps) %in% mm10.viruses$V1,]
inf_reps <- colSums(inf_reps)/colSums(data)
inf_genes <- melt(inf_genes)
inf_genes$cat <- "genes"
inf_genes$sample <- rownames(inf_genes)
inf_virus <- melt(inf_virus)
inf_virus$cat <- "virus"
inf_virus$sample <- rownames(inf_virus)
inf_reps <- melt(inf_reps)
inf_reps$cat <- "repeats"
inf_reps$sample <- rownames(inf_reps)
inf_data <- rbind(inf_genes,inf_virus,inf_reps)
inf_data$replicates <- inf_data$sample
inf_data$replicates <- gsub("_1","",inf_data$replicates)
inf_data$replicates <- gsub("_2","",inf_data$replicates)
inf_data$replicates <- gsub("_3","",inf_data$replicates)
inf_data$replicates <- gsub("_4","",inf_data$replicates)
#inf_data$replicates <- gsub("infected","inf",inf_data$replicates)
#inf_data$replicates <- gsub("control","mock",inf_data$replicates)
inf_data$replicates <- paste(inf_data$cat,inf_data$replicates,sep='_')
test <- inf_data %>% group_by(replicates) %>% summarize(avg = mean(value))
test$cat <- c('genes','genes','genes','genes','genes','genes','genes','genes','repeats','repeats','repeats','repeats','repeats','repeats','repeats','repeats','virus','virus','virus','virus','virus','virus','virus','virus')
test$status <- rep(c("inf_AECII","inf_AM","inf_CC","inf_CD103","mock_AECII","mock_AM","mock_CC","mock_CD103"),3)
library(ggplot2)
fill <- c("#5F9EA0", "#E1B378","#FF0033")
p4 <- ggplot()+geom_bar(aes(y=avg*100,x=status,fill=cat),data=test,stat="identity")+scale_fill_manual(values=fill)
p4



###  Breakdown of reads mapping to gene and virus

inf_virus <- data[rownames(data) %in% mm10.viruses$V1,]
inf_genes <- subset(data,grepl("^EN",rownames(data)))
inf_virus <- colSums(inf_virus)/(colSums(inf_virus)+colSums(inf_genes))
inf_genes <- colSums(inf_genes)/(colSums(data[rownames(data) %in% mm10.viruses$V1,])+colSums(inf_genes))
#inf_reps <- subset(data,!grepl("^EN",rownames(data)))
#inf_reps <- inf_reps[!rownames(inf_reps) %in% mm10.viruses$V1,]
#inf_reps <- colSums(inf_reps)/colSums(data)
inf_genes <- melt(inf_genes)
inf_genes$cat <- "genes"
inf_genes$sample <- rownames(inf_genes)
inf_virus <- melt(inf_virus)
inf_virus$cat <- "virus"
inf_virus$sample <- rownames(inf_virus)
#inf_reps <- melt(inf_reps)
#inf_reps$cat <- "repeats"
#inf_reps$sample <- rownames(inf_reps)
inf_data <- rbind(inf_genes,inf_virus)
inf_data$replicates <- inf_data$sample
inf_data$replicates <- gsub("_1","",inf_data$replicates)
inf_data$replicates <- gsub("_2","",inf_data$replicates)
inf_data$replicates <- gsub("_3","",inf_data$replicates)
inf_data$replicates <- gsub("_4","",inf_data$replicates)
inf_data$replicates <- paste(inf_data$cat,inf_data$replicates,sep='_')
test <- inf_data %>% group_by(replicates) %>% summarize(avg = mean(value))
test$cat <- c('genes','genes','genes','genes','genes','genes','genes','genes','virus','virus','virus','virus','virus','virus','virus','virus')
test$status <- rep(c("inf_AECII","inf_AM","inf_CC","inf_CD103","mock_AECII","mock_AM","mock_CC","mock_CD103"),2)
library(ggplot2)
# reverse last two colors so that virus is red
fill <- c("#5F9EA0","#FF0033","#E1B378")
p4 <- ggplot()+geom_bar(aes(y=avg*100,x=status,fill=cat),data=test,stat="identity")+scale_fill_manual(values=fill)
p4



###  Breakdown of reads mapping to gene and repeat types

# repeats
data2 <- data[!rownames(data) %in% mm10.viruses$V1,]
#data2 <- data[!rownames(data) %in% mm10.viruses$V1 & !rownames(data) %in% SRs,] # obsolete
ERVs <- rownames(rep_table[rep_table$class == 'LTR' | rep_table$class == 'LTR?',])
LINEs <- rownames(rep_table[rep_table$class == 'LINE',])
SINEs <- rownames(rep_table[rep_table$class == 'SINE',])
DNA <- rownames(rep_table[rep_table$class == 'DNA' | rep_table$class == 'DNA?',])
rRNA <- rownames(rep_table[rep_table$class == 'rRNA',])
inf_ERV <- data2[rownames(data2) %in% ERVs,]
inf_LINE <- data2[rownames(data2) %in% LINEs,]
inf_SINE <- data2[rownames(data2) %in% SINEs,]
inf_DNA <- data2[rownames(data2) %in% DNA,]
inf_rRNA <- data2[rownames(data2) %in% rRNA,]
inf_ERV <- colSums(inf_ERV)/colSums(data2)
inf_LINE <- colSums(inf_LINE)/colSums(data2)
inf_SINE <- colSums(inf_SINE)/colSums(data2)
inf_DNA <- colSums(inf_DNA)/colSums(data2)
inf_rRNA <- colSums(inf_rRNA)/colSums(data2)
inf_reps <- subset(data2,!grepl("^EN",rownames(data2)))
inf_other <- inf_reps[!rownames(inf_reps) %in% mm10.viruses$V1 & !rownames(inf_reps) %in% ERVs & !rownames(inf_reps) %in% LINEs & !rownames(inf_reps) %in% SINEs & !rownames(inf_reps) %in% DNA & !rownames(inf_reps) %in% rRNA,]
inf_other <- colSums(inf_other)/colSums(data2)
# genes
inf_genes <- subset(data2,grepl("^EN",rownames(data2)))
inf_genes <- colSums(inf_genes)/colSums(data2)


inf_genes <- melt(inf_genes)
inf_genes$cat <- "genes"
inf_genes$sample <- rownames(inf_genes)
inf_ERV <- melt(inf_ERV)
inf_ERV$cat <- "ERV"
inf_ERV$sample <- rownames(inf_ERV)
inf_LINE <- melt(inf_LINE)
inf_LINE$cat <- "LINE"
inf_LINE$sample <- rownames(inf_LINE)
inf_SINE <- melt(inf_SINE)
inf_SINE$cat <- "SINE"
inf_SINE$sample <- rownames(inf_SINE)
inf_DNA <- melt(inf_DNA)
inf_DNA$cat <- "DNA"
inf_DNA$sample <- rownames(inf_DNA)
inf_rRNA <- melt(inf_rRNA)
inf_rRNA$cat <- "rRNA"
inf_rRNA$sample <- rownames(inf_rRNA)
inf_other <- melt(inf_other)
inf_other$cat <- "other"
inf_other$sample <- rownames(inf_other)

inf_data <- rbind(inf_genes,inf_ERV,inf_LINE,inf_SINE,inf_DNA,inf_rRNA,inf_other)

inf_data$replicates <- inf_data$sample
inf_data$replicates <- gsub("_1","",inf_data$replicates)
inf_data$replicates <- gsub("_2","",inf_data$replicates)
inf_data$replicates <- gsub("_3","",inf_data$replicates)
inf_data$replicates <- gsub("infected","inf",inf_data$replicates)
inf_data$replicates <- gsub("control","mock",inf_data$replicates)
inf_data$replicates <- paste(inf_data$cat,inf_data$replicates,sep='-')

test <- inf_data %>% group_by(replicates) %>% summarize(avg = mean(value))
test2 <- do.call(rbind,strsplit(as.character(test$replicates),"-"))
test$cat <- test2[,1]
test$status <- test2[,2]
library(ggplot2)
library("RColorBrewer")
library(wesanderson)
#fill=terrain.colors(7)
fill=brewer.pal(n = 7, name = "RdYlBu")
#fill=brewer.pal(n = 7, name = "BrBG")
#fill=brewer.pal(n = 7, name = "Paired")
p4 <- ggplot()+geom_bar(aes(y=avg*100,x=status,fill=cat),colour="black",data=test,stat="identity")+scale_fill_manual(values=fill)
p4


library(edgeR)
# create counts matrices for comparisons
data <- data[!rownames(data) %in% mm10.viruses$V1,]
AM_naive <- data[,c(4,5,10,16,24,30,31,32)]
AECII_naive <- data[,c(7,8,9,14,17,22,22,29)]
CC_naive <- data[,c(1,2,3,11,19,21,26,28)]
CD103_naive <- data[,c(6,12,13,15,18,20,25,27)]

# create normalized count matrix
GS_cpm <- cpm(data)
GS_tpm <- tpm(data,master_gene)

### edgeR pipeline - AM vs. naive comparison
y <- DGEList(counts=AM_naive,group=c(1,1,1,1,2,2,2,2))
keep <- rowSums(cpm(y)>1) >= 4
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(AM_naive)[1])
tt <- tt$table
test <- subset(tt,!grepl("^EN",rownames(tt)))
test3 <- test[test$FDR < 0.05 & abs(test$logFC) > 1,]
#test <- test[test$FDR < 0.05 & abs(test$logFC) > 2,]
INF_AM_naive2 <- test3
upINF_AM_naive2 <- test3[test3$logFC < 0,]
dim(upINF_AM_naive2)
downINF_AM_naive2 <- test3[test3$logFC > 0,]
dim(downINF_AM_naive2)

data_reps <- high_naive[rownames(high_naive) %in% rownames(test),]
INF_AM_naive_cpm_de <- GS_cpm[rownames(GS_cpm) %in% rownames(test),]
INF_AM_naive_tpm_de <- GS_tpm[rownames(GS_tpm) %in% rownames(test),]
INF_AM_naive <- test3

upINF_AM_naive <- INF_AM_naive[INF_AM_naive$logFC < 0 & !rownames(INF_AM_naive) %in% mm10.viruses$V1,]
dim(upINF_AM_naive)
downINF_AM_naive<- INF_AM_naive[INF_AM_naive$logFC > 0 & !rownames(INF_AM_naive) %in% mm10.viruses$V1,]
dim(downINF_AM_naive)

upINF_AM_naive <- do.call(rbind,strsplit(as.character(rownames(upINF_AM_naive)),":"))
downINF_AM_naive <- do.call(rbind,strsplit(as.character(rownames(downINF_AM_naive)),":"))
table(upINF_AM_naive[,3])
table(downINF_AM_naive[,3])

test2 <- subset(tt,grepl("^EN",rownames(tt)))
test2 <- test2[test2$FDR < 0.05 & abs(test2$logFC) > 2,]
INF_AM_naive_g <- test2
upINF_AM_naive_g <- test2[test2$logFC < 0 & !rownames(test2) %in% mm10.viruses$V1,]
dim(upINF_AM_naive_g)
downINF_AM_naive_g<- test2[test2$logFC > 0 & !rownames(test2) %in% mm10.viruses$V1,]
dim(downINF_AM_naive_g)



### edgeR pipeline - AECII vs. naive comparison
y <- DGEList(counts=AECII_naive,group=c(1,1,1,1,2,2,2,2))
keep <- rowSums(cpm(y)>1) >= 4
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(AECII_naive)[1])
tt <- tt$table
test <- subset(tt,!grepl("^EN",rownames(tt)))
test3 <- test[test$FDR < 0.05 & abs(test$logFC) > 1,]
test <- test[test$FDR < 0.05 & abs(test$logFC) > 2,]
INF_AECII_naive2 <- test3
upINF_AECII_naive2 <- test3[test3$logFC < 0,]
dim(upINF_AECII_naive2)
downINF_AECII_naive2 <- test3[test3$logFC > 0,]
dim(downINF_AECII_naive2)

data_reps <- high_naive[rownames(high_naive) %in% rownames(test),]
INF_AECII_naive_cpm_de <- GS_cpm[rownames(GS_cpm) %in% rownames(test),]
INF_AECII_naive_tpm_de <- GS_tpm[rownames(GS_tpm) %in% rownames(test),]
INF_AECII_naive <- test3

upINF_AECII_naive <- INF_AECII_naive[INF_AECII_naive$logFC < 0 & !rownames(INF_AECII_naive) %in% mm10.viruses$V1,]
dim(upINF_AECII_naive)
downINF_AECII_naive<- INF_AECII_naive[INF_AECII_naive$logFC > 0 & !rownames(INF_AECII_naive) %in% mm10.viruses$V1,]
dim(downINF_AECII_naive)

upINF_AECII_naive <- do.call(rbind,strsplit(as.character(rownames(upINF_AECII_naive)),":"))
downINF_AECII_naive <- do.call(rbind,strsplit(as.character(rownames(downINF_AECII_naive)),":"))
table(upINF_AECII_naive[,3])
table(downINF_AECII_naive[,3])

test2 <- subset(tt,grepl("^EN",rownames(tt)))
test2 <- test2[test2$FDR < 0.05 & abs(test2$logFC) > 2,]
INF_AECII_naive_g <- test2
upINF_AECII_naive_g <- test2[test2$logFC < 0 & !rownames(test2) %in% mm10.viruses$V1,]
dim(upINF_AECII_naive_g)
downINF_AECII_naive_g<- test2[test2$logFC > 0 & !rownames(test2) %in% mm10.viruses$V1,]
dim(downINF_AECII_naive_g)



### edgeR pipeline - CC vs. naive comparison
y <- DGEList(counts=CC_naive,group=c(1,1,1,1,2,2,2,2))
keep <- rowSums(cpm(y)>1) >= 4
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(CC_naive)[1])
tt <- tt$table
test <- subset(tt,!grepl("^EN",rownames(tt)))
test3 <- test[test$FDR < 0.05 & abs(test$logFC) > 1,]
test <- test[test$FDR < 0.05 & abs(test$logFC) > 2,]
upINF_CC_naive2 <- test3[test3$logFC < 0,]
dim(upINF_CC_naive2)
downINF_CC_naive2 <- test3[test3$logFC > 0,]
dim(downINF_CC_naive2)

data_reps <- high_naive[rownames(high_naive) %in% rownames(test),]
INF_CC_naive_cpm_de <- GS_cpm[rownames(GS_cpm) %in% rownames(test),]
INF_CC_naive_tpm_de <- GS_tpm[rownames(GS_tpm) %in% rownames(test),]
INF_CC_naive <- test3

upINF_CC_naive <- INF_CC_naive[INF_CC_naive$logFC < 0 & !rownames(INF_CC_naive) %in% mm10.viruses$V1,]
dim(upINF_CC_naive)
downINF_CC_naive<- INF_CC_naive[INF_CC_naive$logFC > 0 & !rownames(INF_CC_naive) %in% mm10.viruses$V1,]
dim(downINF_CC_naive)

upINF_CC_naive <- do.call(rbind,strsplit(as.character(rownames(upINF_CC_naive)),":"))
downINF_CC_naive <- do.call(rbind,strsplit(as.character(rownames(downINF_CC_naive)),":"))
table(upINF_CC_naive[,3])
table(downINF_CC_naive[,3])

test2 <- subset(tt,grepl("^EN",rownames(tt)))
test2 <- test2[test2$FDR < 0.05 & abs(test2$logFC) > 2,]
INF_CC_naive_g <- test2
upINF_CC_naive_g <- test2[test2$logFC < 0 & !rownames(test2) %in% mm10.viruses$V1,]
dim(upINF_CC_naive_g)
downINF_CC_naive_g<- test2[test2$logFC > 0 & !rownames(test2) %in% mm10.viruses$V1,]
dim(downINF_CC_naive_g)




### edgeR pipeline - CD103 vs. naive comparison
y <- DGEList(counts=CD103_naive,group=c(1,1,1,1,2,2,2,2))
keep <- rowSums(cpm(y)>1) >= 4
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y)
y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)
et <- exactTest(y)
tt <- topTags(et,n=dim(CD103_naive)[1])
tt <- tt$table
test <- subset(tt,!grepl("^EN",rownames(tt)))
test3 <- test[test$FDR < 0.05 & abs(test$logFC) > 1,]
test <- test[test$FDR < 0.05 & abs(test$logFC) > 2,]
upINF_CD103_naive2 <- test3[test3$logFC < 0,]
dim(upINF_CD103_naive2)
downINF_CD103_naive2 <- test3[test3$logFC > 0,]
dim(downINF_CD103_naive2)

INF_CD103_naive_cpm_de <- GS_cpm[rownames(GS_cpm) %in% rownames(test),]
INF_CD103_naive_tpm_de <- GS_tpm[rownames(GS_tpm) %in% rownames(test),]
INF_CD103_naive <- test3

upINF_CD103_naive <- INF_CD103_naive[INF_CD103_naive$logFC < 0 & !rownames(INF_CD103_naive) %in% mm10.viruses$V1,]
dim(upINF_CD103_naive)
downINF_CD103_naive<- INF_CD103_naive[INF_CD103_naive$logFC > 0 & !rownames(INF_CD103_naive) %in% mm10.viruses$V1,]
dim(downINF_CD103_naive)

upINF_CD103_naive <- do.call(rbind,strsplit(as.character(rownames(upINF_CD103_naive)),":"))
downINF_CD103_naive <- do.call(rbind,strsplit(as.character(rownames(downINF_CD103_naive)),":"))
table(upINF_CD103_naive[,3])
table(downINF_CD103_naive[,3])

test2 <- subset(tt,grepl("^EN",rownames(tt)))
test2 <- test2[test2$FDR < 0.05 & abs(test2$logFC) > 2,]
INF_CD103_naive_g <- test2 
upINF_CD103_naive_g <- test2[test2$logFC < 0 & !rownames(test2) %in% mm10.viruses$V1,]
dim(upINF_CD103_naive_g)
downINF_CD103_naive_g<- test2[test2$logFC > 0 & !rownames(test2) %in% mm10.viruses$V1,]
dim(downINF_CD103_naive_g)

