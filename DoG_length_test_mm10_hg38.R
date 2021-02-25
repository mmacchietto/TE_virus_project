##########################################################################################
##########################################################################################
##
## DoG analysis HUMAN -- look at lengths of DoGs and compare infected versus mock in HUMAN
##
##########################################################################################
##########################################################################################


setwd("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/DoG/")
df <- data.frame()
x <- list.files("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/DoG/")
for (i in x){
	y <- list.files(paste("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/DoG/",i,sep=""),pattern=".bed")
	for (j in y){
		print(c(i,y))
		data <- read.delim(paste("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/DoG",i,j,sep='/'),header=FALSE)
		len <- data$V3-data$V2
		test <- cbind(i,j,len)
		df <- rbind(df,test)
	}
}

write.table(df,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/hg38/DoG/DoG.lengths.all.viruses.txt",quote=FALSE,sep="\t")
library(ggplot2)
df$len <- as.numeric(as.character(df$len))
df$ij <- paste(df$i,df$j,sep="_")

df$status <- "inf"
df[grepl("mock",df$j),]$status <- "mock"
df[grepl("ctrl",df$j),]$status <- "mock"
df$statusp <- paste(df$i,df$status,sep="_")
df$len <- df$len/1000
#ggplot(df,aes(x=ij,y=log2(len)))+ geom_boxplot()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
#ggplot(df,aes(x=ij,y=len))+ geom_boxplot()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))

pdf("hg38_DoG_lengths.pdf",width=12,height=8)
ggplot(df,aes(x=statusp,y=log2(len)))+ geom_boxplot()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
ggplot(df,aes(x=statusp,y=len))+ geom_boxplot()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))


test <- df[df$i == "DFV",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "EBOV_CD4Tcells",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "EBOV_MDM",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "EBOVarpe",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "EBV",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "HCMV",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "HCV",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "HIVact",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "HIVrest",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "HPV",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "HRSV",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "HSV-1",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "IAV",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "KSHVhuvec",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "KSHVmc116",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)


test <- df[df$i == "MeV",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "ORFV",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "RESTV",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "RVFV",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "sendai",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "VZV",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "ZIKV",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)



##########################################################################################
##########################################################################################
##
## DoG analysis MOUSE -- look at lengths of DoGs and compare infected versus mock in MOUSE
##
##########################################################################################
##########################################################################################


setwd("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/DoGs/")
df <- data.frame()
x <- list.files("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/DoGs/")
for (i in x){
	y <- list.files(paste("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/DoGs/",i,sep=""),pattern=".bed")
	for (j in y){
		print(c(i,y))
		data <- read.delim(paste("/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/DoGs",i,j,sep='/'),header=FALSE)
		len <- data$V3-data$V2
		test <- cbind(i,j,len)
		df <- rbind(df,test)
	}
}

write.table(df,"/Users/mmacchie/Desktop/TE_project/mouse_encode/TEtranscripts/viruses/mm10/DoGs/DoG.lengths.all.viruses.txt",quote=FALSE,sep="\t")
library(ggplot2)
df$len <- as.numeric(as.character(df$len))
df$ij <- paste(df$i,df$j,sep="_")

df$status <- "inf"
df[grepl("mock",df$j),]$status <- "mock"

df$statusp <- paste(df$i,df$status,sep="_")
df$len <- df$len/1000
#ggplot(df,aes(x=ij,y=log2(len)))+ geom_boxplot()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
#ggplot(df,aes(x=ij,y=len))+ geom_boxplot()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))

pdf("mm10_DoG_lengths.pdf",width=15,height=10)
ggplot(df,aes(x=statusp,y=log2(len)))+ geom_boxplot()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
ggplot(df,aes(x=statusp,y=len))+ geom_boxplot()+ theme(axis.text.x = element_text(angle = 90, hjust = 1))


test <- df[df$i == "H5N1",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "H5N8",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "HSV-1_DRG",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "HSV-1_LIM",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "IAVaecii",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "IAVam",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "IAVcc",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "IAVcd103",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "IBV",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "LCMV",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "MCMV",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "MHV-A59",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "MHV68",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "WNV",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)

test <- df[df$i == "WNV12hpi",]
test_inf <- test[test$status == "inf",]
test_mock <- test[test$status == "mock",]
t.test(test_inf$len,test_mock$len)
