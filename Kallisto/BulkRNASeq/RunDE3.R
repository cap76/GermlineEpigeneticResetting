library(DESeq2)

cts <- as.matrix(read.csv("WolframDTGM.txt",sep="\t",row.names="Geneid"))
coldata <- read.csv('sampleAn3.txt', sep='\t', row.names=1)

#Do the thing
coldata2 <- data.frame(condition = coldata$condition)
rownames(coldata2) <- colnames(cts)

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata2, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)

res1<-results(dds, contrast=c("condition","hPGC_wk7_F","hPGC_wk9_F"))
write.csv(as.data.frame(res1), file = "PGC_F_9_vs_7_GM.csv")

res2<-results(dds, contrast=c("condition","hPGC_wk7_F","hSoma_wk7_F"))
write.csv(as.data.frame(res2), file = "PGC_Som_F_7_GM.csv")

res3<-results(dds, contrast=c("condition","hPGC_wk7_M","hSoma_wk7_M"))
write.csv(as.data.frame(res3), file = "PGC_Som_M_7_GM.csv")


library(DESeq2)
#cts <- as.matrix(read.csv("CountMatrix.txt",sep="\t",row.names="Genes"))
cts <- as.matrix(read.csv("WolframDTGM.txt",sep="\t",row.names="Geneid"))
coldata <- read.csv('sampleAn3.txt', sep='\t', row.names=1)


coldata2 <- data.frame(condition = coldata$condition2)
rownames(coldata2) <- colnames(cts)

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata2, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)

res4<-results(dds, contrast=c("condition","hPGC_wk5_M","hPGC_wk7_M"))
write.csv(as.data.frame(res4), file = "PGC_M_5.5_vs_7_9_GM.csv")


#res2<-results(dds, contrast=c("condition2","PGC_5.5_M","PGC_7_M"))
#write.csv(as.data.frame(res2), file = "PGC_M_5.5_vs_7_9.csv")

