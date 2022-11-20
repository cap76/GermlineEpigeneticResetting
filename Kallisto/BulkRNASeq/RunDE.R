library(DESeq2)

cts <- as.matrix(read.csv("CountMatrix.txt",sep="\t",row.names="Genes"))
coldata <- read.csv('sampleAn.txt', sep='\t', row.names=1)

#Do the thing
coldata2 <- data.frame(condition = coldata$condition)
rownames(coldata2) <- colnames(cts)

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata2, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)

res1<-results(dds, contrast=c("condition","PGC_7_F","PGC_9_F"))
write.csv(as.data.frame(res1), file = "PGC_F_9_vs_7.csv")

res2<-results(dds, contrast=c("condition","PGC_7_F","Som_7_F"))
write.csv(as.data.frame(res2), file = "PGC_Som_F_7.csv")

res3<-results(dds, contrast=c("condition","PGC_7_M","Som_7_M"))
write.csv(as.data.frame(res3), file = "PGC_Som_M_7.csv")


library(DESeq2)
cts <- as.matrix(read.csv("CountMatrix.txt",sep="\t",row.names="Genes"))

coldata <- read.csv('sampAn2.txt', sep='\t', row.names=1)
coldata2 <- data.frame(condition = coldata$condition2)
rownames(coldata2) <- colnames(cts)

dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata2, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)

res4<-results(dds, contrast=c("condition","PGC_5.5_M","PGC_7_M"))
write.csv(as.data.frame(res4), file = "PGC_M_5.5_vs_7_9.csv")


#res2<-results(dds, contrast=c("condition2","PGC_5.5_M","PGC_7_M"))
#write.csv(as.data.frame(res2), file = "PGC_M_5.5_vs_7_9.csv")

