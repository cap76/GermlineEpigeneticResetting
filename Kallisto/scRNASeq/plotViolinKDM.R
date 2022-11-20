

D1 <- read.csv("/Volumes/GoogleDrive/My\ Drive/Wolfram/Kallisto/KDMA1.csv",header = TRUE,sep=",",row.names = 1)

Set1 <- grep("F_PGC_8W", colnames(D1))
Set2 <- grep("M_PGC_7W", colnames(D1))

DataF <- D1[,Set1]
DataM <- D1[,Set2]



D <- data.frame(x=unlist(c(DataF[1,],x=DataM[1,])), y = as.factor(c(matrix("F", ncol = length(Set1), nrow = length(1)),matrix("M", ncol = length(Set2), nrow = length(1))) ) )
p <- ggplot(D, aes(x=y, y=x)) + geom_violin() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_bw() + ylab("TPM") + ggtitle("KDM6A")
ggsave(filename=paste("BOX_KDM6A_TF.pdf",sep=""),width = 13, height = 13, plot = p)

D <- data.frame(x=unlist(c(DataF[2,],x=DataM[2,])), y = as.factor(c(matrix("F", ncol = length(Set1), nrow = length(1)),matrix("M", ncol = length(Set2), nrow = length(1))) ) )
p <- ggplot(D, aes(x=y, y=x)) + geom_violin() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_bw() + ylab("TPM") + ggtitle("KDM6B")
ggsave(filename=paste("BOX_KDM6B_TF.pdf",sep=""),width = 13, height = 13, plot = p)


D <- data.frame(x=unlist(c(DataF[1,],x=DataM[1,])), y = as.factor(c(matrix("F", ncol = length(Set1), nrow = length(1)),matrix("M", ncol = length(Set2), nrow = length(1))) ) )
p <- ggplot(D, aes(x=y, y=x)) + geom_violin() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_bw() + ylab("TPM") + ggtitle("DPPA5")
ggsave(filename=paste("BOX_DPPA5_TF.pdf",sep=""),width = 13, height = 13, plot = p)


D <- data.frame(x=unlist(c(DataF[2,],x=DataM[2,])), y = as.factor(c(matrix("F", ncol = length(Set1), nrow = length(1)),matrix("M", ncol = length(Set2), nrow = length(1))) ) )
p <- ggplot(D, aes(x=y, y=x)) + geom_violin() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_bw() + ylab("TPM") + ggtitle("PDPN")
ggsave(filename=paste("BOX_PDPN_TF.pdf",sep=""),width = 13, height = 13, plot = p)


D <- data.frame(x=unlist(c(DataF[3,],x=DataM[3,])), y = as.factor(c(matrix("F", ncol = length(Set1), nrow = length(1)),matrix("M", ncol = length(Set2), nrow = length(1))) ) )
p <- ggplot(D, aes(x=y, y=x)) + geom_violin() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_bw() + ylab("TPM") + ggtitle("PIWIL1")
ggsave(filename=paste("BOX_piWIL1_TF.pdf",sep=""),width = 13, height = 13, plot = p)

D <- data.frame(x=unlist(c(DataF[4,],x=DataM[4,])), y = as.factor(c(matrix("F", ncol = length(Set1), nrow = length(1)),matrix("M", ncol = length(Set2), nrow = length(1))) ) )
p <- ggplot(D, aes(x=y, y=x)) + geom_violin() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_bw() + ylab("TPM") + ggtitle("PIWIL2")
ggsave(filename=paste("BOX_PIWIL2_TF.pdf",sep=""),width = 13, height = 13, plot = p)


D <- data.frame(x=unlist(c(DataF[5,],x=DataM[5,])), y = as.factor(c(matrix("F", ncol = length(Set1), nrow = length(1)),matrix("M", ncol = length(Set2), nrow = length(1))) ) )
p <- ggplot(D, aes(x=y, y=x)) + geom_violin() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_bw() + ylab("TPM") + ggtitle("OR5AN1")
ggsave(filename=paste("BOX_OR5AN1_TF.pdf",sep=""),width = 13, height = 13, plot = p)
