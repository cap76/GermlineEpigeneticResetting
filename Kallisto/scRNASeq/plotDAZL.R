
library(ggplot2)
D1 <- read.csv("/Volumes/GoogleDrive/My\ Drive/Wolfram/Kallisto/DE_gene_level/DAZL.csv",header = TRUE,sep=",",row.names = 1)


l1 <- c("F_PGC_8W","F_PGC_17W")
l2 <- c("M_PGC_7W","M_PGC_19W")

#Set1 <- grep("F_PGC_8W", colnames(D1))
#Set2 <- grep("M_PGC_7W", colnames(D1))


Set1 <- grep(l1[1], colnames(D1))
Set2 <- grep(l1[2], colnames(D1))
#Set3 <- grep(l1[3], colnames(D1))
#Set4 <- grep(l1[4], colnames(D1))
#Set5 <- grep(l1[5], colnames(D1))

Data1 <- D1[,Set1]
Data2 <- D1[,Set2]
#Data3 <- D1[,Set3]
#Data4 <- D1[,Set4]
#Data5 <- D1[,Set5]




D <- data.frame(
  
  x=unlist(c(Data1[1,],Data2[1,])), 
  y = as.factor(c(
    matrix("T1", ncol = length(Set1), nrow = length(1)),matrix("T2", ncol = length(Set2), nrow = length(1))
    ) ) 
  )
p <- ggplot(D, aes(x=y, y=x)) + geom_violin() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_bw() + ylab("TPM") + ggtitle("DAZL") + stat_summary(fun.data = "mean_cl_boot", geom = "crossbar",
                                                                                                                                                                 colour = "red", width = 0.2)
ggsave(filename=paste("BOX_DAZL_F_TF.pdf",sep=""),width = 13, height = 13, plot = p, useDingbats=FALSE)




Set1 <- grep(l2[1], colnames(D1))
Set2 <- grep(l2[2], colnames(D1))
#Set3 <- grep(l2[3], colnames(D1))
#Set4 <- grep(l2[4], colnames(D1))
#Set5 <- grep(l2[5], colnames(D1))

Data1 <- D1[,Set1]
Data2 <- D1[,Set2]
#Data3 <- D1[,Set3]
#Data4 <- D1[,Set4]
#Data5 <- D1[,Set5]


D <- data.frame(
  
  x=unlist(c(Data1[1,],Data2[1,])), 
  y = as.factor(c(
    matrix("T1", ncol = length(Set1), nrow = length(1)),matrix("T2", ncol = length(Set2), nrow = length(1))
   ) ) 
)
p <- ggplot(D, aes(x=y, y=x)) + geom_violin() + geom_jitter(shape=16, position=position_jitter(0.2)) + theme_bw() + ylab("TPM") + ggtitle("DAZL") + stat_summary(fun.data = "mean_cl_boot", geom = "crossbar",
                                                                                                                                                                 colour = "red", width = 0.2)
ggsave(filename=paste("BOX_DAZL_M_TF.pdf",sep=""),width = 13, height = 13, plot = p, useDingbats=FALSE)

