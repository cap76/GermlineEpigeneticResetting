library(rstatix)

#we are missing H3K27me3_soma_wk9 in spreadsheet

#Matching list of marks in PGC and in Soma
Set1 <- c("H2AK119u1_hPGC",
"H3K27ac_hPGC",
"H3K27me3_wk9.hPGC",
"H3K27me3_wk7.hPGC",
"H3K4me3_hPGC",
"H3K9me3_hPGC",
"H3K4me1_hPGC",
"ATAC_hPGC",
"expr_wk7.hPGC",
"meth_wk7.hPGC")

Set2 <- c("H2AK119u1_Soma",
"H3K27ac_Soma",
"H3K27me3_wk9.Soma",
"H3K27me3_wk9.hPGC",
"H3K4me3_Soma",
"H3K9me3_Soma",
"H3K4me1_Soma",
"ATAC_Soma.M",
"expr.wk7.Soma",
"meth.wk7.Soma")

#Read the data table
Data <- read.table("TE_SOM_input_with_nodes.csv",sep=",", header = TRUE)

#Get unique list of nodes
nodes <- unique(Data$Nodes)

#Empty matrix for storing p value and effect size
PV <- matrix(, ncol = length(Set1), nrow = length(nodes))
ES <- matrix(, ncol = length(Set1), nrow = length(nodes))

#Loop though all nodes and conditions
for (i in 1:length(nodes)) {
  for (j in 1:length(Set1)) {
    
  #Get data for ith node and jth condition
    curNode <- nodes[i]
    D1 <- Data[which(Data$Nodes==curNode),c(Set1[j])]
    D2 <- Data[which(Data$Nodes==curNode),c(Set2[j])]
    
    # Do Wilcox test, paired (two-sideds?)
    p <- wilcox.test(D1, D2, paired = TRUE, alternative = "two.sided")
    PV[i,j] <- p$p.value #Store p val
    #p <- wilcox.test(D1, D2, paired = TRUE, alternative = "greater")

    #Get effect size 
    d <- wilcox_effsize(data.frame(x=c(D1,D2),y=as.factor(c(matrix("PGC", ncol = length(D1), nrow = length(1)),matrix("Soma", ncol = length(D1), nrow = length(1))))), x ~ y, paired = TRUE)
    ES[i,j] <- d$effsize # store eff size
  }
}


saveRDS(PV,file="PVTE.rds")
saveRDS(ES,file="ESTE.rds")
saveRDS(nodes,file="NodesTE.rds")
saveRDS(Set1,file="Set1TE.rds")
saveRDS(Set2,file="Set2TE.rds")

PV <- as.data.frame(PV)
ES <- as.data.frame(ES)

rownames(PV) <- nodes
colnames(PV) <- paste(Set1,Set2,sep="_")

rownames(ES) <- nodes
colnames(ES) <- paste(Set1,Set2,sep="_")

write.csv(PV, file=paste("./PvalsTE.csv",sep=""),quote=FALSE)
write.csv(ES, file=paste("./EffSizeTE.csv",sep=""),quote=FALSE)
