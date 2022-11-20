library('sleuth')

base_dir <- './'

samples <- c("Kal_SRR2013752",
"Kal_SRR2013753",
"Kal_SRR2013754",
"Kal_SRR2013755",
"Kal_SRR2013756",
"Kal_SRR2013757",
"Kal_SRR2013758",
"Kal_SRR2013759",
"Kal_SRR2013760",
"Kal_SRR2013534",
"Kal_SRR2013543",
"Kal_SRR2013544",
"Kal_SRR2013545",
"Kal_SRR2013546",
"Kal_SRR2013547",
"Kal_SRR2013548",
"Kal_SRR2013549",
"Kal_SRR2013550",
"Kal_SRR2013551",
"Kal_SRR2013552",
"Kal_SRR2013535",
"Kal_SRR2013553",
"Kal_SRR2013554",
"Kal_SRR2013555",
"Kal_SRR2013556",
"Kal_SRR2013557",
"Kal_SRR2013558",
"Kal_SRR2013559",
"Kal_SRR2013560",
"Kal_SRR2013561",
"Kal_SRR2013562",
"Kal_SRR2013536",
"Kal_SRR2013563",
"Kal_SRR2013564",
"Kal_SRR2013537",
"Kal_SRR2013538",
"Kal_SRR2013539",
"Kal_SRR2013540",
"Kal_SRR2013541",
"Kal_SRR2013542",
"Kal_SRR2013565",
"Kal_SRR2013574",
"Kal_SRR2013575",
"Kal_SRR2013576",
"Kal_SRR2013577",
"Kal_SRR2013578",
"Kal_SRR2013579",
"Kal_SRR2013580",
"Kal_SRR2013581",
"Kal_SRR2013582",
"Kal_SRR2013583",
"Kal_SRR2013566",
"Kal_SRR2013584",
"Kal_SRR2013585",
"Kal_SRR2013586",
"Kal_SRR2013587",
"Kal_SRR2013588",
"Kal_SRR2013589",
"Kal_SRR2013590",
"Kal_SRR2013567",
"Kal_SRR2013568",
"Kal_SRR2013569",
"Kal_SRR2013570",
"Kal_SRR2013571",
"Kal_SRR2013572",
"Kal_SRR2013573")

kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))

timepoint <- c("F",
"F",
"F",
"F",
"F",
"F",
"F",
"F",
"F",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M",
"M")



Genes <- readRDS("Mapping1.rds")
names(Genes)[1] <- "target_id" #names(genes)[1]
names(Genes)[2] <- "ext_gene"
Genes <- Genes[2:dim(Genes)[1],]

#gtf <- ensemblGenome()
#read.gtf(gtf, filename="Mus_musculus.GRCm38.96.gtf")
#genes = gtf@ev$genes[ ,c("gene_id","gene_name")]
#genes = unique(gtf@ev$gtf[ ,c("transcript_id","gene_name")])



s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = timepoint, stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes,  aggregation_column = "ext_gene", extra_bootstrap_summary = TRUE)
so <- sleuth_fit(so, ~timepoint, 'full') #Fit full model
so <- sleuth_fit(so, ~1, 'reduced') #And reduced model
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = TRUE)
write.table(sleuth_table_gene, file = "PGC_MvF_17-19W_GeneLevel.csv")

#SNM <- kallisto_table(so, use_filtered = TRUE, normalized = TRUE, include_covariates = TRUE)
#write.table(SNM, file = "PGC_M_Cl2_v_0_6_table.csv")

#SNM <- kallisto_table(so, use_filtered = TRUE, normalized = TRUE, include_covariates = TRUE)
#write.table(SNM, file = "PGC_Soma_FWk7_table.csv")


#s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = timepoint, stringsAsFactors=FALSE)
#so2 <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, extra_bootstrap_summary = TRUE)

#saveRDS(so2,file = "SleuthPr.rds")#
#so3 <- kallisto_table(so2)
#write.table(so3, file = "Test2.csv")

