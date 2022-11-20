library('sleuth')

base_dir <- './'

samples <- c("Kal_SRR2013621",
"Kal_SRR2013630",
"Kal_SRR2013631",
"Kal_SRR2013632",
"Kal_SRR2013633",
"Kal_SRR2013622",
"Kal_SRR2013623",
"Kal_SRR2013624",
"Kal_SRR2013625",
"Kal_SRR2013626",
"Kal_SRR2013627",
"Kal_SRR2013628",
"Kal_SRR2013629",
"Kal_SRR2013634",
"Kal_SRR2013643",
"Kal_SRR2013635",
"Kal_SRR2013636",
"Kal_SRR2013637",
"Kal_SRR2013638",
"Kal_SRR2013639",
"Kal_SRR2013640",
"Kal_SRR2013641",
"Kal_SRR2013642",
"Kal_SRR2013487",
"Kal_SRR2013496",
"Kal_SRR2013497",
"Kal_SRR2013498",
"Kal_SRR2013499",
"Kal_SRR2013500",
"Kal_SRR2013501",
"Kal_SRR2013502",
"Kal_SRR2013503",
"Kal_SRR2013504",
"Kal_SRR2013505",
"Kal_SRR2013488",
"Kal_SRR2013506",
"Kal_SRR2013489",
"Kal_SRR2013490",
"Kal_SRR2013491",
"Kal_SRR2013492",
"Kal_SRR2013493",
"Kal_SRR2013494",
"Kal_SRR2013495",
"Kal_SRR2013507",
"Kal_SRR2013516",
"Kal_SRR2013517",
"Kal_SRR2013518",
"Kal_SRR2013508",
"Kal_SRR2013509",
"Kal_SRR2013510",
"Kal_SRR2013511",
"Kal_SRR2013512",
"Kal_SRR2013513",
"Kal_SRR2013514",
"Kal_SRR2013515",
"Kal_SRR2013519",
"Kal_SRR2013528",
"Kal_SRR2013529",
"Kal_SRR2013530",
"Kal_SRR2013531",
"Kal_SRR2013532",
"Kal_SRR2013533",
"Kal_SRR2013520",
"Kal_SRR2013521",
"Kal_SRR2013522",
"Kal_SRR2013523",
"Kal_SRR2013524",
"Kal_SRR2013525",
"Kal_SRR2013526",
"Kal_SRR2013527")


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
"F",
"F",
"F",
"F",
"F",
"F",
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
write.table(sleuth_table_gene, file = "PGC_MvF_10-11W_GeneLevel.csv")

#SNM <- kallisto_table(so, use_filtered = TRUE, normalized = TRUE, include_covariates = TRUE)
#write.table(SNM, file = "PGC_M_Cl2_v_0_6_table.csv")

#SNM <- kallisto_table(so, use_filtered = TRUE, normalized = TRUE, include_covariates = TRUE)
#write.table(SNM, file = "PGC_Soma_FWk7_table.csv")


#s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = timepoint, stringsAsFactors=FALSE)
#so2 <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, extra_bootstrap_summary = TRUE)

#saveRDS(so2,file = "SleuthPr.rds")#
#so3 <- kallisto_table(so2)
#write.table(so3, file = "Test2.csv")

