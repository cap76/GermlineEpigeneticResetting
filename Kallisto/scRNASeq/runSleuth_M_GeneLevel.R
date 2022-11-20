library('sleuth')

base_dir <- './'

samples <- c("Kal_SRR2013591",
"Kal_SRR2013592",
"Kal_SRR2013593",
"Kal_SRR2013602",
"Kal_SRR2013594",
"Kal_SRR2013595",
"Kal_SRR2013596",
"Kal_SRR2013597",
"Kal_SRR2013598",
"Kal_SRR2013599",
"Kal_SRR2013600",
"Kal_SRR2013601",
"Kal_SRR2013442",
"Kal_SRR2013443",
"Kal_SRR2013444",
"Kal_SRR2013445",
"Kal_SRR2013446",
"Kal_SRR2013447")

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
write.table(sleuth_table_gene, file = "PGC_MvF_4W_GeneLevel.csv")

#SNM <- kallisto_table(so, use_filtered = TRUE, normalized = TRUE, include_covariates = TRUE)
#write.table(SNM, file = "PGC_M_Cl2_v_0_6_table.csv")

#SNM <- kallisto_table(so, use_filtered = TRUE, normalized = TRUE, include_covariates = TRUE)
#write.table(SNM, file = "PGC_Soma_FWk7_table.csv")


#s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = timepoint, stringsAsFactors=FALSE)
#so2 <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, extra_bootstrap_summary = TRUE)

#saveRDS(so2,file = "SleuthPr.rds")#
#so3 <- kallisto_table(so2)
#write.table(so3, file = "Test2.csv")

