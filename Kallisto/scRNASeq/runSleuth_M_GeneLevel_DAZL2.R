library('sleuth')

base_dir <- './'

samples <- c("Kal_SRR2013603",
"Kal_SRR2013604",
"Kal_SRR2013605",
"Kal_SRR2013606",
"Kal_SRR2013607",
"Kal_SRR2013608",
"Kal_SRR2013609",
"Kal_SRR2013610",
"Kal_SRR2013611",
"Kal_SRR2013612",
"Kal_SRR2013613",
"Kal_SRR2013614",
"Kal_SRR2013615",
"Kal_SRR2013616",
"Kal_SRR2013617",
"Kal_SRR2013618",
"Kal_SRR2013619",
"Kal_SRR2013620",
"Kal_SRR2013644",
"Kal_SRR2013645",
"Kal_SRR2013646",
"Kal_SRR2013647",
"Kal_SRR2013648",
"Kal_SRR2013649",
"Kal_SRR2013650",
"Kal_SRR2013651",
"Kal_SRR2013652",
"Kal_SRR2013653",
"Kal_SRR2013654",
"Kal_SRR2013655",
"Kal_SRR2013656",
"Kal_SRR2013657",
"Kal_SRR2013658",
"Kal_SRR2013659",
"Kal_SRR2013660",
"Kal_SRR2013661",
"Kal_SRR2013662",
"Kal_SRR2013663",
"Kal_SRR2013664",
"Kal_SRR2013665",
"Kal_SRR2013666",
"Kal_SRR2013667",
"Kal_SRR2013668",
"Kal_SRR2013669",
"Kal_SRR2013670",
"Kal_SRR2013671",
"Kal_SRR2013672",
"Kal_SRR2013673",
"Kal_SRR2013674")

kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))

timepoint <- c("8W",
"8W",
"8W",
"8W",
"8W",
"8W",
"8W",
"8W",
"8W",
"8W",
"8W",
"8W",
"8W",
"8W",
"8W",
"8W",
"8W",
"8W",
"17W",
"17W",
"17W",
"17W",
"17W",
"17W",
"17W",
"17W",
"17W",
"17W",
"17W",
"17W",
"17W",
"17W",
"17W",
"17W",
"17W",
"17W",
"17W",
"17W",
"17W",
"17W",
"17W",
"17W",
"17W",
"17W",
"17W",
"17W",
"17W",
"17W",
"17W")




Genes <- readRDS("Mapping1.rds")
names(Genes)[1] <- "target_id" #names(genes)[1]
names(Genes)[2] <- "ext_gene"
Genes <- Genes[2:dim(Genes)[1],]

#gtf <- ensemblGenome()
#read.gtf(gtf, filename="Mus_musculus.GRCm38.96.gtf")
#genes = gtf@ev$genes[ ,c("gene_id","gene_name")]
#genes = unique(gtf@ev$gtf[ ,c("transcript_id","gene_name")])



s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = timepoint, stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, target_mapping = Genes, aggregation_column = "ext_gene", extra_bootstrap_summary = TRUE)
so <- sleuth_fit(so, ~timepoint, 'full') #Fit full model
so <- sleuth_fit(so, ~1, 'reduced') #And reduced model
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = TRUE)
write.table(sleuth_table_gene, file = "PGC_M_8-17W_GeneLevel.csv")

saveRDS(so,file="SleuthWk8-17.rds")
#sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)

#plot_bootstrap(so, "PDPN", units = "est_counts", color_by = "condition")
#ggsave(filename=paste("BOX_PDPN_TF.pdf",sep=""),width = 13, height = 13, plot = p1)

#plot_bootstrap(so, "PIWIL2", units = "est_counts", color_by = "condition")
#ggsave(filename=paste("BOX_PIWIL2_TF.pdf",sep=""),width = 13, height = 13, plot = p1)#

#plot_bootstrap(so, "OR5AN1", units = "est_counts", color_by = "condition")
#ggsave(filename=paste("BOX_OR5AN1_TF.pdf",sep=""),width = 13, height = 13, plot = p1)

#plot_bootstrap(so, "DPPA5", units = "est_counts", color_by = "condition")
#ggsave(filename=paste("BOX_DPPA5_TF.pdf",sep=""),width = 13, height = 13, plot = p1)

#plot_bootstrap(so, "PIWIL1", units = "est_counts", color_by = "condition")
#ggsave(filename=paste("BOX_PIWIL1_TF.pdf",sep=""),width = 13, height = 13, plot = p1)


#SNM <- kallisto_table(so, use_filtered = TRUE, normalized = TRUE, include_covariates = TRUE)
#write.table(SNM, file = "PGC_M_Cl2_v_0_6_table.csv")

#SNM <- kallisto_table(so, use_filtered = TRUE, normalized = TRUE, include_covariates = TRUE)
#write.table(SNM, file = "PGC_Soma_FWk7_table.csv")


#s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = timepoint, stringsAsFactors=FALSE)
#so2 <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, extra_bootstrap_summary = TRUE)

#saveRDS(so2,file = "SleuthPr.rds")#
#so3 <- kallisto_table(so2)
#write.table(so3, file = "Test2.csv")

