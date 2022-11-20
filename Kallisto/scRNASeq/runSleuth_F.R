library('sleuth')

base_dir <- './'

samples <- c("Kal_SRR2013637",
"Kal_SRR2013640",
"Kal_SRR2013591",
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
"Kal_SRR2013603",
"Kal_SRR2013614",
"Kal_SRR2013615",
"Kal_SRR2013616",
"Kal_SRR2013617",
"Kal_SRR2013618",
"Kal_SRR2013619",
"Kal_SRR2013620",
"Kal_SRR2013621",
"Kal_SRR2013622",
"Kal_SRR2013624",
"Kal_SRR2013634",
"Kal_SRR2013643",
"Kal_SRR2013635",
"Kal_SRR2013636",
"Kal_SRR2013638",
"Kal_SRR2013639",
"Kal_SRR2013641",
"Kal_SRR2013647",
"Kal_SRR2013604",
"Kal_SRR2013605",
"Kal_SRR2013606",
"Kal_SRR2013607",
"Kal_SRR2013608",
"Kal_SRR2013609",
"Kal_SRR2013631",
"Kal_SRR2013632",
"Kal_SRR2013633",
"Kal_SRR2013625",
"Kal_SRR2013626",
"Kal_SRR2013627",
"Kal_SRR2013642",
"Kal_SRR2013659",
"Kal_SRR2013612",
"Kal_SRR2013613",
"Kal_SRR2013610",
"Kal_SRR2013611")

kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))

timepoint <- c("Early",
"Early",
"Early",
"Early",
"Early",
"Early",
"Early",
"Early",
"Early",
"Early",
"Early",
"Early",
"Early",
"Early",
"Early",
"Early",
"Early",
"Early",
"Early",
"Early",
"Early",
"Early",
"Late",
"Late",
"Late",
"Late",
"Late",
"Late",
"Late",
"Late",
"Late",
"Late",
"Late",
"Late",
"Late",
"Late",
"Late",
"Late",
"Late",
"Late",
"Late",
"Late",
"Late",
"Late",
"Late",
"Late",
"Late",
"Late",
"Late",
"Late",
"Late")


s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = timepoint, stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, extra_bootstrap_summary = TRUE)
so <- sleuth_fit(so, ~timepoint, 'full') #Fit full model
so <- sleuth_fit(so, ~1, 'reduced') #And reduced model
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = TRUE)
write.table(sleuth_table_gene, file = "PGC_F_Cl2_v_0_6.csv")

SNM <- kallisto_table(so, use_filtered = TRUE, normalized = TRUE, include_covariates = TRUE)
write.table(SNM, file = "PGC_F_Cl2_v_0_6_table.csv")



#s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = timepoint, stringsAsFactors=FALSE)
#so2 <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, extra_bootstrap_summary = TRUE)

#saveRDS(so2,file = "SleuthPr.rds")#
#so3 <- kallisto_table(so2)
#write.table(so3, file = "Test2.csv")

