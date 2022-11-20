library('sleuth')

base_dir <- './'

#F wk 7 vs 9
samples <- c("KalSRR1777318","KalSRR1777319","KalSRR1777315","KalSRR1777316","KalSRR1777317")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("Wk7","Wk7","Wk9","Wk9","Wk9"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, extra_bootstrap_summary = TRUE)
so <- sleuth_fit(so, ~timepoint, 'full') #Fit full model
so <- sleuth_fit(so, ~1, 'reduced') #And reduced model
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "PGC_F_Wk7-9.csv")

SNM <- kallisto_table(so, use_filtered = TRUE, normalized = TRUE, include_covariates = TRUE)
write.table(SNM, file = "PGC_F_Wk7-9_table.csv")


#M wk 5.5 vs 9
samples <- c("KalSRR1777320","KalSRR1537296","KalSRR1537297","KalSRR1777321","KalSRR1777322")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("Wk7","Wk7","Wk7","Wk9","Wk9"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, extra_bootstrap_summary = TRUE)
so <- sleuth_fit(so, ~timepoint, 'full') #Fit full model
so <- sleuth_fit(so, ~1, 'reduced') #And reduced model
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "PGC_M_Wk5.5and7-9.csv")


SNM <- kallisto_table(so, use_filtered = TRUE, normalized = TRUE, include_covariates = TRUE)
write.table(SNM, file = "PGC_M_Wk5.5and7-9_table.csv")


#M wk 5.5 vs 9
#samples <- c("KalSRR1777320","KalSRR1537296","KalSRR1537297")
#kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
#s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("Wk7","Wk7","Wk7"), stringsAsFactors=FALSE)
#so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, extra_bootstrap_summary = TRUE)
#so <- sleuth_fit(so, ~timepoint, 'full') #Fit full model
#so <- sleuth_fit(so, ~1, 'reduced') #And reduced model
#so <- sleuth_lrt(so, 'reduced', 'full')
#sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
#write.table(sleuth_table_gene, file = "PGC_M_Wk5.5-7.csv")


#F soma vs PGC
samples <- c("KalSRR1777315","KalSRR1777316","KalSRR1777317","KalSRR1777312","KalSRR1777313","KalSRR1777314")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("PGC","PGC","PGC","Soma","Soma","Soma"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, extra_bootstrap_summary = TRUE)
so <- sleuth_fit(so, ~timepoint, 'full') #Fit full model
so <- sleuth_fit(so, ~1, 'reduced') #And reduced model
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "PGC_Soma_FWk7.csv")

SNM <- kallisto_table(so, use_filtered = TRUE, normalized = TRUE, include_covariates = TRUE)
write.table(SNM, file = "PGC_Soma_FWk7_table.csv")


#M soma vs PGC
samples <- c("KalSRR1537296","KalSRR1537297","KalSRR1537298","KalSRR1537299")
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))
s2c <- data.frame(path=kal_dirs, sample=samples, timepoint = c("PGC","PGC","Soma","Soma"), stringsAsFactors=FALSE)
so <- sleuth_prep(s2c,read_bootstrap_tpm = TRUE, extra_bootstrap_summary = TRUE)
so <- sleuth_fit(so, ~timepoint, 'full') #Fit full model
so <- sleuth_fit(so, ~1, 'reduced') #And reduced model
so <- sleuth_lrt(so, 'reduced', 'full')
sleuth_table_gene <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
write.table(sleuth_table_gene, file = "PGC_Soma_MWk7.csv")

SNM <- kallisto_table(so, use_filtered = TRUE, normalized = TRUE, include_covariates = TRUE)
write.table(SNM, file = "PGC_Soma_MWk7_table.csv")

