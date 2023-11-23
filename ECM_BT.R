setwd("R:/TBT/Bioinformatics/ECM_Proteins_zero/")

zero_data_TPM <- read.delim("R:/PM/Work_In_Progress/Molecular_Profiling/RNA_Seq/ExpressionAnalysis/GeneExpression_TPM_Counts.txt",
                            sep="\t", stringsAsFactors = F)

zero_data_counts <- read.delim("R:/PM/Work_In_Progress/Molecular_Profiling/RNA_Seq/ExpressionAnalysis/GeneExpression_rawCounts.txt",
                            sep="\t", stringsAsFactors = F)

patient_metadata <- read.delim("R:/PM/Work_In_Progress/Molecular_Profiling/RNA_Seq/ExpressionAnalysis/Patients_Diagnosis.txt",
                               sep="\t", stringsAsFactors = F)

ECM_gene_list <- read.delim("ECM_genelist_maayanlab.txt",sep="\t", stringsAsFactors = F)

patient_metadata_BT <- patient_metadata[which(patient_metadata$Category == "BT"),]

#colnames(zero_data_TPM) <- gsub("\\.", "-", colnames(zero_data_TPM))

colnames(zero_data_counts) <- gsub("\\.", "-", colnames(zero_data_counts))


# index <- which(colnames(zero_data_TPM) %in% patient_metadata_BT$Patient.ID)
# data_TPM_matched.BT <- zero_data_TPM[,c(1,2,index)]


index <- which(colnames(zero_data_counts) %in% patient_metadata_BT$Patient.ID)
data_counts_matched.BT <- zero_data_counts[,c(1,2,index)]

# index.1 <- which(data_TPM_matched.BT$gene_id %in% ECM_gene_list$Symbol)
# data_TPM_matched.BT.ECM <- data_TPM_matched.BT[index.1,]

index.1 <- which(data_counts_matched.BT$gene_id %in% ECM_gene_list$Symbol)
data_coounts_matched.BT.ECM <- data_counts_matched.BT[index.1,]


write.table(data_TPM_matched.BT.ECM, "Zero_ECM_Markers_BT.txt", sep="\t", row.names = F, quote = F)

write.table(data_coounts_matched.BT.ECM, "Zero_ECM_Markers_BT_counts.txt", sep="\t", row.names = F, quote = F)

write.table(patient_metadata_BT, "BT_patient_metadata.txt", sep="\t", row.names = F, quote = F)


index.2 <- which(ECM_gene_list$Symbol %in% data_TPM_matched.BT.ECM$gene_id)
ECM_gene_list.missing <- ECM_gene_list[-index.2,]

