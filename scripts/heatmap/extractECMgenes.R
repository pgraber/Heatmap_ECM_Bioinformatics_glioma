setwd("R:/TBT/Bioinformatics/ECM_Proteins_zero/")

goi<-read.delim("ECM_genelist_maayanlab.txt",sep="\t",header=T,stringsAsFactors=F)
tpm<-read.delim("R:/PM/Work_In_Progress/Molecular_Profiling/RNA_Seq/ExpressionAnalysis/GeneExpression_TPM_Counts.txt",sep="\t",header=T,stringsAsFactors=F)
colnames(tpm)<-gsub(pattern="\\.",replacement="\\-",colnames(tpm))
dx<-read.delim("R:/PM/Work_In_Progress/Molecular_Profiling/RNA_Seq/ExpressionAnalysis/Patients_Diagnosis.txt",sep="\t",header=T,stringsAsFactors=F)

subdx<-dx[which(dx$Diagnosis=="NBL" | dx$Diagnosis=="EWS" | dx$Diagnosis=="OST" | dx$Diagnosis=="RMS FN" | dx$Diagnosis=="RMS FP"),]
index<-which(colnames(tpm) %in% subdx$Patient.ID)
subtpm<-tpm[,c(1,index)]

index<-which(subtpm$gene_id %in% goi$Symbol)

goiTpm<-subtpm[index,]
rownames(goiTpm)<-goiTpm$gene_id
goiTpm<-goiTpm[,-1]

goiTpm<-as.data.frame(t(goiTpm))
goiTpm$Patient.ID<-rownames(goiTpm)

combdata<-merge(subdx,goiTpm,by=c("Patient.ID"))

write.table(combdata,"ZERO_ECMmarkers_NBL_EWS_OST_RMS.txt",sep="\t",row.names=F,quote=F)
