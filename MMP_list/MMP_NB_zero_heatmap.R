setwd("R:/KCA/Projects/ECM_bioprinting/")
library(readxl)
library(gplots)
library(ComplexHeatmap)
library(fastcluster)
library(circlize)
library(ggplot2)
library(ggpubr)

sheet_nam <- excel_sheets("ZERO_ECMmarkers_NBL_EWS_OST_RMS_MJ.xlsx")

##heatmap of all samples

data_ECM_all <- read_excel("ZERO_ECMmarkers_NBL_EWS_OST_RMS_MJ.xlsx", sheet = "ZERO_ECMmarkers_NBL_EWS_OST_RMS")
metadata <- data_ECM_all[,1:4]
unique(metadata$Diagnosis)

metadata_NBL <- metadata[which(metadata$Diagnosis == "NBL"),]

data_ECM_all_NBL <- data_ECM_all[which(data_ECM_all$Diagnosis == "NBL"),]

index <- grep("MMP", colnames(data_ECM_all_NBL))
data_ECM_all_NBL_MMP <- data_ECM_all_NBL[,c(1:4,index)]

data_tpm_NBL_MMP <- data_ECM_all_NBL_MMP[,c(1,5:ncol(data_ECM_all_NBL_MMP))]
rownames(data_tpm_NBL_MMP) <- data_ECM_all_NBL_MMP$Patient.ID

data_tpm_t <- as.data.frame(t(data_tpm_NBL_MMP[,-1]))
colnames(data_tpm_t) <- data_tpm_NBL_MMP$Patient.ID

mat <- as.matrix(data_tpm_t)
mat[mat==0] <- 0.001

tpm<-as.matrix(mat)

logTPM <- log(tpm)

log_scale <- t(scale(t(logTPM)))
col_fun1 <- colorRamp2(c(-1, 0, 1), c("darkorange", "white", "deeppink"))
col_fun2 <- colorRamp2(c(1, 5, 10), c("seagreen1", "royalblue", "orangered"))
fh = function(x) fastcluster::hclust(dist(x))
metadata_order <- metadata_NBL[order(colnames(log_scale)),]
metadata_order_1 <- as.data.frame(metadata_order[,2])
#col_at <- HeatmapAnnotation(df = metadata_order_1, col = list(Diagnosis = SubtypeCol))

pdf("R:/TBT/Bioinformatics/ECM_Proteins_zero/MMP_list/MMP_genes_zero_NBcohort.pdf",width = 5, height = 5)
Heatmap(as.matrix(log_scale), name = "logTPM", show_row_names = T,show_column_names = F,
        cluster_rows = fh, cluster_columns = fh,
        row_names_gp = grid::gpar(fontsize = 7))
dev.off()




data_tpm_t$avg_exp <- rowMeans(data_tpm_t)
data_tpm_t_avg_order <- data_tpm_t[order(data_tpm_t$avg_exp,decreasing = T),]
mat <- as.matrix(data_tpm_t_avg_order[,-ncol(data_tpm_t_avg_order)])
mat[mat==0] <- 0.001

tpm<-as.matrix(mat)
logTPM <- log(tpm)
log_scale <- t(scale(t(logTPM)))
col_fun1 <- colorRamp2(c(-1, 0, 1), c("darkorange", "white", "deeppink"))
col_fun2 <- colorRamp2(c(1, 5, 10), c("seagreen1", "royalblue", "orangered"))
fh = function(x) fastcluster::hclust(dist(x))
metadata_order <- metadata_NBL[order(colnames(log_scale)),]
metadata_order_1 <- as.data.frame(metadata_order[,2])
#col_at <- HeatmapAnnotation(df = metadata_order_1)

pdf("R:/TBT/Bioinformatics/ECM_Proteins_zero/MMP_list/MMP_genes_zero_NBcohort_ordered.pdf",width = 5, height = 5)
Heatmap(as.matrix(log_scale), name = "logTPM", show_row_names = T,show_column_names = F,
        cluster_rows = F, cluster_columns = fh,row_names_gp = grid::gpar(fontsize = 7)  )
dev.off()

data_tpm <- data_tpm_NBL_MMP[,-1]
rownames(data_tpm) <- data_tpm_NBL_MMP$Patient.ID
MMP_top <- c("MMP2", "MMP9", "MMP14")
for (i in 1:length(MMP_top)) {
  message(i)
  png(paste("R:/TBT/Bioinformatics/ECM_Proteins_zero/MMP_list/boxplot_MMP_vs_other_MMP.png", sep=""))
  index <- which(colnames(data_tpm) %in% MMP_top)
  data_tpm_rowmean <- data_tpm
  data_tpm_rowmean$Other_MMP_genes <-  rowMeans(data_tpm_rowmean[,-c(index)])
  data_tpm_row_gene <- as.data.frame(data_tpm_rowmean[,c(index)])
  #colnames(data_tpm_row_gene)[1] <- MMP_top[i]
  data_tpm_row <- as.data.frame(data_tpm_rowmean[,ncol(data_tpm_rowmean)])
  data_tpm_row.1 <- cbind(data_tpm_row, data_tpm_row_gene)
  x <- reshape2::melt(data_tpm_row.1)
  #data_tpm_row_gene$gene <- rep(colnames(data_tpm_row_gene)[1],nrow(data_tpm_row_gene))
  #data_tpm_row$gene <- rep("Other MMP Genes",nrow(data_tpm_row))
  colnames(x)[2] <- "TPM"
  colnames(x)[1] <- "gene"
  
  #data_tpm_row_b <- rbind(data_tpm_row, data_tpm_row_gene)
  x$logTpm <- x$TPM
  x$logTpm[x$logTpm == 0] <- 0.0001
  x$logTpm <- log(x$logTpm)
  #x$geneInt <- rep(MMP_top[i], nrow(data_tpm_row_b))
  p <- ggboxplot(x, x = "gene", y = "logTpm" ,
                 palette = "jco",
                 short.panel.labs = FALSE)+
    xlab(label = "")+
    #ggtitle(label = MMP_top[i])+
    theme(plot.title = element_text(hjust = 0.5))
  
  # Use only p.format as label. Remove method name.
  #p<- p + stat_compare_means(label = "p.format",method = "t.test",paired = T,
  #                           aes(x = factor(gene, levels = c("Other ECM genes", data_tpm_avg_Top_NB_top10$Name[i])), y = logTpm ), data = data_tpm_row_b)
  
  p<- p + stat_compare_means(method = "t.test",paired = T, label = "p.format",
                             comparisons = list(c("MMP2", "Other_MMP_genes"),
                                                c("MMP14", "Other_MMP_genes"),
                                                c("MMP9","Other_MMP_genes")))
  
  #,
  #                           aes(x = factor(gene, levels = c("Other ECM genes", data_tpm_avg_Top_NB_top10$Name[i])), y = logTpm ), data = data_tpm_row_b)
  
  print(p)
  dev.off()
}
