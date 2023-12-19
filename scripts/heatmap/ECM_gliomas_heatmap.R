setwd("R:/TBT/Bioinformatics/ECM_Proteins_zero/")
library(readxl)
library(gplots)
library(ComplexHeatmap)
library(fastcluster)
library(circlize)

sheet_nam <- excel_sheets("ZERO_ECMmarkers BT_patient_metadata.xlsx")

##heatmap of all samples

data_ECM_Gliomas <- read_excel("ZERO_ECMmarkers BT_patient_metadata.xlsx", sheet = "Gliomas")
data_ECM_Gliomas <- data_ECM_Gliomas[-c(nrow(data_ECM_Gliomas),nrow(data_ECM_Gliomas)-1 ),]
metadata <- data_ECM_Gliomas[,1:4]
unique(metadata$Diagnosis)


data_tpm <- data_ECM_Gliomas[,c(1,5:ncol(data_ECM_Gliomas))]
rownames(data_tpm) <- data_tpm$Patient.ID

data_ECM_avg <- read_excel("ZERO_ECMmarkers BT_patient_metadata.xlsx", sheet = "Gliomas Rank")
data_tpm_avg <- data_tpm
data_tpm_avg <- data_tpm_avg[,-1]
data_tpm_avg <- as.data.frame(t(data_tpm_avg))
colnames(data_tpm_avg) <- data_tpm$Patient.ID
data_tpm_avg$avg_exp <- rowMeans(data_tpm_avg)


data_tpm_avg_order <- data_tpm_avg[order(data_tpm_avg$avg_exp,decreasing = T),]
x <- nrow(data_tpm_avg_order)-29
data_tpm_avg_order_top <- data_tpm_avg_order[c(1:30,x:nrow(data_tpm_avg_order)),]

index <- which(colnames(data_tpm) %in% rownames(data_tpm_avg_order_top))
data_tpm_top <- data_tpm[,c(1,index)]

data_tpm_t <- t(data_tpm_top[,-1])
colnames(data_tpm_t) <- data_tpm$Patient.ID

mat <- data_tpm_t

mat[mat==0] <- 0.001

tpm<-as.matrix(mat)

logTPM <- log(tpm)

###complex heatmap


##scale the matrix before using Complexheatmap 
log_scale <- t(scale(t(logTPM)))


col_fun1 <- colorRamp2(c(-1, 0, 1), c("darkorange", "white", "deeppink"))
col_fun2 <- colorRamp2(c(1, 5, 10), c("seagreen1", "royalblue", "orangered"))

fh = function(x) fastcluster::hclust(dist(x))

#"Glioma other" "HGG"          "DMG"  

# SubtypeCol <- c("EWS" = "goldenrod", "NBL" = "darkgreen", "OST" ="lightskyblue1",
#                 "RMS FN"= "purple", "RMS FP"=  "red")
SubtypeCol <- c("Glioma other" = "goldenrod", "HGG" = "darkgreen", "DMG" ="lightskyblue1")

m <- match(colnames(log_scale), metadata$Patient.ID)
metadata_order <- metadata[m,]


metadata_order_1 <- as.data.frame(metadata_order[,2])
#col_at <- HeatmapAnnotation(df = KPT_aoc_log_1, col = list(AOC = col_fun2,Model  = colGrp ))
col_at <- HeatmapAnnotation(df = metadata_order_1, col = list(Diagnosis = SubtypeCol))





###### 
# png("heatmap_all_samples_diag_high_v1.png", height=800, width=900,res=100)
# Heatmap(as.matrix(log_scale), name = "logTPM", show_row_names = F,show_column_names = F, 
#         top_annotation = col_at  )
# dev.off()
# 
# png("heatmap_all_samples_diag_high_v2.png", height=800, width=900,res=100)
# Heatmap(as.matrix(log_scale), name = "logTPM", show_row_names = T,show_column_names = F, 
#         top_annotation = col_at  )
# dev.off()

#####
png("heatmap_gliomas_top_30_cluster.png", height=800, width=900,res=100)
Heatmap(as.matrix(log_scale), name = "logTPM", show_row_names = T,show_column_names = F,
        top_annotation = col_at , row_names_gp = gpar(fontsize = 5), cluster_rows = fh, 
        cluster_columns = fh )
dev.off()



####
data_tpm_avg_order_top_1 <- data_tpm_avg_order_top[,-ncol(data_tpm_avg_order_top)]

mat <- data_tpm_avg_order_top_1

mat[mat==0] <- 0.001

tpm<-as.matrix(mat)

logTPM <- log(tpm)

##scale the matrix before using Complexheatmap 
log_scale <- t(scale(t(logTPM)))



#"Glioma other" "HGG"          "DMG"  
SubtypeCol <- c("Glioma other" = "goldenrod", "HGG" = "darkgreen", "DMG" ="lightskyblue1")

m <- match(colnames(log_scale), metadata$Patient.ID)
metadata_order <- metadata[m,]


metadata_order_1 <- as.data.frame(metadata_order[,2])
#col_at <- HeatmapAnnotation(df = KPT_aoc_log_1, col = list(AOC = col_fun2,Model  = colGrp ))
col_at <- HeatmapAnnotation(df = metadata_order_1, col = list(Diagnosis = SubtypeCol))

png("heatmap_gliomas_top_30_ordered.png", height=800, width=900,res=100)
Heatmap(as.matrix(log_scale), name = "logTPM", show_row_names = T,show_column_names = F,
        top_annotation = col_at , row_names_gp = gpar(fontsize = 5), cluster_rows = F, cluster_columns = fh )
dev.off()

#######filtering low counts
data_counts <- read.delim("Zero_ECM_Markers_BT_counts.txt",sep="\t", stringsAsFactors = F)
library(edgeR)
keep<-rowSums(cpm(data_counts[,3:ncol(data_counts)])>1)>=(0.5*ncol(data_counts[,3:ncol(data_counts)])+1)
counts<-data_counts[keep,]
counts_left <- data_counts[!keep,]

data_ECM_Gliomas <- read_excel("Gliomas/ZERO_ECMmarkers BT_patient_metadata.xlsx", sheet = "Gliomas")
data_ECM_Gliomas <- data_ECM_Gliomas[-c(nrow(data_ECM_Gliomas),nrow(data_ECM_Gliomas)-1 ),]

index <- which(colnames(data_ECM_Gliomas) %in% counts$gene_id)
data_ECM_Gliomas.1 <- data_ECM_Gliomas[,c(1:4,index)]

data_tpm_out <- data_ECM_Gliomas[,-index]
writeLines(colnames(data_tpm_out)[-c(1:4)],"Gliomas/filtered_geneList.txt")

metadata <- data_ECM_Gliomas[,1:4]
unique(metadata$Diagnosis)


data_tpm <- data_ECM_Gliomas[,c(1,5:ncol(data_ECM_Gliomas))]
rownames(data_tpm) <- data_tpm$Patient.ID

#data_ECM_avg <- read_excel("ZERO_ECMmarkers BT_patient_metadata.xlsx", sheet = "Gliomas Rank")
data_tpm_avg <- data_tpm
data_tpm_avg <- data_tpm_avg[,-1]
data_tpm_avg <- as.data.frame(t(data_tpm_avg))
colnames(data_tpm_avg) <- data_tpm$Patient.ID
data_tpm_avg$avg_exp <- rowMeans(data_tpm_avg)


data_tpm_avg_order <- data_tpm_avg[order(data_tpm_avg$avg_exp,decreasing = T),]
x <- nrow(data_tpm_avg_order)-29
data_tpm_avg_order_top <- data_tpm_avg_order[c(1:30,x:nrow(data_tpm_avg_order)),]

index <- which(colnames(data_tpm) %in% rownames(data_tpm_avg_order_top))
data_tpm_top <- data_tpm[,c(1,index)]

data_tpm_t <- t(data_tpm_top[,-1])
colnames(data_tpm_t) <- data_tpm$Patient.ID

mat <- data_tpm_t

mat[mat==0] <- 0.001

tpm<-as.matrix(mat)

logTPM <- log(tpm)

###complex heatmap


##scale the matrix before using Complexheatmap 
log_scale <- t(scale(t(logTPM)))


col_fun1 <- colorRamp2(c(-1, 0, 1), c("darkorange", "white", "deeppink"))
col_fun2 <- colorRamp2(c(1, 5, 10), c("seagreen1", "royalblue", "orangered"))

fh = function(x) fastcluster::hclust(dist(x))

#"Glioma other" "HGG"          "DMG"  

# SubtypeCol <- c("EWS" = "goldenrod", "NBL" = "darkgreen", "OST" ="lightskyblue1",
#                 "RMS FN"= "purple", "RMS FP"=  "red")
SubtypeCol <- c("Glioma other" = "goldenrod", "HGG" = "darkgreen", "DMG" ="lightskyblue1")

m <- match(colnames(log_scale), metadata$Patient.ID)
metadata_order <- metadata[m,]


metadata_order_1 <- as.data.frame(metadata_order[,2])
#col_at <- HeatmapAnnotation(df = KPT_aoc_log_1, col = list(AOC = col_fun2,Model  = colGrp ))
col_at <- HeatmapAnnotation(df = metadata_order_1, col = list(Diagnosis = SubtypeCol))





###### 
# png("heatmap_all_samples_diag_high_v1.png", height=800, width=900,res=100)
# Heatmap(as.matrix(log_scale), name = "logTPM", show_row_names = F,show_column_names = F, 
#         top_annotation = col_at  )
# dev.off()
# 
# png("heatmap_all_samples_diag_high_v2.png", height=800, width=900,res=100)
# Heatmap(as.matrix(log_scale), name = "logTPM", show_row_names = T,show_column_names = F, 
#         top_annotation = col_at  )
# dev.off()

#####
png("heatmap_gliomas_top_30_cluster_filtered.png", height=800, width=900,res=100)
Heatmap(as.matrix(log_scale), name = "logTPM", show_row_names = T,show_column_names = F,
        top_annotation = col_at , row_names_gp = gpar(fontsize = 5), cluster_rows = fh, 
        cluster_columns = fh )
dev.off()


data_tpm_avg_order_top_1 <- data_tpm_avg_order_top[,-ncol(data_tpm_avg_order_top)]

mat <- data_tpm_avg_order_top_1

mat[mat==0] <- 0.001

tpm<-as.matrix(mat)

logTPM <- log(tpm)

##scale the matrix before using Complexheatmap 
log_scale <- t(scale(t(logTPM)))



#"Glioma other" "HGG"          "DMG"  
SubtypeCol <- c("Glioma other" = "goldenrod", "HGG" = "darkgreen", "DMG" ="lightskyblue1")

m <- match(colnames(log_scale), metadata$Patient.ID)
metadata_order <- metadata[m,]


metadata_order_1 <- as.data.frame(metadata_order[,2])
#col_at <- HeatmapAnnotation(df = KPT_aoc_log_1, col = list(AOC = col_fun2,Model  = colGrp ))
col_at <- HeatmapAnnotation(df = metadata_order_1, col = list(Diagnosis = SubtypeCol))

png("heatmap_gliomas_top_30_ordered_fil.png", height=800, width=900,res=100)
Heatmap(as.matrix(log_scale), name = "logTPM", show_row_names = T,show_column_names = F,
        top_annotation = col_at , row_names_gp = gpar(fontsize = 5), cluster_rows = F, cluster_columns = fh )
dev.off()

write.table(data_tpm_avg_order, "Gliomas/ZERO_ECM_glioma_fil_genelist.txt", sep="\t", quote = F)

########box plot###########
data_ECM_Gliomas <- read_excel("Gliomas/ZERO_ECMmarkers BT_patient_metadata.xlsx", sheet = "Gliomas")
data_ECM_Gliomas <- data_ECM_Gliomas[-c(nrow(data_ECM_Gliomas),nrow(data_ECM_Gliomas)-1 ),]

metadata <- data_ECM_Gliomas[,1:4]

data_tpm <- as.data.frame(data_ECM_Gliomas[,c(1,5:ncol(data_ECM_Gliomas))])
rownames(data_tpm) <- data_tpm$Patient.ID
library(ggpubr)

data_tpm_avg <- data_tpm
data_tpm_avg <- data_tpm_avg[,-1]
data_tpm_avg <- as.data.frame(t(data_tpm_avg))
colnames(data_tpm_avg) <- data_tpm$Patient.ID
data_tpm_avg$avg_exp <- rowMeans(data_tpm_avg)


data_tpm_avg_order <- data_tpm_avg[order(data_tpm_avg$avg_exp,decreasing = T),]

data_tpm_avg_order_top_30 <- data_tpm_avg_order[c(1:30),]


dir.create("Gliomas/Boxplot")

data_tpm_avg_Top_NB_top10 <- data_tpm_avg_order_top_30[1:30,]
data_tpm_avg_Top_NB_top10 <- cbind(rownames(data_tpm_avg_Top_NB_top10), data_tpm_avg_Top_NB_top10)
colnames(data_tpm_avg_Top_NB_top10)[1] <- "Name" 
i <- 1
for (i in 1:nrow(data_tpm_avg_Top_NB_top10)) {
  message(i)
  png(paste("Gliomas/Boxplot/Gliomas_top10_",data_tpm_avg_Top_NB_top10$Name[i],"_all.png", sep=""))
  index <- which(colnames(data_tpm) %in% data_tpm_avg_Top_NB_top10$Name[i])
  data_tpm_rowmean <- data_tpm
  data_tpm_rowmean$Other_ECM_genes <-  rowMeans(data_tpm_rowmean[,-c(1,index)])
  data_tpm_row_gene <- as.data.frame(data_tpm_rowmean[,c(index)])
  colnames(data_tpm_row_gene)[1] <- data_tpm_avg_Top_NB_top10$Name[i]
  data_tpm_row <- as.data.frame(data_tpm_rowmean[,ncol(data_tpm_rowmean)])
  data_tpm_row_gene$gene <- rep(colnames(data_tpm_row_gene)[1],nrow(data_tpm_row_gene))
  data_tpm_row$gene <- rep("Other ECM Genes",nrow(data_tpm_row))
  colnames(data_tpm_row_gene)[1] <- "TPM"
  colnames(data_tpm_row)[1] <- "TPM"
  
  data_tpm_row_b <- rbind(data_tpm_row, data_tpm_row_gene)
  data_tpm_row_b$logTpm <- data_tpm_row_b$TPM
  data_tpm_row_b$logTpm[data_tpm_row_b$logTpm == 0] <- 0.0001
  data_tpm_row_b$logTpm <- log(data_tpm_row_b$logTpm)
  data_tpm_row_b$geneInt <- rep(data_tpm_avg_Top_NB_top10$Name[i], nrow(data_tpm_row_b))
  p <- ggboxplot(data_tpm_row_b, x = "gene", y = "logTpm" ,
                 palette = "jco",
                 short.panel.labs = FALSE)+
    xlab(label = "")+
    ggtitle(label = data_tpm_avg_Top_NB_top10$Name[i])+
    theme(plot.title = element_text(hjust = 0.5))
  
  # Use only p.format as label. Remove method name.
  #p<- p + stat_compare_means(label = "p.format",method = "t.test",paired = T,
  #                           aes(x = factor(gene, levels = c("Other ECM genes", data_tpm_avg_Top_NB_top10$Name[i])), y = logTpm ), data = data_tpm_row_b)
  
  p<- p + stat_compare_means(method = "t.test",paired = T, label = "p.format")
  #,
  #                           aes(x = factor(gene, levels = c("Other ECM genes", data_tpm_avg_Top_NB_top10$Name[i])), y = logTpm ), data = data_tpm_row_b)
  
  print(p)
  dev.off()
}



list.files("Gliomas/", full.names = F)

