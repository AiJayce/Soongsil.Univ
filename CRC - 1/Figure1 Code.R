set.seed(2860)

options(verbose = FALSE) # Seurat 함수들이 실행될 때 로그 메시지를 표시하지 않습니다.
options(tidyverse.quiet = TRUE) # tidyverse 패키지의 로그 메시지가 출력되지 않습니다.
options(warn=-1)
options(future.rng.onMisuse = "ignore")
options(future.globals.maxSize = 1e11)

library(viridis)
library(stringr)
library(ggplot2)
library(patchwork)
require(cluster)
library(colorRamp2)
library(gridExtra )
library(scater)
library(Seurat)
library(tidyverse)
library(cowplot)
library(dplyr)
library(reshape2)
library(S4Vectors)
library(pheatmap)
library(RColorBrewer)
library(data.table)
library(ComplexHeatmap)
require(circlize)
library(scRepertoire )


total_CRC <- readRDS("/home/Data_Drive_8TB_3/Coloncancer/new_metadata_rds/H1-31_PC25_Seurat_v4.2.2_Final_2.rds") #stromal cell meta


Immune_CRC <- readRDS("/home/Data_Drive_8TB_3/Coloncancer/new_metadata_rds/H1-31_CRC_immune.rds")



T_CRC_3 <- readRDS("/home/Data_Drive_8TB_3/Coloncancer/new_metadata_rds/Re_Step5_Final_annotation.scale.Tcell_nofiltered.rds")

T_CRC_3_filtered <- readRDS("/home/Data_Drive_8TB_3/Coloncancer/new_metadata_rds/Re_Step5_Final_annotation.scale.Tcells.rds")
source("~/Code/Seurat_CRC/CRC_Metadata_Calling.R")

Tonly_CRC_filtered <- subset(T_CRC_3_filtered, Annotation_T_tmp != "NK")

CD4T_CRC <- subset(T_CRC_3_filtered, Annotation_T_Fig1B  %in% c("CD4 Tnaive","CD4 Tcm","CD4 Tfh","CD4 Th17","CD4 Treg")   )

CD8T_CRC <- subset(T_CRC_3_filtered, Annotation_T_Fig1B  %in% c("CD8 Tex","CD8 Trm","CD8 Temra","CD8 MAIT","CD8 Tem") )




# Epi_CRC <- readRDS("/home/Data_Drive_8TB_3/Coloncancer/new_metadata_rds/H1-31_Epi_Scoring_and_Scaling_epi_with.batchcorrection_PC19_res0.3_exclude_nUMI.untilUMAP_2.rds")
Epi_Cyling_CRC <- readRDS("/home/Data_Drive_8TB_3/Coloncancer/new_metadata_rds/H1-31_Epi_CyclingCell.rds")



nPatients_list <- paste0("H",seq(1:31)[-10])





as.character(T_CRC_3$seurat_clusters)-> T_CRC_3$Annotation_total_lv4

as.data.frame(T_CRC_3$Annotation_total_lv4) -> tmp_meta
colnames(tmp_meta) <- "Celltype"

# T_CRC$Annotation_T_lv1 

as.data.frame(T_CRC_3_filtered$Annotation_T_Fig1B) -> tmp_meta_T
colnames(tmp_meta_T) <- "Celltype"

tmp_meta[rownames(tmp_meta_T), "Celltype"] <- tmp_meta_T[, 1]

T_CRC_3$Annotation_T_Fig1B_upload <- tmp_meta

T_CRC_3$Annotation_T_Fig1B_upload[T_CRC_3$Annotation_T_Fig1B_upload %in% c("19" ,"20","24")] <- "Unused"



as.character(Tonly_CRC_filtered$seurat_clusters)-> Tonly_CRC_filtered$Annotation_total_lv4

as.data.frame(Tonly_CRC_filtered$Annotation_total_lv4) -> tmp_meta
colnames(tmp_meta) <- "Celltype"

# T_CRC$Annotation_T_lv1 

as.data.frame(T_CRC_3_filtered$Annotation_T_Fig1B) -> tmp_meta_T
colnames(tmp_meta_T) <- "Celltype"

tmp_meta[rownames(tmp_meta_T), "Celltype"] <- tmp_meta_T[, 1]

Tonly_CRC_filtered$Annotation_T_Fig1B_upload <- tmp_meta





Idents(T_CRC_3) <- T_CRC_3$Annotation_T_Fig1B_upload
clusters <- T_CRC_3$Annotation_T_Fig1B_upload
cluster_names <- rep("Unknown", length(clusters))
# cluster_names <- clusters

cluster_names[clusters %in% c("CD4 Tnaive","CD4 Tcm","CD4 Tfh","CD4 Th17","CD4 Treg")] <- "CD4_T"
cluster_names[clusters %in% c("CD8 Tex","CD8 Trm","CD8 Temra","CD8 MAIT","CD8 Tem")] <- "CD8_T"

cluster_names[clusters %in% c("NK")] <- "NK"
cluster_names[clusters %in% c("Unknown(MT+)")] <- "Unknown"

cluster_names[clusters %in% c("Unused")] <- "Unused"

T_CRC_3 <- AddMetaData(T_CRC_3, metadata = cluster_names, col.name = "Annotation_T_tmp_rough")




Idents(Tonly_CRC_filtered) <- Tonly_CRC_filtered$Annotation_T_Fig1B_upload
clusters <- Tonly_CRC_filtered$Annotation_T_Fig1B_upload
cluster_names <- rep("Unknown", length(clusters))
# cluster_names <- clusters

cluster_names[clusters %in% c("CD4 Tnaive","CD4 Tcm","CD4 Tfh","CD4 Th17","CD4 Treg")] <- "CD4_T"
cluster_names[clusters %in% c("CD8 Tex","CD8 Trm","CD8 Temra","CD8 MAIT","CD8 Tem")] <- "CD8_T"

cluster_names[clusters %in% c("NK")] <- "NK"
cluster_names[clusters %in% c("Unknown(MT+)")] <- "Unknown"

cluster_names[clusters %in% c("Unknown")] <- "Unknown"

Tonly_CRC_filtered <- AddMetaData(Tonly_CRC_filtered, metadata = cluster_names, col.name = "Annotation_T_tmp_rough")



as.character(total_CRC$annotation)-> total_CRC$Annotation_total


as.data.frame(total_CRC$Annotation_total) -> tmp_meta
colnames(tmp_meta) <- "Celltype"

# T_CRC$Annotation_T_lv1 

as.data.frame(T_CRC_3$Annotation_T_Fig1B_upload) -> tmp_meta_T
colnames(tmp_meta_T) <- "Celltype"

tmp_meta[rownames(tmp_meta_T), "Celltype"] <- tmp_meta_T[, 1]

total_CRC$Annotation_total_Tsub <- tmp_meta

# DimPlot(total_CRC, reduction = "umap.unintegrated", group.by = "Annotation_total_lv3_T", label=T)



as.character(total_CRC$annotation)-> total_CRC$Annotation_total_lv4

as.data.frame(total_CRC$Annotation_total_lv4) -> tmp_meta
colnames(tmp_meta) <- "Celltype"

# T_CRC$Annotation_T_lv1 

as.data.frame(T_CRC_3$Annotation_T_tmp_rough) -> tmp_meta_T
colnames(tmp_meta_T) <- "Celltype"

tmp_meta[rownames(tmp_meta_T), "Celltype"] <- tmp_meta_T[, 1]

total_CRC$Annotation_total_Trough <- tmp_meta

# DimPlot(total_CRC, reduction = "umap.unintegrated", group.by = "Annotation_total_lv3_Tsub", label=T)






as.character(Immune_CRC$annotation) -> Immune_CRC$Annotation_total_lv4

as.data.frame(Immune_CRC$Annotation_total_lv4) -> tmp_meta
colnames(tmp_meta) <- "Celltype"

# T_CRC$Annotation_T_lv1 

as.data.frame(T_CRC_3$Annotation_T_Fig1B_upload) -> tmp_meta_T
colnames(tmp_meta_T) <- "Celltype"

tmp_meta[rownames(tmp_meta_T), "Celltype"] <- tmp_meta_T[, 1]

Immune_CRC$Annotation_total_Tsub <- tmp_meta

# DimPlot(Immune_CRC, reduction = "umap.unintegrated", group.by = "Annotation_total_lv3_T", label=T)



as.character(Immune_CRC$annotation) -> Immune_CRC$Annotation_total_lv4

as.data.frame(Immune_CRC$Annotation_total_lv4) -> tmp_meta
colnames(tmp_meta) <- "Celltype"

# T_CRC$Annotation_T_lv1 

as.data.frame(T_CRC_3$Annotation_T_tmp_rough) -> tmp_meta_T
colnames(tmp_meta_T) <- "Celltype"

tmp_meta[rownames(tmp_meta_T), "Celltype"] <- as.character(tmp_meta_T[, 1])

Immune_CRC$Annotation_total_Trough <- tmp_meta

# DimPlot(Immune_CRC, reduction = "umap.unintegrated", group.by = "Annotation_total_lv3_Tsub", label=T)


# Cell type Proportion ##########################
total_CRC -> total_CRC_subT
total_CRC -> total_CRC_subTrough

# total_CRC -> total_CRC_sub_T_sub

total_CRC$annotation_tmp <- total_CRC$annotation
total_CRC_subT$annotation_tmp <- total_CRC$Annotation_total_Tsub
total_CRC_subTrough$annotation_tmp <- total_CRC$Annotation_total_Trough

# total_CRC_sub_T_sub$annotation_tmp <- total_CRC_sub_T_sub$Annotation_total_lv3_Tsub


# tmp_tcrc <-subset(T_CRC , Annotation_T_lv1 == "CD4T"|Annotation_T_lv1 == "CD8T"|Annotation_T_lv1 == "Prolifer_T")
Tonly_CRC_filtered -> Tonly_CRC_filtered_2
# Tonly_CRC_filtered
# "epithelial" -> Epi_Cyling_CRC$annotation_tmp

Immune_CRC -> Immune_CRC_2
Immune_CRC -> Immune_CRC_3


Immune_CRC$annotation_tmp <- Immune_CRC$annotation
Immune_CRC_2$annotation_tmp <- Immune_CRC_2$Annotation_total_Tsub    
Immune_CRC_3$annotation_tmp <- Immune_CRC_3$Annotation_total_Trough

Tonly_CRC_filtered$annotation_tmp <- Tonly_CRC_filtered$Annotation_T_Fig1B_upload
Tonly_CRC_filtered_2$annotation_tmp <- Tonly_CRC_filtered$Annotation_T_tmp_rough

CD4T_CRC$annotation_tmp <- CD4T_CRC$Annotation_T_Fig1B
CD8T_CRC$annotation_tmp <- CD8T_CRC$Annotation_T_Fig1B



# CD8TNK_CRC$annotation.CD8TNK -> CD8TNK_CRC$annotation_tmp
# Myeloid_CRC$annotation.mye -> Myeloid_CRC$annotation_tmp
# B_CRC$annotation.B -> B_CRC$annotation_tmp
# Stromal_CRC$annotation.stromal -> Stromal_CRC$annotation_tmp





# CellType_CRClist <- list(total_CRC ,Epi_Cyling_CRC, CD4T_CRC  ,CD8TNK_CRC, Myeloid_CRC , B_CRC  ,Stromal_CRC)
CellType_CRClist <- list(total_CRC , total_CRC_subT , total_CRC_subTrough, 
                         Immune_CRC ,Immune_CRC_2 ,Immune_CRC_3 , 
                         Tonly_CRC_filtered ,Tonly_CRC_filtered_2  , CD4T_CRC  ,CD8T_CRC )
CellType_CRCnamelist <- c("total" , "total_subT", "total_sub_Trough" , 
                          "Immune", "Immune_T_sub","Immune_Trough",
                          "T_subtype", "T_rough", "CD4T"  ,"CD8T")


for (i in 1:length(CellType_CRClist) ){
  CellType_CRC <- CellType_CRClist[[i]]
  
  # 데이터프레임 만들기
  data_df <- data.frame(
    Patient = CellType_CRC$patients,     # 환자 정보
    CellType = CellType_CRC$annotation_tmp  # 세포 타입 정보
  )
  
  mat_temp <- as.data.frame(table(data_df) )
  
  mat <- as.data.frame(matrix(mat_temp$Freq, ncol = length(unique(mat_temp$Patient)), byrow = TRUE ,dimnames = list(colnames( table(data_df)  ),rownames(table(data_df) ) )  )  )
  rm(data_df)
  rm(mat_temp)
  
  
  mat <- mat[, nPatients_list]
  head(mat)
  
  nCell_patients <- as.data.frame(t(as.data.frame(as.matrix(table(CellType_CRC$patients)))[1]) )
  rownames(nCell_patients) <- "nCells"
  
  
  mat_percent <- mat
  
  for (col in 1:ncol(mat)) {
    col_name <- colnames(mat)[col]
    mat_percent[[col_name]] <- (mat[[col_name]] / nCell_patients[[col_name]]) * 100
  }
  
  assign( sprintf("matpercent_%s", CellType_CRCnamelist[i]), mat_percent )
  
}

# CellType_pct_list <- list(pct_total = matpercent_total , pct_CD4T = matpercent_CD4T  ,pct_CD8TNK =  matpercent_CD8TNK, pct_Myeloid =  matpercent_Myeloid , pct_B =  matpercent_B  ,pct_Stromal =  matpercent_Stromal)
CellType_pct_list <- list(pct_total = matpercent_total , pct_total_subT = matpercent_total_subT ,  pct_total_subTrough = matpercent_total_sub_Trough ,
                          pct_Immune = matpercent_Immune , pct_Immune_Tsub = matpercent_Immune_T_sub, pct_Immune_Trough = matpercent_Immune_Trough ,
                          pct_Tsub = matpercent_T_subtype  , pct_Trough = matpercent_T_rough ,pct_CD4T = matpercent_CD4T  ,pct_CD8T =  matpercent_CD8T
)

Epi_CRC <- NULL
Immune_CRC <- NULL

# CD4T_CRC <- NULL
# CD8T_CRC <- NULL
rm(Epi_CRC )
rm(Immune_CRC )
# rm(CD4T_CRC )
# rm(CD8T_CRC )


rm(total_CRC_subT )
rm(total_CRC_subTrough )
rm(Immune_CRC_2 )
rm(Immune_CRC_3 )
rm(Tonly_CRC_filtered_2 )
rm(CellType_CRClist )
rm(CellType_CRCnamelist )
rm(CellType_CRC )


rm(tmp_meta )
rm(tmp_meta_T )

rm(mat )
rm(mat_percent )
rm(matpercent_CD4T )

rm(matpercent_CD8T )
rm(matpercent_Immune )
rm(matpercent_Immune_T_sub )
rm(matpercent_Immune_Trough )
rm(matpercent_T_rough )
rm(matpercent_T_subtype )
rm(matpercent_total )
rm(matpercent_total_subT )
rm(matpercent_total_sub_Trough )



# Fig1A - Total Cell Profiling #################################


Idents(total_CRC) <- total_CRC$seurat_clusters

# sum(table(combined.integrated@meta.data$SCT_snn_res.0.2))
clusters <- total_CRC$seurat_clusters
cluster_names <- rep("Unknown", length(clusters))

cluster_names[clusters %in% c(1,4)] <- "CD4T"
cluster_names[clusters %in% c(2)] <- "CD8T"
cluster_names[clusters %in% c(7)] <- "CD8T/NK"
cluster_names[clusters %in% c(0, 8)] <- "Epithelial"
cluster_names[clusters %in% c(9)] <- "Fibroblast"   # stromal
cluster_names[clusters %in% c(10)] <- "Endothelial"
cluster_names[clusters %in% c(11)] <- "Mast"
cluster_names[clusters %in% c(6)] <- "Plasma"
cluster_names[clusters %in% c(5)] <- "B"
cluster_names[clusters %in% c(3,12,13)] <- "Myeloid"

total_CRC <- AddMetaData(total_CRC, metadata = cluster_names, col.name = "Celltype_Fig1")


total_CRC$Annotation_Fig1 <- total_CRC$Annotation_total_Trough 
total_CRC$Annotation_Fig1[total_CRC$Annotation_Fig1 == "Cycling cells"] <- "Epithelial"
total_CRC$Annotation_Fig1[total_CRC$Annotation_Fig1 == "Unused"] <- "Unknown"
total_CRC$Annotation_Fig1[total_CRC$Annotation_Fig1 == "Unknown(MT+)"] <- "Unknown"
total_CRC$Annotation_Fig1[total_CRC$Annotation_Fig1 == "CD4_T"] <- "CD4T"
total_CRC$Annotation_Fig1[total_CRC$Annotation_Fig1 == "CD8_T"] <- "CD8T"


total_CRC$Annotation_Fig1A <- total_CRC$Annotation_total_Trough 
total_CRC$Annotation_Fig1A[total_CRC$Annotation_Fig1A == "Cycling cells"] <- "Epithelial"
total_CRC$Annotation_Fig1A[total_CRC$Annotation_Fig1A == "Unused"] <- "T/NK"
total_CRC$Annotation_Fig1A[total_CRC$Annotation_Fig1A == "Unknown"] <- "T/NK"
total_CRC$Annotation_Fig1A[total_CRC$Annotation_Fig1A == "Unknown(MT+)"] <- "T/NK"
total_CRC$Annotation_Fig1A[total_CRC$Annotation_Fig1A == "CD4_T"] <- "T/NK"
total_CRC$Annotation_Fig1A[total_CRC$Annotation_Fig1A == "CD8_T"] <- "T/NK"
total_CRC$Annotation_Fig1A[total_CRC$Annotation_Fig1A == "NK"] <- "T/NK"
total_CRC$Annotation_Fig1A[total_CRC$Annotation_Fig1A == "CD19B"] <- "B"
total_CRC$Annotation_Fig1A[total_CRC$Annotation_Fig1A == "plasmaB"] <- "Plasma"

# saveRDS(total_CRC, "/home/Data_Drive_8TB_3/Coloncancer/Temp_Rds/CRC_Total_Seurat.rds")

Idents(total_CRC) <- total_CRC$Annotation_Fig1A
total_CRC@active.ident <- factor(total_CRC@active.ident , levels = c("Epithelial",  "T/NK"  , "Myeloid",   "B",    "Plasma"  ,   "Fibroblasts",     "Endothelial",    "Mast"         ))
DimPlot(total_CRC ,label=T , label.size = 7 , )   + NoAxes()


# `Color Code
# Epi - "#FFB6C1", 
# T/NK - "#87CEEB"  ( cd4t "#4169E1"  cd8t "#6A5ACD"    Nk"#2E8B57" ) 
# 
# Myeloid - "#FDFD96",
# B - "#FFB347",  Plasma - "#FFCC99"  ,  
# Fibroblasts - "#77DD77",    Endothelial -  "#BFFF00" , 
# Mast - "#FF6961"`


DimPlot(total_CRC ,label=T , label.size = 7,
        cols = c("Epithelial" = "#FFB6C1",  "T/NK" = "#87CEEB", "Myeloid"= "#BDFCC9",   "B"= "#FFB347",    "Plasma"= "#FFCC99"  ,   "Fibroblasts"= "#77DD77",     "Endothelial" = "#BFFF00" ,    "Mast" = "#FF6961"
                 ))  + theme(legend.text = element_text(size= 14) )+ NoAxes()  + xlim(-14, 18) + ylim(-18, 14)

## Fig1A-1 #######
# cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure1A1_total_UMAP.pdf" , width = 6, height = 4)
# DimPlot(total_CRC ,label=T , label.size = 6) + ggtitle(NULL) + xlim(-14, 18) + ylim(-18, 14)  + NoAxes()
# dev.off()
# 
# cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure1A1_total_UMAP_labelF.pdf" , width = 6, height = 4)
# DimPlot(total_CRC ,label=F , label.size = 6) + ggtitle(NULL) + xlim(-14, 18) + ylim(-18, 14)  + NoAxes()
# dev.off()

cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure1A1_total_UMAP_labelT_margin55.pdf" , width = 6, height = 4)
DimPlot(total_CRC ,label=T , label.size = 6, 
        cols = c("Epithelial" = "#FFB6C1",  "T/NK" = "#87CEEB", "Myeloid"= "#BDFCC9",   "B"= "#FFB347",    "Plasma"= "#FFCC99"  ,   "Fibroblasts"= "#77DD77",     "Endothelial" = "#BFFF00" ,    "Mast" = "#FF6961"))+ 
  ggtitle(NULL) + xlim(-14, 18) + ylim(-18, 14)  + NoAxes() +theme(legend.text = element_text(margin = margin(t = 5, b = 5)) )
dev.off()




# Idents(total_CRC) <- total_CRC$seurat_clusters
# 
# total_CRC@active.ident <- factor(total_CRC@active.ident, 
#                                  levels=c("0","8","1","4","2","7","3","12","13","5","6","9","10","11") )
# 
# 
# 
# 
# DotPlot(total_CRC, features =  list("Epithelial" = c("EPCAM"), "CD4T/CD8T/NK" = c("CD3D","CD3E", "CD4", "CD8A", "CD8B", "GZMK" ,"NKG7"),  "Myeloid" = c("CD68","S100A8", "FLT3"),
#                                     "B" = c("CD79A", "MS4A1", "IGHG1","IGHA2"), "Fibroblast" = c("IGFBP7","DCN") , "Endothelial" = c("PECAM1","CLDN5"),
#                                      "Mast" = c("KIT")    ) , dot.scale = 5  )  +   scale_colour_viridis(option="viridis")  +  theme_bw() + theme(axis.text.x = element_text(angle = 90))   + xlab('') +  ylab('')


# cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/figure1_DotPlot.pdf" , width = 12, height = 4)
# 
# DotPlot(total_CRC, features =  list("Epithelial" = c("EPCAM"), "CD4T/CD8T/NK" = c("CD3D","CD3E", "CD4", "CD8A", "CD8B", "GZMK" ,"NKG7"),  "Myeloid" = c("CD68","S100A8", "FLT3"),
#                                     "B" = c("CD79A", "MS4A1", "IGHG1","IGHA2"), "Fibroblast" = c("IGFBP7","DCN") , "Endothelial" = c("PECAM1","CLDN5"),
#                                     "Mast" = c("KIT")    ) , dot.scale = 8  )  +   scale_colour_viridis(option="viridis")  +  theme_bw() + 
#   theme(axis.text.x = element_text(angle = 90, family="Arial", face='bold' , size = 10 ), axis.text.y = element_text(angle = 90 , size = 8) )   + xlab('') +  ylab('')
# 
# 
# dev.off()




Idents(total_CRC) <- total_CRC$Annotation_Fig1A
table(total_CRC$Annotation_Fig1A)
total_CRC@active.ident <- factor(total_CRC@active.ident, 
                                 levels=c("Mast" ,"Endothelial" , "Fibroblasts",  "Plasma",  "B", "Myeloid", "T/NK",   "Epithelial") )


## Fig1A-2 ##############
cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure1A2_total_DotPlot.pdf" , width = 8, height = 4)

DotPlot(total_CRC, features =  list("Epithelial" = c("EPCAM"), "T/NK" = c("CD3D","CD3E","NKG7"),  "Myeloid" = c("CD68","S100A8", "FLT3"),
                                    "B" = c("CD79A", "MS4A1"),"Plasma" = c("IGHG1","IGHA2"), "Fibroblast" = c("IGFBP7","DCN") , "Endothelial" = c("PECAM1","CLDN5"),
                                    "Mast" = c("KIT")    ) , dot.scale = 8 ,  )  +   scale_colour_viridis(option="viridis")  +  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, family="Arial", face='bold' , size = 10 ), 
        axis.text.y = element_text(angle = 0 , size = 10, hjust = 0 ,  face='bold'), 
        strip.text.x = element_blank()
  )   + xlab('') +  ylab('')

dev.off()





## Fig1A-2 Sup Cluster UMAP##############

total_CRC@active.ident <- factor(total_CRC@active.ident, 
                                 levels= c(seq(0,13) )) 


cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure1A1Sup_total_Cluster_UMAP.pdf" , width = 6, height = 4)
DimPlot(total_CRC ,label=T , label.size = 6 )+ 
  ggtitle(NULL) + xlim(-14, 18) + ylim(-18, 14)  + NoAxes() +theme(legend.text = element_text(margin = margin(t = 0, b = 0)) )
dev.off()


## Fig1A-2 Sup ##############
DimPlot(total_CRC, group.by = "seurat_clusters", label=T)

Idents(total_CRC) <- total_CRC$seurat_clusters

total_CRC@active.ident <- factor(total_CRC@active.ident, 
                                 levels= rev(c(seq(0,13) )) )

total_CRC_marker <- FindAllMarkers(total_CRC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
total_CRC_marker <- total_CRC_marker[total_CRC_marker$p_val_adj <= 0.05, ]

total_CRC_marker %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 2) %>%
  ungroup() -> total_CRC_marker_topgene

total_CRC@active.ident <- factor(total_CRC@active.ident, 
                                 levels= rev(c(seq(0,13) )) )

### Fig1A-2 Sup clusterDEG ##############
cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure1A2Sup_total_DotPlot_clusterDEG_increasing.pdf" , width = 10, height = 4)

DotPlot(total_CRC, features =  rev(unique(total_CRC_marker_topgene$gene))     , dot.scale = 8 ,  )  +   scale_colour_viridis(option="viridis")  +  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, family="Arial", face='bold' , size = 10 ), 
        axis.text.y = element_text(angle = 0 , size = 10, hjust = 0 ,  face='bold'), 
        strip.text.x = element_blank()
  )   + xlab('') +  ylab('')

dev.off()


Idents(total_CRC) <- total_CRC$seurat_clusters
total_CRC@active.ident <- factor(total_CRC@active.ident, 
                                 levels= rev(c(0,8 , 1,2,4,7 , 3,12,13 , 5, 6, 9 , 10 , 11 )) )

total_CRC_marker <- FindAllMarkers(total_CRC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
total_CRC_marker <- total_CRC_marker[total_CRC_marker$p_val_adj <= 0.05, ]

total_CRC_marker %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 2) %>%
  ungroup() -> total_CRC_marker_topgene


total_CRC@active.ident <- factor(total_CRC@active.ident, 
                                 levels= rev(c(0,8 , 1,2,4,7 , 3,12,13 , 5, 6, 9 , 10 , 11 )) )


cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure1A2Sup_total_DotPlot_clusterDEG_Stepsorted.pdf" , width = 10, height = 4)

DotPlot(total_CRC, features = rev(unique(total_CRC_marker_topgene$gene) )    , dot.scale = 8 ,  )  +   scale_colour_viridis(option="viridis")  +  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, family="Arial", face='bold' , size = 10 ), 
        axis.text.y = element_text(angle = 0 , size = 10, hjust = 0 ,  face='bold'), 
        strip.text.x = element_blank()
  )   + xlab('') +  ylab('')

dev.off()

 





Idents(total_CRC) <- total_CRC$seurat_clusters
total_CRC@active.ident <- factor(total_CRC@active.ident, 
                                 levels= rev(c(0,8 , 1,2,4,7 , 3,12,13 , 5, 6, 9 , 10 , 11 )) )


cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure1A2Sup_total_DotPlot_clustermarker.pdf" , width = 10, height = 4)    

DotPlot(total_CRC, features =  list("Epithelial" = c("EPCAM"), "T/NK" = c("CD3D","CD3E","NKG7"),  "Myeloid" = c("CD68","S100A8", "FLT3"),
                                    "B" = c("CD79A", "MS4A1"),"Plasma" = c("IGHG1","IGHA2"), "Fibroblast" = c("IGFBP7","DCN") , "Endothelial" = c("PECAM1","CLDN5"),
                                    "Mast" = c("KIT")    ) , dot.scale = 8 ,  )  +   scale_colour_viridis(option="viridis")  +  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, family="Arial", face='bold' , size = 10 ), 
        axis.text.y = element_text(angle = 0 , size = 10, hjust = 0 ,  face='bold'), 
        strip.text.x = element_blank()
  )   + xlab('') +  ylab('')

dev.off()




DotPlot(total_CRC, features = rev(unique(total_CRC_marker_topgene$gene) )    , dot.scale = 8  )  +   scale_colour_viridis(option="viridis")  +  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, family="Arial", face='bold' , size = 10 ), 
        axis.text.y = element_text(angle = 0 , size = 10, hjust = 0 ,  face='bold'), 
        strip.text.x = element_blank()
  )   + xlab('') +  ylab('')



library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(circlize)
library(ComplexHeatmap)
library(grid)
library(viridis)
library(Polychrome)



Idents(total_CRC) <- total_CRC$seurat_clusters
total_CRC@active.ident <- factor(total_CRC@active.ident, 
                                 levels= rev(c(0,8 , 1,2,4,7 , 3,12,13 , 5, 6, 9 , 10 , 11 )) )



# DotPlot 생성
p <- DotPlot(total_CRC, features =  c("EPCAM","CD3D","CD3E","NKG7","CD68","S100A8", "FLT3","CD79A", "MS4A1","IGHG1","IGHA2","IGFBP7","DCN","PECAM1","CLDN5","KIT")
             , dot.scale = 8 ,  )  +   scale_colour_viridis(option="viridis")  +  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, family="Arial", face='bold' , size = 10 ), 
        axis.text.y = element_text(angle = 0 , size = 10, hjust = 0 ,  face='bold'), 
        strip.text.x = element_blank()
  )   + xlab('') +  ylab('')


df <- p$data
head(df)

# exp_mat 전치
exp_mat <- df %>%
  dplyr::select(-pct.exp, -avg.exp) %>%
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>%
  as.data.frame()

row.names(exp_mat) <- exp_mat$features.plot
exp_mat <- exp_mat[,-1] %>% as.matrix()
exp_mat <- as.data.frame(t(exp_mat)) # 전치
exp_mat <- exp_mat[nrow(exp_mat):1, ]

head(exp_mat)

# percent_mat 전치
percent_mat <- df %>%
  dplyr::select(-avg.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = id, values_from = pct.exp) %>%
  as.data.frame()

row.names(percent_mat) <- percent_mat$features.plot
percent_mat <- percent_mat[,-1] %>% as.matrix()
percent_mat <- as.data.frame(t(percent_mat)) # 전치
percent_mat <- percent_mat[nrow(percent_mat):1, ]

head(percent_mat)

breaks <- seq(min(df$avg.exp.scaled), max(df$avg.exp.scaled), length.out = 30)
colors <- viridis(length(breaks))
col_fun <- circlize::colorRamp2(breaks, colors)

# 그룹 정의
row_groups <- factor(c("Epithelial", "Epithelial", 
                       "T/NK", "T/NK", "T/NK", "T/NK",
                       "Myeloid", "Myeloid", "Myeloid",
                       "B",
                       "Plasma",
                       "Fibroblasts",
                       "Endothelial",
                       "Mast"), 
                     levels = c("Epithelial", "T/NK", "Myeloid", "B", "Plasma", "Fibroblasts", "Endothelial", "Mast"))

# 행 annotation 생성
row_anno <- rowAnnotation(Group = row_groups, 
                          col = list(Group = c("Epithelial" = "#FFB6C1", 
                                               "T/NK" = "#87CEEB", "Myeloid"= "#BDFCC9", 
                                               "B"= "#FFB347",    "Plasma"= "#FFCC99"  ,  
                                               "Fibroblasts"= "#77DD77",     "Endothelial" = "#BFFF00" ,  
                                               "Mast" = "#FF6961")),
                          annotation_label = "CellType", show_legend=F,
                          show_annotation_name = T)


column_groups <- factor(c(rep("A", length.out = 1), 
                          rep("B", length.out = 3), 
                          rep("C", length.out = 3), 
                          rep("D", length.out = 4), 
                          rep("E", length.out = 2), 
                          rep("F", length.out = 2), 
                          
                          rep("G", length.out = 1)  ))



# Modified cell_fun with grid lines behind circles
cell_fun <- function(j, i, x, y, w, h, fill){
  # Draw horizontal lines
  grid.lines(x = unit(c(x - w/2, x + w/2), "npc"),
             y = unit(c(y, y), "npc"), 
             gp = gpar(col = "#EBEBEB", lty = "solid"))
  
  # Draw vertical lines
  # Only draw once per row
  grid.lines(x = unit(c(x, x), "npc"),
             y = unit(c(y - h/2, y + h/2), "npc"),
             gp = gpar(col = "#EBEBEB", lty = "solid"))
  
  # Draw the circle for each cell
  grid.circle(x = x, y = y, 
              r = percent_mat[i, j]/100 * max(unit.c(w, h)) * 0.2,
              gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))
}
ht <- Heatmap(as.matrix(exp_mat), 
              cluster_rows = FALSE, cluster_columns = FALSE,
              heatmap_legend_param = list(title = "Average Expression"),
              show_row_names = TRUE, row_names_side = "left", row_names_rot = 0, row_names_gp = gpar(fontsize = 10, fontfamily = "Arial"),
              show_column_names = TRUE, column_names_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Arial"),
              row_title = NULL,
              col = col_fun, 
              cell_fun = cell_fun,
              rect_gp = gpar(type = "none"),
              border = "black",
              show_heatmap_legend = FALSE,
              left_annotation = row_anno, 
              column_split = column_groups, 
              column_title = NULL)

# Corrected legends
annotation_legend <- Legend(labels = c("Epithelial", "T/NK", "Myeloid", "B", "Plasma", "Fibroblasts", "Endothelial", "Mast"),
                            title = "CellType",
                            legend_gp = gpar(fill =c( "#FFB6C1", 
                                                      "#87CEEB",
                                                      "#BDFCC9", 
                                                      "#FFB347",  "#FFCC99"  ,  
                                                      "#77DD77",  "#BFFF00" ,   
                                                      "#FF6961" 
                                                      
                            )))

percent_expressed_legend <- Legend(at = c(0, 25, 50, 75), 
                                   title = "Percent Expressed",
                                   legend_gp = gpar(col = "black"),
                                   type = "points", 
                                   size = unit(c(0.5, 1.25, 2.5, 5), "mm"),
                                   background = NULL) # Corrected

average_expression_legend <- Legend(at = c(-1, 0, 1, 2, 3),
                                    col_fun = col_fun,
                                    title = "Average Expression")

combined_legend <- packLegend(annotation_legend, average_expression_legend, percent_expressed_legend,  direction = "vertical")

# Draw heatmap with corrected legends
draw(ht, annotation_legend_side = "right", annotation_legend_list = combined_legend, padding = unit(c(2, 2, 2, 2), "mm"))




### Fig1A-2 Sup Marker Annotated Dotplot ##############
cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure1A2Sup_total_DotPlot_clustermarker_Annotation.pdf" , width = 10, height = 4)
draw(ht, annotation_legend_side = "right", annotation_legend_list = combined_legend, padding = unit(c(3, 3, 3, 3), "mm"))

dev.off()





# Fig1Sup - CEBPB analysis #####################

## CEBPB cor Epithelial Proportion #####################
cor.test (colSums(CellType_pct_list$pct_total[c("Epithelial","Cycling cells"),]) , CEBPB_MeanExp_Vector, method = "pearson" )
cor.test (colSums(CellType_pct_list$pct_total[c("Epithelial"),]) , CEBPB_MeanExp_Vector, method = "pearson" )

cor.test (colSums(CellType_pct_list$pct_total[c("Epithelial","Cycling cells"),]) , CEBPB_MeanExp_Vector, method = "spearman" )
cor.test (colSums(CellType_pct_list$pct_total[c("Epithelial"),]) , CEBPB_MeanExp_Vector, method = "spearman" )


colSums(CellType_pct_list$pct_total[c("Epithelial","Cycling cells"),]) -> epi_pct


ScatterFrame_0 <- data.frame( CEBPB_Epi_Exp = CEBPB_MeanExp_Vector ,
                              epi_pct = epi_pct
)

# theme_minimal()
# theme_void()
# theme_classic()


### Fig1_sup #######
cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure1SUP_CEBPB_Epi_Cor.pdf" , width = 4, height = 4)

ggplot( ScatterFrame_0, aes(x = CEBPB_MeanExp_Vector)) +
  geom_smooth(aes(y = epi_pct, color = "CEBPB_Epi_Exp"), method = "lm", se = TRUE , linetype="longdash") +
  
  geom_point(aes(y = epi_pct, color = "CEBPB_Epi_Exp"), size = 1 , alpha = 0.5) +
  
  # geom_point(aes(y = CEBPB_Hpct, color = "CEBPB_Hpct"), size = 1 , alpha = 0.5) +
  # geom_smooth(aes(y = CEBPB_Hpct, color = "CEBPB_Hpct"), method = "lm", se = FALSE , linetype="longdash") +
  
  labs(y = "Epithelial Proportion",
       x = "Normalized CEBPB(Epi) Expr", color = "") +
  theme_classic() + 
  ggtitle("")+
  scale_color_manual( values = c("CEBPB_Epi_Exp" = "#F8766D"),
                      labels = c("CEBPB_Epi_Exp" = sprintf("Cor : %s \n P : %s", round(as.numeric(cor.test( CEBPB_MeanExp_Vector  , epi_pct)$estimate),2) ,round(cor.test( CEBPB_MeanExp_Vector  , epi_pct)$p.value,2) )
                                 
                      ) 
  ) + 
  theme(legend.text = element_text(size = 10), legend.position = c(0.2, 0.9))

dev.off()




## GMM CEBPB ###############
library(Seurat)
library(mclust)
library(ggplot2)
library(dplyr)


# temp_seurat_epi <- readRDS("/home/Data_Drive_8TB_2/ColonCancer/human1-31/Epi_script/final_version/Scoring_and_Scaling_epi_with.batchcorrection_PC19_res0.3_exclude_nUMI.untilUMAP_2.rds")



gene_counts <- Epi_Cyling_CRC@assays$RNA@data["CEBPB", ]
gene_counts <- gene_counts[gene_counts > 0]


max(gene_counts)
min(gene_counts)

# gene_counts <- round(gene_counts, digits = 3)
gene_counts <- as.data.frame(gene_counts)

ggplot() +
  geom_histogram(data = gene_counts, aes(x = gene_counts), fill = "black", alpha = 0.5, bins = 100)


ggplot() +
  geom_histogram(data = gene_counts, aes(x = gene_counts, y= ..density..), fill = "black", alpha = 0.5, bins = 1000)

hist(gene_counts$gene_counts, main = "CEBPB", xlab = "Expression", ylab = "Frequency", breaks= 100)
hist_plot <- hist(gene_counts$gene_counts, main = "CEBPB", xlab = "Expression", ylab = "Frequency", breaks= 100)

# max(log2(gene_counts) )
# min(log2(gene_counts) ) 
# 
# hist(log2(gene_counts), main = "CEBPB", xlab = "Count", ylab = "Frequency")

sum(is.na(gene_counts))


data <- gene_counts
min(data)
max(data)


# GMM 모델 적합
fit <- Mclust(data, G = 2)

# 적합 결과 출력
summary(fit)

# GMM 모델 시각화
plot(fit, what = "classification")
max(fit$data[fit$classification == 1])
min(fit$data[fit$classification == 2])

fit$parameters

sig_value <- sqrt( fit$parameters$variance$sigmasq[1] )
mu_value<- fit$parameters$mean[1]

curve(dnorm(x, mean = mu_value, sd = sig_value) , 
      from = 0 , to = 5,
      main = "Fitted Gaussian Distribution",
      xlab = "Value", ylab = "Density")

curve(dnorm(x, mean = fit$parameters$mean[2], sd = fit$parameters$variance$sigmasq[2]) , 
      from = 0 , to = 5,
      main = "Fitted Gaussian Distribution",
      xlab = "Value", ylab = "Density")



# View(gene_counts)


# 필요한 패키지 불러오기
library(ggplot2)
library(dplyr)

# 정규 분포 A와 B의 파라미터 설정
mean_A <- mu_value
sd_A <- sig_value

mean_B <- fit$parameters$mean[2]
sd_B <- sqrt(fit$parameters$variance$sigmasq[2])

# 정규 분포 A와 B의 확률 밀도 함수(PDF) 정의
pdf_A <- function(x) dnorm(x, mean = mean_A, sd = sd_A)
pdf_B <- function(x) dnorm(x, mean = mean_B, sd = sd_B)

# 그래프 그리기
x <- seq(0.01, 5, by = 0.0001)
data <- data.frame(x = x, A = pdf_A(x), B = pdf_B(x))

ggplot(data, aes(x)) +
  geom_line(aes(y = A), color = "blue", linetype = "dashed") +
  geom_line(aes(y = B), color = "red", linetype = "dashed") +
  labs(title = "CEBPB GMM model fitting",x="Expression", y = "count") +
  theme(plot.title = element_text(hjust = 0.5)) 

# 그래프가 만나는 지점 찾기
intersection <- data %>%
  filter(x >= 0.5, abs(A - B) < 0.01) %>%
  head(1)

intersection$x  # 그래프가 만나는 첫 번째 x 값 (0.5 이상인 값)
intersection$A  # 정규 분포 A에서의 해당 x 값에 대한 y 값
intersection$B  # 정규 분포 B에서의 해당 x 값에 대한 y 값


ggplot() +
  geom_histogram(data = gene_counts, aes(x = gene_counts, y= (..density..)*2 ), fill = "black", alpha = 0.5, bins = 100) +
  geom_line(data = data , aes(x=x ,y = A), color = "blue", linetype = "dashed") +
  geom_line(data= data , aes(x=x, y = B), color = "red", linetype = "dashed") +
  geom_vline(xintercept = intersection$x , linetype = "solid", color = "#DAA520") +
  labs(title = "CEBPB GMM model fitting",x="Expression", y = "count") +
  theme(plot.title = element_text(hjust = 0.5)) 

ggplot() +
  geom_histogram(data = gene_counts, aes(x = gene_counts, y = (..density..)*2), fill = "black", alpha = 0.5, bins = 100) +
  geom_line(data = data, aes(x = x, y = A), color = "blue", linetype = "dashed") +
  geom_line(data = data, aes(x = x, y = B), color = "red", linetype = "dashed") +
  geom_vline(xintercept = intersection$x, linetype = "solid", color = "#DAA520") +
  annotate("text", x = intersection$x, y = 0, label = intersection$x, vjust = -0.25, color = "#DAA520", size=0) +
  labs(title = "CEBPB GMM model fitting", x = "Expression", y = "count") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 15))



### Fig1_sup  ###########
cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure1SUP_CEBPB_GMM.pdf" , width = 3, height = 2)

ggplot() +
  geom_histogram(data = gene_counts, aes(x = gene_counts, y= (..density..)*2 ), fill = "black", alpha = 0.5, bins = 100) +
  geom_line(data = data , aes(x=x ,y = A), color = "blue", linetype = "dashed") +
  geom_line(data= data , aes(x=x, y = B), color = "red", linetype = "dashed") +
  geom_vline(xintercept = intersection$x , linetype = "solid", color = "#DAA520") +
  labs(title = "CEBPB GMM model fitting",x="Normalized Expr", y = "density") +
  theme_classic()+ 
  xlim(0,3.9)+
  theme(plot.title = element_text(hjust = 0.5)) 

dev.off()







# Fig1B - T cell Profiling ###################

# https://ashpublications.org/blood/article/105/7/2877/20218/CCR6-expression-defines-regulatory-effector-memory
# https://www.researchgate.net/figure/Markers-of-memory-T-cell-subsets-and-their-precursors-Terminally-differentiated-TCM_fig3_339483745
# https://arthritis-research.biomedcentral.com/articles/10.1186/s13075-017-1343-8#Fig1

# IL6ST
# https://www.biocompare.com/Editorial-Articles/597618-A-Guide-to-Naive-T-Cell-Markers/
# https://www.rndsystems.com/product-highlights/antibodies-memory-t-cell-subset-identification

DimPlot(T_CRC_3_filtered, group.by = "seurat_clusters", label=T)

Idents(T_CRC_3_filtered) <- T_CRC_3_filtered$seurat_clusters

# sum(table(combined.integrated@meta.data$SCT_snn_res.0.2))
clusters <- T_CRC_3_filtered$seurat_clusters
cluster_names <- rep("Unknown", length(clusters))

cluster_names[clusters %in% c(0)] <- "0"
cluster_names[clusters %in% c(1)] <- "1"
cluster_names[clusters %in% c(2)] <- "2"
cluster_names[clusters %in% c(3)] <- "3"
cluster_names[clusters %in% c(4)] <- "4"
cluster_names[clusters %in% c(5)] <- "5"
cluster_names[clusters %in% c(6)] <- "6"
cluster_names[clusters %in% c(7)] <- "7"
cluster_names[clusters %in% c(8)] <- "8"
cluster_names[clusters %in% c(9)] <- "9"
cluster_names[clusters %in% c(10)] <- "10"
cluster_names[clusters %in% c(11)] <- "11"
cluster_names[clusters %in% c(12)] <- "12"
cluster_names[clusters %in% c(13)] <- "13"
cluster_names[clusters %in% c(14)] <- "14"
cluster_names[clusters %in% c(15)] <- "15"
cluster_names[clusters %in% c(16)] <- "16"
cluster_names[clusters %in% c(17)] <- "17"
cluster_names[clusters %in% c(18)] <- "18"
cluster_names[clusters %in% c(21)] <- "19"
cluster_names[clusters %in% c(22)] <- "20"
cluster_names[clusters %in% c(23)] <- "21"
cluster_names[clusters %in% c(25)] <- "22"
cluster_names[clusters %in% c(26)] <- "23"
cluster_names[clusters %in% c(27)] <- "24"
cluster_names[clusters %in% c(28)] <- "25"
cluster_names[clusters %in% c(29)] <- "26"

T_CRC_3_filtered <- AddMetaData(T_CRC_3_filtered, metadata = cluster_names, col.name = "Fig1B_2Sup_clusters")

Idents(T_CRC_3_filtered) <- T_CRC_3_filtered$Fig1B_2Sup_clusters
T_CRC_3_filtered@active.ident <- factor(T_CRC_3_filtered@active.ident , levels = seq(0,26)   )



## Fig1B-2 Sup Cluster UMAP ##############

T_CRC_3_filtered@active.ident <- factor(T_CRC_3_filtered@active.ident, 
                                        levels = seq(0,26) ) 


cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure1B1Sup_total_Cluster_UMAP.pdf" , width = 6, height = 4)
DimPlot(T_CRC_3_filtered ,label=T , label.size = 6 )+ 
  ggtitle(NULL)  + NoAxes() + theme(legend.text = element_text(margin = margin(t = 0, b = 0)) )
dev.off()



T_CRC_3_filtered_clusterDEG <- FindAllMarkers(T_CRC_3_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
T_CRC_3_filtered_clusterDEG <- T_CRC_3_filtered_clusterDEG[T_CRC_3_filtered_clusterDEG$p_val_adj <= 0.05, ]

T_CRC_3_filtered_clusterDEG %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 2) %>%
  ungroup() -> T_CRC_3_filtered_marker_topgene

T_CRC_3_filtered@active.ident <- factor(T_CRC_3_filtered@active.ident , levels = rev(seq(0,26))   )

DimPlot(T_CRC_3_filtered, group.by = "Fig1B_2Sup_clusters", label = T)

write.table(T_CRC_3_filtered_clusterDEG, "/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/Figure1_Gene_list/Filtered_T_seuratcluster_DEG_minpct0.25_logfc0.5_adjP0.05_comm27.txt", sep = "\t",quote = F)


## Fig1B-2 Sup clusterDEG ##############
cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure1B2Sup_t_DotPlot_clusterDEG_increasing.pdf" , width = 12, height = 6)

DotPlot(T_CRC_3_filtered, features =  unique(T_CRC_3_filtered_marker_topgene$gene)     , dot.scale = 8 ,  )  +   scale_colour_viridis(option="viridis")  +  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, family="Arial", face='bold' , size = 10 ), 
        axis.text.y = element_text(angle = 0 , size = 10, hjust = 0 ,  face='bold'), 
        strip.text.x = element_blank()
  )   + xlab('') +  ylab('')
dev.off()



# VlnPlot(T_CRC_3_filtered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# VlnPlot(T_CRC_3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = "seurat_clusters")
# DimPlot(T_CRC_3 , group.by = "seurat_clusters" , label=T)






Idents(T_CRC_3_filtered) <- T_CRC_3_filtered$Fig1B_2Sup_clusters
clusters <- T_CRC_3_filtered$Fig1B_2Sup_clusters
cluster_names <- T_CRC_3_filtered$Fig1B_2Sup_clusters
# cluster_names <- rep("Unknown", length(clusters))
# cluster_names <- clusters

cluster_names[clusters %in% c("6","7","11","20","24")] <- "CD4 Treg"

cluster_names[clusters %in% c("5")] <- "CD4 Th17"
cluster_names[clusters %in% c("12","4")] <- "CD4 Tcm"
cluster_names[clusters %in% c("8")] <- "CD4 Tfh"
cluster_names[clusters %in% c("2")] <- "CD4 Tnaive"

cluster_names[clusters %in% c("16")] <- "CD8 MAIT"
cluster_names[clusters %in% c("0","1","18")] <- "CD8 Tem"
cluster_names[clusters %in% c("21")] <- "CD8 Temra"
cluster_names[clusters %in% c("3","14","23")] <- "CD8 Tex"
cluster_names[clusters %in% c("15")] <- "CD8 Trm"


cluster_names[clusters %in% c("10","13")] <- "NK"


cluster_names[clusters %in% c("9","17","19","22")] <- "Unknown"
cluster_names[clusters %in% c("25","26")] <- "Unknown"


T_CRC_3_filtered <- AddMetaData(T_CRC_3_filtered, metadata = cluster_names, col.name = "Annotation_T_Fig1B")

# T_CRC_3_filtered$Annotation_T_Fig1B -> T_CRC_3_filtered$Annotation_T_Fig1B_v2
# T_CRC_3_filtered$Annotation_T_Fig1B_v2[ T_CRC_3_filtered$Annotation_T_Fig1B_v2 == "CD4 Treg" ] <- "CD4\nTreg"
# T_CRC_3_filtered$Annotation_T_Fig1B_v2[ T_CRC_3_filtered$Annotation_T_Fig1B_v2 == "CD4 Th17" ] <- "CD4\nTh17"
# T_CRC_3_filtered$Annotation_T_Fig1B_v2[ T_CRC_3_filtered$Annotation_T_Fig1B_v2 == "CD4 Tfh" ] <- "CD4\nTfh"
# T_CRC_3_filtered$Annotation_T_Fig1B_v2[ T_CRC_3_filtered$Annotation_T_Fig1B_v2 == "CD4 Tnaive" ] <- "CD4\nTnaive"
# T_CRC_3_filtered$Annotation_T_Fig1B_v2[ T_CRC_3_filtered$Annotation_T_Fig1B_v2 == "CD4 Tcm" ] <- "CD4\nTcm"
# 
# T_CRC_3_filtered$Annotation_T_Fig1B_v2[ T_CRC_3_filtered$Annotation_T_Fig1B_v2 == "CD8 MAIT" ] <-  "CD8\nMAIT" 
# T_CRC_3_filtered$Annotation_T_Fig1B_v2[ T_CRC_3_filtered$Annotation_T_Fig1B_v2 == "CD8 Tem" ] <-  "CD8\nTem" 
# T_CRC_3_filtered$Annotation_T_Fig1B_v2[ T_CRC_3_filtered$Annotation_T_Fig1B_v2 == "CD8 Temra" ] <-  "CD8\nTemra" 
# T_CRC_3_filtered$Annotation_T_Fig1B_v2[ T_CRC_3_filtered$Annotation_T_Fig1B_v2 == "CD8 Tex" ] <-  "CD8\nTex" 
# T_CRC_3_filtered$Annotation_T_Fig1B_v2[ T_CRC_3_filtered$Annotation_T_Fig1B_v2 == "CD8 Trm" ] <- "CD8\nTrm"

# Idents(T_CRC_3_filtered) <- T_CRC_3_filtered$Annotation_T_Fig1B
# T_CRC_3_filtered@active.ident <- factor(T_CRC_3_filtered@active.ident , levels = rev(c("Unknown", "NK", 
#                                                                                        "CD8\nTrm" ,"CD8\nTemra" , "CD8\nTem" , "CD8\nMAIT" ,  "CD8\nTex" , 
#                                                                                        "CD4\nTcm" ,"CD4\nTnaive", "CD4\nTfh" , "CD4\nTh17", "CD4\nTreg"   ))    )


Idents(T_CRC_3_filtered)  <- T_CRC_3_filtered$Annotation_T_Fig1B
T_CRC_3_filtered@active.ident <- factor(T_CRC_3_filtered@active.ident, 
                                        levels=rev(c("Unknown", "NK", 
                                                 "CD8 Trm" ,"CD8 Temra" , "CD8 Tem" , "CD8 MAIT" ,  "CD8 Tex" , 
                                                 "CD4 Tcm" ,"CD4 Tnaive", "CD4 Tfh" , "CD4 Th17", "CD4 Treg" )  ) )


DimPlot(T_CRC_3_filtered ,label=T , label.size = 7 , cols =     c("#FFD700", "#00CED1", "#DC143C", "#7B68EE", "#228B22", "#FF4500", "#32CD32", "#663399", "#FF1493", "#00FA9A" , "#4682B4", "#F961DD")     )   + NoAxes()

# rev(c("Unknown", "NK","CD8 Tem" ,"CD8 Temra" , "CD8 MAIT" , "CD8 Trm" , "CD8 Tex" , "CD4 Tfh" ,"CD4 Tem", "CD4 Tcm" ,   "CD4 Tnaive", "CD4 Treg"   ))

## Fig1B-1 ###################
cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure1B1_T_umap.pdf" , width = 6, height = 4)
DimPlot(T_CRC_3_filtered ,label=F , label.size = 5 , cols =     c("#FFD700", "#00CED1", "#DC143C", "#7B68EE", "#228B22", "#FF4500", "#32CD32", "#663399", "#FF1493", "#00FA9A" , "#4682B4", "#F961DD")     )   + NoAxes()

dev.off()



Idents(T_CRC_3_filtered) <- T_CRC_3_filtered$Annotation_T_Fig1B
T_CRC_3_filtered_Annotgroupmarker <- FindAllMarkers(T_CRC_3_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
T_CRC_3_filtered_Annotgroupmarker <- T_CRC_3_filtered_Annotgroupmarker[T_CRC_3_filtered_Annotgroupmarker$p_val_adj <= 0.05, ]

T_CRC_3_filtered_Annotgroupmarker %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 10) %>%
  ungroup() -> T_topgene

write.table(T_CRC_3_filtered_Annotgroupmarker, "/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/Figure1_Gene_list/Filtered_T_AnnotationGroup_DEG_minpct0.25_logfc0.5_adjP0.05_comm27.txt", sep = "\t",quote = F)





Idents(T_CRC_3_filtered)  <- T_CRC_3_filtered$Annotation_T_Fig1B
T_CRC_3_filtered@active.ident <- factor(T_CRC_3_filtered@active.ident, 
                                        levels=c("Unknown", "NK", 
                                                 "CD8 Trm" ,"CD8 Temra" , "CD8 Tem" , "CD8 MAIT" ,  "CD8 Tex" , 
                                                 "CD4 Tcm" ,"CD4 Tnaive", "CD4 Tfh" , "CD4 Th17", "CD4 Treg"   ) )
## Fig1B-2 ###################
# CD4 Treg # 6 7 11 20 24	FOXP3 IL2RA
# CD4 naive # 2	CCR7 SELL TCF7 LEF1
# CD4 Tcm # 4 12	FAS IL2RB
# CD4 Th17 # 5	CXCR6 KLRB1 RORA CCR6 FURIN CD6


# CD8 Tex # 3 14 23 	LAG3 HAVCR2 PDCD1 TIGIT 
# CD8 MAIT # 16	 TRAV1-2
# CD8 Tem # 0 1 18 19	GZMK
# CD8 Trm # 15	 ITGAE HOPX  XCL1+ CD8 T cell
# CD8 Temra # 21	 ITGB1 FGFBP2

# 9 17 22	MT+ T cell
# 25 26 Unknown

DotPlot(T_CRC_3_filtered, features = list( "CD4T" = c("CD4","FOXP3","KLRB1","RORA","CCR6","MAF","PDCD1","CCR7", "SELL","IL6ST" ,"FAS"), 
                                           "CD8T" = c("CD8A","HAVCR2","LAG3","TRAV1-2","GZMK","ITGB1","FGFBP2","ITGAE","HOPX" ), 
                                           "NK" = c("TYROBP","FCGR3A")
                                           
)  , 
dot.scale = 5, dot.min = 0.01  )  +   
  scale_colour_viridis(option="viridis")  +  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, family="Arial", face='bold' , size = 10 ), 
        axis.text.y = element_text(angle = 0 , size = 10, hjust = 0 ,  face='bold')
  )   + xlab('') +  ylab('')





cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure1B2_T_DotPlot.pdf" , width = 8, height = 4)

DotPlot(T_CRC_3_filtered, features = list( "CD4T" = c("CD4","FOXP3","KLRB1","RORA","MAF","PDCD1","CCR7", "SELL","IL6ST" ,"FAS"), 
                                           "CD8T" = c("CD8A","HAVCR2","LAG3","TRAV1-2","GZMK","ITGB1","FGFBP2","ITGAE","HOPX" ), 
                                           "NK" = c("TYROBP","FCGR3A")
                                           
)  , 
dot.scale = 5, dot.min = 0.01  )  +   
  scale_colour_viridis(option="viridis")  +  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, family="Arial", face='bold' , size = 10 ), 
        axis.text.y = element_text(angle = 0 , size = 10, hjust = 0 ,  face='bold')
  )   + xlab('') +  ylab('')

dev.off()


cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure1B2_T_DotPlot_simplever.pdf" , width = 7, height = 3.5)

DotPlot(T_CRC_3_filtered, features = list( "CD4T" = c("CD4","FOXP3","KLRB1","RORA","MAF","CCR7","IL6ST" ,"FAS"), 
                                           "CD8T" = c("CD8A","HAVCR2","LAG3","TRAV1-2","GZMK","FGFBP2","ITGAE" ), 
                                           "NK" = c("TYROBP","FCGR3A")
                                           
)  , 
dot.scale = 5, dot.min = 0.01  )  +   
  scale_colour_viridis(option="viridis")  +  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, family="Arial", face='bold' , size = 8 ), 
        axis.text.y = element_text(angle = 0 , size = 10, hjust = 0 ,  face='bold')
  )   + xlab('') +  ylab('')

dev.off()




## Fig1B-2 Sup###################

# T_CRC_3_filtered@active.ident <- factor(T_CRC_3_filtered@active.ident, 
#                                         levels=c("Unknown", "NK",
#                                                  "CD8 Trm" ,"CD8 Temra" , "CD8 Tem" , "CD8 MAIT" ,  "CD8 Tex" , 
#                                                  "CD4 Tcm" ,"CD4 Tnaive", "CD4 Tfh" , "CD4 Th17", "CD4 Treg"   ) )

Idents(T_CRC_3_filtered) <- T_CRC_3_filtered$Fig1B_2Sup_clusters
T_CRC_3_filtered@active.ident <- factor(T_CRC_3_filtered@active.ident , levels = c("6","7","11","20","24",
                                                                                   "5",
                                                                                   "8",
                                                                                   "2",
                                                                                   "4","12",
                                                                                   
                                                                                   "3","14","23",
                                                                                   "16",
                                                                                   "0","1","18",
                                                                                   "21",
                                                                                   "15",
                                                                                   
                                                                                   "10","13",
                                                                                   
                                                                                   "9","17","19","22",

                                                                                   "25","26"
)   )


T_CRC_3_filtered_clusterDEG <- FindAllMarkers(T_CRC_3_filtered, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
T_CRC_3_filtered_clusterDEG <- T_CRC_3_filtered_clusterDEG[T_CRC_3_filtered_clusterDEG$p_val_adj <= 0.05, ]

T_CRC_3_filtered_clusterDEG %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 2) %>%
  ungroup() -> T_CRC_3_filtered_clusterDEG_topgene





T_CRC_3_filtered@active.ident <- factor(T_CRC_3_filtered@active.ident , levels = rev( c("6","7","11","20","24",
                                                                                        "5",
                                                                                        "8",
                                                                                        "2",
                                                                                        "4","12",
                                                                                        
                                                                                        
                                                                                        "3","14","23",
                                                                                        "16",
                                                                                        "0","1","18",
                                                                                        "21",
                                                                                        "15",

                                                                                        "10","13",
                                                                                        
                                                                                        "9","17","19","22",
                                                                                        
                                                                                        "25","26")
)   )



DotPlot(T_CRC_3_filtered, features = unique(T_CRC_3_filtered_clusterDEG_topgene$gene)     , dot.scale = 8 ,  )  +   scale_colour_viridis(option="viridis")  +  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, family="Arial", face='bold' , size = 10 ), 
        axis.text.y = element_text(angle = 0 , size = 10, hjust = 0 ,  face='bold'), 
        strip.text.x = element_blank()
  )   + xlab('') +  ylab('')


### Fig1B-2Sup_clusterMarker ########################
cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure1B2Sup_t_DotPlot_clusterDEG_Stepsorted.pdf" , width = 12, height = 6)

DotPlot(T_CRC_3_filtered, features = unique(T_CRC_3_filtered_clusterDEG_topgene$gene)     , dot.scale = 8 ,  )  +   scale_colour_viridis(option="viridis")  +  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, family="Arial", face='bold' , size = 10 ), 
        axis.text.y = element_text(angle = 0 , size = 10, hjust = 0 ,  face='bold'), 
        strip.text.x = element_blank()
  )   + xlab('') +  ylab('')



dev.off()



library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(circlize)
library(ComplexHeatmap)
library(grid)
library(viridis)
library(Polychrome)


# DotPlot 생성
p2 <- DotPlot(T_CRC_3_filtered, features = unique(T_CRC_3_filtered_clusterDEG_topgene$gene)     , dot.scale = 8 ,  )  +   scale_colour_viridis(option="viridis")  +  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, family="Arial", face='bold' , size = 10 ), 
        axis.text.y = element_text(angle = 0 , size = 10, hjust = 0 ,  face='bold'), 
        strip.text.x = element_blank()
  )   + xlab('') +  ylab('')


df2 <- p2$data
head(df2)

# exp_mat 전치
exp_mat <- df2 %>%
  dplyr::select(-pct.exp, -avg.exp) %>%
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>%
  as.data.frame()

row.names(exp_mat) <- exp_mat$features.plot
exp_mat <- exp_mat[,-1] %>% as.matrix()
exp_mat <- as.data.frame(t(exp_mat)) # 전치
exp_mat <- exp_mat[nrow(exp_mat):1, ]

head(exp_mat)

# percent_mat 전치
percent_mat <- df2 %>%
  dplyr::select(-avg.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = id, values_from = pct.exp) %>%
  as.data.frame()

row.names(percent_mat) <- percent_mat$features.plot
percent_mat <- percent_mat[,-1] %>% as.matrix()
percent_mat <- as.data.frame(t(percent_mat)) # 전치
percent_mat <- percent_mat[nrow(percent_mat):1, ]

head(percent_mat)

breaks <- seq(min(df2$avg.exp.scaled), max(df2$avg.exp.scaled), length.out = 30)
colors <- viridis(length(breaks))
col_fun <- circlize::colorRamp2(breaks, colors)




# 그룹 정의
row_groups <- factor(c("CD4 Treg", "CD4 Treg","CD4 Treg","CD4 Treg",  "CD4 Treg" ,
                       "CD4 Th17",
                       "CD4 Tfh" ,
                       "CD4 Tnaive",
                       "CD4 Tcm" ,"CD4 Tcm" ,
                       
                       "CD8 Tex","CD8 Tex","CD8 Tex",
                       "CD8 MAIT",
                       "CD8 Tem" ,"CD8 Tem" ,"CD8 Tem" ,
                       "CD8 Temra" ,
                       "CD8 Trm" ,
                       
                       
                       "NK","NK",
                       
                       
                       "Unknown", "Unknown", "Unknown", "Unknown", 
                       
                       "Unknown", "Unknown")
                     
                     ,  levels = rev(c("Unknown", "NK",
                                       "CD8 Trm" ,"CD8 Temra" , "CD8 Tem" , "CD8 MAIT" ,  "CD8 Tex" , 
                                       "CD4 Tcm" ,"CD4 Tnaive", "CD4 Tfh" , "CD4 Th17", "CD4 Treg"   ))
)






# 행 annotation 생성
row_anno <- rowAnnotation(Group = row_groups, 
                          col = list(Group = c("CD4 Treg" = "#FFD700",  
                                               "CD4 Th17" = "#00CED1", 
                                               "CD4 Tfh" = "#DC143C", 
                                               "CD4 Tnaive" = "#7B68EE",
                                               "CD4 Tcm" = "#228B22",
                                               
                                               "CD8 Tex" = "#FF4500",
                                               "CD8 MAIT" = "#32CD32",
                                               "CD8 Tem" = "#663399",
                                               "CD8 Temra"= "#FF1493",
                                               "CD8 Trm"= "#00FA9A",
                                               
                                               "NK"= "#4682B4",
                                               "Unknown" = "#F961DD")),
                          
                          
                          
                          annotation_label = "CellType", show_legend=F,
                          show_annotation_name = T)

# Modified cell_fun with grid lines behind circles
cell_fun <- function(j, i, x, y, w, h, fill){
  # Draw horizontal lines
  grid.lines(x = unit(c(x - w/2, x + w/2), "npc"),
             y = unit(c(y, y), "npc"), 
             gp = gpar(col = "#EBEBEB", lty = "solid"))
  
  # Draw vertical lines
  # Only draw once per row
  grid.lines(x = unit(c(x, x), "npc"),
             y = unit(c(y - h/2, y + h/2), "npc"),
             gp = gpar(col = "#EBEBEB", lty = "solid"))
  
  # Draw the circle for each cell
  grid.circle(x = x, y = y, 
              r = percent_mat[i, j]/100 * min(unit.c(w, h)) * 0.75,
              gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))
}
ht <- Heatmap(as.matrix(exp_mat), 
              cluster_rows = FALSE, cluster_columns = FALSE,
              heatmap_legend_param = list(title = "Average Expression"),
              show_row_names = TRUE, row_names_side = "left", row_names_rot = 0, row_names_gp = gpar(fontsize = 10, fontfamily = "Arial"),
              show_column_names = TRUE, column_names_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Arial"),
              row_title = NULL,
              col = col_fun, 
              cell_fun = cell_fun,
              rect_gp = gpar(type = "none"),
              border = "black",
              show_heatmap_legend = FALSE,
              left_annotation = row_anno)

# Corrected legends
annotation_legend <- Legend(labels = rev( c("Unknown", "NK", 
                                            "CD8 Trm" ,"CD8 Temra" , "CD8 Tem" , "CD8 MAIT" ,  "CD8 Tex" , 
                                            "CD4 Tcm" ,"CD4 Tnaive", "CD4 Tfh" , "CD4 Th17", "CD4 Treg"  ) ),
                            title = "CellType",
                            legend_gp = gpar(fill = c("#FFD700", "#00CED1", "#DC143C", "#7B68EE", "#228B22", "#FF4500", "#32CD32", "#663399", "#FF1493", "#00FA9A" , "#4682B4", "#F961DD")      ))

percent_expressed_legend <- Legend(at = c(0, 25, 50, 75), 
                                   title = "Percent Expressed",
                                   legend_gp = gpar(col = "black"),
                                   type = "points", 
                                   size = unit(c(0.5, 1.25, 2.5, 5), "mm"),
                                   background = NULL) # Corrected

average_expression_legend <- Legend(at = c(-1, 0, 1, 2, 3),
                                    col_fun = col_fun,
                                    title = "Average Expression")

combined_legend <- packLegend(annotation_legend, average_expression_legend, percent_expressed_legend,  direction = "vertical")

# Draw heatmap with corrected legends
draw(ht, annotation_legend_side = "right", annotation_legend_list = combined_legend, padding = unit(c(2, 2, 2, 2), "mm"))


### Fig1B-2Sup_clusterDEG_Annotation ########################
cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure1B2Sup_t_DotPlot_clusterDEG_Annotation.pdf" , width = 12, height = 6)

draw(ht, annotation_legend_side = "right", annotation_legend_list = combined_legend, padding = unit(c(3, 3, 3, 3), "mm"))

dev.off()







### Fig1B-2Sup_clusterMarker ########################
Idents(T_CRC_3_filtered) <- T_CRC_3_filtered$Fig1B_2Sup_clusters
T_CRC_3_filtered@active.ident <- factor(T_CRC_3_filtered@active.ident , levels = rev(c("6","7","11","20","24",
                                                                                   "5",
                                                                                   "8",
                                                                                   "2",
                                                                                   "4","12",
                                                                                   
                                                                                   "3","14","23",
                                                                                   "16",
                                                                                   "0","1","18",
                                                                                   "21",
                                                                                   "15",
                                                                                   
                                                                                   "10","13",
                                                                                   
                                                                                   "9","17","19","22",
                                                                                   
                                                                                   "25","26")
)   )


cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure1B2Sup_t_DotPlot_clustermarker_mergeTnTcm.pdf" , width = 12, height = 6)    

DotPlot(T_CRC_3_filtered, features = list( "CD4T" = c("CD4"),
                                           "CD4 Treg" = c("FOXP3","CTLA4"),
                                           "CD4 Th17" = c("KLRB1","RORA","CCR6"),
                                           "CD4 Tfh" = c("MAF","PDCD1"),
                                           "CD4 Tnaive/Tcm" = c("CCR7", "SELL","IL6ST","IFNGR1","FAS","TIMP1"),
                                           
                                           
                                           "CD8T" = c("CD8A","CD8B"),
                                           "CD8 Tex" = c("HAVCR2","LAG3"),
                                           "CD8 MAIT" = c("TRAV1-2","SLC4A10"),
                                           "CD8 Tem" = c("GZMK","GZMM"),
                                           "CD8 Temra" = c("ITGB1","FGFBP2"),
                                           "CD8 Trm" = c("ITGAE","HOPX" ),
                                           

                                           "NK" = c("TYROBP","FCGR3A") 
                                           
)  , dot.scale = 5, dot.min = 0.01  )  +   
  scale_colour_viridis(option="viridis")  +  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, family="Arial", face='bold' , size = 10 ), 
        axis.text.y = element_text(angle = 0 , size = 10, hjust = 0 ,  face='bold')
  )   + xlab('') +  ylab('')


dev.off()



cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure1B2Sup_t_DotPlot_clustermarker_seperTnTcm.pdf" , width = 12, height = 6)    

DotPlot(T_CRC_3_filtered, features = list( "CD4T" = c("CD4"),
                                           "CD4 Treg" = c("FOXP3","CTLA4"),
                                           "CD4 Th17" = c("KLRB1","RORA","CCR6"),
                                           "CD4 Tfh" = c("MAF","PDCD1"),
                                           "CD4 Tnaive" = c("CCR7", "SELL","IL6ST" ),
                                           "CD4 Tcm" = c("IFNGR1","FAS","TIMP1"), 
                                           
                                           
                                           
                                           "CD8T" = c("CD8A","CD8B"),
                                           "CD8 Tex" = c("HAVCR2","LAG3"),
                                           "CD8 MAIT" = c("TRAV1-2","SLC4A10"),
                                           "CD8 Tem" = c("GZMK","GZMM"),
                                           "CD8 Temra" = c("ITGB1","FGFBP2"),
                                           "CD8 Trm" = c("ITGAE","HOPX" ),
                                           
                                           "NK" = c("TYROBP","FCGR3A")
                                           
)  , 
dot.scale = 5, dot.min = 0.01  )  +   
  scale_colour_viridis(option="viridis")  +  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, family="Arial", face='bold' , size = 10 ), 
        axis.text.y = element_text(angle = 0 , size = 10, hjust = 0 ,  face='bold')
  )   + xlab('') +  ylab('')


dev.off()




library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(circlize)
library(ComplexHeatmap)
library(grid)
library(viridis)
library(Polychrome)


# DotPlot 생성
p2 <- DotPlot(T_CRC_3_filtered, features =  c("CD4","FOXP3","CTLA4","KLRB1","RORA","CCR6","MAF","PDCD1","CCR7", "SELL","IL6ST" ,"IFNGR1","FAS","TIMP1",
                                              "CD8A","CD8B","HAVCR2","LAG3","TRAV1-2","SLC4A10", "GZMK","GZMM","ITGB1","FGFBP2","ITGAE","HOPX",
                                              "TYROBP","FCGR3A")
              , dot.scale = 8 ,  )  +   scale_colour_viridis(option="viridis")  +  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, family="Arial", face='bold' , size = 10 ), 
        axis.text.y = element_text(angle = 0 , size = 10, hjust = 0 ,  face='bold'), 
        strip.text.x = element_blank()
  )   + xlab('') +  ylab('')


df2 <- p2$data
head(df2)

# exp_mat 전치
exp_mat <- df2 %>%
  dplyr::select(-pct.exp, -avg.exp) %>%
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>%
  as.data.frame()

row.names(exp_mat) <- exp_mat$features.plot
exp_mat <- exp_mat[,-1] %>% as.matrix()
exp_mat <- as.data.frame(t(exp_mat)) # 전치
exp_mat <- exp_mat[nrow(exp_mat):1, ]

head(exp_mat)

# percent_mat 전치
percent_mat <- df2 %>%
  dplyr::select(-avg.exp, -avg.exp.scaled) %>%
  pivot_wider(names_from = id, values_from = pct.exp) %>%
  as.data.frame()

row.names(percent_mat) <- percent_mat$features.plot
percent_mat <- percent_mat[,-1] %>% as.matrix()
percent_mat <- as.data.frame(t(percent_mat)) # 전치
percent_mat <- percent_mat[nrow(percent_mat):1, ]

head(percent_mat)

breaks <- seq(min(df2$avg.exp.scaled), max(df2$avg.exp.scaled), length.out = 30)
colors <- viridis(length(breaks))
col_fun <- circlize::colorRamp2(breaks, colors)




row_groups <- factor(c("CD4 Treg", "CD4 Treg","CD4 Treg","CD4 Treg",  "CD4 Treg" ,
                       "CD4 Th17",
                       "CD4 Tfh" ,
                       "CD4 Tnaive",
                       "CD4 Tcm" ,"CD4 Tcm" ,
                       
                       "CD8 Tex","CD8 Tex","CD8 Tex",
                       "CD8 MAIT",
                       "CD8 Tem" ,"CD8 Tem" ,"CD8 Tem" ,
                       "CD8 Temra" ,
                       "CD8 Trm" ,
                       

                       "NK","NK",
                       
                       "Unknown", "Unknown", "Unknown", "Unknown",
                       "Unknown", "Unknown")
                     
                     ,  levels = rev(c("Unknown", "NK",
                                       "CD8 Trm" ,"CD8 Temra" , "CD8 Tem" , "CD8 MAIT" ,  "CD8 Tex" , 
                                       "CD4 Tcm" ,"CD4 Tnaive", "CD4 Tfh" , "CD4 Th17", "CD4 Treg"   ))
)



# 행 annotation 생성
row_anno <- rowAnnotation(Group = row_groups, 
                          col = list(Group = c("CD4 Treg" = "#FFD700",  
                                               "CD4 Th17" = "#00CED1", 
                                               "CD4 Tfh" = "#DC143C", 
                                               "CD4 Tnaive" = "#7B68EE",
                                               "CD4 Tcm" = "#228B22",
                                               
                                               "CD8 Tex" = "#FF4500",
                                               "CD8 MAIT" = "#32CD32",
                                               "CD8 Tem" = "#663399",
                                               "CD8 Temra"= "#FF1493",
                                               "CD8 Trm"= "#00FA9A",
                                               
                                               "NK"= "#4682B4",
                                               "Unknown" = "#F961DD")),
                          
                          
                          annotation_label = "CellType", show_legend=F,
                          show_annotation_name = T)


column_groups <- factor(c(rep("A", length.out = 1), 
                          rep("B", length.out = 2), 
                          rep("C", length.out = 3), 
                          rep("D", length.out = 2), 
                          rep("E", length.out = 3), 
                          rep("F", length.out = 3), 
                          
                          rep("G", length.out = 2), 
                          rep("H", length.out = 2), 
                          rep("I", length.out = 2), 
                          rep("J", length.out = 2), 
                          rep("K", length.out = 2), 
                          rep("L", length.out = 2), 
                          rep("M", length.out = 2)
                          
                          ))



# Modified cell_fun with grid lines behind circles
cell_fun <- function(j, i, x, y, w, h, fill){
  # Draw horizontal lines
  grid.lines(x = unit(c(x - w/2, x + w/2), "npc"),
             y = unit(c(y, y), "npc"), 
             gp = gpar(col = "#EBEBEB", lty = "solid"))
  
  # Draw vertical lines
  # Only draw once per row
  grid.lines(x = unit(c(x, x), "npc"),
             y = unit(c(y - h/2, y + h/2), "npc"),
             gp = gpar(col = "#EBEBEB", lty = "solid"))
  
  # Draw the circle for each cell
  grid.circle(x = x, y = y, 
              r = percent_mat[i, j]/100 * max(unit.c(w, h)) * 0.3,
              gp = gpar(fill = col_fun(exp_mat[i, j]), col = NA))
}

ht <- Heatmap(as.matrix(exp_mat), 
              cluster_rows = FALSE, cluster_columns = FALSE,
              heatmap_legend_param = list(title = "Average Expression"),
              show_row_names = TRUE, row_names_side = "left", row_names_rot = 0, row_names_gp = gpar(fontsize = 10, fontfamily = "Arial"),
              show_column_names = TRUE, column_names_gp = gpar(fontsize = 10, fontface = "bold", fontfamily = "Arial"),
              row_title = NULL,
              col = col_fun, 
              cell_fun = cell_fun,
              rect_gp = gpar(type = "none"),
              border = "black",
              show_heatmap_legend = FALSE,
              left_annotation = row_anno, 
              column_split = column_groups, 
              column_title = NULL)

# Corrected legends
annotation_legend <- Legend(labels = rev( c("Unknown", "NK", 
                                            "CD8 Trm" ,"CD8 Temra" , "CD8 Tem" , "CD8 MAIT" ,  "CD8 Tex" , 
                                            "CD4 Tcm" ,"CD4 Tnaive", "CD4 Tfh" , "CD4 Th17", "CD4 Treg"  ) ),
                            title = "CellType",
                            legend_gp = gpar(fill = c("#FFD700", "#00CED1", "#DC143C", "#7B68EE", "#228B22", "#FF4500", "#32CD32", "#663399", "#FF1493", "#00FA9A" , "#4682B4", "#F961DD")        ))




percent_expressed_legend <- Legend(at = c(0, 25, 50, 75), 
                                   title = "Percent Expressed",
                                   legend_gp = gpar(col = "black"),
                                   type = "points", 
                                   size = unit(c(0.5, 1.25, 2.5, 5), "mm"),
                                   background = NULL) # Corrected

average_expression_legend <- Legend(at = c(-1, 0, 1, 2, 3),
                                    col_fun = col_fun,
                                    title = "Average Expression")

combined_legend <- packLegend(annotation_legend, average_expression_legend, percent_expressed_legend,  direction = "vertical")

# Draw heatmap with corrected legends
draw(ht, annotation_legend_side = "right", annotation_legend_list = combined_legend, padding = unit(c(2, 2, 2, 2), "mm"))


### Fig1B-2Sup_clusterMarker_Annotation ########################
cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure1B2Sup_t_DotPlot_clustermarker_Annotation.pdf" , width = 12, height = 6)

draw(ht, annotation_legend_side = "right", annotation_legend_list = combined_legend, padding = unit(c(3, 3, 3, 3), "mm"))

dev.off()







##### T pex  marker ##############
DotPlot(T_CRC_3_filtered, features = list( "CD4T" = c("CD4"),
                                           "CD4 Treg" = c("FOXP3","CTLA4"),
                                           "CD4 Th17" = c("KLRB1","RORA","CCR6"),
                                           "CD4 Tfh" = c("MAF","PDCD1"),
                                           "CD4 Tnaive" = c("CCR7", "SELL","IL6ST" ),
                                           "CD4 Tcm" = c("IFNGR1","FAS","TIMP1"), 
                                           
                                           
                                           "Prolifer" = c("STMN1","MKI67","CCNA"),
                                           
                                           "CD8T" = c("CD8A","CD8B"),
                                           "CD8 Tpex" = c("TCF7","EOMES","CXCR5", "CD28", "IL7R"),
                                           "CD8 Tex" = c("HAVCR2","LAG3"),
                                           "CD8 MAIT" = c("TRAV1-2","SLC4A10"),
                                           "CD8 Tem" = c("GZMK","GZMM"),
                                           "CD8 Temra" = c("ITGB1","FGFBP2"),
                                           "CD8 Trm" = c("ITGAE","HOPX" ),
                                           
                                           "NK" = c("TYROBP","FCGR3A")
                                           
)  , 
dot.scale = 5, dot.min = 0.01  )  +   
  scale_colour_viridis(option="viridis")  +  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, family="Arial", face='bold' , size = 10 ), 
        axis.text.y = element_text(angle = 0 , size = 10, hjust = 0 ,  face='bold')
  )   + xlab('') +  ylab('')


DimPlot(T_CRC_3_filtered, label = T)
FeaturePlot(T_CRC_3_filtered, features = c("STMN1","MKI67","TCF7","LAG3","HAVCR2","CTLA4"), ncol=3  )

























# Fig1C - CRC total Proportion profiling  ###############################################################

# total_CRC$Annotation_Fig1 <- total_CRC$Annotation_total_Trough
# total_CRC$Annotation_Fig1[total_CRC$Annotation_Fig1 == "Cycling cells"] <- "Epithelial"
# total_CRC$Annotation_Fig1[total_CRC$Annotation_Fig1 == "Unused"] <- "Unknown"
# total_CRC$Annotation_Fig1[total_CRC$Annotation_Fig1 == "CD4_T"] <- "CD4T"
# total_CRC$Annotation_Fig1[total_CRC$Annotation_Fig1 == "CD8_T"] <- "CD8T"

total_CRC$Annotation_Fig1 -> total_CRC$Annotation_Fig1C
total_CRC$Annotation_Fig1C[total_CRC$Annotation_Fig1C == "plasmaB"] <- "Plasma"
total_CRC$Annotation_Fig1C[total_CRC$Annotation_Fig1C == "CD19B"] <- "B"



table(total_CRC$Annotation_Fig1C)


CellType_CRC <- total_CRC

# 데이터프레임 만들기
data_df <- data.frame(
  Patient = CellType_CRC$patients,     # 환자 정보
  CellType = CellType_CRC$Annotation_Fig1C  # 세포 타입 정보
)

mat_temp <- as.data.frame(table(data_df) )

mat <- as.data.frame(matrix(mat_temp$Freq, ncol = length(unique(mat_temp$Patient)), byrow = TRUE ,dimnames = list(colnames( table(data_df)  ),rownames(table(data_df) ) )  )  )
rm(data_df)
rm(mat_temp)



mat <- mat[, nPatients_list]
head(mat)
colSums(mat)

table(CellType_CRC$patients)


nCell_patients <- as.data.frame(t(as.data.frame(as.matrix(table(CellType_CRC$patients)))[1]) )
rownames(nCell_patients) <- "nCell"




# # CellType_CRC
# {
#   total_CRC
#   max(table(total_CRC$patients)) -> tmp_Pati_maxCell
#   
#   Vln_matrix <- matrix(0, nrow = tmp_Pati_maxCell , ncol = 30)
#   
#   
#   colnames(Vln_matrix) <- nPatients_list
#   as.data.frame(Vln_matrix)  -> Vln_matrix
#   
#   for (i in 1:length(nPatients_list) ){
#     Patients_n <- nPatients_list[i]
#     
#     tmp_Pati_cebpb = unlist( total_CRC@assays$RNA$data["CEBPB",][total_CRC$patients == Patients_n] ,use.names = F)
#     
#     Vln_matrix[i] <- c( tmp_Pati_cebpb, rep(NA, tmp_Pati_maxCell - length( tmp_Pati_cebpb )  ) )
#     
#   }
#   
#   
#   # draw(columnAnnotation(CEBPB = anno_density( Vln_matrix , type = "violin", which = "column"))  )
#   # p3 <- columnAnnotation(CEBPB = anno_density( Vln_matrix , type = "violin", which = "column"))
# }

mat_percent <- mat

for (col in 1:ncol(mat)) {
  col_name <- colnames(mat)[col]
  mat_percent[[col_name]] <- (mat[[col_name]] / nCell_patients[[col_name]]) * 100
}

# write.table(mat_percent, "/home/Data_Drive_8TB_3/ch2860/ColonCancer/Basic_seurat/Proportion_heatmap/Tseg_EpiCyc_percent.txt", quote = F , sep = "\t")

nCell_celltype <- data.frame(Celltype_nCell = rowSums(mat) )

head(nPatients_list)








colfunc <- colorRampPalette(c("white", "blue"))
colors <- colfunc(24)
lower_colors <- colors[5:20]

colfunc <- colorRampPalette(c("white", "red"))
colors <- colfunc(21)
upper_colors <- colors[4:17]


tmp_fr <- data.frame(index = patients_newnum , patient = nPatients_list, TMB = TMB)
tmp_fr[order(tmp_fr$TMB), ]
tmp_fr[order(tmp_fr$index), ]


tmp_fr <- cbind (tmp_fr[order(tmp_fr$TMB), ], "Color" = c(lower_colors,upper_colors))
tmp_fr
# tmp_fr[order(tmp_fr$index), ]["Color"]
# tmp_fr <- as.vector(tmp_fr[order(tmp_fr$index), ]["Color"])$Color




# 각 열의 벡터를 데이터 프레임으로 저장( 가장 위에 쓴 것이 가장 본 히트맵에 가깝게)
annotation_col <- data.frame(
  # "nCell`(Patients)`" = as.numeric(nCell_patients[,nPatients_list] ),
  # "CS_Ratio" = log2_CS_Ratio ,
  
  "Nstage" = N_stages,
  "Tstage" = T_stages ,
  "Location" = Locations  ,
  
  "RAS" = RAS_mutation_vector,
  "TP53" = TP53_mutation_vector,
  "APC" = APC_vector,
  "MSI" = MSS_vector 
  

  
  # "TMB_HL" = TMB_vector,
  # "TMB" = TMB
  # "Percent_TP53" = Percent_TP53,
  # "Percent_CEBPB" = Percent_CEBPB,
  # "CD8T_CTLA4_H_Cell" = as.numeric(CD8T_CTLA4_Patients_Mat[,nPatients_list]["H_pct",] ),
  # "CD4T_CTLA4_H_Cell" = as.numeric(CD4T_CTLA4_Patients_Mat[,nPatients_list]["H_pct",] ),
  # "CEBPB_H_Cell" = CEBPB_pct
  
  # KRAS_mutation = KRAS_mutation_vector,
  # NRAS_mutation = NRAS_mutation_vector,
  # BRAF_mutation = BRAF_mutation_vector
)

# colnames(annotation_col)[1] <- "nCell(Patients)"



# Vln_matrix ->ttttmmmppp 

row_anno_colors <- list(# CEBPB_H_Cell = c("white", "skyblue"),
  # CD4T_CTLA4_H_Cell = c("white" ,"#77AA77"),
  # CD8T_CTLA4_H_Cell = c("white" ,"#77AA77"),
  
  # Percent_CEBPB =  c("white", "#005D67"), 
  # Percent_TP53 = c("white","#A79FE1","#312271"),
  # TMB_HL = c("TMB_H" ="#2062B0", "TMB_L" = "#8098B0"),
  # "TMB" = setNames(as.character(tmp_fr$Color), tmp_fr$TMB), 
  
  
  
  
  "MSI" = c("MSS" = "darkgrey", "MSI-H" = "darkred"), 
  "RAS" = c("RAS mutation" = "purple", "WT" = "darkgrey"),
  "TP53" = c("TP53 WT" = "darkgrey", "TP53 mutation" = "orange"),
  "APC" = c("APC mutation" = "#B06D2B", "APC WT" = "darkgrey"),
  # CS_Ratio = c("#FF03FF", "#FF00C0" ,"white","#00FCC0"),
  # "nCell(Patients)" = c("white", "yellow")
  
  Location=c("D-colon" = "#98D8AA", "T-colon" = "#F3E99F", "A-colon" = "#F7D060", "S-colon" = "#FF6D60",
             "Rectum" = "#F6BA6F", "Cecum" = "#6DA9E4", "Rectosigmoid" = "blue"),
  "Tstage" = c("T2" = "#FFD93D", "T3" = "#FF8400", "T4" = "#4F200D"),
  "Nstage" = c("N0" = "#A7C04A", "N1" = "#539165", "N2" = "#3F497F")
)




col_fun = colorRamp2(c(0, 10, 25, 40, 50, 65, 80),  c("#FFEEE7", "#FFBEA9", "#FF7E67" ,"#E52833" ,"#A90335" ,"#80002C" ,"#6D0026" ))
lgd_Propor = Legend(col_fun = col_fun, title = "Proportion(%)", at = c(0, 25, 50, 75), 
                    labels = c("0", "25","50","75"))

col_fun = colorRamp2(c(0,30), c("white","skyblue") )
lgd_CEBPBH = Legend(col_fun = col_fun, title = "CEBPB-H(%)", at = c(0, 10, 20, 30),
                    labels = c("0", "10","20","30"))

col_fun = colorRamp2(c(1,2.5,2.51,4), c("#D2D2FF", "#2C2CFF", "#FFD8D8", "#FF3232"))
lgd_TMB = Legend(col_fun = col_fun, title = "TMB", at = c(1, 1.75, 2.5 ,3.25, 4), 
                 labels = c("0", "Low", "15","High","91"))

at_binary = seq(0, 1, by = 1)
col_fun_binary = colorRamp2(c(0, 1), c("darkgrey", "darkred"))
lgd_MSS = Legend(at = at_binary, title = "MSI", labels = c("MSS","MSI-H"),legend_gp = gpar(fill = col_fun_binary(at_binary)))

at_binary = seq(0, 1, by = 1)
col_fun_binary = colorRamp2(c(0, 1), c("darkgrey", "orange"))
lgd_TP53 = Legend(at = at_binary, title = "TP53", labels = c("TP53 WT","TP53 MT"),legend_gp = gpar(fill = col_fun_binary(at_binary)))

at_binary = seq(0, 2, by = 1)
col_fun_binary = colorRamp2(c(0, 1), c("darkgrey", "purple"))
lgd_RAS = Legend(at = at_binary, title = "RAS", labels = c("RAS WT","RAS MT"),legend_gp = gpar(fill = col_fun_binary(at_binary)))

at_binary = seq(0, 1, by = 1)
col_fun_binary = colorRamp2(c(0, 1), c("darkgrey", "#B06D2B"))
lgd_APC =Legend(at = at_binary, title = "APC", labels = c("APC WT","APC MT"),legend_gp = gpar(fill = col_fun_binary(at_binary)))



at_binary = seq(0, 6, by = 1)
col_fun_binary = colorRamp2(c(0, 1, 2, 3, 4, 5, 6), c("#98D8AA", "#F3E99F", "#F7D060", "#FF6D60",
                                        "#F6BA6F",  "#6DA9E4", "blue"))
lgd_Loc =Legend(at = at_binary, title = "Location", labels = c("D-colon" , "T-colon" , "A-colon" , "S-colon" ,
                                                               "Rectum" , "Cecum" , "Rectosigmoid"),legend_gp = gpar(fill = col_fun_binary(at_binary)))


at_binary = seq(0, 2, by = 1)
col_fun_binary = colorRamp2(c(0,1, 2), c("#FFD93D", "#FF8400", "#4F200D"))
lgd_StageT =Legend(at = at_binary, title = "Stage T", labels = c("T2","T3","T4"),legend_gp = gpar(fill = col_fun_binary(at_binary)))


at_binary = seq(0, 2, by = 1)
col_fun_binary = colorRamp2(c(0,1, 2), c("#A7C04A", "#539165", "#3F497F"))
lgd_StageN =Legend(at = at_binary, title = "Stage N", labels = c("N0","N1","N2"),legend_gp = gpar(fill = col_fun_binary(at_binary)))



col_fun = colorRamp2(c(-10,0,5), c("#FF03FF","white","#00FCC0") )
lgd_CS_logRatio = Legend(col_fun = col_fun, title = "CS_Ratio(log)")

col_fun = colorRamp2(c(0,11000), c("white","yellow") )
lgd_nCellP  = Legend(col_fun = col_fun, title = "nCell(Patients)", at = c(0, 10000),  legend_height = unit(1.5, "cm"),
                     labels = c("0", "10000"))



# pd = packLegend(list = list(lgd_Propor, lgd_CEBPBH, lgd_TMB, lgd_MSS, lgd_TP53, lgd_RAS, lgd_APC ,lgd_CS_logRatio, lgd_nCellP ))
pd = packLegend(lgd_MSS, lgd_APC , lgd_TP53, lgd_RAS, lgd_Propor,
                max_height = unit(18, "cm"), 
                column_gap = unit(1, "cm") )

pd = packLegend(lgd_MSS, lgd_APC , lgd_TP53, lgd_RAS, lgd_Propor, lgd_Loc, lgd_StageT, lgd_StageN ,
                max_height = unit(15, "cm"),
                column_gap = unit(1, "cm") )




# col_fun = colorRamp2(c(min(nCell_patients), max(nCell_patients)), c("white", "red"))
# 
# limmax <- floor( (max(nCell_patients)+ 1000) / 1000) * 1000
# lim_mid <- floor( (limmax/2) / 1000) * 1000
# 
# p0 = HeatmapAnnotation(nCells = unlist(nCell_patients,use.names = F) ,   col =  list(nCells  = col_fun) ,annotation_name_side = "left" ,annotation_label = "nCells", name = "nCells",show_annotation_name = T ,
#                        annotation_legend_param = list(title = "nCells-Patients",  at = c(0,lim_mid,limmax ), labels = c("0", as.character(lim_mid), as.character(limmax)))   )


p2 <- rowAnnotation(nCell = anno_barplot( nCell_celltype[1][c("Epithelial","CD4T","CD8T","Myeloid","B","Plasma","NK","Fibroblasts","Endothelial","Mast","Unknown"),], add_numbers = T, width = unit(30, "mm"), gp = gpar(border =NA, fill="gray" ),which = "row" ) )





quantiles_df <- function(df, probs = seq(0, 1, length.out = 1000)) {
  result <- as.data.frame(lapply(df, function(column) {
    quantile(column, probs = probs, na.rm = TRUE)
  }))
  return(result)
}




# p4 <- HeatmapAnnotation(CEBPB_PCT = anno_barplot( CEBPB_pct , bar_width = 1),
#                         CEBPB_EXP = anno_density( Vln_matrix , type = "violin", which = "column", width = unit(0, "cm"), height = unit(1.2, "cm") ,heatmap_colors = c("#FFFFFF","#BBBBBB")  )  ,
#                         
#                         show_annotation_name = T  ,annotation_name_side = "left" ,annotation_label = c("CEBPB(%)","CEBPB(exp)") , annotation_name_rot = 0 , annotation_name_gp = gpar(fontsize = 10),
#                         show_legend = F , annotation_legend_param = list(title = "CEBPB Exp level",  at = c(0,4 ), labels = c("0", "4") )
#                         
#                         
#                         
# )  

# p5 <- HeatmapAnnotation(CEBPB = anno_density( Vln_matrix , type = "heatmap", which = "column", width = unit(4, "cm"), height = unit(2, "cm")  ,joyplot_scale = 0.9 ) ,annotation_name_side = "left" , show_annotation_name = T  
#                         ,annotation_label = "CEBPB" , annotation_legend_param = list(title = "CEBPB Exp level",  at = c(0,2 ), labels = c("0", "2") ) ,  )  

# rownames(mat_percent) <- c("B", "CD4T", "CD8T" ,  "Endothelial", "Epithelial" ,"Fibroblasts" ,"Mast"  , "Myeloid" , "NK" ,  "Plasma"  ,   "Unknown"  )
# 
# mat_percent_tmp <- mat_percent[c("Epithelial","CD4T","CD8T","Myeloid","B","Plasma","NK","Fibroblasts","Endothelial","Mast","Unknown"),]
# 
# p1 <- pheatmap(as.matrix(mat_percent_tmp[,]), annotation_col = annotation_col, annotation_colors = row_anno_colors,  #scale = "col", 
#                # annotation_row = celltype_ncell[Index_row,][2],  
#                bottom_annotation = p4, 
#                annotation_legend = F, 
#                
#                
#                
#                # column_split = annotation_col$MSS, 
#                show_colnames = T,  column_title = "Patients",  column_title_side ="bottom", column_title_gp = gpar(fontsize = 12),
#                show_rownames = T,  row_names_side = 'left',  fontsize_row = 10 ,
#                
#                cluster_rows = F  , cluster_cols = F,
#                clustering_distance_cols = "correlation", # cutree_cols = 3 , #clustering_distance_cols = "euclidean",
#                
#                color = c("#FFEEE7", "#FFBEA9", "#FF7E67" ,"#E52833" ,"#A90335" ,"#80002C" ,"#6D0026" )
#                
#                # color = rev(hcl.colors(10, "Reds") )
#                
#                , name= "Proportion(%)"
#                , border_color = "gray" 
# )
# 
# 
# draw(p1 + p2 ,  annotation_legend_list= pd)
# 
# 
# # Legend(col_fun = c("#FFEEE7", "#FFBEA9", "#FF7E67" ,"#E52833" ,"#A90335" ,"#80002C" ,"#6D0026" ), title = "Proportion(%)",  at = c(0, 20 , 40, 60, 80))
# 
# VlnPlot(Epi_Cyling_CRC, features = "CEBPB" , group.by = "patients_num", pt.size = 0 , cols = rep("gray",30)) ->aaa


# png("/home/Data_Drive_8TB_3/Coloncancer/SequencingData/tmpfolder/figure_1000_600_res100_in.png", width = 10, height = 5.5 , res = 120,
#     type = "cairo" , units = "in" , antialias = "subpixel")
# draw(p1 + p2 ,  annotation_legend_list= pd)
# dev.off()
# 
# ggsave("/home/Data_Drive_8TB_3/Coloncancer/SequencingData/tmpfolder/figure_dpi50.png",
#        plot = draw(p1 + p2 ,  annotation_legend_list= pd), 
#        dpi=50, dev='png', height=4.5, width=8.5, units="in")
# dev.off()


# cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/figure1_Main_Proportion.pdf" , width = 10, height = 5.5)
# draw(p1 + p2 ,  annotation_legend_list= pd)
# dev.off()
# 
# 
# cairo_pdf("/home/Data_Drive_8TB_3/Coloncancer/SequencingData/tmpfolder/figure1_sub.pdf" , width = 13, height = 2)
# VlnPlot(Epi_Cyling_CRC, features = "CEBPB" , group.by = "patients_num", pt.size = 0 , cols = rep("gray",30))
# dev.off()
# 



rownames(mat_percent) <- c("B", "CD4T", "CD8T" ,  "Endothelial", "Epithelial" ,"Fibroblasts" ,"Mast"  , "Myeloid" , "NK" ,  "Plasma"  ,   "Unknown"  )

mat_percent_tmp <- mat_percent[c("Epithelial","CD4T","CD8T","Myeloid","B","Plasma","NK","Fibroblasts","Endothelial","Mast","Unknown"),]


p1 <- pheatmap(as.matrix(mat_percent_tmp[,]), annotation_col = annotation_col, annotation_colors = row_anno_colors,  #scale = "col", 
               # annotation_row = celltype_ncell[Index_row,][2],  
               # bottom_annotation = p4, 
               annotation_legend = F, legend = F, 
               
               
               
               # column_split = annotation_col$MSS, 
               show_colnames = T,  column_title = "Patients",  column_title_side ="bottom", column_title_gp = gpar(fontsize = 12),
               show_rownames = T,  row_names_side = 'left',  fontsize_row = 10 ,
               
               cluster_rows = F  , cluster_cols = F,
               clustering_distance_cols = "correlation", # cutree_cols = 3 , #clustering_distance_cols = "euclidean",
               
            
               color = c("#FFEEE7", "#FFBEA9", "#FF7E67" ,"#E52833" ,"#A90335" ,"#80002C" ,"#6D0026" )
               
               # color = rev(hcl.colors(10, "Reds") )
               
               , name= "Proportion(%)"
               , border_color = "gray" 
               
)

draw(p1 + p2 ,  annotation_legend_list= pd)


## Fig1C ##########
cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure1C_total_Proportion_NoCEBPB_updateMutation.pdf" , width = 10, height = 4)

draw(p1 + p2 ,  annotation_legend_list= pd)
dev.off()

# 
## Fig1C ##########
cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure1C_total_Proportion_NoCEBPB_updateMutation_for.pdf" , width = 10, height = 4.5)

draw(p1 + p2 ,  annotation_legend_list= pd)
dev.off()


# 





library(ggtext)

# Fig1D - Cell Proportion(total) wilcox test ##############
CellType_CRC <- total_CRC

# table(total_CRC$Annotation_Fig1)


# CellType_CRC <- subset(total_CRC , Annotation_Fig1 != "Unknown")

# 데이터프레임 만들기
data_df <- data.frame(
  Patient = CellType_CRC$patients,     # 환자 정보
  CellType = CellType_CRC$Annotation_Fig1C  # 세포 타입 정보
)

mat_temp <- as.data.frame(table(data_df) )

mat <- as.data.frame(matrix(mat_temp$Freq, ncol = length(unique(mat_temp$Patient)), byrow = TRUE ,dimnames = list(colnames( table(data_df)  ),rownames(table(data_df) ) )  )  )
rm(data_df)
rm(mat_temp)



mat <- mat[, nPatients_list]
head(mat)
colSums(mat)

table(CellType_CRC$patients)


nCell_patients <- as.data.frame(t(as.data.frame(as.matrix(table(CellType_CRC$patients)))[1]) )
rownames(nCell_patients) <- "nCells"


mat_percent <- mat

for (col in 1:ncol(mat)) {
  col_name <- colnames(mat)[col]
  mat_percent[[col_name]] <- (mat[[col_name]] / nCell_patients[[col_name]]) * 100
}

# write.table(mat_percent, "/home/Data_Drive_8TB_3/ch2860/ColonCancer/Basic_seurat/Proportion_heatmap/Tseg_EpiCyc_percent.txt", quote = F , sep = "\t")

nCell_celltype <- data.frame(Celltype_nCell = rowSums(mat) )

head(nPatients_list)


mat_percent_for_temp <- mat_percent





total_CRC$Annotation_Fig1_Immune <- total_CRC$Annotation_total_Trough 
total_CRC$Annotation_Fig1_Immune[total_CRC$Annotation_Fig1_Immune == "Cycling cells"] <- "Epithelial"
total_CRC$Annotation_Fig1_Immune[total_CRC$Annotation_Fig1_Immune == "Unused"] <- "Unknown"
total_CRC$Annotation_Fig1_Immune[total_CRC$Annotation_Fig1_Immune == "CD4_T"] <- "CD4T"
total_CRC$Annotation_Fig1_Immune[total_CRC$Annotation_Fig1_Immune == "CD8_T"] <- "CD8T"


total_CRC$Annotation_Fig1_Immune[total_CRC$Annotation_Fig1_Immune == "B"] <- "Immune Cell"
total_CRC$Annotation_Fig1_Immune[total_CRC$Annotation_Fig1_Immune == "CD4T"] <- "Immune Cell"
total_CRC$Annotation_Fig1_Immune[total_CRC$Annotation_Fig1_Immune == "CD8T"] <- "Immune Cell"
total_CRC$Annotation_Fig1_Immune[total_CRC$Annotation_Fig1_Immune == "Mast"] <- "Immune Cell"
total_CRC$Annotation_Fig1_Immune[total_CRC$Annotation_Fig1_Immune == "Myeloid"] <- "Immune Cell"
total_CRC$Annotation_Fig1_Immune[total_CRC$Annotation_Fig1_Immune == "NK"] <- "Immune Cell"
total_CRC$Annotation_Fig1_Immune[total_CRC$Annotation_Fig1_Immune == "Plasma"] <- "Immune Cell"
total_CRC$Annotation_Fig1_Immune[total_CRC$Annotation_Fig1_Immune == "Unknown"] <- "Immune Cell"

table(total_CRC$Annotation_Fig1_Immune )



CellType_CRC <- total_CRC

# 데이터프레임 만들기
data_df <- data.frame(
  Patient = CellType_CRC$patients,     # 환자 정보
  CellType = CellType_CRC$Annotation_Fig1_Immune  # 세포 타입 정보
)

mat_temp <- as.data.frame(table(data_df) )

mat <- as.data.frame(matrix(mat_temp$Freq, ncol = length(unique(mat_temp$Patient)), byrow = TRUE ,dimnames = list(colnames( table(data_df)  ),rownames(table(data_df) ) )  )  )
rm(data_df)
rm(mat_temp)



mat <- mat[, nPatients_list]
head(mat)
colSums(mat)

table(CellType_CRC$patients)


nCell_patients <- as.data.frame(t(as.data.frame(as.matrix(table(CellType_CRC$patients)))[1]) )
rownames(nCell_patients) <- "nCells"


mat_percent <- mat

for (col in 1:ncol(mat)) {
  col_name <- colnames(mat)[col]
  mat_percent[[col_name]] <- (mat[[col_name]] / nCell_patients[[col_name]]) * 100
}



mat_percent[c("Immune"),] -> mat_percent_for_temp[c("Immune"),]

# c("Epithelial","CD4T","CD8T","Myeloid","B","Plasma","NK","Fibroblasts","Endothelial","Mast","Unknown"),]

mat_percent_for_temp <- mat_percent_for_temp[c("Epithelial","CD4T","CD8T","Myeloid","B","Plasma","NK","Fibroblasts","Endothelial","Mast","Unknown","Immune"),]
rownames(mat_percent_for_temp) <- c("Epithelial","CD4T","CD8T","Myeloid","B","Plasma","NK","Fibroblasts","Endothelial","Mast","Unknown","Immune")


# Vec_TMB <- TMB_vector
Vec_MSS <- MSS_vector
Vec_APC <- ifelse(APC_vector == "APC mutation", "APC MT", APC_vector)
Vec_TP53 <- ifelse(TP53_mutation_vector == "TP53 mutation", "TP53 MT", TP53_mutation_vector)
Vec_RAS <- ifelse(RAS_mutation_vector == "RAS mutation", "RAS MT", RAS_mutation_vector)
Vec_RAS <- ifelse(Vec_RAS == "WT", "RAS WT", Vec_RAS)

# Split_info <- list(Vec_TMB, Vec_MSS, Vec_TP53, Vec_RAS, Vec_APC)  
# Split_Col <- c("TMB","MSS","TP53","RAS", "APC")
Split_info <- list( Vec_MSS, Vec_APC, Vec_TP53, Vec_RAS)  
Split_Col <- c("MSS", "APC","TP53","RAS")


tpval_result_frame <- data.frame(row.names = rownames(mat_percent_for_temp) ) #, MSS= c() ,TP53=c() ,  RAS= c()  )        


for (i in  1:length(Split_info))  {
  colnames(mat_percent_for_temp) <- Split_info[[i]]
  value_1 <- levels(as.factor(Split_info[[i]] ))[1]
  value_2 <- levels(as.factor(Split_info[[i]] ))[2]
  
  for (j in 1:length(rownames(mat_percent_for_temp)) ){
    vector_1 <- as.numeric(mat_percent_for_temp[colnames(mat_percent_for_temp) == value_1][j,])
    vector_2 <- as.numeric(mat_percent_for_temp[colnames(mat_percent_for_temp) == value_2][j,])
    
    # t.test(vector_1, vector_2)$p.value 
    tpval_result_frame[rownames(mat_percent_for_temp)[j], Split_Col[i]] <- wilcox.test(vector_1, vector_2)$p.value
    
  }
}


rev(names(table( Split_info[[1]] )))
names(table( Split_info[[2]] ))
names(table( Split_info[[3]] ))
names(table( Split_info[[4]] ))



mat_percent_for_temp
tpval_result_frame
# MSS_split_table_pct

# ggplot 객체들을 저장하는 리스트 생성
Vln_list <- list()
n_col = nrow(mat_percent_for_temp)

# 예시 데이터와 분할 정보를 가정한 상태에서 반복문을 통해 Vln_list를 채움
for (i in  1:length(Split_info))  {
  colnames(mat_percent_for_temp) <- Split_info[[i]]
  
  for (j in 1:length(rownames(mat_percent_for_temp)) ){
    
    if (j == 1 ){
      
      Vln_frame <- data.frame(
        Condition = colnames(mat_percent_for_temp) ,
        Proportion = unlist(mat_percent_for_temp[j,], use.names = F)
      )
      
      Vln_frame$Condition <- factor(Vln_frame$Condition, levels = rev(names(table( Split_info[[i]] ))) )
      
      # 첫 번째 행인지 확인
      is_first_row <- ((i-1) * length(rownames(mat_percent_for_temp)) + j) <= n_col
      
      # 첫 번째 행에만 제목과 x축 레이블을 추가
      title_text <- if(is_first_row) sprintf("%s", rownames(mat_percent_for_temp[j,])) else ""
      
      # p-value 계산
      p_value <- round(tpval_result_frame[rownames(mat_percent_for_temp[j,]), Split_Col[i]], digits = 3)
      
      # p-value가 0.05 미만이면 숫자를 빨갛게
      if (p_value < 0.05) {
        x_text <- sprintf("P = <span style='color:red;'>%s</span>", p_value)
      } else {
        x_text <- sprintf("P = %s", p_value)
      }
      
      
      p <- ggplot(Vln_frame, aes(x = Condition, y = Proportion, fill = Condition)) +
        geom_violin(trim = F) +
        labs(
          title = title_text,
          x = x_text,
          # x = "Condition",
          
          y = "Proportion(%)"
        ) + 
        ylim(0, 100) +
        theme_minimal() +  
        scale_fill_manual(values = c( "TP53 WT" = "skyblue",  "TP53 MT"= "salmon" ,  "RAS WT" = "skyblue",  "RAS MT"= "salmon" , 
                                      "APC WT" = "skyblue",  "APC MT"= "salmon" ,   "MSS"= "skyblue" , "MSI-H" = "salmon")) +       
        
        
        # coord_fixed(ratio = 0.02 )+ 
        
        theme(legend.position = "none",
              axis.text.x=element_blank(),
              plot.title = element_text(size = 20,hjust = 0.5),
              axis.title = element_text(size = 15),
              axis.title.x = element_markdown(size = 20),
              axis.title.y =  element_text(size = 18),
              axis.text.y = element_text(size = 15)
        )+ 
        geom_jitter(shape=16, position=position_jitter(0.0001))+ 
        stat_summary(fun=median, geom="point", size=2.5, color="red")
      
      
      Vln_list[[((i-1) * length(rownames(mat_percent_for_temp)) + j)]] <- p
    }
    
    else if (j == length(rownames(mat_percent_for_temp)) ){
      
      Vln_frame <- data.frame(
        Condition = colnames(mat_percent_for_temp) ,
        Proportion = unlist(mat_percent_for_temp[j,], use.names = F)
      )
      
      Vln_frame$Condition <- factor(Vln_frame$Condition, levels = rev(names(table( Split_info[[i]] ))) )
      
      
      # 첫 번째 행인지 확인
      is_first_row <- ((i-1) * length(rownames(mat_percent_for_temp)) + j) <= n_col
      
      # 첫 번째 행에만 제목과 x축 레이블을 추가
      title_text <- if(is_first_row) sprintf("%s", rownames(mat_percent_for_temp[j,])) else ""

      # p-value 계산
      p_value <- round(tpval_result_frame[rownames(mat_percent_for_temp[j,]), Split_Col[i]], digits = 3)
      
      # p-value가 0.05 미만이면 숫자를 빨갛게
      if (p_value < 0.05) {
        x_text <- sprintf("P = <span style='color:red;'>%s</span>", p_value)
      } else {
        x_text <- sprintf("P = %s", p_value)
      }
      
      
      p <- ggplot(Vln_frame, aes(x = Condition, y = Proportion, fill = Condition)) +
        geom_violin(trim = FALSE) +
        labs(title = title_text, x = x_text, y = "") + 
        ylim(0, 100) +
        theme_minimal() +  
        scale_fill_manual(values = c( "TP53 WT" = "skyblue",  "TP53 MT"= "salmon" ,  "RAS WT" = "skyblue",  "RAS MT"= "salmon" , 
                                      "APC WT" = "skyblue",  "APC MT"= "salmon" ,   "MSS"= "skyblue" , "MSI-H" = "salmon")) +             
        
        # coord_fixed(ratio = 0.02 )+ 
        theme(legend.position = "none",
              axis.text.x = element_blank(),
              plot.title = element_text(size = 20, hjust = 0.5),
              axis.title = element_text(size = 15),
              axis.title.x = element_markdown(size = 20),
              axis.title.y =  element_text(size = 0),
              axis.text.y = element_text(size = 15)
        ) + 
        geom_jitter(shape = 16, position = position_jitter(0.0001)) + 
        stat_summary(fun = median, geom = "point", size = 2.5, color = "red")
      
      Vln_list[[((i-1) * length(rownames(mat_percent_for_temp)) + j)]] <- p
    }
    
    else if ( j != 1) {
      
      Vln_frame <- data.frame(
        Condition = colnames(mat_percent_for_temp) ,
        Proportion = unlist(mat_percent_for_temp[j,], use.names = F)
      )
      Vln_frame$Condition <- factor(Vln_frame$Condition, levels = rev(  names(table( Split_info[[i]] ))  ) )
      
      # 첫 번째 행인지 확인
      is_first_row <- ((i-1) * length(rownames(mat_percent_for_temp)) + j) <= n_col
      
      # 첫 번째 행에만 제목과 x축 레이블을 추가
      title_text <- if(is_first_row) sprintf("%s", rownames(mat_percent_for_temp[j,])) else ""

      
      # p-value 계산
      p_value <- round(tpval_result_frame[rownames(mat_percent_for_temp[j,]), Split_Col[i]], digits = 3)
      
      # p-value가 0.05 미만이면 숫자를 빨갛게
      if (p_value < 0.05) {
        x_text <- sprintf("P = <span style='color:red;'>%s</span>", p_value)
      } else {
        x_text <- sprintf("P = %s", p_value)
      }
      
      p <- ggplot(Vln_frame, aes(x = Condition, y = Proportion, fill = Condition)) +
        geom_violin(trim = FALSE) +
        labs(title = title_text, x = x_text, y = "") + 
        ylim(0, NA) +
        theme_minimal() +  
        scale_fill_manual(values = c( "TP53 WT" = "skyblue",  "TP53 MT"= "salmon" ,  "RAS WT" = "skyblue",  "RAS MT"= "salmon" , 
                                      "APC WT" = "skyblue",  "APC MT"= "salmon" ,   "MSS"= "skyblue" , "MSI-H" = "salmon")) +      
        # coord_fixed(ratio = 0.02 )+ 
        theme(legend.position = "none",
              axis.text.x = element_blank(),
              plot.title = element_text(size = 20, hjust = 0.5),
              axis.title = element_text(size = 15),
              axis.title.x = element_markdown(size = 20),
              axis.title.y =  element_text(size = 0),
              axis.text.y = element_text(size = 15)
        ) + 
        geom_jitter(shape = 16, position = position_jitter(0.0001)) + 
        stat_summary(fun = median, geom = "point", size = 2.5, color = "red")
      
      Vln_list[[((i-1) * length(rownames(mat_percent_for_temp)) + j)]] <- p
    }
    
    
  }
}

Vln_list[[1]]
Vln_list[[12]]
grid.arrange(grobs = Vln_list, ncol = n_col) # 리스트의 각 ggplot 객체들을 그리드 형태로 표시

table(T_CRC_3_filtered$Annotation_T_Fig1B)

setwd("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/")
total_width <- n_col * 170 # 각 플롯의 가로 크기를 150픽셀로 가정
total_height <- ceiling(length(Vln_list) / n_col) * 190 # 각 플롯의 세로 크기를 160픽셀로 가정
#png(filename = "/home/Data_Drive_8TB_3/Coloncancer/SequencingData/Proportion_stat_test/Immune_lineage.png", width = total_width, height = total_height, units = "in", res = 300)

## Fig1D ##########
cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure1D_Total_propor_wilcox_uploadMutation.pdf" , width = 24, height = 10 )
grid.arrange(grobs = Vln_list, ncol = n_col) # 리스트의 각 ggplot 객체들을 그리드 형태로 표시
dev.off()






# ggplot 객체들을 저장하는 리스트 생성
Vln_list_lgd <- list()

# Vec_TMB <- TMB_vector
Vec_MSS <- MSS_vector
Vec_APC <- ifelse(APC_vector == "APC mutation", "APC MT", APC_vector)
Vec_TP53 <- ifelse(TP53_mutation_vector == "TP53 mutation", "TP53 MT", TP53_mutation_vector)
Vec_RAS <- ifelse(RAS_mutation_vector == "RAS mutation", "RAS MT", RAS_mutation_vector)
Vec_RAS <- ifelse(Vec_RAS == "WT", "RAS WT", Vec_RAS)

# Split_info <- list(Vec_TMB, Vec_MSS, Vec_TP53, Vec_RAS, Vec_APC)  
# Split_Col <- c("TMB","MSS","TP53","RAS", "APC")
Split_info <- list( Vec_MSS, Vec_APC, Vec_TP53, Vec_RAS)  
Split_Col <- c("MSS", "APC","TP53","RAS")


mat_percent_for_temp[1,] -> mat_percent_for_temp_forlegend

# 예시 데이터와 분할 정보를 가정한 상태에서 반복문을 통해 Vln_list를 채움
for (i in  1:length(Split_info))  {
  colnames(mat_percent_for_temp_forlegend) <- Split_info[[i]]
  
  for (j in 1:length(rownames(mat_percent_for_temp_forlegend)) ){
    Vln_frame <- data.frame(
      Condition = colnames(mat_percent_for_temp_forlegend) ,
      Proportion = unlist(mat_percent_for_temp_forlegend[j,], use.names = F)
    )
    Vln_frame$Condition <- factor(Vln_frame$Condition, levels = rev(names(table( Split_info[[i]] ))) )
    
    # 첫 번째 행인지 확인
    is_first_row <- ((i-1) * length(rownames(mat_percent_for_temp_forlegend)) + j) <= 1
    
    # 첫 번째 행에만 제목과 x축 레이블을 추가
    title_text <- if(is_first_row) sprintf("%s", rownames(mat_percent_for_temp_forlegend[j,])) else ""

    # p-value 계산
    p_value <- round(tpval_result_frame[rownames(mat_percent_for_temp[j,]), Split_Col[i]], digits = 3)
    
    # p-value가 0.05 미만이면 숫자를 빨갛게
    if (p_value < 0.05) {
      x_text <- sprintf("P = <span style='color:red;'>%s</span>", p_value)
    } else {
      x_text <- sprintf("P = %s", p_value)
    }
    
    
    p <- ggplot(Vln_frame, aes(x = Condition, y = Proportion, fill = Condition)) +
      geom_violin(trim=FALSE, width = 0.5) +
      labs(title = title_text, x = x_text) +
      labs(
        title = title_text,
        x = x_text,
        # x = "Condition",
        
        y = "Proportion(%)"
      ) + 
      ylim(0, NA) +
      theme_minimal() +  
      scale_fill_manual(values = c( "TP53 WT" = "skyblue",  "TP53 MT"= "salmon" ,  "RAS WT" = "skyblue",  "RAS MT"= "salmon" , 
                                    "APC WT" = "skyblue",  "APC MT"= "salmon" ,   "MSS"= "skyblue" , "MSI-H" = "salmon")) +  
      theme(
        axis.text.x=element_blank(),
        plot.title = element_text(size = 18,hjust = 0.5),
        axis.title = element_text(size = 15),
        axis.title.x = element_markdown(size = 20),
        axis.title.y =  element_text(size = 12),
        axis.text.y = element_text(size = 12)
      )+ 
      geom_jitter(shape=16, position=position_jitter(0.0001))+ 
      stat_summary(fun=median, geom="point", size=2.5, color="red")
    
    Vln_list_lgd[[((i-1) * length(rownames(mat_percent_for_temp_forlegend)) + j)]] <- p
  }
}

Vln_list_lgd[[1]]

setwd("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/")
total_width <- 1 * 150 # 각 플롯의 가로 크기를 150픽셀로 가정
total_height <- ceiling(length(Vln_list_lgd) / 1) * 160 # 각 플롯의 세로 크기를 160픽셀로 가정
grid.arrange(grobs = Vln_list_lgd, ncol = 1) # 리스트의 각 ggplot 객체들을 그리드 형태로 표시


## Fig1D legend ##########

cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure1D_Total_propor_wilcox_legend.pdf" , width = 3, height = 10 )
grid.arrange(grobs = Vln_list_lgd, ncol = 1) # 리스트의 각 ggplot 객체들을 그리드 형태로 표시
dev.off()





# Fig1D Addtional T(total) cell Proportion Analysis ########################
CellType_CRC <- total_CRC


# 데이터프레임 만들기
data_df <- data.frame(
  Patient = CellType_CRC$patients,     # 환자 정보
  CellType = CellType_CRC$Annotation_total_Tsub  # 세포 타입 정보
)

mat_temp <- as.data.frame(table(data_df) )

mat <- as.data.frame(matrix(mat_temp$Freq, ncol = length(unique(mat_temp$Patient)), byrow = TRUE ,dimnames = list(colnames( table(data_df)  ),rownames(table(data_df) ) )  )  )
rm(data_df)
rm(mat_temp)



mat <- mat[, nPatients_list]
head(mat)
colSums(mat)

table(CellType_CRC$patients)


nCell_patients <- as.data.frame(t(as.data.frame(as.matrix(table(CellType_CRC$patients)))[1]) )
rownames(nCell_patients) <- "nCells"


mat_percent <- mat

for (col in 1:ncol(mat)) {
  col_name <- colnames(mat)[col]
  mat_percent[[col_name]] <- (mat[[col_name]] / nCell_patients[[col_name]]) * 100
}

# write.table(mat_percent, "/home/Data_Drive_8TB_3/ch2860/ColonCancer/Basic_seurat/Proportion_heatmap/Tseg_EpiCyc_percent.txt", quote = F , sep = "\t")

nCell_celltype <- data.frame(Celltype_nCell = rowSums(mat) )

head(nPatients_list)


mat_percent_for_temp <- mat_percent[c("CD4 Treg","CD4 Th17","CD4 Tfh","CD4 Tnaive","CD4 Tcm","CD8 Tex","CD8 MAIT","CD8 Tem","CD8 Temra","CD8 Trm","Unknown"),]


# mat_percent_for_temp["Unknown", ] + mat_percent_for_temp["Unknown(MT+)", ] -> mat_percent_for_temp["Unknown", ]
# 
# mat_percent_for_temp <- mat_percent_for_temp[!rownames(mat_percent_for_temp) %in% "Unknown(MT+)", ]
# 
# # 원하는 행들만 선택
# mat_percent_for_temp <- mat_percent_for_temp[c("CD4 Treg", "CD4 Th17", "CD4 Tfh", "CD4 Tnaive", "CD4 Tcm", "CD8 Tex", "CD8 MAIT", "CD8 Tem", "CD8 Temra", "CD8 Trm", "Unknown"), ]



# Vec_TMB <- TMB_vector
Vec_MSS <- MSS_vector
Vec_APC <- ifelse(APC_vector == "APC mutation", "APC MT", APC_vector)
Vec_TP53 <- ifelse(TP53_mutation_vector == "TP53 mutation", "TP53 MT", TP53_mutation_vector)
Vec_RAS <- ifelse(RAS_mutation_vector == "RAS mutation", "RAS MT", RAS_mutation_vector)
Vec_RAS <- ifelse(Vec_RAS == "WT", "RAS WT", Vec_RAS)

# Split_info <- list(Vec_TMB, Vec_MSS, Vec_TP53, Vec_RAS, Vec_APC)  
# Split_Col <- c("TMB","MSS","TP53","RAS", "APC")
Split_info <- list( Vec_MSS, Vec_APC, Vec_TP53, Vec_RAS)  
Split_Col <- c("MSS", "APC","TP53","RAS")


tpval_result_frame <- data.frame(row.names = rownames(mat_percent_for_temp) ) #, MSS= c() ,TP53=c() ,  RAS= c()  )        


for (i in  1:length(Split_info))  {
  colnames(mat_percent_for_temp) <- Split_info[[i]]
  value_1 <- levels(as.factor(Split_info[[i]] ))[1]
  value_2 <- levels(as.factor(Split_info[[i]] ))[2]
  
  for (j in 1:length(rownames(mat_percent_for_temp)) ){
    vector_1 <- as.numeric(mat_percent_for_temp[colnames(mat_percent_for_temp) == value_1][j,])
    vector_2 <- as.numeric(mat_percent_for_temp[colnames(mat_percent_for_temp) == value_2][j,])
    
    # t.test(vector_1, vector_2)$p.value 
    tpval_result_frame[rownames(mat_percent_for_temp)[j], Split_Col[i]] <- wilcox.test(vector_1, vector_2)$p.value
    
  }
}

mat_percent_for_temp
tpval_result_frame
# MSS_split_table_pct

# ggplot 객체들을 저장하는 리스트 생성
Vln_list <- list()
n_col = nrow(mat_percent_for_temp)

# 예시 데이터와 분할 정보를 가정한 상태에서 반복문을 통해 Vln_list를 채움
for (i in  1:length(Split_info))  {
  colnames(mat_percent_for_temp) <- Split_info[[i]]
  
  for (j in 1:length(rownames(mat_percent_for_temp)) ){
    
    if (j == 1 ){
      
      Vln_frame <- data.frame(
        Condition = colnames(mat_percent_for_temp) ,
        Proportion = unlist(mat_percent_for_temp[j,], use.names = F)
      )
      Vln_frame$Condition <- factor(Vln_frame$Condition, levels = rev(names(table( Split_info[[i]] ))) )
      
      # 첫 번째 행인지 확인
      is_first_row <- ((i-1) * length(rownames(mat_percent_for_temp)) + j) <= n_col
      
      # 첫 번째 행에만 제목과 x축 레이블을 추가
      title_text <- if(is_first_row) sprintf("%s", rownames(mat_percent_for_temp[j,])) else ""
      
      # p-value 계산
      p_value <- round(tpval_result_frame[rownames(mat_percent_for_temp[j,]), Split_Col[i]], digits = 3)
      
      # p-value가 0.05 미만이면 숫자를 빨갛게
      if (p_value < 0.05) {
        x_text <- sprintf("P = <span style='color:red;'>%s</span>", p_value)
      } else {
        x_text <- sprintf("P = %s", p_value)
      }
      
      
      p <- ggplot(Vln_frame, aes(x = Condition, y = Proportion, fill = Condition)) +
        geom_violin(trim = F) +
        labs(
          title = title_text,
          x = x_text,
          # x = "Condition",
          
          y = "Proportion(%)"
        ) + 
        # ylim(0, 100) +
        theme_minimal() +  
        scale_fill_manual(values = c( "TP53 WT" = "skyblue",  "TP53 MT"= "salmon" ,  "RAS WT" = "skyblue",  "RAS MT"= "salmon" , 
                                      "APC WT" = "skyblue",  "APC MT"= "salmon" ,   "MSS"= "skyblue" , "MSI-H" = "salmon")) +  
        # coord_fixed(ratio = 0.02 )+ 
        
        theme(legend.position = "none",
              axis.text.x=element_blank(),
              plot.title = element_text(size = 20,hjust = 0.5),
              axis.title = element_text(size = 15),
              axis.title.x = element_markdown(size = 20),
              axis.title.y =  element_text(size = 18),
              axis.text.y = element_text(size = 15)
        )+ 
        geom_jitter(shape=16, position=position_jitter(0.0001))+ 
        stat_summary(fun=median, geom="point", size=2.5, color="red")
      
      
      Vln_list[[((i-1) * length(rownames(mat_percent_for_temp)) + j)]] <- p
    }
    
  
    else  {
      
      Vln_frame <- data.frame(
        Condition = colnames(mat_percent_for_temp) ,
        Proportion = unlist(mat_percent_for_temp[j,], use.names = F)
      )
      
      Vln_frame$Condition <- factor(Vln_frame$Condition, levels = rev(names(table( Split_info[[i]] ))) )
      
      # 첫 번째 행인지 확인
      is_first_row <- ((i-1) * length(rownames(mat_percent_for_temp)) + j) <= n_col
      
      # 첫 번째 행에만 제목과 x축 레이블을 추가
      title_text <- if(is_first_row) sprintf("%s", rownames(mat_percent_for_temp[j,])) else ""

      # p-value 계산
      p_value <- round(tpval_result_frame[rownames(mat_percent_for_temp[j,]), Split_Col[i]], digits = 3)
      
      # p-value가 0.05 미만이면 숫자를 빨갛게
      if (p_value < 0.05) {
        x_text <- sprintf("P = <span style='color:red;'>%s</span>", p_value)
      } else {
        x_text <- sprintf("P = %s", p_value)
      }
      
      
      p <- ggplot(Vln_frame, aes(x = Condition, y = Proportion, fill = Condition)) +
        geom_violin(trim = FALSE) +
        labs(title = title_text, x = x_text, y = "") + 
        ylim(0, NA) +
        theme_minimal() +  
        scale_fill_manual(values = c( "TP53 WT" = "skyblue",  "TP53 MT"= "salmon" ,  "RAS WT" = "skyblue",  "RAS MT"= "salmon" , 
                                      "APC WT" = "skyblue",  "APC MT"= "salmon" ,   "MSS"= "skyblue" , "MSI-H" = "salmon")) +    
        # coord_fixed(ratio = 0.02 )+ 
        theme(legend.position = "none",
              axis.text.x = element_blank(),
              plot.title = element_text(size = 20, hjust = 0.5),
              axis.title = element_text(size = 15),
              axis.title.x = element_markdown(size = 20),
              axis.title.y =  element_text(size = 0),
              axis.text.y = element_text(size = 15)
        ) + 
        geom_jitter(shape = 16, position = position_jitter(0.0001)) + 
        stat_summary(fun = median, geom = "point", size = 2.5, color = "red")
      
      Vln_list[[((i-1) * length(rownames(mat_percent_for_temp)) + j)]] <- p
    }
    
    
  }
}

grid.arrange(grobs = Vln_list, ncol = n_col) # 리스트의 각 ggplot 객체들을 그리드 형태로 표시


setwd("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/")
total_width <- n_col * 170 # 각 플롯의 가로 크기를 150픽셀로 가정
total_height <- ceiling(length(Vln_list) / n_col) * 190 # 각 플롯의 세로 크기를 160픽셀로 가정
#png(filename = "/home/Data_Drive_8TB_3/Coloncancer/SequencingData/Proportion_stat_test/Immune_lineage.png", width = total_width, height = total_height, units = "in", res = 300)

## Fig1D Add - T(total) ##########
cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure1Dadd_Tintotal_propor_wilcox_uploadMutation.pdf" , width = 24, height = 10 )
grid.arrange(grobs = Vln_list, ncol = n_col) # 리스트의 각 ggplot 객체들을 그리드 형태로 표시
dev.off()







# Fig1D Addtional T cell Proportion Analysis ########################
CellType_CRC <- Tonly_CRC_filtered

# table(total_CRC$Annotation_Fig1)


# CellType_CRC <- subset(total_CRC , Annotation_Fig1 != "Unknown")

# 데이터프레임 만들기
data_df <- data.frame(
  Patient = CellType_CRC$patients,     # 환자 정보
  CellType = CellType_CRC$Annotation_T_Fig1B  # 세포 타입 정보
)

mat_temp <- as.data.frame(table(data_df) )

mat <- as.data.frame(matrix(mat_temp$Freq, ncol = length(unique(mat_temp$Patient)), byrow = TRUE ,dimnames = list(colnames( table(data_df)  ),rownames(table(data_df) ) )  )  )
rm(data_df)
rm(mat_temp)



mat <- mat[, nPatients_list]
head(mat)
colSums(mat)

table(CellType_CRC$patients)


nCell_patients <- as.data.frame(t(as.data.frame(as.matrix(table(CellType_CRC$patients)))[1]) )
rownames(nCell_patients) <- "nCells"


mat_percent <- mat

for (col in 1:ncol(mat)) {
  col_name <- colnames(mat)[col]
  mat_percent[[col_name]] <- (mat[[col_name]] / nCell_patients[[col_name]]) * 100
}

# write.table(mat_percent, "/home/Data_Drive_8TB_3/ch2860/ColonCancer/Basic_seurat/Proportion_heatmap/Tseg_EpiCyc_percent.txt", quote = F , sep = "\t")

nCell_celltype <- data.frame(Celltype_nCell = rowSums(mat) )

head(nPatients_list)


mat_percent_for_temp <- mat_percent


# c("Epithelial","CD4T","CD8T","Myeloid","B","Plasma","NK","Fibroblasts","Endothelial","Mast","Unknown"),]

mat_percent_for_temp <- mat_percent_for_temp[c("CD4 Treg","CD4 Th17","CD4 Tfh","CD4 Tnaive","CD4 Tcm","CD8 Tex","CD8 MAIT","CD8 Tem","CD8 Temra","CD8 Trm","Unknown"),]
# rownames(mat_percent_for_temp) <- c("Epithelial","CD4T","CD8T","Myeloid","B","Plasma","NK","Fibroblasts","Endothelial","Mast","Unknown","Immune")


# Vec_TMB <- TMB_vector
Vec_MSS <- MSS_vector
Vec_APC <- ifelse(APC_vector == "APC mutation", "APC MT", APC_vector)
Vec_TP53 <- ifelse(TP53_mutation_vector == "TP53 mutation", "TP53 MT", TP53_mutation_vector)
Vec_RAS <- ifelse(RAS_mutation_vector == "RAS mutation", "RAS MT", RAS_mutation_vector)
Vec_RAS <- ifelse(Vec_RAS == "WT", "RAS WT", Vec_RAS)

# Split_info <- list(Vec_TMB, Vec_MSS, Vec_TP53, Vec_RAS, Vec_APC)  
# Split_Col <- c("TMB","MSS","TP53","RAS", "APC")
Split_info <- list( Vec_MSS, Vec_APC, Vec_TP53, Vec_RAS)  
Split_Col <- c("MSS", "APC","TP53","RAS")


tpval_result_frame <- data.frame(row.names = rownames(mat_percent_for_temp) ) #, MSS= c() ,TP53=c() ,  RAS= c()  )        


for (i in  1:length(Split_info))  {
  colnames(mat_percent_for_temp) <- Split_info[[i]]
  value_1 <- levels(as.factor(Split_info[[i]] ))[1]
  value_2 <- levels(as.factor(Split_info[[i]] ))[2]
  
  for (j in 1:length(rownames(mat_percent_for_temp)) ){
    vector_1 <- as.numeric(mat_percent_for_temp[colnames(mat_percent_for_temp) == value_1][j,])
    vector_2 <- as.numeric(mat_percent_for_temp[colnames(mat_percent_for_temp) == value_2][j,])
    
    # t.test(vector_1, vector_2)$p.value 
    tpval_result_frame[rownames(mat_percent_for_temp)[j], Split_Col[i]] <- wilcox.test(vector_1, vector_2)$p.value
    
  }
}

mat_percent_for_temp
tpval_result_frame
# MSS_split_table_pct

# ggplot 객체들을 저장하는 리스트 생성
Vln_list <- list()
n_col = nrow(mat_percent_for_temp)

# 예시 데이터와 분할 정보를 가정한 상태에서 반복문을 통해 Vln_list를 채움
for (i in  1:length(Split_info))  {
  colnames(mat_percent_for_temp) <- Split_info[[i]]
  
  for (j in 1:length(rownames(mat_percent_for_temp)) ){
    
    if (j == 1 ){
      
      Vln_frame <- data.frame(
        Condition = colnames(mat_percent_for_temp) ,
        Proportion = unlist(mat_percent_for_temp[j,], use.names = F)
      )
      Vln_frame$Condition <- factor(Vln_frame$Condition, levels = rev(names(table( Split_info[[i]] ))) )
      
      # 첫 번째 행인지 확인
      is_first_row <- ((i-1) * length(rownames(mat_percent_for_temp)) + j) <= n_col
      
      # 첫 번째 행에만 제목과 x축 레이블을 추가
      title_text <- if(is_first_row) sprintf("%s", rownames(mat_percent_for_temp[j,])) else ""
      x_text <- sprintf("P = %s", round(tpval_result_frame[rownames(mat_percent_for_temp[j,]), Split_Col[i] ],digits = 3) )
      
      x_text_style <- if (round(tpval_result_frame[rownames(mat_percent_for_temp[j,]), Split_Col[i] ],digits = 3) < 0.05) {
        element_text(size = 20, color = "red", face = "bold")
      } else {
        element_text(size = 20)
      }
      
      
      p <- ggplot(Vln_frame, aes(x = Condition, y = Proportion, fill = Condition)) +
        geom_violin(trim = F) +
        labs(
          title = title_text,
          x = x_text,
          # x = "Condition",
          
          y = "Proportion(%)"
        ) + 
        # ylim(0, 100) +
        theme_minimal() +  
        scale_fill_manual(values = c( "TP53 WT" = "skyblue",  "TP53 MT"= "salmon" ,  "RAS WT" = "skyblue",  "RAS MT"= "salmon" , 
                                      "APC WT" = "skyblue",  "APC MT"= "salmon" ,   "MSS"= "skyblue" , "MSI-H" = "salmon")) +    
        # coord_fixed(ratio = 0.02 )+ 
        
        theme(legend.position = "none",
              axis.text.x=element_blank(),
              plot.title = element_text(size = 20,hjust = 0.5),
              axis.title = element_text(size = 15),
              axis.title.x = x_text_style,
              axis.title.y =  element_text(size = 18),
              axis.text.y = element_text(size = 15)
        )+ 
        geom_jitter(shape=16, position=position_jitter(0.0001))+ 
        stat_summary(fun=median, geom="point", size=2.5, color="red")
      
      
      Vln_list[[((i-1) * length(rownames(mat_percent_for_temp)) + j)]] <- p
    }
    
    
    else  {
      
      Vln_frame <- data.frame(
        Condition = colnames(mat_percent_for_temp) ,
        Proportion = unlist(mat_percent_for_temp[j,], use.names = F)
      )
      Vln_frame$Condition <- factor(Vln_frame$Condition, levels = rev(names(table( Split_info[[i]] ))) )
      
      # 첫 번째 행인지 확인
      is_first_row <- ((i-1) * length(rownames(mat_percent_for_temp)) + j) <= n_col
      
      # 첫 번째 행에만 제목과 x축 레이블을 추가
      title_text <- if(is_first_row) sprintf("%s", rownames(mat_percent_for_temp[j,])) else ""
      x_text <- sprintf("P = %s", round(tpval_result_frame[rownames(mat_percent_for_temp[j,]), Split_Col[i] ],digits = 3) )
      
      x_text_style <- if (round(tpval_result_frame[rownames(mat_percent_for_temp[j,]), Split_Col[i] ],digits = 3) < 0.05) {
        element_text(size = 20, color = "red", face = "bold")
      } else {
        element_text(size = 20)
      }
      
      
      p <- ggplot(Vln_frame, aes(x = Condition, y = Proportion, fill = Condition)) +
        geom_violin(trim = FALSE) +
        labs(title = title_text, x = x_text, y = "") + 
        ylim(0, NA) +
        theme_minimal() +  
        scale_fill_manual(values = c( "TP53 WT" = "skyblue",  "TP53 MT"= "salmon" ,  "RAS WT" = "skyblue",  "RAS MT"= "salmon" , 
                                      "APC WT" = "skyblue",  "APC MT"= "salmon" ,   "MSS"= "skyblue" , "MSI-H" = "salmon")) +    
        # coord_fixed(ratio = 0.02 )+ 
        theme(legend.position = "none",
              axis.text.x = element_blank(),
              plot.title = element_text(size = 20, hjust = 0.5),
              axis.title = element_text(size = 15),
              axis.title.x = x_text_style,
              axis.title.y =  element_text(size = 0),
              axis.text.y = element_text(size = 15)
        ) + 
        geom_jitter(shape = 16, position = position_jitter(0.0001)) + 
        stat_summary(fun = median, geom = "point", size = 2.5, color = "red")
      
      Vln_list[[((i-1) * length(rownames(mat_percent_for_temp)) + j)]] <- p
    }
    
    
  }
}

grid.arrange(grobs = Vln_list, ncol = n_col) # 리스트의 각 ggplot 객체들을 그리드 형태로 표시


setwd("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/")
total_width <- n_col * 170 # 각 플롯의 가로 크기를 150픽셀로 가정
total_height <- ceiling(length(Vln_list) / n_col) * 190 # 각 플롯의 세로 크기를 160픽셀로 가정
#png(filename = "/home/Data_Drive_8TB_3/Coloncancer/SequencingData/Proportion_stat_test/Immune_lineage.png", width = total_width, height = total_height, units = "in", res = 300)

## Fig1D Add - T ##########
cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure1Dadd_T_propor_wilcox_uploadMutation.pdf" , width = 24, height = 10 )
grid.arrange(grobs = Vln_list, ncol = n_col) # 리스트의 각 ggplot 객체들을 그리드 형태로 표시
dev.off()






# Fig1D Addtional CD4T cell Proportion Analysis ########################
CellType_CRC <- CD4T_CRC

# table(total_CRC$Annotation_Fig1)


# CellType_CRC <- subset(total_CRC , Annotation_Fig1 != "Unknown")

# 데이터프레임 만들기
data_df <- data.frame(
  Patient = CellType_CRC$patients,     # 환자 정보
  CellType = CellType_CRC$Annotation_T_Fig1B  # 세포 타입 정보
)

mat_temp <- as.data.frame(table(data_df) )

mat <- as.data.frame(matrix(mat_temp$Freq, ncol = length(unique(mat_temp$Patient)), byrow = TRUE ,dimnames = list(colnames( table(data_df)  ),rownames(table(data_df) ) )  )  )
rm(data_df)
rm(mat_temp)



mat <- mat[, nPatients_list]
head(mat)
colSums(mat)

table(CellType_CRC$patients)


nCell_patients <- as.data.frame(t(as.data.frame(as.matrix(table(CellType_CRC$patients)))[1]) )
rownames(nCell_patients) <- "nCells"


mat_percent <- mat

for (col in 1:ncol(mat)) {
  col_name <- colnames(mat)[col]
  mat_percent[[col_name]] <- (mat[[col_name]] / nCell_patients[[col_name]]) * 100
}

# write.table(mat_percent, "/home/Data_Drive_8TB_3/ch2860/ColonCancer/Basic_seurat/Proportion_heatmap/Tseg_EpiCyc_percent.txt", quote = F , sep = "\t")

nCell_celltype <- data.frame(Celltype_nCell = rowSums(mat) )

head(nPatients_list)


mat_percent_for_temp <- mat_percent


# c("Epithelial","CD4T","CD8T","Myeloid","B","Plasma","NK","Fibroblasts","Endothelial","Mast","Unknown"),]

mat_percent_for_temp <- mat_percent_for_temp[c("CD4 Treg","CD4 Th17","CD4 Tfh","CD4 Tnaive","CD4 Tcm"),]
# rownames(mat_percent_for_temp) <- c("Epithelial","CD4T","CD8T","Myeloid","B","Plasma","NK","Fibroblasts","Endothelial","Mast","Unknown","Immune")


# Vec_TMB <- TMB_vector
Vec_MSS <- MSS_vector
Vec_APC <- ifelse(APC_vector == "APC mutation", "APC MT", APC_vector)
Vec_TP53 <- ifelse(TP53_mutation_vector == "TP53 mutation", "TP53 MT", TP53_mutation_vector)
Vec_RAS <- ifelse(RAS_mutation_vector == "RAS mutation", "RAS MT", RAS_mutation_vector)
Vec_RAS <- ifelse(Vec_RAS == "WT", "RAS WT", Vec_RAS)

# Split_info <- list(Vec_TMB, Vec_MSS, Vec_TP53, Vec_RAS, Vec_APC)  
# Split_Col <- c("TMB","MSS","TP53","RAS", "APC")
Split_info <- list( Vec_MSS, Vec_APC, Vec_TP53, Vec_RAS)  
Split_Col <- c("MSS", "APC","TP53","RAS")

tpval_result_frame <- data.frame(row.names = rownames(mat_percent_for_temp) ) #, MSS= c() ,TP53=c() ,  RAS= c()  )        


for (i in  1:length(Split_info))  {
  colnames(mat_percent_for_temp) <- Split_info[[i]]
  value_1 <- levels(as.factor(Split_info[[i]] ))[1]
  value_2 <- levels(as.factor(Split_info[[i]] ))[2]
  
  for (j in 1:length(rownames(mat_percent_for_temp)) ){
    vector_1 <- as.numeric(mat_percent_for_temp[colnames(mat_percent_for_temp) == value_1][j,])
    vector_2 <- as.numeric(mat_percent_for_temp[colnames(mat_percent_for_temp) == value_2][j,])
    
    # t.test(vector_1, vector_2)$p.value 
    tpval_result_frame[rownames(mat_percent_for_temp)[j], Split_Col[i]] <- wilcox.test(vector_1, vector_2)$p.value
    
  }
}

mat_percent_for_temp
tpval_result_frame
# MSS_split_table_pct

# ggplot 객체들을 저장하는 리스트 생성
Vln_list <- list()
n_col = nrow(mat_percent_for_temp)

# 예시 데이터와 분할 정보를 가정한 상태에서 반복문을 통해 Vln_list를 채움
for (i in  1:length(Split_info))  {
  colnames(mat_percent_for_temp) <- Split_info[[i]]
  
  for (j in 1:length(rownames(mat_percent_for_temp)) ){
    
    if (j == 1 ){
      
      Vln_frame <- data.frame(
        Condition = colnames(mat_percent_for_temp) ,
        Proportion = unlist(mat_percent_for_temp[j,], use.names = F)
      )
      Vln_frame$Condition <- factor(Vln_frame$Condition, levels = rev(names(table( Split_info[[i]] ))) )
      
      # 첫 번째 행인지 확인
      is_first_row <- ((i-1) * length(rownames(mat_percent_for_temp)) + j) <= n_col
      
      # 첫 번째 행에만 제목과 x축 레이블을 추가
      title_text <- if(is_first_row) sprintf("%s", rownames(mat_percent_for_temp[j,])) else ""
      x_text <- sprintf("P = %s", round(tpval_result_frame[rownames(mat_percent_for_temp[j,]), Split_Col[i] ],digits = 3) )
      
      x_text_style <- if (round(tpval_result_frame[rownames(mat_percent_for_temp[j,]), Split_Col[i] ],digits = 3) < 0.05) {
        element_text(size = 20, color = "red", face = "bold")
      } else {
        element_text(size = 20)
      }
      
      
      p <- ggplot(Vln_frame, aes(x = Condition, y = Proportion, fill = Condition)) +
        geom_violin(trim = F) +
        labs(
          title = title_text,
          x = x_text,
          # x = "Condition",
          
          y = "Proportion(%)"
        ) + 
        # ylim(0, 100) +
        theme_minimal() +  
        scale_fill_manual(values = c( "TP53 WT" = "skyblue",  "TP53 MT"= "salmon" ,  "RAS WT" = "skyblue",  "RAS MT"= "salmon" , 
                                      "APC WT" = "skyblue",  "APC MT"= "salmon" ,   "MSS"= "skyblue" , "MSI-H" = "salmon")) +    
        # coord_fixed(ratio = 0.02 )+ 
        
        theme(legend.position = "none",
              axis.text.x=element_blank(),
              plot.title = element_text(size = 20,hjust = 0.5),
              axis.title = element_text(size = 15),
              axis.title.x = x_text_style,
              axis.title.y =  element_text(size = 18),
              axis.text.y = element_text(size = 15)
        )+ 
        geom_jitter(shape=16, position=position_jitter(0.0001))+ 
        stat_summary(fun=median, geom="point", size=2.5, color="red")
      
      
      Vln_list[[((i-1) * length(rownames(mat_percent_for_temp)) + j)]] <- p
    }
    
    
    else  {
      
      Vln_frame <- data.frame(
        Condition = colnames(mat_percent_for_temp) ,
        Proportion = unlist(mat_percent_for_temp[j,], use.names = F)
      )
      Vln_frame$Condition <- factor(Vln_frame$Condition, levels = rev(names(table( Split_info[[i]] ))) )
      
      # 첫 번째 행인지 확인
      is_first_row <- ((i-1) * length(rownames(mat_percent_for_temp)) + j) <= n_col
      
      # 첫 번째 행에만 제목과 x축 레이블을 추가
      title_text <- if(is_first_row) sprintf("%s", rownames(mat_percent_for_temp[j,])) else ""
      x_text <- sprintf("P = %s", round(tpval_result_frame[rownames(mat_percent_for_temp[j,]), Split_Col[i] ],digits = 3) )
      
      x_text_style <- if (round(tpval_result_frame[rownames(mat_percent_for_temp[j,]), Split_Col[i] ],digits = 3) < 0.05) {
        element_text(size = 20, color = "red", face = "bold")
      } else {
        element_text(size = 20)
      }
      
      
      p <- ggplot(Vln_frame, aes(x = Condition, y = Proportion, fill = Condition)) +
        geom_violin(trim = FALSE) +
        labs(title = title_text, x = x_text, y = "") + 
        ylim(0, NA) +
        theme_minimal() +  
        scale_fill_manual(values = c( "TP53 WT" = "skyblue",  "TP53 MT"= "salmon" ,  "RAS WT" = "skyblue",  "RAS MT"= "salmon" , 
                                      "APC WT" = "skyblue",  "APC MT"= "salmon" ,   "MSS"= "skyblue" , "MSI-H" = "salmon")) +   
        # coord_fixed(ratio = 0.02 )+ 
        theme(legend.position = "none",
              axis.text.x = element_blank(),
              plot.title = element_text(size = 20, hjust = 0.5),
              axis.title = element_text(size = 15),
              axis.title.x = x_text_style,
              axis.title.y =  element_text(size = 0),
              axis.text.y = element_text(size = 15)
        ) + 
        geom_jitter(shape = 16, position = position_jitter(0.0001)) + 
        stat_summary(fun = median, geom = "point", size = 2.5, color = "red")
      
      Vln_list[[((i-1) * length(rownames(mat_percent_for_temp)) + j)]] <- p
    }
    
    
  }
}

grid.arrange(grobs = Vln_list, ncol = n_col) # 리스트의 각 ggplot 객체들을 그리드 형태로 표시

table(T_CRC_3_filtered$Annotation_T_Fig1B)

setwd("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/")
total_width <- n_col * 170 # 각 플롯의 가로 크기를 150픽셀로 가정
total_height <- ceiling(length(Vln_list) / n_col) * 190 # 각 플롯의 세로 크기를 160픽셀로 가정
#png(filename = "/home/Data_Drive_8TB_3/Coloncancer/SequencingData/Proportion_stat_test/Immune_lineage.png", width = total_width, height = total_height, units = "in", res = 300)

## Fig1D Add - CD4T ##########
cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure1Dadd_CD4T_propor_wilcox_uploadMutation.pdf" , width = 10, height = 10 )
grid.arrange(grobs = Vln_list, ncol = n_col) # 리스트의 각 ggplot 객체들을 그리드 형태로 표시
dev.off()




# Fig1D Addtional CD8T cell Proportion Analysis ########################
CellType_CRC <- CD8T_CRC

# table(total_CRC$Annotation_Fig1)


# CellType_CRC <- subset(total_CRC , Annotation_Fig1 != "Unknown")

# 데이터프레임 만들기
data_df <- data.frame(
  Patient = CellType_CRC$patients,     # 환자 정보
  CellType = CellType_CRC$Annotation_T_Fig1B  # 세포 타입 정보
)

mat_temp <- as.data.frame(table(data_df) )

mat <- as.data.frame(matrix(mat_temp$Freq, ncol = length(unique(mat_temp$Patient)), byrow = TRUE ,dimnames = list(colnames( table(data_df)  ),rownames(table(data_df) ) )  )  )
rm(data_df)
rm(mat_temp)



mat <- mat[, nPatients_list]
head(mat)
colSums(mat)

table(CellType_CRC$patients)


nCell_patients <- as.data.frame(t(as.data.frame(as.matrix(table(CellType_CRC$patients)))[1]) )
rownames(nCell_patients) <- "nCells"


mat_percent <- mat

for (col in 1:ncol(mat)) {
  col_name <- colnames(mat)[col]
  mat_percent[[col_name]] <- (mat[[col_name]] / nCell_patients[[col_name]]) * 100
}

# write.table(mat_percent, "/home/Data_Drive_8TB_3/ch2860/ColonCancer/Basic_seurat/Proportion_heatmap/Tseg_EpiCyc_percent.txt", quote = F , sep = "\t")

nCell_celltype <- data.frame(Celltype_nCell = rowSums(mat) )

head(nPatients_list)


mat_percent_for_temp <- mat_percent


# c("Epithelial","CD4T","CD8T","Myeloid","B","Plasma","NK","Fibroblasts","Endothelial","Mast","Unknown"),]

mat_percent_for_temp <- mat_percent_for_temp[c("CD8 Tex","CD8 MAIT","CD8 Tem","CD8 Temra","CD8 Trm"),]


# Vec_TMB <- TMB_vector
Vec_MSS <- MSS_vector
Vec_APC <- ifelse(APC_vector == "APC mutation", "APC MT", APC_vector)
Vec_TP53 <- ifelse(TP53_mutation_vector == "TP53 mutation", "TP53 MT", TP53_mutation_vector)
Vec_RAS <- ifelse(RAS_mutation_vector == "RAS mutation", "RAS MT", RAS_mutation_vector)
Vec_RAS <- ifelse(Vec_RAS == "WT", "RAS WT", Vec_RAS)

# Split_info <- list(Vec_TMB, Vec_MSS, Vec_TP53, Vec_RAS, Vec_APC)  
# Split_Col <- c("TMB","MSS","TP53","RAS", "APC")
Split_info <- list( Vec_MSS, Vec_APC, Vec_TP53, Vec_RAS)  
Split_Col <- c("MSS", "APC","TP53","RAS")

tpval_result_frame <- data.frame(row.names = rownames(mat_percent_for_temp) ) #, MSS= c() ,TP53=c() ,  RAS= c()  )        


for (i in  1:length(Split_info))  {
  colnames(mat_percent_for_temp) <- Split_info[[i]]
  value_1 <- levels(as.factor(Split_info[[i]] ))[1]
  value_2 <- levels(as.factor(Split_info[[i]] ))[2]
  
  for (j in 1:length(rownames(mat_percent_for_temp)) ){
    vector_1 <- as.numeric(mat_percent_for_temp[colnames(mat_percent_for_temp) == value_1][j,])
    vector_2 <- as.numeric(mat_percent_for_temp[colnames(mat_percent_for_temp) == value_2][j,])
    
    # t.test(vector_1, vector_2)$p.value 
    tpval_result_frame[rownames(mat_percent_for_temp)[j], Split_Col[i]] <- wilcox.test(vector_1, vector_2)$p.value
    
  }
}

mat_percent_for_temp
tpval_result_frame
# MSS_split_table_pct

# ggplot 객체들을 저장하는 리스트 생성
Vln_list <- list()
n_col = nrow(mat_percent_for_temp)

# 예시 데이터와 분할 정보를 가정한 상태에서 반복문을 통해 Vln_list를 채움
for (i in  1:length(Split_info))  {
  colnames(mat_percent_for_temp) <- Split_info[[i]]
  
  for (j in 1:length(rownames(mat_percent_for_temp)) ){
    
    if (j == 1 ){
      
      Vln_frame <- data.frame(
        Condition = colnames(mat_percent_for_temp) ,
        Proportion = unlist(mat_percent_for_temp[j,], use.names = F)
      )
      Vln_frame$Condition <- factor(Vln_frame$Condition, levels = rev(names(table( Split_info[[i]] ))) )
      
      # 첫 번째 행인지 확인
      is_first_row <- ((i-1) * length(rownames(mat_percent_for_temp)) + j) <= n_col
      
      # 첫 번째 행에만 제목과 x축 레이블을 추가
      title_text <- if(is_first_row) sprintf("%s", rownames(mat_percent_for_temp[j,])) else ""
      x_text <- sprintf("P = %s", round(tpval_result_frame[rownames(mat_percent_for_temp[j,]), Split_Col[i] ],digits = 3) )
      
      x_text_style <- if (round(tpval_result_frame[rownames(mat_percent_for_temp[j,]), Split_Col[i] ],digits = 3) < 0.05) {
        element_text(size = 20, color = "red", face = "bold")
      } else {
        element_text(size = 20)
      }
      
      
      p <- ggplot(Vln_frame, aes(x = Condition, y = Proportion, fill = Condition)) +
        geom_violin(trim = F) +
        labs(
          title = title_text,
          x = x_text,
          # x = "Condition",
          
          y = "Proportion(%)"
        ) + 
        # ylim(0, 100) +
        theme_minimal() +  
        scale_fill_manual(values = c( "TP53 WT" = "skyblue",  "TP53 MT"= "salmon" ,  "RAS WT" = "skyblue",  "RAS MT"= "salmon" , 
                                      "APC WT" = "skyblue",  "APC MT"= "salmon" ,   "MSS"= "skyblue" , "MSI-H" = "salmon")) +   
        # coord_fixed(ratio = 0.02 )+ 
        
        theme(legend.position = "none",
              axis.text.x=element_blank(),
              plot.title = element_text(size = 20,hjust = 0.5),
              axis.title = element_text(size = 15),
              axis.title.x = x_text_style,
              axis.title.y =  element_text(size = 18),
              axis.text.y = element_text(size = 15)
        )+ 
        geom_jitter(shape=16, position=position_jitter(0.0001))+ 
        stat_summary(fun=median, geom="point", size=2.5, color="red")
      
      
      Vln_list[[((i-1) * length(rownames(mat_percent_for_temp)) + j)]] <- p
    }
    
    
    else  {
      
      Vln_frame <- data.frame(
        Condition = colnames(mat_percent_for_temp) ,
        Proportion = unlist(mat_percent_for_temp[j,], use.names = F)
      )
      Vln_frame$Condition <- factor(Vln_frame$Condition, levels = rev(names(table( Split_info[[i]] ))) )
      
      # 첫 번째 행인지 확인
      is_first_row <- ((i-1) * length(rownames(mat_percent_for_temp)) + j) <= n_col
      
      # 첫 번째 행에만 제목과 x축 레이블을 추가
      title_text <- if(is_first_row) sprintf("%s", rownames(mat_percent_for_temp[j,])) else ""
      x_text <- sprintf("P = %s", round(tpval_result_frame[rownames(mat_percent_for_temp[j,]), Split_Col[i] ],digits = 3) )
      
      x_text_style <- if (round(tpval_result_frame[rownames(mat_percent_for_temp[j,]), Split_Col[i] ],digits = 3) < 0.05) {
        element_text(size = 20, color = "red", face = "bold")
      } else {
        element_text(size = 20)
      }
      
      
      p <- ggplot(Vln_frame, aes(x = Condition, y = Proportion, fill = Condition)) +
        geom_violin(trim = FALSE) +
        labs(title = title_text, x = x_text, y = "") + 
        ylim(0, NA) +
        theme_minimal() +  
        scale_fill_manual(values = c( "TP53 WT" = "skyblue",  "TP53 MT"= "salmon" ,  "RAS WT" = "skyblue",  "RAS MT"= "salmon" , 
                                      "APC WT" = "skyblue",  "APC MT"= "salmon" ,   "MSS"= "skyblue" , "MSI-H" = "salmon")) +  
        
        # coord_fixed(ratio = 0.02 )+ 
        theme(legend.position = "none",
              axis.text.x = element_blank(),
              plot.title = element_text(size = 20, hjust = 0.5),
              axis.title = element_text(size = 15),
              axis.title.x = x_text_style,
              axis.title.y =  element_text(size = 0),
              axis.text.y = element_text(size = 15)
        ) + 
        geom_jitter(shape = 16, position = position_jitter(0.0001)) + 
        stat_summary(fun = median, geom = "point", size = 2.5, color = "red")
      
      Vln_list[[((i-1) * length(rownames(mat_percent_for_temp)) + j)]] <- p
    }
    
    
  }
}


grid.arrange(grobs = Vln_list, ncol = n_col) # 리스트의 각 ggplot 객체들을 그리드 형태로 표시

table(T_CRC_3_filtered$Annotation_T_Fig1B)

setwd("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/")
total_width <- n_col * 170 # 각 플롯의 가로 크기를 150픽셀로 가정
total_height <- ceiling(length(Vln_list) / n_col) * 190 # 각 플롯의 세로 크기를 160픽셀로 가정
#png(filename = "/home/Data_Drive_8TB_3/Coloncancer/SequencingData/Proportion_stat_test/Immune_lineage.png", width = total_width, height = total_height, units = "in", res = 300)

## Fig1D Add - CD8T ##########
cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure1Dadd_CD8T_propor_wilcox_uploadMutation.pdf" , width = 10, height = 10 )
grid.arrange(grobs = Vln_list, ncol = n_col) # 리스트의 각 ggplot 객체들을 그리드 형태로 표시
dev.off()













# monocyte, T cell  ####################################
Myeloid_CRC <- readRDS("/home/Data_Drive_8TB_2/ColonCancer/human1-31/Myeloid/Step5_H1-31_Myeloid_Seurat_v4.2.2_annotation_Final_2.rds")
table(Myeloid_CRC$annotation.mye)


clusters <- Myeloid_CRC$patients
cluster_names <- rep("CEBPB_GMM_Low", length(clusters))
cluster_names[clusters %in% CEBPBH_indices_GMM_Human] <- "CEBPB_GMM_High"
Myeloid_CRC <- AddMetaData(Myeloid_CRC, metadata = as.factor(cluster_names), col.name = "CEBPB_GMM")










##  Myeloid in  total ##############

clusters <- total_CRC$patients
cluster_names <- rep("CEBPB_GMM_Low", length(clusters))
cluster_names[clusters %in% CEBPBH_indices_GMM_Human] <- "CEBPB_GMM_High"
total_CRC <- AddMetaData(total_CRC, metadata = as.factor(cluster_names), col.name = "CEBPB_GMM")


as.character(total_CRC$Annotation_Fig1 )-> total_CRC$Annotation_total_lv4

as.data.frame(total_CRC$Annotation_total_lv4) -> tmp_meta
colnames(tmp_meta) <- "Celltype"

# T_CRC$Annotation_T_lv1 

as.data.frame(Myeloid_CRC$annotation.mye) -> tmp_meta_mye
colnames(tmp_meta_mye) <- "Celltype"

tmp_meta[rownames(tmp_meta_mye), "Celltype"] <- as.character(tmp_meta_mye[, 1])

total_CRC$Annotation_total_Myerough <- tmp_meta
table(total_CRC$Annotation_total_Myerough)


cluster_counts <- table(total_CRC$Annotation_total_Myerough)
# 먼저 "unassigned-CD74+" 와 "unassigned-KRT8+" 의 값을 추출합니다.
unassigned_values <- cluster_counts[c("unassigned-CD74+", "unassigned-KRT8+", "Unknown")]

# 그 다음 unassigned_values 를 제외한 나머지 값들을 정렬합니다.
sorted_clusters <- c(names(sort(cluster_counts[!names(cluster_counts) %in% c("unassigned-CD74+", "unassigned-KRT8+","Unknown")], decreasing = TRUE)), "unassigned-CD74+", "unassigned-KRT8+")# Seurat 객체의 factor 레벨을 재설정


total_CRC$Annotation_total_Myerough <- factor(total_CRC$Annotation_total_Myerough, levels = rev(sorted_clusters))
Idents(total_CRC) <- total_CRC$Annotation_total_Myerough

mT_annoted_facotr <- sorted_clusters


total_CRC_TP53_WT <- subset( total_CRC , TP53_Mut == "TP53WT" )
total_CRC_TP53_MT <- subset( total_CRC , TP53_Mut == "TP53mutation"  )


round(table(total_CRC_TP53_WT$Annotation_total_Myerough)[ rev(mT_annoted_facotr) ] / sum(table(total_CRC_TP53_WT$Annotation_total_Myerough) ) * 100 , 2) -> H_Mono_TP53WT_pct

round(table(total_CRC_TP53_MT$Annotation_total_Myerough)[ rev(mT_annoted_facotr) ] / sum(table(total_CRC_TP53_MT$Annotation_total_Myerough) ) * 100 , 2) -> H_Mono_TP53MT_pct


sorted_clusters

library(ggplot2)
# 예시 데이터 생성
fig3b2_category <-  mT_annoted_facotr 
fig3b2_wt <- as.numeric(  rev( H_Mono_TP53WT_pct))   # 첫 번째 데이터 값
fig3b2_tp53 <- as.numeric( rev(H_Mono_TP53MT_pct))   # 두 번째 데이터 값

# 데이터를 데이터프레임으로 변환
df <- data.frame(fig3b2_category, fig3b2_wt, fig3b2_tp53)

# 데이터를 long format으로 변환
df_long <- df %>%
  pivot_longer(cols = c(fig3b2_wt, fig3b2_tp53), names_to = "variable", values_to = "value") %>%
  mutate(variable = factor(variable, levels = c("fig3b2_wt", "fig3b2_tp53")))

# fig3b_category를 factor로 변환하여 순서를 정의
df_long$fig3b2_category <- factor(df_long$fig3b2_category, levels = fig3b2_category)

# 그래프 그리기
ggplot(df_long, aes(x = variable, y = value, fill = fig3b2_category)) +
  geom_bar(stat = "identity", position = "stack") +
  
  # scale_fill_manual(name = "Cell Type",
  #                   values = rev(c("#FFFF55", "#DDAA11", 
  #                                  "#BB3191", "#CCFF55", 
  #                                  "#6161FF", "#1111AA", "#A131FF",
  #                                  "#AFAF55", "#FF5151",  "#BB3131" ))
  #                   
  # ) +
  theme_classic() +
  labs(title = "Cell Composition by TP53 Mutation", fill = "Cell Type", 
       x = "",  # x 축 값 변경
       y = "Proportion (%)") +
  scale_x_discrete(labels = c("fig3b2_wt" = "Trp53 WT", "fig3b2_tp53" = "Trp53 MT")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 8), 
        axis.title.y = element_text(hjust = 0.5, size = 11), 
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(labels = scales::number_format(scale = 1, suffix = "%"))





CellType_CRC <-  total_CRC


# 데이터프레임 만들기
data_df <- data.frame(
  Patient = CellType_CRC$patients,     # 환자 정보
  CellType = CellType_CRC$Annotation_total_Myerough  # 세포 타입 정보
)

mat_temp <- as.data.frame(table(data_df) )

mat <- as.data.frame(matrix(mat_temp$Freq, ncol = length(unique(mat_temp$Patient)), byrow = TRUE ,dimnames = list(colnames( table(data_df)  ),rownames(table(data_df) ) )  )  )
rm(data_df)
rm(mat_temp)

mat <- mat[, nPatients_list]


nCell_patients <- as.data.frame(t(as.data.frame(as.matrix(table(CellType_CRC$patients)))[1]) )
rownames(nCell_patients) <- "nCells"

mat_percent <- mat

for (col in 1:ncol(mat)) {
  col_name <- colnames(mat)[col]
  mat_percent[[col_name]] <- (mat[[col_name]] / nCell_patients[[col_name]]) * 100
}

mat_percent_Myeloid <- mat_percent


TP53_mutation_vector -> TP53_mutation_vector_tmp
TP53_mutation_vector_tmp[TP53_mutation_vector_tmp == "TP53 mutation" ] <- "TP53 MT"
TP53_mutation_vector_tmp[TP53_mutation_vector_tmp == "TP53WT" ] <- "TP53 WT"


data <- data.frame(
  value = as.numeric(mat_percent_Myeloid["monocytes", ]),
  group = TP53_mutation_vector_tmp
)

# data$value[data$group == "TP53 WT"]

round(wilcox.test(data$value[data$group == "TP53 WT"] , data$value[data$group == "TP53 MT"])$p.value, 30) -> custom_p_value


ggviolin(data, x = "group", y = "value",
         color = "group", palette = my_colors, title = "", 
         xlab ="" , ylab = "Monocyte Proportion(Total)", ylim = c(0, round_up_to_tens(max(data$value))+10 )) + 
  geom_jitter(alpha = 0.2, width = 0.02 ,aes( color=group)) + 
  theme(legend.position="none",
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 10), 
        axis.title.y = element_text(hjust = 0.5, size = 11), 
        plot.title = element_text(hjust = 0.5)) +
  # theme_minimal( )+  theme(legend.position = "none") 
  stat_compare_means(label =  "p.format", label.x = 1.5 , label.y = round_up_to_tens(max(data$value))+ 5 ,  comparisons = list( c("TP53 WT","TP53 MT")  ) , size = 5 , show.legend = F  , tip.length = 0.1) # +  annotate("text", x = 1.5, y = round_up_to_tens(max(data$value)) +33 , label = custom_p_value, size = 5)





Myeloid_CEBPB_low <- subset( Myeloid_CRC , CEBPB_GMM == "CEBPB_GMM_Low"  )
Myeloid_CEBPB_high <- subset( Myeloid_CRC , CEBPB_GMM == "CEBPB_GMM_High" )


round(table(Myeloid_CEBPB_low$annotation.mye)[ rev(mT_annoted_facotr) ] / sum(table(Myeloid_CEBPB_low$annotation.mye) ) * 100 , 2) -> H_Mono_CEBPBlow_pct

round(table(Myeloid_CEBPB_high$annotation.mye)[ rev(mT_annoted_facotr) ] / sum(table(Myeloid_CEBPB_high$annotation.mye) ) * 100 , 2) -> H_Mono_CEBPBhigh_pct


sorted_clusters

library(ggplot2)
# 예시 데이터 생성
fig3b2_category <-  mT_annoted_facotr 
fig3b2_wt <- as.numeric(  rev( H_Mono_CEBPBlow_pct))   # 첫 번째 데이터 값
fig3b2_tp53 <- as.numeric( rev(H_Mono_CEBPBhigh_pct))   # 두 번째 데이터 값

# 데이터를 데이터프레임으로 변환
df <- data.frame(fig3b2_category, fig3b2_wt, fig3b2_tp53)

# 데이터를 long format으로 변환
df_long <- df %>%
  pivot_longer(cols = c(fig3b2_wt, fig3b2_tp53), names_to = "variable", values_to = "value") %>%
  mutate(variable = factor(variable, levels = c("fig3b2_wt", "fig3b2_tp53")))

# fig3b_category를 factor로 변환하여 순서를 정의
df_long$fig3b2_category <- factor(df_long$fig3b2_category, levels = fig3b2_category)

# 그래프 그리기
ggplot(df_long, aes(x = variable, y = value, fill = fig3b2_category)) +
  geom_bar(stat = "identity", position = "stack") +
  # scale_fill_manual(name = "Cell Type",
  #                   values = rev(c("#FFFF55", "#DDAA11", 
  #                                  "#BB3191", "#CCFF55", 
  #                                  "#6161FF", "#1111AA", "#A131FF",
  #                                  "#AFAF55", "#FF5151",  "#BB3131" ))
  #                   
  # ) +
  theme_classic() +
  labs(title = "Cell Composition by CEBPB H/L", fill = "Cell Type",
       x = "",  # x 축 값 변경
       y = "Proportion (%)") +
  scale_x_discrete(labels = c("fig3b2_wt" = "Trp53 WT", "fig3b2_tp53" = "Trp53 MT")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 8), 
        axis.title.y = element_text(hjust = 0.5, size = 11), 
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(labels = scales::number_format(scale = 1, suffix = "%"))





CellType_CRC <-  total_CRC


# 데이터프레임 만들기
data_df <- data.frame(
  Patient = CellType_CRC$patients,     # 환자 정보
  CellType = CellType_CRC$Annotation_total_Myerough  # 세포 타입 정보
)

mat_temp <- as.data.frame(table(data_df) )

mat <- as.data.frame(matrix(mat_temp$Freq, ncol = length(unique(mat_temp$Patient)), byrow = TRUE ,dimnames = list(colnames( table(data_df)  ),rownames(table(data_df) ) )  )  )
rm(data_df)
rm(mat_temp)

mat <- mat[, nPatients_list]


nCell_patients <- as.data.frame(t(as.data.frame(as.matrix(table(CellType_CRC$patients)))[1]) )
rownames(nCell_patients) <- "nCells"

mat_percent <- mat

for (col in 1:ncol(mat)) {
  col_name <- colnames(mat)[col]
  mat_percent[[col_name]] <- (mat[[col_name]] / nCell_patients[[col_name]]) * 100
}

mat_percent_Myeloid <- mat_percent


CEBPB_Patients_index$CEBPB.GMM -> CEBPB_Vector
CEBPB_Vector[CEBPB_Vector == "Low" ] <- "CEBPB Low"
CEBPB_Vector[CEBPB_Vector == "High" ] <- "CEBPB High"



data <- data.frame(
  value = as.numeric(mat_percent_Myeloid["monocytes", ]),
  group = CEBPB_Vector
)

# data$value[data$group == "TP53 WT"]

round(wilcox.test(data$value[data$group == "CEBPB Low"] , data$value[data$group == "CEBPB High"])$p.value, 30) -> custom_p_value


ggviolin(data, x = "group", y = "value",
         color = "group", palette = c("CEBPB Low"= "skyblue", "CEBPB High" = "salmon"), title = "", 
         xlab ="" , ylab = "Monocyte Proportion(Total)", ylim = c(0, round_up_to_tens(max(data$value))+10 )) + 
  geom_jitter(alpha = 0.2, width = 0.02 ,aes( color=group)) + 
  theme(legend.position="none",
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 10), 
        axis.title.y = element_text(hjust = 0.5, size = 11), 
        plot.title = element_text(hjust = 0.5)) +
  # theme_minimal( )+  theme(legend.position = "none") 
  stat_compare_means(label =  "p.format", label.x = 1.5 , label.y = round_up_to_tens(max(data$value))+ 5 ,  comparisons = list( c("CEBPB Low","CEBPB High")  ) , size = 5 , show.legend = F  , tip.length = 0.1) # +  annotate("text", x = 1.5, y = round_up_to_tens(max(data$value)) +33 , label = custom_p_value, size = 5)






















  


## Myeloid ##############

cluster_counts <- table(Myeloid_CRC$annotation.mye)
# 먼저 "unassigned-CD74+" 와 "unassigned-KRT8+" 의 값을 추출합니다.
unassigned_values <- cluster_counts[c("unassigned-CD74+", "unassigned-KRT8+")]

# 그 다음 unassigned_values 를 제외한 나머지 값들을 정렬합니다.
sorted_clusters <- c(names(sort(cluster_counts[!names(cluster_counts) %in% c("unassigned-CD74+", "unassigned-KRT8+")], decreasing = TRUE)), "unassigned-CD74+", "unassigned-KRT8+")# Seurat 객체의 factor 레벨을 재설정


Myeloid_CRC$annotation.mye <- factor(Myeloid_CRC$annotation.mye, levels = rev(sorted_clusters))
Idents(Myeloid_CRC) <- Myeloid_CRC$annotation.mye

mT_annoted_facotr <- sorted_clusters


Myeloid_TP53_WT <- subset( Myeloid_CRC , TP53_Mut == "TP53WT" )
Myeloid_TP53_MT <- subset( Myeloid_CRC , TP53_Mut == "TP53mutation"  )


round(table(Myeloid_TP53_WT$annotation.mye)[ rev(mT_annoted_facotr) ] / sum(table(Myeloid_TP53_WT$annotation.mye) ) * 100 , 2) -> H_Mono_TP53WT_pct

round(table(Myeloid_TP53_MT$annotation.mye)[ rev(mT_annoted_facotr) ] / sum(table(Myeloid_TP53_MT$annotation.mye) ) * 100 , 2) -> H_Mono_TP53MT_pct


sorted_clusters

library(ggplot2)
# 예시 데이터 생성
fig3b2_category <-  mT_annoted_facotr 
fig3b2_wt <- as.numeric(  rev( H_Mono_TP53WT_pct))   # 첫 번째 데이터 값
fig3b2_tp53 <- as.numeric( rev(H_Mono_TP53MT_pct))   # 두 번째 데이터 값

# 데이터를 데이터프레임으로 변환
df <- data.frame(fig3b2_category, fig3b2_wt, fig3b2_tp53)

# 데이터를 long format으로 변환
df_long <- df %>%
  pivot_longer(cols = c(fig3b2_wt, fig3b2_tp53), names_to = "variable", values_to = "value") %>%
  mutate(variable = factor(variable, levels = c("fig3b2_wt", "fig3b2_tp53")))

# fig3b_category를 factor로 변환하여 순서를 정의
df_long$fig3b2_category <- factor(df_long$fig3b2_category, levels = fig3b2_category)

# 그래프 그리기
ggplot(df_long, aes(x = variable, y = value, fill = fig3b2_category)) +
  geom_bar(stat = "identity", position = "stack") +
  
  # scale_fill_manual(name = "Cell Type",
  #                   values = rev(c("#FFFF55", "#DDAA11", 
  #                                  "#BB3191", "#CCFF55", 
  #                                  "#6161FF", "#1111AA", "#A131FF",
  #                                  "#AFAF55", "#FF5151",  "#BB3131" ))
  #                   
  # ) +
  theme_classic() +
  labs(title = "Cell Composition by TP53 Mutation", fill = "Cell Type", 
       x = "",  # x 축 값 변경
       y = "Proportion (%)") +
  scale_x_discrete(labels = c("fig3b2_wt" = "TP53 WT", "fig3b2_tp53" = "TP53 MT")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 8), 
        axis.title.y = element_text(hjust = 0.5, size = 11), 
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(labels = scales::number_format(scale = 1, suffix = "%"))





CellType_CRC <-  Myeloid_CRC


# 데이터프레임 만들기
data_df <- data.frame(
  Patient = CellType_CRC$patients,     # 환자 정보
  CellType = CellType_CRC$annotation.mye  # 세포 타입 정보
)

mat_temp <- as.data.frame(table(data_df) )

mat <- as.data.frame(matrix(mat_temp$Freq, ncol = length(unique(mat_temp$Patient)), byrow = TRUE ,dimnames = list(colnames( table(data_df)  ),rownames(table(data_df) ) )  )  )
rm(data_df)
rm(mat_temp)

mat <- mat[, nPatients_list]


nCell_patients <- as.data.frame(t(as.data.frame(as.matrix(table(CellType_CRC$patients)))[1]) )
rownames(nCell_patients) <- "nCells"

mat_percent <- mat

for (col in 1:ncol(mat)) {
  col_name <- colnames(mat)[col]
  mat_percent[[col_name]] <- (mat[[col_name]] / nCell_patients[[col_name]]) * 100
}

mat_percent_Myeloid <- mat_percent


TP53_mutation_vector -> TP53_mutation_vector_tmp
TP53_mutation_vector_tmp[TP53_mutation_vector_tmp == "TP53 mutation" ] <- "TP53 MT"
TP53_mutation_vector_tmp[TP53_mutation_vector_tmp == "TP53WT" ] <- "TP53 WT"


data <- data.frame(
  value = as.numeric(mat_percent_Myeloid["monocytes", ]),
  group = TP53_mutation_vector_tmp
)

# data$value[data$group == "TP53 WT"]

round(wilcox.test(data$value[data$group == "TP53 WT"] , data$value[data$group == "TP53 MT"])$p.value, 30) -> custom_p_value


ggviolin(data, x = "group", y = "value",
         color = "group", palette = my_colors, title = "", 
         xlab ="" , ylab = "Monocyte Proportion(Myeloid)", ylim = c(0, round_up_to_tens(max(data$value))+40 )) + 
  geom_jitter(alpha = 0.2, width = 0.02 ,aes( color=group)) + 
  theme(legend.position="none",
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 10), 
        axis.title.y = element_text(hjust = 0.5, size = 11), 
        plot.title = element_text(hjust = 0.5)) +
  # theme_minimal( )+  theme(legend.position = "none") 
  stat_compare_means(label =  "p.format", label.x = 1.5 , label.y = round_up_to_tens(max(data$value))+ 35 ,  comparisons = list( c("TP53 WT","TP53 MT")  ) , size = 5 , show.legend = F  , tip.length = 0.1) # +  annotate("text", x = 1.5, y = round_up_to_tens(max(data$value)) +33 , label = custom_p_value, size = 5)





Myeloid_CEBPB_low <- subset( Myeloid_CRC , CEBPB_GMM == "CEBPB_GMM_Low"  )
Myeloid_CEBPB_high <- subset( Myeloid_CRC , CEBPB_GMM == "CEBPB_GMM_High" )


round(table(Myeloid_CEBPB_low$annotation.mye)[ rev(mT_annoted_facotr) ] / sum(table(Myeloid_CEBPB_low$annotation.mye) ) * 100 , 2) -> H_Mono_CEBPBlow_pct

round(table(Myeloid_CEBPB_high$annotation.mye)[ rev(mT_annoted_facotr) ] / sum(table(Myeloid_CEBPB_high$annotation.mye) ) * 100 , 2) -> H_Mono_CEBPBhigh_pct


sorted_clusters

library(ggplot2)
# 예시 데이터 생성
fig3b2_category <-  mT_annoted_facotr 
fig3b2_wt <- as.numeric(  rev( H_Mono_CEBPBlow_pct))   # 첫 번째 데이터 값
fig3b2_tp53 <- as.numeric( rev(H_Mono_CEBPBhigh_pct))   # 두 번째 데이터 값

# 데이터를 데이터프레임으로 변환
df <- data.frame(fig3b2_category, fig3b2_wt, fig3b2_tp53)

# 데이터를 long format으로 변환
df_long <- df %>%
  pivot_longer(cols = c(fig3b2_wt, fig3b2_tp53), names_to = "variable", values_to = "value") %>%
  mutate(variable = factor(variable, levels = c("fig3b2_wt", "fig3b2_tp53")))

# fig3b_category를 factor로 변환하여 순서를 정의
df_long$fig3b2_category <- factor(df_long$fig3b2_category, levels = fig3b2_category)

# 그래프 그리기
ggplot(df_long, aes(x = variable, y = value, fill = fig3b2_category)) +
  geom_bar(stat = "identity", position = "stack") +
  # scale_fill_manual(name = "Cell Type",
  #                   values = rev(c("#FFFF55", "#DDAA11", 
  #                                  "#BB3191", "#CCFF55", 
  #                                  "#6161FF", "#1111AA", "#A131FF",
  #                                  "#AFAF55", "#FF5151",  "#BB3131" ))
  #                   
  # ) +
  theme_classic() +
  labs(title = "Cell Composition by CEBPB H/L", fill = "Cell Type",
       x = "",  # x 축 값 변경
       y = "Proportion (%)") +
  scale_x_discrete(labels = c("fig3b2_wt" = "CEBPB Low", "fig3b2_tp53" = "CEBPB High")) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 8), 
        axis.title.y = element_text(hjust = 0.5, size = 11), 
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(labels = scales::number_format(scale = 1, suffix = "%"))





CellType_CRC <-  Myeloid_CRC


# 데이터프레임 만들기
data_df <- data.frame(
  Patient = CellType_CRC$patients,     # 환자 정보
  CellType = CellType_CRC$annotation.mye  # 세포 타입 정보
)

mat_temp <- as.data.frame(table(data_df) )

mat <- as.data.frame(matrix(mat_temp$Freq, ncol = length(unique(mat_temp$Patient)), byrow = TRUE ,dimnames = list(colnames( table(data_df)  ),rownames(table(data_df) ) )  )  )
rm(data_df)
rm(mat_temp)

mat <- mat[, nPatients_list]


nCell_patients <- as.data.frame(t(as.data.frame(as.matrix(table(CellType_CRC$patients)))[1]) )
rownames(nCell_patients) <- "nCells"

mat_percent <- mat

for (col in 1:ncol(mat)) {
  col_name <- colnames(mat)[col]
  mat_percent[[col_name]] <- (mat[[col_name]] / nCell_patients[[col_name]]) * 100
}

mat_percent_Myeloid <- mat_percent


CEBPB_Patients_index$CEBPB.GMM -> CEBPB_Vector
CEBPB_Vector[CEBPB_Vector == "Low" ] <- "CEBPB Low"
CEBPB_Vector[CEBPB_Vector == "High" ] <- "CEBPB High"



data <- data.frame(
  value = as.numeric(mat_percent_Myeloid["monocytes", ]),
  group = CEBPB_Vector
)

# data$value[data$group == "TP53 WT"]

round(wilcox.test(data$value[data$group == "CEBPB Low"] , data$value[data$group == "CEBPB High"])$p.value, 30) -> custom_p_value


ggviolin(data, x = "group", y = "value",
         color = "group", palette = c("CEBPB Low"= "skyblue", "CEBPB High" = "salmon"), title = "", 
         xlab ="" , ylab = "Monocyte Proportion(Myeloid)", ylim = c(0, round_up_to_tens(max(data$value))+40 )) + 
  geom_jitter(alpha = 0.2, width = 0.02 ,aes( color=group)) + 
  theme(legend.position="none",
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 10), 
        axis.title.y = element_text(hjust = 0.5, size = 11), 
        plot.title = element_text(hjust = 0.5)) +
  # theme_minimal( )+  theme(legend.position = "none") 
  stat_compare_means(label =  "p.format", label.x = 1.5 , label.y = round_up_to_tens(max(data$value))+ 35 ,  comparisons = list( c("CEBPB Low","CEBPB High")  ) , size = 5 , show.legend = F  , tip.length = 0.1) # +  annotate("text", x = 1.5, y = round_up_to_tens(max(data$value)) +33 , label = custom_p_value, size = 5)





























Myeloid_monocyte$TP53_Fig3_sampletype <- Myeloid_monocyte$TP53_Mut 

Myeloid_monocyte$TP53_Mut[Myeloid_monocyte$TP53_Fig3_sampletype == "TP53mutation" ] <- "TP53 MT"
Myeloid_monocyte$TP53_Mut[Myeloid_monocyte$TP53_Fig3_sampletype == "TP53WT" ] <- "TP53 WT"
table(Myeloid_monocyte$TP53_Mut)


data <- data.frame(
  value = Myeloid_monocyte@assays$RNA$data["CEBPB",] ,
  group = Myeloid_monocyte$TP53_Fig3_sampletype  
)

# data$value[data$group == "TP53 WT"]

round(wilcox.test(data$value[data$group == "TP53 WT"] , data$value[data$group == "TP53 MT"])$p.value, 30) -> custom_p_value


ggviolin(data, x = "group", y = "value",
         color = "group", palette = my_colors, title = "", 
         xlab ="" , ylab = "Normalized CEBPB(Epi) Expr", ylim = c(0, round_up_to_tens(max(data$value))-2 ))  +  
  theme(legend.position="none",
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
        axis.text.y = element_text(angle = 0, hjust = 0.5, size = 10), 
        axis.title.y = element_text(hjust = 0.5, size = 11), 
        plot.title = element_text(hjust = 0.5)) +
  # theme_minimal( )+  theme(legend.position = "none") 
  stat_compare_means(label =  "p.signif", label.x = 1.5 , label.y = round_up_to_tens(max(data$value))-4 ,  comparisons = list( c("TP53 WT","TP53 MT")  ) , size = 8 , show.legend = F  , tip.length = 0.1) + 
  annotate("text", x = 1.5, y = round_up_to_tens(max(data$value)) - 4, label = custom_p_value, size = 5)













































