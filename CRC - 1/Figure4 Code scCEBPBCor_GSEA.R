# Setting #################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(fgsea)
library(gridExtra)
set.seed(2860)
options(digits = 10)

# Treg Cells Promote Immunosuppression in Cancer Immune Escape
# Treg-Cell Differentiation
# Effector T-cell Inactivation in Cancer Immune Escape
# (Elsevier_Pathway_Collection CD4T, CD8T )
# 
# 
# CTLA4 Inhibitory Signaling
# (Reactome CD4T)
# 
# 
# Regulation Of Regulatory T Cell Differentiation
# Positive Regulation Of Programmed Cell Death
# (GO biological process, CD4T)
# 
# PD-L1 expression and PD-1 checkpoint pathway in cancer
# (KEGG, CD4T, CD8T)

# 
# Regulation_Of_Regulatory_T_Cell_Differentiation <- c("LAG3", "KAT2A", "DUSP10", "SOCS1", "BCL6", "KAT5", "CTLA4", "FOXO3", "CD46")
# 
# Positive_Regulation_Of_Programmed_Cell_Death <- c("FAF1", "HMGB1", "UBE2Z", "CASP8", "CTLA4", "EMILIN2", "DNM1L", "PHLDA1", "MAP3K5", "PHLDA3", "UNC13B", "DAPK2", "SIRT1", "DNM2", "LATS1", "IFNG", "ITGA6", "ARHGEF7", "SLC27A4", "HTRA1", "FOXO3", "STK4", "FOXO1", "SAV1", "BCLAF1", "BCL2L11", "E2F1", "PMAIP1", "ECT2", "MARK4", "GBP2", "ATG7", "MCL1", "GBP3", "GADD45A", "SIAH1", "DAB2IP", "TNFRSF10A", "CFLAR", "GADD45G", "REST", "BCL6", "STK17A", "STK17B", "JMY", "BCL2", "NF1", "CTNNB1", "BAX", "PDCD2")
# 
# PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer <- c("CD274", "PTEN", "CD3G", "PIK3R1", "HIF1A", "RASGRP1", "EGFR", "PPP3CA", "PPP3CB", "PPP3R1", "PPP3CC", "AKT3", "AKT1", "MAPK1", "PLCG1", "JAK2", "HRAS", "JAK1", "MAPK3", "MAP3K3", "BATF3", "MAP2K1", "MAP2K2", "STAT1", "CSNK2A2", "NFATC3", "NFATC2", "PTPN11", "MAPK14", "NFKB1", "MAPK11", "IFNG", "PIK3CA", "PRKCQ", "RAF1", "NFKBIE", "TLR4")
# 
# PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer_2 <- c("PTEN", "PIK3CB", "HIF1A", "RASGRP1", "PPP3CA", "PPP3CB", "PPP3R1", "PPP3CC", "AKT1", "MAPK1", "PLCG1", "JAK2", "HRAS", "JAK1", "MAPK3", "MAP3K3", "MAP2K1", "MAP2K2", "STAT1", "CSNK2A2", "STAT3", "PTPN11", "MAPK14", "NFKB1", "MAPK11", "NFKBIE")
# 
# CTLA4_Inhibitory_Signaling <- c("LYN", "PPP2CB", "PPP2R1B", "PPP2R5E", "AKT3", "PPP2R5A", "AKT1", "CTLA4", "PTPN11", "FYN")
# 
# Treg_Cells_Promote_Immunosuppression_in_Cancer_Immune_Escape <- c("MAP2K1", "MAP2K2", "SMAD3", "TGFB1", "ITPR1", "GZMB", "NFATC2", "AHR", "HIF1A", "TGFBR1", "TGFBR2", "IRF4", "AKT1", "CTLA4", "MAPK1", "GRB2", "PLCG1", "RAF1", "JAK3", "MAP3K14", "JAK1", "MAPK3")
# 
# Treg_Cell_Differentiation <- c("ITPR1", "CD3G", "MALT1", "MAPK9", "AKT1", "CTLA4", "MAPK1", "FYN", "PLCG1", "JAK3", "JAK1", "MAPK3", "MAP2K1", "MAP2K2", "SMAD3", "TGFB1", "MAPK14", "TGFBR1", "TGFBR2", "PTPRC", "IRF4", "PRKCQ", "GRB2", "RAF1", "MAP3K14", "CARD11")
# 
# Effector_T_cell_Inactivation_in_Cancer_Immune_Escape <- c("CD274", "CIITA", "LAG3", "NFATC2", "PTPN11", "TNFRSF10A", "PIK3R1", "CEACAM1", "CASP8", "IFNG", "BCL2", "AKT1", "CTLA4", "PRKCQ", "FYN", "GRB2", "PLCG1", "MAP3K5")




library(fgsea)

# gene_set <- list(
#   "Regulation_Of_Regulatory_T_Cell_Differentiation" = c("LAG3", "KAT2A", "DUSP10", "SOCS1", "BCL6", "KAT5", "CTLA4", "FOXO3", "CD46"),
#   
#   "Positive_Regulation_Of_Programmed_Cell_Death" = c("FAF1", "HMGB1", "UBE2Z", "CASP8", "CTLA4", "EMILIN2", "DNM1L", "PHLDA1", "MAP3K5", "PHLDA3", "UNC13B", "DAPK2", "SIRT1", "DNM2", "LATS1", "IFNG", "ITGA6", "ARHGEF7", "SLC27A4", "HTRA1", "FOXO3", "STK4", "FOXO1", "SAV1", "BCLAF1", "BCL2L11", "E2F1", "PMAIP1", "ECT2", "MARK4", "GBP2", "ATG7", "MCL1", "GBP3", "GADD45A", "SIAH1", "DAB2IP", "TNFRSF10A", "CFLAR", "GADD45G", "REST", "BCL6", "STK17A", "STK17B", "JMY", "BCL2", "NF1", "CTNNB1", "BAX", "PDCD2"),
#   # (GO:0043068)
#   
#   "PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer(CD4T)" = c("CD274", "PTEN", "CD3G", "PIK3R1", "HIF1A", "RASGRP1", "EGFR", "PPP3CA", "PPP3CB", "PPP3R1", "PPP3CC", "AKT3", "AKT1", "MAPK1", "PLCG1", "JAK2", "HRAS", "JAK1", "MAPK3", "MAP3K3", "BATF3", "MAP2K1", "MAP2K2", "STAT1", "CSNK2A2", "NFATC3", "NFATC2", "PTPN11", "MAPK14", "NFKB1", "MAPK11", "IFNG", "PIK3CA", "PRKCQ", "RAF1", "NFKBIE", "TLR4"),
#   
#   "PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer(CD8T)" = c("PTEN", "PIK3CB", "HIF1A", "RASGRP1", "PPP3CA", "PPP3CB", "PPP3R1", "PPP3CC", "AKT1", "MAPK1", "PLCG1", "JAK2", "HRAS", "JAK1", "MAPK3", "MAP3K3", "MAP2K1", "MAP2K2", "STAT1", "CSNK2A2", "STAT3", "PTPN11", "MAPK14", "NFKB1", "MAPK11", "NFKBIE"),
#   
#   "CTLA4_Inhibitory_Signaling" = c("LYN", "PPP2CB", "PPP2R1B", "PPP2R5E", "AKT3", "PPP2R5A", "AKT1", "CTLA4", "PTPN11", "FYN"),
#   
#   "Treg_Cells_Promote_Immunosuppression_in_Cancer_Immune_Escape" = c("MAP2K1", "MAP2K2", "SMAD3", "TGFB1", "ITPR1", "GZMB", "NFATC2", "AHR", "HIF1A", "TGFBR1", "TGFBR2", "IRF4", "AKT1", "CTLA4", "MAPK1", "GRB2", "PLCG1", "RAF1", "JAK3", "MAP3K14", "JAK1", "MAPK3"),
#   
#   "Treg_Cell_Differentiation" = c("ITPR1", "CD3G", "MALT1", "MAPK9", "AKT1", "CTLA4", "MAPK1", "FYN", "PLCG1", "JAK3", "JAK1", "MAPK3", "MAP2K1", "MAP2K2", "SMAD3", "TGFB1", "MAPK14", "TGFBR1", "TGFBR2", "PTPRC", "IRF4", "PRKCQ", "GRB2", "RAF1", "MAP3K14", "CARD11"),
#   
#   "Effector_T_cell_Inactivation_in_Cancer_Immune_Escape" = c("CD274", "CIITA", "LAG3", "NFATC2", "PTPN11", "TNFRSF10A", "PIK3R1", "CEACAM1", "CASP8", "IFNG", "BCL2", "AKT1", "CTLA4", "PRKCQ", "FYN", "GRB2", "PLCG1", "MAP3K5")
# 
# )



gene_set_fullgene <- list(
  "Regulation_Of_Regulatory_T_Cell_Differentiation" =  c("CR1", "LAG3", "FANCA", "LILRB2", "PSG4", "SOX12", "HLA-G", "LILRB4", "FOXO3", "VSIR", "FOXP3", 
                                                         "IL4I1", "IRF1", "DUSP10", "MDK", "SOCS1", "KAT2A", "BCL6", "KAT5", "TNFSF4", "CTLA4", "CD46"),
  
  "Positive_Regulation_Of_Programmed_Cell_Death" = c("UNC13B", "BAD", "CRADD", "PNMA5", "PNMA3", "SCRIB", "CCAR1", "SFRP4", "SFRP2", "APC", "NEURL1", 
                                                     "HIP1R", "BAX", "TP53", "SFRP1", "CCAR2", "SPDEF", "TFPT", "E2F1", "APBB1", "CYP1B1", "GBP1", "GBP3",
                                                     "GBP2", "GBP5", "SNCA", "ZNF268", "HOXA5", "KLF11", "SIAH1", "GPLD1", "WNT10B", "GADD45B", "GADD45A", 
                                                     "GADD45G", "FZD9", "LATS2", "LATS1", "MLLT11", "BIN1", "KCNMA1", "RNPS1", "ATF6", "ATF4", "FANK1", 
                                                     "BCL10", "DCUN1D3", "PRF1", "IFIT2", "SLC6A2", "SOX4", "RASSF2", "LGALS2", "LGALS1", "ARL6IP5", "TNFSF10", 
                                                     "ATG7", "BARD1", "DUSP6", "JUN", "NLRP2B", "TNFRSF12A", "SPHK2", "LILRB1", "SIRT1", "PAWR", "RHOB", "TFAP4", 
                                                     "ATM", "EEF1E1", "FCER2", "ARHGEF7", "HRK", "PXT1", "CSRNP3", "TOP2A", "OMA1", "TMEM164", "BCLAF1", "CLEC7A", 
                                                     "C1QBP", "HTATIP2", "ANKRD1", "CTLA4", "FADD", "PLEKHN1", "MAP3K20", "RIPK1", "CTSD", "NEUROD1", "PRMT2", "SYK", 
                                                     "DKKL1", "LAPTM5", "RRP1B", "SEPTIN4", "ANO6", "DAB2IP", "CLIP3", "RYBP", "STK17A", "STK17B", "RNF122", "FOXA1",
                                                     "FAF1", "RARG", "CTNNB1", "HTT", "ING5", "ING4", "AIFM2", "PNMA2", "PNMA1", "ADAMTSL4", "TCTN3", "AIFM1", "SMPD1",
                                                     "RACK1", "BID", "PHLDA1", "TGFB3", "RPS6", "POU4F2", "OLFM1", "MELK", "SHQ1", "UBE2Z", "IAPP", "FAS", "PDCD5", 
                                                     "PDCD2", "SQSTM1", "PPID", "CLU", "TNF", "FOXO3", "FOXO1", "PTPA", "PYCARD", "BCL2L11", "FRZB", "SAP18", "SLIT2",
                                                     "PHLDA3", "MAP3K5", "UTP11", "BNIP3L", "IGFBP3", "TNFRSF10C", "TNFRSF10B", "TNFRSF10A", "MAP3K10", "MAP3K11", 
                                                     "SLC27A4", "HTRA4", "HTRA2", "FASLG", "HTRA1", "TGM2", "CALHM2", "FAM162A", "UBD", "TNFRSF8", "PMAIP1", "BMF", 
                                                     "MARK4", "MOAP1", "DNM1L", "STPG1", "RFPL1", "APAF1", "CFLAR", "NRG1", "FOXL2", "DNM2", "BMP4", "BMP2", "IL6", 
                                                     "JMY", "BCL6", "IRF5", "IFNG", "PANO1", "BCL2", "NF1", "ITGA6", "QRICH1", "DDX3X", "DIABLO", "BCL2A1", "TNFAIP8",
                                                     "LRRK2", "TRADD", "HMGB1", "NLRC4", "SAV1", "STK3", "KNG1", "STK4", "CASP8", "SCIN", "CASP9", "CASP6", "PIP5KL1", 
                                                     "GPER1", "CASP2", "CASP3", "EMILIN2", "ECT2", "EMILIN1", "BOK", "IP6K2", "ZNF622", "CTNNBL1", "MCL1", "PIDD1", 
                                                     "NTRK1", "DAPK2", "DAPK1", "DAPK3", "GZMA", "SOD1", "DNAJA1", "REST", "ERCC3", "MNDA", "TRIM35", "MTCH1", "HDAC6",
                                                     "ITGB1", "MTCH2", "USP27X", "PRR7", "DDX20", "ACVR1C", "HSPD1", "FLCN", "RBM5", "MAPK8", "RPS6KA2", "BNIP3", 
                                                     "IGF2BP1", "BAK1", "NOS2", "C3ORF38") ,  
  
  "PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer" = c("MAP3K3", "EGF", "JUN", "NFATC3", "PTPN11", "NFATC2", "NFATC1", "TICAM2", "TICAM1", "BATF", "MAPK13",
                                                               "EML4", "MAPK14", "MAPK11", "NFKBIB", "TLR2", "MAPK12", "NFKBIA", "CD4", "PIK3CB", "PIK3CA", "IFNG", 
                                                               "CD28", "TLR9", "CD247", "PDCD1", "NFKBIE", "PTPN6", "TLR4", "ALK", "LAT", "RAF1", "RELA", "PIK3R3", 
                                                               "CD3G", "PIK3R2", "CD3E", "EGFR", "CD3D", "PPP3CA", "NRAS", "PPP3CB", "PPP3CC", "AKT3", "AKT1", "AKT2",
                                                               "BATF2", "BATF3", "MAP2K2", "MAP2K3", "CSNK2A2", "IFNGR1", "STAT1", "CHUK", "CSNK2A1", "IFNGR2", "MAP2K1",
                                                               "STAT3", "CSNK2A3", "FOS", "TIRAP", "MTOR", "NFKB1", "MYD88", "LCK", "ZAP70", "PIK3R1", "RPS6KB2", 
                                                               "RPS6KB1", "CSNK2B", "TRAF6", "KRAS", "PRKCQ", "CD274", "PTEN", "PIK3CD", "HIF1A", "IKBKB", "PPP3R1", 
                                                               "PPP3R2", "RASGRP1", "MAPK3", "HRAS", "MAPK1", "PLCG1", "IKBKG", "JAK1", "MAP2K6", "JAK2"),
  
  "CTLA4_Inhibitory_Signaling" = c("YES1", "CD86", "LYN", "SRC", "CD80", "PPP2R5C", "PPP2R5B", "PTPN11", "PPP2R1A", "PPP2R5E", "PPP2R5D", "PPP2R5A", "PPP2CA", "PPP2CB", "LCK", "PPP2R1B", "AKT3", "AKT1", "AKT2", "CTLA4", "FYN"),
  
  "Treg_Cells_Promote_Immunosuppression_in_Cancer_Immune_Escape" = c("SMAD3", "TGFB1", "ENTPD1", "MAP2K2", "IL10", "JUN", "MAP2K1", "GZMB", "PDPK1", "FOS", "NFATC2", "TGFBR1", "FOXP3", "TGFBR2", "IL2", "NFKBIA", "ZAP70", "IRF4", 
                                                                     "IL2RA", "IL2RB", "GRB2", "MAP3K14", "LAT", "RAF1", "SOS1", "PRF1", "AHR", "IL2RG", "HIF1A", "ITPR1", "IKBKB", "NT5E", "LGALS1", "AKT1", "MAPK3", "PTK2B", 
                                                                     "CTLA4", "MAPK1", "PLCG1", "TNFSF10", "JAK3", "MAP3K8", "JAK1", "STAT5A"),
  
  "Treg_Cell_Differentiation" = c("TGFB1", "IL10", "TGFBR1", "FOXP3", "TGFBR2", "IL2", "MAPK14", "NFKBIA", "PTPRC", "CD4", "IRF4", "CD28", "CD247", "PTPN6", "LAT", "RAF1", "SOS1", "CD86", "TRB", "TRA", "CD80", "CD3G",
                                  "IL2RG", "CD3E", "CD3D", "ITPR1", "AKT1", "GRAP2", "CTLA4", "FCER1G", "MAP3K7", "MAP3K8", "ICOS", "SMAD3", "CARD11", "MAP2K2", "MAP2K1", "PDPK1", "LCK", "ZAP70", 
                                  "IL2RA", "IL2RB", "TRAF6", "GRB2", "PRKCQ", "CD40LG", "MAP3K14", "BCL10", "MALT1", "MAPK9", "MAPK3", "MAPK1", "PLCG1", "FYN", "JAK3", "JAK1", "MAP2K7", "STAT5A"),
  
  "Effector_T_cell_Inactivation_in_Cancer_Immune_Escape" = c("LAG3", "BAD", "PTPN11", "NFATC2", "PDCD1LG2", "FOXP3", "TNFRSF1A", "NFKBIA", "CEACAM1", "FAS", "IFNG", "BCL2", "CD28", "PDCD1", "CLEC4G", "PTPN6", 
                                                             "TNFRSF25", "CD86", "VTCN1", "CD80", "TRADD", "TNF", "CASP8", "CASP9", "CASP10", "AKT1", "CASP2", "CASP3", "CTLA4", "FADD", 
                                                             "RIPK1", "MAP3K7", "CIITA", "MAP3K5", "DAXX", "PDPK1", "TNFRSF10B", "TNFRSF10A", "TRAF2", "ESR1", "LCK", "ZAP70", "PIK3R1", "BIRC5", 
                                                             "GRB2", "PRKCQ", "SIVA1", "EBAG9", "CD276", "CD274", "FASLG", "IKBKB", "PTK2B", "BTLA", "PLCG1", "FYN", "PLCG2", "MAP2K7"),
  

  "T_Cell_Receptor_Signaling" = c("RHOA", "IL2", "MAPK14", "VAV3", "IL5", "IL4", "NFKBIA", "PTPRC", "CD4", "IFNG", "CD28", "LCP2", 
                                  "CD247", "PDCD1", "PTPN6", "LAT", "RAF1", "SOS1", "CD86", "TRB", "TRA", "CD80", "CD3G", "CD3E", 
                                  "CD3D", "ITPR1", "PAK1", "AKT1", "GRAP2", "CTLA4", "FCER1G", "MAP3K7", "MAP3K8", "ICOS", "CARD11", 
                                  "MAP2K2", "MAP2K1", "PDPK1", "LCK", "ZAP70", "TRAF6", "GRB2", "PRKCQ", "CD40LG", "MAP3K14", "BCL10", 
                                  "MALT1", "MAPK9", "MAPK3", "MAPK1", "PLCG1", "FYN", "NCK1", "MAP2K7"),
  
  "T_Cell_Receptor_signaling_pathway" =  c("BUB1B-PAK6", "MAPK13", "MAPK14", "MAPK11", "TEC", "MAPK12", "CD4", "MAPK10", 
                                           "PDCD1", "PTPN6", "RAF1", "RELA", "CSF2", "TNF", "CBLB", "NRAS", "PAK2", "PAK1", 
                                           "AKT3", "AKT1", "AKT2", "GRAP2", "MAP3K7", "PAK4", "MAP3K8", "PAK3", "PAK6", 
                                           "PAK5", "CARD11", "MAP2K2", "MAP2K1", "CD8B2", "NFKB1", "DLG1", "PIK3R1", "CD8B", 
                                           "CD8A", "GRB2", "PRKCQ", "MAP3K14", "BCL10", "PPP3R1", "PPP3R2", "RASGRP1", 
                                           "PLCG1", "MAP2K7", "RHOA", "IL10", "JUN", "NFATC3", "NFATC2", "NFATC1", "VAV1", 
                                           "VAV2", "IL2", "VAV3", "IL5", "NFKBIB", "IL4", "NFKBIA", "PTPRC", "PIK3CB", 
                                           "PIK3CA", "CDK4", "IFNG", "CD28", "LCP2", "CD247", "SOS2", "NFKBIE", "ITK", 
                                           "LAT", "SOS1", "GSK3B", "PIK3R3", "CD3G", "PIK3R2", "CD3E", "CD3D", "PPP3CA", 
                                           "PPP3CB", "PPP3CC", "CTLA4", "ICOS", "CHUK", "PDPK1", "FOS", "CDC42", "LCK", 
                                           "ZAP70", "KRAS", "CD40LG", "PIK3CD", "MALT1", "IKBKB", "MAPK9", "MAPK8", "MAPK3", 
                                           "HRAS", "MAPK1", "FYN", "IKBKG", "NCK1", "NCK2")
  
  
  )

SC_CEBPB_PCC_COR_ENRICHR_GENE <- readRDS("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/Figure4_Gene_list/SC_CEBPB_PCC_COR_ENRICHR_GENE.rds")
sc_cor_frame <- readRDS("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/Figure4_Gene_list/sc_cor_frame.rds")

SC_CEBPB_PCC_COR_ENRICHR_GENE$CD4T_pos[which(SC_CEBPB_PCC_COR_ENRICHR_GENE$CD4T_pos != 0)] -> CD4T_posgene
SC_CEBPB_PCC_COR_ENRICHR_GENE$CD8T_pos[which(SC_CEBPB_PCC_COR_ENRICHR_GENE$CD8T_pos != 0)] -> CD8T_posgene

SC_CEBPB_PCC_COR_ENRICHR_GENE$CD4T_neg[which(SC_CEBPB_PCC_COR_ENRICHR_GENE$CD4T_neg != 0)] -> CD4T_neggene
SC_CEBPB_PCC_COR_ENRICHR_GENE$CD8T_neg[which(SC_CEBPB_PCC_COR_ENRICHR_GENE$CD8T_neg != 0)] -> CD8T_neggene


# Gene_define ####################################
## Correlation Whole Gene in CD4T #################
sc_cor_frame[2] -> genelist_tmp
genelist_tmp$CD4T_Cor -> geneList
rownames(genelist_tmp) -> names(geneList)
sort(geneList, decreasing = T) -> geneList
length(geneList)


## Correlation Sig Gene in CD4T  #################
sc_cor_frame[ rownames(sc_cor_frame)[rownames(sc_cor_frame) %in% c(CD4T_posgene, CD4T_neggene ) ]  , ][2] -> genelist_tmp
genelist_tmp$CD4T_Cor -> geneList
rownames(genelist_tmp) -> names(geneList)
sort(geneList, decreasing = T) -> geneList
length(geneList)


## Correlation Pos Sig Gene in CD4T  #################
sc_cor_frame[ rownames(sc_cor_frame)[rownames(sc_cor_frame) %in% CD4T_posgene]  , ][2] -> genelist_tmp
genelist_tmp$CD4T_Cor -> geneList
rownames(genelist_tmp) -> names(geneList)
sort(geneList, decreasing = T) -> geneList
length(geneList)





## Correlation Whole Gene in CD8T #################
sc_cor_frame[3] -> genelist_tmp
genelist_tmp$CD8T_Cor -> geneList
rownames(genelist_tmp) -> names(geneList)
sort(geneList, decreasing = T) -> geneList
length(geneList)


## Correlation Sig Gene in CD8T  #################
sc_cor_frame[ rownames(sc_cor_frame)[rownames(sc_cor_frame) %in% c(CD8T_posgene, CD8T_neggene ) ]  , ][3] -> genelist_tmp
genelist_tmp$CD8T_Cor -> geneList
rownames(genelist_tmp) -> names(geneList)
sort(geneList, decreasing = T) -> geneList
length(geneList)


## Correlation Pos Sig Gene in CD8T  #################
sc_cor_frame[ rownames(sc_cor_frame)[rownames(sc_cor_frame) %in% CD8T_posgene]  , ][3] -> genelist_tmp
genelist_tmp$CD8T_Cor -> geneList
rownames(genelist_tmp) -> names(geneList)
sort(geneList, decreasing = T) -> geneList
length(geneList)





fgsea_results <- fgsea(pathways = gene_set_fullgene, stats = geneList , nperm = 1000 )

# 결과를 데이터프레임으로 변환합니다.
fgsea_results_df <- as.data.frame(fgsea_results)







plotEnrichment(gene_set_fullgene[["Regulation_Of_Regulatory_T_Cell_Differentiation"]], geneList) +
  labs(title = "GSEA Plot: Regulation Of Regulatory T Cell Differentiation")


plotEnrichment(gene_set_fullgene[["Positive_Regulation_Of_Programmed_Cell_Death"]], geneList) +
  labs(title = "GSEA Plot: Positive Regulation Of Programmed Cell Death")


plotEnrichment(gene_set_fullgene[["PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer(CD4T)"]], geneList) +
  labs(title = "GSEA Plot: PD L1 expression and PD 1 checkpoint pathway in cancer")


plotEnrichment(gene_set_fullgene[["PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer(CD4T)"]], geneList)  -> aaaa
aaaa$data$rank
geneList[940]

plotEnrichment(gene_set_fullgene[["PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer(CD8T)"]], geneList) +
  labs(title = "GSEA Plot: PD L1 expression and PD 1 checkpoint pathway in cancer(CD8T)")


plotEnrichment(gene_set_fullgene[["CTLA4_Inhibitory_Signaling"]], geneList) +
  labs(title = "GSEA Plot: CTLA4 Inhibitory Signaling")

plotEnrichment(gene_set_fullgene[["CTLA4_Inhibitory_Signaling"]], geneList) -> bbbb
bbbb$data$rank
geneList[2216]


plotEnrichment(gene_set_fullgene[["Treg_Cells_Promote_Immunosuppression_in_Cancer_Immune_Escape"]], geneList) +
  labs(title = "GSEA Plot: Treg Cells Promote Immunosuppression in Cancer Immune Escape")

plotEnrichment(gene_set_fullgene[["Treg_Cell_Differentiation"]], geneList) +
  labs(title = "GSEA Plot: Treg Cell Differentiation")

plotEnrichment(gene_set_fullgene[["Effector_T_cell_Inactivation_in_Cancer_Immune_Escape"]], geneList) +
  labs(title = "GSEA Plot: Effector T cell Inactivation in Cancer Immune Escape")







# For Example #########################
## Correlation Whole Gene in CD4T #################
sc_cor_frame[2] -> genelist_tmp
genelist_tmp$CD4T_Cor -> geneList
rownames(genelist_tmp) -> names(geneList)
sort(geneList, decreasing = T) -> geneList
length(geneList)



fgsea_results <- fgsea(pathways = gene_set_fullgene, stats = geneList  )

fgsea_results_CD4T <- fgsea_results

# 결과를 데이터프레임으로 변환합니다.
fgsea_results_df <- as.data.frame(fgsea_results)


ranked_list <- geneList


ranked_df <- data.frame(
  Gene = names(ranked_list),
  Correlation = ranked_list
)

# Bar plot 생성
bar_plot <- ggplot(ranked_df, aes(x = seq_along(Correlation), y = 1, fill = Correlation)) +
  geom_tile(height = 0.1) +  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(y = "")+
  theme(axis.title.x = element_blank(),
        plot.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6, colour = "white"),
        
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        plot.margin = unit(c(0, 1, 0, 1), "lines"))


# scatter_plot <- ggplot(ranked_df, aes(x = seq_along(Correlation), y = Correlation)) +
#   geom_point(size = 0.1) +
#   theme_minimal() +
#   labs(x = "Rank", y = "Correlation")+
#   theme(plot.margin = unit(c(0, 1, 1, 1), "lines"))


scatter_plot <- ggplot(ranked_df, aes(x = seq_along(Correlation), y = Correlation)) +
  geom_hline(yintercept = c(-0.5,0,  0.5), color = "#EEEEEE", linetype = "solid") +
  geom_vline(xintercept = c(0,10000,20000), color = "#EEEEEE", linetype = "solid") +
  
  geom_bar(stat = "identity", size = 0.1, color = "lightgray") +
  # geom_text(aes(x = 17000, y = 0, label = "Zero"), size = 3, color = "lightgray", family = "Arial") +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.5)) +
  scale_x_continuous(breaks = seq(0, 30000, by = 10000)) +
  labs(x = "Rank", y = "Correlation") +
  
  theme_void()+
  theme(axis.text.y = element_text(size = 8,  family = "arial", margin = margin(0,4,0,0 ) ),
        axis.text.x = element_text(size = 8,  family = "arial"),
        axis.title.x = element_text(size = 11, margin = margin(3,1,1,1), family = "arial" , vjust = 0.5),
        axis.title.y = element_text(size = 11, margin = margin(1,5,1,1), angle = 90, family = "arial" , hjust = 0.5 ),
        
        
        plot.margin = unit(c(0, 1, 1, 1), "lines"))






cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/CEBPBcorGSEAPlot/GOBP_CD4T_TregDiffere.pdf" , width = 6, height = 4)


# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[["Regulation_Of_Regulatory_T_Cell_Differentiation"]])

# Enrichment plot 생성
enrichment_plot <- plotEnrichment(gene_set_fullgene[["Regulation_Of_Regulatory_T_Cell_Differentiation"]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "GSEA Plot: Regulation Of Regulatory T Cell Differentiation" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title.x = element_blank(),   
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 8),
        
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) 



grid.arrange(enrichment_plot, bar_plot, scatter_plot, ncol = 1, heights = c(3, 0.5, 2))




dev.off()




cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/CEBPBcorGSEAPlot/GOBP_CD4T_ProgrammedCellDeath.pdf" , width = 6, height = 4)

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[["Positive_Regulation_Of_Programmed_Cell_Death"]])

# Enrichment plot 생성
enrichment_plot <- plotEnrichment(gene_set_fullgene[["Positive_Regulation_Of_Programmed_Cell_Death"]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "GSEA Plot: Positive Regulation Of Programmed Cell Death" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 8),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) 


grid.arrange(enrichment_plot, bar_plot, scatter_plot, ncol = 1, heights = c(3, 0.5, 2))

dev.off()


# plotEnrichment(gene_set_fullgene[["PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer(CD4T)"]], geneList) +
#   labs(title = "GSEA Plot: PD L1 expression and PD 1 checkpoint pathway in cancer")



cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/CEBPBcorGSEAPlot/Kegg_CD4T_PDL1.pdf" , width = 6, height = 4)

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[["PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer"]])

# Enrichment plot 생성
enrichment_plot <- plotEnrichment(gene_set_fullgene[["PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer"]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "GSEA Plot: PD L1 expression and PD 1 checkpoint pathway in cancer" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 8),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) 

grid.arrange(enrichment_plot, bar_plot, scatter_plot, ncol = 1, heights = c(3, 0.5, 2))


dev.off()





cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/CEBPBcorGSEAPlot/Reactome_CTLA4inhi.pdf" , width = 6, height = 4)

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[["CTLA4_Inhibitory_Signaling"]])

# Enrichment plot 생성
enrichment_plot <- plotEnrichment(gene_set_fullgene[["CTLA4_Inhibitory_Signaling"]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "GSEA Plot: CTLA4 Inhibitory Signaling" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 8),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) 



grid.arrange(enrichment_plot, bar_plot, scatter_plot, ncol = 1, heights = c(3, 0.5, 2))


dev.off()





cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/CEBPBcorGSEAPlot/Eles_TregImmuneEscape.pdf" , width = 6, height = 4)

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[["Treg_Cells_Promote_Immunosuppression_in_Cancer_Immune_Escape"]])

# Enrichment plot 생성
enrichment_plot <- plotEnrichment(gene_set_fullgene[["Treg_Cells_Promote_Immunosuppression_in_Cancer_Immune_Escape"]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "GSEA Plot: Treg Cells Promote Immunosuppression in Cancer Immune Escape" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 8),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) 



grid.arrange(enrichment_plot, bar_plot, scatter_plot, ncol = 1, heights = c(3, 0.5, 2))


dev.off()




cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/CEBPBcorGSEAPlot/Eles_Treg_Cell_Differentiation.pdf" , width = 6, height = 4)

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[["Treg_Cell_Differentiation"]])

# Enrichment plot 생성
enrichment_plot <- plotEnrichment(gene_set_fullgene[["Treg_Cell_Differentiation"]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "GSEA Plot: Treg Cell Differentiation" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 8),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) 



grid.arrange(enrichment_plot, bar_plot, scatter_plot, ncol = 1, heights = c(3, 0.5, 2))


dev.off()




cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/CEBPBcorGSEAPlot/Eles_Effector_CD4T_cell_Inactivation_in_Cancer_Immune_Escape.pdf" , width = 6, height = 4)

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[["Effector_T_cell_Inactivation_in_Cancer_Immune_Escape"]])

# Enrichment plot 생성
enrichment_plot <- plotEnrichment(gene_set_fullgene[["Effector_T_cell_Inactivation_in_Cancer_Immune_Escape"]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "GSEA Plot: Effector T cell Inactivation in Cancer Immune Escape" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 8),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) 


grid.arrange(enrichment_plot, bar_plot, scatter_plot, ncol = 1, heights = c(3, 0.5, 2))


dev.off()










## Correlation Whole Gene in CD8T #################
sc_cor_frame[3] -> genelist_tmp
genelist_tmp$CD8T_Cor -> geneList
rownames(genelist_tmp) -> names(geneList)
sort(geneList, decreasing = T) -> geneList
length(geneList)


fgsea_results <- fgsea(pathways = gene_set_fullgene, stats = geneList , nperm = 10000 )

fgsea_results_CD8T <- fgsea_results
fgsea_results_CD8T$padj

# 결과를 데이터프레임으로 변환합니다.
fgsea_results_df <- as.data.frame(fgsea_results_CD8T)


ranked_list <- geneList



ranked_df <- data.frame(
  Gene = names(ranked_list),
  Correlation = ranked_list
)

# Bar plot 생성
bar_plot <- ggplot(ranked_df, aes(x = seq_along(Correlation), y = 1, fill = Correlation)) +
  geom_tile(height = 0.1) +  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(y = "")+
  theme(axis.title.x = element_blank(),
        plot.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6, colour = "white"),
        
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        plot.margin = unit(c(0, 1, 0, 1), "lines"))



# scatter_plot <- ggplot(ranked_df, aes(x = seq_along(Correlation), y = Correlation)) +
#   geom_point(size = 0.1) +
#   theme_minimal() +
#   labs(x = "Rank", y = "Correlation")+
#   theme(plot.margin = unit(c(0, 1, 1, 1), "lines"))


scatter_plot <- ggplot(ranked_df, aes(x = seq_along(Correlation), y = Correlation)) +
  geom_hline(yintercept = c(-0.5,0,  0.5), color = "#EEEEEE", linetype = "solid") +
  geom_vline(xintercept = c(0,10000,20000), color = "#EEEEEE", linetype = "solid") +
  
  geom_bar(stat = "identity", size = 0.1, color = "lightgray") +
  # geom_text(aes(x = 17000, y = 0, label = "Zero"), size = 3, color = "lightgray", family = "Arial") +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.5)) +
  scale_x_continuous(breaks = seq(0, 30000, by = 10000)) +
  labs(x = "Rank", y = "Correlation") +
  
  theme_void()+
  theme(axis.text.y = element_text(size = 8,  family = "arial", margin = margin(0,4,0,0 ) ),
        axis.text.x = element_text(size = 8,  family = "arial"),
        axis.title.x = element_text(size = 11, margin = margin(3,1,1,1), family = "arial" , vjust = 0.5),
        axis.title.y = element_text(size = 11, margin = margin(1,5,1,1), angle = 90, family = "arial" , hjust = 0.5 ),
        
        
        plot.margin = unit(c(0, 1, 1, 1), "lines"))


scatter_plot_b <- ggplot(ranked_df, aes(x = seq_along(Correlation), y = Correlation)) +
  geom_hline(yintercept = c(-0.5,0,  0.5), color = "#EEEEEE", linetype = "solid") +
  geom_vline(xintercept = c(0,10000,20000), color = "#EEEEEE", linetype = "solid") +
  
  geom_bar(stat = "identity", size = 0.1, color = "lightgray") +
  # geom_text(aes(x = 17000, y = 0, label = "Zero"), size = 3, color = "lightgray", family = "Arial") +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.5)) +
  scale_x_continuous(breaks = seq(0, 30000, by = 10000)) +
  labs(x = "Rank", y = "Correlation") +
  
  theme_void()+
  theme(axis.text.y = element_text(size = 8,  family = "arial", margin = margin(0,4,0,0 ) ),
        axis.text.x = element_text(size = 8,  family = "arial"),
        axis.title.x = element_text(size = 11, margin = margin(3,1,1,1), family = "arial" , vjust = 0.5),
        axis.title.y = element_text(size = 11, margin = margin(1,5,1,1), angle = 90, family = "arial" , hjust = 0.5 , color = "white"),
        
        
        plot.margin = unit(c(0, 1, 1, 1), "lines"))



# plotEnrichment(gene_set[["PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer(CD8T)"]], geneList) +
#   labs(title = "GSEA Plot: PD L1 expression and PD 1 checkpoint pathway in cancer(CD8T)")


cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/CEBPBcorGSEAPlot/Kegg_CD8T_PDL1.pdf" , width = 5, height = 4)

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[["PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer"]])

# Enrichment plot 생성
enrichment_plot_cd8ta <- plotEnrichment(gene_set_fullgene[["PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer"]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "PD L1 expression and PD 1 checkpoint pathway in cancer" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 8),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) 

grid.arrange(enrichment_plot_a, bar_plot, scatter_plot, ncol = 1, heights = c(3, 0.5, 2))


dev.off()



cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/CEBPBcorGSEAPlot/Eles_Effector_CD8T_cell_Inactivation_in_Cancer_Immune_Escape.pdf" , width = 6, height = 4)

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[["Effector_T_cell_Inactivation_in_Cancer_Immune_Escape"]])

# Enrichment plot 생성
enrichment_plot_cd8tb <- plotEnrichment(gene_set_fullgene[["Effector_T_cell_Inactivation_in_Cancer_Immune_Escape"]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Effector T cell Inactivation in Cancer Immune Escape" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 8),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) 


grid.arrange(enrichment_plot_b, bar_plot, scatter_plot, ncol = 1, heights = c(3, 0.5, 2))


dev.off()




# 원하는 항목들만 새로운 리스트에 저장
gene_set_fullgene_tmp <- list(
  Treg_Cell_Differentiation = gene_set_fullgene$Treg_Cell_Differentiation,
  PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer = gene_set_fullgene$PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer,
  Treg_Cells_Promote_Immunosuppression_in_Cancer_Immune_Escape = gene_set_fullgene$Treg_Cells_Promote_Immunosuppression_in_Cancer_Immune_Escape,
  T_Cell_Receptor_Signaling = gene_set_fullgene$T_Cell_Receptor_Signaling,
  Effector_T_cell_Inactivation_in_Cancer_Immune_Escape = gene_set_fullgene$Effector_T_cell_Inactivation_in_Cancer_Immune_Escape,
  CTLA4_Inhibitory_Signaling = gene_set_fullgene$CTLA4_Inhibitory_Signaling,
  Regulation_Of_Regulatory_T_Cell_Differentiation = gene_set_fullgene$Regulation_Of_Regulatory_T_Cell_Differentiation,
  Positive_Regulation_Of_Programmed_Cell_Death = gene_set_fullgene$Positive_Regulation_Of_Programmed_Cell_Death
)

# Fig4G  ###################
## Correlation Whole Gene in CD4T #################
sc_cor_frame[2] -> genelist_tmp
genelist_tmp$CD4T_Cor -> geneList
rownames(genelist_tmp) -> names(geneList)
sort(geneList, decreasing = T) -> geneList
length(geneList)



# fgsea_results <- fgsea(pathways = gene_set_fullgene_tmp, stats = geneList, gseaParam = 1000 )
# fgsea_results <- fgsea(pathways = gene_set_fullgene_tmp, stats = geneList )

# fgsea_results_CD4T <- fgsea_results
# fgsea_results_CD4T

# saveRDS(as.data.frame(fgsea_results_CD4T), "/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/Figure4_Gene_list/fgsea_results_CD4T.rds")


# 결과를 데이터프레임으로 변환합니다.
# fgsea_results_df <- as.data.frame(fgsea_results)
# fgsea_results_df

fgsea_results_CD4T <- readRDS( "/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/Figure4_Gene_list/fgsea_results_CD4T.rds")

ranked_list <- geneList


ranked_df <- data.frame(
  Gene = names(ranked_list),
  Correlation = ranked_list
)





# Bar plot 생성
bar_plot <- ggplot(ranked_df, aes(x = seq_along(Correlation), y = 1, fill = Correlation)) +
  geom_tile(height = 0.1) +  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(y = "")+
  theme(axis.title.x = element_blank(),
        plot.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 9, colour = "white" , margin = margin(0,3,0,3)),
        
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        plot.margin = unit(c(0, 1, 0, 1), "lines"))


scatter_plot <- ggplot(ranked_df, aes(x = seq_along(Correlation), y = Correlation)) +
  geom_hline(yintercept = c(-0.5,0,  0.5), color = "#EEEEEE", linetype = "solid") +
  geom_vline(xintercept = c(0,10000,20000), color = "#EEEEEE", linetype = "solid") +
  
  geom_bar(stat = "identity", size = 0.1, color = "lightgray") +
  # geom_text(aes(x = 17000, y = 0, label = "Zero"), size = 3, color = "lightgray", family = "Arial") +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.5)) +
  scale_x_continuous(breaks = seq(0, 30000, by = 10000)) +
  labs(x = "Rank", y = "Correlation") +
  
  theme_void()+
  theme(axis.text.y = element_text(size = 12,  family = "arial", margin = margin(0,4,0,2 ) ),
        axis.text.x = element_text(size = 12,  family = "arial"),
        axis.title.x = element_text(size = 14, margin = margin(3,1,1,1), family = "arial" , vjust = 0.5),
        axis.title.y = element_text(size = 14, margin = margin(1,5,1,1), angle = 90, family = "arial" , hjust = 0.5 ),
        
        plot.margin = unit(c(0, 1, 1, 1), "lines"))



# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[["PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer"]])

# Enrichment plot 생성
enrichment_plot_cd4ta <- plotEnrichment(gene_set_fullgene[["PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer"]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "PD L1 expression and PD 1 checkpoint\npathway in cancer" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == "PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer" ] ,3 )* 0.95 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == "PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer" ] ,6 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == "PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer" ] ,1 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")



# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[["Treg_Cell_Differentiation"]])

# Enrichment plot 생성
enrichment_plot_cd4tb <- plotEnrichment(gene_set_fullgene[["Treg_Cell_Differentiation"]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Treg Cell Differentiation\n" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == "Treg_Cell_Differentiation" ] ,3 )* 0.95 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == "Treg_Cell_Differentiation" ] ,6 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == "Treg_Cell_Differentiation" ] ,1 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")


grid.arrange(enrichment_plot_cd4ta, enrichment_plot_cd4tb,  bar_plot, bar_plot,  scatter_plot, scatter_plot, ncol = 2, heights = c(2.5, 0.3, 1.5))


### Fig4G  - CD4T###################

cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/figure4G_GSEA_CD4T.pdf" , width = 10 , height = 4)

grid.arrange(enrichment_plot_cd4ta, enrichment_plot_cd4tb,  bar_plot, bar_plot,  scatter_plot, scatter_plot, ncol = 2, heights = c(2.5, 0.3, 1.5) )


dev.off()






gene_set_fullgene$Treg_Cells_Promote_Immunosuppression_in_Cancer_Immune_Escape
gene_set_fullgene$T_Cell_Receptor_Signaling
gene_set_fullgene$Effector_T_cell_Inactivation_in_Cancer_Immune_Escape
gene_set_fullgene$CTLA4_Inhibitory_Signaling
gene_set_fullgene$Regulation_Of_Regulatory_T_Cell_Differentiation
gene_set_fullgene$Positive_Regulation_Of_Programmed_Cell_Death


termname <- "Treg_Cells_Promote_Immunosuppression_in_Cancer_Immune_Escape"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd4t_11 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Treg Cells Promote Immunosuppression\nin Cancer Immune Escape" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.95,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               # round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,6 ) ,
                                                                               formatC(fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ], format = "e", digits = 1),
                                                                               
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,1 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")




termname <- "T_Cell_Receptor_Signaling"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd4t_12 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "T Cell Receptor Signaling\n" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.95 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               # round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,4 ) ,
                                                                               formatC(fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ], format = "e", digits = 1),
                                                                               
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,1 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")



termname <- "Effector_T_cell_Inactivation_in_Cancer_Immune_Escape"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd4t_21 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Effector T cell Inactivation\nin Cancer Immune Escape" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 ) * 0.95,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               # round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,4 ) ,
                                                                               formatC(fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ], format = "e", digits = 1),
                                                                               
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,1 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")







termname <- "CTLA4_Inhibitory_Signaling"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd4t_22 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "CTLA4 Inhibitory Signaling\n" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 ) * 0.95,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               # round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,4 ) ,
                                                                               formatC(fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ], format = "e", digits = 1),
                                                                               
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,1 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")



termname <- "Regulation_Of_Regulatory_T_Cell_Differentiation"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd4t_31 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Regulation Of Regulatory T Cell \nDifferentiation" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.95,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               # round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,3 ) ,
                                                                               round(fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ],  digits = 3),
                                                                               
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,1 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")




termname <- "Positive_Regulation_Of_Programmed_Cell_Death"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd4t_32 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Positive Regulation Of Programmed \nCell Death" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.95 ,
                                                               label = sprintf("FDR : %s\nNES : %s",
                                                                               # round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,4 ) ,
                                                                               formatC(fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ], format = "e", digits = 1),
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,1 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")









### Fig4Gsup  - CD4T###################

cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/figure4Gsup1_GSEA_CD4T.pdf" , width = 10 , height = 4)

grid.arrange(enrichment_plot_cd4t_11, enrichment_plot_cd4t_12,  bar_plot, bar_plot,  scatter_plot, scatter_plot, ncol = 2, heights = c(2.5, 0.3, 1.5) )


dev.off()


cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/figure4Gsup2_GSEA_CD4T.pdf" , width = 10 , height = 4)

grid.arrange(enrichment_plot_cd4t_21, enrichment_plot_cd4t_22,  bar_plot, bar_plot,  scatter_plot, scatter_plot, ncol = 2, heights = c(2.5, 0.3, 1.5) )


dev.off()


cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/figure4Gsup3_GSEA_CD4T.pdf" , width = 10 , height = 4)

grid.arrange(enrichment_plot_cd4t_31, enrichment_plot_cd4t_32,  bar_plot, bar_plot,  scatter_plot, scatter_plot, ncol = 2, heights = c(2.5, 0.3, 1.5) )


dev.off()








# # 원하는 항목들만 새로운 리스트에 저장
# gene_set_fullgene_tmp_cd8t <- list(
#   PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer = gene_set_fullgene$PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer,
#   Effector_T_cell_Inactivation_in_Cancer_Immune_Escape = gene_set_fullgene$Effector_T_cell_Inactivation_in_Cancer_Immune_Escape
# )

# Fig4G  ###################
## Correlation Whole Gene in CD8T #################
sc_cor_frame[3] -> genelist_tmp
genelist_tmp$CD8T_Cor -> geneList
rownames(genelist_tmp) -> names(geneList)
sort(geneList, decreasing = T) -> geneList
length(geneList)

# fgsea_results <- fgsea(pathways = gene_set_fullgene, stats = geneList  )

# fgsea_results <- fgsea(pathways = gene_set_fullgene_tmp, stats = geneList , nperm = 10000 )
# fgsea_results <- fgsea(pathways = gene_set_fullgene_tmp_cd8t, stats = geneList , nperm = 1000 )
# 
# 
# fgsea_results_CD8T <- fgsea_results
# fgsea_results_CD8T$padj


fgsea_results_CD8T <- readRDS( "/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/Figure4_Gene_list/fgsea_results_CD8T.rds")

# 결과를 데이터프레임으로 변환합니다.
fgsea_results_df <- as.data.frame(fgsea_results_CD8T)

ranked_list <- geneList

ranked_df <- data.frame(
  Gene = names(ranked_list),
  Correlation = ranked_list
)

# Bar plot 생성
bar_plot <- ggplot(ranked_df, aes(x = seq_along(Correlation), y = 1, fill = Correlation)) +
  geom_tile(height = 0.1) +  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(y = "")+
  theme(axis.title.x = element_blank(),
        plot.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 9, colour = "white" , margin = margin(0,3,0,3)),
        
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        plot.margin = unit(c(0, 1, 0, 1), "lines"))


scatter_plot <- ggplot(ranked_df, aes(x = seq_along(Correlation), y = Correlation)) +
  geom_hline(yintercept = c(-0.5,0,  0.5), color = "#EEEEEE", linetype = "solid") +
  geom_vline(xintercept = c(0,10000,20000), color = "#EEEEEE", linetype = "solid") +
  
  geom_bar(stat = "identity", size = 0.1, color = "lightgray") +
  # geom_text(aes(x = 17000, y = 0, label = "Zero"), size = 3, color = "lightgray", family = "Arial") +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.5)) +
  scale_x_continuous(breaks = seq(0, 30000, by = 10000)) +
  labs(x = "Rank", y = "Correlation") +
  
  theme_void()+
  theme(axis.text.y = element_text(size = 12,  family = "arial", margin = margin(0,4,0,2 ) ),
        axis.text.x = element_text(size = 12,  family = "arial"),
        axis.title.x = element_text(size = 14, margin = margin(3,1,1,1), family = "arial" , vjust = 0.5),
        axis.title.y = element_text(size = 14, margin = margin(1,5,1,1), angle = 90, family = "arial" , hjust = 0.5 ),
        
        plot.margin = unit(c(0, 1, 1, 1), "lines"))





# plotEnrichment(gene_set[["PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer(CD8T)"]], geneList) +
#   labs(title = "GSEA Plot: PD L1 expression and PD 1 checkpoint pathway in cancer(CD8T)")


# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[["PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer"]])

# Enrichment plot 생성
enrichment_plot_cd8ta <- plotEnrichment(gene_set_fullgene[["PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer"]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "PD L1 expression and PD 1 checkpoint\npathway in cancer" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD8T$ES[fgsea_results_CD8T$pathway == "PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer" ] ,3 )* 0.95 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               formatC(fgsea_results_CD8T$padj[fgsea_results_CD8T$pathway == "PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer"], format = "e", digits = 1),
                                                                               # round( fgsea_results_CD8T$padj[fgsea_results_CD8T$pathway == "PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer" ] ,5 ) ,
                                                                               round( fgsea_results_CD8T$NES[fgsea_results_CD8T$pathway == "PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer" ] ,1 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")




# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[["Effector_T_cell_Inactivation_in_Cancer_Immune_Escape"]])

# Enrichment plot 생성
enrichment_plot_cd8tb <- plotEnrichment(gene_set_fullgene[["Effector_T_cell_Inactivation_in_Cancer_Immune_Escape"]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Effector T cell Inactivation\nin Cancer Immune Escape" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf, y = round( fgsea_results_CD8T$ES[fgsea_results_CD8T$pathway == "Effector_T_cell_Inactivation_in_Cancer_Immune_Escape" ] ,3 )* 0.95 , 
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               # round( fgsea_results_CD8T$padj[fgsea_results_CD8T$pathway == "Effector_T_cell_Inactivation_in_Cancer_Immune_Escape" ] ,4 ) ,
                                                                               formatC(fgsea_results_CD8T$padj[fgsea_results_CD8T$pathway == "Effector_T_cell_Inactivation_in_Cancer_Immune_Escape"], format = "f", digits = 3),
                                                                               
                                                                               round( fgsea_results_CD8T$NES[fgsea_results_CD8T$pathway == "Effector_T_cell_Inactivation_in_Cancer_Immune_Escape" ] ,1 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")




grid.arrange(enrichment_plot_cd8ta, enrichment_plot_cd8tb,  bar_plot, bar_plot,  scatter_plot, scatter_plot, ncol = 2, heights = c(2.5, 0.3, 1.5))

### Fig4G - CD8T ###################
cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/figure4G_GSEA_CD8T.pdf" , width = 10 , height = 4)

grid.arrange(enrichment_plot_cd8ta, enrichment_plot_cd8tb,  bar_plot, bar_plot,  scatter_plot, scatter_plot, ncol = 2, heights = c(2.5, 0.3, 1.5))


dev.off()






gene_set_fullgene$Treg_Cells_Promote_Immunosuppression_in_Cancer_Immune_Escape
gene_set_fullgene$T_Cell_Receptor_Signaling
gene_set_fullgene$Effector_T_cell_Inactivation_in_Cancer_Immune_Escape
gene_set_fullgene$CTLA4_Inhibitory_Signaling
gene_set_fullgene$Regulation_Of_Regulatory_T_Cell_Differentiation
gene_set_fullgene$Positive_Regulation_Of_Programmed_Cell_Death


termname <- "Treg_Cells_Promote_Immunosuppression_in_Cancer_Immune_Escape"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd8t_11 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Treg Cells Promote Immunosuppression\nin Cancer Immune Escape" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,4 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")




termname <- "T_Cell_Receptor_Signaling"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd8t_12 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "T Cell Receptor Signaling" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,4 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")



termname <- "Effector_T_cell_Inactivation_in_Cancer_Immune_Escape"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd8t_21 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Effector T cell Inactivation in Cancer Immune Escape" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,4 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")







termname <- "CTLA4_Inhibitory_Signaling"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd8t_22 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "CTLA4 Inhibitory Signaling" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,4 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")



termname <- "Regulation_Of_Regulatory_T_Cell_Differentiation"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd8t_31 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Regulation_Of_Regulatory_T_Cell_Differentiation" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,4 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")




termname <- "Positive_Regulation_Of_Programmed_Cell_Death"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd8t_32 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Positive Regulation Of Programmed Cell Death" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,4 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")









### Fig4Gsup  - CD8T###################

cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/figure4Gsup1_GSEA_CD8T.pdf" , width = 10 , height = 4)

grid.arrange(enrichment_plot_cd8t_11, enrichment_plot_cd8t_12,  bar_plot, bar_plot,  scatter_plot, scatter_plot, ncol = 2, heights = c(2.5, 0.3, 1.5) )


dev.off()


cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/figure4Gsup2_GSEA_CD8T.pdf" , width = 10 , height = 4)

grid.arrange(enrichment_plot_cd8t_21, enrichment_plot_cd8t_22,  bar_plot, bar_plot,  scatter_plot, scatter_plot, ncol = 2, heights = c(2.5, 0.3, 1.5) )


dev.off()


cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/figure4Gsup3_GSEA_CD8T.pdf" , width = 10 , height = 4)

grid.arrange(enrichment_plot_cd8t_31, enrichment_plot_cd8t_32,  bar_plot, bar_plot,  scatter_plot, scatter_plot, ncol = 2, heights = c(2.5, 0.3, 1.5) )


dev.off()












# p-value 계산
p_value <- round(tpval_result_frame[rownames(mat_percent_for_temp[j,]), Split_Col[i]], digits = 4)

if (p_value < 0.01) {
  # 0.01 미만이면 지수 표기법 사용
  formatted_p_value <- formatC(p_value, format = "e", digits = 1)
  parts <- strsplit(formatted_p_value, "e")[[1]]
  mantissa <- parts[1]
  exponent <- as.numeric(parts[2])
  x_text <- sprintf("<i>P</i> = %s × 10<sup>%d</sup>", mantissa, exponent)
} else if (p_value < 0.1) {
  # 0.010 미만이면 지수 표기법 사용
  x_text <- paste("<i>P</i> =", formatC(p_value, format = "f", digits = 3))        
} else if (p_value == 1) {
  # p_value가 정확히 1인 경우에는 1.0으로 표시
  x_text <-  paste("<i>P</i> = 1.0")
} else {
  # 0.01 이상인 경우 항상 두 자리로 표시
  x_text <- paste("<i>P</i> =", formatC(p_value, format = "f", digits = 2))
}





























# Additional Term ################################################################
## GSEA_result_excel_create ######################################################

library(readxl)

# 엑셀 파일 경로 설정
file_path <- "/home/Data_Drive_8TB_3/Coloncancer/CorGSEA/CRC_GSEA_Total_filtered.xlsx"

# 엑셀 파일의 모든 시트 이름 가져오기
sheet_names <- excel_sheets(file_path)

# 각 시트별 데이터를 처리하여 리스트로 저장
sheet_gene_list <- lapply(sheet_names, function(sheet) {
  # 시트 읽기
  data <- read_excel(file_path, sheet = sheet, col_names = FALSE)
  
  # A열의 Term
  terms <- data[[1]]
  
  # C열부터 마지막 열까지 유전자들
  genes <- data[, 3:ncol(data)]
  
  # Term 별로 유전자 리스트 생성
  gene_list <- lapply(seq_along(terms), function(i) {
    term <- terms[i]
    term_genes <- genes[i, ]
    term_genes <- unlist(term_genes, use.names = FALSE)
    term_genes <- term_genes[!is.na(term_genes) & term_genes != ""] # NA 및 빈 문자열 제거
    return(term_genes)
  })
  
  names(gene_list) <- terms
  
  # 유전자 개수가 10개 미만인 Term 제외
  gene_list <- gene_list[sapply(gene_list, length) >= 10]
  
  return(gene_list)
})

# 시트 이름을 리스트의 이름으로 설정
names(sheet_gene_list) <- sheet_names

# 결과 출력
sheet_gene_list



# library(writexl)
# output_file <- "/home/Data_Drive_8TB_3/Coloncancer/CorGSEA/CRC_GSEA_Total_filtered_output.xlsx"
# 
# # 각 시트 데이터를 변환 및 저장
# sheet_data <- lapply(sheet_gene_list, function(gene_list) {
#   # Term과 유전자 리스트를 데이터프레임으로 변환
#   do.call(rbind, lapply(names(gene_list), function(term) {
#     data.frame(Term = term, Genes = paste(gene_list[[term]], collapse = ", "), stringsAsFactors = FALSE)
#   }))
# })
# 
# # 엑셀 파일로 저장
# write_xlsx(sheet_data, path = output_file)
# 
# cat("Data has been successfully saved to:", output_file, "\n")




sc_cor_frame[2] -> genelist_tmp
genelist_tmp$CD4T_Cor -> geneList
rownames(genelist_tmp) -> names(geneList)
sort(geneList, decreasing = T) -> geneList
length(geneList)
ranked_list <- geneList
geneList->CD4TgeneList

ranked_CD4T_df <- data.frame(
  Gene = names(ranked_list),
  Correlation = ranked_list
)




sc_cor_frame[3] -> genelist_tmp
genelist_tmp$CD8T_Cor -> geneList
rownames(genelist_tmp) -> names(geneList)
sort(geneList, decreasing = T) -> geneList
length(geneList)
ranked_list <- geneList
geneList->CD8TgeneList


ranked_CD8T_df <- data.frame(
  Gene = names(ranked_list),
  Correlation = ranked_list
)










# ranked_CD4T_df의 Gene 순서를 기준으로 정렬
ranked_genes <- ranked_CD4T_df$Gene

# sheet_gene_list의 각 Term에서 유전자 순서를 재정렬하며, 유전자 개수가 10개 미만인 Term은 제외
sheet_gene_list_CD4T_sorted <- lapply(sheet_gene_list, function(gene_list) {
  filtered_gene_list <- lapply(gene_list, function(genes) {
    # ranked_genes에 따라 genes를 정렬
    sorted_genes <- genes[genes %in% ranked_genes] # ranked_genes에 없는 유전자는 제외
    sorted_genes <- ranked_genes[ranked_genes %in% sorted_genes] # ranked_genes의 순서로 정렬
    return(sorted_genes)
  })
  
  # 유전자 개수가 10개 이상인 Term만 필터링
  filtered_gene_list <- filtered_gene_list[sapply(filtered_gene_list, length) >= 15]
  return(filtered_gene_list)
})

# 결과 출력
sheet_gene_list_CD4T_sorted




library(openxlsx)

# PosCorr, NegCorr, NotSig 기준 값
pos_corr_genes <- ranked_CD4T_df$Gene[ranked_CD4T_df$Correlation > 0.25]
neg_corr_genes <- ranked_CD4T_df$Gene[ranked_CD4T_df$Correlation < -0.25]
not_sig_genes <- ranked_CD4T_df$Gene[ranked_CD4T_df$Correlation >= -0.25 & ranked_CD4T_df$Correlation <= 0.25]

# 리스트를 데이터프레임 형식으로 변환
sheet_data <- lapply(sheet_gene_list_CD4T_sorted, function(gene_list) {
  do.call(rbind, lapply(names(gene_list), function(term) {
    genes <- paste(gene_list[[term]], collapse = ", ")
    gene_count <- length(gene_list[[term]])
    
    pos_corr <- intersect(gene_list[[term]], pos_corr_genes)
    pos_corr_genes_str <- paste(pos_corr, collapse = ", ")
    pos_corr_count <- length(pos_corr)
    
    not_sig <- intersect(gene_list[[term]], not_sig_genes)
    not_sig_genes_str <- paste(not_sig, collapse = ", ")
    not_sig_count <- length(not_sig)
    
    neg_corr <- intersect(gene_list[[term]], neg_corr_genes)
    neg_corr_genes_str <- paste(neg_corr, collapse = ", ")
    neg_corr_count <- length(neg_corr)
    
    data.frame(
      Term = term,
      Genes = genes,
      GeneCount = gene_count,
      PosCorr = pos_corr_genes_str,
      PosCorrCount = pos_corr_count,
      NotSig = not_sig_genes_str,
      NotSigCount = not_sig_count,
      NegCorr = neg_corr_genes_str,
      NegCorrCount = neg_corr_count,
      stringsAsFactors = FALSE
    )
  }))
})

# 엑셀 파일 생성
output_file <- "/home/Data_Drive_8TB_3/Coloncancer/CorGSEA/CRC_GSEA_Total_filtered_CD4T_sorted.xlsx"
wb <- createWorkbook()

# 시트별 데이터 추가
for (sheet_name in names(sheet_data)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, sheet_data[[sheet_name]])
  
  # 중앙 정렬 스타일 생성
  center_style <- createStyle(halign = "center", valign = "center")
  
  # 중앙 정렬 적용 (GeneCount, PosCorrCount, NotSigCount, NegCorrCount 열)
  cols_to_align <- c("GeneCount", "PosCorrCount", "NotSigCount", "NegCorrCount")
  for (col_name in cols_to_align) {
    col_index <- which(names(sheet_data[[sheet_name]]) == col_name)
    if (length(col_index) > 0) {
      addStyle(
        wb,
        sheet = sheet_name,
        style = center_style,
        cols = col_index,
        rows = 2:(nrow(sheet_data[[sheet_name]]) + 1),
        gridExpand = TRUE
      )
    }
  }
}

# 엑셀 파일 저장
saveWorkbook(wb, file = output_file, overwrite = TRUE)

cat("Data with center-aligned counts has been saved to:", output_file, "\n")



# PosCorr, NegCorr, NotSig 기준 값
pos_corr_genes <- ranked_CD4T_df$Gene[ranked_CD4T_df$Correlation > 0.25]
neg_corr_genes <- ranked_CD4T_df$Gene[ranked_CD4T_df$Correlation < -0.25]
not_sig_genes <- ranked_CD4T_df$Gene[ranked_CD4T_df$Correlation >= -0.25 & ranked_CD4T_df$Correlation <= 0.25]

# 리스트를 데이터프레임 형식으로 변환
sheet_data <- lapply(sheet_gene_list_CD4T_sorted, function(gene_list) {
  do.call(rbind, lapply(names(gene_list), function(term) {
    # B 열: 모든 유전자
    genes <- paste(gene_list[[term]], collapse = ", ")
    gene_count <- length(gene_list[[term]])
    
    # C 열: Correlation > 0.25 (PosCorr)
    pos_corr <- intersect(gene_list[[term]], pos_corr_genes)
    pos_corr_genes_str <- paste(pos_corr, collapse = ", ")
    pos_corr_count <- length(pos_corr)
    
    # D 열: -0.25 <= Correlation <= 0.25 (NotSig)
    not_sig <- intersect(gene_list[[term]], not_sig_genes)
    not_sig_genes_str <- paste(not_sig, collapse = ", ")
    not_sig_count <- length(not_sig)
    
    # E 열: Correlation < -0.25 (NegCorr)
    neg_corr <- intersect(gene_list[[term]], neg_corr_genes)
    neg_corr_genes_str <- paste(neg_corr, collapse = ", ")
    neg_corr_count <- length(neg_corr)
    
    # 데이터프레임으로 변환
    data.frame(
      Term = term,
      Genes = genes,
      GeneCount = gene_count,
      PosCorr = pos_corr_genes_str,
      PosCorrCount = pos_corr_count,
      NotSig = not_sig_genes_str,
      NotSigCount = not_sig_count,
      NegCorr = neg_corr_genes_str,
      NegCorrCount = neg_corr_count,
      stringsAsFactors = FALSE
    )
  }))
})

# Term 필터링 (PosCorrCount > NegCorrCount인 Term들만 포함)
filtered_sheet_data <- lapply(sheet_data, function(data) {
  subset(data, PosCorrCount > NegCorrCount)
})

filtered_sheet_data -> filtered_sheet_data_tmp 



for (i in 1:length(filtered_sheet_data)){
  
  tmp_fgsea <- fgsea(pathways = sheet_gene_list_CD4T_sorted[[i]], stats = geneList , nperm = 10000 )
  
  # fgsea 결과를 데이터프레임으로 변환
  tmp_fgsea_df <- data.frame(
    Term = tmp_fgsea$pathway,
    NES = tmp_fgsea$NES,
    padj = tmp_fgsea$padj,
    stringsAsFactors = FALSE
  )
  
  # 원본 데이터와 fgsea 결과 병합
  filtered_sheet_data[[i]] <- filtered_sheet_data[[i]] %>%
    left_join(tmp_fgsea_df, by = "Term")
  
  # 소수점 처리
  filtered_sheet_data[[i]]$NES <- round(filtered_sheet_data[[i]]$NES, 2)
  filtered_sheet_data[[i]]$padj <- signif(filtered_sheet_data[[i]]$padj, 3)

  
}
  

  
# 엑셀 파일 생성
output_file <- "/home/Data_Drive_8TB_3/Coloncancer/CorGSEA/CRC_GSEA_Total_filtered_CD4T_Final.xlsx"
wb <- createWorkbook()

# 필터링된 데이터를 시트별로 추가
for (sheet_name in names(filtered_sheet_data)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, filtered_sheet_data[[sheet_name]])
  
  # 중앙 정렬 스타일 생성
  center_style <- createStyle(halign = "center", valign = "center")
  
  # 중앙 정렬 적용 (GeneCount, PosCorrCount, NotSigCount, NegCorrCount 열)
  cols_to_align <- c("GeneCount", "PosCorrCount", "NotSigCount", "NegCorrCount")
  for (col_name in cols_to_align) {
    col_index <- which(names(filtered_sheet_data[[sheet_name]]) == col_name)
    if (length(col_index) > 0) {
      addStyle(
        wb,
        sheet = sheet_name,
        style = center_style,
        cols = col_index,
        rows = 2:(nrow(filtered_sheet_data[[sheet_name]]) + 1),
        gridExpand = TRUE
      )
    }
  }
}

# 엑셀 파일 저장
saveWorkbook(wb, file = output_file, overwrite = TRUE)

cat("Filtered data has been saved to:", output_file, "\n")










# ranked_CD8T_df Gene 순서를 기준으로 정렬
ranked_genes <- ranked_CD8T_df$Gene

# sheet_gene_list의 각 Term에서 유전자 순서를 재정렬하며, 유전자 개수가 10개 미만인 Term은 제외
sheet_gene_list_CD8T_sorted <- lapply(sheet_gene_list, function(gene_list) {
  filtered_gene_list <- lapply(gene_list, function(genes) {
    # ranked_genes에 따라 genes를 정렬
    sorted_genes <- genes[genes %in% ranked_genes] # ranked_genes에 없는 유전자는 제외
    sorted_genes <- ranked_genes[ranked_genes %in% sorted_genes] # ranked_genes의 순서로 정렬
    return(sorted_genes)
  })
  
  # 유전자 개수가 10개 이상인 Term만 필터링
  filtered_gene_list <- filtered_gene_list[sapply(filtered_gene_list, length) >= 15]
  return(filtered_gene_list)
})

# 결과 출력
sheet_gene_list_CD8T_sorted



# pos_corr_genes <- ranked_CD8T_df$Gene[ranked_CD8T_df$Correlation > 0.25]
# neg_corr_genes <- ranked_CD8T_df$Gene[ranked_CD8T_df$Correlation < -0.25]
# not_sig_genes <- ranked_CD8T_df$Gene[ranked_CD8T_df$Correlation >= -0.25 & ranked_CD8T_df$Correlation <= 0.25]
# 
# # 리스트를 데이터프레임 형식으로 변환하여 저장
# output_file <- "/home/Data_Drive_8TB_3/Coloncancer/CorGSEA/CRC_GSEA_Total_filtered_CD8T_sorted.xlsx"
# 
# sheet_data <- lapply(sheet_gene_list_CD8T_sorted, function(gene_list) {
#   do.call(rbind, lapply(names(gene_list), function(term) {
#     # B 열: 모든 유전자
#     genes <- paste(gene_list[[term]], collapse = ", ")
#     
#     # C 열: Correlation > 0.25 (PosCorr)
#     pos_corr <- paste(intersect(gene_list[[term]], pos_corr_genes), collapse = ", ")
#     
#     # E 열: Correlation < -0.25 (NegCorr)
#     neg_corr <- paste(intersect(gene_list[[term]], neg_corr_genes), collapse = ", ")
#     
#     # D 열: -0.25 <= Correlation <= 0.25 (NotSig)
#     not_sig <- paste(intersect(gene_list[[term]], not_sig_genes), collapse = ", ")
#     
#     # 데이터프레임으로 변환
#     data.frame(
#       Term = term,
#       Genes = genes,
#       PosCorr = pos_corr,
#       NotSig = not_sig,
#       NegCorr = neg_corr,
#       stringsAsFactors = FALSE
#     )
#   }))
# })
# 
# # 엑셀 파일로 저장
# library(writexl)
# write_xlsx(sheet_data, path = output_file)







# PosCorr, NegCorr, NotSig 기준 값
pos_corr_genes <- ranked_CD8T_df$Gene[ranked_CD8T_df$Correlation > 0.25]
neg_corr_genes <- ranked_CD8T_df$Gene[ranked_CD8T_df$Correlation < -0.25]
not_sig_genes <- ranked_CD8T_df$Gene[ranked_CD8T_df$Correlation >= -0.25 & ranked_CD8T_df$Correlation <= 0.25]

# 리스트를 데이터프레임 형식으로 변환
sheet_data <- lapply(sheet_gene_list_CD8T_sorted, function(gene_list) {
  do.call(rbind, lapply(names(gene_list), function(term) {
    # B 열: 모든 유전자
    genes <- paste(gene_list[[term]], collapse = ", ")
    gene_count <- length(gene_list[[term]])
    
    # C 열: Correlation > 0.25 (PosCorr)
    pos_corr <- intersect(gene_list[[term]], pos_corr_genes)
    pos_corr_genes_str <- paste(pos_corr, collapse = ", ")
    pos_corr_count <- length(pos_corr)
    
    # D 열: -0.25 <= Correlation <= 0.25 (NotSig)
    not_sig <- intersect(gene_list[[term]], not_sig_genes)
    not_sig_genes_str <- paste(not_sig, collapse = ", ")
    not_sig_count <- length(not_sig)
    
    # E 열: Correlation < -0.25 (NegCorr)
    neg_corr <- intersect(gene_list[[term]], neg_corr_genes)
    neg_corr_genes_str <- paste(neg_corr, collapse = ", ")
    neg_corr_count <- length(neg_corr)
    
    # 데이터프레임으로 변환
    data.frame(
      Term = term,
      Genes = genes,
      GeneCount = gene_count,
      PosCorr = pos_corr_genes_str,
      PosCorrCount = pos_corr_count,
      NotSig = not_sig_genes_str,
      NotSigCount = not_sig_count,
      NegCorr = neg_corr_genes_str,
      NegCorrCount = neg_corr_count,
      stringsAsFactors = FALSE
    )
  }))
})

# 엑셀 파일 생성
output_file <- "/home/Data_Drive_8TB_3/Coloncancer/CorGSEA/CRC_GSEA_Total_filtered_CD8T_sorted.xlsx"
wb <- createWorkbook()

# 시트별 데이터 추가
for (sheet_name in names(sheet_data)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, sheet_data[[sheet_name]])
  
  # 중앙 정렬 스타일 생성
  center_style <- createStyle(halign = "center", valign = "center")
  
  # 중앙 정렬 적용 (GeneCount, PosCorrCount, NotSigCount, NegCorrCount 열)
  cols_to_align <- c("GeneCount", "PosCorrCount", "NotSigCount", "NegCorrCount")
  for (col_name in cols_to_align) {
    col_index <- which(names(sheet_data[[sheet_name]]) == col_name)
    if (length(col_index) > 0) {
      addStyle(
        wb,
        sheet = sheet_name,
        style = center_style,
        cols = col_index,
        rows = 2:(nrow(sheet_data[[sheet_name]]) + 1),
        gridExpand = TRUE
      )
    }
  }
}

# 엑셀 파일 저장
saveWorkbook(wb, file = output_file, overwrite = TRUE)

cat("Data with center-aligned counts has been saved to:", output_file, "\n")








# PosCorr, NegCorr, NotSig 기준 값
pos_corr_genes <- ranked_CD8T_df$Gene[ranked_CD8T_df$Correlation > 0.25]
neg_corr_genes <- ranked_CD8T_df$Gene[ranked_CD8T_df$Correlation < -0.25]
not_sig_genes <- ranked_CD8T_df$Gene[ranked_CD8T_df$Correlation >= -0.25 & ranked_CD8T_df$Correlation <= 0.25]

# 리스트를 데이터프레임 형식으로 변환
sheet_data <- lapply(sheet_gene_list_CD8T_sorted, function(gene_list) {
  do.call(rbind, lapply(names(gene_list), function(term) {
    # B 열: 모든 유전자
    genes <- paste(gene_list[[term]], collapse = ", ")
    gene_count <- length(gene_list[[term]])
    
    # C 열: Correlation > 0.25 (PosCorr)
    pos_corr <- intersect(gene_list[[term]], pos_corr_genes)
    pos_corr_genes_str <- paste(pos_corr, collapse = ", ")
    pos_corr_count <- length(pos_corr)
    
    # D 열: -0.25 <= Correlation <= 0.25 (NotSig)
    not_sig <- intersect(gene_list[[term]], not_sig_genes)
    not_sig_genes_str <- paste(not_sig, collapse = ", ")
    not_sig_count <- length(not_sig)
    
    # E 열: Correlation < -0.25 (NegCorr)
    neg_corr <- intersect(gene_list[[term]], neg_corr_genes)
    neg_corr_genes_str <- paste(neg_corr, collapse = ", ")
    neg_corr_count <- length(neg_corr)
    
    # 데이터프레임으로 변환
    data.frame(
      Term = term,
      Genes = genes,
      GeneCount = gene_count,
      PosCorr = pos_corr_genes_str,
      PosCorrCount = pos_corr_count,
      NotSig = not_sig_genes_str,
      NotSigCount = not_sig_count,
      NegCorr = neg_corr_genes_str,
      NegCorrCount = neg_corr_count,
      stringsAsFactors = FALSE
    )
  }))
})

# Term 필터링 (PosCorrCount > NegCorrCount인 Term들만 포함)
filtered_sheet_data <- lapply(sheet_data, function(data) {
  subset(data, PosCorrCount > NegCorrCount)
})


filtered_sheet_data -> filtered_sheet_data_tmp 



for (i in 1:length(filtered_sheet_data)){
  
  tmp_fgsea <- fgsea(pathways = sheet_gene_list_CD8T_sorted[[i]], stats = geneList , nperm = 10000 )
  
  # fgsea 결과를 데이터프레임으로 변환
  tmp_fgsea_df <- data.frame(
    Term = tmp_fgsea$pathway,
    NES = tmp_fgsea$NES,
    padj = tmp_fgsea$padj,
    stringsAsFactors = FALSE
  )
  
  # 원본 데이터와 fgsea 결과 병합
  filtered_sheet_data[[i]] <- filtered_sheet_data[[i]] %>%
    left_join(tmp_fgsea_df, by = "Term")
  
  # 소수점 처리
  filtered_sheet_data[[i]]$NES <- round(filtered_sheet_data[[i]]$NES, 2)
  filtered_sheet_data[[i]]$padj <- signif(filtered_sheet_data[[i]]$padj, 3)
  
  
}



# 엑셀 파일 생성
output_file <- "/home/Data_Drive_8TB_3/Coloncancer/CorGSEA/CRC_GSEA_Total_filtered_CD8T_Final.xlsx"
wb <- createWorkbook()

# 필터링된 데이터를 시트별로 추가
for (sheet_name in names(filtered_sheet_data)) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, filtered_sheet_data[[sheet_name]])
  
  # 중앙 정렬 스타일 생성
  center_style <- createStyle(halign = "center", valign = "center")
  
  # 중앙 정렬 적용 (GeneCount, PosCorrCount, NotSigCount, NegCorrCount 열)
  cols_to_align <- c("GeneCount", "PosCorrCount", "NotSigCount", "NegCorrCount")
  for (col_name in cols_to_align) {
    col_index <- which(names(filtered_sheet_data[[sheet_name]]) == col_name)
    if (length(col_index) > 0) {
      addStyle(
        wb,
        sheet = sheet_name,
        style = center_style,
        cols = col_index,
        rows = 2:(nrow(filtered_sheet_data[[sheet_name]]) + 1),
        gridExpand = TRUE
      )
    }
  }
}

# 엑셀 파일 저장
saveWorkbook(wb, file = output_file, overwrite = TRUE)

cat("Filtered data has been saved to:", output_file, "\n")






CEBPBL_indices_GMM_index <- which(CEBPB_MeanExp_Vector < 1.2 )
CEBPBH_indices_GMM_index <- which(CEBPB_MeanExp_Vector > 1.2 )
CEBPBL_indices_GMM_Human <-  paste0("H",  ifelse(CEBPBL_indices_GMM_index >= 10, CEBPBL_indices_GMM_index + 1, CEBPBL_indices_GMM_index))
CEBPBH_indices_GMM_Human <-  paste0("H",  ifelse(CEBPBH_indices_GMM_index >= 10, CEBPBH_indices_GMM_index + 1, CEBPBH_indices_GMM_index))


Idents(CD8T_CRC) <- CD8T_CRC$patients
clusters <- CD8T_CRC$patients
cluster_names <- CD8T_CRC$patients
# cluster_names <- rep("Unknown", length(clusters))
# cluster_names <- clusters

cluster_names[clusters %in% c(CEBPBL_indices_GMM_Human)] <- "CEBPB L"
cluster_names[clusters %in% c(CEBPBH_indices_GMM_Human)] <- "CEBPB H"

CD8T_CRC <- AddMetaData(CD8T_CRC, metadata = cluster_names, col.name = "CEBPB_GMM_HL")

library(patchwork)
VlnPlot(Epi_Cyling_CRC, features = "CEBPB", pt.size = 0)+ theme(legend.position = "none") + ggtitle("CEBPB in Epithelail")-> a
VlnPlot(CD8T_CRC, features = "ENTPD1", pt.size = 0) + theme(legend.position = "none") + ggtitle("ENTPD1 in CD8T") -> b

combined_plot <- a / b + plot_layout(heights = c(1, 1)) 
combined_plot


ggsave("/home/Data_Drive_8TB_3/Coloncancer/CorGSEA/ENTPD1.png", plot = combined_plot , height = 5 ,width =  10,units = "in" )


CEBPBL_indices_median_index <- which(CEBPB_MeanExp_Vector < median(CEBPB_MeanExp_Vector) )
CEBPBH_indices_median_index <- which(CEBPB_MeanExp_Vector > median(CEBPB_MeanExp_Vector) )

CEBPBL_indices_median_Human <-  paste0("H",  ifelse(CEBPBL_indices_median_index >= 10, CEBPBL_indices_median_index + 1, CEBPBL_indices_median_index))
CEBPBH_indices_median_Human <-  paste0("H",  ifelse(CEBPBH_indices_median_index >= 10, CEBPBH_indices_median_index + 1, CEBPBH_indices_median_index))




AverageExpression(CD8T_CRC , features = "ENTPD1", assays = "RNA")$RNA


mean(AverageExpression(CD8T_CRC , features = "ENTPD1", assays = "RNA")$RNA[,CEBPBH_indices_median_Human])
mean(AverageExpression(CD8T_CRC , features = "ENTPD1", assays = "RNA")$RNA[,CEBPBL_indices_median_Human])

t.test( as.numeric(AverageExpression(CD8T_CRC , features = "ENTPD1", assays = "RNA")$RNA[,CEBPBH_indices_median_Human]), 
        as.numeric( AverageExpression(CD8T_CRC , features = "ENTPD1", assays = "RNA")$RNA[,CEBPBL_indices_median_Human] )
        )

group_high <- as.numeric(AverageExpression(CD8T_CRC, features = "ENTPD1", assays = "RNA")$RNA[, CEBPBH_indices_median_Human])
group_low <- as.numeric(AverageExpression(CD8T_CRC, features = "ENTPD1", assays = "RNA")$RNA[, CEBPBL_indices_median_Human])

# 데이터프레임 생성
vln_data <- data.frame(
  Group = c(rep("High", length(group_high)), rep("Low", length(group_low))),
  Expression = c(group_high, group_low)
)

# VlnPlot 생성
vln_plot <- ggplot(vln_data, aes(x = Group, y = Expression, fill = Group)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_classic() +
  labs(title = "High vs Low Groups(Median) in CD8T", y = "ENTPD1 Expression", x = "Group") +
  theme(legend.position = "none")

# 플롯 출력
print(vln_plot)

ggsave("/home/Data_Drive_8TB_3/Coloncancer/CorGSEA/ENTPD1_median_Vln.png", plot = vln_plot , height = 4 ,width =  4 ,units = "in" )







AverageExpression(CD8T_CRC , features = "ENTPD1", assays = "RNA")$RNA
CEBPB_MeanExp_Vector


expr_vector1 <- as.numeric(AverageExpression(CD8T_CRC, features = "ENTPD1", assays = "RNA")$RNA)
expr_vector2 <- as.numeric(CEBPB_MeanExp_Vector)

# 데이터프레임 생성
correlation_data <- data.frame(
  Expression1 = expr_vector1,
  Expression2 = expr_vector2
)

# 상관 계수 계산 (유효숫자 2자리)
correlation <- cor(expr_vector1, expr_vector2, use = "complete.obs", method = "pearson")
correlation_text <- sprintf("Correlation: %.2f", correlation)

# Scatter plot 생성
scatter_plot <- ggplot(correlation_data, aes(x = Expression1, y = Expression2)) +
  geom_point(color = "blue", size = 3, alpha = 0.7) + # 점 추가
  geom_smooth(method = "lm", color = "red", se = TRUE) + # 추세선 추가
  theme_classic() +
  ylim(0, 3)+
  labs(
    title = "",
    x = "ENTPD1 in CD8T",
    y = "CEBPB in epithelial"
  ) +
  annotate("text", x = 0.2 , y = 2.7 , label = correlation_text, size = 5, hjust = 0)

# 플롯 출력
print(scatter_plot)



ggsave("/home/Data_Drive_8TB_3/Coloncancer/CorGSEA/ENTPD1_Corr.png", plot = scatter_plot , height = 4 ,width =  4 ,units = "in" )






#  TCR  ##################################


gene_set_fullgene_TCR <- list(
  "T_Cell_Receptor_Signaling_Pathway_(GO:0050852)" = c("BTNL10P", "BTNL9", "BTN3A2", "BTN3A1", "MOG", "BTN1A1", "BTNL2", "NFKBIZ", "RC3H1", "RC3H2", 
                                          "BTN3A3", "BTNL3", "THEMIS", "BTNL8", "TEC", "TNFRSF21", "PTPRC", "NCOR2", "THEMIS2", 
                                          "LOC102723996", "PIK3CA", "GNAO1", "CD28", "RBCK1", "LCP2", "CD247", "SKAP1", "ITK", 
                                          "NFKBID", "UBE2N", "KHDRBS1", "STOML2", "VTCN1", "DENND1B", "PTPN22", "CD3G", "THY1", 
                                          "CD3E", "PTPRJ", "CD3D", "RNF31", "TRBC1", "FCHO1", "ORC3", "TRBC2", "PDE4B", "CTLA4", 
                                          "RIPK2", "TRDC", "MAP3K7", "FYB1", "EIF2B5", "EIF2B4", "EIF2B3", "FYB2", "BTN2A2", 
                                          "BTN2A1", "RFTN1", "PDE4D", "HLA-A", "LCK", "ZAP70", "WNK1", "CD8B", "CD8A", "TRAF6", 
                                          "BTK", "PRKD2", "HLA-DRB3", "HLA-DRB1", "CD276", "HHLA2", "TXK", "BCL10", "PIK3CD", 
                                          "GATA3", "SPPL3", "MALT1", "PRAM1", "TRGC1", "IKBKB", "INPP5D", "TRGC2", "BTN2A3P", 
                                          "ZC3H12A", "HLA-DPA1", "PLCG1", "CSK", "FYN", "PLCG2", "EIF2B2", "ERMAP", "IKBKG", 
                                          "EIF2B1", "ICOSLG"),
  
  "Modulators_Of_TCR_Signaling_And_T_Cell_Activation_WP5072" = c("RHOH", "VAV1", "NFKBIA", "CD5", "RASA2", "AGO1", "CD28", "REL", "LCP2", "CD247", 
                                                                 "PTPN6", "ITK", "LAT", "RELA", "SMARCB1", "DGKA", "NDUFB10", "TNFAIP3", "ZFP36L1", 
                                                                 "RNF7", "CBLB", "CD3G", "PIK3R2", "CD3E", "CD3D", "SOCS1", "AKT1", "GRAP2", 
                                                                 "UBASH3A", "RPRD1B", "MAP3K8", "SH2B3", "CARD11", "MEF2D", "MAP4K1", "CHUK", 
                                                                 "PDPK1", "SH2D1A", "NFKB1", "DGKZ", "FIBP", "LCK", "ZAP70", "PIK3R1", "CD8A", 
                                                                 "TRAF6", "GRB2", "PRKCQ", "MAP3K14", "CDKN1B", "CUL5", "BCL10", "MALT1", 
                                                                 "TMEM222", "IKBKB", "PCBP2", "GNA13", "ELOB", "PLCG1", "ARIH2", "IKBKG"),
  
  "Downstream_TCR_Signaling_R-HSA-202424" = c("FBXW11", "PSMB9", "NFKBIA", "PSMC6", "CD4", "PIK3CB", "PSMC4", "PIK3CA", "PSMC5", "PSMC2", 
                                              "PSMC3", "HLA-DRA", "PSMC1", "CD247", "HLA-DQB2", "UBE2N", "PSMD11", "PSMD10", "RELA", 
                                              "PSMD13", "PSMD12", "PSMD14", "UBE2D1", "UBE2D2", "PIK3R2", "CD3E", "TARP", "CD3D", 
                                              "PSMB11", "PSMB10", "PSMD9", "PSMD7", "PSMD8", "PSMD5", "PSMD6", "PSMD3", "PSMD4", "PSMD2", 
                                              "RIPK2", "MAP3K7", "BTRC", "RPS27A", "HLA-DQA2", "HLA-DQA1", "CARD11", "HLA-DRB5", "CHUK", 
                                              "PDPK1", "NFKB1", "PSMA6", "PSMA7", "CDC34", "LCK", "PSMA4", "PSMA5", "SEM1", "PIK3R1", 
                                              "PSMA2", "PSMA3", "PSME4", "PSMA1", "TRAF6", "PSME2", "PSME3", "UBE2V1", "PRKCQ", "PSME1", 
                                              "TAB2", "HLA-DRB4", "UBA52", "HLA-DRB3", "HLA-DRB1", "BCL10", "CUL1", "PTEN", "PSMA8", 
                                              "MALT1", "PSMB7", "IKBKB", "INPP5D", "PSMB8", "PSMB5", "PSMB6", "PSMB3", "UBC", "PSMB4", 
                                              "UBB", "PSMB1", "PSMB2", "HLA-DPA1", "PSMF1", "TRAT1", "IKBKG", "SKP1"),
  
  "TCR_Signaling_R-HSA-202403" = c("FBXW11", "PSMB9", "PSMC6", "CD4", "PSMC4", "PSMC5", "PSMC2", "PSMC3", "PSMC1", "SPNS1", "UBE2N", "RELA", 
                                   "UBE2D1", "UBE2D2", "PSMB11", "PSMB10", "PSMD9", "PSMD7", "PSMD8", "PSMD5", "PSMD6", "PSMD3", "PSMD4", 
                                   "GRAP2", "PSMD2", "MAP3K7", "BTRC", "FYB1", "CARD11", "HLA-DRB5", "NFKB1", "VASP", "SEM1", "PIK3R1", 
                                   "PSME4", "TRAF6", "PSME2", "PSME3", "PRKCQ", "PSME1", "TAB2", "HLA-DRB4", "HLA-DRB3", "HLA-DRB1", "BCL10", 
                                   "CUL1", "INPP5D", "UBC", "UBB", "PLCG1", "PSMF1", "PLCG2", "SKP1", "ENAH", "NFKBIA", "PTPRC", "PIK3CB", 
                                   "PIK3CA", "HLA-DRA", "LCP2", "CD247", "ITK", "ARHGEF7", "HLA-DQB2", "PSMD11", "PSMD10", "PSMD13", 
                                   "PSMD12", "PSMD14", "PTPN22", "PIK3R2", "CD3E", "PTPRJ", "TARP", "CD3D", "RIPK2", "RPS27A", "HLA-DQA2", 
                                   "HLA-DQA1", "PAG1", "CHUK", "PDPK1", "PSMA6", "PSMA7", "CDC34", "LCK", "PSMA4", "ZAP70", "PSMA5", "PSMA2", 
                                   "PSMA3", "PSMA1", "PKN1", "UBE2V1", "EVL", "PKN2", "UBA52", "CD101", "PTEN", "WAS", "PSMA8", "MALT1", 
                                   "PSMB7", "IKBKB", "PSMB8", "PSMB5", "PSMB6", "PSMB3", "PSMB4", "PSMB1", "PSMB2", "HLA-DPA1", "TRAT1", 
                                   "CSK", "IKBKG", "NCK1"),
  
  "Phosphorylation_Of_CD3_And_TCR_Zeta_Chains_R-HSA-202427" = c("HLA-DRB5", "PTPN22", "CD3E", "PTPRJ", "TARP", "CD3D", "LCK", "PTPRC", 
                                                                "CD4", "HLA-DPA1", "HLA-DRA", "CSK", "CD247", "HLA-DRB4", "HLA-DRB3", 
                                                                "HLA-DQA2", "HLA-DQA1", "HLA-DQB2", "HLA-DRB1", "PAG1"),
  
  
  
  
  "T_Cell_Receptor_Signaling" = c("RHOA", "IL2", "MAPK14", "VAV3", "IL5", "IL4", "NFKBIA", "PTPRC", "CD4", "IFNG", "CD28", "LCP2", 
                                  "CD247", "PDCD1", "PTPN6", "LAT", "RAF1", "SOS1", "CD86", "TRB", "TRA", "CD80", "CD3G", "CD3E", 
                                  "CD3D", "ITPR1", "PAK1", "AKT1", "GRAP2", "CTLA4", "FCER1G", "MAP3K7", "MAP3K8", "ICOS", "CARD11", 
                                  "MAP2K2", "MAP2K1", "PDPK1", "LCK", "ZAP70", "TRAF6", "GRB2", "PRKCQ", "CD40LG", "MAP3K14", "BCL10", 
                                  "MALT1", "MAPK9", "MAPK3", "MAPK1", "PLCG1", "FYN", "NCK1", "MAP2K7"),
  
  "T_Cell_Receptor_to_AP1_Signaling" = c("CARD11", "MAP2K2", "MAP2K1", "LCK", "PTPRC", "ZAP70", "CD4", "TRAF6", "GRB2", "PRKCQ", 
                                         "CD247", "PDCD1", "PTPN6", "LAT", "RAF1", "SOS1", "TRB", "TRA", "BCL10", "CD3G", "CD3E", 
                                         "CD3D", "MALT1", "MAPK9", "MAPK3", "GRAP2", "CTLA4", "MAPK1", "PLCG1", "FYN", "FCER1G", 
                                         "MAP3K7", "MAP2K7"),
  
  "T_Cell_Receptor_to_ATF_CREB_Signaling" = c("MAP2K4", "MAP3K1", "MAP2K3", "CD72", "IL16", "VAV1", "CDC42", "MAPK14", "CREB1", 
                                              "LCK", "PTPRC", "ZAP70", "CD4", "CD8B", "CD8A", "MAPKAPK2", "MAPKAPK3", "CD28", 
                                              "MAP3K11", "LCP2", "PDCD1", "PTPN6", "CD22", "ATF2", "CD86", "ATF1", "CD80", 
                                              "RPS6KA5", "CTLA4", "RAC1", "FYN", "MAP2K6"),
  
  "T_Cell_Receptor_to_CREBBP_Signaling" = c("CREBBP", "TRB", "TRA", "CD3G", "CAMKK1", "CD3E", "CAMKK2", "CD3D", "ITPR1", "LCK", 
                                            "PTPRC", "ZAP70", "CD4", "CAMK4", "GRAP2", "CTLA4", "PLCG1", "FYN", "GRB2", "CD247", 
                                            "FCER1G", "PDCD1", "PTPN6", "LAT"),
  
  "T_Cell_Receptor_to_NFATC_Signaling" = c("NFATC4", "NFATC3", "NFATC2", "NFATC1", "MAPK14", "LCK", "PTPRC", "ZAP70", "CD4", "GRB2", 
                                           "CD247", "PDCD1", "PTPN6", "LAT", "NFAT5", "TRB", "TRA", "CD3G", "CD3E", "CD3D", 
                                           "ITPR1", "GRAP2", "CTLA4", "PLCG1", "FYN", "FCER1G"),
  
  "T_Cell_Receptor_to_NFkB_Signaling" = c("CARD11", "LYN", "CHUK", "IL16", "PDPK1", "PRKCB", "VAV1", "NFKBIA", "LCK", "ZAP70", 
                                          "CD4", "TRAF6", "CD28", "GRB2", "LCP2", "PRKCQ", "CD247", "MAP3K14", "PTPN6", "ITK", 
                                          "LAT", "CD86", "TRB", "TRA", "CD80", "BCL10", "CD3G", "CD3E", "CD3D", "ITPR1", 
                                          "MALT1", "IKBKB", "AKT1", "GRAP2", "CTLA4", "PLCG1", "FYN", "FCER1G", "MAP3K7", 
                                          "MAP3K8", "ICOS"),
  
  "T_Cell_Receptor_to_STAT_Signaling" = c("CD86", "TRB", "TRA", "STAT1", "IL16", "STAT3", "CD3G", "CD3E", "CD3D", "LCK", "ZAP70", 
                                          "CD4", "STAT6", "CTLA4", "FYN", "CD247", "FCER1G", "JAK3", "PDCD1", "PTPN6", "JAK2"),
  
  
  "T_Cell_Receptor_signaling_pathway" =  c("BUB1B-PAK6", "MAPK13", "MAPK14", "MAPK11", "TEC", "MAPK12", "CD4", "MAPK10", 
                                          "PDCD1", "PTPN6", "RAF1", "RELA", "CSF2", "TNF", "CBLB", "NRAS", "PAK2", "PAK1", 
                                          "AKT3", "AKT1", "AKT2", "GRAP2", "MAP3K7", "PAK4", "MAP3K8", "PAK3", "PAK6", 
                                          "PAK5", "CARD11", "MAP2K2", "MAP2K1", "CD8B2", "NFKB1", "DLG1", "PIK3R1", "CD8B", 
                                          "CD8A", "GRB2", "PRKCQ", "MAP3K14", "BCL10", "PPP3R1", "PPP3R2", "RASGRP1", 
                                          "PLCG1", "MAP2K7", "RHOA", "IL10", "JUN", "NFATC3", "NFATC2", "NFATC1", "VAV1", 
                                          "VAV2", "IL2", "VAV3", "IL5", "NFKBIB", "IL4", "NFKBIA", "PTPRC", "PIK3CB", 
                                          "PIK3CA", "CDK4", "IFNG", "CD28", "LCP2", "CD247", "SOS2", "NFKBIE", "ITK", 
                                          "LAT", "SOS1", "GSK3B", "PIK3R3", "CD3G", "PIK3R2", "CD3E", "CD3D", "PPP3CA", 
                                          "PPP3CB", "PPP3CC", "CTLA4", "ICOS", "CHUK", "PDPK1", "FOS", "CDC42", "LCK", 
                                          "ZAP70", "KRAS", "CD40LG", "PIK3CD", "MALT1", "IKBKB", "MAPK9", "MAPK8", "MAPK3", 
                                          "HRAS", "MAPK1", "FYN", "IKBKG", "NCK1", "NCK2"),
  
  "T_cell_receptor_downstream_signaling" = c(
    "CARD11", "RELA", "CHUK", "PDPK1", "BCL10", "PTEN", "NFKB1", "MALT1", 
    "IKBKB", "INPP5D", "NFKBIA", "LCK", "CD4", "TRAF6", "RIPK2", "TRAT1", 
    "UBE2V1", "PRKCQ", "TAB2", "IKBKG", "MAP3K7", "UBE2N"
  ),
  
  "T_cell_receptor_naive_CD4_T_cells" =  c(
    "GAB2", "PTPN11", "VAV1", "PTPRC", "CD4", "STIM1", "TRPV6", "HLA-DRA", "CD28", 
    "ORAI1", "LCP2", "CD247", "PTPN6", "ITK", "LAT", "SOS1", "CD86", "TRB", "TRA", 
    "SHC1", "CD80", "CD3G", "CD3E", "CD3D", "NRAS", "RAP1A", "SH3BP2", "STK39", 
    "AKT1", "FLNA", "GRAP2", "MAP3K8", "PAG1", "CARD11", "MAP4K1", "DBNL", "CHUK", 
    "PDPK1", "PRKCE", "PRKCB", "PRKCA", "CDC42", "LCK", "ZAP70", "TRAF6", "GRB2", 
    "KRAS", "PRKCQ", "MAP3K14", "HLA-DRB1", "BCL10", "PTEN", "WAS", "CBL", "RASGRP2", 
    "FYB", "MALT1", "SLA2", "IKBKB", "RASGRP1", "RASSF5", "HRAS", "PLCG1", "CSK", 
    "FYN", "IKBKG", "NCK1"
  ),
  
  "T_cell_receptor_naive_CD8_T_cells" = c(
    "VAV1", "PTPRC", "STIM1", "TRPV6", "CD28", "ORAI1", "LCP2", "CD247", "PTPN6", 
    "LAT", "SOS1", "CD86", "TRB", "TRA", "SHC1", "CD80", "CD3G", "CD3E", "CD3D", 
    "NRAS", "RAP1A", "AKT1", "GRAP2", "B2M", "MAP3K8", "PAG1", "CARD11", "CHUK", 
    "PDPK1", "PRKCE", "HLA-A", "PRKCB", "PRKCA", "LCK", "ZAP70", "CD8B", "CD8A", 
    "TRAF6", "GRB2", "KRAS", "PRKCQ", "MAP3K14", "BCL10", "PRF1", "CBL", "RASGRP2", 
    "MALT1", "IKBKB", "RASGRP1", "RASSF5", "HRAS", "PLCG1", "CSK", "FYN", "IKBKG"
  ),
  
  
  "T_cell_receptor_signaling_pathway" = c(
    "RALBP1", "MAP3K1", "MAPK13", "MAPK14", "MAPK11", "TEC", "MAPK12", "CD4", 
    "PTPN7", "CALM2", "PDCD1", "CALM3", "PTPN6", "CALM1", "RAF1", "RELA", "CSF2", 
    "ASAP2", "CBLC", "TNF", "CBLB", "NRAS", "PAK2", "PAK1", "ASAP1", "AKT3", 
    "AKT1", "AKT2", "CHP2", "GRAP2", "PAK7", "CHP1", "MAP3K7", "PAK4", "MAP3K8", 
    "PAK3", "PAK6", "MAP2K4", "CARD11", "MAP2K2", "MAP2K1", "PRKCB", "PRKCA", 
    "NFKB1", "ARFGAP3", "ARFGAP1", "DLG1", "PIK3R1", "CD8B", "CD8A", "GRB2", 
    "PRKCQ", "MAP3K14", "HLA-DRB1", "NCF2", "BCL10", "ARHGAP1", "ARHGAP6", "ELK1", 
    "ARHGAP4", "ARHGAP5", "PPP3R1", "PPP3R2", "RASGRP1", "PLCG1", "PDK1", "MAP2K7", 
    "RHOA", "IL10", "JUN", "NFATC4", "NFATC3", "ARAP1", "NFATC2", "ARAP2", 
    "NFATC1", "VAV1", "VAV2", "IL2", "VAV3", "IL5", "NFKBIB", "IL4", "NFKBIA", 
    "RASA1", "PTPRC", "PIK3CB", "PIK3CA", "CDK4", "IFNG", "HLA-DRA", "CD28", 
    "LCP2", "CD247", "SOS2", "NFKBIE", "ITK", "LAT", "SOS1", "GSK3B", "NFAT5", 
    "TRB", "TRA", "SHC1", "PIK3R5", "PIK3R3", "CD3G", "PIK3R2", "CD3E", "CD3D", 
    "PPP3CA", "PPP3CB", "PPP3CC", "CTLA4", "ICOS", "CHUK", "FOS", "CDC42", "LCK", 
    "ZAP70", "KRAS", "CD40LG", "CAMK2B", "PIK3CD", "CBL", "MALT1", "PIK3CG", 
    "IKBKB", "MAPK9", "MAPK8", "CHN1", "MAPK3", "HRAS", "MAPK1", "RAC1", "FYN", 
    "IKBKG", "NCK1", "NCK2"
  ),
  
  "T_cell_receptor_JNK_pathway" =c(
    "MAP2K4", "MAP3K1", "JUN", "MAP4K1", "DBNL", "PRKCB", "CRKL", "MAPK8", "GRAP2", 
    "CRK", "LCP2", "MAP3K7", "MAP3K8", "LAT"
  ),
  
  "T_cell_receptor_Ras_pathway" = c(
    "MAP2K1", "FOS", "PRKCB", "BRAF", "PRKCA", "ELK1", "NRAS", "MAPK3", "HRAS", 
    "MAPK1", "PTPN7", "KRAS", "MAP3K8", "RAF1"
  ),
  
  
  "T_cell_signal_transduction" = c(
    "MAP3K3", "MAP3K4", "MAP3K1", "MAP3K2", "NFATC4", "NFATC3", "NFATC2", "NFATC1", 
    "MAPK15", "VAV1", "MAPK13", "VAV2", "MAPK14", "VAV3", "MAPK11", "LAT2", 
    "MAPK12", "PTPRC", "MAPK10", "CD28", "LCP2", "SOS2", "ITK", "LAT", "RAF1", 
    "SOS1", "NFAT5", "TRB", "TRA", "TRD", "TRG", "CBLC", "CBLB", "CD3G", "CD3E", 
    "CD3D", "PPP3CA", "NRAS", "PAK1", "GRAP2", "CTLA4", "MAP3K5", "MAP2K4", 
    "PRKCH", "MAP2K5", "PRKCG", "MAP2K2", "MAP2K3", "PRKCI", "PRKCD", "MAP2K1", 
    "PRKCE", "PRKCB", "PRKCA", "NFKB1", "NFKB2", "PAdVCgp06", "LCK", "PRKCZ", 
    "ZAP70", "GRB2", "KRAS", "PRKCQ", "RASGRP3", "CBL", "RASGRP2", "RASGRP4", 
    "PPP3R1", "MAPK9", "RASGRP1", "MAPK8", "PRDX1", "MAPK3", "HRAS", "MAPK1", 
    "PLCG1", "CSK", "MAPK7", "MAPK6", "MAP2K6", "NCK1", "MAPK4", "NCK2"
  )
  

  
  
  
  
)





# Correlation Whole Gene in CD4T #################
sc_cor_frame[2] -> genelist_tmp
genelist_tmp$CD4T_Cor -> geneList
rownames(genelist_tmp) -> names(geneList)
sort(geneList, decreasing = T) -> geneList

fgsea_results_CD4T_TCR <- fgsea(pathways = gene_set_fullgene_TCR , stats = geneList , nperm = 10000 )

# 결과를 데이터프레임으로 변환합니다.
fgsea_results_CD4T_TCR <- as.data.frame(fgsea_results_CD4T_TCR)


ranked_list <- geneList


ranked_df <- data.frame(
  Gene = names(ranked_list),
  Correlation = ranked_list
)

# Bar plot 생성
bar_plot <- ggplot(ranked_df, aes(x = seq_along(Correlation), y = 1, fill = Correlation)) +
  geom_tile(height = 0.1) +  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(y = "")+
  theme(axis.title.x = element_blank(),
        plot.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6, colour = "white"),
        
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        plot.margin = unit(c(0, 1, 0, 1), "lines"))


# scatter_plot <- ggplot(ranked_df, aes(x = seq_along(Correlation), y = Correlation)) +
#   geom_point(size = 0.1) +
#   theme_minimal() +
#   labs(x = "Rank", y = "Correlation")+
#   theme(plot.margin = unit(c(0, 1, 1, 1), "lines"))


scatter_plot <- ggplot(ranked_df, aes(x = seq_along(Correlation), y = Correlation)) +
  geom_hline(yintercept = c(-0.5,0,  0.5), color = "#EEEEEE", linetype = "solid") +
  geom_vline(xintercept = c(0,10000,20000), color = "#EEEEEE", linetype = "solid") +
  
  geom_bar(stat = "identity", size = 0.1, color = "lightgray") +
  # geom_text(aes(x = 17000, y = 0, label = "Zero"), size = 3, color = "lightgray", family = "Arial") +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.5)) +
  scale_x_continuous(breaks = seq(0, 30000, by = 10000)) +
  labs(x = "Rank", y = "Correlation") +
  
  theme_void()+
  theme(axis.text.y = element_text(size = 8,  family = "arial", margin = margin(0,4,0,0 ) ),
        axis.text.x = element_text(size = 8,  family = "arial"),
        axis.title.x = element_text(size = 11, margin = margin(3,1,1,1), family = "arial" , vjust = 0.5),
        axis.title.y = element_text(size = 11, margin = margin(1,5,1,1), angle = 90, family = "arial" , hjust = 0.5 ),
        
        
        plot.margin = unit(c(0, 1, 1, 1), "lines"))



# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene_TCR[["Regulation_Of_Regulatory_T_Cell_Differentiation"]])

# Enrichment plot 생성
enrichment_plot <- plotEnrichment(gene_set_fullgene_TCR[["Regulation_Of_Regulatory_T_Cell_Differentiation"]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "GSEA Plot: Regulation Of Regulatory T Cell Differentiation" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title.x = element_blank(),   
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 8),
        
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) 



grid.arrange(enrichment_plot, bar_plot, scatter_plot, ncol = 1, heights = c(3, 0.5, 2))





names(gene_set_fullgene_TCR)[1] -> tcr_pathwayname

plotEnrichment(gene_set_fullgene_TCR[[tcr_pathwayname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title =  gsub("_", " ", tcr_pathwayname)  ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf, y = round( fgsea_results_CD4T_TCR$ES[fgsea_results_CD4T_TCR$pathway == tcr_pathwayname ] ,3 )* 0.9 , 
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T_TCR$padj[fgsea_results_CD4T_TCR$pathway == tcr_pathwayname ] ,4 ) ,
                                                                               round( fgsea_results_CD4T_TCR$NES[fgsea_results_CD4T_TCR$pathway == tcr_pathwayname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")




list() -> cd4tcr_list 
for (i in 1:length(names(gene_set_fullgene_TCR)) ) {

names(gene_set_fullgene_TCR)[i] -> tcr_pathwayname

assign( paste0("CD4tcr_", i), plotEnrichment(gene_set_fullgene_TCR[[tcr_pathwayname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title =  gsub("_", " ", tcr_pathwayname)  ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf, y = round( fgsea_results_CD4T_TCR$ES[fgsea_results_CD4T_TCR$pathway == tcr_pathwayname ] ,3 )* 0.9 , 
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T_TCR$padj[fgsea_results_CD4T_TCR$pathway == tcr_pathwayname ] ,4 ) ,
                                                                               round( fgsea_results_CD4T_TCR$NES[fgsea_results_CD4T_TCR$pathway == tcr_pathwayname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")  )

cd4tcr_list[[i]] <- get(paste0("CD4tcr_", i))


}




### Fig4G - CD8T ###################
cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/Fig4G_GSEA_CD4T_tcr.pdf" , width = 30 , height = 20)

grid.arrange(   CD4tcr_1, CD4tcr_2, CD4tcr_3, CD4tcr_4, CD4tcr_5, 
                CD4tcr_6, CD4tcr_7, CD4tcr_8, CD4tcr_9, CD4tcr_10, 
                CD4tcr_11, CD4tcr_12, CD4tcr_13, CD4tcr_14, CD4tcr_15, 
                CD4tcr_16, CD4tcr_17, CD4tcr_18, CD4tcr_19, CD4tcr_20 , ncol = 4 )
dev.off()











# Correlation Whole Gene in CD8T #################
sc_cor_frame[3] -> genelist_tmp
genelist_tmp$CD8T_Cor -> geneList
rownames(genelist_tmp) -> names(geneList)
sort(geneList, decreasing = T) -> geneList
length(geneList)


fgsea_results_CD8T_TCR <- fgsea(pathways = gene_set_fullgene_TCR , stats = geneList , nperm = 10000 )

# 결과를 데이터프레임으로 변환합니다.
fgsea_results_CD8T_TCR <- as.data.frame(fgsea_results_CD8T_TCR)


ranked_list <- geneList



ranked_df <- data.frame(
  Gene = names(ranked_list),
  Correlation = ranked_list
)

# Bar plot 생성
bar_plot <- ggplot(ranked_df, aes(x = seq_along(Correlation), y = 1, fill = Correlation)) +
  geom_tile(height = 0.1) +  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(y = "")+
  theme(axis.title.x = element_blank(),
        plot.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 6, colour = "white"),
        
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        plot.margin = unit(c(0, 1, 0, 1), "lines"))


# scatter_plot <- ggplot(ranked_df, aes(x = seq_along(Correlation), y = Correlation)) +
#   geom_point(size = 0.1) +
#   theme_minimal() +
#   labs(x = "Rank", y = "Correlation")+
#   theme(plot.margin = unit(c(0, 1, 1, 1), "lines"))


scatter_plot <- ggplot(ranked_df, aes(x = seq_along(Correlation), y = Correlation)) +
  geom_hline(yintercept = c(-0.5,0,  0.5), color = "#EEEEEE", linetype = "solid") +
  geom_vline(xintercept = c(0,10000,20000), color = "#EEEEEE", linetype = "solid") +
  
  geom_bar(stat = "identity", size = 0.1, color = "lightgray") +
  # geom_text(aes(x = 17000, y = 0, label = "Zero"), size = 3, color = "lightgray", family = "Arial") +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.5)) +
  scale_x_continuous(breaks = seq(0, 30000, by = 10000)) +
  labs(x = "Rank", y = "Correlation") +
  
  theme_void()+
  theme(axis.text.y = element_text(size = 8,  family = "arial", margin = margin(0,4,0,0 ) ),
        axis.text.x = element_text(size = 8,  family = "arial"),
        axis.title.x = element_text(size = 11, margin = margin(3,1,1,1), family = "arial" , vjust = 0.5),
        axis.title.y = element_text(size = 11, margin = margin(1,5,1,1), angle = 90, family = "arial" , hjust = 0.5 ),
        
        
        plot.margin = unit(c(0, 1, 1, 1), "lines"))






# plotEnrichment(gene_set[["PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer(CD8T)"]], geneList) +
#   labs(title = "GSEA Plot: PD L1 expression and PD 1 checkpoint pathway in cancer(CD8T)")


# cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/CEBPBcorGSEAPlot/Kegg_CD8T_PDL1.pdf" , width = 6, height = 4)

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[["PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer"]])

# Enrichment plot 생성
enrichment_plot <- plotEnrichment(gene_set_fullgene[["PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer"]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "GSEA Plot: PD L1 expression and PD 1 checkpoint pathway in cancer" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 11),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 8),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) 

grid.arrange(enrichment_plot, bar_plot, scatter_plot, ncol = 1, heights = c(3, 0.5, 2))


# dev.off()


names(gene_set_fullgene_TCR)[1] -> tcr_pathwayname

plotEnrichment(gene_set_fullgene_TCR[[tcr_pathwayname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title =  gsub("_", " ", tcr_pathwayname)  ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf, y = round( fgsea_results_CD8T_TCR$ES[fgsea_results_CD8T_TCR$pathway == tcr_pathwayname ] ,3 )* 0.9 , 
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD8T_TCR$padj[fgsea_results_CD8T_TCR$pathway == tcr_pathwayname ] ,4 ) ,
                                                                               round( fgsea_results_CD8T_TCR$NES[fgsea_results_CD8T_TCR$pathway == tcr_pathwayname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")




list() -> cd8tcr_list 
for (i in 1:length(names(gene_set_fullgene_TCR)) ) {
  
  names(gene_set_fullgene_TCR)[i] -> tcr_pathwayname
  
  assign( paste0("CD8tcr_", i), plotEnrichment(gene_set_fullgene_TCR[[tcr_pathwayname]], ranked_list) +
            scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
            
            labs(title =  gsub("_", " ", tcr_pathwayname)  ,y= "Enrichment Score") + 
            theme(plot.title = element_text(hjust = 0.5, size = 18),
                  axis.title.y = element_text(size= 14),
                  axis.title.x = element_blank(),        
                  axis.text.x = element_blank(),
                  axis.text.y = element_text(size= 12),
                  axis.ticks.x = element_blank(),
                  axis.ticks.y = element_blank(),
                  plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf, y = round( fgsea_results_CD8T_TCR$ES[fgsea_results_CD8T_TCR$pathway == tcr_pathwayname ] ,3 )* 0.9 , 
                                                                         label = sprintf("FDR : %s\nNES : %s", 
                                                                                         round( fgsea_results_CD8T_TCR$padj[fgsea_results_CD8T_TCR$pathway == tcr_pathwayname ] ,4 ) ,
                                                                                         round( fgsea_results_CD8T_TCR$NES[fgsea_results_CD8T_TCR$pathway == tcr_pathwayname ] ,3 )),
                                                                         hjust = 1, vjust = 1, size = 5, color = "black") )
  
  cd8tcr_list[[i]] <- get(paste0("CD8tcr_", i))
  
  
}


names(gene_set_fullgene_TCR)[3] -> tcr_pathwayname

CD8tcr_3 <- plotEnrichment(gene_set_fullgene_TCR[[tcr_pathwayname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title =  gsub("_", " ", tcr_pathwayname)  ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf, y = 0.18  , 
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD8T_TCR$padj[fgsea_results_CD8T_TCR$pathway == tcr_pathwayname ] ,4 ) ,
                                                                               round( fgsea_results_CD8T_TCR$NES[fgsea_results_CD8T_TCR$pathway == tcr_pathwayname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")



### Fig4G - CD8T ###################
cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/Fig4G_GSEA_CD8T_tcr.pdf" , width = 30 , height = 20)

grid.arrange(   CD8tcr_1, CD8tcr_2, CD8tcr_3, CD8tcr_4, CD8tcr_5, 
                CD8tcr_6, CD8tcr_7, CD8tcr_8, CD8tcr_9, CD8tcr_10, 
                CD8tcr_11, CD8tcr_12, CD8tcr_13, CD8tcr_14, CD8tcr_15, 
                CD8tcr_16, CD8tcr_17, CD8tcr_18, CD8tcr_19, CD8tcr_20 , ncol = 4 )
dev.off()













































# Fig4G Sig ###################
## ################
### Correlation Sig Gene in CD4T #################
sc_cor_frame[ rownames(sc_cor_frame)[rownames(sc_cor_frame) %in% c(CD4T_posgene, CD4T_neggene ) ]  , ][2] -> genelist_tmp
genelist_tmp$CD4T_Cor -> geneList
rownames(genelist_tmp) -> names(geneList)
sort(geneList, decreasing = T) -> geneList
length(geneList)



# fgsea_results <- fgsea(pathways = gene_set_fullgene, stats = geneList, gseaParam = 100000 )
fgsea_results <- fgsea(pathways = gene_set_fullgene, stats = geneList )

fgsea_results_CD4T <- fgsea_results

# 결과를 데이터프레임으로 변환합니다.
fgsea_results_df <- as.data.frame(fgsea_results)
fgsea_results_df

ranked_list <- geneList


ranked_df <- data.frame(
  Gene = names(ranked_list),
  Correlation = ranked_list
)



length(geneList)

# Bar plot 생성
bar_plot <- ggplot(ranked_df, aes(x = seq_along(Correlation), y = 1, fill = Correlation)) +
  geom_tile(height = 0.1) +  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(y = "")+
  theme(axis.title.x = element_blank(),
        plot.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 9, colour = "white" , margin = margin(0,3,0,3)),
        
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        plot.margin = unit(c(0, 1, 0, 1), "lines"))


scatter_plot <- ggplot(ranked_df, aes(x = seq_along(Correlation), y = Correlation)) +
  geom_hline(yintercept = c(-0.5,0,  0.5), color = "#EEEEEE", linetype = "solid") +
  geom_vline(xintercept = c(0,2000), color = "#EEEEEE", linetype = "solid") +
  
  geom_bar(stat = "identity", size = 0.1, color = "lightgray") +
  # geom_text(aes(x = 17000, y = 0, label = "Zero"), size = 3, color = "lightgray", family = "Arial") +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.5)) +
  scale_x_continuous(breaks = seq(0, length(geneList), by = 2000)) +
  labs(x = "Rank", y = "Correlation") +
  
  theme_void()+
  theme(axis.text.y = element_text(size = 12,  family = "arial", margin = margin(0,4,0,2 ) ),
        axis.text.x = element_text(size = 12,  family = "arial"),
        axis.title.x = element_text(size = 14, margin = margin(3,1,1,1), family = "arial" , vjust = 0.5),
        axis.title.y = element_text(size = 14, margin = margin(1,5,1,1), angle = 90, family = "arial" , hjust = 0.5 ),
        
        plot.margin = unit(c(0, 1, 1, 1), "lines"))



# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[["PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer"]])

# Enrichment plot 생성
enrichment_plot_cd4ta <- plotEnrichment(gene_set_fullgene[["PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer"]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "PD L1 expression and PD 1 checkpoint\npathway in cancer" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == "PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer" ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == "PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer" ] ,6 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == "PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer" ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")



# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[["Treg_Cell_Differentiation"]])

# Enrichment plot 생성
enrichment_plot_cd4tb <- plotEnrichment(gene_set_fullgene[["Treg_Cell_Differentiation"]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Treg Cell Differentiation\n" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == "Treg_Cell_Differentiation" ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == "Treg_Cell_Differentiation" ] ,5 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == "Treg_Cell_Differentiation" ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")


grid.arrange(enrichment_plot_cd4ta, enrichment_plot_cd4tb,  bar_plot, bar_plot,  scatter_plot, scatter_plot, ncol = 2, heights = c(2.5, 0.3, 1.5))


#### Fig4G  - CD4T###################

cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/figure4G_GSEA_CD4T.pdf" , width = 10 , height = 4)

grid.arrange(enrichment_plot_cd4ta, enrichment_plot_cd4tb,  bar_plot, bar_plot,  scatter_plot, scatter_plot, ncol = 2, heights = c(2.5, 0.3, 1.5) )


dev.off()



if (p_value < 0.01) {
  # 0.01 미만이면 지수 표기법 사용
  formatted_p_value <- formatC(p_value, format = "e", digits = 1)
  parts <- strsplit(formatted_p_value, "e")[[1]]
  mantissa <- parts[1]
  exponent <- as.numeric(parts[2])
  x_text <- sprintf("<i>P</i> = %s × 10<sup>%d</sup>", mantissa, exponent)
} else if (p_value < 0.1) {
  # 0.010 미만이면 지수 표기법 사용
  x_text <- paste("<i>P</i> =", formatC(p_value, format = "f", digits = 3))        
} else if (p_value == 1) {
  # p_value가 정확히 1인 경우에는 1.0으로 표시
  x_text <-  paste("<i>P</i> = 1.0")
} else {
  # 0.01 이상인 경우 항상 두 자리로 표시
  x_text <- paste("<i>P</i> =", formatC(p_value, format = "f", digits = 2))
}















gene_set_fullgene$Treg_Cells_Promote_Immunosuppression_in_Cancer_Immune_Escape
gene_set_fullgene$T_Cell_Receptor_Signaling
gene_set_fullgene$Effector_T_cell_Inactivation_in_Cancer_Immune_Escape
gene_set_fullgene$CTLA4_Inhibitory_Signaling
gene_set_fullgene$Regulation_Of_Regulatory_T_Cell_Differentiation
gene_set_fullgene$Positive_Regulation_Of_Programmed_Cell_Death


termname <- "Treg_Cells_Promote_Immunosuppression_in_Cancer_Immune_Escape"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd4t_11 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Treg Cells Promote Immunosuppression\nin Cancer Immune Escape" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,5 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")




termname <- "T_Cell_Receptor_Signaling"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd4t_12 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "T Cell Receptor Signaling\n" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,4 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")



termname <- "Effector_T_cell_Inactivation_in_Cancer_Immune_Escape"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd4t_21 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Effector T cell Inactivation\nin Cancer Immune Escape" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,3 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")







termname <- "CTLA4_Inhibitory_Signaling"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd4t_22 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "CTLA4 Inhibitory Signaling\n" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,3 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")



termname <- "Regulation_Of_Regulatory_T_Cell_Differentiation"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd4t_31 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Regulation Of Regulatory T Cell \nDifferentiation" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,2 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")




termname <- "Positive_Regulation_Of_Programmed_Cell_Death"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd4t_32 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Positive Regulation Of Programmed \nCell Death" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,3 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")









#### Fig4Gsup  - CD4T###################

cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/figure4Gsup1_GSEA_CD4T.pdf" , width = 10 , height = 4)

grid.arrange(enrichment_plot_cd4t_11, enrichment_plot_cd4t_12,  bar_plot, bar_plot,  scatter_plot, scatter_plot, ncol = 2, heights = c(2.5, 0.3, 1.5) )


dev.off()


cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/figure4Gsup2_GSEA_CD4T.pdf" , width = 10 , height = 4)

grid.arrange(enrichment_plot_cd4t_21, enrichment_plot_cd4t_22,  bar_plot, bar_plot,  scatter_plot, scatter_plot, ncol = 2, heights = c(2.5, 0.3, 1.5) )


dev.off()


cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/figure4Gsup3_GSEA_CD4T.pdf" , width = 10 , height = 4)

grid.arrange(enrichment_plot_cd4t_31, enrichment_plot_cd4t_32,  bar_plot, bar_plot,  scatter_plot, scatter_plot, ncol = 2, heights = c(2.5, 0.3, 1.5) )


dev.off()









### Correlation Sig Gene in CD8T #################
sc_cor_frame[ rownames(sc_cor_frame)[rownames(sc_cor_frame) %in% c(CD4T_posgene, CD4T_neggene ) ]  , ][3] -> genelist_tmp
genelist_tmp$CD8T_Cor -> geneList
rownames(genelist_tmp) -> names(geneList)
sort(geneList, decreasing = T) -> geneList
length(geneList)


fgsea_results <- fgsea(pathways = gene_set_fullgene, stats = geneList , nperm = 10000 )

fgsea_results_CD8T <- fgsea_results
fgsea_results_CD8T$padj

# 결과를 데이터프레임으로 변환합니다.
fgsea_results_df <- as.data.frame(fgsea_results_CD8T)

ranked_list <- geneList

ranked_df <- data.frame(
  Gene = names(ranked_list),
  Correlation = ranked_list
)

# Bar plot 생성
bar_plot <- ggplot(ranked_df, aes(x = seq_along(Correlation), y = 1, fill = Correlation)) +
  geom_tile(height = 0.1) +  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(y = "")+
  theme(axis.title.x = element_blank(),
        plot.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 9, colour = "white" , margin = margin(0,3,0,3)),
        
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        plot.margin = unit(c(0, 1, 0, 1), "lines"))


scatter_plot <- ggplot(ranked_df, aes(x = seq_along(Correlation), y = Correlation)) +
  geom_hline(yintercept = c(-0.5,0,  0.5), color = "#EEEEEE", linetype = "solid") +
  geom_vline(xintercept = c(0,2000), color = "#EEEEEE", linetype = "solid") +
  
  geom_bar(stat = "identity", size = 0.1, color = "lightgray") +
  # geom_text(aes(x = 17000, y = 0, label = "Zero"), size = 3, color = "lightgray", family = "Arial") +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.5)) +
  scale_x_continuous(breaks = seq(0, length(geneList), by = 2000)) +
  labs(x = "Rank", y = "Correlation") +
  
  theme_void()+
  theme(axis.text.y = element_text(size = 12,  family = "arial", margin = margin(0,4,0,2 ) ),
        axis.text.x = element_text(size = 12,  family = "arial"),
        axis.title.x = element_text(size = 14, margin = margin(3,1,1,1), family = "arial" , vjust = 0.5),
        axis.title.y = element_text(size = 14, margin = margin(1,5,1,1), angle = 90, family = "arial" , hjust = 0.5 ),
        
        plot.margin = unit(c(0, 1, 1, 1), "lines"))





# plotEnrichment(gene_set[["PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer(CD8T)"]], geneList) +
#   labs(title = "GSEA Plot: PD L1 expression and PD 1 checkpoint pathway in cancer(CD8T)")


# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[["PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer"]])

# Enrichment plot 생성
enrichment_plot_cd8ta <- plotEnrichment(gene_set_fullgene[["PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer"]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "PD L1 expression and PD 1 checkpoint\npathway in cancer" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD8T$ES[fgsea_results_CD8T$pathway == "PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer" ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD8T$padj[fgsea_results_CD8T$pathway == "PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer" ] ,4 ) ,
                                                                               round( fgsea_results_CD8T$NES[fgsea_results_CD8T$pathway == "PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer" ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")




# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[["Effector_T_cell_Inactivation_in_Cancer_Immune_Escape"]])

# Enrichment plot 생성
enrichment_plot_cd8tb <- plotEnrichment(gene_set_fullgene[["Effector_T_cell_Inactivation_in_Cancer_Immune_Escape"]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Effector T cell Inactivation\nin Cancer Immune Escape" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf, y = round( fgsea_results_CD8T$ES[fgsea_results_CD8T$pathway == "Effector_T_cell_Inactivation_in_Cancer_Immune_Escape" ] ,3 )* 0.9 , 
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD8T$padj[fgsea_results_CD8T$pathway == "Effector_T_cell_Inactivation_in_Cancer_Immune_Escape" ] ,4 ) ,
                                                                               round( fgsea_results_CD8T$NES[fgsea_results_CD8T$pathway == "Effector_T_cell_Inactivation_in_Cancer_Immune_Escape" ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")



#### Fig4G - CD8T ###################
cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/figure4G_GSEA_CD8T.pdf" , width = 10 , height = 4)

grid.arrange(enrichment_plot_cd8ta, enrichment_plot_cd8tb,  bar_plot, bar_plot,  scatter_plot, scatter_plot, ncol = 2, heights = c(2.5, 0.3, 1.5))


dev.off()






gene_set_fullgene$Treg_Cells_Promote_Immunosuppression_in_Cancer_Immune_Escape
gene_set_fullgene$T_Cell_Receptor_Signaling
gene_set_fullgene$Effector_T_cell_Inactivation_in_Cancer_Immune_Escape
gene_set_fullgene$CTLA4_Inhibitory_Signaling
gene_set_fullgene$Regulation_Of_Regulatory_T_Cell_Differentiation
gene_set_fullgene$Positive_Regulation_Of_Programmed_Cell_Death


termname <- "Treg_Cells_Promote_Immunosuppression_in_Cancer_Immune_Escape"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd8t_11 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Treg Cells Promote Immunosuppression\nin Cancer Immune Escape" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,4 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")




termname <- "T_Cell_Receptor_Signaling"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd8t_12 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "T Cell Receptor Signaling" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,4 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")



termname <- "Effector_T_cell_Inactivation_in_Cancer_Immune_Escape"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd8t_21 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Effector T cell Inactivation in Cancer Immune Escape" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,4 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")







termname <- "CTLA4_Inhibitory_Signaling"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd8t_22 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "CTLA4 Inhibitory Signaling" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,4 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")



termname <- "Regulation_Of_Regulatory_T_Cell_Differentiation"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd8t_31 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Regulation_Of_Regulatory_T_Cell_Differentiation" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,4 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")




termname <- "Positive_Regulation_Of_Programmed_Cell_Death"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd8t_32 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Positive Regulation Of Programmed Cell Death" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,4 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")









#### Fig4Gsup  - CD8T###################

cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/figure4Gsup1_GSEA_CD8T.pdf" , width = 10 , height = 4)

grid.arrange(enrichment_plot_cd8t_11, enrichment_plot_cd8t_12,  bar_plot, bar_plot,  scatter_plot, scatter_plot, ncol = 2, heights = c(2.5, 0.3, 1.5) )


dev.off()


cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/figure4Gsup2_GSEA_CD8T.pdf" , width = 10 , height = 4)

grid.arrange(enrichment_plot_cd8t_21, enrichment_plot_cd8t_22,  bar_plot, bar_plot,  scatter_plot, scatter_plot, ncol = 2, heights = c(2.5, 0.3, 1.5) )


dev.off()


cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/figure4Gsup3_GSEA_CD8T.pdf" , width = 10 , height = 4)

grid.arrange(enrichment_plot_cd8t_31, enrichment_plot_cd8t_32,  bar_plot, bar_plot,  scatter_plot, scatter_plot, ncol = 2, heights = c(2.5, 0.3, 1.5) )


dev.off()










# X Fig4G ZeroCorGenedrop ###################
## ################
### Correlation ZeroCorGenedrop Gene in CD4T #################
sc_cor_frame[ rownames(sc_cor_frame)[rownames(sc_cor_frame) %in% c(CD4T_posgene, CD4T_neggene ) ]  , ][2] -> genelist_tmp
subset(sc_cor_frame, CD4T_Cor != 0) -> genelist_tmp

genelist_tmp$CD4T_Cor -> geneList
rownames(genelist_tmp) -> names(geneList)
sort(geneList, decreasing = T) -> geneList
length(geneList)
table(geneList)[table(geneList) > 1]
sum(table(geneList)[table(geneList) > 1])

# fgsea_results <- fgsea(pathways = gene_set_fullgene, stats = geneList, gseaParam = 100000 )
fgsea_results <- fgsea(pathways = gene_set_fullgene, stats = geneList )

fgsea_results_CD4T <- fgsea_results

# 결과를 데이터프레임으로 변환합니다.
fgsea_results_df <- as.data.frame(fgsea_results)
fgsea_results_df

ranked_list <- geneList


ranked_df <- data.frame(
  Gene = names(ranked_list),
  Correlation = ranked_list
)



length(geneList)

# Bar plot 생성
bar_plot <- ggplot(ranked_df, aes(x = seq_along(Correlation), y = 1, fill = Correlation)) +
  geom_tile(height = 0.1) +  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(y = "")+
  theme(axis.title.x = element_blank(),
        plot.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 9, colour = "white" , margin = margin(0,3,0,3)),
        
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        plot.margin = unit(c(0, 1, 0, 1), "lines"))


scatter_plot <- ggplot(ranked_df, aes(x = seq_along(Correlation), y = Correlation)) +
  geom_hline(yintercept = c(-0.5,0,  0.5), color = "#EEEEEE", linetype = "solid") +
  geom_vline(xintercept = c(0,2000), color = "#EEEEEE", linetype = "solid") +
  
  geom_bar(stat = "identity", size = 0.1, color = "lightgray") +
  # geom_text(aes(x = 17000, y = 0, label = "Zero"), size = 3, color = "lightgray", family = "Arial") +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.5)) +
  scale_x_continuous(breaks = seq(0, length(geneList), by = 2000)) +
  labs(x = "Rank", y = "Correlation") +
  
  theme_void()+
  theme(axis.text.y = element_text(size = 12,  family = "arial", margin = margin(0,4,0,2 ) ),
        axis.text.x = element_text(size = 12,  family = "arial"),
        axis.title.x = element_text(size = 14, margin = margin(3,1,1,1), family = "arial" , vjust = 0.5),
        axis.title.y = element_text(size = 14, margin = margin(1,5,1,1), angle = 90, family = "arial" , hjust = 0.5 ),
        
        plot.margin = unit(c(0, 1, 1, 1), "lines"))



# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[["PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer"]])

# Enrichment plot 생성
enrichment_plot_cd4ta <- plotEnrichment(gene_set_fullgene[["PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer"]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "PD L1 expression and PD 1 checkpoint\npathway in cancer" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == "PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer" ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == "PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer" ] ,6 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == "PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer" ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")



# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[["Treg_Cell_Differentiation"]])

# Enrichment plot 생성
enrichment_plot_cd4tb <- plotEnrichment(gene_set_fullgene[["Treg_Cell_Differentiation"]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Treg Cell Differentiation\n" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == "Treg_Cell_Differentiation" ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == "Treg_Cell_Differentiation" ] ,5 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == "Treg_Cell_Differentiation" ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")


grid.arrange(enrichment_plot_cd4ta, enrichment_plot_cd4tb,  bar_plot, bar_plot,  scatter_plot, scatter_plot, ncol = 2, heights = c(2.5, 0.3, 1.5))


#### Fig4G  - CD4T###################

cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/figure4G_GSEA_CD4T.pdf" , width = 10 , height = 4)

grid.arrange(enrichment_plot_cd4ta, enrichment_plot_cd4tb,  bar_plot, bar_plot,  scatter_plot, scatter_plot, ncol = 2, heights = c(2.5, 0.3, 1.5) )


dev.off()






gene_set_fullgene$Treg_Cells_Promote_Immunosuppression_in_Cancer_Immune_Escape
gene_set_fullgene$T_Cell_Receptor_Signaling
gene_set_fullgene$Effector_T_cell_Inactivation_in_Cancer_Immune_Escape
gene_set_fullgene$CTLA4_Inhibitory_Signaling
gene_set_fullgene$Regulation_Of_Regulatory_T_Cell_Differentiation
gene_set_fullgene$Positive_Regulation_Of_Programmed_Cell_Death


termname <- "Treg_Cells_Promote_Immunosuppression_in_Cancer_Immune_Escape"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd4t_11 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Treg Cells Promote Immunosuppression\nin Cancer Immune Escape" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,5 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")




termname <- "T_Cell_Receptor_Signaling"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd4t_12 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "T Cell Receptor Signaling\n" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,4 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")



termname <- "Effector_T_cell_Inactivation_in_Cancer_Immune_Escape"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd4t_21 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Effector T cell Inactivation\nin Cancer Immune Escape" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,3 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")







termname <- "CTLA4_Inhibitory_Signaling"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd4t_22 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "CTLA4 Inhibitory Signaling\n" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,3 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")



termname <- "Regulation_Of_Regulatory_T_Cell_Differentiation"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd4t_31 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Regulation Of Regulatory T Cell \nDifferentiation" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,2 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")




termname <- "Positive_Regulation_Of_Programmed_Cell_Death"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd4t_32 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Positive Regulation Of Programmed \nCell Death" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,3 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")









#### Fig4Gsup  - CD4T###################

cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/figure4Gsup1_GSEA_CD4T.pdf" , width = 10 , height = 4)

grid.arrange(enrichment_plot_cd4t_11, enrichment_plot_cd4t_12,  bar_plot, bar_plot,  scatter_plot, scatter_plot, ncol = 2, heights = c(2.5, 0.3, 1.5) )


dev.off()


cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/figure4Gsup2_GSEA_CD4T.pdf" , width = 10 , height = 4)

grid.arrange(enrichment_plot_cd4t_21, enrichment_plot_cd4t_22,  bar_plot, bar_plot,  scatter_plot, scatter_plot, ncol = 2, heights = c(2.5, 0.3, 1.5) )


dev.off()


cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/figure4Gsup3_GSEA_CD4T.pdf" , width = 10 , height = 4)

grid.arrange(enrichment_plot_cd4t_31, enrichment_plot_cd4t_32,  bar_plot, bar_plot,  scatter_plot, scatter_plot, ncol = 2, heights = c(2.5, 0.3, 1.5) )


dev.off()









### Correlation ZeroCorGenedrop Gene in CD8T #################
sc_cor_frame[ rownames(sc_cor_frame)[rownames(sc_cor_frame) %in% c(CD4T_posgene, CD4T_neggene ) ]  , ][3] -> genelist_tmp
genelist_tmp$CD8T_Cor -> geneList
rownames(genelist_tmp) -> names(geneList)
sort(geneList, decreasing = T) -> geneList
length(geneList)


fgsea_results <- fgsea(pathways = gene_set_fullgene, stats = geneList , nperm = 10000 )

fgsea_results_CD8T <- fgsea_results
fgsea_results_CD8T$padj

# 결과를 데이터프레임으로 변환합니다.
fgsea_results_df <- as.data.frame(fgsea_results_CD8T)

ranked_list <- geneList

ranked_df <- data.frame(
  Gene = names(ranked_list),
  Correlation = ranked_list
)

# Bar plot 생성
bar_plot <- ggplot(ranked_df, aes(x = seq_along(Correlation), y = 1, fill = Correlation)) +
  geom_tile(height = 0.1) +  
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  labs(y = "")+
  theme(axis.title.x = element_blank(),
        plot.title = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 9, colour = "white" , margin = margin(0,3,0,3)),
        
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        panel.grid = element_blank(),
        plot.margin = unit(c(0, 1, 0, 1), "lines"))


scatter_plot <- ggplot(ranked_df, aes(x = seq_along(Correlation), y = Correlation)) +
  geom_hline(yintercept = c(-0.5,0,  0.5), color = "#EEEEEE", linetype = "solid") +
  geom_vline(xintercept = c(0,2000), color = "#EEEEEE", linetype = "solid") +
  
  geom_bar(stat = "identity", size = 0.1, color = "lightgray") +
  # geom_text(aes(x = 17000, y = 0, label = "Zero"), size = 3, color = "lightgray", family = "Arial") +
  scale_y_continuous(breaks = seq(-1, 1, by = 0.5)) +
  scale_x_continuous(breaks = seq(0, length(geneList), by = 2000)) +
  labs(x = "Rank", y = "Correlation") +
  
  theme_void()+
  theme(axis.text.y = element_text(size = 12,  family = "arial", margin = margin(0,4,0,2 ) ),
        axis.text.x = element_text(size = 12,  family = "arial"),
        axis.title.x = element_text(size = 14, margin = margin(3,1,1,1), family = "arial" , vjust = 0.5),
        axis.title.y = element_text(size = 14, margin = margin(1,5,1,1), angle = 90, family = "arial" , hjust = 0.5 ),
        
        plot.margin = unit(c(0, 1, 1, 1), "lines"))





# plotEnrichment(gene_set[["PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer(CD8T)"]], geneList) +
#   labs(title = "GSEA Plot: PD L1 expression and PD 1 checkpoint pathway in cancer(CD8T)")


# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[["PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer"]])

# Enrichment plot 생성
enrichment_plot_cd8ta <- plotEnrichment(gene_set_fullgene[["PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer"]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "PD L1 expression and PD 1 checkpoint\npathway in cancer" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD8T$ES[fgsea_results_CD8T$pathway == "PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer" ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD8T$padj[fgsea_results_CD8T$pathway == "PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer" ] ,4 ) ,
                                                                               round( fgsea_results_CD8T$NES[fgsea_results_CD8T$pathway == "PD_L1_expression_and_PD_1_checkpoint_pathway_in_cancer" ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")




# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[["Effector_T_cell_Inactivation_in_Cancer_Immune_Escape"]])

# Enrichment plot 생성
enrichment_plot_cd8tb <- plotEnrichment(gene_set_fullgene[["Effector_T_cell_Inactivation_in_Cancer_Immune_Escape"]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Effector T cell Inactivation\nin Cancer Immune Escape" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf, y = round( fgsea_results_CD8T$ES[fgsea_results_CD8T$pathway == "Effector_T_cell_Inactivation_in_Cancer_Immune_Escape" ] ,3 )* 0.9 , 
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD8T$padj[fgsea_results_CD8T$pathway == "Effector_T_cell_Inactivation_in_Cancer_Immune_Escape" ] ,4 ) ,
                                                                               round( fgsea_results_CD8T$NES[fgsea_results_CD8T$pathway == "Effector_T_cell_Inactivation_in_Cancer_Immune_Escape" ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")



#### Fig4G - CD8T ###################
cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/figure4G_GSEA_CD8T.pdf" , width = 10 , height = 4)

grid.arrange(enrichment_plot_cd8ta, enrichment_plot_cd8tb,  bar_plot, bar_plot,  scatter_plot, scatter_plot, ncol = 2, heights = c(2.5, 0.3, 1.5))


dev.off()






gene_set_fullgene$Treg_Cells_Promote_Immunosuppression_in_Cancer_Immune_Escape
gene_set_fullgene$T_Cell_Receptor_Signaling
gene_set_fullgene$Effector_T_cell_Inactivation_in_Cancer_Immune_Escape
gene_set_fullgene$CTLA4_Inhibitory_Signaling
gene_set_fullgene$Regulation_Of_Regulatory_T_Cell_Differentiation
gene_set_fullgene$Positive_Regulation_Of_Programmed_Cell_Death


termname <- "Treg_Cells_Promote_Immunosuppression_in_Cancer_Immune_Escape"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd8t_11 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Treg Cells Promote Immunosuppression\nin Cancer Immune Escape" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,4 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")




termname <- "T_Cell_Receptor_Signaling"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd8t_12 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "T Cell Receptor Signaling" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,4 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")



termname <- "Effector_T_cell_Inactivation_in_Cancer_Immune_Escape"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd8t_21 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Effector T cell Inactivation in Cancer Immune Escape" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,4 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")







termname <- "CTLA4_Inhibitory_Signaling"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd8t_22 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "CTLA4 Inhibitory Signaling" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,4 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")



termname <- "Regulation_Of_Regulatory_T_Cell_Differentiation"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd8t_31 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Regulation_Of_Regulatory_T_Cell_Differentiation" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,4 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")




termname <- "Positive_Regulation_Of_Programmed_Cell_Death"

# 각 유전자의 위치 찾기
hits <- which(names(ranked_list) %in% gene_set_fullgene[[termname]])

# Enrichment plot 생성
enrichment_plot_cd8t_32 <- plotEnrichment(gene_set_fullgene[[termname]], ranked_list) +
  scale_y_continuous(breaks = seq(0, 0.9, by = 0.25)) +
  
  labs(title = "Positive Regulation Of Programmed Cell Death" ,y= "Enrichment Score") + 
  theme(plot.title = element_text(hjust = 0.5, size = 18),
        axis.title.y = element_text(size= 14),
        axis.title.x = element_blank(),        
        axis.text.x = element_blank(),
        axis.text.y = element_text(size= 12),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin = unit(c(1, 1, 0, 1), "lines")) + annotate("text", x = Inf,  y = round( fgsea_results_CD4T$ES[fgsea_results_CD4T$pathway == termname ] ,3 )* 0.9 ,
                                                               label = sprintf("FDR : %s\nNES : %s", 
                                                                               round( fgsea_results_CD4T$padj[fgsea_results_CD4T$pathway == termname ] ,4 ) ,
                                                                               round( fgsea_results_CD4T$NES[fgsea_results_CD4T$pathway == termname ] ,3 )),
                                                               hjust = 1, vjust = 1, size = 5, color = "black")









#### Fig4Gsup  - CD8T###################

cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/figure4Gsup1_GSEA_CD8T.pdf" , width = 10 , height = 4)

grid.arrange(enrichment_plot_cd8t_11, enrichment_plot_cd8t_12,  bar_plot, bar_plot,  scatter_plot, scatter_plot, ncol = 2, heights = c(2.5, 0.3, 1.5) )


dev.off()


cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/figure4Gsup2_GSEA_CD8T.pdf" , width = 10 , height = 4)

grid.arrange(enrichment_plot_cd8t_21, enrichment_plot_cd8t_22,  bar_plot, bar_plot,  scatter_plot, scatter_plot, ncol = 2, heights = c(2.5, 0.3, 1.5) )


dev.off()


cairo_pdf("/home/ch2860/Code/Seurat_CRC/Submit_figure_Code/Fig/figure4/figure4Gsup3_GSEA_CD8T.pdf" , width = 10 , height = 4)

grid.arrange(enrichment_plot_cd8t_31, enrichment_plot_cd8t_32,  bar_plot, bar_plot,  scatter_plot, scatter_plot, ncol = 2, heights = c(2.5, 0.3, 1.5) )


dev.off()
