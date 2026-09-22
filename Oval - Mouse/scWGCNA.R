setwd("/home/Data_Drive_8TB/kykim/scWGCNA/")
renv::init()
source("/home/Data_Drive_8TB/kykim/scWGCNA/renv/activate.R")
renv::project()

####
save.image("my_session.RData")
save(my.small_MmLimbE155, file = "my.small_MmLimbE155.RData")
save(MmLimbE155.scWGCNA, file = "MmLimbE155.scWGCNA.RData")
load("my_session.RData")
my.small_MmLimbE155 <- get(load("my.small_MmLimbE155.RData"))
MmLimbE155.scWGCNA <- get(load("MmLimbE155.scWGCNA.RData"))
####

if(!requireNamespace("devtools", quietly = TRUE))
  install.packages("devtools")
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("escape", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("escape")
}


install.packages("reticulate")
install.packages("Cairo")
install.packages("Seurat")
install.packages("dplyr")
install.packages("tibble")
install.packages("WGCNA")
install.packages("remotes")
install.packages("intergraph")
install.packages("VennDiagram")
install.packages("pheatmap")

BiocManager::install(c("GO.db", "impute", "preprocessCore"))
BiocManager::install('biomaRt')
devtools::install_github("cferegrino/scWGCNA", ref="main")
reticulate::py_install("anndata", pip=TRUE)
remotes::install_version("Seurat", version = "4.3.0")
renv::install("Seurat@4.3.0")

library(Seurat)
library(reticulate)
library(Cairo)
library(devtools)
library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(data.table)
library(tibble)
library(WGCNA)
library(scWGCNA)
library(igraph)
library(network)
library(GGally)
library(intergraph)
library(biomaRt)
library(VennDiagram)
library(pheatmap)
library(matrixStats)

######
remove.packages("Seurat")
packageVersion("Seurat")
######
###### Basal -> Luminal Epithelium #####
ad <- import("anndata")
data <- ad$read_h5ad("/home/Data_Drive_8TB/kykim/1. Oval/Monocle/Luminal_adata.h5ad")
cell_names <- py_to_r(data$obs_names$tolist())
gene_names <- py_to_r(data$var_names$tolist())
expr_matrix <- t(as.matrix(data$X))
rownames(expr_matrix) <- gene_names
colnames(expr_matrix) <- cell_names
data <- CreateSeuratObject(counts = expr_matrix, meta.data = data$obs)
data <- subset(data, subset = `clinical information` != "Glandular Epithelium")
data <- NormalizeData(data, normalization.method = "LogNormalize", scale.factor = 10000)

######
DefaultAssay(data) <- "RNA"
data <- ScaleData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 5000)

muc_genes <- grep("Muc", rownames(data), value = TRUE)
vf <- VariableFeatures(data)
vf <- unique(c(vf, muc_genes))
VariableFeatures(data) <- vf

data <- subset(data, features = VariableFeatures(data))
data <- RunPCA(data, features = VariableFeatures(data))

my.small_MmLimbE155 <- data
MmLimbE155.pcells = calculate.pseudocells(s.cells = my.small_MmLimbE155, # Single cells in Seurat object
                                          seeds=0.1, # Fraction of cells to use as seeds to aggregate pseudocells
                                          nn = 10, # Number of neighbors to aggregate
                                          reduction = "pca", # Reduction to use
                                          dims = 1:20) # The dimensions to use

CairoPNG("/home/Data_Drive_8TB/kykim/1. Oval/WGCNA/Scale_Independence.png", width = 1500, height = 1200, res = 200)
MmLimbE155.scWGCNA = run.scWGCNA(p.cells = MmLimbE155.pcells, # Pseudocells (recommended), or Seurat single cells
                                 s.cells = my.small_MmLimbE155, # single cells in Seurat format
                                 is.pseudocell = T, # We are using single cells twice this time
                                 features = VariableFeatures(my.small_MmLimbE155))
dev.off()

## Module - Module correlation ##
MEs <- MmLimbE155.scWGCNA$MEs[, 1:10]
colnames(MEs) <- c("module1","module2","module3","module4","module5","module6","module7","module8","module9","module10")
module_cor <- cor(MEs, use = "pairwise.complete.obs")
CairoPNG("/home/Data_Drive_8TB/kykim/1. Oval/WGCNA/Correlation_heatmap.png", width = 1500, height = 1200, res = 200)
pheatmap(module_cor,
         color = colorRampPalette(c("blue","white","red"))(100),
         main = "Module–Module Correlation",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         treeheight_row = 0,  # 가지 안보이게
         treeheight_col = 0,  # 가지 안보이게
         border_color = NA,
         display_numbers = TRUE,
         number_format = "%.2f")
dev.off()


## Dendroidgram ##
CairoPNG("/home/Data_Drive_8TB/kykim/1. Oval/WGCNA/Total_dendroidgram.png", width = 1500, height = 1200, res = 150)
scW.p.dendro(scWGCNA.data = MmLimbE155.scWGCNA)
dev.off()

CairoPNG("/home/Data_Drive_8TB/kykim/1. Oval/WGCNA/TSNE_plot.png", width = 1500, height = 1200, res = 200)
scW.p.expression(s.cells = my.small_MmLimbE155, # Single cells in Seurat format
                 scWGCNA.data = MmLimbE155.scWGCNA, # scWGCNA list dataset
                 modules = "all", # Which modules to plot?
                 reduction = "pca", # Which reduction to plot?
                 ncol=3) # How many columns to use, in case we're plotting several?
dev.off()

CairoPNG("/home/Data_Drive_8TB/kykim/1. Oval/WGCNA/module4_TSNE.png", width = 1500, height = 1200, res = 150)
scW.p.expression(s.cells = my.small_MmLimbE155, scWGCNA.data = MmLimbE155.scWGCNA, reduction = "pca",  modules = 4)
dev.off()


names(MmLimbE155.scWGCNA$modules)
my.module = MmLimbE155.scWGCNA$modules$`5_greenyellow`
my.gnames = my.small_MmLimbE155@misc$gnames
my.module$gname = my.gnames[rownames(my.module), "Gene.name"]


# Total Network (Each model)
MmLimbE155.scWGCNA = scWGCNA.networks(scWGCNA.data = MmLimbE155.scWGCNA)

CairoPNG("/home/Data_Drive_8TB/kykim/1. Oval/WGCNA/model2_network.png", width = 1500, height = 1200, res = 80)
scW.p.network(MmLimbE155.scWGCNA, module=2)
dev.off()
CairoPNG("/home/Data_Drive_8TB/kykim/1. Oval/WGCNA/model7_network.png", width = 1500, height = 1200, res = 100)
scW.p.network(MmLimbE155.scWGCNA, module=7)
dev.off()

module_genes <- MmLimbE155.scWGCNA$module.genes[[2]]
adj <- MmLimbE155.scWGCNA$TOM[module_genes, module_genes]
g <- graph_from_adjacency_matrix(adj, mode="undirected", weighted=TRUE, diag=FALSE)
selected_genes <- c('Tgfbr2','Smad3','Hif1a','Plau','Egfr','Ceacam1','Epcam',
                        'Foxo1', 'Muc4')
sub_nodes <- intersect(selected_genes, V(g)$name)
g_sub <- induced_subgraph(g, sub_nodes)
CairoPNG("/home/Data_Drive_8TB/kykim/1. Oval/WGCNA/model2_network_gene_pick.png",width = 1500, height = 1200, res = 150)
plot(g_sub, vertex.label = V(g_sub)$name,
     vertex.color = "lightblue",
     vertex.size = 25,
     edge.width = E(g_sub)$weight * 100)  # 굵기 증가
dev.off()

module_genes <- MmLimbE155.scWGCNA$module.genes[[7]]
adj <- MmLimbE155.scWGCNA$TOM[module_genes, module_genes]
g <- graph_from_adjacency_matrix(adj, mode="undirected", weighted=TRUE, diag=FALSE)
sub_nodes <- intersect(selected_genes, V(g)$name)
g_sub <- induced_subgraph(g, sub_nodes)
CairoPNG("/home/Data_Drive_8TB/kykim/1. Oval/WGCNA/model7_network_gene_pick.png",
         width = 1500, height = 1200, res = 150)
plot(g_sub, vertex.label = V(g_sub)$name,
     vertex.color = "lightblue",
     vertex.size = 25,
     edge.width = E(g_sub)$weight * 100)  # 굵기 증가
dev.off()





Gg.ps=calculate.pseudocells(my.small_GgLimbHH29, dims = 1:10)
my.ortho = merge(my.small_MmLimbE155@misc$gnames,my.small_GgLimbHH29@misc$gnames, by = "Gene.name")
my.ortho=my.ortho[,2:3]

### Data setting
expr <- t(MmLimbE155.scWGCNA$expr)
MmLimbE155.scWGCNA$expr <- expr
gene_map <- my.small_MmLimbE155@misc$gnames
common <- intersect(rownames(MmLimbE155.scWGCNA$expr), gene_map$Gene.name)
rownames(MmLimbE155.scWGCNA$expr)[match(common, rownames(MmLimbE155.scWGCNA$expr))] <-  gene_map$Gene.stable.ID[match(common, gene_map$Gene.name)]


MmvGg.comparative = scWGNA.compare(scWGCNA.data = MmLimbE155.scWGCNA,
                                   test.list = list(Gg.ps),
                                   test.names = c("Gg"),
                                   ortho = my.ortho, # not needed unless reference and tests have different gene names
                                   ortho.sp = c(2))


CairoPNG("/home/Data_Drive_8TB/kykim/1. Oval/WGCNA/model_gene_percentage.png", width = 1500, height = 1200, res = 150)
scW.p.modulefrac(MmvGg.comparative)
dev.off()


CairoPNG("/home/Data_Drive_8TB/kykim/1. Oval/WGCNA/Zsummary_median.png", width = 1500, height = 1200, res = 150)
scW.p.preservation(scWGCNA.comp.data = MmvGg.comparative,
                   to.plot=c("preservation", "median.rank"))
dev.off()


CairoPNG("/home/Data_Drive_8TB/kykim/1. Oval/WGCNA/density_conn.png", width = 1500, height = 1200, res = 150)
scW.p.preservation(scWGCNA.comp.data = MmvGg.comparative, to.plot=c("density", "connectivity"))
dev.off()









###############
names(MmLimbE155.scWGCNA$modules)

module1 <- MmLimbE155.scWGCNA$module.genes[[1]]
module2 <- MmLimbE155.scWGCNA$module.genes[[2]]
module3 <- MmLimbE155.scWGCNA$module.genes[[3]]
module4 <- MmLimbE155.scWGCNA$module.genes[[4]]
module5 <- MmLimbE155.scWGCNA$module.genes[[5]]
module6 <- MmLimbE155.scWGCNA$module.genes[[6]]
module7 <- MmLimbE155.scWGCNA$module.genes[[7]]
module8 <- MmLimbE155.scWGCNA$module.genes[[8]]
module9 <- MmLimbE155.scWGCNA$module.genes[[9]]
module10 <- MmLimbE155.scWGCNA$module.genes[[10]]
#module11 <- MmLimbE155.scWGCNA$module.genes[[11]]
#module12 <- MmLimbE155.scWGCNA$module.genes[[12]]
#module13 <- MmLimbE155.scWGCNA$module.genes[[13]]

modules <- list(module1, module2, module3,module4, module5, module6, module7, module8)#, module9, module10)#, module11, module12, module13)
muc4_in_module <- sapply(modules, function(x) "Muc1" %in% x)
which(muc4_in_module)


plot_module_gene_counts <- function(scWGCNA_object, 
                                    main_title = "Gene Counts per Module",
                                    colors = NULL) {
  module_counts <- sapply(scWGCNA_object$module.genes, length)
  module_names <- paste0("Module ", seq_along(module_counts))
  if (is.null(colors)) {
    colors <- rainbow(length(module_counts))
  }
  pie(module_counts,
      labels = paste0(module_names, "\n(n=", module_counts, ")"),
      col = colors,
      main = main_title,
      cex = 0.6)
  legend("topright", 
         legend = module_names,
         fill = colors,
         cex = 1,
         bty = "n",
         title = "Modules")
  result_df <- data.frame(
    Module = module_names,
    Gene_Count = module_counts
  )
  return(invisible(result_df))
}

CairoPNG("/home/Data_Drive_8TB/kykim/1. Oval/WGCNA/pie_plot.png",  width = 1500, height = 1200, res = 150)
plot_module_gene_counts(MmLimbE155.scWGCNA,     main_title = "Gene Counts (Module)")
dev.off()

##############
out_dir <- "/home/Data_Drive_8TB/kykim/1. Oval/WGCNA/"
write.table(module1,file = paste0(out_dir, "module1.txt"),  quote = FALSE,  row.names = FALSE,col.names = FALSE)
write.table(module2,file = paste0(out_dir, "module2.txt"),  quote = FALSE,  row.names = FALSE,col.names = FALSE)
write.table(module3,file = paste0(out_dir, "module3.txt"),  quote = FALSE,  row.names = FALSE,col.names = FALSE)
write.table(module4,file = paste0(out_dir, "module4.txt"),  quote = FALSE,  row.names = FALSE,col.names = FALSE)
write.table(module5,file = paste0(out_dir, "module5.txt"),  quote = FALSE,  row.names = FALSE,col.names = FALSE)
write.table(module6,file = paste0(out_dir, "module6.txt"),  quote = FALSE,  row.names = FALSE,col.names = FALSE)
write.table(module7,file = paste0(out_dir, "module7.txt"),  quote = FALSE,  row.names = FALSE,col.names = FALSE)
write.table(module8,file = paste0(out_dir, "module8.txt"),  quote = FALSE,  row.names = FALSE,col.names = FALSE)
write.table(module9,file = paste0(out_dir, "module9.txt"),  quote = FALSE,  row.names = FALSE,col.names = FALSE)
write.table(module10,file = paste0(out_dir, "module10.txt"),  quote = FALSE,  row.names = FALSE,col.names = FALSE)





###############
module_sizes <- sapply(MmLimbE155.scWGCNA$module.genes, length)
total_genes <- sum(module_sizes)
module_percent <- round(module_sizes / total_genes * 100, 2)

result <- data.frame(
  module = paste0("module", seq_along(module_sizes)),
  gene_count = module_sizes,
  percent = module_percent
)

print(result)
cat("Total genes:", total_genes, "\n")
