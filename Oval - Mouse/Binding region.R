setwd("/home/Data_Drive_8TB/kykim/rtracklayer/")
renv::init()
source("/home/Data_Drive_8TB/kykim/rtracklayer/renv/activate.R")
renv::project()
.libPaths(c("/home/Data_Drive_8TB/kykim/R_libs", .libPaths()))
my_lib <- "/home/Data_Drive_8TB/kykim/R_libs"
.libPaths(c(my_lib, .libPaths()))

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

options(timeout = 3600)

BiocManager::install(c("BiocGenerics", "S4Vectors", "IRanges", "GenomeInfoDb", "GenomicRanges", "Biostrings", "BSgenome"))
BiocManager::install(c("motifmatchr", "JASPAR2022"), lib = my_lib, update = FALSE)
BiocManager::install("TFBSTools")
BiocManager::install("TFMPvalue")
install.packages("/home/Data_Drive_8TB/kykim/BSgenome.Hsapiens.UCSC.hg38_1.4.5.tar.gz", 
                 repos = NULL, 
                 type = "source")
install.packages("/home/Data_Drive_8TB/kykim/TFMPvalue_0.0.9.tar.gz", 
                 repos = NULL, 
                 type = "source",
                 lib = "/home/Data_Drive_8TB/kykim/R_libs")


library(motifmatchr)
library(TFBSTools)
library(JASPAR2022)
library(BSgenome.Hsapiens.UCSC.hg38)
library(rtracklayer)
library(GenomicRanges)
library(TFMPvalue)
library(Gviz)

bed_path <- "/home/Data_Drive_8TB/kykim/1. Oval/Save_data/Human/ATAC/uteru_tissue_female_adult_46_years/ENCFF547ONJ.bed"
genome <- BSgenome.Hsapiens.UCSC.hg38

chr_name <- "chr3"

#### Promoter ####
#plot_start <- 195811613
#plot_end <- 195812135

#### Enhencer ####
plot_start <- 195741765
plot_end <- 195854515

#plot_start <- 195760383
#plot_end <- 195760400

target_region <- GRanges(chr_name, IRanges(plot_start, plot_end))

all_peaks <- import(bed_path, format = "narrowPeak")
target_peaks <- subsetByOverlaps(all_peaks, target_region)

opts <- list(name = "MAX", species = 9606, all_versions = FALSE)
mxd1_pwm <- getMatrixSet(JASPAR2022, opts)

# 모티프 매칭 실행 (p.cutoff 0.0001)
motif_ix <- matchMotifs(mxd1_pwm, target_peaks, genome = genome, out = "matches")
motif_pos <- matchMotifs(mxd1_pwm, target_peaks, genome = genome, out = "positions")

peak_results <- as.data.frame(target_peaks)
peak_results$has_MXD1 <- as.logical(motifMatches(motif_ix))

mxd1_locations <- as.data.frame(unlist(motif_pos))

# 결과 출력
cat("--- Analysis Summary ---\n")
cat("Total peaks in region:", length(target_peaks), "\n")
cat("Peaks containing MXD1 motif:", sum(peak_results$has_MXD1), "\n\n")

if(sum(peak_results$has_MXD1) > 0) {
  cat("--- Detected MXD1 Binding Sites ---\n")
  print(mxd1_locations[, c("seqnames", "start", "end", "strand", "score")])
} else {
  cat("No MXD1 motifs found in the target peaks under current thresholds.\n")
}


#################### 실제 시퀀스 ######################

peak_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, target_peaks)
my_pattern <- "CACGTG"

match_results <- vmatchPattern(my_pattern, peak_seqs)

match_counts <- elementNROWS(match_results)
peak_results_exact <- as.data.frame(target_peaks)
peak_results_exact$sequence_match_count <- match_counts

print(subset(peak_results_exact, sequence_match_count > 0))

################## PLOT #####################

getSeq(Hsapiens, "chr3", start=195747890, end=195747915)

chr <- "chr3"
reg_start <- 195820338
reg_end   <- 195820358

peak_seq <- getSeq(Hsapiens, chr, start=reg_start, end=reg_end)
seq_vec <- strsplit(as.character(peak_seq), "")[[1]]
base_colors <- c(A="#109618", C="#3366CC", G="#FF9900", T="#DC3912")
cols <- base_colors[seq_vec]

sTrack <- AnnotationTrack(start = reg_start:reg_end,
                          width = 0,
                          chromosome = chr,
                          genome = "hg38",
                          name = "Sequence",
                          fill = cols,
                          col = NA,
                          id = seq_vec)

match_pos <- matchPattern("CACGTG", peak_seq)
motif_starts <- reg_start + start(match_pos) - 1
motif_ends   <- reg_start + end(match_pos) - 1
highlight_gr <- GRanges(chr, IRanges(motif_starts, motif_ends))

htrack <- HighlightTrack(trackList = list(sTrack), 
                         range = highlight_gr, 
                         fill = "#FFE0E0", col = "red", lwd = 2)

CairoPNG(filename = "/home/Data_Drive_8TB/kykim/1. Oval/Save_data/Human/ATAC/uteru_tissue_female_adult_46_years/MXD1_Actual_Match_Plot.png", 
         width = 2400, height = 500, res = 200)

plotTracks(list(GenomeAxisTrack(add53 = TRUE), htrack),
           from = reg_start, 
           to = reg_end,
           chromosome = chr,
           main = "MXD1 Transcription factor binding site",
           cex.main = 1.5,
           showId = FALSE,
           showFeatureId = TRUE,
           fontcolor.feature = "white",
           cex.feature = 2,
           cex.title = 1.5,
           featureAnnotation = "id",
           just.feature = "center",
           stacking = "dense")

dev.off()

