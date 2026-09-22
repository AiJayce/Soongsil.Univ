setwd("/home/Data_Drive_8TB/kykim/Chipseq/")
renv::init()
source("/home/Data_Drive_8TB/kykim/monocle/renv/activate.R")
renv::project()

options(verbose = FALSE, tidyverse.quiet = TRUE, warn = -1,
        future.rng.onMisuse = "ignore", future.globals.maxSize = 1e11,
        timeout = 10000)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c(
  "GenomicRanges",
  "IRanges",
  "rtracklayer",
  "ChIPseeker",
  "BSgenome.Hsapiens.UCSC.hg38",
  "Gviz",
  "biomaRt"
))
install.packages("Cairo")

library(GenomicRanges)                          # GRanges, IRanges 등 핵심 자료구조
library(IRanges)
library(rtracklayer)                            # bedGraph import
library(ChIPseeker)                             # covplot(커버리지 플롯)
library(BSgenome.Hsapiens.UCSC.hg38)            # hg38 염색체 길이 정보(Seqinfo)
library(Gviz)                                   # 트랙 기반 시각화
library(biomaRt)                                # 유전자 주석 트랙(BiomartGeneRegionTrack)
library(Cairo)
# ---- 작업 디렉토리 및 데이터 불러오기 -----------------------------------------
setwd("/home/Data_Drive_8TB/kykim/2. CRC/Chip-seq/HCT116/")

# bedGraph → GRanges로 로드 (score 메타컬럼
GSM1917772 <- rtracklayer::import("GSM1917772.bedGraph", format = "bedGraph")
GSM1917773 <- rtracklayer::import("GSM1917773.bedGraph", format = "bedGraph")
#GSM1417250 <- rtracklayer::import("GSM1417250.bedGraph", format = "bedGraph")
# ---- 커버리지 대략 확인 (선택) -------------------------------------------------
# 특정 염색체에서 bedGraph의 score를 가중치로 사용해 커버리지 형태를 빠르게 점검  
# (필요 없으면 주석 처리)
#covplot(GSM1917773, chrs = c("chr20"), weightCol = "score")

# ---- 프래그먼트 길이로 구간 확장-------------------------------
# 주의: 원래는 "리드" 단위 GRanges에 쓰는 패턴
prepareChIPseq = function(reads, frag.len = 200) {
  cat(paste0('Using user-defined fragment size: ', frag.len, "\n"))
  reads.extended = resize(reads, width = frag.len)
  return(trim(reads.extended))
}
GSM1917772 = prepareChIPseq(GSM1917772, frag.len = 200)
GSM1917773 = prepareChIPseq(GSM1917773, frag.len = 200)
#GSM1417250 = prepareChIPseq(GSM1417250, frag.len = 200)
# ---- 타일(빈) 생성: chr20 특정 구간 시각화를 위한 준비 ------------------------
# hg38의 염색체 길이 정보 로드 및 관심 염색체만 선택
si <- seqinfo(BSgenome.Hsapiens.UCSC.hg38)
si <- si[paste0("chr", c(1:22, "X", "Y"))]

# 관심 구간 설정 (CEBPB 위치 근처 예시: chr20:50,180,000-50,200,000)
chrnum <- "chr20"
#startposition <- 50189000
#endposition   <- 50191000

startposition <- 50180000
endposition   <- 50200000

# 일정 크기의 타일(빈)으로 구간을 분할 (여기서는 200bp)
binsize <- 200
bins <- tileGenome(si[chrnum], tilewidth = binsize, cut.last.tile.in.chrom = TRUE)

# ---- 빈 단위로 시그널 집계------------------------------------
# 단순 카운트(빈별 겹침 개수)
BinChIPseq = function(reads, bins) {
  mcols(bins)$score = countOverlaps(bins, reads)
  return(bins)
}

# 라이브러리 크기 보정(RPM) – 샘플 간 비교 용이
BinChIPseq_RPM = function(reads, bins) {
  raw.counts = countOverlaps(bins, reads)
  rpm.counts = raw.counts / sum(raw.counts) * 1e6
  mcols(bins)$score = rpm.counts
  return(bins)
}

# (1) 원시 카운트 스케일 빈
GSM1917772.200bins <- BinChIPseq(GSM1917772, bins)
GSM1917773.200bins <- BinChIPseq(GSM1917773, bins)
#GSM1417250.200bins <- BinChIPseq(GSM1417250, bins)
# ---- Gviz 시각화(원시 카운트 버전) --------------------------------------------
# 데이터 트랙 생성: 히스토그램 형태로 빈의 score 시각화
GSM1917772.track <- DataTrack(GSM1917772.200bins, type = "histogram",
                              strand = "*", genome = "hg38",
                              col.histogram = "gray", fill.histogram = "black",
                              name = "Camptothecin treated\n(GSM1917772)", col.axis = "black",
                              cex.axis = 0.4, ylim = c(0, 150))

GSM1917773.track <- DataTrack(GSM1917773.200bins, type = "histogram",
                              strand = "*", genome = "hg38",
                              col.histogram = "gray", fill.histogram = "black",
                              name = "TP53 Knockout\n(GSM1917773)", col.axis = "black",
                              cex.axis = 0.4, ylim = c(0, 150))

# 유전자 주석 트랙: Ensembl에서 해당 구간 유전자 정보 로드
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", mirror = "www")
bmTrack <- BiomartGeneRegionTrack(genome = "hg38",
                                  chromosome = chrnum,
                                  start = startposition, end = endposition,
                                  name = "RefSeq Genes",
                                  transcriptAnnotation = "symbol",
                                  biomart = mart)

# 좌표축 트랙
AT <- GenomeAxisTrack()


# Promoter 1
promoter1 <- GRanges(seqnames = "chr20",
                     ranges = IRanges(start = 50189189, end = 50190844),
                     name = "Promotor 1")

promoter1Track <- AnnotationTrack(promoter1,
                                  genome = "hg38",
                                  chromosome = "chr20",
                                  name = "Promotor1",
                                  fill = "red",
                                  col = "black",
                                  stacking = "dense",
                                  just.group = "above")

# Promoter 3
promoter3 <- GRanges(seqnames = "chr20",
                     ranges = IRanges(start = 50190208, end = 50190913),
                     name = "Promotor 3")

promoter3Track <- AnnotationTrack(promoter3,
                                  genome = "hg38",
                                  chromosome = "chr20",
                                  name = "Promotor3",
                                  fill = "blue",
                                  col = "black",
                                  stacking = "dense",
                                  just.group = "above")

displayPars(GSM1917772.track) <- list(cex.title = 0.88, col.title = "black",
                                      cex.axis = 0.8, ylim = c(0,100))
displayPars(GSM1917773.track) <- list(cex.title = 0.92, col.title = "black",
                                      cex.axis = 0.8, ylim = c(0,100))
displayPars(promoter1Track)   <- list(cex.title = 0.65, col.title = "black", rotation.title = 0)
displayPars(promoter3Track)   <- list(cex.title = 0.65, col.title = "black", rotation.title = 0)
displayPars(bmTrack)          <- list(cex.title = 1.2, col.title = "black", fontsize = 10)
displayPars(AT)               <- list(cex.title = 0.8, col.title = "black", cex.axis = 0.6)

# PNG 출력
CairoPNG("/home/Data_Drive_8TB/kykim/2. CRC/Chip-seq/HCT116(width)_raw.png", width = 1600, height = 1200, res = 150)

plotTracks(c(GSM1917772.track, GSM1917773.track ,  bmTrack, AT,
             promoter1Track, promoter3Track),
           from = startposition, to = endposition,
           transcriptAnnotation = "symbol",
           type = "histogram")
#cex.title = 0.2, fontsize = 1)
#
dev.off()

# ---- RPM 정규화 후 다시 시각화 -------------------------------------------------
# (2) RPM 스케일 빈
GSM1917772.200bins <- BinChIPseq_RPM(GSM1917772, bins)
GSM1917773.200bins <- BinChIPseq_RPM(GSM1917773, bins)

# RPM 스케일 데이터 트랙
GSM1917772.track <- DataTrack(GSM1917772.200bins, type = "histogram",
                              strand = "*", genome = "hg38",
                              col.histogram = "gray", fill.histogram = "black",
                              name = "GSM1917772 (RPM)", col.axis = "black",
                              cex.axis = 0.4, ylim = c(0, 150))

GSM1917773.track <- DataTrack(GSM1917773.200bins, type = "histogram",
                              strand = "*", genome = "hg38",
                              col.histogram = "gray", fill.histogram = "black",
                              name = "GSM1917773 (RPM)", col.axis = "black",
                              cex.axis = 0.4, ylim = c(0, 150))

# 화면 출력(RPM)
plotTracks(c(GSM1917772.track, GSM1917773.track, bmTrack, AT),
           from = startposition, to = endposition,
           transcriptAnnotation = "symbol",
           type = "histogram", cex.title = 0.7, fontsize = 10)

# 파일로 저장(RPM)
png("CEBPB_tracks_highres_signorm.png", width = 5, height = 5, units = "in", res = 300)
plotTracks(c(GSM1917772.track, GSM1917773.track, bmTrack, AT),
           from = startposition, to = endposition,
           transcriptAnnotation = "symbol",
           type = "histogram", cex.title = 0.7, fontsize = 10)
dev.off()