
pacman::p_load(magrittr, tidyverse, ENCODExplorer,
  rtracklayer,
  ChIPseeker, TxDb.Hsapiens.UCSC.hg19.knownGene,
  EnsDb.Hsapiens.v75, clusterProfiler, AnnotationDbi,
  org.Hs.eg.db, ChIPpeakAnno)

data(TSS.human.GRCh37)

yy1_biosamples <-
  qs::qread(here::here("data", "qs", "YY1_peaks_per_biosample.qs"))

peaks <- purrr::pluck(yy1_biosamples, "peaks", 1) %>% unlist()
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

chr9_peaks <- subset(peaks, seqnames == "chr9")

# get notch1 symbol
human <- org.Hs.eg.db 
keys(human, keytype = "SYMBOL", pattern = "NOTCH")
c
notch1 <- AnnotationDbi::select(human,
  keys = "NOTCH1", columns = c("ENTREZID", "ENSEMBL"),
  keytype = "SYMBOL")

notch1_gene <- ChIPpeakAnno::toGRanges(txdb)[notch1$ENTREZID]
notch1_tss <- TSS.human.GRCh37[notch1$ENSEMBL]
notch1_tss <- with(notch1_tss,
  GRanges(
    seqnames = "chr9",
    ranges = IRanges(start = start, end = end),
    strand = strand))

resize_mid <- function(gr) {

  start(gr) <- floor(.5 * (start(gr) + end(gr)))
  width(gr) <- 1
  gr
  
}


summits <- peaks
summits <- resize_mid(summits)

# get bw files
data_files <- list.files(here::here("data", "tracks"), full.names = TRUE,
  recursive = TRUE)

# remove K562
data_files <- str_subset(data_files, "K562", negate = TRUE)
data_files <- c(
  data_files %>% str_subset("bigwig"),
  data_files %>% str_subset("bw"))
data_files <- data_files[str_detect(data_files, "YY1")]
data_files %<>% str_subset("RUNX1", negate = TRUE)
data_files %<>% str_subset("IKFZ", negate = TRUE)

bw <- map(data_files, rtracklayer::import.bw)


start(notch1_gene) <- end(notch1_gene)
start(notch1_tss) <- end(notch1_tss)

summit_window <- resize(summits, 2e3, fix = "center")
summit_window <- summit_window[nearest(notch1_tss, summits)]

# pick max
bw <- map(bw, subsetByOverlaps, summit_window)
bw <- map(bw, ~ .[which.max(.$score)])
bw <- map(bw, resize_mid)

map(bw, ~ distanceToNearest(notch1_tss, .))

distanceToNearest(notch1_gene, summits)
distanceToNearest(notch1_tss, summits)

> data_files
[1] "/z/Comp/uwcccseq/rwelch/Xuan_Pan/YY1_chipseq_tracks/data/tracks/CD4+/CD4_YY1_rep0.bigwig"                  
62836

[2] "/z/Comp/uwcccseq/rwelch/Xuan_Pan/YY1_chipseq_tracks/data/tracks/Jurkat/YY1/yy1_chip.bw"
62838

[4] "/z/Comp/uwcccseq/rwelch/Xuan_Pan/YY1_chipseq_tracks/data/tracks/PDX/GSM4321119_YY1_ChIPseq_treat_pileup.bw"
62868