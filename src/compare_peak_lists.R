
pacman::p_load(magrittr, tidyverse, ENCODExplorer,
  GenomicRanges, rtracklayer,
  ChIPseeker, TxDb.Hsapiens.UCSC.hg19.knownGene,
  EnsDb.Hsapiens.v75, clusterProfiler, AnnotationDbi,
  org.Hs.eg.db)

tf_data <-
  qs::qread(here::here("data", "qs", "multiple_TFs_bed_K562_peaks.qs"))

yy1_peaks <- tf_data %>%
  dplyr::filter(target == "YY1") %>%
  purrr::pluck("peaks", 1) %>%
  set_names(NULL) %>%
  unlist()

other_tfs <- tf_data %>%
  dplyr::filter(target != "YY1") %>%
  dplyr::select(target, peaks) %>%
  dplyr::mutate(
    peaks = purrr::map(peaks, set_names, NULL),
    peaks = purrr::map(peaks, unlist))

# compute distances between TFs and YY1
other_tfs %<>%
  dplyr::mutate(
    distances = purrr::map(peaks,
      ~ IRanges::distanceToNearest(yy1_peaks, .)))

yy1_dist <- other_tfs %>%
  dplyr::group_by(target) %>%
  dplyr::mutate(
    target = stringr::str_c(target, seq_len(n()), sep = "_"),
    distances = purrr::map(distances, tibble::as_tibble)) %>%
  dplyr::ungroup() %>%
  dplyr::select(target, distances) %>%
  unnest(cols = c(distances)) %>%
  dplyr::select(-subjectHits) %>%
  tidyr::pivot_wider(names_from = target, values_from = distance)

mcols(yy1_peaks) <- yy1_dist %>%
  dplyr::select(-queryHits) %>%
  DataFrame()

ww1 <- 5e3
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
annot_peaks <- ChIPseeker::annotatePeak(yy1_peaks, TxDb = txdb,
  tssRegion = c(-ww1, ww1), verbose = TRUE)

# Get the entrez IDs
entrez_annot <- annot_peaks@anno$geneId

# Return the gene symbol for the set of Entrez IDs
annotations_edb <- AnnotationDbi::select(
  EnsDb.Hsapiens.v75,
  keys = as.character(entrez_annot),
  columns = c("GENENAME"),
  keytype = "ENTREZID") %>%
  tibble::as_tibble() %>%
  dplyr::mutate(ENTREZID = as.character(ENTREZID)) %>%
  dplyr::distinct() %>%
  dplyr::group_by(ENTREZID) %>%
  dplyr::summarize(symbols = str_c(GENENAME, collapse = ","),
    .groups = "drop")

peaks_table <- annot_peaks@anno %>%
  as.data.frame() %>%
  tibble::as_tibble() %>%
  dplyr::rename(chrom = seqnames) %>%
  dplyr::select(-width, -strand) %>%
  dplyr::left_join(annotations_edb, by = c(geneId = "ENTREZID"))

symbols <- dplyr::pull(peaks_table, symbols)

mcols(yy1_peaks)[["genes"]] <- symbols

create_link <- function(chrom, start, end, ww = 0) {

  glue::glue(
    "https://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=welch16&hgS_otherUserSessionName=YY1_peaks&position={chr}%3A{start}%2D{end}",
    chr = chrom,
    start = start - ww,
    end = end + ww)
}

ww <- 10e3
candidates <- subset(yy1_peaks, between(SPI1_1, 10e3, 20e3) &
  (between(IKZF1_1, 10e3, 20e3) | between(IKZF1_2, 10e3, 20e3))) %>%
  as.data.frame() %>%
  tibble::as_tibble() %>%
  dplyr::rename(chrom = seqnames) %>%
  dplyr::select(-width, -strand) %>%
  dplyr::mutate(
    chrom = as.character(chrom),
    link = create_link(chrom, start, end, ww))
    
candidates %>%
  dplyr::filter(chrom == "chr19") %>% pul


candidates %>%
  dplyr::select(-genes) %>%
  readr::write_csv(here::here("data", "out",
    "K562_peaks_YY1_dist10KB_SPI1_IKZF1.csv"))

# NOTCH1 enhancer area
create_link("chr9", 139364317, 139459993)
 