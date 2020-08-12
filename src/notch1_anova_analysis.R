
pacman::p_load(magrittr, tidyverse, ENCODExplorer,
  GenomicRanges, rtracklayer,
  ChIPseeker, TxDb.Hsapiens.UCSC.hg19.knownGene,
  EnsDb.Hsapiens.v75, clusterProfiler, AnnotationDbi,
  org.Hs.eg.db)

tf_data <- qs::qread(here::here("data", "qs", "TFs_bigwig_K562.qs"))

bw_dir <- here::here("data", "tracks", "K562")

tf_data %<>%
  dplyr::select(target, file_accession) %>%
  dplyr::mutate(
    tracks = purrr::map(
      glue::glue("{dir}/{file}.bigWig",
      dir = bw_dir, file = file_accession),
      rtracklayer::import.bw))

notch1_locus <- GenomicRanges::GRanges(
  seqnames = "chr9",
  ranges = IRanges(start = 139375500, end = 139450500))

# chr9:139,375,587-139,450,524
tf_data %<>%
  mutate(
    reads = map(tracks, subsetByOverlaps, notch1_locus),
    locus_coverage = map(reads,
      ~ coverage(., weight = .$score)[notch1_locus]),
    signal = map(locus_coverage, ~ as(.[["chr9"]], "DFrame")))

# find the local maxima coordinates

coords <- seq(start(notch1_locus), end(notch1_locus))

get_peak_signal1 <- function(reads, all_peaks) {

  cover <- IRanges::coverage(reads, weight = reads$score)
  out <- all_peaks %>%
    as.data.frame() %>%
    select(-width, -strand) %>%
    mutate(peak = str_c("peak", seq_len(n())))

  bind_cols(
    out,
    tibble::tibble(
      max_signal = map_dbl(seq_along(all_peaks),
        ~ max(cover[all_peaks[.]]))))
}


get_peak_signal <- function(tf_data, ww) {

  tf_data %<>%
    mutate(
      peaks = map(signal, tibble::as_tibble),
      peaks = map(peaks, pull, X),
      max_signal = map_dbl(peaks, max),
      peaks = map(peaks, splus2R::peaks, span = ww + 1),
      ranges = map(peaks,
        ~ GRanges(seqnames = "chr9",
          ranges = IRanges(coords[.]))))

  all_peaks <- do.call(c, tf_data$ranges)
  all_peaks <- reduce(
    resize(all_peaks, width = ww, fix = "center"))

  tf_data %<>%
    mutate(
      peak_signal = map(reads, get_peak_signal1, all_peaks))

  tf_data %>%
    select(target, peak_signal) %>%
    unnest(cols = c(peak_signal)) %>%
    group_by(target, peak, seqnames, start, end) %>%
    summarize(max_signal = max(max_signal), .groups = "drop")

}

peak_signal <- tibble::tibble(
  windows = c(100, 200, 500, 1000, 2000, 3500, 5000, 7500, 10000))
peak_signal %<>%
  mutate(
    data = map(windows, ~ get_peak_signal(tf_data, .)))

notch1_locus %>%
  qs::qsave(here::here("data", "qs", "notch1_locus.qs"))

peak_signal %>%
  qs::qsave(here::here("data", "qs", "notch1_peak_signals.qs"))
