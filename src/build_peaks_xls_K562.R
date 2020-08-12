
pacman::p_load(magrittr, tidyverse, ENCODExplorer,
  ChIPseeker, TxDb.Hsapiens.UCSC.hg19.knownGene,
  EnsDb.Hsapiens.v75, clusterProfiler, AnnotationDbi,
  org.Hs.eg.db)


yy1_biosamples <-
  qs::qread(here::here("data", "qs", "YY1_peaks_per_biosample.qs"))

peaks <- purrr::pluck(yy1_biosamples, "peaks", 1) %>% unlist()
ww <- 10e3
ww1 <- 1e3

names(peaks) <- NULL

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

annot_peaks <- ChIPseeker::annotatePeak(peaks, TxDb = txdb,
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
    dplyr::mutate(ENTREZID = as.character(ENTREZID))

create_link <- function(chrom, start, end, ww = 0) {

  glue::glue(
    "https://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=welch16&hgS_otherUserSessionName=YY1_peaks&position={chr}%3A{start}%2D{end}",
    chr = chrom,
    start = start - ww,
    end = end + ww)
}

peaks_table <- annot_peaks@anno %>%
  as.data.frame() %>%
  tibble::as_tibble() %>%
  dplyr::rename(chrom = seqnames) %>%
  dplyr::select(-width, -strand) %>%
  dplyr::left_join(annotations_edb, by = c(geneId = "ENTREZID")) %>%
  dplyr::rename(symbol = GENENAME) %>%
  dplyr::select(-starts_with("gene")) %>%
  dplyr::mutate(
    chrom = as.character(chrom),
    link = create_link(chrom, start, end, ww))



peaks_table %>%
  readr::write_csv(here::here("data", "out", "K562_peaks.csv"))

# NOTCH1 enhancer area
create_link("chr9", 139364317, 139459993)
 