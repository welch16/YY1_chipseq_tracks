
pacman::p_load(magrittr, tidyverse, ENCODExplorer)


yy1_biosamples <-
  qs::qread(here::here("data", "qs", "YY1_peaks_per_biosample.qs"))

peaks <- purrr::pluck(yy1_biosamples, "peaks", 1) %>% unlist()
ww <- 10e3


peaks_table <- tibble::tibble(
  chrom = as.character(GenomicRanges::seqnames(peaks)),
  start = GenomicRanges::start(peaks),
  end =   GenomicRanges::end(peaks)) %>%
  dplyr::mutate(
    link = glue::glue(
  "https://genome.ucsc.edu/cgi-bin/hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=welch16&hgS_otherUserSessionName=YY1_peaks&position={chr}%3A{start}%2D{end}",
    chr = chrom,
    start = start - ww,
    end = end + ww))


# /hgTracks?hgS_doOtherUser=submit&hgS_otherUserName=welch16&hgS_otherUserSessionName=YY1_peaks

# https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr1%3A120454176%2D120612317&hgsid=865638117_Ixb7iWdoDxIgi12KiWrjnklDeRJP

peaks_table %>%
  readr::write_csv(here::here("data", "out", "K562_peaks.csv"))