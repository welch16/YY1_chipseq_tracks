
pacman::p_load(magrittr, tidyverse, GenomicRanges)


yy1_peaks <- here::here("data/peaks/Jurkat/YY1_peaks.txt") %>%
  readr::read_delim(delim = "\t", col_names = FALSE, skip = 1)

peaks <- GRanges(
  seqnames = yy1_peaks$X1,
  ranges = IRanges(start = yy1_peaks$X2, end = yy1_peaks$X3))


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
  readr::write_csv(here::here("data", "out", "Jurkat_peaks.csv"))