
pacman::p_load(magrittr, tidyverse)

jurkat <- here::here("data", "out", "Jurkat_peaks.csv") %>%
  readr::read_csv() %>%
  dplyr::select(-link)

k562 <- here::here("data", "out", "K562_peaks.csv") %>%
  readr::read_csv() %>%
  dplyr::select(-link)

jurkat %<>% dplyr::filter(chrom != "chrM")
k562 %<>% dplyr::filter(chrom != "chrM")


jurkat %<>%
  dplyr::mutate_if(is.numeric, list( ~ . - 1)) %>%
  readr::write_tsv(here::here("data", "figs", "Jurkat", "Jurkat.bed"))

k562 %<>%
  dplyr::mutate_if(is.numeric, list( ~ . - 1)) %>%
  readr::write_tsv(here::here("data", "figs", "K562", "K562.bed"))  