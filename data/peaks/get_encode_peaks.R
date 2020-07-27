
pacman::p_load(magrittr, tidyverse, ENCODExplorer)
encode_df <- get_encode_df()


# explore the possible data
peak_data <- encode_df %>%
  dplyr::filter(str_detect(output_type, "peaks"))

yy1_biosamples <- peak_data %>%
  dplyr::filter(target == "YY1") %>%
  dplyr::count(biosample_name, output_type) %>%
  dplyr::filter(stringr::str_detect(output_type, "IDR"))

safe_query_consensus <- function(name, tf, assembly = "GRCh38") {
  ENCODExplorer::queryConsensusPeaks(
    biosample_name = name,
    assembly = assembly,
    target = tf)
}

yy1_biosamples %<>%
  tibble::as_tibble() %>%
  dplyr::select(biosample_name) %>%
  dplyr::distinct()
  
yy1_biosamples %<>%
  dplyr::mutate(
    peaks = purrr::map(
      biosample_name, safely(safe_query_consensus), "YY1",
      assembly = "GRCh38"))

other_tfs <- peak_data %>%
  dplyr::filter(stringr::str_detect(output_type, "IDR")) %>%
  dplyr::filter(
    biosample_name %in% yy1_biosamples[["biosample_name"]]) %>%
  dplyr::filter(target != "YY1") %>%
  dplyr::distinct(biosample_name, target) %>%
  dplyr::mutate(
    peaks = purrr::map2(
      biosample_name, target, safely(safe_query_consensus),
        assembly = "GRCh38"))

yy1_biosamples %<>%
  dplyr::filter(purrr::map_lgl(peaks, ~ is.null(.$error))) %>%
  dplyr::mutate(
    peaks = purrr::map(peaks, "result"),
    peaks = purrr::map(peaks, ENCODExplorer::peaks))

yy1_biosamples %>%
  qs::qsave(here::here("data", "qs", "YY1_peaks_per_biosample.qs"))

other_tfs %<>%
  dplyr::filter(purrr::map_lgl(peaks, ~ is.null(.$error))) %>%
  dplyr::mutate(
    peaks = purrr::map(peaks, "result"),
    peaks = purrr::map(peaks, ENCODExplorer::peaks))

other_tfs %>%
  tibble::as_tibble() %>%
  qs::qsave(here::here("data", "qs", "TFs_peaks_per_biosample.qs"))