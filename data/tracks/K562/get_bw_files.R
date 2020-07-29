
pacman::p_load(magrittr, tidyverse, ENCODExplorer)
encode_df <- get_encode_df()


# explore the possible data
tf_data <- encode_df %>%
  dplyr::filter(biosample_name == "K562") %>%
  dplyr::filter(target %in% c("YY1", "SPI1", "IKZF1", "RUNX1")) %>%
  dplyr::filter(assembly == "hg19") %>%
  dplyr::filter(stringr::str_detect(output_type, "fold change")) %>%
  dplyr::filter(technical_replicates == "1_1; 2_1") %>%
  dplyr::filter(file_status == "released") %>%
  dplyr::filter(! stringr::str_detect(dataset_description, "v04")) %>%
  dplyr::filter(dataset_description != "ChIP-Seq on K562-human")

ENCODExplorer::downloadEncode(tf_data, format = "bigWig")

tf_data %>%
  tibble::as_tibble() %>%
  qs::qsave(here::here("data", "qs", "TFs_bigwig_K562.qs"))