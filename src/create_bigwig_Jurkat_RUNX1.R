
pacman::p_load(GenomicAlignments)

ga <- readGAlignments("data/tracks/Jurkat/RUNX1/Jurkat_RUNX1_rep1.bam")

gr <- as(ga, "GRanges")
gr_ext <- resize(gr, 200)

cover <- coverage(gr_ext)

rtracklayer::export.bw(
  cover, "data/tracks/Jurkat/RUNX1/Jurkat_RUNX1_rep1.bigwig")