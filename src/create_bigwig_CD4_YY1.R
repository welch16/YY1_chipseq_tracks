
pacman::p_load(liftOver, vroom, GenomicAlignments,
  rtracklayer, GenomicRanges)

# read aligned reads from bed file
aligned_reads <- vroom::vroom(here::here("data",
  "tracks", "CD4+", "GSM630810_YY1.aligned.bed"),
  col_names = c("chr", "start", "end", "x", "y", "strand"))
gr <- with(aligned_reads,
  GRanges(seqnames = chr,
    ranges = IRanges(start = start, end = end),
    strand = strand))

# convert to hg19 from hg18
chain_file <- here::here("data", "tracks",
  "CD4+", "hg18ToHg19.over.chain")
chain <- import.chain(chain_file)
gr <- liftOver(gr, chain)

# extend reads (for a better figure, would have to
# estimate the fragment length with SCC)
gr_ext <- resize(gr, 200)
cover <- coverage(gr_ext)

rtracklayer::export.bw(
  cover,
  here::here("data", "tracks", "CD4+", "CD4_YY1_rep0.bigwig"))
