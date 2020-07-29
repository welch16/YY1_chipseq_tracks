
workdr=/ua/rwelch/Documents/CISR/YY1_chipseq_tracks/data/tracks/Jurkat/IKFZ1
codedr=/ua/rwelch/Documents/CISR/YY1_chipseq_tracks/tools

chromdr=/ua/rwelch/Documents/CISR/YY1_chipseq_tracks/data/tracks

bg2bw="$codedr/bedGraphToBigWig"
bedsort="$codedr/bedSort"

$bedsort "$workdr/GSM2759908_Jurkat_Helios_Rep1.bedGraph" \
  "$workdr/Jurkat_IKFZ1_rep1.bg"

$bedsort "$workdr/GSM2759909_Jurkat_Helios_Rep2.bedGraph" \
  "$workdr/Jurkat_IKFZ1_rep2.bg"

$bg2bw "$workdr/Jurkat_IKFZ1_rep1.bg" \
  "$chromdr/hg19.chrom.sizes" \
  "$workdr/Jurkat_IKFZ1_rep1.bigwig"

$bg2bw "$workdr/Jurkat_IKFZ1_rep2.bg" \
  "$chromdr/hg19.chrom.sizes" \
  "$workdr/Jurkat_IKFZ1_rep2.bigwig"
  