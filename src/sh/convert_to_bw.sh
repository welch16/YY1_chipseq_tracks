
workdr=/ua/rwelch/Documents/CISR/YY1_chipseq_tracks

chr_sizes=hg19.chrom.sizes
chip_file=GSM2773998_20170825_6177.8590.rpm.WIG
input_file=GSM2773999_20170825_6173.8592.rpm.WIG

wig2bw="$workdr/tools/wigToBigWig"

$wig2bw "$workdr/data/tracks/Jurkat/YY1/$chip_file" \
  "$workdr/data/tracks/$chr_sizes" \
  "$workdr/data/tracks/Jurkat/YY1/yy1_chip.bw" 

$wig2bw "$workdr/data/tracks/Jurkat/YY1/$input_file" \
  "$workdr/data/tracks/$chr_sizes" \
  "$workdr/data/tracks/Jurkat/YY1/yy1_input.bw"
