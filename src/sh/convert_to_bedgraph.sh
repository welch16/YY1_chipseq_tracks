
workdr=/ua/rwelch/Documents/CISR/YY1_chipseq_tracks

chip_file=yy1_chip.bw
input_file=yy1_input.bw

bw2bg="$workdr/tools/bigWigToBedGraph"

$bw2bg "$workdr/data/tracks/Jurkat/YY1/$chip_file" \
  "$workdr/data/tracks/Jurkat/YY1/yy1_chip.bedGraph"

$bw2bg "$workdr/data/tracks/Jurkat/YY1/$input_file" \
  "$workdr/data/tracks/Jurkat/YY1/yy1_input.bedGraph"
