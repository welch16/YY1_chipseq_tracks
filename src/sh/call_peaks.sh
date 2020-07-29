
workdr=/ua/rwelch/Documents/CISR/YY1_chipseq_tracks

chip_file=yy1_chip.bedGraph
input_file=yy1_input.bedGraph

macs="$workdr/tools/macs2"

$macs bdgpeakcall -i "$workdr/data/tracks/Jurkat/YY1/$chip_file" \
  -l 150 \
  -o "$workdr/data/peaks/Jurkat/YY1_peaks.txt"

  # -g hs \
  # --call-summits \
  # --nomodel --extsize 147 \
  # --mfold 1 50 \
  # -n YY1 \
  # -p .1 \