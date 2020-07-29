
workdr=/ua/rwelch/Documents/CISR/YY1_chipseq_tracks/

bedfile="$workdr/data/tracks/Jurkat/RUNX1/GSM1045363_RUNX1_rep1_alignment.bed"

bamfile="$workdr/data/tracks/Jurkat/RUNX1/Jurkat_RUNX1_rep1.bam"

codedr=/ua/rwelch/Documents/CISR/YY1_chipseq_tracks/tools

bedtools="$codedr/bedtools"
chr_sizes=hg19.chrom.sizes

$bedtools bedtobam -i $bedfile -g "$workdr/data/tracks/$chr_sizes" > \
  "$bamfile"
