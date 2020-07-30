
workdr=/z/Comp/uwcccseq/rwelch/zz_other_collab/YY1_chipseq_tracks/data/figs

dirj="$workdr/Jurkat"
dirk="$workdr/K562"
w=10000
wk="10K"
dpi=1000


pal=Blues

computeMatrix scale-regions -S "$dirj/YY1_Jurkat_rep0.bigWig" \
                 "$dirj/RUNX1_Jurkat_rep0.bigWig" \
                 "$dirj/IKFZ1_Jurkat_rep1.bigWig" \
                 "$dirj/IKFZ1_Jurkat_rep2.bigWig" \
              -R "$dirj/Jurkat.bed" \
              --regionBodyLength "$w" \
              --skipZeros -o "$dirj/matrix.mat.gz"

plotHeatmap -m "$dirj/matrix.mat.gz" \
     -out "$dirj/Jurkat_heatmap.png" \
     --colorMap "$pal" \
     --whatToShow 'heatmap and colorbar' \
     --zMin 0 --zMax 10 \
     --startLabel "s. - $wk" \
     --endLabel "s. + $wk" \
     --xAxisLabel "distance to summit (bp)" \
     --dpi $dpi

computeMatrix scale-regions -S "$dirk/YY1_K562_rep0.bigWig" \
                               "$dirk/RUNX1_K562_rep0.bigWig" \
                               "$dirk/SPI1_K562_rep0.bigWig" \
                               "$dirk/IKFZ1_K562_rep1.bigWig" \
                               "$dirk/IKFZ1_K562_rep2.bigWig" \
              -R "$dirk/K562.bed" \
              --regionBodyLength "$w" \
              --skipZeros -o "$dirk/matrix.mat.gz"

plotHeatmap -m "$dirk/matrix.mat.gz" \
     -out "$dirk/K562_heatmap.png" \
     --colorMap "$pal" \
     --whatToShow 'heatmap and colorbar' \
     --startLabel "s. - $wk" \
     --endLabel "s. + $wk" \
     --zMin 0 --zMax 10 \
     --xAxisLabel "distance to summit (bp)" \
     --kmeans 4 \
     --dpi $dpi

w=500
wk="500"
pal=Blues

computeMatrix scale-regions -S "$dirk/YY1_K562_rep0.bigWig" \
                               "$dirk/RUNX1_K562_rep0.bigWig" \
                               "$dirk/SPI1_K562_rep0.bigWig" \
                               "$dirk/IKFZ1_K562_rep1.bigWig" \
                               "$dirk/IKFZ1_K562_rep2.bigWig" \
              -R "$dirk/K562.bed" \
              --regionBodyLength "$w" \
              --skipZeros -o "$dirk/matrix_mini.mat.gz"

plotHeatmap -m "$dirk/matrix_mini.mat.gz" \
     -out "$dirk/K562_heatmap_500.png" \
     --colorMap "$pal" \
     --whatToShow 'heatmap and colorbar' \
     --startLabel "s. - $wk" \
     --endLabel "s. + $wk" \
     --zMin 0 --zMax 10 \
     --xAxisLabel "distance to summit (bp)" \
     --kmeans 4 \
     --dpi $dpi


computeMatrix scale-regions -S "$dirj/YY1_Jurkat_rep0.bigWig" \
                 "$dirj/RUNX1_Jurkat_rep0.bigWig" \
                 "$dirj/IKFZ1_Jurkat_rep1.bigWig" \
                 "$dirj/IKFZ1_Jurkat_rep2.bigWig" \
              -R "$dirj/Jurkat.bed" \
              --regionBodyLength "$w" \
              --skipZeros -o "$dirj/matrix_mini.mat.gz"

plotHeatmap -m "$dirj/matrix_mini.mat.gz" \
     -out "$dirj/Jurkat_heatmap_500.png" \
     --colorMap "$pal" \
     --whatToShow 'heatmap and colorbar' \
     --zMin 0 --zMax 10 \
     --startLabel "s. - $wk" \
     --endLabel "s. + $wk" \
     --xAxisLabel "distance to summit (bp)" \
     --dpi $dpi
