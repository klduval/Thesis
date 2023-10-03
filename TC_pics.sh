#!/bin/bash
#SBATCH --job-name=K9_timecourse
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=60gb
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kld57880@uga.edu

BASEDIR="/scratch/kld57880/TC_final"
WORK="/work/mglab/kld"

###lets make these bedgraphs into bigwigs for data visualization
module load ucsc/359
mkdir $BASEDIR/bws

for infile in $BASEDIR/bdgrphs/*norm.bga
do
  base=$(basename ${infile} .norm.bga)
  bedSort $infile $infile
  bedGraphToBigWig $infile $BASEDIR/genome/chrNameLength.txt $BASEDIR/bws/$base.bw
done

###lets do some broad comparisons to see what our data looks like before moving on
module load deepTools/3.3.1-intel-2019b-Python-3.7.4

multiBigwigSummary bins -b $BASEDIR/bws/*.bw -o $BASEDIR/bw_summ.npz -p 24
plotCorrelation -in $BASEDIR/bw_summ.npz -c spearman -p heatmap -o $BASEDIR/timecourse_bw_summ_heatmap.pdf
plotPCA -in $BASEDIR/bw_summ.npz -o $BASEDIR/timecourse_bw_summ_PCA.pdf

###I want to make merged/average bigwig files
bigwigCompare -b1 $BASEDIR/bws/2hpf_K9_1.bw -b2 $BASEDIR/bws/2hpf_K9_2.bw --operation mean -bs 10 -p 20 -o $BASEDIR/bws/2hpf_K9_AVG.bw
bigwigCompare -b1 $BASEDIR/bws/2.5hpf_K9_1.bw -b2 $BASEDIR/bws/2.5hpf_K9_2.bw --operation add --scaleFactors 0.333:0.333 -bs 10 -p 20 -o $BASEDIR/bws/2.5hpf_K9_rep1rep2.bw
bigwigCompare -b1 $BASEDIR/bws/2.5hpf_K9_rep1rep2.bw -b2 $BASEDIR/bws/2.5hpf_K9_3.bw --operation add --scaleFactors 1:0.333 -bs 10 -p 20 -o $BASEDIR/bws/2.5hpf_K9_AVG.bw
bigwigCompare -b1 $BASEDIR/bws/3hpf_K9_1.bw -b2 $BASEDIR/bws/3hpf_K9_2.bw --operation add --scaleFactors 0.333:0.333 -bs 10 -p 20 -o $BASEDIR/bws/3hpf_K9_rep1rep2.bw
bigwigCompare -b1 $BASEDIR/bws/3hpf_K9_rep1rep2.bw -b2 $BASEDIR/bws/3hpf_K9_3.bw --operation add --scaleFactors 1:0.333 -bs 10 -p 20 -o $BASEDIR/bws/3hpf_K9_AVG.bw
bigwigCompare -b1 $BASEDIR/bws/3.5hpf_K9_1.bw -b2 $BASEDIR/bws/3.5hpf_K9_2.bw --operation add --scaleFactors 0.333:0.333 -bs 10 -p 20 -o $BASEDIR/bws/3.5hpf_K9_rep1rep2.bw
bigwigCompare -b1 $BASEDIR/bws/3.5hpf_K9_rep1rep2.bw -b2 $BASEDIR/bws/3.5hpf_K9_3.bw --operation add --scaleFactors 1:0.333 -bs 10 -p 20 -o $BASEDIR/bws/3.5hpf_K9_AVG.bw
bigwigCompare -b1 $BASEDIR/bws/4hpf_K9_1.bw -b2 $BASEDIR/bws/4hpf_K9_2.bw --operation add --scaleFactors 0.333:0.333 -bs 10 -p 20 -o $BASEDIR/bws/4hpf_K9_rep1rep2.bw
bigwigCompare -b1 $BASEDIR/bws/4hpf_K9_rep1rep2.bw -b2 $BASEDIR/bws/4hpf_K9_3.bw --operation add --scaleFactors 1:0.333 -bs 10 -p 20 -o $BASEDIR/bws/4hpf_K9_AVG.bw
bigwigCompare -b1 $BASEDIR/bws/4.5hpf_K9_1.bw -b2 $BASEDIR/bws/4.5hpf_K9_2.bw --operation add --scaleFactors 0.333:0.333 -bs 10 -p 20 -o $BASEDIR/bws/4.5hpf_K9_rep1rep2.bw
bigwigCompare -b1 $BASEDIR/bws/4.5hpf_K9_rep1rep2.bw -b2 $BASEDIR/bws/4.5hpf_K9_3.bw --operation add --scaleFactors 1:0.333 -bs 10 -p 20 -o $BASEDIR/bws/4.5hpf_K9_AVG.bw

rm $BASEDIR/bws/*rep1rep2*

###making pictures now
computeMatrix reference-point -S $BASEDIR/bws/2hpf_K9_AVG.bw $BASEDIR/bws/2.5hpf_K9_AVG.bw $BASEDIR/bws/3hpf_K9_AVG.bw $BASEDIR/bws/3.5hpf_K9_AVG.bw $BASEDIR/bws/4hpf_K9_AVG.bw $BASEDIR/bws/4.5hpf_K9_AVG.bw -R $BASEDIR/peaks/4.5hpf_K9_final.bed --referencePoint center -p max -a 10000 -b 10000 -bs 10 --missingDataAsZero -o $BASEDIR/matrices/TC_4.5peaks.gz
computeMatrix reference-point -S $BASEDIR/bws/3.5hpf_K9_AVG.bw $BASEDIR/bws/4hpf_K9_AVG.bw $BASEDIR/bws/4.5hpf_K9_AVG.bw -R $BASEDIR/peaks/4.5hpf_K9_final.bed --referencePoint center -p max -a 10000 -b 10000 -bs 10 --missingDataAsZero -o $BASEDIR/matrices/lateTC_4.5peaks.gz
computeMatrix reference-point -S $BASEDIR/bws/2hpf_K9_AVG.bw $BASEDIR/bws/2.5hpf_K9_AVG.bw $BASEDIR/bws/3hpf_K9_AVG.bw -R $BASEDIR/peaks/4.5hpf_K9_final.bed --referencePoint center -p max -a 10000 -b 10000 -bs 10 --missingDataAsZero -o $BASEDIR/matrices/earlyTC_4.5peaks.gz

for infile in $BASEDIR/matrices/*.gz
do
  base=$(basename ${infile} .gz)
  plotHeatmap -m $infile --colorMap Blues --legendLocation none --regionsLabel "Peaks" -o $BASEDIR/figs/"$base"_heatmap.pdf
done

###going to make matrices so I can compute peak "height" for a simple line chart
computeMatrix scale-regions -S $BASEDIR/bws/2hpf_K9_AVG.bw -R $BASEDIR/peaks/2hpf_K9_final.bed --missingDataAsZero -p max --outFileNameMatrix $BASEDIR/2hpf_peakmatrix.tab -o $BASEDIR/matrices/2hpf.gz
computeMatrix scale-regions -S $BASEDIR/bws/2.5hpf_K9_AVG.bw -R $BASEDIR/peaks/2.5hpf_K9_final.bed --missingDataAsZero -p max --outFileNameMatrix $BASEDIR/2.5hpf_peakmatrix.tab -o $BASEDIR/matrices/2.5hpf.gz
computeMatrix scale-regions -S $BASEDIR/bws/3hpf_K9_AVG.bw -R $BASEDIR/peaks/3hpf_K9_final.bed --missingDataAsZero -p max --outFileNameMatrix $BASEDIR/3hpf_peakmatrix.tab -o $BASEDIR/matrices/3hpf.gz
computeMatrix scale-regions -S $BASEDIR/bws/3.5hpf_K9_AVG.bw -R $BASEDIR/peaks/3.5hpf_K9_final.bed --missingDataAsZero -p max --outFileNameMatrix $BASEDIR/3.5hpf_peakmatrix.tab -o $BASEDIR/matrices/3.5hpf.gz
computeMatrix scale-regions -S $BASEDIR/bws/4hpf_K9_AVG.bw -R $BASEDIR/peaks/4hpf_K9_final.bed --missingDataAsZero -p max --outFileNameMatrix $BASEDIR/4hpf_peakmatrix.tab -o $BASEDIR/matrices/4hpf.gz
computeMatrix scale-regions -S $BASEDIR/bws/4.5hpf_K9_AVG.bw -R $BASEDIR/peaks/4.5hpf_K9_final.bed --missingDataAsZero -p max --outFileNameMatrix $BASEDIR/4.5hpf_peakmatrix.tab  -o $BASEDIR/matrices/4.5hpf.gz
