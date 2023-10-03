#!/bin/bash
#SBATCH --job-name=MO_cut&run
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=80gb
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kld57880@uga.edu

BASEDIR="/scratch/kld57880/MO_CnR2"
TOOLDIR='/home/kld57880/Git2/toolbox'


###trimming, multiQC, and aligning to danio genome
for file in $BASEDIR/raw/*_R*.fastq.gz;
do
 if [[ $prefix ]]; then
       base=$(basename ${first} _L001_R1_001.fastq.gz)
       sh $TOOLDIR/PE_trim_and_star.sh -o $BASEDIR -n $base -m one $first $file
       prefix=
   else
       first=$file
       prefix=${file%%_*}
   fi
done

###Remove PCR duplicates####
ml picard/2.27.4-Java-13.0.2

for infile in $BASEDIR/bams/*q1.bam
do
 base=$(basename ${infile} _q1.bam)
 java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I $infile -M $BASEDIR/bams/"$base"_dupmetrics.txt -O $BASEDIR/bams/"$base"_nodups.bam --REMOVE_DUPLICATES true
done

####spike in normalization####
mkdir $BASEDIR/bdgrphs

for file in $BASEDIR/bams/*nodups.bam;
do
   base=$(basename ${file} _nodups.bam)
   sh $TOOLDIR/kmet_spike.sh $BASEDIR/bdgrphs $base $BASEDIR/trimmed/"$base"_L001_R1_001_val_1.fq.gz \
   $BASEDIR/trimmed/"$base"_L001_R2_001_val_2.fq.gz $file bga $BASEDIR/genome/chrNameLength.txt
done

mkdir $BASEDIR/figs
module load SAMtools/1.16.1-GCC-11.3.0

###lets make  bedgraphs into bigwigs for data visualization###
module load ucsc/443
mkdir $BASEDIR/bws

for infile in $BASEDIR/bdgrphs/*kmet.bga
do
 base=$(basename ${infile} _kmet.bga)
 bedSort $infile $infile
 bedGraphToBigWig $infile $BASEDIR/genome/chrNameLength.txt $BASEDIR/bws/$base.bw
done

module load deepTools/3.5.2-foss-2022a

multiBigwigSummary bins -b $BASEDIR/bws/*.bw -o $BASEDIR/bw_summ.npz -p 24
plotCorrelation -in $BASEDIR/bw_summ.npz -c spearman -p heatmap -o $BASEDIR/figs/MO_bw_summ_heatmap.pdf
plotPCA -in $BASEDIR/bw_summ.npz -o $BASEDIR/figs/MO_bw_summ_PCA.pdf

###peak calling###
module load Homer/4.11-foss-2022a
mkdir $BASEDIR/peaks

for infile in $BASEDIR/bdgrphs/6*kmet.bga
 do base=$(basename ${infile} _kmet.bga)
 cat $infile | awk '{print $1 "\t" $2 "\t" $3 "\t" "+" "\t" "+" "\t" "+"}' > $BASEDIR/peaks/$base.bgato.bed
done

for infile in $BASEDIR/peaks/6*bgato.bed
 do base=$(basename ${infile} .bgato.bed)
 makeTagDirectory $BASEDIR/peaks/$base.BtB.tagdir $infile -format bed
done

for infile in $BASEDIR/peaks/6*K9*.tagdir
 do base=$(basename ${infile} .BtB.tagdir)
 findPeaks $infile -style histone -minDist 1000 -gsize 1.5e9 -F 4 -i $BASEDIR/peaks/*IgG*.tagdir -o $BASEDIR/peaks/$base.txt
done

for infile in $BASEDIR/peaks/6*.txt
do
 base=$(basename ${infile} .txt)
 sed '/^#/d' $infile | awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $8 "\t" $5 "\t" $6 "\t" $12 "\t" "-1"}' | sed 's/\.000000//g' > $BASEDIR/peaks/$base.peaks.bed
done

###ChIPr time!###
module load ChIP-R
chipr -i $BASEDIR/peaks/1_suvAB_MO_K9_4_5h_1_S1.peaks.bed $BASEDIR/peaks/2_suvAB_MO_K9_4_5h_2_S2.peaks.bed $BASEDIR/peaks/3_suvAB_MO_K9_4_5h_3_S3.peaks.bed -m 2 -o $BASEDIR/peaks/suvAB_K9_4.5h_repPeaks
chipr -i $BASEDIR/peaks/4_setAB_MO_K9_4_5h_1_S4.peaks.bed $BASEDIR/peaks/5_setAB_MO_K9_4_5h_2_S5.peaks.bed $BASEDIR/peaks/6_setAB_MO_K9_4_5h_3_S6.peaks.bed -m 2 -o $BASEDIR/peaks/setAB_K9_4.5h_repPeaks
chipr -i $BASEDIR/peaks/7_gfp_MO_K9_4_5h_1_S7.peaks.bed $BASEDIR/peaks/8_gfp_MO_K9_4_5h_2_S8.peaks.bed $BASEDIR/peaks/9_gfp_MO_K9_4_5h_3_S9.peaks.bed -m 2 -o $BASEDIR/peaks/gfpMO_K9_4.5h_repPeaks

###differential peaks###
module load BEDTools

cat $BASEDIR/peaks/*_all.bed | bedtools sort -i stdin | bedtools merge -i stdin > $BASEDIR/peaks/MO_MERGEDpeaks.bed
bedtools intersect -a $BASEDIR/peaks/MO_MERGEDpeaks.bed -b $BASEDIR/peaks/gfpMO_K9_4.5h_repPeaks_all.bed $BASEDIR/peaks/setAB_K9_4.5h_repPeaks_all.bed $BASEDIR/peaks/suvAB_K9_4.5h_repPeaks_all.bed \
-wa -wb -names gfp setAB suvAB > $BASEDIR/peaks/MO_interPeaks.bed
sh $TOOLDIR/multi_inter.sh $BASEDIR/peaks/MO_interPeaks.bed $BASEDIR/peaks

####peak annotation####
curl -s ftp://ftp.ensembl.org/pub/release-98/gtf/danio_rerio/Danio_rerio.GRCz11.98.gtf.gz | gunzip -c > $BASEDIR/refann.gtf
mkdir $BASEDIR/peaks/ann

for infile in $BASEDIR/peaks/*all.bed
do
 base=$( basename ${infile} _repPeaks_all.bed)
 annotatePeaks.pl $infile danRer11 -gtf $BASEDIR/refann.gtf > $BASEDIR/peaks/ann/$base.maskann.txt
done

for infile in $BASEDIR/peaks/ann/*maskann.txt
do
 base=$(basename ${infile} .maskann.txt)
 awk -F'\t' 'sqrt($10*$10) <=1000' $infile > $BASEDIR/peaks/ann/$base.1000bp_ann.txt
done

for infile in $BASEDIR/peaks/ann/*maskann.txt
do
 base=$(basename ${infile} .maskann.txt)
 awk -F'\t' 'sqrt($10*$10) >=1000' $infile | awk '{print $2 "\t" $3 "\t" $4 }' > $BASEDIR/peaks/ann/$base.MOREthan1000bp.bed
done

for infile in $BASEDIR/peaks/ann/*.MOREthan1000bp.bed
do
 base=$( basename ${infile} .MOREthan1000bp.bed)
 bedtools intersect -a $infile -b /work/mglab/kld/TEfiles_Chang_etal/TEann.gtf -F 0.50 -wo > $BASEDIR/peaks/ann/$base.TEann.txt
done

for infile in $BASEDIR/peaks/ann/*TEann.txt
do
	base=$(basename ${infile} .txt)
 sed 's:;:\t:g' $infile > $BASEDIR/peaks/ann/$base.trim.txt
done

for infile in $BASEDIR/peaks/ann/*TEann.trim.txt
do
	base=$(basename ${infile} .trim.txt)
	sed 's:"::g' $infile | sed 's:gene_id::g' | sed 's:transcript_id::g' | sed 's:family_id::g' | sed 's:class_id::g' | sed 's:Cluster::g' | sed 's:TE_classification::g' > $BASEDIR/peaks/ann/$base.txt
done

###averaged bws###
bigwigCompare -b1 $BASEDIR/bws/1_suvAB_MO_K9_4_5h_1_S1.bw -b2 $BASEDIR/bws/2_suvAB_MO_K9_4_5h_2_S2.bw --operation add --scaleFactors 0.333:0.333 -bs 10 -p 20 -o $BASEDIR/bws/suvAB_4.5h_K9_rep1rep2.bw
bigwigCompare -b1 $BASEDIR/bws/suvAB_4.5h_K9_rep1rep2.bw -b2 $BASEDIR/bws/3_suvAB_MO_K9_4_5h_3_S3.bw --operation add --scaleFactors 1:0.333 -bs 10 -p 20 -o $BASEDIR/bws/suvAB_4.5h_K9_AVG.bw

bigwigCompare -b1 $BASEDIR/bws/4_setAB_MO_K9_4_5h_1_S4.bw -b2 $BASEDIR/bws/5_setAB_MO_K9_4_5h_2_S5.bw --operation add --scaleFactors 0.333:0.333 -bs 10 -p 20 -o $BASEDIR/bws/setAB_4.5h_K9_rep1rep2.bw
bigwigCompare -b1 $BASEDIR/bws/setAB_4.5h_K9_rep1rep2.bw -b2 $BASEDIR/bws/6_setAB_MO_K9_4_5h_3_S6.bw --operation add --scaleFactors 1:0.333 -bs 10 -p 20 -o $BASEDIR/bws/setAB_4.5h_K9_AVG.bw

bigwigCompare -b1 $BASEDIR/bws/7_gfp_MO_K9_4_5h_1_S7.bw -b2 $BASEDIR/bws/8_gfp_MO_K9_4_5h_2_S8.bw --operation add --scaleFactors 0.333:0.333 -bs 10 -p 20 -o $BASEDIR/bws/gfp_4.5h_K9_rep1rep2.bw
bigwigCompare -b1 $BASEDIR/bws/gfp_4.5h_K9_rep1rep2.bw -b2 $BASEDIR/bws/9_gfp_MO_K9_4_5h_3_S9.bw --operation add --scaleFactors 1:0.333 -bs 10 -p 20 -o $BASEDIR/bws/gfp_4.5h_K9_AVG.bw

###make a couple pictures###
mkdir $BASEDIR/matrices

computeMatrix reference-point -S $BASEDIR/bws/*AVG*.bw -R $BASEDIR/peaks/gfpMO_K9_4.5h_repPeaks_all.bed --referencePoint center -p max -a 10000 -b 10000 -bs 10 --missingDataAsZero -o $BASEDIR/matrices/K9AVG_gpfPeaks.gz
computeMatrix reference-point -S $BASEDIR/bws/*AVG*.bw -R $BASEDIR/peaks/suvAB_K9_4.5h_repPeaks_all.bed --referencePoint center -p max -a 10000 -b 10000 -bs 10 --missingDataAsZero -o $BASEDIR/matrices/K9AVG_suvPeaks.gz
computeMatrix reference-point -S $BASEDIR/bws/*AVG*.bw -R $BASEDIR/peaks/setMO_K9_4.5h_repPeaks_all.bed --referencePoint center -p max -a 10000 -b 10000 -bs 10 --missingDataAsZero -o $BASEDIR/matrices/K9AVG_setPeaks.gz

for infile in $BASEDIR/matrices/*.gz
do
  base=$(basename ${infile} .gz)
  plotHeatmap -m $infile --colorMap Blues -o $BASEDIR/figs/$base.heatmap.pdf
  plotHeatmap -m $infile --colorMap Blues --kmeans 2 -o $BASEDIR/figs/$base.kmeans2.heatmap.pdf
  plotHeatmap -m $infile --colorMap Blues --kmeans 3 -o $BASEDIR/figs/$base.kmeans3.heatmap.pdf
done

###calculate SAT1 read depth in each treatment
ml deepTools

computeMatrix scale-regions -S $BASEDIR/bws/gfp_4.5h_K9_AVG.bw -R /scratch/kld57880/TC_final/peric/sat1_trimmed.bed --missingDataAsZero -p max --outFileNameMatrix $BASEDIR/matrices/gfpMO_sat1matrix.tab -o $BASEDIR/matrices/gfpMO_sat1.gz
computeMatrix scale-regions -S $BASEDIR/bws/setAB_4.5h_K9_AVG.bw -R /scratch/kld57880/TC_final/peric/sat1_trimmed.bed --missingDataAsZero -p max --outFileNameMatrix $BASEDIR/matrices/setMO_sat1matrix.tab -o $BASEDIR/matrices/setMO_sat1.gz
computeMatrix scale-regions -S $BASEDIR/bws/suvAB_4.5h_K9_AVG.bw -R /scratch/kld57880/TC_final/peric/sat1_trimmed.bed --missingDataAsZero -p max --outFileNameMatrix $BASEDIR/matrices/suvMO_sat1matrix.tab -o $BASEDIR/matrices/suvMO_sat1.gz

####peak annotation####
curl -s ftp://ftp.ensembl.org/pub/release-98/gtf/danio_rerio/Danio_rerio.GRCz11.98.gtf.gz | gunzip -c > $BASEDIR/refann.gtf
mkdir $BASEDIR/peaks/ann

for infile in $BASEDIR/peaks/diff/*.bed
do
 base=$( basename ${infile} .bed)
 annotatePeaks.pl $infile danRer11 -gtf $BASEDIR/refann.gtf > $BASEDIR/peaks/diff/$base.maskann.txt
done

for infile in $BASEDIR/peaks/diff/*maskann.txt
do
 base=$(basename ${infile} .maskann.txt)
 awk -F'\t' 'sqrt($10*$10) <=1000' $infile > $BASEDIR/peaks/diff/$base.1000bp_ann.txt
done

for infile in $BASEDIR/peaks/diff/*maskann.txt
do
 base=$(basename ${infile} .maskann.txt)
 awk -F'\t' 'sqrt($10*$10) >=1000' $infile | awk '{print $2 "\t" $3 "\t" $4 }' > $BASEDIR/peaks/diff/$base.MOREthan1000bp.bed
done

for infile in $BASEDIR/peaks/diff/*.MOREthan1000bp.bed
do
 base=$( basename ${infile} .MOREthan1000bp.bed)
 bedtools intersect -a $infile -b /work/mglab/kld/TEfiles_Chang_etal/TEann.gtf -F 0.50 -wo > $BASEDIR/peaks/diff/$base.TEann.txt
done

for infile in $BASEDIR/peaks/diff/*TEann.txt
do
	base=$(basename ${infile} .txt)
 sed 's:;:\t:g' $infile > $BASEDIR/peaks/diff/$base.trim.txt
done

for infile in $BASEDIR/peaks/diff/*TEann.trim.txt
do
	base=$(basename ${infile} .trim.txt)
	sed 's:"::g' $infile | sed 's:gene_id::g' | sed 's:transcript_id::g' | sed 's:family_id::g' | sed 's:class_id::g' | sed 's:Cluster::g' | sed 's:TE_classification::g' | awk '{print $12 "\t" $16 "\t" $17}' > $BASEDIR/peaks/diff/$base.txt
done
