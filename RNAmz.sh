#!/bin/bash
#SBATCH --job-name=MZ_RNA
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=160gb
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kld57880@uga.edu

BASEDIR="/scratch/kld57880/RNAseq_suvMZ_6.2023"
TOOLDIR='/home/kld57880/Git2/toolbox'

###trimming, multiQC, and aligning to  genome
for file in $BASEDIR/raw/*L001*001.fastq.gz;
do
 if [[ $prefix ]]; then
       base=$(basename ${first} _L001_R1_001.fastq.gz)
       sh $TOOLDIR/RNApe_trim_n_star.sh -o $BASEDIR -n $base -m one $first $file
       prefix=
   else
       first=$file
       prefix=${file%%_*}
   fi
done

####Remove PCR duplicates
ml picard/2.27.4-Java-13.0.2

for infile in $BASEDIR/bams/*q1.bam
do
 base=$(basename ${infile} _q1.bam)
 java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I $infile -M $BASEDIR/bams/"$base"_dupmetrics.txt -O $BASEDIR/bams/"$base"_nodups.bam --REMOVE_DUPLICATES true
done

###use featureCounts to generate count files for genes and repeats
module load Subread
curl -s ftp://ftp.ensembl.org/pub/release-98/gtf/danio_rerio/Danio_rerio.GRCz11.98.gtf.gz | gunzip -c > $BASEDIR/refann.gtf

mkdir $BASEDIR/counts
featureCounts -M -B -p -T 24 -g gene_name --minOverlap 20 -a $BASEDIR/refann.gtf -o $BASEDIR/counts/MZrna_refAnn_counts1.txt $BASEDIR/bams/*_nodups.bam
featureCounts -M -B -p -T 24 -g gene_id --minOverlap 20 -a $BASEDIR/refann.gtf -o $BASEDIR/counts/MZrna_refAnn_counts.txt $BASEDIR/bams/*_nodups.bam
featureCounts -M -B -p -T 24 -g gene_id --minOverlap 20 -a /work/mglab/kld/TEfiles_Chang_etal/TEann.gtf -o $BASEDIR/counts/MZrna_repeats_counts.txt $BASEDIR/bams/*_nodups.bam

###make bigwigs
module load SAMtools/1.10-iccifort-2019.5.281
module load deepTools
mkdir $BASEDIR/bws

for infile in $BASEDIR/bams/*_nodups.bam
do
 base=$(basename ${infile} _nodups.bam)
 samtools index -@ 24 $infile $BASEDIR/bams/"$base".bai
 bamCoverage -b $infile -bs 25 -p 24 --normalizeUsing BPM -o $BASEDIR/bws/"$base".TPMnorm.bw
done

###averaged bw files
bigwigCompare -b1 $BASEDIR/bws/1_AB_2h_TC_1_S1.TPM.bw -b2 $BASEDIR/bws/2_AB_2h_TC_2_S2.TPM.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/AB_2h_TC_rep1.bw
bigwigCompare -b1 $BASEDIR/bws/3_AB_2h_TC_3_S3.TPM.bw -b2 $BASEDIR/bws/4_AB_2h_TC_4_S4.TPM.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/AB_2h_TC_rep2.bw
bigwigCompare -b1 $BASEDIR/bws/AB_2h_TC_rep1.bw -b2 $BASEDIR/bws/AB_2h_TC_rep2.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/AB_2h_TC_AVG.bw

bigwigCompare -b1 $BASEDIR/bws/5_AB_2_5h_TC_1_S5.TPM.bw -b2 $BASEDIR/bws/6_AB_2_5h_TC_2_S6.TPM.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/AB_2.5h_TC_rep1.bw
bigwigCompare -b1 $BASEDIR/bws/7_AB_2_5h_TC_3_S7.TPM.bw -b2 $BASEDIR/bws/8_AB_2_5h_TC_4_S8.TPM.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/AB_2.5h_TC_rep2.bw
bigwigCompare -b1 $BASEDIR/bws/AB_2.5h_TC_rep1.bw -b2 $BASEDIR/bws/AB_2.5h_TC_rep2.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/AB_2.5h_TC_AVG.bw

bigwigCompare -b1 $BASEDIR/bws/9_AB_3_5h_TC_1_S9.TPM.bw -b2 $BASEDIR/bws/10_AB_3_5h_TC_2_S10.TPM.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/AB_3.5h_TC_rep1.bw
bigwigCompare -b1 $BASEDIR/bws/11_AB_3_5h_TC_3_S11.TPM.bw -b2 $BASEDIR/bws/12_AB_3_5h_TC_4_S12.TPM.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/AB_3.5h_TC_rep2.bw
bigwigCompare -b1 $BASEDIR/bws/AB_3.5h_TC_rep1.bw -b2 $BASEDIR/bws/AB_3.5h_TC_rep2.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/AB_3.5h_TC_AVG.bw

bigwigCompare -b1 $BASEDIR/bws/13_AB_4_5h_TC_1_S13.TPM.bw -b2 $BASEDIR/bws/14_AB_4_5h_TC_2_S14.TPM.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/AB_4.5h_TC_rep1.bw
bigwigCompare -b1 $BASEDIR/bws/15_AB_4_5h_TC_3_S15.TPM.bw -b2 $BASEDIR/bws/16_AB_4_5h_TC_4_S16.TPM.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/AB_4.5h_TC_rep2.bw
bigwigCompare -b1 $BASEDIR/bws/AB_4.5h_TC_rep1.bw -b2 $BASEDIR/bws/AB_4.5h_TC_rep2.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/AB_4.5h_TC_AVG.bw

bigwigCompare -b1 $BASEDIR/bws/17_AB_6h_TC_1_S17.TPM.bw -b2 $BASEDIR/bws/18_AB_6h_TC_2_S18.TPM.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/AB_6h_TC_rep1.bw
bigwigCompare -b1 $BASEDIR/bws/19_AB_6h_TC_3_S19.TPM.bw -b2 $BASEDIR/bws/20_AB_6h_TC_4_S20.TPM.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/AB_6h_TC_rep2.bw
bigwigCompare -b1 $BASEDIR/bws/AB_6h_TC_rep1.bw -b2 $BASEDIR/bws/AB_6h_TC_rep2.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/AB_6h_TC_AVG.bw

bigwigCompare -b1 $BASEDIR/bws/24_aMZ_2h_TC_1_S24.TPM.bw -b2 $BASEDIR/bws/25_aMZ_2h_TC_2_S25.TPM.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/aMZ_2h_TC_rep1.bw
bigwigCompare -b1 $BASEDIR/bws/26_aMZ_2h_TC_3_S26.TPM.bw -b2 $BASEDIR/bws/27_aMZ_2h_TC_4_S27.TPM.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/aMZ_2h_TC_rep2.bw
bigwigCompare -b1 $BASEDIR/bws/aMZ_2h_TC_rep1.bw -b2 $BASEDIR/bws/aMZ_2h_TC_rep2.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/aMZ_2h_TC_AVG.bw

bigwigCompare -b1 $BASEDIR/bws/28_aMZ_2_5h_TC_1_S28.TPM.bw -b2 $BASEDIR/bws/29_aMZ_2_5h_TC_2_S29.TPM.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/aMZ_2.5h_TC_rep1.bw
bigwigCompare -b1 $BASEDIR/bws/30_aMZ_2_5h_TC_3_S30.TPM.bw -b2 $BASEDIR/bws/31_aMZ_2_5h_TC_4_S31.TPM.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/aMZ_2.5h_TC_rep2.bw
bigwigCompare -b1 $BASEDIR/bws/aMZ_2.5h_TC_rep1.bw -b2 $BASEDIR/bws/aMZ_2.5h_TC_rep2.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/aMZ_2.5h_TC_AVG.bw

bigwigCompare -b1 $BASEDIR/bws/32_aMZ_3_5h_TC_1_S32.TPM.bw -b2 $BASEDIR/bws/33_aMZ_3_5h_TC_2_S33.TPM.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/aMZ_3.5h_TC_rep1.bw
bigwigCompare -b1 $BASEDIR/bws/34_aMZ_3_5h_TC_3_S34.TPM.bw -b2 $BASEDIR/bws/35_aMZ_3_5h_TC_4_S35.TPM.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/aMZ_3.5h_TC_rep2.bw
bigwigCompare -b1 $BASEDIR/bws/aMZ_3.5h_TC_rep1.bw -b2 $BASEDIR/bws/aMZ_3.5h_TC_rep2.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/aMZ_3.5h_TC_AVG.bw

bigwigCompare -b1 $BASEDIR/bws/*aMZ_4_5h_TC_1*.TPM.bw -b2 $BASEDIR/bws/*aMZ_4_5h_TC_2*.TPM.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/aMZ_4.5h_TC_rep1.bw
bigwigCompare -b1 $BASEDIR/bws/*aMZ_4_5h_TC_3*.TPM.bw -b2 $BASEDIR/bws/*aMZ_4_5h_TC_4*.TPM.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/aMZ_4.5h_TC_rep2.bw
bigwigCompare -b1 $BASEDIR/bws/aMZ_4.5h_TC_rep1.bw -b2 $BASEDIR/bws/aMZ_4.5h_TC_rep2.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/aMZ_4.5h_TC_AVG.bw

bigwigCompare -b1 $BASEDIR/bws/*aMZ_6h_TC_1*.TPM.bw -b2 $BASEDIR/bws/*aMZ_6h_TC_2*.TPM.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/aMZ_6h_TC_rep1.bw
bigwigCompare -b1 $BASEDIR/bws/*aMZ_6h_TC_3*.TPM.bw -b2 $BASEDIR/bws/*aMZ_6h_TC_4*.TPM.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/aMZ_6h_TC_rep2.bw
bigwigCompare -b1 $BASEDIR/bws/aMZ_6h_TC_rep1.bw -b2 $BASEDIR/bws/aMZ_6h_TC_rep2.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/aMZ_6h_TC_AVG.bw

bigwigCompare -b1 $BASEDIR/bws/*bMZ_6h_TC_1*.TPM.bw -b2 $BASEDIR/bws/*bMZ_6h_TC_2*.TPM.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/bMZ_6h_TC_rep1.bw
bigwigCompare -b1 $BASEDIR/bws/*bMZ_6h_TC_3*.TPM.bw -b2 $BASEDIR/bws/*bMZ_6h_TC_4*.TPM.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/bMZ_6h_TC_rep2.bw
bigwigCompare -b1 $BASEDIR/bws/bMZ_6h_TC_rep1.bw -b2 $BASEDIR/bws/bMZ_6h_TC_rep2.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/bMZ_6h_TC_AVG.bw

bigwigCompare -b1 $BASEDIR/bws/*bMZ_2h_TC_1*.TPM.bw -b2 $BASEDIR/bws/*bMZ_2h_TC_2*.TPM.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/bMZ_2h_TC_rep1.bw
bigwigCompare -b1 $BASEDIR/bws/*bMZ_2h_TC_3*.TPM.bw -b2 $BASEDIR/bws/*bMZ_2h_TC_4*.TPM.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/bMZ_2h_TC_rep2.bw
bigwigCompare -b1 $BASEDIR/bws/bMZ_2h_TC_rep1.bw -b2 $BASEDIR/bws/bMZ_2h_TC_rep2.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/bMZ_2h_TC_AVG.bw

bigwigCompare -b1 $BASEDIR/bws/*bMZ_2_5h_TC_1*.TPM.bw -b2 $BASEDIR/bws/*bMZ_2_5h_TC_2*.TPM.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/bMZ_2.5h_TC_rep1.bw
bigwigCompare -b1 $BASEDIR/bws/*bMZ_2_5h_TC_3*.TPM.bw -b2 $BASEDIR/bws/*bMZ_2_5h_TC_4*.TPM.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/bMZ_2.5h_TC_rep2.bw
bigwigCompare -b1 $BASEDIR/bws/bMZ_2.5h_TC_rep1.bw -b2 $BASEDIR/bws/bMZ_2.5h_TC_rep2.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/bMZ_2.5h_TC_AVG.bw

bigwigCompare -b1 $BASEDIR/bws/*bMZ_3_5h_TC_1*.TPM.bw -b2 $BASEDIR/bws/*bMZ_3_5h_TC_2*.TPM.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/bMZ_3.5h_TC_rep1.bw
bigwigCompare -b1 $BASEDIR/bws/*bMZ_3_5h_TC_3*.TPM.bw -b2 $BASEDIR/bws/*bMZ_3_5h_TC_4*.TPM.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/bMZ_3.5h_TC_rep2.bw
bigwigCompare -b1 $BASEDIR/bws/bMZ_3.5h_TC_rep1.bw -b2 $BASEDIR/bws/bMZ_3.5h_TC_rep2.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/bMZ_3.5h_TC_AVG.bw

bigwigCompare -b1 $BASEDIR/bws/*bMZ_4_5h_TC_1*.TPM.bw -b2 $BASEDIR/bws/*bMZ_4_5h_TC_2*.TPM.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/bMZ_4.5h_TC_rep1.bw
bigwigCompare -b1 $BASEDIR/bws/*bMZ_4_5h_TC_3*.TPM.bw -b2 $BASEDIR/bws/*bMZ_4_5h_TC_4*.TPM.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/bMZ_4.5h_TC_rep2.bw
bigwigCompare -b1 $BASEDIR/bws/bMZ_4.5h_TC_rep1.bw -b2 $BASEDIR/bws/bMZ_4.5h_TC_rep2.bw --operation mean -bs 10 -p 24 -o $BASEDIR/bws/bMZ_4.5h_TC_AVG.bw

####make some pictures just about overall transcription

computeMatrix reference-point -S $BASEDIR/bws/*2h*AVG*.bw -R $BASEDIR/genes.bed --missingDataAsZero --referencePoint TSS -p max -a 1000 -b 1000 -bs 10 -o $BASEDIR/matrices/2h_AVG_geneTSS.gz
computeMatrix scale-regions -S $BASEDIR/bws/*2h*AVG*.bw -R $BASEDIR/genes.bed --missingDataAsZero -p max -a 1000 -b 1000 -bs 10 -o $BASEDIR/matrices/2h_AVG_geneBody.gz

computeMatrix reference-point -S $BASEDIR/bws/*2.5h*AVG*.bw -R $BASEDIR/genes.bed --missingDataAsZero --referencePoint TSS -p max -a 1000 -b 1000 -bs 10 -o $BASEDIR/matrices/2.5h_AVG_geneTSS.gz
computeMatrix scale-regions -S $BASEDIR/bws/*2.5h*AVG*.bw -R $BASEDIR/genes.bed --missingDataAsZero -p max -a 1000 -b 1000 -bs 10 -o $BASEDIR/matrices/2.5h_AVG_geneBody.gz

computeMatrix reference-point -S $BASEDIR/bws/*3.5h*AVG*.bw -R $BASEDIR/genes.bed --missingDataAsZero --referencePoint TSS -p max -a 1000 -b 1000 -bs 10 -o $BASEDIR/matrices/3.5h_AVG_geneTSS.gz
computeMatrix scale-regions -S $BASEDIR/bws/*3.5h*AVG*.bw -R $BASEDIR/genes.bed --missingDataAsZero -p max -a 1000 -b 1000 -bs 10 -o $BASEDIR/matrices/3.5h_AVG_geneBody.gz

computeMatrix reference-point -S $BASEDIR/bws/*4.5h*AVG*.bw -R $BASEDIR/genes.bed --missingDataAsZero --referencePoint TSS -p max -a 1000 -b 1000 -bs 10 -o $BASEDIR/matrices/4.5h_AVG_geneTSS.gz
computeMatrix scale-regions -S $BASEDIR/bws/*4.5h*AVG*.bw -R $BASEDIR/genes.bed --missingDataAsZero -p max -a 1000 -b 1000 -bs 10 -o $BASEDIR/matrices/4.5h_AVG_geneBody.gz

computeMatrix reference-point -S $BASEDIR/bws/*6h*AVG*.bw -R $BASEDIR/genes.bed --missingDataAsZero --referencePoint TSS -p max -a 1000 -b 1000 -bs 10 -o $BASEDIR/matrices/6h_AVG_geneTSS.gz
computeMatrix scale-regions -S $BASEDIR/bws/*6h*AVG*.bw -R $BASEDIR/genes.bed --missingDataAsZero -p max -a 1000 -b 1000 -bs 10 -o $BASEDIR/matrices/6h_AVG_geneBody.gz

for infile in $BASEDIR/matrices/*.gz
do
  base=$(basename ${infile} .gz)
  plotHeatmap -m $infile --colorMap Reds -o $BASEDIR/figs/$base.heatmap.pdf
done
