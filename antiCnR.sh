#!/bin/bash
#SBATCH --job-name=anti_cut&run
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=80gb
#SBATCH --time=24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kld57880@uga.edu

BASEDIR="/scratch/kld57880/antiK9_5.2022"
TOOLDIR='/home/kld57880/Git2/toolbox'

###trimming, multiQC, and aligning to danio genome
module load STAR/2.7.2b-GCC-8.3.0

for file in $BASEDIR/raw/*_R*.fastq.gz;
do
  if [[ $prefix ]]; then
        base=$(basename ${first} _R1.fastq.gz)
        sh $TOOLDIR/PE_trim_and_star.sh -o $BASEDIR -n $base -m one $first $file
        prefix=
    else
        first=$file
        prefix=${file%%_*}
    fi
done

mkdir $BASEDIR/figs
module load SAMtools/1.10-iccifort-2019.5.281

###aligning to ecoli genome
curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz | gunzip -c > $BASEDIR/ecoli_refseq.fa
###note here that STAR suggests SAindex = 10 but that makes the alignment FAIL, do 8 instead
mkdir $BASEDIR/ecoli_genome
STAR --runThreadN 20 --genomeSAindexNbases 8 --runMode genomeGenerate --genomeDir $BASEDIR/ecoli_genome --genomeFastaFiles $BASEDIR/ecoli_refseq.fa

for file in $BASEDIR/trimmed/*_val_*.fq.gz;
do
  if [[ $prefix ]]; then
        base=$(basename ${first} _R1_val_1.fq.gz)
        STAR --runThreadN 20 --genomeDir $BASEDIR/ecoli_genome --outFileNamePrefix $BASEDIR/bams/"$base"_ecoli \
        --readFilesCommand zcat --readFilesIn "$first" "$file" --outSAMtype BAM SortedByCoordinate \
        --outSAMmultNmax 1 --alignEndsType EndToEnd --alignIntronMax 1 --alignMatesGapMax 2000
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

for infile in $BASEDIR/bams/*ecoli*.bam
do
  base=$(basename ${infile} Aligned.sortedByCoord.out.bam)
  java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I $infile -M $BASEDIR/bams/"$base"_dupmetrics.txt -O $BASEDIR/bams/"$base"_nodups.bam --REMOVE_DUPLICATES true
done

samtools merge -f $BASEDIR/bams/mIgG_nodups.bam $BASEDIR/bams/*IgG*hpf_nodups.bam
samtools merge -f $BASEDIR/bams/mIgG_ecoli_nodups.bam $BASEDIR/bams/*IgG*ecoli*nodups.bam

###Now we need to extract all the aligned reads in preperation for spike in normalization
module load BEDTools/2.29.2-GCC-8.3.0

for infile in $BASEDIR/bams/*nodups.bam
do
  base=$(basename ${infile} .bam)
  bedtools bamtobed -i $infile | awk -v OFS='\t' '{len = $3 - $2; print $0, len }' > $BASEDIR/bams/$base.btb.bed
done

mkdir $BASEDIR/bdgrphs

for file in $BASEDIR/bams/*nodups.btb.bed;
do
  if [[ $prefix ]]; then
        base=$(basename ${file} _nodups.btb.bed)
        sh $TOOLDIR/DNA_spike.kd.sh $file $first \
        100000 bga $BASEDIR/genome/chrNameLength.txt 1 1000 $BASEDIR/bdgrphs/"$base".norm.bga
        prefix=
    else
        first=$file
        prefix=${file%%_*}
    fi
done

###lets make  bedgraphs into bigwigs for data visualization
module load ucsc/359
mkdir $BASEDIR/bws

for infile in $BASEDIR/bdgrphs/*.norm.bga
do
  base=$(basename ${infile} .norm.bga)
  bedSort $infile $infile
  bedGraphToBigWig $infile $BASEDIR/genome/chrNameLength.txt $BASEDIR/bws/$base.bw
done

module load deepTools/3.3.1-intel-2019b-Python-3.7.4

multiBigwigSummary bins -b $BASEDIR/bws/*.bw -o $BASEDIR/bw_summ.npz -p 24
plotCorrelation -in $BASEDIR/bw_summ.npz -c spearman -p heatmap -o $BASEDIR/antiCnR_bw_summ_heatmap.pdf
plotPCA -in $BASEDIR/bw_summ.npz -o $BASEDIR/antiCnR_bw_summ_PCA.pdf

###peak calling
module load Homer/4.11-foss-2019b
mkdir $BASEDIR/peaks

for infile in $BASEDIR/bdgrphs/*norm.bga
  do base=$(basename ${infile} .norm.bga)
  cat $infile | awk '{print $1 "\t" $2 "\t" $3 "\t" "+" "\t" "+" "\t" "+"}' > $BASEDIR/peaks/$base.bgato.bed
done

for infile in $BASEDIR/peaks/*bgato.bed
  do base=$(basename ${infile} .bgato.bed)
  makeTagDirectory $BASEDIR/peaks/$base.BtB.tagdir $infile -format bed
done

for infile in $BASEDIR/peaks/*K9*.tagdir
  do base=$(basename ${infile} .BtB.tagdir)
  findPeaks $infile -style histone -minDist 1000 -gsize 1.5e9 -F 4 -i $BASEDIR/peaks/*mIgG*.tagdir -o $BASEDIR/peaks/$base.txt
done

for infile in $BASEDIR/peaks/*.txt
do
  base=$(basename ${infile} .txt)
  sed '/^#/d' $infile | awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" $8 "\t" $5 "\t" $6 "\t" $12 "\t" "-1"}' | sed 's/\.000000//g' > $BASEDIR/peaks/$base.peaks.bed
done

###ChIPr
module load ChIP-R/1.1.0-foss-2019b-Python-3.7.4

chipr -i $BASEDIR/peaks/K9abcam_24hpf*.peaks.bed -m 2 -o $BASEDIR/peaks/K9abcam_24hpfrepPeaks
chipr -i $BASEDIR/peaks/K9active_24hpf*.peaks.bed -m 2 -o $BASEDIR/peaks/K9active_24hpfrepPeaks
chipr -i $BASEDIR/peaks/K9diag_24hpf*.peaks.bed -m 2 -o $BASEDIR/peaks/K9diag_24hpfrepPeaks

chipr -i $BASEDIR/peaks/K9abcam_2.5hpf*.peaks.bed -m 2 -o $BASEDIR/peaks/K9abcam_2.5hpfrepPeaks
chipr -i $BASEDIR/peaks/K9active_2.5hpf*.peaks.bed -m 2 -o $BASEDIR/peaks/K9active_2.5hpfrepPeaks
chipr -i $BASEDIR/peaks/K9dia_2.5hpf*.peaks.bed -m 2 -o $BASEDIR/peaks/K9diag_2.5hpfrepPeaks

chipr -i $BASEDIR/peaks/K9abcam_4.5hpf*.peaks.bed -m 2 -o $BASEDIR/peaks/K9abcam_4.5hpfrepPeaks
chipr -i $BASEDIR/peaks/K9active_4.5hpf*.peaks.bed -m 2 -o $BASEDIR/peaks/K9active_4.5hpfrepPeaks
chipr -i $BASEDIR/peaks/K9diag_4.5hpf*.peaks.bed -m 2 -o $BASEDIR/peaks/K9diag_4.5hpfrepPeaks

###differenetial peaks
cat $BASEDIR/peaks/*_24hpfrepPeaks_all.bed | bedtools sort -i stdin | bedtools merge -i stdin > $BASEDIR/peaks/24hpf_MERGEDpeaks.bed
bedtools intersect -a $BASEDIR/peaks/24hpf_MERGEDpeaks.bed -b $BASEDIR/peaks/K9abcam_24hpfrepPeaks_all.bed $BASEDIR/peaks/K9active_24hpfrepPeaks_all.bed $BASEDIR/peaks/K9diag_24hpfrepPeaks_all.bed \
-wa -wb -names abcam active diag > $BASEDIR/peaks/24hpf_interPeaks.bed
sh $TOOLDIR/multi_inter.sh $BASEDIR/peaks/24hpf_interPeaks.bed $BASEDIR/peaks 24hpf

cat $BASEDIR/peaks/*_4.5hpfrepPeaks_all.bed | bedtools sort -i stdin | bedtools merge -i stdin > $BASEDIR/peaks/4.5hpf_MERGEDpeaks.bed
bedtools intersect -a $BASEDIR/peaks/4.5hpf_MERGEDpeaks.bed -b $BASEDIR/peaks/K9abcam_4.5hpfrepPeaks_all.bed $BASEDIR/peaks/K9active_4.5hpfrepPeaks_all.bed $BASEDIR/peaks/K9diag_4.5hpfrepPeaks_all.bed \
-wa -wb -names abcam active diag > $BASEDIR/peaks/4.5hpf_interPeaks.bed
sh $TOOLDIR/multi_inter.sh $BASEDIR/peaks/4.5hpf_interPeaks.bed $BASEDIR/peaks 4.5hpf

cat $BASEDIR/peaks/*_2.5hpfrepPeaks_all.bed | bedtools sort -i stdin | bedtools merge -i stdin > $BASEDIR/peaks/2.5hpf_MERGEDpeaks.bed
bedtools intersect -a $BASEDIR/peaks/2.5hpf_MERGEDpeaks.bed -b $BASEDIR/peaks/K9abcam_2.5hpfrepPeaks_all.bed $BASEDIR/peaks/K9active_2.5hpfrepPeaks_all.bed $BASEDIR/peaks/K9diag_2.5hpfrepPeaks_all.bed \
-wa -wb -names abcam active diag > $BASEDIR/peaks/2.5hpf_interPeaks.bed
sh $TOOLDIR/multi_inter.sh $BASEDIR/peaks/2.5hpf_interPeaks.bed $BASEDIR/peaks 2.5hpf

####peak annotation####
curl -s ftp://ftp.ensembl.org/pub/release-98/gtf/danio_rerio/Danio_rerio.GRCz11.98.gtf.gz | gunzip -c > $BASEDIR/refann.gtf
mkdir $BASEDIR/peaks/ann

for infile in $BASEDIR/peaks/*all.bed
do
  base=$( basename ${infile} repPeaks_all.bed)
  annotatePeaks.pl $infile danRer11 -gtf $BASEDIR/refann.gtf > $BASEDIR/peaks/ann/$base.maskann.txt
done

for infile in $BASEDIR/*.bed
do
  base=$( basename ${infile} .bed)
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

###average bw files
module load deepTools/3.3.1-intel-2019b-Python-3.7.4

bigwigCompare -b1 $BASEDIR/bws/K9abcam_24hpf_1.bw -b2 $BASEDIR/bws/K9abcam_24hpf_2.bw --operation add --scaleFactors 0.333:0.333 -bs 10 -p 20 -o $BASEDIR/bws/K9abcam_24hpf_rep1rep2.bw
bigwigCompare -b1 $BASEDIR/bws/K9abcam_24hpf_rep1rep2.bw -b2 $BASEDIR/bws/K9abcam_24hpf_3.bw --operation add --scaleFactors 1:0.333 -bs 10 -p 20 -o $BASEDIR/bws/K9abcam_24hpf_AVG.bw
bigwigCompare -b1 $BASEDIR/bws/K9active_24hpf_1.bw -b2 $BASEDIR/bws/K9active_24hpf_2.bw --operation add --scaleFactors 0.333:0.333 -bs 10 -p 20 -o $BASEDIR/bws/K9active_24hpf_rep1rep2.bw
bigwigCompare -b1 $BASEDIR/bws/K9active_24hpf_rep1rep2.bw -b2 $BASEDIR/bws/K9active_24hpf_3.bw --operation add --scaleFactors 1:0.333 -bs 10 -p 20 -o $BASEDIR/bws/K9active_24hpf_AVG.bw
bigwigCompare -b1 $BASEDIR/bws/K9diag_24hpf_1.bw -b2 $BASEDIR/bws/K9diag_24hpf_2.bw --operation add --scaleFactors 0.333:0.333 -bs 10 -p 20 -o $BASEDIR/bws/K9diag_24hpf_rep1rep2.bw
bigwigCompare -b1 $BASEDIR/bws/K9diag_24hpf_rep1rep2.bw -b2 $BASEDIR/bws/K9diag_24hpf_3.bw --operation add --scaleFactors 1:0.333 -bs 10 -p 20 -o $BASEDIR/bws/K9diag_24hpf_AVG.bw

bigwigCompare -b1 $BASEDIR/bws/K9abcam_4.5hpf_1.bw -b2 $BASEDIR/bws/K9abcam_4.5hpf_2.bw --operation add --scaleFactors 0.333:0.333 -bs 10 -p 20 -o $BASEDIR/bws/K9abcam_4.5hpf_rep1rep2.bw
bigwigCompare -b1 $BASEDIR/bws/K9abcam_4.5hpf_rep1rep2.bw -b2 $BASEDIR/bws/K9abcam_4.5hpf_3.bw --operation add --scaleFactors 1:0.333 -bs 10 -p 20 -o $BASEDIR/bws/K9abcam_4.5hpf_AVG.bw
bigwigCompare -b1 $BASEDIR/bws/K9active_4.5hpf_1.bw -b2 $BASEDIR/bws/K9active_4.5hpf_2.bw --operation add --scaleFactors 0.333:0.333 -bs 10 -p 20 -o $BASEDIR/bws/K9active_4.5hpf_rep1rep2.bw
bigwigCompare -b1 $BASEDIR/bws/K9active_4.5hpf_rep1rep2.bw -b2 $BASEDIR/bws/K9active_4.5hpf_3.bw --operation add --scaleFactors 1:0.333 -bs 10 -p 20 -o $BASEDIR/bws/K9active_4.5hpf_AVG.bw
bigwigCompare -b1 $BASEDIR/bws/K9diag_4.5hpf_1.bw -b2 $BASEDIR/bws/K9diag_4.5hpf_2.bw --operation add --scaleFactors 0.333:0.333 -bs 10 -p 20 -o $BASEDIR/bws/K9diag_4.5hpf_rep1rep2.bw
bigwigCompare -b1 $BASEDIR/bws/K9diag_4.5hpf_rep1rep2.bw -b2 $BASEDIR/bws/K9diag_4.5hpf_3.bw --operation add --scaleFactors 1:0.333 -bs 10 -p 20 -o $BASEDIR/bws/K9diag_4.5hpf_AVG.bw

bigwigCompare -b1 $BASEDIR/bws/K9abcam_2.5hpf_1.bw -b2 $BASEDIR/bws/K9abcam_2.5hpf_2.bw --operation add --scaleFactors 0.333:0.333 -bs 10 -p 20 -o $BASEDIR/bws/K9abcam_2.5hpf_rep1rep2.bw
bigwigCompare -b1 $BASEDIR/bws/K9abcam_2.5hpf_rep1rep2.bw -b2 $BASEDIR/bws/K9abcam_2.5hpf_3.bw --operation add --scaleFactors 1:0.333 -bs 10 -p 20 -o $BASEDIR/bws/K9abcam_2.5hpf_AVG.bw
bigwigCompare -b1 $BASEDIR/bws/K9active_2.5hpf_1.bw -b2 $BASEDIR/bws/K9active_2.5hpf_2.bw --operation add --scaleFactors 0.333:0.333 -bs 10 -p 20 -o $BASEDIR/bws/K9active_2.5hpf_rep1rep2.bw
bigwigCompare -b1 $BASEDIR/bws/K9active_2.5hpf_rep1rep2.bw -b2 $BASEDIR/bws/K9active_2.5hpf_3.bw --operation add --scaleFactors 1:0.333 -bs 10 -p 20 -o $BASEDIR/bws/K9active_2.5hpf_AVG.bw
bigwigCompare -b1 $BASEDIR/bws/K9dia_2.5hpf_1.bw -b2 $BASEDIR/bws/K9dia_2.5hpf_2.bw --operation add --scaleFactors 0.333:0.333 -bs 10 -p 20 -o $BASEDIR/bws/K9diag_2.5hpf_rep1rep2.bw
bigwigCompare -b1 $BASEDIR/bws/K9diag_2.5hpf_rep1rep2.bw -b2 $BASEDIR/bws/K9dia_2.5hpf_3.bw --operation add --scaleFactors 1:0.333 -bs 10 -p 20 -o $BASEDIR/bws/K9diag_2.5hpf_AVG.bw

###blacklist intersect
for infile in $BASEDIR/peaks/*all.bed
do
  base=$( basename ${infile} repPeaks_all.bed)
  bedtools intersect -a $infile -b /scratch/kld57880/TC_final/peaks/blacklist.bed -v > $BASEDIR/peaks/"$base"_final.bed
done

###going to do a couple peak intersections for some analysis
ml BEDTools

bedtools intersect -a $BASEDIR/peaks/K9abcam_24hpf_final.bed -b $BASEDIR/peaks/K9abcam_4.5hpf_final.bed -wa > $BASEDIR/peaks/K9abcam_24_4.5_peaks.bed
bedtools intersect -a $BASEDIR/peaks/K9abcam_24_4.5_peaks.bed -b $BASEDIR/peaks/K9abcam_2.5hpf_final.bed -wa > $BASEDIR/peaks/K9abcam_24_4.5_2.5_peaks.bed
bedtools intersect -a $BASEDIR/peaks/K9abcam_4.5hpf_final.bed -b $BASEDIR/peaks/K9abcam_2.5hpf_final.bed -wa > $BASEDIR/peaks/K9abcam_4.5_2.5_peaks.bed

bedtools intersect -a $BASEDIR/peaks/K9abcam_2.5hpf_final.bed -b $BASEDIR/peaks/K9abcam_24_4.5_2.5_peaks.bed  -v > $BASEDIR/peaks/K9abcam_2.5only_peaks.bed
bedtools intersect -a $BASEDIR/peaks/K9abcam_4.5hpf_final.bed -b $BASEDIR/peaks/K9abcam_4.5_2.5_peaks.bed -v > $BASEDIR/peaks/K9abcam_4.5only_peaks.bed
bedtools intersect -a $BASEDIR/peaks/K9abcam_24hpf_final.bed -b $BASEDIR/peaks/K9abcam_24_4.5_peaks.bed -v > $BASEDIR/peaks/K9abcam_24only_peaks.bed

###make some pics
mkdir $BASEDIR/matrices

computeMatrix reference-point -S $BASEDIR/bws/K9abcam_24hpf_AVG.bw -R $BASEDIR/peaks/K9abcam_4.5hpfrepPeaks_all.bed --referencePoint center -p max -a 10000 -b 10000 -bs 10 --missingDataAsZero -o $BASEDIR/matrices/abcam24hpf.gz
computeMatrix reference-point -S $BASEDIR/bws/K9abcam_4.5hpf_AVG.bw -R $BASEDIR/peaks/K9abcam_4.5hpfrepPeaks_all.bed --referencePoint center -p max -a 10000 -b 10000 -bs 10 --missingDataAsZero -o $BASEDIR/matrices/abcam4.5hpf.gz
computeMatrix reference-point -S $BASEDIR/bws/K9abcam_2.5hpf_AVG.bw -R $BASEDIR/peaks/K9abcam_2.5hpfrepPeaks_all.bed --referencePoint center -p max -a 10000 -b 10000 -bs 10 --missingDataAsZero -o $BASEDIR/matrices/abcam2.5hpf.gz
computeMatrix reference-point -S $BASEDIR/bws/K9abcam_24hpf_AVG.bw $BASEDIR/bws/K9active_24hpf_AVG.bw $BASEDIR/bws/K9diag_24hpf_AVG.bw -R $BASEDIR/peaks/24hpf_a*.bed $BASEDIR/peaks/24hpf_d*.bed --referencePoint center -p max -a 10000 -b 10000 -bs 10 --missingDataAsZero -o $BASEDIR/matrices/24hpfPeaks.gz
computeMatrix reference-point -S $BASEDIR/bws/K9abcam_4.5hpf_AVG.bw $BASEDIR/bws/K9active_4.5hpf_AVG.bw $BASEDIR/bws/K9diag_4.5hpf_AVG.bw -R $BASEDIR/peaks/4.5hpf_a*.bed $BASEDIR/peaks/4.5hpf_d*.bed --referencePoint center -p max -a 10000 -b 10000 -bs 10 --missingDataAsZero -o $BASEDIR/matrices/4.5hpfPeaks.gz
computeMatrix reference-point -S $BASEDIR/bws/K9abcam_2.5hpf_AVG.bw $BASEDIR/bws/K9active_2.5hpf_AVG.bw $BASEDIR/bws/K9diag_2.5hpf_AVG.bw -R $BASEDIR/peaks/2.5hpf_a*.bed $BASEDIR/peaks/2.5hpf_d*.bed --referencePoint center -p max -a 10000 -b 10000 -bs 10 --missingDataAsZero -o $BASEDIR/matrices/2.5hpfPeaks.gz

for infile in $BASEDIR/matrices/*.gz
do
  base=$(basename ${infile} .gz)
  plotHeatmap -m $infile --colorMap Blues --legendLocation none --whatToShow "plot and heatmap" --regionsLabel "$base" -o $BASEDIR/figs/"$base"_heatmap.pdf
done

###pericentromere intersect
mkdir $BASEDIR/peaks/peric2

for infile in $BASEDIR/peaks/K9abcam*.bed
do
  base=$(basename ${infile} .bed)
  bedtools intersect -a $infile -b /scratch/kld57880/TC_final/peric/centromeres_sloppy2.bed -wa > $BASEDIR/peaks/peric2/$base.peric.bed
  bedtools intersect -a $infile -b /scratch/kld57880/TC_final/peric/centromeres_sloppy2.bed -v > $BASEDIR/peaks/peric2/$base.non_peric.bed
done

for infile in /scratch/kld57880/TC_final/peric/null/NULLcoords2*.bed
do
  base=$(basename ${infile} .bed)
  bedtools intersect -a $BASEDIR/peaks/int/pre_ZGA_final.bed -b $infile -wa > $BASEDIR/peaks/peric2/pre_ZGA.$base.bed
  bedtools intersect -a $BASEDIR/peaks/int/post_ZGA_final.bed -b $infile -wa > $BASEDIR/peaks/peric2/post_ZGA.$base.bed
done

##making pics to show abcam is better
computeMatrix reference-point -S $BASEDIR/bws/K9abcam_2.5hpf_AVG.bw $BASEDIR/bws/K9active_2.5hpf_AVG.bw $BASEDIR/bws/K9diag_2.5hpf_AVG.bw -R $BASEDIR/peaks/other_anti/2.5hpf_abcam_active_diag.bed $BASEDIR/peaks/2.5hpf_abcam_only.bed --referencePoint center -p max -a 10000 -b 10000 -bs 10 --missingDataAsZero -o $BASEDIR/matrices/2.5hpf_anti_comp2.gz
computeMatrix reference-point -S $BASEDIR/bws/K9abcam_4.5hpf_AVG.bw $BASEDIR/bws/K9active_4.5hpf_AVG.bw $BASEDIR/bws/K9diag_4.5hpf_AVG.bw -R $BASEDIR/peaks/other_anti/4.5hpf_abcam_active_diag.bed $BASEDIR/peaks/4.5hpf_abcam_only.bed --referencePoint center -p max -a 10000 -b 10000 -bs 10 --missingDataAsZero -o $BASEDIR/matrices/4.5hpf_anti_comp2.gz
computeMatrix reference-point -S $BASEDIR/bws/K9abcam_24hpf_AVG.bw $BASEDIR/bws/K9active_24hpf_AVG.bw $BASEDIR/bws/K9diag_24hpf_AVG.bw -R $BASEDIR/peaks/other_anti/24hpf_abcam_active_diag.bed $BASEDIR/peaks/24hpf_abcam_only.bed --referencePoint center -p max -a 10000 -b 10000 -bs 10 --missingDataAsZero -o $BASEDIR/matrices/24hpf_anti_comp2.gz

for infile in $BASEDIR/matrices/*comp2.gz
do
  base=$(basename ${infile} .gz)
  plotHeatmap -m $infile --colorMap Blues -o $BASEDIR/figs/"$base"_heatmap.pdf
done
