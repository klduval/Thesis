#!/bin/bash
#SBATCH --job-name=K9_timecourse
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=100gb
#SBATCH --time=12:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kld57880@uga.edu

BASEDIR="/scratch/kld57880/TC_final"

###peak calling
module load Homer/4.11-foss-2019b
mkdir $BASEDIR/peaks

for infile in $BASEDIR/bdgrphs/*.norm.bga
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

module load ChIP-R/1.1.0-foss-2019b-Python-3.7.4

chipr -i $BASEDIR/peaks/2hpf_K9_1.peaks.bed $BASEDIR/peaks/2hpf_K9_2.peaks.bed -m 2 -o $BASEDIR/peaks/2hpf_K9_repPeaks
chipr -i $BASEDIR/peaks/2.5hpf_K9_1.peaks.bed $BASEDIR/peaks/2.5hpf_K9_2.peaks.bed $BASEDIR/peaks/2.5hpf_K9_3.peaks.bed -m 2 -o $BASEDIR/peaks/2.5hpf_K9_repPeaks
chipr -i $BASEDIR/peaks/3hpf_K9_1.peaks.bed $BASEDIR/peaks/3hpf_K9_2.peaks.bed $BASEDIR/peaks/3hpf_K9_3.peaks.bed -m 2 -o $BASEDIR/peaks/3hpf_K9_repPeaks
chipr -i $BASEDIR/peaks/3.5hpf_K9_1.peaks.bed $BASEDIR/peaks/3.5hpf_K9_2.peaks.bed $BASEDIR/peaks/3.5hpf_K9_3.peaks.bed -m 2 -o $BASEDIR/peaks/3.5hpf_K9_repPeaks
chipr -i $BASEDIR/peaks/4hpf_K9_1.peaks.bed $BASEDIR/peaks/4hpf_K9_2.peaks.bed $BASEDIR/peaks/4hpf_K9_3.peaks.bed -m 2 -o $BASEDIR/peaks/4hpf_K9_repPeaks
chipr -i $BASEDIR/peaks/4.5hpf_K9_1.peaks.bed $BASEDIR/peaks/4.5hpf_K9_2.peaks.bed $BASEDIR/peaks/4.5hpf_K9_3.peaks.bed -m 2 -o $BASEDIR/peaks/4.5hpf_K9_repPeaks

###make a blacklist file
findPeaks $BASEDIR/peaks/tagdirs/IgG.meh.tagdir -style factor -o $BASEDIR/peaks/IgG.txt
sed '/^#/d' $BASEDIR/peaks/IgG.txt | awk '{print $2 "\t" $3 "\t" $4 "\t" $1 "\t" "1" "\t" $5 "\t" $6 "\t" $12 "\t" "-1"}' > $BASEDIR/peaks/blacklist.bed

ml BEDTools

###intersect the peaks with the blacklist file to make sure we aren't looking at sticky regions before this step
for infile in $BASEDIR/peaks/*all.bed
do
  base=$( basename ${infile} _repPeaks_all.bed)
  bedtools intersect -a $infile -b $BASEDIR/peaks/blacklist.bed -v > $BASEDIR/peaks/"$base"_final.bed
done

####peak annotation####
curl -s ftp://ftp.ensembl.org/pub/release-98/gtf/danio_rerio/Danio_rerio.GRCz11.98.gtf.gz | gunzip -c > $BASEDIR/refann.gtf
mkdir $BASEDIR/peaks/ann

for infile in $BASEDIR/peaks/*final.bed
do
  base=$( basename ${infile} final.bed)
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
	base=$(basename ${infile} _.TEann.txt)
	sed 's:;:\t:g' $infile > $BASEDIR/peaks/ann/$base.trim.txt
done

for infile in $BASEDIR/peaks/ann/*.trim.txt
do
	base=$(basename ${infile} .trim.txt)
	sed 's:"::g' $infile | sed 's:gene_id::g' | sed 's:transcript_id::g' | sed 's:family_id::g' | sed 's:class_id::g' | sed 's:Cluster::g' | sed 's:TE_classification::g' > $BASEDIR/peaks/ann/$base.txt
done

#####want to split peaks with pre-ZGA seeding and those without
bedtools intersect -a $BASEDIR/peaks/4.5hpf_K9_final.bed -b $BASEDIR/peaks/2.5hpf_K9_final.bed -v > $BASEDIR/peaks/4.5hpf_NOseed_peaks.bed
bedtools intersect -a $BASEDIR/peaks/4.5hpf_K9_final.bed -b $BASEDIR/peaks/2.5hpf_K9_final.bed -wa > $BASEDIR/peaks/4.5hpf_seed_peaks.bed

for infile in $BASEDIR/peaks/*seed*.bed
do
  base=$( basename ${infile} final.bed)
  annotatePeaks.pl $infile danRer11 -gtf $BASEDIR/refann.gtf > $BASEDIR/peaks/ann/$base.maskann.txt
done

for infile in $BASEDIR/peaks/ann/4.5*seed*maskann.txt
do
  base=$(basename ${infile} .maskann.txt)
  awk -F'\t' 'sqrt($10*$10) <=1000' $infile > $BASEDIR/peaks/ann/$base.1000bp_ann.txt
done

for infile in $BASEDIR/peaks/ann/4.5*seed*maskann.txt
do
  base=$(basename ${infile} .maskann.txt)
  awk -F'\t' 'sqrt($10*$10) >=1000' $infile | awk '{print $2 "\t" $3 "\t" $4 }' > $BASEDIR/peaks/ann/$base.MOREthan1000bp.bed
done

for infile in $BASEDIR/peaks/ann/4.5*seed*.MOREthan1000bp.bed
do
  base=$( basename ${infile} .MOREthan1000bp.bed)
  bedtools intersect -a $infile -b /work/mglab/kld/TEfiles_Chang_etal/TEann.gtf -F 0.50 -wo > $BASEDIR/peaks/ann/$base.TEann.txt
done

for infile in $BASEDIR/peaks/ann/4.5*seed*TEann.txt
do
	base=$(basename ${infile} _.TEann.txt)
	sed 's:;:\t:g' $infile > $BASEDIR/peaks/ann/$base.trim.txt
done

for infile in $BASEDIR/peaks/ann/4.5*seed*.trim.txt
do
	base=$(basename ${infile} .trim.txt)
	sed 's:"::g' $infile | sed 's:gene_id::g' | sed 's:transcript_id::g' | sed 's:family_id::g' | sed 's:class_id::g' | sed 's:Cluster::g' | sed 's:TE_classification::g' > $BASEDIR/peaks/ann/$base.txt
done

####now going to intersect peaks with the TC_pericentromeres
mkdir $BASEDIR/peric/peak_int

for infile in $BASEDIR/peaks/*final.bed
do
  base=$(basename ${infile} _final.bed)
  bedtools intersect -a $infile -b $BASEDIR/peric/centromeres_sloppy.bed -wa > $BASEDIR/peric/peak_int/$base.peric.bed
  bedtools intersect -a $infile -b $BASEDIR/peric/centNULL_total.bed -wa > $BASEDIR/peric/peak_int/$base.null.bed
done

for infile in $BASEDIR/peric/peak_int/*.bed
do
  base=$( basename ${infile} final.bed)
  annotatePeaks.pl $infile danRer11 -gtf $BASEDIR/refann.gtf > $BASEDIR/peric/peak_int/$base.maskann.txt
done

for infile in $BASEDIR/peric/peak_int/*maskann.txt
do
  base=$(basename ${infile} .maskann.txt)
  awk -F'\t' 'sqrt($10*$10) <=1000' $infile > $BASEDIR/peric/peak_int/$base.1000bp_ann.txt
done

for infile in $BASEDIR/peric/peak_int/*maskann.txt
do
  base=$(basename ${infile} .maskann.txt)
  awk -F'\t' 'sqrt($10*$10) >=1000' $infile | awk '{print $2 "\t" $3 "\t" $4 }' > $BASEDIR/peric/peak_int/$base.MOREthan1000bp.bed
done

for infile in $BASEDIR/peric/peak_int/*.MOREthan1000bp.bed
do
  base=$( basename ${infile} .MOREthan1000bp.bed)
  bedtools intersect -a $infile -b /work/mglab/kld/TEfiles_Chang_etal/TEann.gtf -F 0.50 -wo > $BASEDIR/peric/peak_int/$base.TEann.txt
done

for infile in $BASEDIR/peric/peak_int/*TEann.txt
do
	base=$(basename ${infile} _.TEann.txt)
	sed 's:;:\t:g' $infile > $BASEDIR/peric/peak_int/$base.trim.txt
done

for infile in $BASEDIR/peric/peak_int/*.trim.txt
do
	base=$(basename ${infile} .trim.txt)
	sed 's:"::g' $infile | sed 's:gene_id::g' | sed 's:transcript_id::g' | sed 's:family_id::g' | sed 's:class_id::g' | sed 's:Cluster::g' | sed 's:TE_classification::g' > $BASEDIR/peric/peak_int/$base.txt
done

for infile in $BASEDIR/peaks/*seed*.bed
do
  base=$(basename ${infile} .bed)
  bedtools intersect -a $infile -b $BASEDIR/peric/centromeres_sloppy.bed -wa > $BASEDIR/peric/peak_int/$base.peric.bed
  bedtools intersect -a $infile -b $BASEDIR/peric/centNULL_total.bed -wa > $BASEDIR/peric/peak_int/$base.null.bed
done
