#!/bin/bash
#SBATCH --job-name=K9_timecourse
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=120gb
#SBATCH --time=48:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kld57880@uga.edu

BASEDIR="/scratch/kld57880/TC_final"
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

###aligning to ecoli genome
curl -s https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz | gunzip -c > $BASEDIR/ecoli_refseq.fa
###note here that STAR suggests SAindex = 10 but that makes the alignment FAIL, do 8 instead
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
module load SAMtools/1.10-iccifort-2019.5.281

for infile in $BASEDIR/bams/*q1.bam
do
  base=$(basename ${infile} _q1.bam)
  java -jar $EBROOTPICARD/picard.jar MarkDuplicates -I $infile -M $BASEDIR/bams/"$base"_dupmetrics.txt -O $BASEDIR/bams/"$base"_nodups.bam --REMOVE_DUPLICATES true
done

###going to merge all my IgG samples together to create uniformity in peak calling later
###since I'm about to normalize for library size that should account for increase in read numbers
samtools merge -f $BASEDIR/bams/mIgG_nodups.bam $BASEDIR/bams/*IgG*[1-3]_nodups.bam
samtools merge -f $BASEDIR/bams/mIgG_ecoli_nodups.bam $BASEDIR/bams/*IgG*ecoli*nodups.bam

#Now we need to extract all the aligned reads in preperation for spike in normalization
module load BEDTools/2.29.2-GCC-8.3.0

for infile in $BASEDIR/bams/*nodups.bam
do
  base=$(basename ${infile} .bam)
  bedtools bamtobed -i $infile | awk -v OFS='\t' '{len = $3 - $2; print $0, len }' > $BASEDIR/bams/$base.btb.bed
done

##spike in normalization
mkdir $BASEDIR/bdgrphs

for file in $BASEDIR/bams/*.btb.bed;
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
