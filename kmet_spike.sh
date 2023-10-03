#!/bin/bash
####incorporating kmet spike-in panel from epicypher
####to the our version of the Henikoff spike-in normaliation
###original spike-in script from henikoff lab:
#https://raw.githubusercontent.com/Henikoff/Cut-and-Run/master/spike_in_calibration.csh.

##########################################
Help()
{
  #Display Help
  echo "This is a script to normalize CUT&RUN sequencing that uses the kmet spike in from Epicypher"
  echo "If you didn't use the kmet spike-in when performing CUT&RUN, this won't do anything!"
  echo
  echo "USAGE: $0 output_dir base_name trimmed_1.fq.gz trimmed_2.fq.gz alignment.bam output_format(bg|bga|d) genome_chr_lens"
  echo
  echo "~~~~~~~~"
  echo "output_dir       output directory"
  echo "base_nam         base name for output files, required"
  echo "trimmed_1        first trimmed sequencing file corresponding to the sample in question, required"
  echo "trimmed_2        second trimmed sequencing file corresponding to the sample in question, required"
  echo "alignment.bam    alignment bam file we are trying to normalize, required"
  echo "output_format(bg|bga|d)"
  echo "                         bg = bedgraph format "
  echo "                         bga = bedgraph format, including zeros (reccomended) "
  echo "                         d = reports the depth at each genome position with 1-based coordinates. "
  echo "genome_chr_lens   direction to file listing chromosome names and lengths (usually chrNameLength.txt in STAR genome file)"
  echo
  echo "-h      print this help message"
  echo
  echo "good luck with your analysis! go science!"
  echo
}
#########################################

if [ "$#" -ne 7 ]; then
    echo "error!"
    Help
    exit 1
fi

OUTDIR=${1}
base=${2}
trimmed_1=${3}
trimmed_2=${4}
bam=${5}
report=${6}
genome_len=${7}

echo "**Running KMET Spike-in Normalization**"
echo
echo "  $OUTDIR is the output directory"
echo "  Normalizing $base sample"
echo

echo "Checking if $base trimmed read files are real and Non-empty"
if [ -s $trimmed_1 ] && [ -s $trimmed_2 ]; then
    echo "  $trimmed_1 and $trimmed_2 found and not empty"
else
    echo "  $trimmed_1 and $trimmed_2 was found to either not exist or be empty"
    echo "exiting"
    exit 0
fi

echo
echo "Counting KMET barcodes"

###first need to count the kmet barcodes
for barcode in TTCGCGCGTAACGACGTACCGT CGCGATACGACCGCGTTACGCG CGACGTTAACGCGTTTCGTACG CGCGACTATCGCGCGTAACGCG CCGTACGTCGTGTCGAACGACG CGATACGCGTTGGTACGCGTAA TAGTTCGCGACACCGTTCGTCG TCGACGCGTAAACGGTACGTCG TTATCGCGTCGCGACGGACGTA CGATCGTACGATAGCGTACCGA CGCATATCGCGTCGTACGACCG ACGTTCGACCGCGGTCGTACGA ACGATTCGACGATCGTCGACGA CGATAGTCGCGTCGCACGATCG CGCCGATTACGTGTCGCGCGTA ATCGTACCGCGCGTATCGGTCG CGTTCGAACGTTCGTCGACGAT TCGCGATTACGATGTCGCGCGA ACGCGAATCGTCGACGCGTATA CGCGATATCACTCGACGCGATA CGCGAAATTCGTATACGCGTCG CGCGATCGGTATCGGTACGCGC GTGATATCGCGTTAACGTCGCG TATCGCGCGAAACGACCGTTCG CCGCGCGTAATGCGCGACGTTA CCGCGATACGACTCGTTCGTCG GTCGCGAACTATCGTCGATTCG CCGCGCGTATAGTCCGAGCGTA CGATACGCCGATCGATCGTCGG CCGCGCGATAAGACGCGTAACG CGATTCGACGGTCGCGACCGTA TTTCGACGCGTCGATTCGGCGA ;
	do
	zgrep -c $barcode $trimmed_1 >> $OUTDIR/kmet_1."$base".tmp.txt
done

kmet_1=`awk '{s+=$1}END{print s}' $OUTDIR/kmet_1."$base".tmp.txt`
echo "  kmet count 1 is $kmet_1"

for barcode in TTCGCGCGTAACGACGTACCGT CGCGATACGACCGCGTTACGCG CGACGTTAACGCGTTTCGTACG CGCGACTATCGCGCGTAACGCG CCGTACGTCGTGTCGAACGACG CGATACGCGTTGGTACGCGTAA TAGTTCGCGACACCGTTCGTCG TCGACGCGTAAACGGTACGTCG TTATCGCGTCGCGACGGACGTA CGATCGTACGATAGCGTACCGA CGCATATCGCGTCGTACGACCG ACGTTCGACCGCGGTCGTACGA ACGATTCGACGATCGTCGACGA CGATAGTCGCGTCGCACGATCG CGCCGATTACGTGTCGCGCGTA ATCGTACCGCGCGTATCGGTCG CGTTCGAACGTTCGTCGACGAT TCGCGATTACGATGTCGCGCGA ACGCGAATCGTCGACGCGTATA CGCGATATCACTCGACGCGATA CGCGAAATTCGTATACGCGTCG CGCGATCGGTATCGGTACGCGC GTGATATCGCGTTAACGTCGCG TATCGCGCGAAACGACCGTTCG CCGCGCGTAATGCGCGACGTTA CCGCGATACGACTCGTTCGTCG GTCGCGAACTATCGTCGATTCG CCGCGCGTATAGTCCGAGCGTA CGATACGCCGATCGATCGTCGG CCGCGCGATAAGACGCGTAACG CGATTCGACGGTCGCGACCGTA TTTCGACGCGTCGATTCGGCGA ;
	do
  zgrep -c $barcode $trimmed_2 >> $OUTDIR/kmet_2."$base".tmp.txt
done

kmet_2=`awk '{s+=$1}END{print s}' $OUTDIR/kmet_2."$base".tmp.txt`
echo "  kmet count 2 is $kmet_2"

kmet_total=`echo "$kmet_1 + $kmet_2" | bc -l`

echo "KMET total barcode count is $kmet_total"

if [ -s $OUTDIR/kmet_stats.txt ]; then
    echo "$base   $kmet_total" >> $OUTDIR/kmet_stats.txt
else
    echo "Sample  KMET_Count" >> $OUTDIR/kmet_stats.txt
    echo "                  " >> $OUTDIR/kmet_stats.txt
    echo "$base   $kmet_total" >> $OUTDIR/kmet_stats.txt
fi

echo "KMET counts being saved in $OUTDIR/kmet_stats.txt"

###making btb files
echo
echo "Extracting aligned reads for normalization"
echo "Loading BEDTools"
module load BEDTools

bedtools bamtobed -i $bam | awk -v OFS='\t' '{len = $3 - $2; print $0, len }' > $OUTDIR/$base.btb.bed
genome_bed=`echo $OUTDIR/"$base".btb.bed`

echo "Extracted reads saved in $genome_bed"

######now for the spike in
echo
echo "Checking if $base bam-to-bed file is real and Non-empty"
if [ -s $genome_bed ] ; then
    echo "${genome_bed} found and not empty"
else
    echo "${genome_bed} was found to either not exist or be empty"
    exit 0
fi

echo "Checking if $genome_len file is real and Non-empty"
if [ -s $genome_len ] ; then
    echo "$genome_len found and not empty"
else
    echo "${genome_len} was found to either not exist or be empty"
    exit 0
fi

echo
echo "Generating Scale Factor"
#MPR is 1000000/total danio reads so total danio reads is getting multiplied by spike in reads
if [ -s $genome_bed ] ; then
        read_count=`wc -l $genome_bed | awk '{print $1}'`
        echo "Read count is $read_count"
        MPR=`echo "1000000/$read_count" | bc -l`
        echo "MPR is $MPR"
        kmet_scale=`echo "1000000 / $kmet_total" | bc -l`
        echo "KMET scale is $kmet_scale"
        echo "Multiplying MPR by KMET scale"
        scale_factor=`echo "$MPR * $kmet_scale" | bc -l`
        echo "Scale Factor is $scale_factor"
else
              echo "No Genome File given"
              echo "exiting"
              exit 0
fi

#Select fragments within the length range, assumes fragment length is in column 6 of the bed file
echo
echo "Generating Fragment Lengths"
cat $genome_bed | awk -v min=10 -v max=1000 '$7 >= min &&  $7 <= max {print $0}' > $OUTDIR/frag_length.temp.bed

echo "Generating Genome Coverage"
bedtools genomecov -${report} -scale ${scale_factor} -i $OUTDIR/frag_length.temp.bed -g ${genome_len} > ${OUTDIR}/${base}_kmet.${report}

echo "SCRIPT DONE"

rm $OUTDIR/frag_length.temp.bed
rm $OUTDIR/kmet*"$base"*.tmp.txt
