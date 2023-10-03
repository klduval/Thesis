#!/bin/bash
##########################################
##Help
##########################################
Help()
{
  #Display Help
  echo "TC_peaks is used to compare peaks across timepoints/samples"
  echo
  echo "SYNTAX: $0 [options] peakfile1.bed peakfile2.bed . . . peakfile6.bed"
  echo
  echo "OPTIONS:"
  echo "-o      output directory, required"
  echo
  echo "-m      determines mode of analysis"
  echo "  -m all    performs analysis on all peaks (default)"
  echo "  -m merged performs analysis on peaks merged across timepoints"
  echo
  echo "-n      {non-functional rn} indicates that you want non-linear output in addition to linear"
  echo "        for example, peaks that are in timepoints 1 and 3, but not 2"
  echo "        default is NOT to provide this information"
  echo
  echo "-t      timepoint/sample labels in order, no spaces"
  echo "        IMPORTANT - each timepoint label goes after it's own flag"
  echo "        for example: -t 2hpf -t 3hpf -t 4hpf"
  echo "        If no timepoint labels are provided, times will be labelled 1 2 3  ... 6"
  echo
  echo "INPUT FILES: at least 3 BED files are listed after options"
  echo "      Currently only works on ChIPr output formatted files"
  echo "      Currently only works for up to 6 timepoints"
  echo
  echo "-h      print this help message"
  echo
}
#########################################
###this is me trying to write a script for timecourse peak analysis
#####Process the input options
if [ "$#" == 0 ]; then
    echo "Error"
    Help
    exit
fi

while getopts ":ho:m:t:" option; do
  case "$option" in
    h) # display Help
      Help
      exit
      ;;
    o) #setting output directory
      OUTDIR=$OPTARG
      ;;
    :) echo "Option -$OPTARG requires an argument"
      echo
      Help
      exit
      ;;
    m) #setting mode
      mode=$OPTARG
      ;;
    t) times+=("$OPTARG");;
      #...
    \?) # invalid option
      echo "Invalid option -$OPTARG"
      echo
      echo "Help:"
      Help
      exit
      ;;
  esac
done

peak_file1=${@:$OPTIND:1}
peak_file2=${@:$OPTIND+1:1}
peak_file3=${@:$OPTIND+2:1}
peak_file4=${@:$OPTIND+3:1}
peak_file5=${@:$OPTIND+4:1}
peak_file6=${@:$OPTIND+5:1}
peak_file7=${@:$OPTIND+6:1}

echo
echo "OUTDIR is $OUTDIR"

echo
echo "PEAK FILES GIVEN"
for infile in $peak_file1 $peak_file2 $peak_file3 $peakfile_4 $peak_file5 $peak_file6
do
  if [ -f "$infile" ] ; then
  echo $(basename $infile)
  fi
done

if [ -f "$peak_file1" ] && [ -f "$peak_file2" ] && [ -f "$peak_file3" ] && [ -f "$peak_file4" ] && [ -f "$peak_file5" ] && [ -f "$peak_file6" ] ; then
  num_timepoints=6
  echo "...There are $num_timepoints timepoints"
elif [ -f "$peak_file1" ] && [ -f "$peak_file2" ] && [ -f "$peak_file3" ] && [ -f "$peak_file4" ] && [ -f "$peak_file5" ] ; then
  num_timepoints=5
  echo "...There are $num_timepoints timepoints"
elif [ -f "$peak_file1" ] && [ -f "$peak_file2" ] && [ -f "$peak_file3" ] && [ -f "$peak_file4" ] ; then
  num_timepoints=4
  echo "...There are $num_timepoints timepoints"
elif [ -f "$peak_file1" ] && [ -f "$peak_file2" ] && [ -f "$peak_file3" ] ; then
  num_timepoints=3
  echo "...There are $num_timepoints timepoints"
elif [ -f "$peak_file7"] ; then
  echo "...Error: Too many timepoints or file input error. Quitting"
  exit
else
  echo "...Error: Not enough timepoints. Quitting"
  exit
fi

###Detecting the names of timepoints based on input labels with the -t function
echo
if [ "${#times[@]}" -gt 0 ]; then
  echo "TIMEPOINTS GIVEN"
  echo "There are ${#times[@]} labels given:"
  for val in "${times[@]}" ; do
    echo " - $val"
  done
else
  echo "No timepoint labels given"
  echo "Timepoints will be labelled 1 to $num_timepoints"
fi
echo

if [ "${#times[@]}" -gt 0 ] && [ "${#times[@]}" != "$num_timepoints" ]; then
echo "Error: number of peak files and timepoint labels do not match. Quitting"
exit
fi

###Assining our timepoint labels to variables to make this easier later
if [ "${#times[@]}" -gt 0 ]; then
  pf1_time="${times[0]}"
  pf2_time="${times[1]}"
  pf3_time="${times[2]}"
  pf4_time="${times[3]}"
  pf5_time="${times[4]}"
  pf6_time="${times[5]}"
else
  pf1_time=1
  pf2_time=2
  pf3_time=3
  pf4_time=4
  pf5_time=5
  pf6_time=6
fi

###Printing labels assigned to peak files
echo "TIMEPOINTS ASSIGNED TO PEAK FILES"
if [ -f "$peak_file1" ] ; then
echo "$(basename $peak_file1) = $pf1_time"
fi
if [ -f "$peak_file2" ] ; then
echo "$(basename $peak_file2) = $pf2_time"
fi
if [ -f "$peak_file3" ] ; then
echo "$(basename $peak_file3) = $pf3_time"
fi
if [ -f "$peak_file4" ] ; then
echo "$(basename $peak_file4) = $pf4_time"
fi
if [ -f "$peak_file5" ] ; then
echo "$(basename $peak_file5) = $pf5_time"
fi
if [ -f "$peak_file6" ] ; then
echo "$(basename $peak_file6) = $pf6_time"
fi
echo

echo "Determining Analysis Mode"
if [ -z "$mode" ] ; then
  analysis_mode=ALL
  echo " Analysis Mode is $analysis_mode by default"
elif [ $mode == 'merged' ] ; then
  analysis_mode=MERGED
  echo "  Analysis Mode is $analysis_mode"
elif [ $mode == 'all' ] ; then
  analysis_mode=ALL
  echo "  Analysis Mode is $analysis_mode"
else
  analysis_mode=ALL
  echo "  Analysis Mode could not be determined so Analysis Mode is $analysis_mode"
fi

echo
echo "All data input checks complete."
echo
echo "Starting Analysis"
echo "...Loading BEDtools"
module load BEDTools/2.29.2-GCC-8.3.0
results=results_"$analysis_mode"_TCpeaks
echo "...Making results summary file in $OUTDIR/$results.txt"
touch $OUTDIR/"$results".txt
sumfile=$OUTDIR/"$results".txt

if [ $analysis_mode == 'ALL' ] ; then
cat $peak_file1 $peak_file2 $peak_file3 $peak_file4 $peak_file5 $peak_file6 | bedtools sort -i stdin > $OUTDIR/ALLpeaks.bed
echo "The number of ALL peaks is" >> $sumfile
wc -l $OUTDIR/ALLpeaks.bed >> $sumfile
echo "...ALL peaks have been added to one file"
elif [ $analysis_mode == 'MERGED' ] ; then
cat $peak_file1 $peak_file2 $peak_file3 $peak_file4 $peak_file5 $peak_file6 | bedtools sort -i stdin | bedtools merge -i stdin > $OUTDIR/MERGEDpeaks.bed
echo "The number of MERGED peaks is" >> $sumfile
wc -l $OUTDIR/MERGEDpeaks.bed >> $sumfile
echo "...Peaks have been MERGED"
else
echo "...Error: peak add/merge failed. Quitting"
exit
fi

###now going to intersect these peaks with each timepoint peaks
if [ $analysis_mode == 'ALL' ] ; then
  if [ -f "$peak_file1" ] ; then
    bedtools intersect -a $OUTDIR/ALLpeaks.bed -b $peak_file1 -wa  > $OUTDIR/aPeaks_int.$pf1_time.bed
  fi
  if [ -f "$peak_file2" ] ; then
    bedtools intersect -a $OUTDIR/ALLpeaks.bed -b $peak_file2 -wa  > $OUTDIR/aPeaks_int.$pf2_time.bed
  fi
  if [ -f "$peak_file3" ] ; then
    bedtools intersect -a $OUTDIR/ALLpeaks.bed -b $peak_file3 -wa  > $OUTDIR/aPeaks_int.$pf3_time.bed
  fi
  if [ -f "$peak_file4" ] ; then
    bedtools intersect -a $OUTDIR/ALLpeaks.bed -b $peak_file4 -wa  > $OUTDIR/aPeaks_int.$pf4_time.bed
  fi
  if [ -f "$peak_file5" ] ; then
    bedtools intersect -a $OUTDIR/ALLpeaks.bed -b $peak_file5 -wa  > $OUTDIR/aPeaks_int.$pf5_time.bed
  fi
  if [ -f "$peak_file6" ] ; then
    bedtools intersect -a $OUTDIR/ALLpeaks.bed -b $peak_file6 -wa  > $OUTDIR/aPeaks_int.$pf6_time.bed
  fi
  echo "...Peaks have been intersected"
elif [ $analysis_mode == 'MERGED' ] ; then
  if [ -f "$peak_file1" ] ; then
    bedtools intersect -a $OUTDIR/MERGEDpeaks.bed -b $peak_file1 -wa  > $OUTDIR/mPeaks_int.$pf1_time.bed
  fi
  if [ -f "$peak_file2" ] ; then
    bedtools intersect -a $OUTDIR/MERGEDpeaks.bed -b $peak_file2 -wa  > $OUTDIR/mPeaks_int.$pf2_time.bed
  fi
  if [ -f "$peak_file3" ] ; then
    bedtools intersect -a $OUTDIR/MERGEDpeaks.bed -b $peak_file3 -wa  > $OUTDIR/mPeaks_int.$pf3_time.bed
  fi
  if [ -f "$peak_file4" ] ; then
    bedtools intersect -a $OUTDIR/MERGEDpeaks.bed -b $peak_file4 -wa  > $OUTDIR/mPeaks_int.$pf4_time.bed
  fi
  if [ -f "$peak_file5" ] ; then
    bedtools intersect -a $OUTDIR/MERGEDpeaks.bed -b $peak_file5 -wa  > $OUTDIR/mPeaks_int.$pf5_time.bed
  fi
  if [ -f "$peak_file6" ] ; then
    bedtools intersect -a $OUTDIR/MERGEDpeaks.bed -b $peak_file6 -wa  > $OUTDIR/mPeaks_int.$pf6_time.bed
  fi
  echo "...Peaks have been intersected"
else
  echo "...Error: peak intersection failed. Quitting"
  exit
fi

###now making the individual peak intersection counts files
if [ $analysis_mode == 'MERGED' ] ; then
  for infile in $OUTDIR/mPeaks_int*.bed
  do
    base=$(basename ${infile} .bed)
    bedtools intersect -a $OUTDIR/MERGEDpeaks.bed -b $infile -wa -c -f 1 -r > $OUTDIR/$base.count.temp.bed
  done
  echo "...Peak intersections have been counted"
elif [ $analysis_mode == 'ALL' ] ; then
  for infile in $OUTDIR/aPeaks_int*.bed
  do
    base=$(basename ${infile} .bed)
    bedtools intersect -a $OUTDIR/ALLpeaks.bed -b $infile -wa -c -f 1 -r > $OUTDIR/$base.count.temp.bed
  done
  echo "...Peak intersections have been counted"
else
  echo "...Error: intersection counts failed. Quitting"
  exit
fi

###Now we need to make the summary file
if [ $analysis_mode == 'MERGED' ]; then
  if [ $num_timepoints == 3 ]; then
    cut -f4 $OUTDIR/mPeaks_int."$pf2_time".count.temp.bed | paste $OUTDIR/mPeaks_int."$pf1_time".count.temp.bed - > $OUTDIR/mPeaks_1-2count.temp.bed
    cut -f4 $OUTDIR/mPeaks_int."$pf3_time".count.temp.bed | paste $OUTDIR/mPeaks_1-2count.temp.bed - > $OUTDIR/mPeaks_total.count.bed
    echo "...Summary count file made in $OUTDIR/mPeaks_total.count.bed"
  elif [ $num_timepoints == 4 ]; then
    cut -f4 $OUTDIR/mPeaks_int."$pf2_time".count.temp.bed | paste $OUTDIR/mPeaks_int."$pf1_time".count.temp.bed - > $OUTDIR/mPeaks_1-2count.temp.bed
    cut -f4 $OUTDIR/mPeaks_int."$pf3_time".count.temp.bed | paste $OUTDIR/mPeaks_1-2count.temp.bed - > $OUTDIR/mPeaks_1-3.count.temp.bed
    cut -f4 $OUTDIR/mPeaks_int."$pf4_time".count.temp.bed | paste $OUTDIR/mPeaks_1-3.count.temp.bed - > $OUTDIR/mPeaks_total.count.bed
    echo "...Summary count file made in $OUTDIR/mPeaks_total.count.bed"
  elif [ $num_timepoints == 5 ]; then
    cut -f4 $OUTDIR/mPeaks_int."$pf2_time".count.temp.bed | paste $OUTDIR/mPeaks_int."$pf1_time".count.temp.bed - > $OUTDIR/mPeaks_1-2count.temp.bed
    cut -f4 $OUTDIR/mPeaks_int."$pf3_time".count.temp.bed | paste $OUTDIR/mPeaks_1-2count.temp.bed - > $OUTDIR/mPeaks_1-3.count.temp.bed
    cut -f4 $OUTDIR/mPeaks_int."$pf4_time".count.temp.bed | paste $OUTDIR/mPeaks_1-3.count.temp.bed - > $OUTDIR/mPeaks_1-4.count.temp.bed
    cut -f4 $OUTDIR/mPeaks_int."$pf5_time".count.temp.bed | paste $OUTDIR/mPeaks_1-4.count.temp.bed - > $OUTDIR/mPeaks_total.count.bed
    echo "...Summary count file made in $OUTDIR/mPeaks_total.count.bed"
  elif [ $num_timepoints == 6 ]; then
    cut -f4 $OUTDIR/mPeaks_int."$pf2_time".count.temp.bed | paste $OUTDIR/mPeaks_int."$pf1_time".count.temp.bed - > $OUTDIR/mPeaks_1-2count.temp.bed
    cut -f4 $OUTDIR/mPeaks_int."$pf3_time".count.temp.bed | paste $OUTDIR/mPeaks_1-2count.temp.bed - > $OUTDIR/mPeaks_1-3.count.temp.bed
    cut -f4 $OUTDIR/mPeaks_int."$pf4_time".count.temp.bed | paste $OUTDIR/mPeaks_1-3.count.temp.bed - > $OUTDIR/mPeaks_1-4.count.temp.bed
    cut -f4 $OUTDIR/mPeaks_int."$pf5_time".count.temp.bed | paste $OUTDIR/mPeaks_1-4.count.temp.bed - > $OUTDIR/mPeaks_1-5.count.temp.bed
    cut -f4 $OUTDIR/mPeaks_int."$pf6_time".count.temp.bed | paste $OUTDIR/mPeaks_1-5.count.temp.bed - > $OUTDIR/mPeaks_total.count.bed
    echo "...Summary count file made in $OUTDIR/mPeaks_total.count.bed"
  fi
elif [ $analysis_mode == 'ALL' ]; then
  if [ $num_timepoints == 3 ]; then
    cut -f1,2,3,10 $OUTDIR/aPeaks_int."$pf1_time".count.temp.bed | paste - > $OUTDIR/aPeaks_1count.temp.bed
    cut -f10 $OUTDIR/aPeaks_int."$pf2_time".count.temp.bed | paste $OUTDIR/aPeaks_1count.temp.bed - > $OUTDIR/aPeaks_1-2count.temp.bed
    cut -f10 $OUTDIR/aPeaks_int."$pf3_time".count.temp.bed | paste $OUTDIR/aPeaks_1-2count.temp.bed - > $OUTDIR/aPeaks_total.count.bed
    echo "...Summary count file made in $OUTDIR/aPeaks_total.count.bed"
  elif [ $num_timepoints == 4 ]; then
    cut -f1,2,3,10 $OUTDIR/aPeaks_int."$pf1_time".count.temp.bed | paste - > $OUTDIR/aPeaks_1count.temp.bed
    cut -f10 $OUTDIR/aPeaks_int."$pf2_time".count.temp.bed | paste $OUTDIR/aPeaks_1count.temp.bed - > $OUTDIR/aPeaks_1-2count.temp.bed
    cut -f10 $OUTDIR/aPeaks_int."$pf3_time".count.temp.bed | paste $OUTDIR/aPeaks_1-2count.temp.bed - > $OUTDIR/aPeaks_1-3count.temp.bed
    cut -f10 $OUTDIR/aPeaks_int."$pf4_time".count.temp.bed | paste $OUTDIR/aPeaks_1-3count.temp.bed - > $OUTDIR/aPeaks_total.count.bed
    echo "...Summary count file made in $OUTDIR/aPeaks_total.count.bed"
  elif [ $num_timepoints == 5 ]; then
    cut -f1,2,3,10 $OUTDIR/aPeaks_int."$pf1_time".count.temp.bed | paste - > $OUTDIR/aPeaks_1count.temp.bed
    cut -f10 $OUTDIR/aPeaks_int."$pf2_time".count.temp.bed | paste $OUTDIR/aPeaks_1count.temp.bed - > $OUTDIR/aPeaks_1-2count.temp.bed
    cut -f10 $OUTDIR/aPeaks_int."$pf3_time".count.temp.bed | paste $OUTDIR/aPeaks_1-2count.temp.bed - > $OUTDIR/aPeaks_1-3count.temp.bed
    cut -f10 $OUTDIR/aPeaks_int."$pf4_time".count.temp.bed | paste $OUTDIR/aPeaks_1-3count.temp.bed - > $OUTDIR/aPeaks_1-4count.temp.bed
    cut -f10 $OUTDIR/aPeaks_int."$pf5_time".count.temp.bed | paste $OUTDIR/aPeaks_1-4count.temp.bed - > $OUTDIR/aPeaks_total.count.bed
    echo "...Summary count file made in $OUTDIR/aPeaks_total.count.bed"
  elif [ $num_timepoints == 6 ]; then
    cut -f1,2,3,10 $OUTDIR/aPeaks_int."$pf1_time".count.temp.bed | paste - > $OUTDIR/aPeaks_1count.temp.bed
    cut -f10 $OUTDIR/aPeaks_int."$pf2_time".count.temp.bed | paste $OUTDIR/aPeaks_1count.temp.bed - > $OUTDIR/aPeaks_1-2count.temp.bed
    cut -f10 $OUTDIR/aPeaks_int."$pf3_time".count.temp.bed | paste $OUTDIR/aPeaks_1-2count.temp.bed - > $OUTDIR/aPeaks_1-3count.temp.bed
    cut -f10 $OUTDIR/aPeaks_int."$pf4_time".count.temp.bed | paste $OUTDIR/aPeaks_1-3count.temp.bed - > $OUTDIR/aPeaks_1-4count.temp.bed
    cut -f10 $OUTDIR/aPeaks_int."$pf5_time".count.temp.bed | paste $OUTDIR/aPeaks_1-4count.temp.bed - > $OUTDIR/aPeaks_1-5count.temp.bed
    cut -f10 $OUTDIR/aPeaks_int."$pf6_time".count.temp.bed | paste $OUTDIR/aPeaks_1-5count.temp.bed - > $OUTDIR/aPeaks_total.count.bed
    echo "...Summary count file made in $OUTDIR/aPeaks_total.count.bed"
  fi
else
  echo "...Error: summary of counts failed. Quitting"
  exit
fi

echo "...Deleting temporary files"
rm $OUTDIR/*Peaks_int.*bed
rm $OUTDIR/*.temp.bed

if [ $analysis_mode == 'MERGED' ]; then
  totalcountfile="$OUTDIR/mPeaks_total.count.bed"
elif [ $analysis_mode == 'ALL' ]; then
  totalcountfile="$OUTDIR/aPeaks_total.count.bed"
fi

subdir="$analysis_mode"peak_groups
echo "...Making a subdirectory for grouped peak files in $OUTDIR/$subdir"
mkdir $OUTDIR/$subdir

###now lets make the nexted awk commands to count the groups of peaks
if [ $num_timepoints == 6 ]; then
  awk '$4!=0 {print}' $totalcountfile | awk '$5!=0 {print}' | awk '$6!=0 {print}' | awk '$7!=0 {print}' | awk '$8!=0 {print}' | awk '$9!=0 {print}' > $OUTDIR/$subdir/peaks_"$pf1_time"thru"$pf6_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5!=0 {print}' | awk '$6!=0 {print}' | awk '$7!=0 {print}' | awk '$8!=0 {print}' | awk '$9!=0 {print}' > $OUTDIR/$subdir/peaks_"$pf2_time"thru"$pf6_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5==0 {print}' | awk '$6!=0 {print}' | awk '$7!=0 {print}' | awk '$8!=0 {print}' | awk '$9!=0 {print}' > $OUTDIR/$subdir/peaks_"$pf3_time"thru"$pf6_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5==0 {print}' | awk '$6==0 {print}' | awk '$7!=0 {print}' | awk '$8!=0 {print}' | awk '$9!=0 {print}' > $OUTDIR/$subdir/peaks_"$pf4_time"thru"$pf6_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5==0 {print}' | awk '$6==0 {print}' | awk '$7==0 {print}' | awk '$8!=0 {print}' | awk '$9!=0 {print}' > $OUTDIR/$subdir/peaks_"$pf5_time"thru"$pf6_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5==0 {print}' | awk '$6==0 {print}' | awk '$7==0 {print}' | awk '$8==0 {print}' | awk '$9!=0 {print}' > $OUTDIR/$subdir/peaks_"$pf6_time"ONLY.bed
  awk '$4!=0 {print}' $totalcountfile | awk '$5!=0 {print}' | awk '$6!=0 {print}' | awk '$7!=0 {print}' | awk '$8!=0 {print}' | awk '$9==0 {print}' > $OUTDIR/$subdir/peaks_"$pf1_time"thru"$pf5_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5!=0 {print}' | awk '$6!=0 {print}' | awk '$7!=0 {print}' | awk '$8!=0 {print}' | awk '$9==0 {print}' > $OUTDIR/$subdir/peaks_"$pf2_time"thru"$pf5_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5==0 {print}' | awk '$6!=0 {print}' | awk '$7!=0 {print}' | awk '$8!=0 {print}' | awk '$9==0 {print}' > $OUTDIR/$subdir/peaks_"$pf3_time"thru"$pf5_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5==0 {print}' | awk '$6==0 {print}' | awk '$7!=0 {print}' | awk '$8!=0 {print}' | awk '$9==0 {print}' > $OUTDIR/$subdir/peaks_"$pf4_time"thru"$pf5_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5==0 {print}' | awk '$6==0 {print}' | awk '$7==0 {print}' | awk '$8!=0 {print}' | awk '$9==0 {print}' > $OUTDIR/$subdir/peaks_"$pf5_time"ONLY.bed
  awk '$4!=0 {print}' $totalcountfile | awk '$5!=0 {print}' | awk '$6!=0 {print}' | awk '$7!=0 {print}' | awk '$8==0 {print}' | awk '$9==0 {print}' > $OUTDIR/$subdir/peaks_"$pf1_time"thru"$pf4_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5!=0 {print}' | awk '$6!=0 {print}' | awk '$7!=0 {print}' | awk '$8==0 {print}' | awk '$9==0 {print}' > $OUTDIR/$subdir/peaks_"$pf2_time"thru"$pf4_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5==0 {print}' | awk '$6!=0 {print}' | awk '$7!=0 {print}' | awk '$8==0 {print}' | awk '$9==0 {print}' > $OUTDIR/$subdir/peaks_"$pf3_time"thru"$pf4_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5==0 {print}' | awk '$6==0 {print}' | awk '$7!=0 {print}' | awk '$8==0 {print}' | awk '$9==0 {print}' > $OUTDIR/$subdir/peaks_"$pf4_time"ONLY.bed
  awk '$4!=0 {print}' $totalcountfile | awk '$5!=0 {print}' | awk '$6!=0 {print}' | awk '$7==0 {print}' | awk '$8==0 {print}' | awk '$9==0 {print}' > $OUTDIR/$subdir/peaks_"$pf1_time"thru"$pf3_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5!=0 {print}' | awk '$6!=0 {print}' | awk '$7==0 {print}' | awk '$8==0 {print}' | awk '$9==0 {print}' > $OUTDIR/$subdir/peaks_"$pf2_time"thru"$pf3_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5==0 {print}' | awk '$6!=0 {print}' | awk '$7==0 {print}' | awk '$8==0 {print}' | awk '$9==0 {print}' > $OUTDIR/$subdir/peaks_"$pf3_time"ONLY.bed
  awk '$4!=0 {print}' $totalcountfile | awk '$5!=0 {print}' | awk '$6==0 {print}' | awk '$7==0 {print}' | awk '$8==0 {print}' | awk '$9==0 {print}' > $OUTDIR/$subdir/peaks_"$pf1_time"thru"$pf2_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5!=0 {print}' | awk '$6==0 {print}' | awk '$7==0 {print}' | awk '$8==0 {print}' | awk '$9==0 {print}' > $OUTDIR/$subdir/peaks_"$pf2_time"ONLY.bed
  awk '$4!=0 {print}' $totalcountfile | awk '$5==0 {print}' | awk '$6==0 {print}' | awk '$7==0 {print}' | awk '$8==0 {print}' | awk '$9==0 {print}' > $OUTDIR/$subdir/peaks_"$pf1_time"ONLY.bed
  echo "Peak group counts:" >> $sumfile
  wc -l $OUTDIR/$subdir/peaks* >> $sumfile
  echo "...Peak groups counted and summarised in $sumfile"
elif [ $num_timepoints == 5 ]; then
  awk '$4!=0 {print}' $totalcountfile | awk '$5!=0 {print}' | awk '$6!=0 {print}' | awk '$7!=0 {print}' | awk '$8!=0 {print}' > $OUTDIR/$subdir/peaks_"$pf1_time"thru"$pf5_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5!=0 {print}' | awk '$6!=0 {print}' | awk '$7!=0 {print}' | awk '$8!=0 {print}' > $OUTDIR/$subdir/peaks_"$pf2_time"thru"$pf5_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5==0 {print}' | awk '$6!=0 {print}' | awk '$7!=0 {print}' | awk '$8!=0 {print}' > $OUTDIR/$subdir/peaks_"$pf3_time"thru"$pf5_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5==0 {print}' | awk '$6==0 {print}' | awk '$7!=0 {print}' | awk '$8!=0 {print}' > $OUTDIR/$subdir/peaks_"$pf4_time"thru"$pf5_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5==0 {print}' | awk '$6==0 {print}' | awk '$7==0 {print}' | awk '$8!=0 {print}' > $OUTDIR/$subdir/peaks_"$pf5_time"ONLY.bed
  awk '$4!=0 {print}' $totalcountfile | awk '$5!=0 {print}' | awk '$6!=0 {print}' | awk '$7!=0 {print}' | awk '$8==0 {print}' > $OUTDIR/$subdir/peaks_"$pf1_time"thru"$pf4_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5!=0 {print}' | awk '$6!=0 {print}' | awk '$7!=0 {print}' | awk '$8==0 {print}' > $OUTDIR/$subdir/peaks_"$pf2_time"thru"$pf4_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5==0 {print}' | awk '$6!=0 {print}' | awk '$7!=0 {print}' | awk '$8==0 {print}' > $OUTDIR/$subdir/peaks_"$pf3_time"thru"$pf4_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5==0 {print}' | awk '$6==0 {print}' | awk '$7!=0 {print}' | awk '$8==0 {print}' > $OUTDIR/$subdir/peaks_"$pf4_time"ONLY.bed
  awk '$4!=0 {print}' $totalcountfile | awk '$5!=0 {print}' | awk '$6!=0 {print}' | awk '$7==0 {print}' | awk '$8==0 {print}' > $OUTDIR/$subdir/peaks_"$pf1_time"thru"$pf3_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5!=0 {print}' | awk '$6!=0 {print}' | awk '$7==0 {print}' | awk '$8==0 {print}' > $OUTDIR/$subdir/peaks_"$pf2_time"thru"$pf3_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5==0 {print}' | awk '$6!=0 {print}' | awk '$7==0 {print}' | awk '$8==0 {print}' > $OUTDIR/$subdir/peaks_"$pf3_time"ONLY.bed
  awk '$4!=0 {print}' $totalcountfile | awk '$5!=0 {print}' | awk '$6==0 {print}' | awk '$7==0 {print}' | awk '$8==0 {print}' > $OUTDIR/$subdir/peaks_"$pf1_time"thru"$pf2_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5!=0 {print}' | awk '$6==0 {print}' | awk '$7==0 {print}' | awk '$8==0 {print}' > $OUTDIR/$subdir/peaks_"$pf2_time"ONLY.bed
  awk '$4!=0 {print}' $totalcountfile | awk '$5==0 {print}' | awk '$6==0 {print}' | awk '$7==0 {print}' | awk '$8==0 {print}' > $OUTDIR/$subdir/peaks_"$pf1_time"ONLY.bed
  echo "Peak group counts:" >> $sumfile
  wc -l $OUTDIR/$subdir/peaks* >> $sumfile
  echo "...Peak groups counted and summarised in $sumfile"
elif [ $num_timepoints == 4 ]; then
  awk '$4!=0 {print}' $totalcountfile | awk '$5!=0 {print}' | awk '$6!=0 {print}' | awk '$7!=0 {print}' > $OUTDIR/$subdir/peaks_"$pf1_time"thru"$pf4_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5!=0 {print}' | awk '$6!=0 {print}' | awk '$7!=0 {print}' > $OUTDIR/$subdir/peaks_"$pf2_time"thru"$pf4_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5==0 {print}' | awk '$6!=0 {print}' | awk '$7!=0 {print}' > $OUTDIR/$subdir/peaks_"$pf3_time"thru"$pf4_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5==0 {print}' | awk '$6==0 {print}' | awk '$7!=0 {print}' > $OUTDIR/$subdir/peaks_"$pf4_time"ONLY.bed
  awk '$4!=0 {print}' $totalcountfile | awk '$5!=0 {print}' | awk '$6!=0 {print}' | awk '$7==0 {print}' > $OUTDIR/$subdir/peaks_"$pf1_time"thru"$pf3_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5!=0 {print}' | awk '$6!=0 {print}' | awk '$7==0 {print}' > $OUTDIR/$subdir/peaks_"$pf2_time"thru"$pf3_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5==0 {print}' | awk '$6!=0 {print}' | awk '$7==0 {print}' > $OUTDIR/$subdir/peaks_"$pf3_time"ONLY.bed
  awk '$4!=0 {print}' $totalcountfile | awk '$5!=0 {print}' | awk '$6==0 {print}' | awk '$7==0 {print}' > $OUTDIR/$subdir/peaks_"$pf1_time"thru"$pf2_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5!=0 {print}' | awk '$6==0 {print}' | awk '$7==0 {print}' > $OUTDIR/$subdir/peaks_"$pf2_time"ONLY.bed
  awk '$4!=0 {print}' $totalcountfile | awk '$5==0 {print}' | awk '$6==0 {print}' | awk '$7==0 {print}' > $OUTDIR/$subdir/peaks_"$pf1_time"ONLY.bed
  echo "Peak group counts:" >> $sumfile
  wc -l $OUTDIR/$subdir/peaks* >> $sumfile
  echo "...Peak groups counted and summarised in $sumfile"
elif [ $num_timepoints == 3 ]; then
  awk '$4!=0 {print}' $totalcountfile | awk '$5!=0 {print}' | awk '$6!=0 {print}' > $OUTDIR/$subdir/peaks_"$pf1_time"thru"$pf3_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5!=0 {print}' | awk '$6!=0 {print}' > $OUTDIR/$subdir/peaks_"$pf2_time"thru"$pf3_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5==0 {print}' | awk '$6!=0 {print}' > $OUTDIR/$subdir/peaks_"$pf3_time"ONLY.bed
  awk '$4!=0 {print}' $totalcountfile | awk '$5!=0 {print}' | awk '$6==0 {print}' > $OUTDIR/$subdir/peaks_"$pf1_time"thru"$pf2_time".bed
  awk '$4==0 {print}' $totalcountfile | awk '$5!=0 {print}' | awk '$6==0 {print}' > $OUTDIR/$subdir/peaks_"$pf2_time"ONLY.bed
  awk '$4!=0 {print}' $totalcountfile | awk '$5==0 {print}' | awk '$6==0 {print}' > $OUTDIR/$subdir/peaks_"$pf1_time"ONLY.bed
  echo "Peak group counts:" >> $sumfile
  wc -l $OUTDIR/$subdir/peaks* >> $sumfile
  echo "...Peak groups counted and summarised in $sumfile"
else
  echo "...Error: Peak groups could not be determined. Quitting"
  exit
fi


echo
echo "...$0 completed running succesffuly."
