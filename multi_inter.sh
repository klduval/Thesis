#!/bin/bash

# Check if the input filename, output directory, and base name are provided as command-line arguments
if [[ $# -lt 2 ]]; then
  echo "Error: Input filename or output directory is missing."
  echo "Usage: ./multi_inter.sh <input_filename> <output_directory> [base_name]"
  exit 1
fi

# Read the input filename and output directory from the command-line arguments
input_filename="$1"
output_directory="$2"

# Read the optional base name
base_name=""
if [[ $# -eq 3 ]]; then
  base_name="$3"
fi

# Create the output directory if it doesn't exist
mkdir -p "$output_directory"

# Read the input file using awk
awk -F'\t' -v outdir="$output_directory" -v basename="$base_name" '{
  coordinates = $1"\t"$2"\t"$3
  if ($4 != "") {
    antibodies[$4]++
    if (coordinates in interval_dict)
      interval_dict[coordinates] = interval_dict[coordinates] "," $4
    else
      interval_dict[coordinates] = $4
  }
} END {
  for (coordinates in interval_dict) {
    split(interval_dict[coordinates], arr, ",")
    count = length(arr)
    output_filename = ""

    if (count == 1) {
      output_filename = arr[1] "ONLY.bed"
    } else {
      asort(arr)  # Sort the antibodies array to ensure consistent filename order
      output_filename = arr[1]
      for (i = 2; i <= count; i++) {
        if (arr[i] != arr[i-1])  # Add only unique antibodies to the filename
          output_filename = output_filename "_" arr[i]
      }
      output_filename = output_filename ".bed"
    }

    if (output_filename != "") {
      if (basename != "") {
        output_filename = basename "_" output_filename
      }
      output_path = outdir "/" output_filename
      print coordinates >> output_path
      print "Interval", coordinates, "written to", output_path
    }
  }
}' "$input_filename"

# Combine files with matching names into a single file and delete the individual files
for file in "$output_directory"/*ONLY.bed; do
  matching_file="${file%ONLY.bed}.bed"
  combined_file="${file%ONLY.bed}_only.bed"
  cat "$file" "$matching_file" > "$combined_file"
  echo "Combined $file and $matching_file into $combined_file"
  rm "$file" "$matching_file"
  echo "Deleted $file and $matching_file"
done
