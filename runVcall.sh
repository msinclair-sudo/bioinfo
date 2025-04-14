#!/bin/bash
set -e  # Abort on any error

########################################################################
# Setup: Create an intermediates directory and logging functions.
########################################################################
intermdir="./intermediates"
mkdir -p "$intermdir"
LOGFILE="$intermdir/pipeline.log"

log_message() {
  # Log a message with a timestamp (print to stdout and append to log file)
  echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$LOGFILE"
}

usage() {
  echo "Usage: $0 -b <BAM file> -r <reference FASTA file> -o <output VCF file> [-l <BED file>] [-t <threads>]"
  echo "  -b  Path to the unsorted BAM file."
  echo "  -r  Path to the exhaustive reference FASTA file."
  echo "  -o  Path to the output VCF file."
  echo "  -l  (Optional) BED file with regions to call variants."
  echo "  -t  (Optional) Number of threads to use (default: 4)."
  exit 1
}

# Set default threads.
THREADS=4
while getopts ":b:r:o:l:t:" opt; do
  case $opt in
    b) BAM="$OPTARG" ;;
    r) FASTA="$OPTARG" ;;
    o) OUTVCF="$OPTARG" ;;
    l) BED="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    *) usage ;;
  esac
done
if [ -z "$BAM" ] || [ -z "$FASTA" ] || [ -z "$OUTVCF" ]; then
  usage
fi

# Create output directory if it does not exist.
output_dir=$(dirname "$OUTVCF")
if [ ! -d "$output_dir" ]; then
  mkdir -p "$output_dir"
  log_message "Created output directory: $output_dir"
else
  log_message "Output directory already exists: $output_dir"
fi

log_message "Starting pipeline with parameters:"
log_message "  Unsorted BAM: $BAM"
log_message "  Reference FASTA: $FASTA"
log_message "  Output VCF: $OUTVCF"
if [ -n "$BED" ]; then
  log_message "  BED file: $BED"
fi
log_message "  Threads: $THREADS"
log_message "-------------------------------------------"

########################################################################
# Step 1: Process the BAM file: sort, index, and extract first 5 contigs.
########################################################################
log_message "Step 1: Sorting BAM file with $THREADS threads..."
sorted_bam="$intermdir/$(basename "$BAM" .bam).sorted.bam"
if [ ! -s "$sorted_bam" ]; then
  samtools sort -@ "$THREADS" -o "$sorted_bam" "$BAM"
  log_message "Sorted BAM saved as: $sorted_bam"
else
  log_message "Sorted BAM already exists; skipping sorting."
fi

log_message "Indexing sorted BAM file..."
if [ ! -f "${sorted_bam}.bai" ]; then
  samtools index "$sorted_bam"
  log_message "BAM indexing completed."
else
  log_message "BAM index already exists; skipping indexing."
fi

# Extract only the first word (the contig name up to the first whitespace) from the BAM header.
log_message "Extracting first 5 contig names from BAM header..."
bam_contigs=$(samtools view -H "$sorted_bam" | grep '^@SQ' | head -n5 | sed 's/.*SN:\([^[:space:]]*\).*/\1/')
log_message "BAM contigs:"
echo "$bam_contigs" | tee -a "$LOGFILE"

# Determine expected naming: if the first contig starts with "chr", expected_prefix is "chr"; otherwise, numeric.
bam_first=$(echo "$bam_contigs" | head -n1)
if [[ $bam_first == chr* ]]; then
  expected_prefix="chr"
else
  expected_prefix=""
fi
log_message "Expected naming prefix (from BAM): [$expected_prefix]"
log_message "-------------------------------------------"

########################################################################
# Step 2: Process Reference FASTA (Filter, Reformat, and Transform)
########################################################################
FASTA_PROCESSED="$intermdir/$(basename "$FASTA" | sed 's/\.[^.]*$//').processed.fasta"
if [ ! -s "$FASTA_PROCESSED" ]; then
  log_message "Step 2: Processing reference FASTA (Filter, Reformat, Transform)"
  
  # 2a: Filter the reference to include only sequences with "GRCh37.p13 Primary Assembly" (header only)
  # and only those that start with "NC_" (RefSeq primary assemblies).
  filtered_fasta="$intermdir/$(basename "$FASTA" | sed 's/\.[^.]*$//').filtered.fasta"
  log_message "Filtering reference FASTA for primary assemblies..."
  seqkit grep -r -i -n -p "GRCh37\.p13 Primary Assembly" "$FASTA" | seqkit grep -r -n -p "^NC_" > "$filtered_fasta"
  log_message "Filtered FASTA saved as: $filtered_fasta"
  
  # 2b: Reformat the filtered FASTA to a fixed line width.
  reformatted_fasta="$intermdir/$(basename "$filtered_fasta" | sed 's/\.[^.]*$//').reformatted.fasta"
  log_message "Reformatting filtered FASTA to fixed line width..."
  seqtk seq -l60 "$filtered_fasta" > "$reformatted_fasta"
  log_message "Formatted FASTA saved as: $reformatted_fasta"
  
  # 2c: Extract the first 5 contig names from the reformatted FASTA.
  fa_contigs=$(grep '^>' "$reformatted_fasta" | head -n5 | sed 's/^>//; s/ .*//')
  log_message "FASTA contigs (before transformation):"
  echo "$fa_contigs" | tee -a "$LOGFILE"
  
  # 2d: Transform headers if they are in RefSeq NC_ format.
  fast_first=$(echo "$fa_contigs" | head -n1)
  if [[ $fast_first == NC_* ]]; then
    log_message "Detected RefSeq NC_ format in FASTA headers. Converting to desired format."
    transformed_fasta="$intermdir/$(basename "$reformatted_fasta" | sed 's/\.[^.]*$//').transformed.fasta"
    if [[ "$expected_prefix" == "chr" ]]; then
      log_message "Expected prefix is 'chr'. Converting to 'chrN' format."
      sed -E 's/^>NC_[0]*([0-9]+)\..*/>chr\1/; t; s/.*[cC]hromosome[[:space:]]+X.*/>chr23/; s/.*[cC]hromosome[[:space:]]+Y.*/>chr24/' "$reformatted_fasta" > "$transformed_fasta"
    else
      log_message "No expected prefix. Converting to numeric format."
      sed -E 's/^>NC_[0]*([0-9]+)\..*/>\1/; t; s/.*[cC]hromosome[[:space:]]+X.*/>23/; s/.*[cC]hromosome[[:space:]]+Y.*/>24/' "$reformatted_fasta" > "$transformed_fasta"
    fi
    log_message "Transformed FASTA saved as: $transformed_fasta"
    FASTA_TEMP="$transformed_fasta"
  else
    log_message "No RefSeq NC_ format detected; assuming headers are already in desired format."
    FASTA_TEMP="$reformatted_fasta"
  fi
  
  # 2e: Double-check that contig names are extracted after transformation.
  fa_contigs=$(grep '^>' "$FASTA_TEMP" | head -n5 | sed 's/^>//; s/ .*//')
  if [ -z "$fa_contigs" ]; then
    log_message "ERROR: No contig names found in the processed FASTA. Aborting."
    exit 1
  fi
  log_message "FASTA contigs (after transformation):"
  echo "$fa_contigs" | tee -a "$LOGFILE"
  
  # Save the final processed FASTA.
  cp "$FASTA_TEMP" "$FASTA_PROCESSED"
  log_message "Reference processing complete. Final FASTA: $FASTA_PROCESSED"
else
  log_message "Reference processing already complete; using $FASTA_PROCESSED"
fi
log_message "-------------------------------------------"

########################################################################
# Step 3: Process the BED file if provided.
########################################################################
if [ -n "$BED" ]; then
  log_message "Step 3: Processing BED file..."
  bed_contigs=$(grep -v '^#' "$BED" | head -n5 | awk '{print $1}')
  log_message "BED contigs (original):"
  echo "$bed_contigs" | tee -a "$LOGFILE"
  
  bed_first=$(echo "$bed_contigs" | head -n1)
  if [[ $bed_first == chr* ]]; then
    bed_prefix="chr"
  else
    bed_prefix=""
  fi

  if [ "$bed_prefix" != "$expected_prefix" ]; then
    log_message "Inconsistent BED naming detected. Expected prefix: [$expected_prefix] ; BED prefix: [$bed_prefix]."
    modified_bed="$intermdir/$(basename "$BED" | sed 's/\.[^.]*$//').modified.bed"
    if [ "$expected_prefix" = "" ]; then
      log_message "Removing 'chr' prefix from BED file..."
      sed 's/^chr//' "$BED" > "$modified_bed"
    else
      log_message "Adding 'chr' prefix to BED file..."
      sed 's/^\([0-9XYM]\)/chr\1/' "$BED" > "$modified_bed"
    fi
    log_message "Modified BED file saved as: $modified_bed"
    log_message "First 5 lines of modified BED file:"
    head -n 5 "$modified_bed" | tee -a "$LOGFILE"
    BED_FINAL="$modified_bed"
  else
    log_message "BED naming is consistent with BAM naming."
    BED_FINAL="$BED"
  fi

  bed_contigs=$(grep -v '^#' "$BED_FINAL" | head -n5 | awk '{print $1}')
  log_message "BED contigs (final):"
  echo "$bed_contigs" | tee -a "$LOGFILE"
  log_message "-------------------------------------------"
fi

########################################################################
# Step 4: Compare contig names across BAM, Reference FASTA, and BED (if provided)
########################################################################
log_message "Step 4: Comparing contig names across files..."
log_message "BAM contigs:"
echo "$bam_contigs" | tee -a "$LOGFILE"
log_message "Reference FASTA contigs:"
echo "$fa_contigs" | tee -a "$LOGFILE"
if [ -n "$BED" ]; then
  log_message "BED contigs:"
  echo "$bed_contigs" | tee -a "$LOGFILE"
fi

# Updated check_prefix: if expected_prefix is empty, require contig be entirely numeric.
check_prefix() {
  local contigs="$1"
  local prefix="$2"
  while IFS= read -r line; do
    if [ "$prefix" = "chr" ]; then
      if [[ $line != chr* ]]; then
        return 1
      fi
    else
      if ! [[ $line =~ ^[0-9]+$ ]]; then
        return 1
      fi
    fi
  done <<< "$contigs"
  return 0
}

if ! check_prefix "$bam_contigs" "$expected_prefix"; then
  log_message "ERROR: BAM contigs do not follow the expected naming prefix."
  exit 1
fi

if ! check_prefix "$fa_contigs" "$expected_prefix"; then
  log_message "ERROR: Reference FASTA contigs do not follow the expected naming prefix after processing."
  exit 1
fi

if [ -n "$BED" ]; then
  if ! check_prefix "$bed_contigs" "$expected_prefix"; then
    log_message "ERROR: BED contigs do not follow the expected naming prefix after processing."
    exit 1
  fi
fi
log_message "Contig naming is consistent across BAM, Reference FASTA, and BED (if provided)."
log_message "-------------------------------------------"

########################################################################
# Step 5: Filter and reorder the final FASTA based on BAM contigs.
########################################################################
log_message "Step 5: Filtering and reordering FASTA using BAM contig names..."
bam_contigs_full=$(samtools view -H "$sorted_bam" | awk '/^@SQ/ { for(i=1;i<=NF;i++) { if($i ~ /^SN:/) { split($i,a,":"); print a[2] } } }')
echo "$bam_contigs_full" > "$intermdir/contigs.order.txt"
log_message "Contigs for filtering (first 5 lines):"
head -n 5 "$intermdir/contigs.order.txt" | tee -a "$LOGFILE"

ordered_fasta="$intermdir/$(basename "$FASTA_PROCESSED" | sed 's/\.[^.]*$//').ordered.fasta"
seqkit grep -n --threads "$THREADS" -f "$intermdir/contigs.order.txt" "$FASTA_PROCESSED" > "$ordered_fasta"
if [ ! -s "$ordered_fasta" ]; then
  log_message "ERROR: Ordered FASTA file ($ordered_fasta) is empty. Please check your contig names."
  exit 1
fi
log_message "Ordered FASTA created: $ordered_fasta"
log_message "Indexing ordered FASTA..."
samtools faidx "$ordered_fasta"
log_message "-------------------------------------------"

########################################################################
# Step 6: Run the variant calling pipeline.
########################################################################
log_message "Step 6: Running variant calling pipeline with $THREADS threads..."
if [ -n "$BED" ]; then
  bcftools mpileup -q 30 -Q 20 --threads "$THREADS" --targets-file "$BED_FINAL" -f "$ordered_fasta" "$sorted_bam" |
    bcftools call -c -v --ploidy 2 --threads "$THREADS" -o "$OUTVCF"
else
  bcftools mpileup -q 30 -Q 20 --threads "$THREADS" -f "$ordered_fasta" "$sorted_bam" |
    bcftools call -c -v --ploidy 2 --threads "$THREADS" -o "$OUTVCF"
fi
log_message "Variant calling complete! Output saved to: $OUTVCF"

########################################################################
# Step 7: Cleanup intermediate files.
########################################################################
log_message "Cleaning up intermediate files..."
rm -rf "$intermdir"
log_message "Intermediate files have been removed."
