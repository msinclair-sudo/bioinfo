#!/bin/bash
set -e  # abort on error

# Function to display usage information.
usage() {
  echo "Usage: $0 -b <BAM file> -r <reference FASTA file> -o <output VCF file> [-l <bed file>] [-t <threads>]"
  echo "  -b  Path to the unsorted BAM file."
  echo "  -r  Path to the exhaustive reference FASTA file."
  echo "  -o  Path to the output VCF file."
  echo "  -l  (Optional) BED file with regions to call variants."
  echo "  -t  (Optional) Number of threads to use (default: 4)."
  exit 1
}

# Default thread count.
THREADS=4

# Parse command-line options.
while getopts ":b:r:o:l:t:" opt; do
  case $opt in
    b)
      BAM="$OPTARG"
      ;;
    r)
      FASTA="$OPTARG"
      ;;
    o)
      OUTVCF="$OPTARG"
      ;;
    l)
      BED="$OPTARG"
      ;;
    t)
      THREADS="$OPTARG"
      ;;
    *)
      usage
      ;;
  esac
done

# Check required parameters.
if [ -z "$BAM" ] || [ -z "$FASTA" ] || [ -z "$OUTVCF" ]; then
    usage
fi

echo "Parameters:"
echo "  Unsorted BAM: $BAM"
echo "  Reference FASTA: $FASTA"
echo "  Output VCF: $OUTVCF"
if [ -n "$BED" ]; then
  echo "  BED file: $BED"
fi
echo "  Threads: $THREADS"
echo

#####################################################################
# Step 1: Process BAM - sort, index, and extract first five contigs.
#####################################################################
echo "Sorting BAM file with $THREADS threads..."
sorted_bam="$(dirname "$BAM")/$(basename "$BAM" .bam).sorted.bam"
samtools sort -@ "$THREADS" -o "$sorted_bam" "$BAM"
echo "Sorted BAM: $sorted_bam"

echo "Indexing sorted BAM file..."
samtools index "$sorted_bam"
echo

# Extract first five contig names from BAM header
echo "Extracting first five contig names from BAM header..."
bam_contigs=$(samtools view -H "$sorted_bam" | grep '^@SQ' | head -n5 | sed 's/.*SN:\([^ \t]*\).*/\1/')
echo "BAM contigs:"
echo "$bam_contigs"
echo

# Determine expected naming prefix from BAM (e.g. if first contig starts with "chr")
bam_first=$(echo "$bam_contigs" | head -n1)
if [[ $bam_first == chr* ]]; then
  expected_prefix="chr"
else
  expected_prefix=""
fi
echo "Expected naming prefix (from BAM): [$expected_prefix]"
echo

########################################################################
# Step 2: Process FASTA - reformat, check header naming, transform if needed,
#           and extract first five contigs from FASTA.
########################################################################
echo "Reformatting FASTA file to fixed line width..."
formatted_fasta="$(dirname "$FASTA")/$(basename "$FASTA" | sed 's/\.[^.]*$//').reformatted.fasta"
seqtk seq -l60 "$FASTA" > "$formatted_fasta"
echo "Formatted FASTA: $formatted_fasta"
echo

# Extract first five contigs from the reformatted FASTA
fa_contigs=$(grep '^>' "$formatted_fasta" | head -n5 | sed 's/^>//; s/ .*//')
echo "FASTA contigs (before transformation):"
echo "$fa_contigs"
echo

# Check if FASTA naming is consistent with BAM naming
fast_first=$(echo "$fa_contigs" | head -n1)
if [[ $fast_first == chr* ]]; then
  fasta_prefix="chr"
else
  fasta_prefix=""
fi

if [ "$fasta_prefix" != "$expected_prefix" ]; then
  echo "Inconsistent FASTA naming detected."
  echo "Expected prefix: [$expected_prefix] ; FASTA prefix: [$fasta_prefix]"
  transformed_fasta="$(dirname "$formatted_fasta")/$(basename "$formatted_fasta" | sed 's/\.[^.]*$//').transformed.fasta"
  if [ "$expected_prefix" = "" ]; then
    echo "Removing 'chr' prefix from FASTA headers..."
    sed 's/^>chr/>/' "$formatted_fasta" > "$transformed_fasta"
  else
    echo "Adding 'chr' prefix to FASTA headers..."
    sed 's/^>/>chr/' "$formatted_fasta" > "$transformed_fasta"
  fi
  echo "Transformed FASTA file created: $transformed_fasta"
  # Re-extract FASTA contigs from transformed FASTA.
  fa_contigs=$(grep '^>' "$transformed_fasta" | head -n5 | sed 's/^>//; s/ .*//')
  echo "FASTA contigs (after transformation):"
  echo "$fa_contigs"
  FASTA_final="$transformed_fasta"
else
  echo "FASTA naming is consistent with BAM header."
  FASTA_final="$formatted_fasta"
fi
echo

##########################################################################
# Step 3: Process BED (Optional) - check naming and transform if needed.
##########################################################################
if [ -n "$BED" ]; then
  echo "Processing BED file..."
  bed_contigs=$(grep -v '^#' "$BED" | head -n5 | awk '{print $1}')
  echo "BED contigs (original):"
  echo "$bed_contigs"
  
  # Determine BED file prefix based on first contig.
  bed_first=$(echo "$bed_contigs" | head -n1)
  if [[ $bed_first == chr* ]]; then
    bed_prefix="chr"
  else
    bed_prefix=""
  fi

  if [ "$bed_prefix" != "$expected_prefix" ]; then
    echo "Inconsistent BED naming detected."
    echo "Expected prefix: [$expected_prefix] ; BED prefix: [$bed_prefix]"
    modified_bed="$(dirname "$BED")/$(basename "$BED" | sed 's/\.[^.]*$//').modified.bed"
    if [ "$expected_prefix" = "" ]; then
      echo "Removing 'chr' prefix from BED file..."
      sed 's/^chr//' "$BED" > "$modified_bed"
    else
      echo "Adding 'chr' prefix to BED file..."
      sed 's/^\([0-9XYM]\)/chr\1/' "$BED" > "$modified_bed"
    fi
    echo "Modified BED file created: $modified_bed"
    echo "First 5 lines of modified BED file:"
    head -n 5 "$modified_bed"
    BED_final="$modified_bed"
  else
    echo "BED naming is consistent with BAM header."
    BED_final="$BED"
  fi
  echo

  # Re-extract BED contigs for later comparison
  bed_contigs=$(grep -v '^#' "$BED_final" | head -n5 | awk '{print $1}')
  echo "BED contigs (final):"
  echo "$bed_contigs"
  echo
fi

##########################################################################
# Step 4: Compare first five contigs from BAM, FASTA, and BED (if provided)
##########################################################################
echo "Comparing contig names from BAM, FASTA, and BED (if provided)..."
echo "BAM contigs:"
echo "$bam_contigs"
echo "FASTA contigs:"
echo "$fa_contigs"
if [ -n "$BED" ]; then
  echo "BED contigs:"
  echo "$bed_contigs"
fi

# Simple check: ensure that each file's contigs start with the expected prefix.
# If expected_prefix is "chr", then each contig should start with "chr".
# If expected_prefix is empty, then none should start with "chr".

check_prefix() {
  local contigs="$1"
  local prefix="$2"
  while IFS= read -r line; do
    if [ "$prefix" = "chr" ]; then
      if [[ $line != chr* ]]; then
        return 1
      fi
    else
      if [[ $line == chr* ]]; then
        return 1
      fi
    fi
  done <<< "$contigs"
  return 0
}

if ! check_prefix "$bam_contigs" "$expected_prefix"; then
  echo "Error: BAM contigs do not follow the expected naming prefix."
  exit 1
fi

if ! check_prefix "$fa_contigs" "$expected_prefix"; then
  echo "Error: FASTA contigs do not follow the expected naming prefix after transformation."
  exit 1
fi

if [ -n "$BED" ]; then
  if ! check_prefix "$bed_contigs" "$expected_prefix"; then
    echo "Error: BED contigs do not follow the expected naming prefix after transformation."
    exit 1
  fi
fi

echo "Contig naming is consistent across BAM, FASTA, and BED (if provided)."
echo

##########################################################################
# Step 5: Filter and reorder the FASTA to include only BAM contigs.
##########################################################################
echo "Filtering and reordering FASTA using BAM contig names..."
# Re-extract full list of contigs from sorted BAM header for filtering.
samtools view -H "$sorted_bam" | awk '/^@SQ/ { for(i=1;i<=NF;i++) { if($i ~ /^SN:/) { split($i,a,":"); print a[2] } } }' > contigs.order.txt
echo "BAM contigs used for filtering (first 5 shown):"
head -n 5 contigs.order.txt
echo

ordered_fasta="$(dirname "$FASTA_final")/$(basename "$FASTA_final" | sed 's/\.[^.]*$//').ordered.fasta"
seqkit grep -n --threads "$THREADS" -f contigs.order.txt "$FASTA_final" > "$ordered_fasta"
if [ ! -s "$ordered_fasta" ]; then
  echo "Error: Ordered FASTA file ($ordered_fasta) is empty. Please check your contig names."
  exit 1
fi
echo "Ordered FASTA created: $ordered_fasta"
echo "Indexing ordered FASTA..."
samtools faidx "$ordered_fasta"
echo

##########################################################################
# Step 6: Run the variant calling pipeline.
##########################################################################
echo "Running variant calling pipeline with $THREADS threads..."
if [ -n "$BED" ]; then
  bcftools mpileup -q 30 -Q 20 --threads "$THREADS" -l "$BED_final" -f "$ordered_fasta" "$sorted_bam" | \
    bcftools call -c -v --ploidy 2 --threads "$THREADS" -o "$OUTVCF"
else
  bcftools mpileup -q 30 -Q 20 --threads "$THREADS" -f "$ordered_fasta" "$sorted_bam" | \
    bcftools call -c -v --ploidy 2 --threads "$THREADS" -o "$OUTVCF"
fi

echo "Variant calling complete! Output saved to: $OUTVCF"
