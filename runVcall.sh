#!/bin/bash
set -euo pipefail

################################################################################
# GLOBAL SETTINGS & LOGGING
################################################################################
# Create intermediates and output directories.
# Initialize two logging systems:
#   1. Typical log (LOGFILE): contains step‐by‐step messages.
#   2. JSON log (JSON_LOG): records key variables (expected naming format, inputs,
#      and checkpoints) for resume functionality.
################################################################################

intermdir="./intermediates"
mkdir -p "$intermdir"
LOGFILE="$intermdir/pipeline.log"
JSON_LOG="$intermdir/pipeline_status.json"

# Initialize JSON log if not already present.
if [ ! -f "$JSON_LOG" ]; then
  echo "{}" > "$JSON_LOG"
fi

# Typical logging function.
log_message() {
  echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" | tee -a "$LOGFILE"
}

# Update JSON log with a key/value pair (requires jq).
update_json_log() {
  local key="$1"
  local value="$2"
  local tmp
  tmp=$(mktemp)
  jq --arg key "$key" --arg value "$value" '.[$key] = $value' "$JSON_LOG" > "$tmp" && mv "$tmp" "$JSON_LOG"
}

# Check Provided Inputs (n38)
# Verify required parameters exist and record them.
check_provided_inputs() {
  if [ -z "${BAM:-}" ] || [ -z "${FASTA:-}" ] || [ -z "${OUTVCF:-}" ]; then
    echo "Usage: $0 -b <BAM file> -r <reference FASTA file> -o <output VCF file> [-l <BED file>] [-t <threads>] [-m <mapping file>]"
    exit 1
  fi

  test_file_valid "$BAM"
  test_file_valid "$FASTA"
  if [ -n "${BED:-}" ]; then
    test_file_valid "$BED"
  fi
  if [ -n "${MAPPING:-}" ]; then
    test_file_valid "$MAPPING"
  fi

  update_json_log "BAM" "$BAM"
  update_json_log "FASTA" "$FASTA"
  update_json_log "OUTVCF" "$OUTVCF"
  if [ -n "${BED:-}" ]; then
    update_json_log "BED" "$BED"
  fi
  if [ -n "${MAPPING:-}" ]; then
    update_json_log "MAPPING" "$MAPPING"
  fi
  update_json_log "THREADS" "$THREADS"
}

# Resume state (HP2): Reads JSON log for the last checkpoint.
resume_state() {
  if [ -f "$JSON_LOG" ]; then
    local last_checkpoint
    last_checkpoint=$(jq -r '.last_checkpoint // empty' "$JSON_LOG")
    if [ -n "$last_checkpoint" ]; then
      log_message "Resuming from checkpoint: $last_checkpoint"
      LAST_CHECKPOINT="$last_checkpoint"
    fi
  fi
}

################################################################################
# HELPER FUNCTIONS (per HELPERS subgraph)
################################################################################

# HP1: check_prefix
# Verifies that each contig (one per line) matches the expected format.
# Expected format is either "chr" (e.g., "chr1") or numeric-only (e.g., "1").
check_prefix() {
  local contigs="$1"
  local expected="$2"  # Expected value: "chr" or ""
  while IFS= read -r line; do
    if [ "$expected" == "chr" ]; then
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

# n27: Check.Naming.Format – Examine a sample (first 5) of contigs for correct naming.
check_naming_format() {
  local names="$1"
  local expected="$2"
  local sample
  sample=$(echo "$names" | head -n 5)
  if check_prefix "$sample" "$expected"; then
    return 0
  else
    return 1
  fi
}

# n29: Test.File.Valid – Checks that the provided file exists and is non-empty.
test_file_valid() {
  local file="$1"
  if [ ! -s "$file" ]; then
    log_message "ERROR: File '$file' does not exist or is empty. Aborting."
    exit 1
  fi
}

# Subgraph s1: Edit.Names helper function.
#
# Converts a given set of contig/chromosome names to the expected format.
#
# It handles three types of input:
#   n6: "RefSeq format" (e.g., NC_000001) – requires a mapping file.
#   n5: "'chrN' format" (e.g., chr1)
#   n4: "'n' format" (e.g., 1)
#
# n45: It accepts a vector of current names.
# n44: It uses the appropriate method based on the file type.
#
# It echoes a vector of names in the expected format.
edit_names() {
  # Usage: edit_names "current_names" "expected_format"
  # current_names: multi-line string of contig names.
  # expected_format: either "chr" (for chrN) or "" (for numeric-only).
  local names="$1"
  local expected="$2"
  while IFS= read -r name; do
    if [[ $name == NC_* ]]; then
      # RefSeq format detected.
      if [ -z "${MAPPING:-}" ]; then
        log_message "ERROR: Detected RefSeq format ($name) but no mapping file provided. Aborting."
        exit 1
      fi
      # Use the mapping file (CSV with header, comma-separated) to map the RefSeq name.
      mapped=$(awk -v ref="$name" 'BEGIN {FS=","} NR>1 {if($2 == ref) print $1}' "$MAPPING")
      if [ -z "$mapped" ]; then
        log_message "ERROR: Unable to map RefSeq name '$name' using mapping file '$MAPPING'. Aborting."
        exit 1
      fi
      echo "$mapped"
    elif [[ $name == chr* ]]; then
      if [ "$expected" == "" ]; then
        echo "$name" | sed 's/^chr//'
      else
        echo "$name"
      fi
    else
      if [ "$expected" == "chr" ]; then
        echo "chr${name}"
      else
        echo "$name"
      fi
    fi
  done <<< "$names"
}

################################################################################
# STEP A: Parse Command-Line Arguments (Start Node)
################################################################################
usage() {
  echo "Usage: $0 -b <BAM file> -r <reference FASTA file> -o <output VCF file> [-l <BED file>] [-t <threads>] [-m <mapping file>]"
  exit 1
}

# Parse arguments. Note the new -m argument for the mapping file.
THREADS=4
while getopts ":b:r:o:l:t:m:" opt; do
  case $opt in
    b) BAM="$OPTARG" ;;
    r) FASTA="$OPTARG" ;;
    o) OUTVCF="$OPTARG" ;;
    l) BED="$OPTARG" ;;
    t) THREADS="$OPTARG" ;;
    m) MAPPING="$OPTARG" ;;
    *) usage ;;
  esac
done

check_provided_inputs
log_message "Parsed command-line arguments and confirmed required files exist."
update_json_log "last_checkpoint" "A: Inputs provided and checked"

################################################################################
# STEP B & C: Check/Create Output and Intermediates Directories
################################################################################
output_dir=$(dirname "$OUTVCF")
if [ ! -d "$output_dir" ]; then
  mkdir -p "$output_dir"
  log_message "Created output directory: $output_dir"
fi
if [ ! -d "$intermdir" ]; then
  mkdir -p "$intermdir"
  log_message "Created intermediates directory: $intermdir"
fi
update_json_log "last_checkpoint" "B: Directories checked/created"

################################################################################
# RESUME STATE CHECK (n1, HP2)
################################################################################
resume_state

################################################################################
# STEP G: Process BAM File (Step 1)
################################################################################
# Sort the BAM file (I), index it (J), extract all contigs (K) and verify naming format (L).
################################################################################

sorted_bam="$intermdir/$(basename "$BAM" .bam).sorted.bam"
if [ ! -s "$sorted_bam" ]; then
  log_message "Sorting BAM file..."
  samtools sort -@ "$THREADS" -o "$sorted_bam" "$BAM"
else
  log_message "Using existing sorted BAM: $sorted_bam"
fi

if [ ! -f "${sorted_bam}.bai" ]; then
  log_message "Indexing BAM file..."
  samtools index "$sorted_bam"
else
  log_message "BAM index already exists."
fi

# Extract all contigs from the BAM header.
bam_contigs=$(samtools view -H "$sorted_bam" | grep '^@SQ' | sed -E 's/.*SN:([^[:space:]]+).*/\1/')
log_message "BAM contigs (complete list logged to $LOGFILE):"
echo "$bam_contigs" >> "$LOGFILE"
echo "BAM contigs (first 5):"
echo "$bam_contigs" | head -n5

# Determine the expected naming format from the first BAM contig.
bam_first=$(echo "$bam_contigs" | head -n1)
if [[ $bam_first == chr* ]]; then
  expected_prefix="chr"
else
  expected_prefix=""
fi
update_json_log "expected_prefix" "$expected_prefix"
log_message "Expected contig naming format: $( [ "$expected_prefix" == "chr" ] && echo \"'chrN'\" || echo \"'n'\" )."

# Verify naming format for BAM contigs.
if ! check_naming_format "$bam_contigs" "$expected_prefix"; then
  log_message "ERROR: BAM contig names do not match expected format. Aborting."
  exit 1
fi
update_json_log "last_checkpoint" "G: BAM file processed and validated"

################################################################################
# STEP M: Process Reference FASTA (Step 2)
################################################################################
# Reformat FASTA, then filter to include only NC contigs (RefSeq).
################################################################################

FASTA_BASENAME=$(basename "$FASTA" | sed 's/\.[^.]*$//')
# Reformat FASTA: ensure consistent line lengths.
reformatted_fasta="$intermdir/${FASTA_BASENAME}.reformatted.fasta"
seqtk seq -l60 "$FASTA" > "$reformatted_fasta"
test_file_valid "$reformatted_fasta"
log_message "Reformatted FASTA saved."

# Filter: Keep only NC contigs from the reformatted FASTA.
filtered_fasta="$intermdir/${FASTA_BASENAME}.filtered.fasta"
seqkit grep -r -p '^NC_' "$reformatted_fasta" > "$filtered_fasta"
test_file_valid "$filtered_fasta"
log_message "Filtered FASTA to include only NC contigs."

# Extract contig names from filtered FASTA.
fa_contigs=$(grep '^>' "$filtered_fasta" | sed 's/^>//; s/ .*//')
# Remove any empty lines from the list.
fa_contigs=$(echo "$fa_contigs" | grep -v '^\s*$')
log_message "Reference FASTA contigs (filtered, complete list logged), first 5 displayed below:"
echo "$fa_contigs" >> "$LOGFILE"
echo "Reference FASTA contigs (first 5):"
echo "$fa_contigs" | head -n5

# Check naming format; if not matching expected, convert using Edit.Names helper.
if ! check_naming_format "$fa_contigs" "$expected_prefix"; then
  log_message "FASTA contig names not in expected format. Converting using Edit.Names helper."
  fa_contigs_fixed=$(edit_names "$fa_contigs" "$expected_prefix")
  # Create a new FASTA with fixed headers.
  FASTA_PROCESSED="$intermdir/${FASTA_BASENAME}.processed.fasta"
  awk -v fixed_names="$fa_contigs_fixed" '
    BEGIN {
      split(fixed_names, arr, "\n");
      idx = 1;
    }
    /^>/ {
      print ">" arr[idx++];
      next;
    }
    { print }
  ' "$filtered_fasta" > "$FASTA_PROCESSED"
  fa_contigs=$(grep '^>' "$FASTA_PROCESSED" | sed 's/^>//; s/ .*//')
else
  FASTA_PROCESSED="$intermdir/${FASTA_BASENAME}.processed.fasta"
  cp "$filtered_fasta" "$FASTA_PROCESSED"
fi
update_json_log "last_checkpoint" "M: Reference FASTA processed and validated"

################################################################################
# STEP V: Process BED File (Step 3)
################################################################################
# If BED is provided, extract its contig names (X), check naming format (Y),
# and if needed, convert using Edit.Names (Z).
################################################################################
if [ -n "${BED:-}" ]; then
  log_message "Processing BED file: $BED"
  bed_contigs=$(grep -v '^#' "$BED" | awk '{print $1}')
  echo "BED contigs (first 5):"
  echo "$bed_contigs" | head -n5
  echo "$bed_contigs" >> "$LOGFILE"
  
  if ! check_naming_format "$bed_contigs" "$expected_prefix"; then
    log_message "BED contig names not in expected format. Converting using Edit.Names helper."
    modified_bed="$intermdir/$(basename "$BED" | sed 's/\.[^.]*$//').modified.bed"
    fixed_bed_contigs=$(edit_names "$bed_contigs" "$expected_prefix")
    awk -v fixed_names="$fixed_bed_contigs" 'BEGIN {
         split(fixed_names, arr, "\n");
         idx = 1;
       }
       {
         $1 = arr[idx++];
         print;
       }' "$BED" > "$modified_bed"
    echo "Modified BED contigs (first 5):"
    head -n5 "$modified_bed" | tee -a "$LOGFILE"
    BED_FINAL="$modified_bed"
  else
    BED_FINAL="$BED"
  fi
  update_json_log "last_checkpoint" "V: BED file processed and validated"
else
  log_message "No BED file provided; skipping BED processing."
fi

################################################################################
# STEP AD & AE: Compare Contig Names Across Files (Step 4)
################################################################################
# Display the first 5 contigs for BAM and Reference FASTA (and BED if provided)
# without enforcing strict consistency (since the BAM file may contain extra contigs).
################################################################################

log_message "Comparing contig naming across files:"
echo "BAM contigs (first 5):"
echo "$bam_contigs" | head -n5
echo "Reference FASTA contigs (first 5):"
echo "$fa_contigs" | head -n5
if [ -n "${BED:-}" ]; then
  bed_contigs=$(grep -v '^#' "$BED_FINAL" | awk '{print $1}')
  echo "BED contigs (first 5):"
  echo "$bed_contigs" | head -n5
fi
update_json_log "last_checkpoint" "AD/AE: Contig naming compared and logged"

################################################################################
# STEP AH & AJ: Format Reference FASTA (Step 5)
################################################################################
# Extract full BAM contig order (AJ) and reorder FASTA accordingly (AK),
# then index the ordered FASTA.
################################################################################

log_message "Formatting reference FASTA to match BAM contig order..."
bam_contigs_full=$(samtools view -H "$sorted_bam" | \
  awk '/^@SQ/ { for(i=1;i<=NF;i++) { if($i ~ /^SN:/) { split($i,a,":"); print a[2] } } }')
echo "$bam_contigs_full" > "$intermdir/contigs.order.txt"
echo "BAM contig order (first 5):"
head -n5 "$intermdir/contigs.order.txt" | tee -a "$LOGFILE"

ordered_fasta="$intermdir/${FASTA_BASENAME}.ordered.fasta"
seqkit grep -n --threads "$THREADS" -f "$intermdir/contigs.order.txt" "$FASTA_PROCESSED" > "$ordered_fasta"
test_file_valid "$ordered_fasta"
samtools faidx "$ordered_fasta"
update_json_log "last_checkpoint" "AH/AJ/AK: FASTA ordered and indexed"

################################################################################
# STEP AO & AP: Variant Calling (Step 6)
################################################################################
# Run variant calling using bcftools mpileup and bcftools call.
# Instead of piping mpileup output to call, write mpileup output to an
# intermediate BCF file, then call variants on that file.
################################################################################

log_message "Starting variant calling with $THREADS threads..."
update_json_log "last_checkpoint" "AO: Starting variant calling"

mpileup_file="$intermdir/mpileup.bcf"
if [ -n "${BED:-}" ]; then
  bcftools mpileup -q 30 -Q 20 --threads "$THREADS" --targets-file "$BED_FINAL" \
    -f "$ordered_fasta" "$sorted_bam" -Ou -o "$mpileup_file"
else
  bcftools mpileup -q 30 -Q 20 --threads "$THREADS" \
    -f "$ordered_fasta" "$sorted_bam" -Ou -o "$mpileup_file"
fi
test_file_valid "$mpileup_file"

bcftools call -c -v --ploidy 2 --threads "$THREADS" -Ou -o "$OUTVCF" "$mpileup_file"
test_file_valid "$OUTVCF"
update_json_log "last_checkpoint" "AO/AP/AS: Variant calling complete; VCF produced"

################################################################################
# STEP AT & AU: Cleanup (Step 7)
################################################################################
# Clean up intermediates and log final successful checkpoint.
################################################################################

log_message "Cleaning up intermediate files..."
rm -rf "$intermdir"
log_message "Cleanup complete. Pipeline finished successfully."
update_json_log "last_checkpoint" "AT/AU: Pipeline completed successfully"
exit 0
