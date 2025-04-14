#!/bin/bash
set -euo pipefail

################################################################################
# GLOBAL SETTINGS & LOGGING
################################################################################
# (Nodes: A, E, n38, and Check.Provided.Inputs)
#
# - Create intermediates and output directories.
# - Initialize two logging systems:
#   1. Typical log (LOGFILE): contains all steps and messages.
#   2. JSON log (JSON_LOG): records key variables (e.g. expected naming format, input files,
#      and checkpoints) to support a resume mechanism.
################################################################################

intermdir="./intermediates"
mkdir -p "$intermdir"
LOGFILE="$intermdir/pipeline.log"
JSON_LOG="$intermdir/pipeline_status.json"

# Initialize JSON log if not already present.
if [ ! -f "$JSON_LOG" ]; then
  echo "{}" > "$JSON_LOG"
fi

# Typical logging: prints with timestamp and appends to LOGFILE.
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
# Records input variable values to JSON and verifies required files exist.
check_provided_inputs() {
  if [ -z "${BAM:-}" ] || [ -z "${FASTA:-}" ] || [ -z "${OUTVCF:-}" ]; then
    echo "Usage: $0 -b <BAM file> -r <reference FASTA file> -o <output VCF file> [-l <BED file>] [-t <threads>]"
    exit 1
  fi

  # Ensure provided input files exist and are non-empty.
  test_file_valid "$BAM"
  test_file_valid "$FASTA"
  if [ -n "${BED:-}" ]; then
    test_file_valid "$BED"
  fi

  update_json_log "BAM" "$BAM"
  update_json_log "FASTA" "$FASTA"
  update_json_log "OUTVCF" "$OUTVCF"
  if [ -n "${BED:-}" ]; then
    update_json_log "BED" "$BED"
  fi
  update_json_log "THREADS" "$THREADS"
}

# Resume state (HP2)
# Reads JSON log to determine if a checkpoint exists for resume.
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
# HELPER FUNCTIONS (see HELPERS subgraph)
################################################################################

# HP1: check_prefix
# Verifies that each contig in the provided string (one per line) matches the expected format.
# Expected format is either "chr" (e.g. "chr1") or numeric (e.g., "1").
check_prefix() {
  local contigs="$1"
  local expected="$2"  # "chr" or ""
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

# n27: Check.Naming.Format – examine a sample of five contigs for correct naming.
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

# n29: Test.File.Valid – checks that the provided file exists and is non-empty.
test_file_valid() {
  local file="$1"
  if [ ! -s "$file" ]; then
    log_message "ERROR: File '$file' does not exist or is empty. Aborting."
    exit 1
  fi
}

# Subgraph s1: Edit.Names helper function
#
# This function converts a given vector of contig/chromosome names to the expected format.
# It handles three types of input formats:
#   n6: "RefSeq format" (e.g., NC_000001)
#   n5: "'chrN' format" (e.g., chr1)
#   n4: "'n' format" (e.g., 1)
#
# It returns (echoes) a vector of names in the expected format.
#
# n45: The function expects a vector of current names supplied via STDIN.
# n44: The conversion uses an appropriate method (here using sed processing) depending on the format.
edit_names() {
  # Usage: edit_names "current_names" "expected_format"
  # current_names: multi-line string of contig names.
  # expected_format: either "chr" (for chrN) or "" (for numeric-only).
  local names="$1"
  local expected="$2"
  while IFS= read -r name; do
    # If the name is in RefSeq (NC_) format.
    if [[ $name == NC_* ]]; then
      # Remove "NC_", strip any leading zeroes, extract numeric part.
      num=$(echo "$name" | sed -E 's/^NC_0*([0-9]+).*/\1/')
      if [ "$expected" == "chr" ]; then
        echo "chr${num}"
      else
        echo "${num}"
      fi
    elif [[ $name == chr* ]]; then
      # If the name already has "chr" but expected is numeric, remove it.
      if [ "$expected" == "" ]; then
        echo "$name" | sed 's/^chr//'
      else
        echo "$name"
      fi
    else
      # Assume it is numeric.
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
  echo "Usage: $0 -b <BAM file> -r <reference FASTA file> -o <output VCF file> [-l <BED file>] [-t <threads>]"
  exit 1
}

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
#
# G: Process BAM file.
# I: Sort BAM using samtools sort.
# J: Index sorted BAM.
# K: Extract all contigs (display first 5 to stdout, log all).
# L: Check.Naming.Format for BAM contigs.
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

# Extract all contigs from BAM header.
bam_contigs=$(samtools view -H "$sorted_bam" | grep '^@SQ' | sed -E 's/.*SN:([^[:space:]]+).*/\1/')
log_message "BAM contigs (logged complete list in $LOGFILE):"
echo "$bam_contigs" >> "$LOGFILE"
echo "BAM contigs (first 5):"
echo "$bam_contigs" | head -n5

# Determine expected naming format from first contig.
bam_first=$(echo "$bam_contigs" | head -n1)
if [[ $bam_first == chr* ]]; then
  expected_prefix="chr"
else
  expected_prefix=""
fi
update_json_log "expected_prefix" "$expected_prefix"
log_message "Expected contig naming format: $( [ "$expected_prefix" == "chr" ] && echo \"'chrN'\" || echo \"'n'\" )."

# Check if the first five BAM contigs are in the expected format.
if ! check_naming_format "$bam_contigs" "$expected_prefix"; then
  log_message "ERROR: BAM contig names do not match expected format. Aborting."
  exit 1
fi
update_json_log "last_checkpoint" "G: BAM file processed and validated"

################################################################################
# STEP M: Process Reference FASTA (Step 2)
################################################################################
#
# M: Process Reference FASTA.
# P: Extract first 5 contigs from reformatted FASTA.
# n26: Call Check.Naming.Format.
# If not valid, then we call Edit.Names (subgraph s1) to convert.
################################################################################

FASTA_BASENAME=$(basename "$FASTA" | sed 's/\.[^.]*$//')
# Reformat FASTA: Ensure it has a consistent line length.
reformatted_fasta="$intermdir/${FASTA_BASENAME}.reformatted.fasta"
seqtk seq -l60 "$FASTA" > "$reformatted_fasta"
test_file_valid "$reformatted_fasta"
log_message "Reformatted FASTA saved."

# Extract contig names from FASTA header.
fa_contigs=$(grep '^>' "$reformatted_fasta" | sed 's/^>//; s/ .*//')
log_message "Reference FASTA contigs (logged in $LOGFILE), first 5 displayed below:"
echo "$fa_contigs" >> "$LOGFILE"
echo "Reference FASTA contigs (first 5):"
echo "$fa_contigs" | head -n5

# n26: Check naming format for FASTA contigs.
if ! check_naming_format "$fa_contigs" "$expected_prefix"; then
  log_message "FASTAs contig names not in expected format. Converting using Edit.Names helper."
  # Call edit_names: supply current vector and expected prefix.
  fa_contigs_fixed=$(edit_names "$fa_contigs" "$expected_prefix")
  # Save a new, fixed FASTA.
  FASTA_PROCESSED="$intermdir/${FASTA_BASENAME}.processed.fasta"
  # Replace header lines in the FASTA with the fixed names.
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
  ' "$reformatted_fasta" > "$FASTA_PROCESSED"
  fa_contigs=$(grep '^>' "$FASTA_PROCESSED" | sed 's/^>//; s/ .*//')
else
  FASTA_PROCESSED="$intermdir/${FASTA_BASENAME}.processed.fasta"
  cp "$reformatted_fasta" "$FASTA_PROCESSED"
fi
update_json_log "last_checkpoint" "M: Reference FASTA processed and validated"

################################################################################
# STEP V: Process BED File (Step 3)
################################################################################
#
# V: Decision: BED file provided?
# X: Extract first 5 contigs from BED file.
# Y: Call Check.Naming.Format on BED contigs.
# Z: Validate BED format; if inconsistent, call Edit.Names to convert.
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
    # edit_names expects a list of names; here we process each BED contig.
    fixed_bed_contigs=$(edit_names "$bed_contigs" "$expected_prefix")
    # Edit the BED file: replace the first column with corrected names.
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
#
# AD: Compare contig names across files.
# AE: Display the first 5 contigs for BAM, Reference FASTA, and BED (if provided).
################################################################################

log_message "Comparing contig naming across files:"
echo "BAM contigs (first 5):" && echo "$bam_contigs" | head -n5
echo "FASTA contigs (first 5):" && echo "$fa_contigs" | head -n5
if [ -n "${BED:-}" ]; then
  bed_contigs=$(grep -v '^#' "$BED_FINAL" | awk '{print $1}')
  echo "BED contigs (first 5):" && echo "$bed_contigs" | head -n5
fi

# Double-check each file’s contigs conform to expected format.
if ! check_prefix "$bam_contigs" "$expected_prefix" || \
   ! check_prefix "$fa_contigs" "$expected_prefix" || \
   { [ -n "${BED:-}" ] && ! check_prefix "$bed_contigs" "$expected_prefix"; }; then
  log_message "ERROR: Inconsistent contig naming detected across files. Aborting."
  exit 1
fi
update_json_log "last_checkpoint" "AD/AE: Contig comparison successful"

################################################################################
# STEP AH & AJ: Format Reference FASTA (Step 5)
################################################################################
#
# AH: Format reference FASTA.
# AJ: Extract full BAM contig list (ordered) from sorted BAM header.
# AK: Use seqkit grep with BAM contig list to create ordered FASTA.
# n36: Test.File.Valid for ordered FASTA.
# Then index ordered FASTA with samtools faidx.
################################################################################

log_message "Formatting reference FASTA to match BAM contig order..."
bam_contigs_full=$(samtools view -H "$sorted_bam" | awk '/^@SQ/ { for(i=1;i<=NF;i++) { if($i ~ /^SN:/) { split($i,a,":"); print a[2] } } }')
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
#
# AO: Variant calling step.
# AP: Check if BED provided; then run bcftools mpileup using targets if yes.
# n20: Log checkpoint.
# AS: VCF file produced.
################################################################################

log_message "Starting variant calling with $THREADS threads..."
update_json_log "last_checkpoint" "AO: Starting variant calling"
if [ -n "${BED:-}" ]; then
  bcftools mpileup -q 30 -Q 20 --threads "$THREADS" --targets-file "$BED_FINAL" -f "$ordered_fasta" "$sorted_bam" | \
    bcftools call -c -v --ploidy 2 --threads "$THREADS" -o "$OUTVCF"
else
  bcftools mpileup -q 30 -Q 20 --threads "$THREADS" -f "$ordered_fasta" "$sorted_bam" | \
    bcftools call -c -v --ploidy 2 --threads "$THREADS" -o "$OUTVCF"
fi
test_file_valid "$OUTVCF"
update_json_log "last_checkpoint" "AO/AP/AS: Variant calling complete; VCF produced"

################################################################################
# STEP AT & AU: Cleanup (Step 7)
################################################################################
#
# AT: Cleanup intermediate files.
# AU: Remove intermediates directory.
# n21: Log Successful checkpoint to json.
# AV: Exit.
################################################################################

log_message "Cleaning up intermediate files..."
rm -rf "$intermdir"
log_message "Cleanup complete. Pipeline finished successfully."
update_json_log "last_checkpoint" "AT/AU: Pipeline completed successfully"
exit 0
