#!/bin/bash
# check_input.sh
#
# Check input files and append one line per file to a report TSV.
# Always exits 0; main.nf reads the report and stops the run with one clear
# message listing every bad file.
#
#   report columns: file_name  status(OK|FAIL)  reads(- if not counted)  message
#
# Levels
#   quick  milliseconds per file
#          FASTQ: gzip header present (for .gz), not empty, first read parses
#          BAM  : samtools quickcheck -u (header + EOF marker), at least one read
#   full   quick, then read the whole file
#          FASTQ: seqkit stats (catches truncation and broken records anywhere)
#          BAM  : samtools view -c
#   explain  turn a seqkit error log (from the split step) into a FAIL line
#
# Usage
#   check_input.sh --mode fastq|bam --level quick|full --out report.tsv [--threads N] FILE...
#   check_input.sh --mode fastq --level explain --errlog split.err --out report.tsv FILE

set -uo pipefail

mode=""; level=""; out=""; threads=2; errlog=""
files=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode)    mode="$2";    shift 2 ;;
        --level)   level="$2";   shift 2 ;;
        --out)     out="$2";     shift 2 ;;
        --threads) threads="$2"; shift 2 ;;
        --errlog)  errlog="$2";  shift 2 ;;
        -*) echo "Unknown option: $1" >&2; exit 1 ;;
        *)  files+=( "$1" ); shift ;;
    esac
done
[[ "$mode" == "fastq" || "$mode" == "bam" ]] || { echo "--mode must be fastq or bam" >&2; exit 1; }
[[ "$level" == "quick" || "$level" == "full" || "$level" == "explain" ]] || { echo "--level must be quick, full or explain" >&2; exit 1; }
[[ -n "$out" && ${#files[@]} -gt 0 ]] || { echo "Usage: $0 --mode M --level L --out report.tsv FILE..." >&2; exit 1; }

# A missing tool is a problem with the container, not with the input files:
# fail the task loudly instead of reporting the files as bad.
need() { command -v "$1" > /dev/null 2>&1 || { echo "check_input.sh: required tool '$1' not found in this container" >&2; exit 127; }; }
need gzip; need awk; need od
if [[ "$mode" == "bam" ]]; then need samtools; fi
if [[ "$mode" == "fastq" && "$level" == "full" ]]; then need seqkit; fi

REDO="Re-copy or re-download the file."

row() {  # file status reads message
    printf '%s\t%s\t%s\t%s\n' "$(basename "$1")" "$2" "$3" "$4" >> "$out"
    echo "[check_input] $(basename "$1"): $2 $3 $4" >&2
}

# Map seqkit's error text to a plain-language reason
fastq_reason() {  # errlog name
    local reason
    reason=$(tr -d '\033' < "$1" | sed 's/\[[0-9;]*m//g' | grep 'ERRO' | head -n 1 \
             | sed -e 's/^.*ERRO\] *//' -e "s|^.*$2: *||")
    case "$reason" in
        *EOF*)          echo "the file ends unexpectedly - it was probably not copied or downloaded completely" ;;
        *gzip*|*flate*) echo "it is not a valid gzip file" ;;
        *unequal*)      echo "some reads have sequence and quality lines of different lengths" ;;
        *)              echo "it is not a valid FASTQ file${reason:+ ($reason)}" ;;
    esac
}

# --- quick checks: print a FAIL message, or nothing if OK ---------------------
quick_fastq() {
    local f="$1" first
    if [[ "$f" == *.gz ]]; then
        if [[ "$(head -c 2 "$f" | od -An -tx1 | tr -d ' \n')" != "1f8b" ]]; then
            echo "cannot be read: it is not a valid gzip file. $REDO"; return
        fi
        first=$(gzip -cd "$f" 2>/dev/null | head -n 4)
    else
        first=$(head -n 4 "$f")
    fi
    if [[ -z "$first" ]]; then
        echo "contains no reads."; return
    fi
    # first record: @header / sequence / + / quality of the same length
    if ! printf '%s\n' "$first" | awk 'NR==1 && !/^@/ {bad=1}
                                       NR==2 {s=length($0)}
                                       NR==3 && !/^\+/ {bad=1}
                                       NR==4 {q=length($0)}
                                       END {exit (bad || NR<4 || s==0 || s!=q)}'; then
        echo "cannot be read: it is not a valid FASTQ file (the first read is malformed). $REDO"; return
    fi
}

quick_bam() {
    local f="$1" first
    # -u: ONT BAMs are unaligned, so do not require reference targets in the header
    if ! samtools quickcheck -u "$f" 2>/dev/null; then
        echo "cannot be read: it is not a valid BAM file, or it is incomplete (missing end-of-file marker). $REDO"; return
    fi
    first=$(samtools view "$f" 2>/dev/null | head -n 1)
    if [[ -z "$first" ]]; then
        echo "contains no reads."; return
    fi
}

# --- main ---------------------------------------------------------------------
for f in "${files[@]}"; do
    name=$(basename "$f")

    if [[ "$level" == "explain" ]]; then
        row "$f" FAIL - "cannot be read: $(fastq_reason "$errlog" "$name"). $REDO"
        continue
    fi

    if [[ "$mode" == "fastq" ]]; then msg=$(quick_fastq "$f"); else msg=$(quick_bam "$f"); fi
    if [[ -n "$msg" ]]; then row "$f" FAIL - "$msg"; continue; fi

    if [[ "$level" == "quick" ]]; then row "$f" OK - ""; continue; fi

    # full: read the whole file
    if [[ "$mode" == "fastq" ]]; then
        # seqkit stats can exit 0 even on a read error, so also check its stderr
        stats=$(seqkit stats -T -j "$threads" "$f" 2> check.err); rc=$?
        if [[ $rc -ne 0 ]] || grep -q 'ERRO' check.err; then
            row "$f" FAIL - "cannot be read: $(fastq_reason check.err "$name"). $REDO"; continue
        fi
        reads=$(printf '%s\n' "$stats" | awk -F'\t' 'NR==2 {print $4}')
    else
        if ! reads=$(samtools view -c -@ "$threads" "$f" 2> check.err) \
           || grep -qiE 'truncated|error' check.err; then
            row "$f" FAIL - "cannot be read: the BAM file is corrupt or incomplete. $REDO"; continue
        fi
    fi
    if [[ "${reads:-0}" -eq 0 ]]; then row "$f" FAIL 0 "contains no reads."; continue; fi
    row "$f" OK "$reads" ""
done

touch "$out"
exit 0
