#!/bin/bash
# combine_demux.sh
#
# Merge per-input dorado demux output directories into one file per barcode.
#
# Works with both dorado output layouts:
#   - dorado <= 1.1.x : flat files      <dir>/<prefix>_sw_barcode001.fastq
#   - dorado >= 1.3.1 : nested layout   <dir>/.../fastq_pass/<barcode>/<fc>_pass_<barcode>_<run>_0.fastq.gz
#                       (ONT Output Specifications 26.01)
# Files are found recursively and assigned to a barcode from the file name
# first, then from any parent directory name. Plain and gzipped FASTQ are
# both accepted. Barcode numbers are normalised to 3 digits (barcode001) so
# downstream names are unchanged.
#
# Usage:
#   combine_demux.sh --mode fastq --pattern '*_Demux_fastq' --outdir output_fastq
#   combine_demux.sh --mode bam   --pattern '*_Demux_bam'   --outdir output_bam

set -euo pipefail
shopt -s nullglob

usage() {
    echo "Usage: $0 --mode fastq|bam --pattern '<dir glob>' --outdir <dir>" >&2
    exit 1
}

mode=""
pattern=""
outdir=""
while [[ $# -gt 0 ]]; do
    case "$1" in
        --mode)    mode="$2";    shift 2 ;;
        --pattern) pattern="$2"; shift 2 ;;
        --outdir)  outdir="$2";  shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown argument: $1" >&2; usage ;;
    esac
done
[[ "$mode" == "fastq" || "$mode" == "bam" ]] || usage
[[ -n "$pattern" && -n "$outdir" ]] || usage

# shellcheck disable=SC2206
dirs=( $pattern )
(( ${#dirs[@]} )) || { echo "No directories matching '${pattern}' found" >&2; exit 1; }
# Trailing '/' makes find descend into directories that Nextflow staged as
# symlinks (local runs) without needing GNU-only 'find -L'
roots=()
for d in "${dirs[@]}"; do roots+=( "${d%/}/" ); done

mkdir -p "$outdir"
listdir=$(mktemp -d combine_lists.XXXXXX)

# Log the raw dorado layout so a future structure change is easy to diagnose
echo "=== dorado demux output layout (first 200 files) ===" >&2
find "${roots[@]}" -type f | head -n 200 >&2 || true
echo "====================================================" >&2

# Assign one file to a barcode key; prints the key or nothing
classify() {
    local path="$1"
    local base comp n
    base=$(basename "$path")

    # 1) file name
    if [[ "$base" =~ [Bb]arcode0*([0-9]+) ]]; then
        n=$((10#${BASH_REMATCH[1]}))
        printf 'barcode%03d\n' "$n"; return
    fi
    if [[ "$base" == *unclassified* ]]; then
        echo "unclassified"; return
    fi

    # 2) parent directories, nearest first
    local parent
    parent=$(dirname "$path")
    while [[ "$parent" != "." && "$parent" != "/" && -n "$parent" ]]; do
        comp=$(basename "$parent")
        if [[ "$comp" =~ ^([A-Za-z0-9-]+_)?[Bb]arcode0*([0-9]+)$ ]]; then
            n=$((10#${BASH_REMATCH[2]}))
            printf 'barcode%03d\n' "$n"; return
        fi
        if [[ "$comp" == "unclassified" || "$comp" == *_unclassified ]]; then
            echo "unclassified"; return
        fi
        parent=$(dirname "$parent")
    done
}

if [[ "$mode" == "fastq" ]]; then
    find_expr=( \( -name '*.fastq' -o -name '*.fastq.gz' -o -name '*.fq' -o -name '*.fq.gz' \) )
else
    find_expr=( -name '*.bam' )
fi

n_files=0
n_unmatched=0
while IFS= read -r f; do
    n_files=$((n_files + 1))
    key=$(classify "$f")
    if [[ -n "$key" ]]; then
        printf '%s\n' "$f" >> "${listdir}/${key}.list"
    else
        n_unmatched=$((n_unmatched + 1))
        echo "WARNING: could not assign a barcode to: $f" >&2
    fi
done < <(find "${roots[@]}" -type f "${find_expr[@]}")

(( n_files )) || { echo "No ${mode} files found under: ${dirs[*]}" >&2; exit 1; }

# Fail loudly rather than silently dropping reads from the demux report
if (( n_unmatched )); then
    echo "ERROR: ${n_unmatched} of ${n_files} ${mode} file(s) could not be assigned to a barcode or 'unclassified'." >&2
    echo "       The dorado output naming has probably changed - see the layout above." >&2
    exit 1
fi
lists=( "${listdir}"/*.list )

for list in "${lists[@]}"; do
    key=$(basename "$list" .list)
    mapfile -t files < "$list"
    if [[ "$mode" == "fastq" ]]; then
        {
            for f in "${files[@]}"; do
                case "$f" in
                    *.gz) zcat "$f" ;;
                    *)    cat  "$f" ;;
                esac
            done
        } | gzip > "${outdir}/${key}.fastq.gz"
    else
        samtools merge -f "${outdir}/${key}.bam" "${files[@]}"
        samtools index "${outdir}/${key}.bam"
    fi
    echo "${key}: merged ${#files[@]} file(s)" >&2
done

rm -rf "$listdir"
echo "Combined ${n_files} ${mode} file(s) into ${#lists[@]} output(s); ${n_unmatched} unassigned" >&2
