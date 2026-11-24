#!/usr/bin/env bash
# tests/error_tests.sh
#
# Exercise the pipeline's error paths and show the message a customer would see.
#
#   tests/error_tests.sh quick   # input checks only (uses -preview, no Docker, ~1 min)
#   tests/error_tests.sh full    # runtime failures (runs the real pipeline, needs Docker)
#   tests/error_tests.sh all     # both
#
# Run from the repository root. Everything is written under ./error_tests/
# (bad inputs, per-case logs, outdirs and a shared work dir); delete it when done.
#
# Set NXF=/path/to/nextflow to use a specific Nextflow binary.

set -uo pipefail

MODE="${1:-quick}"
NXF="${NXF:-nextflow}"
ROOT="$(pwd)"
T="${ROOT}/error_tests"
FQ_DIR="${ROOT}/test_data/fastq_pass"
BAM_DIR="${ROOT}/test_data/bam_pass"

[[ -f main.nf ]] || { echo "Run this from the pipeline root (where main.nf is)." >&2; exit 1; }
command -v "$NXF" >/dev/null || { echo "Nextflow not found: $NXF" >&2; exit 1; }

mkdir -p "$T"/{logs,work,inputs}
PASS=0; FAIL=0; SHOWN=0
SUMMARY=()

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
# show the part of the console output a customer would read
show_msg() {
    # pipeline's own "... problem:" message if present, else the top of the log
    if grep -q 'problem:' "$1"; then
        sed -n '/problem:/,/^ERROR ~ *$/p' "$1" | grep -v '^ERROR ~ *$' | sed 's/^/    | /'
        return
    fi
    # drop the banner/blank lines, keep the message block
    grep -v -e 'N E X T F L O W' -e '^Launching' -e 'nextflow.io' -e 'is available - Please consider' "$1" \
        | sed '/^[[:space:]]*$/N;/^\n$/D' | head -n 30 | sed 's/^/    | /'
}

# expect_fail <name> <expected text> <nextflow args...>
#   Pipeline must exit non-zero AND print <expected text>
expect_fail() {
    local name="$1" expect="$2"; shift 2
    local log="$T/logs/${name}.log"
    "$NXF" run main.nf -ansi-log false "$@" > "$log" 2>&1
    local rc=$?
    if [[ $rc -ne 0 ]] && grep -qF -- "$expect" "$log"; then
        echo "PASS  $name"; PASS=$((PASS+1)); SUMMARY+=("PASS  $name")
    else
        echo "FAIL  $name  (exit=$rc, expected text: \"$expect\")"; FAIL=$((FAIL+1)); SUMMARY+=("FAIL  $name")
    fi
    show_msg "$log"; echo
}

# observe <name> <nextflow args...>
#   Runtime cases: just run and show what happens, for you to judge the message
observe() {
    local name="$1"; shift
    local log="$T/logs/${name}.log"
    # -ansi-log false: plain log lines; the live progress display uses cursor
    # escapes that overwrite earlier lines when the log is printed back
    "$NXF" run main.nf -ansi-log false -work-dir "$T/work" "$@" > "$log" 2>&1
    local rc=$?
    echo "SHOW  $name  (exit=$rc)"; SHOWN=$((SHOWN+1)); SUMMARY+=("SHOW  $name  exit=$rc")
    # for a failed run, the ERROR block is what matters; for a success, the tail
    if [[ $rc -ne 0 ]]; then
        # first ERROR block only, without the work-dir/container boilerplate
        sed -n '/ERROR ~/,$p' "$log" | sed '/^Work dir:/,$d' | head -n 40 | sed 's/^/    | /'
    else
        tail -n 15 "$log" | sed 's/^/    | /'
    fi
    echo
}

COMMON_FQ=(--data_type fastq --pool_ID errtest)
COMMON_BAM=(--data_type bam --pool_ID errtest)

# ---------------------------------------------------------------------------
# QUICK: input validation (resolveInput) - runs with -preview, no tasks
# ---------------------------------------------------------------------------
run_quick() {
    echo "=== QUICK: input checks (-preview, no Docker) ==="; echo
    local I="$T/inputs"
    rm -rf "$I"/{empty,other,nested,mixed}
    mkdir -p "$I"/empty "$I"/other "$I"/nested/barcode01 "$I"/nested/barcode02
    echo "not reads" > "$I/other/notes.txt"
    cp "$FQ_DIR"/*_3.fastq.gz "$I/nested/barcode01/a.fastq.gz"
    cp "$FQ_DIR"/*_4.fastq.gz "$I/nested/barcode02/b.fastq.gz"
    echo "hello" > "$I/readme.txt"

    local O=(--outdir "$T/out_quick" -preview)

    expect_fail q01_missing_input      "--input is missing"            "${COMMON_FQ[@]}" "${O[@]}"
    expect_fail q02_folder_not_found   "--input not found"             "${COMMON_FQ[@]}" "${O[@]}" --input "$I/does_not_exist/"
    expect_fail q03_file_not_found     "--input not found"             "${COMMON_FQ[@]}" "${O[@]}" --input "$I/nope.fastq.gz"
    expect_fail q04_empty_folder       "(it is empty)"                 "${COMMON_FQ[@]}" "${O[@]}" --input "$I/empty"
    expect_fail q05_no_read_files      "contains no"                   "${COMMON_FQ[@]}" "${O[@]}" --input "$I/other"
    expect_fail q06_bam_given_as_fastq "Did you mean --data_type bam"  "${COMMON_FQ[@]}" "${O[@]}" --input "$BAM_DIR"
    expect_fail q07_fastq_given_as_bam "Did you mean --data_type fastq" "${COMMON_BAM[@]}" "${O[@]}" --input "$FQ_DIR"
    expect_fail q08_single_wrong_type  "Did you mean --data_type bam"  "${COMMON_FQ[@]}" "${O[@]}" --input "$(ls "$BAM_DIR"/*.bam | head -1)"
    expect_fail q09_not_a_read_file    "is not a"                      "${COMMON_FQ[@]}" "${O[@]}" --input "$I/readme.txt"
    expect_fail q10_reads_in_subfolders "sub-folders do"               "${COMMON_FQ[@]}" "${O[@]}" --input "$I/nested"
    expect_fail q11_glob_no_match      "No files match"                "${COMMON_FQ[@]}" "${O[@]}" --input "$I/x*/*.fastq.gz"
    expect_fail q12_bad_data_type      "Invalid --data_type"           --data_type fasta --pool_ID errtest "${O[@]}" --input "$FQ_DIR"

    # start-up checks on file names and the output folder
    rm -rf "$I"/{spaces,dups}; mkdir -p "$I/spaces" "$I/dups"
    cp "$FQ_DIR"/*_3.fastq.gz "$I/spaces/my reads.fastq.gz"
    cp "$FQ_DIR"/*_3.fastq.gz "$I/dups/sampleA.fastq.gz"; gzip -cd "$FQ_DIR"/*_3.fastq.gz > "$I/dups/sampleA.fastq"
    expect_fail q13_space_in_name      "must not contain spaces"       "${COMMON_FQ[@]}" "${O[@]}" --input "$I/spaces"
    expect_fail q14_duplicate_names    "same sample name"              "${COMMON_FQ[@]}" "${O[@]}" --input "$I/dups"
    expect_fail q15_missing_outdir     "--outdir is missing"           "${COMMON_FQ[@]}" -preview --input "$FQ_DIR"
    expect_fail q18_bad_flag_value     "must be true or false"         "${COMMON_FQ[@]}" "${O[@]}" --input "$FQ_DIR" --check_inputs maybe

    # which files the input check reads, by option (must succeed and log it)
    expect_log() {   # expect_log <name> <expected log text> <nextflow args...>
        local name="$1" expect="$2"; shift 2
        local log="$T/logs/${name}.log"
        if "$NXF" run main.nf -ansi-log false "$@" > "$log" 2>&1 && grep -qF -- "$expect" "$log"; then
            echo "PASS  $name"; PASS=$((PASS+1)); SUMMARY+=("PASS  $name")
        else
            echo "FAIL  $name (expected the run to succeed and log: \"$expect\")"
            FAIL=$((FAIL+1)); SUMMARY+=("FAIL  $name"); show_msg "$log"
        fi
    }
    expect_log q17_check_off_default   "Input check off (default)"     "${COMMON_FQ[@]}" "${O[@]}" --input "$FQ_DIR"
    expect_log q19_check_inputs        "Input check: reading"          "${COMMON_FQ[@]}" "${O[@]}" --input "$FQ_DIR" --check_inputs
    expect_log q20_check_all_explicit  "reading all"                   "${COMMON_FQ[@]}" "${O[@]}" --input "$FQ_DIR" --check_inputs false --check_all_inputs true
    echo

    mkdir -p "$I/readonly"; chmod a-w "$I/readonly"
    expect_fail q16_unwritable_outdir  "cannot write to --outdir"      "${COMMON_FQ[@]}" -preview --input "$FQ_DIR" --outdir "$I/readonly/sub"
    chmod u+w "$I/readonly"

    # valid inputs must NOT fail at the input stage
    for ok in "$FQ_DIR" "$FQ_DIR/" "$(ls "$FQ_DIR"/*.fastq.gz | head -1)" "$I/nested/*/*.fastq.gz"; do
        local log="$T/logs/q_ok_$(echo "$ok" | tr -c 'A-Za-z0-9' '_' | tail -c 40).log"
        if "$NXF" run main.nf "${COMMON_FQ[@]}" "${O[@]}" --input "$ok" > "$log" 2>&1; then
            echo "PASS  valid input accepted: $ok"; PASS=$((PASS+1)); SUMMARY+=("PASS  valid: $ok")
        else
            echo "FAIL  valid input rejected: $ok"; FAIL=$((FAIL+1)); SUMMARY+=("FAIL  valid: $ok")
            show_msg "$log"
        fi
    done
    echo
}

# ---------------------------------------------------------------------------
# FULL: runtime failures - runs real tasks in Docker. These show you the
# message; decide whether each one is clear enough for a customer.
# ---------------------------------------------------------------------------
run_full() {
    echo "=== FULL: runtime failures (real run, needs Docker) ==="; echo
    docker info >/dev/null 2>&1 || { echo "Docker is not running - start Docker Desktop first." >&2; return; }
    local R="$T/inputs/runtime"; rm -rf "$R"; mkdir -p "$R"
    local fq1; fq1=$(ls "$FQ_DIR"/*.fastq.gz | head -1)
    local bam1; bam1=$(ls "$BAM_DIR"/*.bam | head -1)

    # r01: not gzip at all, just named .fastq.gz
    mkdir -p "$R/r01"; head -c 5000 /dev/urandom > "$R/r01/garbage.fastq.gz"
    # r02: truncated gzip (download cut off half way)
    mkdir -p "$R/r02"; head -c 200000 "$fq1" > "$R/r02/truncated.fastq.gz"
    # r03: valid gzip, broken FASTQ records (a line missing)
    mkdir -p "$R/r03"; gzip -cd "$fq1" | head -n 4001 | sed '4000d' | gzip > "$R/r03/broken_records.fastq.gz"
    # r04: zero reads (empty gzip)
    mkdir -p "$R/r04"; printf '' | gzip > "$R/r04/empty.fastq.gz"
    # r05: valid reads with no seqWell barcodes (random sequence) -> all unclassified
    mkdir -p "$R/r05"
    python3 - "$R/r05/no_barcodes.fastq.gz" <<'PY'
import gzip, random, sys, uuid
random.seed(1)
with gzip.open(sys.argv[1], "wt") as f:
    for _ in range(200):
        s = "".join(random.choice("ACGT") for _ in range(800))
        f.write(f"@{uuid.uuid4()} runid=x ch=1\n{s}\n+\n{'I'*800}\n")
PY
    # r06: several bad files at once, plus one good one -> all bad ones listed
    mkdir -p "$R/r06"; cp "$fq1" "$R/r06/good.fastq.gz"
    cp "$R/r02/truncated.fastq.gz" "$R/r04/empty.fastq.gz" "$R/r06/"
    # r07: corrupt BAM
    mkdir -p "$R/r07"; head -c 50000 "$bam1" > "$R/r07/truncated.bam"
    # r08: empty BAM (header only)
    mkdir -p "$R/r08"
    docker run --rm --platform linux/amd64 -v "$BAM_DIR":/in:ro -v "$R/r08":/out \
        quay.io/biocontainers/samtools:1.21--h50ea8bc_0 \
        sh -c "samtools view -H /in/$(basename "$bam1") | samtools view -b -o /out/header_only.bam -"

    local W=(-work-dir "$T/work")
    expect_fail r01_not_gzip        "not a valid gzip file"         "${COMMON_FQ[@]}"  "${W[@]}" --input "$R/r01" --outdir "$T/out_r01"
    expect_fail r02_truncated_gzip  "ends unexpectedly"             "${COMMON_FQ[@]}"  "${W[@]}" --input "$R/r02" --outdir "$T/out_r02"
    expect_fail r03_broken_records  "different lengths"             "${COMMON_FQ[@]}"  "${W[@]}" --input "$R/r03" --outdir "$T/out_r03"
    expect_fail r04_zero_reads      "contains no reads"             "${COMMON_FQ[@]}"  "${W[@]}" --input "$R/r04" --outdir "$T/out_r04"
    expect_fail r06_several_bad     "2 input file(s) cannot be used" "${COMMON_FQ[@]}" "${W[@]}" --input "$R/r06" --outdir "$T/out_r06" --check_inputs
    expect_fail r07_corrupt_bam     "not a valid BAM file"          "${COMMON_BAM[@]}" "${W[@]}" --input "$R/r07" --outdir "$T/out_r07"
    expect_fail r08_empty_bam       "contains no reads"             "${COMMON_BAM[@]}" "${W[@]}" --input "$R/r08" --outdir "$T/out_r08"

    # r10: bad file in the MIDDLE of the list - only caught with --check_all_inputs
    mkdir -p "$R/r10"; cp "$fq1" "$R/r10/run_1.fastq.gz"; cp "$R/r02/truncated.fastq.gz" "$R/r10/run_2.fastq.gz"
    cp "$fq1" "$R/r10/run_10.fastq.gz"
    expect_fail r10_middle_file_all "run_2.fastq.gz: cannot be read" "${COMMON_FQ[@]}" "${W[@]}" --input "$R/r10" --outdir "$T/out_r10" --check_all_inputs

    # r05 should SUCCEED with a 0% report and a low-demux warning
    local log="$T/logs/r05_no_barcodes.log"
    "$NXF" run main.nf -ansi-log false "${W[@]}" "${COMMON_FQ[@]}" --input "$R/r05" --outdir "$T/out_r05" > "$log" 2>&1
    local rc=$?
    if [[ $rc -eq 0 ]] && grep -q "of reads were assigned to a barcode" "$log"; then
        echo "PASS  r05_no_barcodes (succeeds, warns about low demux rate)"; PASS=$((PASS+1)); SUMMARY+=("PASS  r05_no_barcodes")
    else
        echo "FAIL  r05_no_barcodes (exit=$rc; expected success + low-demux warning)"; FAIL=$((FAIL+1)); SUMMARY+=("FAIL  r05_no_barcodes")
    fi
    grep -A1 "WARN" "$log" | sed 's/^/    | /'; echo

    # r09: missing barcode FASTA - no pre-check by design; shows Dorado's own message
    observe r09_missing_barcodes "${COMMON_FQ[@]}" --input "$FQ_DIR" --outdir "$T/out_r09" --barcodes "$R/no_such_barcodes.fa"
}

case "$MODE" in
    quick) run_quick ;;
    full)  run_full ;;
    all)   run_quick; run_full ;;
    *) echo "Usage: $0 [quick|full|all]" >&2; exit 1 ;;
esac

echo "=== SUMMARY ==="
printf '%s\n' "${SUMMARY[@]}"
echo "passed=$PASS failed=$FAIL shown=$SHOWN   (full logs: error_tests/logs/)"
[[ $FAIL -eq 0 ]]
