// CUTADAPT_TRIM
//   Step 1: length filter only (--length_filter). Dorado has already
//           trimmed the barcode and its flanking adapter sequences.
//   Step 2: reads that still contain a seqWell adapter (ME or P5, either
//           orientation, anywhere in the read, up to --error_rate) are set
//           aside as ME/P5-tagged. Most of these are chimeras, whose
//           second molecule would otherwise be counted under the wrong
//           barcode (crosstalk).
// The *.ME.tagged.fastq.gz file name is kept for compatibility; it holds
// all ME/P5-tagged reads.
process CUTADAPT_TRIM {
    tag "${pair_id}"
    publishDir path: "${params.outdir}/demuxed_fastq_flat/", mode: 'copy', pattern: "*seqWell*"
    publishDir path: "${params.outdir}/demuxed_fastq/", mode: 'copy', saveAs: { filename ->
        if (filename.contains("seqWell")) "${pair_id}/${filename}"
        else if (filename == "unknown.fastq.gz") filename
        else null
    }
    publishDir path: "${params.outdir}/demuxed_fastq", mode: 'copy', pattern: "unknown.fastq.gz"
    publishDir path: "${params.outdir}/other/ME_tagged_fastq", mode: 'copy', pattern: '*tag*'

    input:
    tuple val(pair_id), path(fq)

    output:
    // optional: outputs with zero reads are deleted in the script below
    path("*.ME.tagged.fastq.gz"),            optional: true
    path("{*seqWell,*.ME.tagged}.fastq.gz"), emit: fq, optional: true

    script:
    """
    # Step 1: Length filter only (no adapter trimming)
    cutadapt \\
        --minimum-length ${params.length_filter} \\
        --report=minimal \\
        -o ${pair_id}.step1.fastq.gz \\
        ${fq} > ${pair_id}.cutadapt_report.step1.txt

    # Step 2: Set aside reads that still contain a seqWell adapter anywhere:
    #   ME (AGATGTGTATAAGAGACAG) and its reverse complement,
    #   P5 (AATGATACGGCGACCACCGAGATCTACAC) and its reverse complement.
    #   Errors allowed = floor(error_rate x matched length): with 0.106,
    #   2 in a full ME, 3 in a full P5. -O 19: a match inside the read must be
    #   full length; a partial match at a read end must be at least 19 nt.
    #   --action=none: nothing is cut, reads are only sorted.
    cutadapt \\
        -b AGATGTGTATAAGAGACAG \\
        -b CTGTCTCTTATACACATCT \\
        -b AATGATACGGCGACCACCGAGATCTACAC \\
        -b GTGTAGATCTCGGTGGTCGCCGTATCATT \\
        -e ${params.error_rate} \\
        -O 19 \\
        --action=none \\
        --report=minimal \\
        --untrimmed-output ${pair_id}.seqWell.fastq.gz \\
        -o ${pair_id}.ME.tagged.fastq.gz \\
        ${pair_id}.step1.fastq.gz > ${pair_id}.cutadapt_report.txt

    rm *step*.fastq.gz

    # Step 3: Drop outputs with no reads, so empty (~20 B) gzip files are
    # never published. A barcode whose reads were all removed by the
    # length filter / adapter check simply has no file.
    for f in ${pair_id}.seqWell.fastq.gz ${pair_id}.ME.tagged.fastq.gz; do
        if [ "\$(gzip -cd "\$f" | head -c 1 | wc -c)" -eq 0 ]; then
            echo "No reads left in \$f - not writing it" >&2
            rm -f "\$f"
        fi
    done
    """
}
