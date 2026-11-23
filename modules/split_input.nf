// Split a single FASTQ (.fastq or .fastq.gz) into n_parts roughly equal chunks.
// seqkit reads plain or gzipped input directly; --extension .gz forces gzipped chunks.
// Output: chunks/<sample_id>.part_001.fastq.gz ...
// Container: quay.io/biocontainers/seqkit
//
// Also serves as the input check for a single split file: the split reads
// every record anyway, so instead of a separate CHECK_FASTQ pass it writes
// <sample_id>.check.tsv (OK/FAIL + plain-language reason). The task itself
// never fails on a bad file; main.nf stops the run from the report.
process SPLIT_FASTQ {
    tag "${sample_id}"
    cpus 4

    input:
    tuple val(sample_id), path(fq)
    val n_parts

    output:
    tuple val(sample_id), path("chunks/*.part_*.gz"), emit: chunks, optional: true
    path("${sample_id}.check.tsv"),                     emit: report

    script:
    """
    report=${sample_id}.check.tsv

    # quick check first (gzip header, not empty, first read parses)
    check_input.sh --mode fastq --level quick --out quick.tsv ${fq}

    if awk -F'\\t' '\$2 == "FAIL" {f = 1} END {exit !f}' quick.tsv; then
        mv quick.tsv \$report
    elif ! seqkit split2 \\
            --by-part ${n_parts} \\
            --threads ${task.cpus} \\
            --extension .gz \\
            --out-dir chunks \\
            ${fq} 2> split.err \\
         || grep -q 'ERRO' split.err; then
        # seqkit can exit 0 on a read error, so its stderr is checked too
        cat split.err >&2
        rm -rf chunks
        check_input.sh --mode fastq --level explain --errlog split.err --out \$report ${fq}
    elif ! ls chunks/*.part_*.gz > /dev/null 2>&1; then
        printf '%s\\tFAIL\\t0\\tcontains no reads.\\n' "${fq}" > \$report
    else
        printf '%s\\tOK\\t-\\t\\n' "${fq}" > \$report
    fi
    """
}


// Split a single unaligned BAM into n_parts chunks, round-robin by record,
// each chunk carrying the full header (RG, PG, etc.). Aux tags (MM/ML) preserved.
// One pass, no read counting; each chunk is written by its own samtools process.
// Output: chunks/<sample_id>.part_000.bam ...
// Container: quay.io/biocontainers/samtools:1.21 (busybox base: uses awk, not GNU split)
//
// Also serves as the input check for a single split file (see SPLIT_FASTQ).
process SPLIT_BAM {
    tag "${sample_id}"
    cpus 4

    input:
    tuple val(sample_id), path(bam)
    val n_parts

    output:
    tuple val(sample_id), path("chunks/*.part_*.bam"), emit: chunks, optional: true
    path("${sample_id}.check.tsv"),                    emit: report

    script:
    """
    set -o pipefail
    report=${sample_id}.check.tsv

    # quick check first (header + EOF marker, at least one read)
    check_input.sh --mode bam --level quick --out quick.tsv ${bam}

    if awk -F'\\t' '\$2 == "FAIL" {f = 1} END {exit !f}' quick.tsv; then
        mv quick.tsv \$report
    else
        mkdir -p chunks
        samtools view -H ${bam} > header.sam

        if samtools view -@ 2 ${bam} 2> split.err \\
            | awk -v n=${n_parts} -v pre="chunks/${sample_id}.part_" '
                {
                    r = (NR - 1) % n
                    if (!(r in cmd)) cmd[r] = sprintf("cat header.sam - | samtools view -b -o %s%03d.bam -", pre, r)
                    print | cmd[r]
                }
                END { for (r in cmd) close(cmd[r]) }' \\
           && ! grep -qiE 'truncated|error' split.err; then
            printf '%s\\tOK\\t-\\t\\n' "${bam}" > \$report
        else
            cat split.err >&2
            rm -rf chunks
            printf '%s\\tFAIL\\t-\\tcannot be read: the BAM file is corrupt or incomplete. Re-copy or re-download the file.\\n' "${bam}" > \$report
        fi
    fi
    """
}
