/*
 * CHECK_BAM
 *
 * Integrity check of input bam files before any demux work: the first
 * and last input file (natural sort) are read end to end by
 * bin/check_input.sh, one task each.
 * Never fails the task itself: it writes a report with OK/FAIL and a
 * plain-language reason per file, and main.nf stops the run with one clear
 * message listing every bad file. Skipped with --skip_input_check.
 */
process CHECK_BAM {
    tag "${files.name}"

    input:
    tuple val(level), path(files)

    output:
    path("*.check.tsv"), emit: report

    script:
    """
    check_input.sh --mode bam --level ${level} --threads ${task.cpus} \\
        --out ${level}_${task.index}.check.tsv \\
        ${files}
    """
}
