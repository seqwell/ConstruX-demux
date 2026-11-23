/*
 * COMBINE_BARCODES / COMBINE_BARCODES_BAM
 *
 * Merge the per-input dorado demux directories into one file per barcode.
 * bin/combine_demux.sh searches each directory recursively, so it handles
 * both the flat layout of dorado 1.1.x and the nested ONT Output
 * Specifications (26.01) layout written by dorado >= 1.3.1 / 2.x,
 * and both plain and gzipped FASTQ.
 *
 * Output names are unchanged: barcode001.fastq.gz ... / unclassified.fastq.gz
 * (or .bam + .bai for the BAM branch).
 */
process COMBINE_BARCODES {
    tag "combine_barcodes"

    // publishDir "${params.outdir}/merged_fastq", mode: 'copy'

    input:
    path dummy  // collected fastq_dir paths

    output:
    path("output_fastq/*"), emit: fastq

    script:
    """
    combine_demux.sh \\
        --mode fastq \\
        --pattern '*_Demux_fastq' \\
        --outdir output_fastq
    """
}

process COMBINE_BARCODES_BAM {
    tag "combine_barcodes_bam"

    input:
    path dummy  // collected bam_dir paths

    output:
    path("output_bam/*"), emit: bam

    script:
    """
    combine_demux.sh \\
        --mode bam \\
        --pattern '*_Demux_bam' \\
        --outdir output_bam
    """
}
