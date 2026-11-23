
process READ_LENGTH {
    tag "$pair_id"
    errorStrategy 'ignore'
    publishDir "${params.outdir}/read_length", pattern: '*weighted*',              mode: 'copy'
    publishDir "${params.outdir}/read_length", pattern: '*.read_length_plot.png',  mode: 'copy'

    input:
    tuple val(pair_id), path(fq)

    output:
    tuple val(pair_id), path("${pair_id}.read_length_matrix.txt"), path("${pair_id}.read_length_plot.png")
    path("${pair_id}.read_length_matrix.txt"), emit: matrix
    path("*bin_table.txt"),                    emit: bin
    path("*weighted*"),                        emit: weighted_plot

    script:
    // One pass over the FASTQ: count reads per length (length<TAB>count).
    // The table has one line per distinct length (a few thousand at most),
    // so only this small table is sorted - not one line per read.
    """
    zcat ${fq} \\
        | awk 'NR % 4 == 2 { c[length(\$0)]++ } END { for (l in c) print l "\\t" c[l] }' \\
        | sort -n -k1,1 > ${pair_id}.length_counts.tsv

    # Fixed 1..16000 bp matrix (length<TAB>count), as before
    awk 'BEGIN { for (i = 1; i <= 16000; i++) c[i] = 0 } { c[\$1] = \$2 } END { for (i = 1; i <= 16000; i++) print i "\\t" c[i] }' \\
        ${pair_id}.length_counts.tsv > ${pair_id}.read_length_matrix.txt

    read_length.py ${pair_id}.length_counts.tsv ${pair_id}.read_length_plot.png
    """
}
