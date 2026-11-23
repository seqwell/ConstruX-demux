#!/usr/bin/env nextflow

include { CHECK_INPUTS         } from './subworkflows/check_inputs.nf'
include { sampleIdOf; asList; asBool } from './functions/input_checks.nf'
include { DORADO_DEMUX         } from './modules/dorado_demux.nf'
include { COMBINE_BARCODES     } from './modules/combine_barcodes.nf'
include { COMBINE_BARCODES_BAM } from './modules/combine_barcodes.nf'
include { BAM_TO_FASTQ         } from './modules/bam_to_fastq.nf'
include { EXTRACT_HEADER       } from './modules/extract_header.nf'
include { REHEADER_READS       } from './modules/reheader_reads.nf'
include { CUTADAPT_TRIM        } from './modules/cutadapt_trim.nf'
include { FILTER_BAM_BY_FASTQ  } from './modules/filter_bam_by_fastq.nf'
include { DEMUX_SUMMARIZE      } from './modules/demux_summarize.nf'
include { NANOSTAT             } from './modules/nanostat.nf'
include { MULTIQC              } from './modules/multiQC.nf'
include { READ_LENGTH          } from './modules/read_length.nf'
include { NANOSTAT_BASES       } from './modules/nanostat_bases.nf'
include { MERGE_DEMUX_BASES    } from './modules/merge_demux_bases.nf'


// Helper functions (sampleIdOf, asList, input checks): functions/input_checks.nf

workflow {

    def wf = workflow
    def p  = params
    def lg = log
    def nf_version = nextflow.version
    
    def all_checked = asBool(params.check_all_inputs, 'check_all_inputs')
    wf.onComplete {

        def ok         = wf.success
        def stats      = wf.stats
        def bar        = '=' * 72
        def outdir     = p.outdir ? file(p.outdir).toUriString() : "${wf.launchDir}/results"
        def cmd        = wf.commandLine
        def resume_cmd = cmd.contains('-resume') ? cmd : "${cmd} -resume"

        def summary = [
            '',
            bar,
            "  seqWell ONT Dorado Demux  —  ${ok ? 'COMPLETED SUCCESSFULLY ✅' : 'FAILED ❌'}",
            bar,
            "  Run name        : ${wf.runName}",
            "  Session ID      : ${wf.sessionId}",
            "  Pipeline rev    : ${wf.revision ?: wf.commitId ?: 'local / no git'}",
            "  Nextflow        : ${nf_version}",
            "  Profile         : ${wf.profile}",
            "  Container       : ${wf.containerEngine ?: 'none'}",
            '',
            "  Started         : ${wf.start}",
            "  Finished        : ${wf.complete}",
            "  Wall time       : ${wf.duration}",
            "  CPU time        : ${stats.computeTimeFmt}",
            '',
            "  Tasks           : ${stats.succeededCount} succeeded, ${stats.cachedCount} cached, " +
                "${stats.failedCount} failed, ${stats.ignoredCount} ignored",
            '',
            "  Input settings",
            "    data_type     : ${p.data_type}",
            "    input         : ${p.input}",
            "    split_n       : ${p.split_n}",
            "    barcodes      : ${p.barcodes}",
            "    arrangement   : ${p.arrangement_toml}",
            '',
            "  Output dir      : ${outdir}",
            "  Work dir        : ${wf.workDir}",
            "  Launch dir      : ${wf.launchDir}",
            '',
            "  Command line    :",
            "    ${cmd}",
        ]

        // Our own input/output checks already printed a complete message;
        // do not repeat it (three times) in the summary
        def own_check_error = (wf.errorMessage ?: '').contains(' problem:')

        if (!ok && own_check_error) {
            summary = summary + [
                '',
                '  Error           : stopped by an input/output check - see the "problem:" message printed with this run',
            ]
        } else if (!ok) {
            def report_lines = wf.errorReport
                ? wf.errorReport.readLines().take(15).collect { l -> "    ${l}" }
                : []
            def error_lines = wf.errorMessage
                ? ["  Error           : ${wf.errorMessage}"]
                : []
            summary = summary + [''] + error_lines + [
                "  Error report    :",
            ] + report_lines + [
                '',
                "  To resume from the last successful step:",
                "    ${resume_cmd}",
            ] + (all_checked ? [] : [
                '',
                "  Tip             : if a step failed while reading an input file, some input",
                "                    files may be corrupt or incomplete. Rerun with --check_all_inputs",
                "                    to read every input file first and list any bad ones.",
            ])
        }

        // Where the debug files are saved (see "Run reports" in nextflow.config)
        def info_dir = "${p.outdir ?: wf.launchDir}/pipeline_info"
        summary = summary + [
            '',
            "  Logs & reports  : ${info_dir}/",
            "                    nextflow.log, run_summary.txt, execution_report.html (CPU/RAM per step),",
            "                    execution_trace.txt, execution_timeline.html" + (ok ? '' : ', error_report.txt'),
            bar, ''
        ]

        // Printed at info level even for a failed run: if the summary were
        // logged as an error, Nextflow 26 would skip printing the error
        // message itself for problems found at start-up.
        lg.info summary.join('\n')

        // ---------------------------------------------------------------
        // Save debug files next to the results, so a customer can send one
        // folder (pipeline_info/) when asking for help. Never let a problem
        // here hide the real result of the run.
        // ---------------------------------------------------------------
        try {
            def info = file(info_dir)
            info.mkdirs()

            // the summary above, as plain text
            info.resolve('run_summary.txt').text = summary.join('\n') + '\n'

            // full error report (all lines, not just the first 15 shown above)
            if (!ok) {
                info.resolve('error_report.txt').text = [
                    "Error message:",
                    wf.errorMessage ?: '(none)',
                    '',
                    "Error report:",
                    wf.errorReport ?: '(none)',
                    ''
                ].join('\n')
            }

            // copy of the Nextflow log (the default .nextflow.log in the launch
            // directory; skipped if the run used a custom -log path)
            def nflog = file("${wf.launchDir}/.nextflow.log")
            if (nflog.exists()) {
                nflog.copyTo(info.resolve('nextflow.log'))
            }
        } catch (Exception e) {
            lg.warn "Could not save debug files to ${info_dir}: ${e.message}"
        }
    }


    // ---------------------------------------------------------------
    // Step 0: Check inputs and outdir, split a single large file
    //   (subworkflows/check_inputs.nf). Emits tuple(sample_id, file)
    //   only after every check has passed.
    // ---------------------------------------------------------------
    CHECK_INPUTS(
        params.input,
        params.data_type,
        params.outdir,
        params.split_n,
        params.check_inputs,
        params.check_all_inputs
    )
    def input_ch = CHECK_INPUTS.out
    def use_bam  = (params.data_type == 'bam')

    def barcode_fasta    = file(params.barcodes)
    def arrangement_toml = file(params.arrangement_toml)


    // ---------------------------------------------------------------
    // Step 1: Demux (one job per input file, or per chunk)
    //   bam input   → dorado emits BAM dir only (--no emit-fastq)
    //   fastq input → dorado emits FASTQ dir only (--emit-fastq)
    // ---------------------------------------------------------------
    DORADO_DEMUX(input_ch, barcode_fasta, arrangement_toml)

    // Set inside whichever branch runs, consumed by the shared QC steps
    def trimmed_fq_ch = channel.empty()   // raw CUTADAPT_TRIM.out.fq
    def demuxed_ch    = channel.empty()   // tuple(sample_id, trimmed_fastq)


    // ================================================================
    //  BAM INPUT BRANCH
    // ================================================================
    if (use_bam) {

        // -----------------------------------------------------------
        // Step 2 (BAM): Combine per-sample BAM demux dirs → per-barcode BAMs
        // -----------------------------------------------------------
        COMBINE_BARCODES_BAM(DORADO_DEMUX.out.bam_dir.collect())

        def bam_barcode_ch = COMBINE_BARCODES_BAM.out.bam
                               .flatMap { bams -> asList(bams) }
                               .filter { bam -> bam.name.endsWith('.bam') }
                               .map { bam -> tuple(bam.baseName, bam) }

        // -----------------------------------------------------------
        // Step 3 (BAM): Convert each per-barcode demuxed BAM → FASTQ
        //   This FASTQ is used as input to cutadapt.
        //   No header extraction is needed for BAM input.
        // -----------------------------------------------------------
        BAM_TO_FASTQ(bam_barcode_ch)

        // -----------------------------------------------------------
        // Step 4 (BAM): Cutadapt trim on the converted FASTQ
        // -----------------------------------------------------------
        CUTADAPT_TRIM(BAM_TO_FASTQ.out.fastq)
        trimmed_fq_ch = CUTADAPT_TRIM.out.fq

        // trimmed_ch keyed by bare barcode ID
        // e.g. "barcode001.seqWell.fastq.gz" → key "barcode001"
        // This must match bam_barcode_ch which is also keyed by "barcode001"
        def base_trimmed_ch = trimmed_fq_ch
                                .flatMap { fqs -> asList(fqs) }
                                .filter { fq -> !fq.name.contains('tagged') && fq.size() > 20 }

        def trimmed_ch = base_trimmed_ch
                           .map { fq -> tuple(fq.baseName.replace('.seqWell.fastq', ''), fq) }

        def trimmed_ch_for_nanostat = base_trimmed_ch
                                        .map { fq -> tuple(fq.baseName.replace('.fastq', ''), fq) }

        // -----------------------------------------------------------
        // Step 5 (BAM): Split the demux-BAM to match cutadapt-FASTQ results
        //   Join on bare barcode ID so the keys match on both sides:
        //     bam_barcode_ch key: "barcode001"
        //     trimmed_ch key:     "barcode001"
        //   Then subset the BAM to only reads whose names appear in the FASTQ,
        //   propagating cutadapt adapter/length/ME-tag filtering onto the BAM.
        // -----------------------------------------------------------
        def bam_fastq_ch = bam_barcode_ch.join(trimmed_ch, by: 0)
                             // produces: tuple(barcode_id, bam, fastq)

        FILTER_BAM_BY_FASTQ(bam_fastq_ch)

        demuxed_ch = trimmed_ch_for_nanostat

    // ================================================================
    //  FASTQ INPUT BRANCH
    // ================================================================
    } else {

        // -----------------------------------------------------------
        // Step 2 (FASTQ): Extract FASTQ headers BEFORE demuxing
        //   awk extracts FASTQ header fields → uuid_tags.tsv
        //   Consumed by REHEADER_READS after demux to restore original headers.
        //   Runs per input file, or per chunk when the input was split.
        // -----------------------------------------------------------
        EXTRACT_HEADER(input_ch)
        def merged_tags = EXTRACT_HEADER.out.collect()

        // -----------------------------------------------------------
        // Step 3 (FASTQ): Combine per-sample FASTQ demux dirs → per-barcode FASTQs
        // -----------------------------------------------------------
        COMBINE_BARCODES(DORADO_DEMUX.out.fastq_dir.collect())

        def fastq_barcode_ch = COMBINE_BARCODES.out.fastq
                                 .flatMap { fqs -> asList(fqs) }
                                 .map { fq -> tuple(fq.baseName.replace('.fastq', ''), fq) }

        // -----------------------------------------------------------
        // Step 4 (FASTQ): Restore original header tags into demuxed FASTQs
        // -----------------------------------------------------------
        REHEADER_READS(fastq_barcode_ch, merged_tags)

        // -----------------------------------------------------------
        // Step 5 (FASTQ): Cutadapt trim — on reheadered FASTQ
        //   No BAM files are created for FASTQ input.
        // -----------------------------------------------------------
        CUTADAPT_TRIM(REHEADER_READS.out.fq)
        trimmed_fq_ch = CUTADAPT_TRIM.out.fq

        demuxed_ch = trimmed_fq_ch
                       .flatMap { fqs -> asList(fqs) }
                       .map { fq -> tuple(fq.baseName.replace('.fastq', ''), fq) }
    }


    // ---------------------------------------------------------------
    // Step 6: Downstream QC — runs on trimmed FASTQ (both branches)
    // ---------------------------------------------------------------
    DEMUX_SUMMARIZE(trimmed_fq_ch.collect())

    // Warn (do not fail) when the overall demux rate looks wrong
    DEMUX_SUMMARIZE.out.subscribe { csv ->
        def line = csv.readLines().find { l -> l.startsWith('Total Demuxed Reads:') }
        def pct  = line ? (line.split(',')[2].replace('%', '').trim() as Double) : 0.0
        if (pct < 50) {
            log.warn "Only ${pct}% of reads were assigned to a barcode (see ${csv.name}). " +
                     'Check that these are seqWell libraries and that the barcode set is correct.'
        }
    }

    def valid_ids_ch = DEMUX_SUMMARIZE.out
                         .splitCsv()
                         .filter { row -> row[0].startsWith('barcode') || row[0] == 'unknown' }
                         .map { row -> tuple(row[0].trim(), true) }

    def filtered_ch = demuxed_ch
                        .join(valid_ids_ch, by: 0)
                        .map { sample_id, fq, _flag -> tuple(sample_id, fq) }

    READ_LENGTH(filtered_ch)
    NANOSTAT(filtered_ch)

    NANOSTAT_BASES(NANOSTAT.out.collect())

    def mqc_config = file("${projectDir}/assets/multiqc_config.yaml")
    MULTIQC(
        NANOSTAT.out.collect().mix(NANOSTAT_BASES.out.q10_bases).collect(),
        mqc_config
    )

    MERGE_DEMUX_BASES(DEMUX_SUMMARIZE.out, NANOSTAT_BASES.out.bases)


    // ===============================================================
    //  Completion summary — runs on success AND failure
    // ===============================================================
}
