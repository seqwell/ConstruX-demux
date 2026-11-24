/*
 * subworkflows/check_inputs.nf
 *
 * CHECK_INPUTS: everything that has to happen before Dorado starts.
 *
 *   1. Start-up checks (run immediately, before any task):
 *        --data_type, --input (folder / file / quoted glob), file names
 *        (no spaces, no duplicate sample names), option combinations,
 *        and that --outdir can be written.
 *   2. Integrity check of the input files, and split of a single file:
 *        One file + split_n > 1:
 *          SPLIT_FASTQ / SPLIT_BAM read the whole file anyway, so they
 *          also act as its check. Chunk IDs look like "<sample>.part_001".
 *        Several files, or one file with split_n 1 (OFF by default):
 *          with --check_inputs the first and last file (natural sort,
 *          run_9 before run_10) are read end to end by CHECK_FASTQ /
 *          CHECK_BAM; with --check_all_inputs every file is.
 *          Without either, the files go straight to Dorado.
 *   3. Any FAIL stops the run with one message listing every bad file.
 *
 * Emits tuple(sample_id, file) (CHECK_INPUTS.out), released only after
 * the check has passed.
 */

include { CHECK_FASTQ } from '../modules/check_fastq.nf'
include { CHECK_BAM   } from '../modules/check_bam.nf'
include { SPLIT_FASTQ } from '../modules/split_input.nf'
include { SPLIT_BAM   } from '../modules/split_input.nf'
include {
    sampleIdOf; asList; resolveInput; checkInputNames; checkOutdir;
    checkDataType; failOnBadInputs; naturalKey; asBool
} from '../functions/input_checks.nf'

workflow CHECK_INPUTS {

    take:
    input              // --input: folder, file or quoted glob
    data_type          // 'fastq' or 'bam'
    outdir             // --outdir
    split_n            // number of chunks for a single input file
    check_inputs       // boolean: read first + last file before Dorado
    check_all_inputs   // boolean: read every file before Dorado

    main:

    // Copy the inputs into local variables before anything else.
    // Nextflow 24 fails to compile ("Variable already defined") when a
    // take: input is used in a statement and then again in a later `def`.
    def input_spec = input
    def dtype      = data_type
    def out_dir    = outdir
    def n_split    = split_n
    def check_on   = asBool(check_inputs, 'check_inputs')
    def check_all  = asBool(check_all_inputs, 'check_all_inputs')

    // ---------------------------------------------------------------
    // 1. Start-up checks
    // ---------------------------------------------------------------
    checkDataType(dtype)
    def use_bam = (dtype == 'bam')

    def input_files = resolveInput(input_spec, use_bam, dtype)
    checkInputNames(input_files)
    checkOutdir(out_dir)

    def n_parts  = n_split as Integer
    def do_split = (input_files.size() == 1 && n_parts > 1)

    if (use_bam) {
        log.info "data_type=bam  →  BAM input: dorado demux (BAM) → BAM-to-FASTQ → cutadapt → filter demux-BAM by cutadapt-FASTQ read IDs"
    } else {
        log.info "data_type=fastq  →  FASTQ input (.fastq or .fastq.gz): extract headers → dorado demux → reheader → cutadapt → FASTQ output only (no BAM created)"
    }
    log.info "Found ${input_files.size()} input file(s)" +
             (do_split ? " → single input, splitting into ${n_parts} chunks for dorado" : " → no splitting")

    def unchecked_ch = channel.fromList(input_files)
                         .map { f -> tuple(sampleIdOf(f), f) }

    // ---------------------------------------------------------------
    // 2. Integrity check (and split for a single file)
    //    files_ch : tuple(sample_id, file) to demultiplex
    //    ok_ch    : emits once, after every check report has passed
    // ---------------------------------------------------------------
    def files_ch = channel.empty()
    def ok_ch    = channel.empty()

    if (do_split) {
        def chunk_out     = channel.empty()
        def split_reports = channel.empty()
        if (use_bam) {
            SPLIT_BAM(unchecked_ch, n_parts)
            chunk_out     = SPLIT_BAM.out.chunks
            split_reports = SPLIT_BAM.out.report
        } else {
            SPLIT_FASTQ(unchecked_ch, n_parts)
            chunk_out     = SPLIT_FASTQ.out.chunks
            split_reports = SPLIT_FASTQ.out.report
        }
        files_ch = chunk_out
                     .flatMap { _id, chunks -> asList(chunks) }
                     .map { c -> tuple(sampleIdOf(c), c) }
        ok_ch    = split_reports.collect().map { reports -> failOnBadInputs(reports) }

    } else if (!check_on && !check_all) {
        log.info "Input check off (default) - add --check_inputs (first + last file) or --check_all_inputs (every file) to read the inputs before demultiplexing"
        files_ch = unchecked_ch
        ok_ch    = channel.of(true)

    } else {
        def ordered    = input_files.toSorted { a, b -> naturalKey(a) <=> naturalKey(b) }
        def full_files = check_all ? ordered : [ordered.first(), ordered.last()].unique()
        if (check_all) {
            log.info "Input check: reading all ${full_files.size()} input files end to end (--check_all_inputs)"
        } else {
            log.info "Input check: reading ${full_files.collect { f -> f.name }.join(' and ')} end to end"
        }

        def check_in = channel.fromList(full_files)
                         .map { f -> tuple('full', f) }
        def check_reports = channel.empty()
        if (use_bam) {
            CHECK_BAM(check_in)
            check_reports = CHECK_BAM.out.report
        } else {
            CHECK_FASTQ(check_in)
            check_reports = CHECK_FASTQ.out.report
        }
        files_ch = unchecked_ch
        ok_ch    = check_reports.collect().map { reports -> failOnBadInputs(reports) }
    }

    emit:
    // 3. files are released only after ok_ch (any FAIL stops the run there)
    files_ch
        .combine(ok_ch)
        .map { id, f, _ok -> tuple(id, f) }
}
