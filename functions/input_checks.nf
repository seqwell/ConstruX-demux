/*
 * functions/input_checks.nf
 *
 * Helper functions for input/output checks, used by
 * subworkflows/check_inputs.nf (and sampleIdOf / asList by main.nf).
 * All values are passed in as arguments; nothing here reads params.
 */

// Strip .bam / .fastq / .fq / .fastq.gz / .fq.gz to get the sample (or chunk) ID
def sampleIdOf(f) {
    return f.name.replaceAll(/\.(bam|fastq|fq)(\.gz)?$/, '')
}

// Normalise a process output that may be a single path or a list of paths
def asList(x) {
    return (x instanceof List) ? x : [x]
}

// ---------------------------------------------------------------
// Input resolution with plain-language errors
//   --input may be a folder (local or S3), a single file, or a quoted glob
// ---------------------------------------------------------------
def describeType(use_bam) {
    return use_bam ? '.bam' : '.fastq / .fq / .fastq.gz / .fq.gz'
}

def inputError(lines) {
    def msg = (['', 'Input problem:'] + lines.collect { l -> '  ' + l } + [
        '',
        'Examples:',
        '  --input /path/to/fastq_pass/            (a folder)',
        '  --input /path/to/reads.fastq.gz         (one file)',
        "  --input '/path/to/run*/*.fastq.gz'      (a glob - keep the quotes)",
        ''
    ]).join('\n')
    error(msg)
}

def resolveInput(input, use_bam, data_type) {
    def want  = use_bam ? /\.bam$/ : /\.(fastq|fq)(\.gz)?$/
    def other = use_bam ? /\.(fastq|fq)(\.gz)?$/ : /\.bam$/
    def otherType = use_bam ? 'fastq' : 'bam'

    if (!input) {
        inputError(['--input is missing. Tell the pipeline where your reads are.'])
    }
    def inStr = input.toString()

    // Quoted glob, e.g. 'reads/*.fastq.gz'
    if (inStr.contains('*') || inStr.contains('?') || inStr.contains('{')) {
        def hits = files(inStr)
        def ok   = hits.findAll { f -> f.name ==~ /.*/ + want }
        if (ok) { return ok }
        if (hits.isEmpty()) {
            inputError(["No files match --input '${inStr}'.",
                        'Check the folder exists and the pattern is spelled correctly.'])
        }
        def nOther = hits.count { f -> f.name ==~ /.*/ + other }
        if (nOther) {
            inputError(["--input '${inStr}' matched ${nOther} ${otherType} file(s), but --data_type is '${data_type}'.",
                        "Did you mean --data_type ${otherType}?"])
        }
        inputError(["--input '${inStr}' matched ${hits.size()} file(s), but none end in ${describeType(use_bam)}."])
    }

    def p = file(inStr)
    if (!p.exists()) {
        inputError(["--input not found: ${inStr}",
                    'Check the path is spelled correctly and that this folder or file exists.'])
    }

    // Single file
    if (!p.isDirectory()) {
        if (p.name ==~ /.*/ + want) { return [p] }
        if (p.name ==~ /.*/ + other) {
            inputError(["--input is a ${otherType} file (${p.name}), but --data_type is '${data_type}'.",
                        "Did you mean --data_type ${otherType}?"])
        }
        inputError(["--input file '${p.name}' is not a ${describeType(use_bam)} file."])
    }

    // Folder: look at files directly inside it
    def entries = p.listFiles() ?: []
    def ok = entries.findAll { f -> !f.isDirectory() && f.name ==~ /.*/ + want }.sort { f -> f.name }
    if (ok) { return ok }

    def nOther = entries.count { f -> !f.isDirectory() && f.name ==~ /.*/ + other }
    if (nOther) {
        inputError(["The folder ${inStr} has ${nOther} ${otherType} file(s) but no ${describeType(use_bam)} files.",
                    "Did you mean --data_type ${otherType}?"])
    }
    def subdirs = entries.findAll { f -> f.isDirectory() }.sort { d -> d.name }
    def nested  = subdirs.findAll { d -> (d.listFiles() ?: []).any { f -> f.name ==~ /.*/ + want } }
    if (nested) {
        def names = nested.take(3).collect { d -> d.name }.join(', ') + (nested.size() > 3 ? ', ...' : '')
        inputError(["The folder ${inStr} has no ${describeType(use_bam)} files directly inside it,",
                    "but its sub-folders do (${names}).",
                    'Point --input at one sub-folder, or use a quoted glob such as',
                    "--input '${inStr.replaceAll('/+$', '')}/*/*${use_bam ? '.bam' : '.fastq.gz'}'"])
    }
    inputError(["The folder ${inStr} exists but contains no ${describeType(use_bam)} files" +
                (entries ? " (it has ${entries.size()} other item(s))." : ' (it is empty).')])
}

// ---------------------------------------------------------------
// Start-up checks on the resolved input files and the output folder
// ---------------------------------------------------------------
def startupError(title, lines) {
    error((['', title] + lines.collect { l -> '  ' + l } + ['']).join('\n'))
}

def checkInputNames(input_files) {
    def spaced = input_files.findAll { f -> f.name.contains(' ') }
    if (spaced) {
        startupError('Input problem: file names must not contain spaces.',
            spaced.take(10).collect { f -> "'${f.name}'" } +
            ['', 'Rename these files (e.g. replace spaces with _) and run again.'])
    }
    def dups = input_files.groupBy { f -> sampleIdOf(f) }.findAll { _id, fs -> fs.size() > 1 }
    if (dups) {
        def lines = dups.collect { id, fs -> "${fs.collect { f -> f.name }.join(' and ')}  ->  both named '${id}'" }
        startupError('Input problem: two or more input files would get the same sample name.',
            lines + ['', 'Remove or rename the duplicates so each file has a unique name.'])
    }
}

def checkOutdir(outdir) {
    if (!outdir) {
        startupError('Output problem: --outdir is missing.',
            ['Tell the pipeline where to write results, e.g. --outdir results/'])
    }
    def od = outdir.toString()
    if (od.contains('://')) {
        return    // cloud storage: permissions are checked when results are written
    }
    def d = file(od)
    def writable = false
    try {
        d.mkdirs()
        def probe = d.resolve('.write_test')
        probe.text = 'ok'
        probe.delete()
        writable = true
    } catch (Exception _e) {
        writable = false
    }
    if (!writable) {
        startupError("Output problem: cannot write to --outdir ${od}",
            ['Check that the folder (or its parent) exists and that you have permission to write there,',
             'or choose a different --outdir.'])
    }
}

// Stop the run if any input-check report (CHECK_FASTQ / CHECK_BAM /
// SPLIT_FASTQ / SPLIT_BAM) says FAIL.  Report columns: file status reads message
def failOnBadInputs(reports) {
    def bad = []
    asList(reports).each { r ->
        r.readLines().each { line ->
            def c = line.split('\t', -1)
            if (c.size() >= 4 && c[1] == 'FAIL') {
                bad.add("${c[0]}: ${c[3]}")
            }
        }
    }
    if (bad) {
        startupError("Input problem: ${bad.size()} input file(s) cannot be used.",
            bad + ['', 'Nothing was demultiplexed. Fix or remove these files and run again.'])
    }
    return true
}

// Natural sort key, so run_9 comes before run_10 (MinKNOW chunk numbering)
def naturalKey(f) {
    return f.name.replaceAll(/\d+/) { m -> m.padLeft(12, '0') }
}

// Reject option combinations that make no sense
def checkOptions(skip_input_check, check_all_inputs) {
    if (skip_input_check && check_all_inputs) {
        startupError('Option problem: --skip_input_check and --check_all_inputs cannot be used together.',
            ['Use --check_all_inputs to read every input file before demultiplexing,',
             'or --skip_input_check to skip the check, but not both.'])
    }
}

// --data_type must be given and be bam or fastq
def checkDataType(data_type) {
    if (!data_type) {
        error "Please specify --data_type [bam|fastq]"
    }
    if (!['bam', 'fastq'].contains(data_type)) {
        error "Invalid --data_type '${data_type}'. Must be 'bam' or 'fastq'."
    }
}

// Command-line flags arrive as booleans on Nextflow 24/25 but as text
// ("false", "true") on Nextflow 26, where the text "false" would count as
// true. Turn any of these into a real true/false, or stop with a message.
def asBool(value, name) {
    if (value instanceof Boolean) {
        return value
    }
    def v = value == null ? '' : value.toString().trim().toLowerCase()
    if (['true', 'yes', '1'].contains(v)) {
        return true
    }
    if (['false', 'no', '0', ''].contains(v)) {
        return false
    }
    startupError("Option problem: --${name} must be true or false (got '${value}').",
        ["Use --${name} true or --${name} false, or leave it out."])
}
