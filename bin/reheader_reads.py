#!/usr/bin/env python3
"""
Restore ONT header tags to a demultiplexed FASTQ using the UUID lookup
TSV(s) written by extract_fastq_tags.py.

Memory: only the tags for reads that are in THIS FASTQ are kept.
  Pass 1 reads the FASTQ and collects its read IDs.
  Pass 2 loads the TSV(s), keeping only lines whose read ID is in that set.
  Pass 3 writes the reheadered FASTQ.
So memory grows with the reads in this barcode, not with the whole run.
The command line and the output are unchanged.
"""
import gzip
import sys
import argparse

def parse_args():
    p = argparse.ArgumentParser(
        description="Restore ONT header tags to demuxed FASTQ using UUID lookup TSV")
    p.add_argument("fastq",          help="Input demuxed FASTQ (gzipped)")
    p.add_argument("output",         help="Output reheadered FASTQ (gzipped)")
    p.add_argument("tsv", nargs="+", help="One or more uuid_tags.tsv files")
    return p.parse_args()

def read_id(header):
    """Read ID from a FASTQ header line, parsed exactly as in reheader():
    text after '@' up to the first whitespace (space OR tab - dorado 2.x
    separates header fields with tabs)."""
    parts = header[1:].rstrip("\n").split(None, 1)
    return parts[0] if parts else ""

def collect_ids(fastq_in):
    """Pass 1: the set of read IDs present in this FASTQ."""
    needed = set()
    with gzip.open(fastq_in, "rt") as inp:
        while True:
            header = inp.readline()
            if not header:
                break           # EOF
            inp.readline(); inp.readline(); inp.readline()   # seq, +, qual
            if header.startswith("@"):
                uuid = read_id(header)
                if uuid:
                    needed.add(uuid)
    print(f"[reheader_reads] {len(needed)} reads in {fastq_in}", file=sys.stderr)
    return needed

def build_lookup(tsv_files, needed):
    """Pass 2: load tags only for read IDs in `needed`."""
    lookup  = {}
    scanned = 0
    for f in tsv_files:
        with open(f) as fh:
            for line in fh:
                line = line.rstrip()
                if not line:
                    continue
                parts = line.split("\t", 1)
                if len(parts) == 2:
                    scanned += 1
                    if parts[0] in needed:
                        lookup[parts[0]] = parts[1]
    print(f"[reheader_reads] Loaded tags for {len(lookup)} of {len(needed)} needed reads "
          f"(scanned {scanned} entries in {len(tsv_files)} TSV(s))", file=sys.stderr)
    return lookup

def reheader(fastq_in, fastq_out, lookup):
    """Pass 3: write the FASTQ with the original header tags restored."""
    missing = 0
    total   = 0
    with gzip.open(fastq_in, "rt") as inp, \
         gzip.open(fastq_out, "wt") as out:
        while True:
            # Read exactly 4 lines per FASTQ record — avoids misidentifying
            # quality score lines that start with '@' as header lines
            header   = inp.readline()
            if not header:
                break           # EOF
            seq      = inp.readline()
            plus     = inp.readline()
            qual     = inp.readline()

            if not header.startswith("@"):
                print(f"[reheader_reads] WARNING: unexpected header line: {header.rstrip()}",
                      file=sys.stderr)
                continue

            uuid = read_id(header)

            if not uuid:
                # Empty header — write as-is and continue
                out.write(header)
                out.write(seq)
                out.write(plus)
                out.write(qual)
                continue

            total += 1
            extra  = lookup.get(uuid)

            if extra:
                out.write(f"@{uuid}{extra}\n")
            else:
                missing += 1
                out.write(header)

            out.write(seq)
            out.write(plus)
            out.write(qual)

    if missing:
        print(f"[reheader_reads] WARNING: {missing}/{total} reads had no tag match",
              file=sys.stderr)
    else:
        print(f"[reheader_reads] {total} reads reheadered -> {fastq_out}",
              file=sys.stderr)

def main():
    args   = parse_args()
    needed = collect_ids(args.fastq)
    lookup = build_lookup(args.tsv, needed)
    del needed                      # not needed for writing; free it early
    reheader(args.fastq, args.output, lookup)

if __name__ == "__main__":
    main()
