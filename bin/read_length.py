#!/usr/bin/env python3
"""
Read-length plots for one barcode, from a per-length count table.

Input : <length_counts.tsv>  two tab-separated columns: read_length  count
        (written by a single awk pass in modules/read_length.nf, so the
        reads never have to be listed or sorted one by one)
Output: <output.png>                     read-length histogram (100 bins)
        <output>_weighted.png            weighted contribution (bin centre x count)
        <output>.bin_table.txt           bin_start  bin_end  count  weighted_sum

The bins, counts and plots are the same as the previous version, which
loaded every read length individually.
"""
import sys
import matplotlib
matplotlib.use("Agg")  # headless mode for Nextflow
import matplotlib.pyplot as plt
import numpy as np

if len(sys.argv) != 3:
    print("Usage: read_length.py <length_counts.tsv> <output.png>")
    sys.exit(1)

counts_file = sys.argv[1]
out_png = sys.argv[2]
out_table = out_png.replace(".png", ".bin_table.txt")
out_weighted_png = out_png.replace(".png", "_weighted.png")

# Load the per-length counts (one line per distinct read length)
lengths, counts = [], []
with open(counts_file) as f:
    for line in f:
        parts = line.split()
        if len(parts) >= 2:
            lengths.append(int(parts[0]))
            counts.append(int(parts[1]))

if not lengths:
    print(f"read_length.py: no reads in {counts_file} - no plots made", file=sys.stderr)
    sys.exit(1)

lengths = np.array(lengths)
counts = np.array(counts)

# Define bins: 100 equal bins from shortest to longest read (as before).
# If every read has the same length, use a 1 bp-wide range instead of
# 101 identical edges, which some numpy/matplotlib versions reject.
lo, hi = lengths.min(), lengths.max()
if lo == hi:
    bins = np.linspace(lo - 0.5, hi + 0.5, 101)
else:
    bins = np.linspace(lo, hi, 101)

# Original histogram (for counts) - weights=counts gives the same bars as
# passing every read length individually
plt.figure(figsize=(10,5))
n, bin_edges, _ = plt.hist(lengths, bins=bins, weights=counts, color='forestgreen', alpha=0.7)
plt.title("Read Length Distribution")
plt.xlabel("Read Length (bp)")
plt.ylabel("Frequency")
plt.grid(axis='y', linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig(out_png, dpi=300)

# Prepare data for weighted bar plot
bin_centers = [(bin_edges[i] + bin_edges[i+1]) / 2 for i in range(len(n))]
weighted_sums = [((bin_edges[i] + bin_edges[i+1]) / 2) * count for i, count in enumerate(n)]

# Save table with weighted_sum column
with open(out_table, "w") as out:
    out.write("bin_start\tbin_end\tcount\tweighted_sum\n")
    for i in range(len(n)):
        out.write(f"{int(bin_edges[i])}\t{int(bin_edges[i+1])}\t{int(n[i])}\t{weighted_sums[i]:.2f}\n")

# Create bar plot using weighted_sum
plt.figure(figsize=(12,5))
plt.bar(bin_centers, weighted_sums, width=(bin_edges[1]-bin_edges[0])*0.9, color='steelblue', alpha=0.8)
plt.title("Weighted Read Length Contribution")
plt.xlabel("Average Read Length (bp)")
plt.ylabel("Weighted Sum (bp * count)")
plt.grid(axis='y', linestyle='--', alpha=0.5)
plt.tight_layout()
plt.savefig(out_weighted_png, dpi=300)
