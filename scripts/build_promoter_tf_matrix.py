#!/usr/bin/env python3

import os
import sys
import glob
import pandas as pd
from collections import defaultdict

if len(sys.argv) != 5:
    print("Usage: build_promoter_tf_matrix.py promoters.fa sarus_dir threshold output.tsv")
    sys.exit(1)

promoter_fasta = sys.argv[1]
sarus_dir = sys.argv[2]
threshold = float(sys.argv[3])
out_tsv = sys.argv[4]


def clean_promoter_id(header):
    header = header.strip()
    if header.startswith(">"):
        header = header[1:]
    header = header.split()[0]
    header = header.split("::")[0]
    return header


def read_promoter_ids_from_fasta(fasta_file):
    promoter_ids = []
    with open(fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                promoter_ids.append(clean_promoter_id(line))
    return promoter_ids


promoter_ids = read_promoter_ids_from_fasta(promoter_fasta)
promoter_set = set(promoter_ids)

counts = defaultdict(lambda: defaultdict(int))

sarus_files = sorted(glob.glob(os.path.join(sarus_dir, "*.sarus.tsv")))
tf_ids = []

for sarus_file in sarus_files:
    tf_id = os.path.basename(sarus_file).replace(".sarus.tsv", "")
    tf_ids.append(tf_id)

    print(f"Parsing {tf_id}", file=sys.stderr)

    current_promoter = None

    with open(sarus_file) as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                current_promoter = clean_promoter_id(line)
                continue

            if current_promoter is None:
                continue

            parts = line.split()

            if len(parts) < 1:
                continue

            try:
                score = float(parts[0])  # score is FIRST column
            except:
                continue

            if current_promoter in promoter_set and score >= threshold:
                counts[current_promoter][tf_id] += 1


df = pd.DataFrame(0, index=promoter_ids, columns=tf_ids, dtype=int)

for promoter_id in promoter_ids:
    for tf_id, value in counts[promoter_id].items():
        df.loc[promoter_id, tf_id] = value

df.insert(0, "promoter_id", df.index)

df.to_csv(out_tsv, sep="\t", index=False)

print(f"Wrote {out_tsv}")
print(f"Shape: {df.shape}")
