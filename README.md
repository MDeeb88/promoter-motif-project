# Promoter × TF Motif-Count Matrix Preparation

## Goal

The goal of this task is to prepare a matrix where:

- rows = promoters of protein-coding genes
- columns = transcription factor motifs
- values = number of motif hits in each promoter

The final output should be two TSV tables:

```text
matrices/promoter_TF_counts_threshold_4.tsv
matrices/promoter_TF_counts_threshold_5.tsv
```

Each value in the matrix is the number of motif hits for a given TF in a given promoter, counted at a selected SARUS score threshold.

The promoter region is defined as:

```text
-500 bp to +500 bp around the TSS
```

Important: use **genomic FASTA**, not RNA FASTA. Promoters are genomic DNA regions.

---

## Input data

We will use the T2T-CHM13v2.0 human genome.

Genome FASTA:

```text
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
```

Genome annotation:

```text
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz
```

HOCOMOCO PWM archive:

```text
https://hocomoco14.autosome.org/final_bundle/hocomoco13/H13CORE/H13CORE_pwm.tar.gz
```

SARUS:

```text
https://github.com/autosome-ru/sarus
```

---

## Expected final files

The final files to return are:

```text
promoters/protein_coding_promoters_500bp.bed
promoters/protein_coding_promoters_500bp.fa
matrices/promoter_TF_counts_threshold_4.tsv
matrices/promoter_TF_counts_threshold_5.tsv
```

Also return a short report with:

```text
Number of promoters:
Number of representative TFs requested:
Number of SARUS files generated:
Number of TFs in final matrix:
Shape of threshold 4 matrix:
Shape of threshold 5 matrix:
Any missing TF motif files:
Any errors:
```

---

## 0. Create project folder

```bash
mkdir -p promoter_motif_project
cd promoter_motif_project

mkdir -p data genome annotation promoters motifs sarus_results matrices scripts
```

Expected folder structure:

```text
promoter_motif_project/
├── data/
├── genome/
├── annotation/
├── promoters/
├── motifs/
├── sarus_results/
├── matrices/
└── scripts/
```

---

## 1. Download genome and annotation

Download the genomic FASTA:

```bash
cd genome

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz

gunzip GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
```

Download the genome annotation:

```bash
cd ../annotation

wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz

gunzip GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz
```

After this, the following files should exist:

```text
genome/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna
annotation/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff
```

Check them:

```bash
cd ..
ls -lh genome/
ls -lh annotation/
```

---

## 2. Make promoter BED file

We need one promoter per protein-coding gene.

For genes on the `+` strand:

```text
TSS = gene start
promoter = TSS - 500 bp to TSS + 500 bp
```

For genes on the `-` strand:

```text
TSS = gene end
promoter = TSS - 500 bp to TSS + 500 bp
```

The GFF file uses 1-based coordinates. BED uses 0-based coordinates. The script below handles this conversion.

Create the script:

```bash
nano scripts/make_promoter_bed.py
```

Paste this code:

```python
#!/usr/bin/env python3

import sys

if len(sys.argv) != 4:
    print("Usage: make_promoter_bed.py annotation.gff output.bed window")
    sys.exit(1)

gff_file = sys.argv[1]
out_bed = sys.argv[2]
window = int(sys.argv[3])


def parse_attributes(attr_string):
    """
    Parse GFF3 attributes column.

    Example:
    ID=gene-ABC;Dbxref=GeneID:123;Name=ABC;gene_biotype=protein_coding
    """
    attrs = {}
    for item in attr_string.strip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = value
    return attrs


with open(gff_file) as fin, open(out_bed, "w") as fout:
    for line in fin:
        if line.startswith("#"):
            continue

        parts = line.rstrip("\n").split("\t")
        if len(parts) < 9:
            continue

        seqid, source, feature_type, start, end, score, strand, phase, attributes = parts

        # Use only gene records.
        if feature_type != "gene":
            continue

        attrs = parse_attributes(attributes)

        # Keep only protein-coding genes.
        # In NCBI GFF files, protein-coding genes usually have gene_biotype=protein_coding.
        gene_biotype = attrs.get("gene_biotype", "")
        if gene_biotype != "protein_coding":
            continue

        start = int(start)
        end = int(end)

        gene_id = attrs.get("ID", "NA")
        gene_name = attrs.get("Name", gene_id)

        # GFF coordinates are 1-based inclusive.
        # BED coordinates are 0-based half-open.
        if strand == "+":
            tss_1based = start
        elif strand == "-":
            tss_1based = end
        else:
            continue

        bed_start = max(0, tss_1based - window - 1)
        bed_end = tss_1based + window

        promoter_id = f"{gene_id}|{gene_name}|{seqid}:{bed_start}-{bed_end}|{strand}"

        # BED6 format:
        # chrom, start, end, name, score, strand
        fout.write(
            f"{seqid}\t{bed_start}\t{bed_end}\t{promoter_id}\t0\t{strand}\n"
        )
```

Run the script:

```bash
python3 scripts/make_promoter_bed.py \
  annotation/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff \
  promoters/protein_coding_promoters_500bp.bed \
  500
```

Check the result:

```bash
head promoters/protein_coding_promoters_500bp.bed
wc -l promoters/protein_coding_promoters_500bp.bed
```

Expected result: approximately around 20,000 protein-coding gene promoters, depending on the exact annotation version.

---

## 3. Extract promoter FASTA sequences

Install `bedtools` if needed.

With conda:

```bash
conda install -c bioconda bedtools
```

Extract promoter sequences:

```bash
bedtools getfasta \
  -fi genome/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna \
  -bed promoters/protein_coding_promoters_500bp.bed \
  -fo promoters/protein_coding_promoters_500bp.fa \
  -name \
  -s
```

Explanation of important flags:

```text
-name = use promoter name from the BED file as the FASTA header
-s    = respect strand information
```

Check the FASTA file:

```bash
grep "^>" promoters/protein_coding_promoters_500bp.fa | head
grep "^>" promoters/protein_coding_promoters_500bp.fa | wc -l
```

The number of FASTA records should match the number of BED records:

```bash
wc -l promoters/protein_coding_promoters_500bp.bed
grep "^>" promoters/protein_coding_promoters_500bp.fa | wc -l
```

---

## 4. Download and test SARUS

Download `sarus.jar` from:

```text
https://github.com/autosome-ru/sarus
```

Put it into the project folder:

```text
promoter_motif_project/sarus.jar
```

Check Java:

```bash
java -version
```

SARUS requires Java. If Java is not installed, install it first.

Test SARUS:

```bash
java -cp sarus.jar ru.autosome.SARUS
```

If SARUS prints help or usage information, it works.

The basic SARUS command has this form:

```bash
java -cp sarus.jar ru.autosome.SARUS <sequences.multifasta> <weight.matrix> <threshold>
```

In this project, we will run SARUS once per TF motif.

---

## 5. Download HOCOMOCO PWM files

Download the PWM archive:

```bash
cd motifs

wget https://hocomoco14.autosome.org/final_bundle/hocomoco13/H13CORE/H13CORE_pwm.tar.gz

tar -xzf H13CORE_pwm.tar.gz

cd ..
```

Check the files:

```bash
find motifs -type f | head
find motifs -type f | wc -l
```

---

## 6. Prepare representative TF list

We need to use only representative TFs from the provided TSV file.

Place the file here:

```text
data/representative_TFs.tsv
```

Check its structure:

```bash
head data/representative_TFs.tsv
```

We need a simple text file with one PWM ID or motif name per line:

```text
data/representative_pwm_ids.txt
```

Example content:

```text
AHR_HUMAN.H13CORE.0.P.B
ARNT_HUMAN.H13CORE.0.P.B
MYC_HUMAN.H13CORE.0.P.B
```

If the first column of the TSV contains PWM IDs, create the file like this:

```bash
cut -f1 data/representative_TFs.tsv | tail -n +2 > data/representative_pwm_ids.txt
```

Check it:

```bash
head data/representative_pwm_ids.txt
wc -l data/representative_pwm_ids.txt
```

Expected: approximately 600 TFs.

Important: if the names in `representative_pwm_ids.txt` do not match the PWM file names in the `motifs/` folder, stop and report this problem.

---

## 7. Run SARUS for all representative TFs

We will run SARUS once per TF.

Input FASTA:

```text
promoters/protein_coding_promoters_500bp.fa
```

Output folder:

```text
sarus_results/
```

We first run SARUS with threshold `1`, then later count hits above thresholds `4` and `5` from the SARUS output.

Create the script:

```bash
nano scripts/run_sarus_all_tfs.sh
```

Paste this code:

```bash
#!/usr/bin/env bash

set -euo pipefail

PROMOTER_FASTA="promoters/protein_coding_promoters_500bp.fa"
PWM_LIST="data/representative_pwm_ids.txt"
PWM_DIR="motifs"
OUT_DIR="sarus_results"
SARUS_JAR="sarus.jar"

mkdir -p "${OUT_DIR}"

while read -r PWM_ID
do
    # Skip empty lines.
    if [[ -z "${PWM_ID}" ]]; then
        continue
    fi

    echo "Processing ${PWM_ID}"

    # Find PWM file.
    # This assumes PWM_ID is part of the PWM file name.
    PWM_FILE=$(find "${PWM_DIR}" -type f -name "*${PWM_ID}*" | head -n 1)

    if [[ -z "${PWM_FILE}" ]]; then
        echo "WARNING: PWM file not found for ${PWM_ID}" >&2
        continue
    fi

    OUT_FILE="${OUT_DIR}/${PWM_ID}.sarus.tsv"

    java -Xmx4G -cp "${SARUS_JAR}" ru.autosome.SARUS \
        "${PROMOTER_FASTA}" \
        "${PWM_FILE}" \
        1 \
        > "${OUT_FILE}"

done < "${PWM_LIST}"
```

Make it executable:

```bash
chmod +x scripts/run_sarus_all_tfs.sh
```

Run it:

```bash
bash scripts/run_sarus_all_tfs.sh
```

Check outputs:

```bash
ls sarus_results | head
ls sarus_results | wc -l
```

Expected: approximately one SARUS output file per representative TF.

---

## 8. Inspect SARUS output format

Before building the matrix, inspect several SARUS output files:

```bash
head sarus_results/*.sarus.tsv | head -n 40
```

You must identify which columns contain:

```text
promoter or sequence name
motif score
position
strand
```

For the final matrix, we only need:

```text
promoter name
score
```

Because the final value is simply the number of hits above a threshold.

If SARUS output contains header lines or comment lines, the parsing script will skip them.

---

## 9. Build promoter × TF matrices

Create the script:

```bash
nano scripts/build_promoter_tf_matrix.py
```

Paste this code:

```python
#!/usr/bin/env python3

import os
import sys
import glob
import pandas as pd
from collections import defaultdict

if len(sys.argv) != 5:
    print(
        "Usage: build_promoter_tf_matrix.py promoters.fa sarus_dir threshold output.tsv"
    )
    sys.exit(1)

promoter_fasta = sys.argv[1]
sarus_dir = sys.argv[2]
threshold = float(sys.argv[3])
out_tsv = sys.argv[4]

# Adjust these after inspecting SARUS output.
# Python uses 0-based column numbers.
#
# Example:
# if SARUS output is:
# sequence_name    position    strand    score
# then:
# PROMOTER_COL = 0
# SCORE_COL = 3

PROMOTER_COL = 0
SCORE_COL = -1


def read_promoter_ids_from_fasta(fasta_file):
    promoter_ids = []
    with open(fasta_file) as f:
        for line in f:
            if line.startswith(">"):
                promoter_id = line[1:].strip().split()[0]
                promoter_ids.append(promoter_id)
    return promoter_ids


promoter_ids = read_promoter_ids_from_fasta(promoter_fasta)

# counts[promoter_id][tf_id] = number of motif hits
counts = defaultdict(lambda: defaultdict(int))

sarus_files = sorted(glob.glob(os.path.join(sarus_dir, "*.sarus.tsv")))

tf_ids = []

for sarus_file in sarus_files:
    tf_id = os.path.basename(sarus_file).replace(".sarus.tsv", "")
    tf_ids.append(tf_id)

    print(f"Parsing {tf_id}", file=sys.stderr)

    with open(sarus_file) as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            # Skip comments or possible header lines.
            if line.startswith("#"):
                continue

            parts = line.split()

            # Skip malformed lines.
            if len(parts) < 2:
                continue

            try:
                promoter_id = parts[PROMOTER_COL]
                score = float(parts[SCORE_COL])
            except Exception:
                # This skips header lines or unexpected rows.
                continue

            if score >= threshold:
                counts[promoter_id][tf_id] += 1


# Build matrix.
df = pd.DataFrame(index=promoter_ids, columns=tf_ids)
df = df.fillna(0)

for promoter_id in promoter_ids:
    for tf_id in tf_ids:
        df.loc[promoter_id, tf_id] = counts[promoter_id].get(tf_id, 0)

df.insert(0, "promoter_id", df.index)

df.to_csv(out_tsv, sep="\t", index=False)
```

Run for threshold 4:

```bash
python3 scripts/build_promoter_tf_matrix.py \
  promoters/protein_coding_promoters_500bp.fa \
  sarus_results \
  4 \
  matrices/promoter_TF_counts_threshold_4.tsv
```

Run for threshold 5:

```bash
python3 scripts/build_promoter_tf_matrix.py \
  promoters/protein_coding_promoters_500bp.fa \
  sarus_results \
  5 \
  matrices/promoter_TF_counts_threshold_5.tsv
```

Check the matrices:

```bash
head matrices/promoter_TF_counts_threshold_4.tsv
head matrices/promoter_TF_counts_threshold_5.tsv
```

Check matrix dimensions:

```bash
python3 - << 'EOF'
import pandas as pd

for f in [
    "matrices/promoter_TF_counts_threshold_4.tsv",
    "matrices/promoter_TF_counts_threshold_5.tsv"
]:
    df = pd.read_csv(f, sep="\t")
    print(f, df.shape)
EOF
```

Expected:

```text
number of rows ≈ number of protein-coding promoters
number of columns ≈ number of TFs + 1 promoter_id column
```

---

## 10. Quality control

Please run all checks below and include the results in the final report.

### A. Number of promoters

```bash
wc -l promoters/protein_coding_promoters_500bp.bed
grep "^>" promoters/protein_coding_promoters_500bp.fa | wc -l
```

These two numbers should match.

---

### B. Number of TF motif files processed

```bash
wc -l data/representative_pwm_ids.txt
ls sarus_results/*.sarus.tsv | wc -l
```

These numbers should be close or identical.

If not identical, report which TFs were missing.

---

### C. Matrix shape

```bash
python3 - << 'EOF'
import pandas as pd

df4 = pd.read_csv("matrices/promoter_TF_counts_threshold_4.tsv", sep="\t")
df5 = pd.read_csv("matrices/promoter_TF_counts_threshold_5.tsv", sep="\t")

print("Threshold 4:", df4.shape)
print("Threshold 5:", df5.shape)
EOF
```

---

### D. Basic motif count summary

```bash
python3 - << 'EOF'
import pandas as pd

for f in [
    "matrices/promoter_TF_counts_threshold_4.tsv",
    "matrices/promoter_TF_counts_threshold_5.tsv"
]:
    df = pd.read_csv(f, sep="\t")
    X = df.drop(columns=["promoter_id"])

    print()
    print(f)
    print("Total motif hits:", X.values.sum())
    print("Mean hits per promoter:", X.sum(axis=1).mean())
    print("Mean hits per TF:", X.sum(axis=0).mean())
    print("Number of all-zero promoters:", (X.sum(axis=1) == 0).sum())
EOF
```

Threshold 5 should usually have fewer hits than threshold 4.

---

## 11. Troubleshooting

### Problem: `bedtools getfasta` gives empty output

Check that chromosome names in the BED file and FASTA file match.

Run:

```bash
head promoters/protein_coding_promoters_500bp.bed
grep ">" genome/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna | head
```

The sequence IDs must match.

---

### Problem: SARUS does not run

Check Java:

```bash
java -version
```

Check that `sarus.jar` is in the project folder:

```bash
ls -lh sarus.jar
```

Try:

```bash
java -cp sarus.jar ru.autosome.SARUS
```

---

### Problem: many PWM files are missing

Check whether IDs in `data/representative_pwm_ids.txt` match the actual motif file names:

```bash
head data/representative_pwm_ids.txt
find motifs -type f | head
```

If the names do not match, do not continue silently. Report the mismatch.

---

### Problem: final matrix contains only zeros

Possible reasons:

1. SARUS output was parsed incorrectly.
2. The wrong score column was used.
3. The threshold is too strict.
4. SARUS output files are empty.

First inspect SARUS output:

```bash
head sarus_results/*.sarus.tsv | head -n 40
```

Then adjust these variables in `scripts/build_promoter_tf_matrix.py`:

```python
PROMOTER_COL = 0
SCORE_COL = -1
```

---

## 12. Important notes

### Use genomic FASTA, not RNA FASTA

Use this file:

```text
GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz
```

Do not use this file for promoter extraction:

```text
GCF_009914755.1_T2T-CHM13v2.0_rna.fna.gz
```

RNA FASTA contains transcript sequences, not genomic promoter regions.

---

### Strand matters

For `+` strand genes, the TSS is at the gene start.

For `-` strand genes, the TSS is at the gene end.

The BED script handles this.

---

### Threshold logic

SARUS is run with threshold `1` to collect many candidate hits.

Then we build two final matrices by filtering hits:

```text
score >= 4
score >= 5
```

This avoids rerunning SARUS separately for every final threshold.

---

### Do not add RNA expression values yet

This task is only for preparing the promoter × TF motif-count matrix.

Expression values will be added later.

Final student matrix should look like this:

```text
promoter_id    TF1    TF2    TF3    ...
```

Later we will merge it with RNA data:

```text
promoter_id    gene_id    TF1    TF2    ...    transcription
```

---

## Final checklist

Before sending the results, confirm that you have:

- downloaded genomic FASTA
- downloaded GFF annotation
- created promoter BED file
- extracted promoter FASTA file
- downloaded HOCOMOCO PWM files
- prepared representative TF list
- run SARUS for all representative TFs
- created threshold 4 matrix
- created threshold 5 matrix
- performed quality control checks
- prepared a short report

