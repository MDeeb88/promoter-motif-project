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

        # Use only gene records
        if feature_type != "gene":
            continue

        attrs = parse_attributes(attributes)

        # Keep only protein-coding genes
        gene_biotype = attrs.get("gene_biotype", "")
        if gene_biotype != "protein_coding":
            continue

        start = int(start)
        end = int(end)

        gene_id = attrs.get("ID", "NA")
        gene_name = attrs.get("Name", gene_id)

        # GFF = 1-based inclusive
        # BED = 0-based half-open
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
