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
    [[ -z "${PWM_ID}" ]] && continue

    PWM_FILE=$(find "${PWM_DIR}" -type f -name "${PWM_ID}.pwm" | head -n 1)

    if [[ -z "${PWM_FILE}" ]]; then
        echo "WARNING: PWM file not found for ${PWM_ID}" >&2
        continue
    fi

    echo "Processing ${PWM_ID}"

    java -Xmx4G -cp "${SARUS_JAR}" ru.autosome.SARUS \
        "${PROMOTER_FASTA}" \
        "${PWM_FILE}" \
        1 \
        > "${OUT_DIR}/${PWM_ID}.sarus.tsv"

done < "${PWM_LIST}"
