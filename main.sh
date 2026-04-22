#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="/mnt/hpc_acegid/home/khadmig/work/data/malaria/smc-impact-under5-malaria-africa"

MALARIAGEN_DIR="/mnt/hpc_acegid/home/khadmig/work/data/malaria/malariaGendata"
SAMPLES_METADATA="${MALARIAGEN_DIR}/Pf8_samples.txt"
BED_FILE="assets/pf3d7_sp_aq_resistance_markers_v1.bed"
BED_FILE2="assets/pf3d7_sp_aq_resistance_markers_v2.bed"
HAPLOTYPE_FILE="assets/pf3d7_sp_aq_haplotypes_panel_v1.tsv"

STEP1_OUTPUT_DIR="${PROJECT_DIR}/results/genomics_smc_malaGen"
STEP1_PIPELINE_PY="scripts/pf8_smc_resistance_haplotype_pipeline_malariaGen.py"

INTERNAL_BASE_DIR="/mnt/hpc_acegid/nfsscratch/khadmig/data/malaria/demux_output_ensemble_v2"
INTERNAL_METADATA_TSV="/mnt/hpc_acegid/nfsscratch/khadmig/data/malaria/demux_input/metadata_with_sample.tsv"
FRAGMENT_MANIFEST="reference/PCR_fragment_complete.tsv"

STEP2_OUTPUT_DIR="${PROJECT_DIR}/results/genomics_internal_smc_2021"
STEP2_PIPELINE_PY="scripts/step2_internal_smc_resistance_extraction_pipeline.py"

STEP3_OUTPUT_DIR="${PROJECT_DIR}/results/genomics_internal_smc_2021_step3"
STEP3_PIPELINE_PY="scripts/step3_internal_smc_mutation_haplotype_cooccurrence_input_pipeline.py"

mkdir -p "${STEP1_OUTPUT_DIR}"
mkdir -p "${STEP2_OUTPUT_DIR}"
mkdir -p "${STEP3_OUTPUT_DIR}"

step1_run_pipeline() {
  python3 "${STEP1_PIPELINE_PY}" \
    --project-dir "${PROJECT_DIR}" \
    --malariagen-dir "${MALARIAGEN_DIR}" \
    --samples-metadata "${SAMPLES_METADATA}" \
    --bed-file "${BED_FILE}" \
    --haplotype-file "${HAPLOTYPE_FILE}" \
    --output-dir "${STEP1_OUTPUT_DIR}" \
    --year-min 2000 \
    --year-max 2026
}

step2_run_pipeline() {
  python3 "${STEP2_PIPELINE_PY}" \
    --base-dir "${INTERNAL_BASE_DIR}" \
    --metadata-tsv "${INTERNAL_METADATA_TSV}" \
    --fragment-manifest "${FRAGMENT_MANIFEST}" \
    --bed-file "${BED_FILE2}" \
    --haplotype-file "${HAPLOTYPE_FILE}" \
    --output-dir "${STEP2_OUTPUT_DIR}" \
    --min-reads 100
}

step3_run_pipeline() {
  python3 "${STEP3_PIPELINE_PY}" \
    --base-dir "${INTERNAL_BASE_DIR}" \
    --metadata-tsv "${INTERNAL_METADATA_TSV}" \
    --fragment-manifest "${FRAGMENT_MANIFEST}" \
    --bed-file "${BED_FILE2}" \
    --haplotype-file "${HAPLOTYPE_FILE}" \
    --output-dir "${STEP3_OUTPUT_DIR}" \
    --min-reads 100
}

step1_run_pipeline
step2_run_pipeline
step3_run_pipeline