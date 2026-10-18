#!/usr/bin/env bash
# Genetics Assistant Internship -- Week 1 pipeline
# Runs the whole analysis end to end and logs every step to results/run_log.txt
set -e
mkdir -p results
{
  echo "=== SIMG1 analysis pipeline ==="
  date
  for step in 01_simulate_dataset 02_qc_and_orfs 03_variant_calling \
              04_msa_conservation 05_local_blast 06_mutation_patterns_markers; do
    echo ""
    echo "--- running scripts/${step}.py ---"
    python3 "scripts/${step}.py"
  done
  echo ""
  echo "=== pipeline complete ==="
} 2>&1 | tee results/run_log.txt
