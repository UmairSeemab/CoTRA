#!/usr/bin/env bash
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_ROOT="${1:-${PWD}/benchmark_results}"
REPS=5
SEED=1234
GENES=20000

BULK_OUT="${OUT_ROOT}/bulk"
WARM_OUT="${OUT_ROOT}/warmup"
LOG_OUT="${OUT_ROOT}/logs"
mkdir -p "$BULK_OUT" "$WARM_OUT" "$LOG_OUT"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

sample_sizes=(6 12 24 48 96)
methods=(DESeq2 edgeR)

for method in "${methods[@]}"; do
  for samples in "${sample_sizes[@]}"; do

    # Preliminary viability run; excluded from measured summaries.
    Rscript --vanilla "${SCRIPT_DIR}/benchmark_bulk.R" \
      --genes "$GENES" --samples "$samples" --method "$method" \
      --rep 0 --seed "$SEED" --outdir "$WARM_OUT" \
      > "${LOG_OUT}/bulk_${method}_${samples}s_warmup.log" 2>&1 || exit 1

    for rep in $(seq 1 "$REPS"); do
      rep2=$(printf "%02d" "$rep")
      id=$(printf "bulk_%s_%05dg_%03ds_rep%s" "$method" "$GENES" "$samples" "$rep2")
      echo "Measured run: $method, $samples samples, repetition $rep"

      /usr/bin/time -v -o "${BULK_OUT}/resources_${id}.txt" \
        Rscript --vanilla "${SCRIPT_DIR}/benchmark_bulk.R" \
          --genes "$GENES" --samples "$samples" --method "$method" \
          --rep "$rep" --seed "$SEED" --outdir "$BULK_OUT" \
        > "${LOG_OUT}/${id}.log" 2>&1 || exit 1
    done
  done
done

echo "Synthetic bulk benchmark completed: $BULK_OUT"
