#!/usr/bin/env bash
set -euo pipefail

BENCH_DIR="/home/umair/Downloads/CoTRA/benchmark"
INPUT="/home/umair/Downloads/CoTRA/data/bulkRNA/raw_gene_counts.tsv"
OUTDIR="${BENCH_DIR}/benchmark_results/bulk_real"
WARMUP="${BENCH_DIR}/benchmark_results/warmup_bulk_real"
LOGDIR="${BENCH_DIR}/benchmark_results/logs"

cd "$BENCH_DIR"

mkdir -p "$OUTDIR" "$WARMUP" "$LOGDIR"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

if [[ ! -f "$INPUT" ]]; then
  echo "ERROR: input file not found: $INPUT" >&2
  exit 1
fi

if [[ ! -f benchmark_bulk_real.R ]]; then
  echo "ERROR: benchmark_bulk_real.R not found in $BENCH_DIR" >&2
  exit 1
fi

echo "Input: $INPUT"
echo "Benchmark directory: $BENCH_DIR"
echo "Contrast: rd10 vs WT"
echo "Measured repetitions: 5"
echo

for method in DESeq2 edgeR; do
  echo "Warm-up: $method"

  /usr/bin/time -v \
    Rscript --vanilla benchmark_bulk_real.R \
      --input "$INPUT" \
      --method "$method" \
      --rep 0 \
      --seed 1234 \
      --reference WT \
      --comparison rd10 \
      --reference-regex '^WT' \
      --comparison-regex '^rd10' \
      --outdir "$WARMUP" \
    > "${LOGDIR}/bulk_real_${method}_warmup.log" \
    2> "${WARMUP}/resources_bulk_real_${method}_warmup.txt"

  for rep in 1 2 3 4 5; do
    echo "Measured run: $method repetition $rep"

    /usr/bin/time -v \
      Rscript --vanilla benchmark_bulk_real.R \
        --input "$INPUT" \
        --method "$method" \
        --rep "$rep" \
        --seed 1234 \
        --reference WT \
        --comparison rd10 \
        --reference-regex '^WT' \
        --comparison-regex '^rd10' \
        --outdir "$OUTDIR" \
      > "${LOGDIR}/bulk_real_${method}_rep${rep}.log" \
      2> "${OUTDIR}/resources_bulk_real_${method}_rep$(printf '%02d' "$rep").txt"
  done
done

echo
echo "Real-data bulk benchmark completed."
echo "Results: $OUTDIR"
