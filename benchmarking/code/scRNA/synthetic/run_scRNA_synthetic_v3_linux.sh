#!/usr/bin/env bash
set -uo pipefail

BENCH_DIR="/home/umair/Downloads/CoTRA/benchmark"
SYN_DIR="${BENCH_DIR}/synthetic_scRNA"
OUTDIR="${BENCH_DIR}/benchmark_results/scrna_synthetic_v3"
WARMUP="${BENCH_DIR}/benchmark_results/warmup_scrna_synthetic_v3"
LOGDIR="${BENCH_DIR}/benchmark_results/logs"

REPS=5
SEED=1234
TRUTH_CLUSTERS=12

cd "$BENCH_DIR" || exit 1
mkdir -p "$OUTDIR" "$WARMUP" "$LOGDIR"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

sizes=(2500 5000 10000 25000 50000)

for cells in "${sizes[@]}"; do
  size6=$(printf "%06d" "$cells")
  input="${SYN_DIR}/counts_synthetic_${size6}c.rds"

  [[ -f "$input" ]] || { echo "Missing $input" >&2; exit 1; }

  echo "Warm-up: $cells cells"
  /usr/bin/time -v \
    Rscript --vanilla benchmark_scRNA_synthetic_v3.R \
      --input "$input" --rep 0 --seed "$SEED" \
      --species mouse --truth-clusters "$TRUTH_CLUSTERS" \
      --outdir "$WARMUP" \
    > "${LOGDIR}/scrna_v3_${size6}c_warmup.log" \
    2> "${WARMUP}/resources_scrna_synthetic_${size6}c_rep00.txt"

  [[ $? -eq 0 ]] || { echo "Warm-up failed for $cells cells" >&2; exit 1; }

  for rep in $(seq 1 "$REPS"); do
    rep2=$(printf "%02d" "$rep")
    id="scrna_synthetic_${size6}c_rep${rep2}"
    echo "Measured run: $cells cells, repetition $rep"

    /usr/bin/time -v \
      Rscript --vanilla benchmark_scRNA_synthetic_v3.R \
        --input "$input" --rep "$rep" --seed "$SEED" \
        --species mouse --truth-clusters "$TRUTH_CLUSTERS" \
        --outdir "$OUTDIR" \
      > "${LOGDIR}/${id}_v3.log" \
      2> "${OUTDIR}/resources_${id}.txt"

    [[ $? -eq 0 ]] || { echo "Run failed: $id" >&2; exit 1; }
  done
done

echo "Synthetic scRNA benchmark v3 completed."
echo "Results: $OUTDIR"
