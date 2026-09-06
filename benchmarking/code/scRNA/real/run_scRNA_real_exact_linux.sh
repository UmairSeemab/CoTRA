#!/usr/bin/env bash
set -uo pipefail

BENCH_DIR="/home/umair/Downloads/CoTRA/benchmark"
REAL_DIR="${BENCH_DIR}/real_scRNA_exact"
OUTDIR="${BENCH_DIR}/benchmark_results/scrna_real_exact"
WARMUP="${BENCH_DIR}/benchmark_results/warmup_scrna_real_exact"
LOGDIR="${BENCH_DIR}/benchmark_results/logs"

REPS=5
SEED=1234
PCS=7

cd "$BENCH_DIR" || exit 1

mkdir -p "$OUTDIR" "$WARMUP" "$LOGDIR"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

if [[ ! -f benchmark_scRNA_real.R ]]; then
  echo "ERROR: benchmark_scRNA_real.R not found in $BENCH_DIR" >&2
  exit 1
fi

sizes=(2500 5000 7500 10000 10430)

echo "CoTRA exact real retinal scRNA benchmark"
echo "Input: exact processed WT + vehicle-rd10 cells reconstructed from raw H5 counts"
echo "PCs used downstream: $PCS"
echo "Measured repetitions per size: $REPS"
echo "Single-thread constrained"
echo

for cells in "${sizes[@]}"; do

  size6=$(printf "%06d" "$cells")
  input="${REAL_DIR}/counts_real_${size6}c.rds"

  if [[ ! -f "$input" ]]; then
    echo "ERROR: missing prepared input: $input" >&2
    exit 1
  fi

  echo "Warm-up: $cells cells"

  /usr/bin/time -v \
    Rscript --vanilla benchmark_scRNA_real.R \
      --input "$input" \
      --rep 0 \
      --seed "$SEED" \
      --pcs "$PCS" \
      --outdir "$WARMUP" \
    > "${LOGDIR}/scrna_real_exact_${size6}c_warmup.log" \
    2> "${WARMUP}/resources_scrna_real_${size6}c_rep00.txt"

  if [[ $? -ne 0 ]]; then
    echo "ERROR: warm-up failed for $cells cells." >&2
    echo "See ${LOGDIR}/scrna_real_exact_${size6}c_warmup.log" >&2
    exit 1
  fi

  for rep in $(seq 1 "$REPS"); do

    rep2=$(printf "%02d" "$rep")
    id="scrna_real_${size6}c_rep${rep2}"

    echo "Measured run: $cells cells, repetition $rep"

    /usr/bin/time -v \
      Rscript --vanilla benchmark_scRNA_real.R \
        --input "$input" \
        --rep "$rep" \
        --seed "$SEED" \
        --pcs "$PCS" \
        --outdir "$OUTDIR" \
      > "${LOGDIR}/${id}_exact.log" \
      2> "${OUTDIR}/resources_${id}.txt"

    if [[ $? -ne 0 ]]; then
      echo "ERROR: measured run failed: $id" >&2
      echo "See ${LOGDIR}/${id}_exact.log" >&2
      exit 1
    fi
  done
done

echo
echo "Exact real retinal scRNA benchmark completed."
echo "Results: $OUTDIR"
