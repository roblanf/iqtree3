#!/bin/bash
IQTREE=./build/iqtree3
ALN=tests/Poeciliopsis8320chars.phy
TREE=tests/concat.treefile
OUTDIR=tests/model_tests
SEED=12345

mkdir -p $OUTDIR

models=(
  "TIM2+F"
  "TIM2+F+I"
  "TIM2+F+G4"
  "TIM2+F+G8"
  "TIM2+F+I+G4"
  "TIM2+F+I+G8"
  "TIM2+F+R2"
  "TIM2+F+R4"
  "TIM2+F+I+R2"
  "TIM2+F+I+R4"
)

for model in "${models[@]}"; do
  # sanitize model name for filenames
  mname=$(echo "$model" | tr '+' '_')

  for method in scf scfl scflg; do
    prefix="${OUTDIR}/${mname}_${method}"
    echo "RUN: $model $method -> $prefix"
    $IQTREE -te $TREE -s $ALN --${method} 1000 -m $model --prefix $prefix -T 1 --seed $SEED --quiet --redo 2>&1 | tail -3
  done

  # Also run scfl with --safe mode
  prefix="${OUTDIR}/${mname}_scfl_safe"
  echo "RUN: $model scfl --safe -> $prefix"
  $IQTREE -te $TREE -s $ALN --scfl 1000 -m $model --prefix $prefix -T 1 --seed $SEED --safe --quiet --redo 2>&1 | tail -3
done
echo "DONE"
