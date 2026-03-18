#!/bin/bash
# Run sCF analyses on all simulated alignments
set -euo pipefail

IQTREE_NEW=./build/iqtree3
IQTREE_OLD=./build/iqtree3-old
SIMDIR=tests/scfl_simulation
ALNDIR=$SIMDIR/alignments
TREEFILE=$SIMDIR/trees/true_tree.nwk
RESDIR=$SIMDIR/results
SEED=12345

mkdir -p $RESDIR

SCF_LEVELS="100 90 80 70 60 50 40 30 20 10 0"
GAP_LEVELS="0 5 10 15 20"
MODELS="GTR GTR+I GTR+I+G8 GTR+I+R4"

total=0
done=0

# Count total runs
for X in $SCF_LEVELS; do
for G in $GAP_LEVELS; do
for model in $MODELS; do
    # New binary: 4 methods
    total=$((total + 4))
    # Old binary: 1 method (scfl only)
    total=$((total + 1))
done; done; done

echo "Total runs: $total"
echo ""

run_scf() {
    local binary=$1 aln=$2 tree=$3 method=$4 model=$5 prefix=$6 extra_args=$7

    if [ -f "${prefix}.cf.stat" ]; then
        return 0
    fi

    $binary -te $tree -s $aln --${method} 1000 -m $model \
        --prefix $prefix -T 1 --seed $SEED $extra_args --quiet --redo 2>/dev/null || true
}

for X in $SCF_LEVELS; do
for G in $GAP_LEVELS; do
    aln=$ALNDIR/scf${X}_gap${G}.phy

    if [ ! -f "$aln" ]; then
        echo "MISSING: $aln"
        continue
    fi

    for model in $MODELS; do
        mname=$(echo "$model" | tr '+' '_')
        tag="scf${X}_gap${G}_${mname}"

        # New binary: scf, scfl, scflg, scfl_safe
        for method in scf scfl scflg; do
            prefix=$RESDIR/new_${tag}_${method}
            run_scf "$IQTREE_NEW" "$aln" "$TREEFILE" "$method" "$model" "$prefix" ""
            done=$((done + 1))
        done

        # scfl --safe
        prefix=$RESDIR/new_${tag}_scfl_safe
        run_scf "$IQTREE_NEW" "$aln" "$TREEFILE" "scfl" "$model" "$prefix" "--safe"
        done=$((done + 1))

        # Old binary: scfl only
        prefix=$RESDIR/old_${tag}_scfl
        run_scf "$IQTREE_OLD" "$aln" "$TREEFILE" "scfl" "$model" "$prefix" ""
        done=$((done + 1))
    done

    echo "Completed: sCF=${X}% gap=${G}% ($done/$total)"
done
done

echo ""
echo "=== All analyses complete ==="
