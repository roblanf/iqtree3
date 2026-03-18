#!/bin/bash
# Simulate alignments with known sCF proportions using alisim + AMAS
set -euo pipefail

IQTREE=./build/iqtree3
SIMDIR=tests/scfl_simulation
ALNDIR=$SIMDIR/alignments
TREEDIR=$SIMDIR/trees

mkdir -p $ALNDIR $TREEDIR

# Write tree files
# True tree: ((A,B),(C,D)) | ((E,F),(G,H))
cat > $TREEDIR/true_tree.nwk << 'EOF'
((A:0.1,B:0.1):0.1,(C:0.1,D:0.1):0.1,((E:0.1,F:0.1):0.1,(G:0.1,H:0.1):0.1):0.1);
EOF

# NNI alternative 1: ((A,B),(E,F)) | ((C,D),(G,H))
cat > $TREEDIR/alt1_tree.nwk << 'EOF'
((A:0.1,B:0.1):0.1,(E:0.1,F:0.1):0.1,((C:0.1,D:0.1):0.1,(G:0.1,H:0.1):0.1):0.1);
EOF

TOTAL=10000

for X in 100 90 80 70 60 50 40 30 20 10 0; do
    n_true=$((TOTAL * X / 100))
    n_alt=$((TOTAL - n_true))

    echo "=== sCF level ${X}%: n_true=$n_true, n_alt=$n_alt ==="

    # Simulate from true tree
    if [ $n_true -gt 0 ]; then
        $IQTREE --alisim $ALNDIR/true_${X} -t $TREEDIR/true_tree.nwk \
            -m JC --length $n_true --seed 42 --quiet 2>/dev/null
    fi

    # Simulate from alt tree
    if [ $n_alt -gt 0 ]; then
        $IQTREE --alisim $ALNDIR/alt_${X} -t $TREEDIR/alt1_tree.nwk \
            -m JC --length $n_alt --seed 43 --quiet 2>/dev/null
    fi

    # Concatenate with AMAS (or just copy if one side is 0)
    if [ $n_true -eq 0 ]; then
        cp $ALNDIR/alt_${X}.phy $ALNDIR/scf${X}_nogaps.phy
    elif [ $n_alt -eq 0 ]; then
        cp $ALNDIR/true_${X}.phy $ALNDIR/scf${X}_nogaps.phy
    else
        AMAS.py concat -i $ALNDIR/true_${X}.phy $ALNDIR/alt_${X}.phy \
            -f phylip -d dna -u phylip -t $ALNDIR/scf${X}_nogaps.phy 2>/dev/null
    fi

    echo "  Created: $ALNDIR/scf${X}_nogaps.phy"
done

# Inject gaps
echo ""
echo "=== Injecting gaps ==="
python3 $SIMDIR/inject_gaps.py

echo ""
echo "=== Simulation complete ==="
echo "Total alignments: 55 (11 sCF levels x 5 gap levels)"
