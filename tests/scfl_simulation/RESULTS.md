# Fix: `--scfl` gap detection for models with rate heterogeneity (Issue #137)

## The bug

`--scfl` gave incorrect results with any model using multiple rate categories (+G, +R, +I+G, +I+R, mixture models). The bug was in `computeSubtreeAncestralState()` in `tree/discordance.cpp`.

**Original code (line 424):**
```cpp
if (params->ancestral_site_concordance == 2 && sum > 1.0) {
    state_best = aln->STATE_UNKNOWN;
}
```

Two problems with the `sum > 1.0` gap detection threshold:

1. **`ncat_mix` scaling**: For internal nodes, partial likelihoods are summed across all `ncat_mix` rate categories. Real data gives `sum ~ ncat_mix`, not `~1.0`, so nearly everything was incorrectly filtered as "all gaps" when `ncat_mix > 1`.

2. **Partial likelihood scaling**: For deep subtrees, partial_lh values are scaled by `2^(256k)` during Felsenstein pruning, making real-data sums astronomically large. However, all-gap subtrees always have `partial_lh = 1.0` (no scaling needed), so `scale_num[ptn] = 0` for gaps and `> 0` for scaled real data.

Leaf nodes were unaffected (they use `tip_partial_lh` directly).

## The fix

```cpp
if (params->ancestral_site_concordance == 2) {
    bool is_gap;
    if (dad_branch->node->isLeaf()) {
        is_gap = sum > 1.0;
    } else {
        bool any_scaled = false;
        if (safe_numeric) {
            for (size_t c = 0; c < ncat_mix && !any_scaled; c++)
                any_scaled = dad_branch->scale_num[ptn * ncat_mix + c] != 0;
        } else {
            any_scaled = dad_branch->scale_num[ptn] != 0;
        }
        is_gap = !any_scaled && sum > (double)ncat_mix;
    }
    if (is_gap) {
        state_best = aln->STATE_UNKNOWN;
    }
}
```

**Key properties:**
- `ncat_mix = 1` and no scaling: threshold is `1.0` -- identical to original behaviour
- Scaled real-data sites: `scale_num > 0` short-circuits, never incorrectly filtered
- All-gap subtrees: `scale_num = 0` and `sum = nstates * ncat_mix > ncat_mix`, always caught
- Handles both normal and `SAFE_NUMERIC` modes (`scale_num` is per-pattern in normal mode, per-pattern-per-category in `SAFE_NUMERIC` mode)

---

## Empirical validation

**Dataset:** Poeciliopsis (8320 characters, from issue #137 reporter)
**Models:** 10 TIM2+F variants (ncat = 1, 2, 4, 5, 8, 9)
**Methods:** `--scf` (parsimony), `--scfl` (likelihood), `--scflg` (lik+gaps), `--scfl --safe`

### sCF is now stable across models (Fig 1)

Before the fix, sCFL values collapsed for multi-category models because nearly all sites were incorrectly filtered as gaps. After the fix, sCFL values are stable regardless of the number of rate categories.

![sCF by number of rate categories](../model_tests/fig1_scf_by_ncat.png)

### sN (informative sites) is now stable across models (Fig 2)

Before the fix, sN dropped to near-zero for multi-category models. After the fix, sN is consistent across all models.

![sN by number of rate categories](../model_tests/fig2_sN_by_ncat.png)

### sN is consistent across models and methods (Fig 3)

Mean sN is comparable across all 10 models and 4 methods.

![Mean sN by model and method](../model_tests/fig4_mean_sN_by_model.png)

### sCFL agrees with sCF (parsimony) across all models (Fig 4)

Pairwise scatter plots (rows = models, columns = method comparisons) show that sCFL, sCFLG, and sCFL --safe all track parsimony sCF along the diagonal.

![Method comparison scatter plots](../model_tests/fig3_method_comparison.png)

### sCFL = sCFL --safe (Fig 5)

The fix correctly handles both normal and `SAFE_NUMERIC` `scale_num` layouts. Max |diff| = 0.01 across all branches and models.

![sCFL vs sCFL --safe](../model_tests/fig5_scfl_vs_safe.png)

---

## Simulation study

### Design

- **Tree:** 8 taxa, 1 target internal branch separating {A,B,C,D} | {E,F,G,H}
- **True tree:** `((A,B),(C,D),((E,F),(G,H)));` (all branches length 0.1)
- **NNI alternative:** `((A,B),(E,F),((C,D),(G,H)));` (swaps C,D and E,F)
- **tree1%:** proportion of 10,000 sites simulated from the true tree (0, 10, 20, ..., 100%)
- **Remaining sites** simulated from the NNI alternative tree
- **Gap injection:** 0, 5, 10, 15, 20% of cells replaced with `?`
- **Models:** GTR, GTR+I, GTR+I+G8, GTR+I+R4
- **Methods:** sCF (parsimony), sCFL, sCFLG, sCFL --safe, sCFL (old/buggy binary)
- **Total:** 11 tree1% levels x 5 gap levels x 4 models x 5 methods = 1,100 runs

### Expectations tested

| Property | Expected |
|---|---|
| sN depends on gap%, NOT tree1% | Tight horizontal bands in Fig 2 |
| sCF depends on tree1%, NOT gap% | Monotonic increase in Fig 4 |
| sCFLG: gaps bias sCF downward | sCFLG lines fan apart with gap% |
| sCFL = sCFL --safe | Identical results (same seed) |
| Old sCFL broken for multi-cat models | Dramatically reduced sN |

### Results

All validation checks pass:

```
1. sCFL vs sCFL --safe max |diff|: 0.01           PASS
2. sCFL sN <= sCFLG sN: 220/220                   PASS
3. GTR+I+G8: old sN=0.7, new sN=589.9             PASS (old broken)
3. GTR+I+R4: old sN=0.7, new sN=589.7             PASS (old broken)
4. GTR (ncat=1): old sCF=41.40, new sCF=41.40      PASS (unchanged)
```

### Fig 1: sN vs Gap% (color = tree1%)

sN decreases with gap% (fewer informative sites as gaps increase). Lines are tightly bundled within each panel, confirming sN is independent of tree1%.

The old buggy sCFL (column 5) shows dramatically reduced sN for GTR+I+G8 and GTR+I+R4 (rows 3-4).

![sN vs Gap%](fig1_sN_vs_gap.png)

### Fig 2: sN vs Tree1% (color = gap%)

Flat horizontal lines confirm sN does not depend on tree topology composition. Each gap level produces a distinct band. sCFLG (column 3) shows the expected reverse pattern: sN *increases* with gaps because it treats gap sites as informative.

Old buggy sCFL (column 5, rows 3-4): sN is near zero regardless of tree1% or gap%.

![sN vs Tree1%](fig2_sN_vs_tree1.png)

### Fig 3: sCF vs Gap% (color = tree1%)

Flat horizontal bands confirm sCF is independent of gap proportion for sCF, sCFL, and sCFL --safe. Each tree1% level produces a distinct band from ~18% (tree1%=0) to ~68% (tree1%=100).

sCFLG (column 3) shows the expected deviation: sCF decreases with gap% because gap sites are treated as informative (diluting signal).

Old buggy sCFL (column 5, rows 3-4): erratic sCF values with very few informative sites.

![sCF vs Gap%](fig3_scf_vs_gap.png)

### Fig 4: sCF vs Tree1% (color = gap%)

sCF increases monotonically from ~18% (tree1%=0, all discordant data) to ~68% (tree1%=100, all concordant data). Gap levels overlap for sCF, sCFL, and sCFL --safe, confirming gap filtering works correctly.

sCFLG lines fan apart at higher gap levels. Old buggy sCFL (rows 3-4) shows erratic behavior.

![sCF vs Tree1%](fig4_scf_vs_tree1.png)
