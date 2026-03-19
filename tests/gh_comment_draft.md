Thanks for reporting this @marianamateos - you've found a real bug.

## The problem

`--scfl` gives incorrect results when used with any model that has multiple rate categories (i.e. `+G`, `+R`, or combinations like `+I+G`, `+I+R`). This includes a lot of models that people would normally use.

The bug is in the site filtering logic. `--scfl` improves on `--scflg` by trying to detect and skip sites where a subtree has no data (all gaps) for a particular site. It does this by checking if the sum of partial likelihoods exceeds 1.0 — which works correctly for leaf nodes (a gap gives sum=4.0 for DNA, a real base gives sum=1.0). But for internal nodes, the partial likelihoods are summed across rate categories, so with e.g. `+G4` the sum is ~4.0 even for perfectly clean data. The filter then incorrectly marks nearly every site at every internal node as "all gaps" and discards them.

The result: for deep branches (where all four subtrees in the quartet are internal), almost no sites survive the filter, and sCFL values become meaningless.

## Evidence

We tested your dataset with 10 models varying from 1 to 9 rate categories, running all three methods (`--scf`, `--scfl`, `--scflg`). The bar chart below shows the mean number of informative sites (sN) used by each method:

[figure]

- **Blue (sCF, parsimony)** and **green (sCFLG)** are stable across all models (~130-140 informative sites per branch on average).
- **Red (sCFL)** drops to ~30 informative sites for any model with gamma or FreeRate categories. For individual branches it can drop as low as 0.2 sites.

The only models where sCFL works correctly are `TIM2+F` (no rate heterogeneity) and `TIM2+F+I` (invariant sites only, which are handled differently in the likelihood kernel).

## What to do now

For now, please use `--scflg` instead of `--scfl`. Despite what the documentation says (and what I said above!!), `--scflg` gives correct results while `--scfl` does not for models with rate heterogeneity.

We're working on a fix for `--scfl`.