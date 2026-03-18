#!/usr/bin/env python3
"""Analyze sCFL simulation results and generate figures."""

import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import numpy as np

SIMDIR = "tests/scfl_simulation"
RESDIR = os.path.join(SIMDIR, "results")

SCF_LEVELS = [100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 0]
GAP_LEVELS = [0, 5, 10, 15, 20]
MODELS = ["GTR", "GTR+I", "GTR+I+G8", "GTR+I+R4"]
NEW_METHODS = ["scf", "scfl", "scflg", "scfl_safe"]
METHOD_LABELS = {
    "scf": "sCF (parsimony)",
    "scfl": "sCFL (likelihood)",
    "scflg": "sCFLG (lik+gaps)",
    "scfl_safe": "sCFL (--safe)",
    "old_scfl": "sCFL (old/buggy)",
}

# ============================================================
# Step 1: Identify target branch
# ============================================================
# The target branch separates {A,B,C,D} from {E,F,G,H}.
# Run one analysis and find the branch ID from the cf.tree file.
# With fixed topology, the branch IDs are consistent across all runs.

def get_clade_taxa_and_label(newick, pos):
    """Starting at '(' in newick at position pos, return (taxa_set, label_str, end_pos).
    label_str is the node label after the closing ')', typically the sCF value."""
    assert newick[pos] == '('
    taxa = set()
    i = pos + 1
    while i < len(newick):
        if newick[i] == '(':
            sub_taxa, _, i = get_clade_taxa_and_label(newick, i)
            taxa |= sub_taxa
        elif newick[i] in '),;:':
            pass  # handled below
        else:
            # Taxon name: read until : , ) or ;
            j = i
            while j < len(newick) and newick[j] not in ':,);':
                j += 1
            name = newick[i:j].strip()
            if name:
                taxa.add(name)
            i = j
            continue
        if newick[i] == ':':
            # Skip branch length
            i += 1
            while i < len(newick) and newick[i] not in ',);':
                i += 1
            continue
        if newick[i] == ',':
            i += 1
            continue
        if newick[i] == ')':
            # Read label after closing paren
            i += 1
            label_start = i
            while i < len(newick) and newick[i] not in ':,);':
                i += 1
            label = newick[label_start:i].strip()
            # Skip branch length if present
            if i < len(newick) and newick[i] == ':':
                i += 1
                while i < len(newick) and newick[i] not in ',);':
                    i += 1
            return taxa, label, i
    return taxa, "", i


def find_target_branch_in_file(cf_tree_file, cf_stat_file,
                                target_taxa=frozenset({"E", "F", "G", "H"})):
    """Parse a cf.tree to find the target bipartition and return matching cf.stat data.

    Branch IDs can change between runs (IQ-TREE renumbers during optimization),
    so we identify the target by its bipartition in the tree, match by sCF label
    to the cf.stat entries, and return the target branch's data dict.
    """
    if not os.path.exists(cf_tree_file) or not os.path.exists(cf_stat_file):
        return None

    with open(cf_tree_file) as f:
        tree_str = f.read().strip().rstrip(';')

    all_taxa = set(re.findall(r'([A-H]):', tree_str))

    # Find all internal clades and their sCF labels
    target_label = None
    for i, ch in enumerate(tree_str):
        if ch == '(':
            taxa, label, _ = get_clade_taxa_and_label(tree_str, i)
            if taxa == target_taxa or taxa == all_taxa - target_taxa:
                target_label = label
                break

    if target_label is None:
        return None

    # Match the target's sCF label to the cf.stat entry
    target_scf = float(target_label)
    branches = parse_stat_file(cf_stat_file)

    # Find branch with sCF closest to the tree label
    best = None
    best_diff = float('inf')
    for b in branches:
        diff = abs(b["sCF"] - target_scf)
        if diff < best_diff:
            best_diff = diff
            best = b

    return best


def parse_stat_file(stat_file):
    """Parse a .cf.stat file and return list of dicts with branch data."""
    rows = []
    if not os.path.exists(stat_file):
        return rows
    with open(stat_file) as f:
        for line in f:
            if line.startswith("#") or line.startswith("ID"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 10:
                continue
            rows.append({
                "branch_id": int(parts[0]),
                "sCF": float(parts[1]),
                "sCF_N": float(parts[2]),
                "sDF1": float(parts[3]),
                "sDF1_N": float(parts[4]),
                "sDF2": float(parts[5]),
                "sDF2_N": float(parts[6]),
                "sN": float(parts[7]),
                "label": parts[8],
                "length": float(parts[9]),
            })
    return rows


# ============================================================
# Step 2: Target branch identification strategy
# ============================================================
# The target branch separates {A,B,C,D} from {E,F,G,H}.
# IMPORTANT: Branch IDs are NOT consistent across runs — IQ-TREE may
# renumber internal nodes when branch lengths change dramatically
# (e.g., when the target branch collapses to near-zero length).
# Therefore, we identify the target branch BY BIPARTITION in each
# cf.tree file, matching its sCF label to the corresponding cf.stat.
print("Target branch: {A,B,C,D}|{E,F,G,H} (identified per-file by bipartition)")

# ============================================================
# Step 3: Parse all results — identify target branch per-file
# ============================================================
rows = []
mismatches = 0

for X in SCF_LEVELS:
    for G in GAP_LEVELS:
        for model in MODELS:
            mname = model.replace("+", "_")
            tag = f"scf{X}_gap{G}_{mname}"

            # New binary methods
            for method in NEW_METHODS:
                tree_file = os.path.join(RESDIR, f"new_{tag}_{method}.cf.tree")
                stat_file = os.path.join(RESDIR, f"new_{tag}_{method}.cf.stat")
                target = find_target_branch_in_file(tree_file, stat_file)
                if target:
                    rows.append({
                        "scf_level": X,
                        "gap_level": G,
                        "model": model,
                        "binary": "new",
                        "method": method,
                        **target,
                    })

            # Old binary: scfl only
            tree_file = os.path.join(RESDIR, f"old_{tag}_scfl.cf.tree")
            stat_file = os.path.join(RESDIR, f"old_{tag}_scfl.cf.stat")
            target = find_target_branch_in_file(tree_file, stat_file)
            if target:
                rows.append({
                    "scf_level": X,
                    "gap_level": G,
                    "model": model,
                    "binary": "old",
                    "method": "old_scfl",
                    **target,
                })

df = pd.DataFrame(rows)

if df.empty:
    print("ERROR: No results found!")
    exit(1)

print(f"\nParsed {len(df)} data points")
print(f"Models: {df['model'].unique().tolist()}")
print(f"Methods: {df['method'].unique().tolist()}")
print(f"True tree proportions: {sorted(df['scf_level'].unique())}")
print(f"Gap levels: {sorted(df['gap_level'].unique())}")

# ============================================================
# Summary table
# ============================================================
print("\n" + "=" * 100)
print("SUMMARY: Mean sCF and sN by model, method, gap level (target branch only)")
print("=" * 100)

summary = df.groupby(["model", "method", "gap_level"]).agg(
    mean_sCF=("sCF", "mean"),
    mean_sN=("sN", "mean"),
).reset_index()

print(f"{'Model':<15s} {'Method':<12s} {'Gap%':>5s} {'mean_sCF':>10s} {'mean_sN':>10s}")
print("-" * 55)
for _, row in summary.sort_values(["model", "method", "gap_level"]).iterrows():
    print(f"{row['model']:<15s} {row['method']:<12s} {row['gap_level']:>5.0f} {row['mean_sCF']:>10.2f} {row['mean_sN']:>10.2f}")

# ============================================================
# Common setup
# ============================================================
ALL_METHODS = ["scf", "scfl", "scflg", "scfl_safe", "old_scfl"]
gap_cmap = plt.cm.viridis(np.linspace(0, 0.9, len(GAP_LEVELS)))
tree1_cmap = plt.cm.plasma(np.linspace(0, 0.9, len(SCF_LEVELS)))

def get_sub(method):
    """Get subset of df for a given method, handling old vs new binary."""
    if method == "old_scfl":
        return df[df["method"] == method]
    return df[(df["method"] == method) & (df["binary"] == "new")]

# Global sN max for consistent y-axes across all sN plots
sn_max = df["sN"].max() * 1.05

# ============================================================
# Figure 1: sN vs Gap% — rows=models, cols=methods, color=tree1%
# Expectation: sN depends on gap%, NOT on tree1%
# ============================================================
fig, axes = plt.subplots(len(MODELS), len(ALL_METHODS), figsize=(30, 24),
                         sharey=True, sharex=True)

for row_idx, model in enumerate(MODELS):
    for col_idx, method in enumerate(ALL_METHODS):
        ax = axes[row_idx, col_idx]
        sub = get_sub(method)
        msub = sub[sub["model"] == model]

        for ti, X in enumerate(SCF_LEVELS):
            xsub = msub[msub["scf_level"] == X].sort_values("gap_level")
            ax.plot(xsub["gap_level"], xsub["sN"], color=tree1_cmap[ti],
                    marker="o", markersize=4, linewidth=1.5, alpha=0.8,
                    label=f"tree1={X}%")

        ax.set_xticks(GAP_LEVELS)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0, sn_max)
        if row_idx == 0:
            ax.set_title(METHOD_LABELS[method], fontsize=11, fontweight="bold")
        if col_idx == 0:
            ax.set_ylabel(f"{model}\nsN (informative sites)", fontsize=10, fontweight="bold")
        if row_idx == len(MODELS) - 1:
            ax.set_xlabel("Gap %", fontsize=10)

# Collect handles from first panel, place single legend at top
handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.965),
           ncol=len(SCF_LEVELS), fontsize=8, frameon=True)
fig.suptitle("sN vs Gap%\n"
             "Rows = models, Columns = methods, Color = tree1%\n"
             "Tight line bundles confirm sN is independent of tree1%",
             fontsize=15, fontweight="bold", y=1.02)
plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig(os.path.join(SIMDIR, "fig1_sN_vs_gap.png"), dpi=150, bbox_inches="tight")
print(f"\nSaved: {SIMDIR}/fig1_sN_vs_gap.png")

# ============================================================
# Figure 2: sN vs Tree1% — rows=models, cols=methods, color=gap%
# Expectation: FLAT lines (sN independent of tree1%)
# ============================================================
fig, axes = plt.subplots(len(MODELS), len(ALL_METHODS), figsize=(30, 24),
                         sharey=True, sharex=True)

for row_idx, model in enumerate(MODELS):
    for col_idx, method in enumerate(ALL_METHODS):
        ax = axes[row_idx, col_idx]
        sub = get_sub(method)
        msub = sub[sub["model"] == model]

        for gi, G in enumerate(GAP_LEVELS):
            gsub = msub[msub["gap_level"] == G].sort_values("scf_level")
            ax.plot(gsub["scf_level"], gsub["sN"], color=gap_cmap[gi],
                    marker="o", markersize=4, linewidth=1.5, alpha=0.8,
                    label=f"gap={G}%")

        ax.set_xticks(SCF_LEVELS)
        ax.tick_params(axis='x', labelsize=7)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0, sn_max)
        if row_idx == 0:
            ax.set_title(METHOD_LABELS[method], fontsize=11, fontweight="bold")
        if col_idx == 0:
            ax.set_ylabel(f"{model}\nsN (informative sites)", fontsize=10, fontweight="bold")
        if row_idx == len(MODELS) - 1:
            ax.set_xlabel("Tree1 %", fontsize=10)

handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.965),
           ncol=len(GAP_LEVELS), fontsize=9, frameon=True)
fig.suptitle("sN vs Tree1%\n"
             "Rows = models, Columns = methods, Color = gap%\n"
             "Flat lines confirm sN is independent of tree topology",
             fontsize=15, fontweight="bold", y=1.02)
plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig(os.path.join(SIMDIR, "fig2_sN_vs_tree1.png"), dpi=150, bbox_inches="tight")
print(f"Saved: {SIMDIR}/fig2_sN_vs_tree1.png")

# ============================================================
# Figure 3: sCF vs Gap% — rows=models, cols=methods, color=tree1%
# Expectation: FLAT horizontal bands (except sCFLG declining)
# ============================================================
fig, axes = plt.subplots(len(MODELS), len(ALL_METHODS), figsize=(30, 24),
                         sharey=True, sharex=True)

for row_idx, model in enumerate(MODELS):
    for col_idx, method in enumerate(ALL_METHODS):
        ax = axes[row_idx, col_idx]
        sub = get_sub(method)
        msub = sub[sub["model"] == model]

        for ti, X in enumerate(SCF_LEVELS):
            xsub = msub[msub["scf_level"] == X].sort_values("gap_level")
            ax.plot(xsub["gap_level"], xsub["sCF"], color=tree1_cmap[ti],
                    marker="o", markersize=4, linewidth=1.5, alpha=0.8,
                    label=f"tree1={X}%")

        ax.set_xticks(GAP_LEVELS)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0, 100)
        if row_idx == 0:
            ax.set_title(METHOD_LABELS[method], fontsize=11, fontweight="bold")
        if col_idx == 0:
            ax.set_ylabel(f"{model}\nsCF (%)", fontsize=10, fontweight="bold")
        if row_idx == len(MODELS) - 1:
            ax.set_xlabel("Gap %", fontsize=10)

handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.965),
           ncol=len(SCF_LEVELS), fontsize=8, frameon=True)
fig.suptitle("sCF vs Gap%\n"
             "Rows = models, Columns = methods, Color = tree1%\n"
             "Flat bands confirm sCF is independent of gaps (except sCFLG)",
             fontsize=15, fontweight="bold", y=1.02)
plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig(os.path.join(SIMDIR, "fig3_scf_vs_gap.png"), dpi=150, bbox_inches="tight")
print(f"Saved: {SIMDIR}/fig3_scf_vs_gap.png")

# ============================================================
# Figure 4: sCF vs Tree1% — rows=models, cols=methods, color=gap%
# Expectation: sCF increases with tree1%, independent of gap% (except sCFLG)
# ============================================================
fig, axes = plt.subplots(len(MODELS), len(ALL_METHODS), figsize=(30, 24),
                         sharey=True, sharex=True)

for row_idx, model in enumerate(MODELS):
    for col_idx, method in enumerate(ALL_METHODS):
        ax = axes[row_idx, col_idx]
        sub = get_sub(method)
        msub = sub[sub["model"] == model]

        for gi, G in enumerate(GAP_LEVELS):
            gsub = msub[msub["gap_level"] == G].sort_values("scf_level")
            ax.plot(gsub["scf_level"], gsub["sCF"], color=gap_cmap[gi],
                    marker="o", markersize=4, linewidth=1.5, alpha=0.8,
                    label=f"gap={G}%")

        ax.set_xticks(SCF_LEVELS)
        ax.tick_params(axis='x', labelsize=7)
        ax.grid(True, alpha=0.3)
        ax.set_ylim(0, 100)
        if row_idx == 0:
            ax.set_title(METHOD_LABELS[method], fontsize=11, fontweight="bold")
        if col_idx == 0:
            ax.set_ylabel(f"{model}\nsCF (%)", fontsize=10, fontweight="bold")
        if row_idx == len(MODELS) - 1:
            ax.set_xlabel("Tree1 %", fontsize=10)

handles, labels = axes[0, 0].get_legend_handles_labels()
fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.965),
           ncol=len(GAP_LEVELS), fontsize=9, frameon=True)
fig.suptitle("sCF vs Tree1%\n"
             "Rows = models, Columns = methods, Color = gap%\n"
             "Overlapping lines confirm sCF is independent of gaps (sCFLG deviates)",
             fontsize=15, fontweight="bold", y=1.02)
plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig(os.path.join(SIMDIR, "fig4_scf_vs_tree1.png"), dpi=150, bbox_inches="tight")
print(f"Saved: {SIMDIR}/fig4_scf_vs_tree1.png")

# ============================================================
# Validation checks
# ============================================================
print("\n" + "=" * 100)
print("VALIDATION CHECKS")
print("=" * 100)

# Check 1: sCFL vs sCFL --safe should be identical (same seed)
scfl_check = df[(df["method"] == "scfl") & (df["binary"] == "new")][
    ["scf_level", "gap_level", "model", "sCF"]].rename(columns={"sCF": "sCFL"})
safe_check = df[(df["method"] == "scfl_safe") & (df["binary"] == "new")][
    ["scf_level", "gap_level", "model", "sCF"]].rename(columns={"sCF": "sCFL_safe"})
scfl_safe_merged = scfl_check.merge(safe_check, on=["scf_level", "gap_level", "model"])
if len(scfl_safe_merged) > 0:
    max_diff = (scfl_safe_merged["sCFL"] - scfl_safe_merged["sCFL_safe"]).abs().max()
    print(f"1. sCFL vs sCFL --safe max |diff|: {max_diff:.4f} {'PASS' if max_diff < 0.1 else 'FAIL'}")

# Check 2: New sCFL sN should be <= sCFLG sN
scfl_sn = df[(df["method"] == "scfl") & (df["binary"] == "new")].groupby(
    ["scf_level", "gap_level", "model"])["sN"].mean().reset_index().rename(columns={"sN": "scfl_sN"})
scflg_sn = df[(df["method"] == "scflg") & (df["binary"] == "new")].groupby(
    ["scf_level", "gap_level", "model"])["sN"].mean().reset_index().rename(columns={"sN": "scflg_sN"})
sn_merged = scfl_sn.merge(scflg_sn, on=["scf_level", "gap_level", "model"])
if len(sn_merged) > 0:
    violations = (sn_merged["scfl_sN"] > sn_merged["scflg_sN"] + 0.1).sum()
    print(f"2. sCFL sN <= sCFLG sN: {len(sn_merged) - violations}/{len(sn_merged)} pass ({violations} violations) {'PASS' if violations == 0 else 'WARN'}")

# Check 3: Old sCFL should be broken for multi-category models
for model in ["GTR+I+G8", "GTR+I+R4"]:
    old_sn = df[(df["model"] == model) & (df["method"] == "old_scfl")]["sN"].mean()
    new_sn = df[(df["model"] == model) & (df["method"] == "scfl") & (df["binary"] == "new")]["sN"].mean()
    print(f"3. {model}: old mean sN={old_sn:.1f}, new mean sN={new_sn:.1f} {'PASS (old broken)' if old_sn < new_sn * 0.5 else 'CHECK'}")

# Check 4: GTR (no rate het) should give same results old vs new
gtr_old = df[(df["model"] == "GTR") & (df["method"] == "old_scfl")]
gtr_new = df[(df["model"] == "GTR") & (df["method"] == "scfl") & (df["binary"] == "new")]
if len(gtr_old) > 0 and len(gtr_new) > 0:
    old_mean = gtr_old["sCF"].mean()
    new_mean = gtr_new["sCF"].mean()
    diff = abs(old_mean - new_mean)
    print(f"4. GTR (ncat=1): old mean sCF={old_mean:.2f}, new mean sCF={new_mean:.2f}, diff={diff:.2f} {'PASS' if diff < 1.0 else 'FAIL'}")

print("\nAnalysis complete.")
