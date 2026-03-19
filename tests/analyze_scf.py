#!/usr/bin/env python3
"""Analyze sCF/sCFL/sCFLG results across models with different rate categories."""

import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')
import numpy as np

# Model definitions: (model_string, sanitized_name, ncat_mix)
# ncat_mix = number of rate/mixture categories that get summed in partial_lh
models = [
    ("TIM2+F",       "TIM2_F",       1),
    ("TIM2+F+I",     "TIM2_F_I",     2),  # invariant + 1 variable
    ("TIM2+F+G4",    "TIM2_F_G4",    4),
    ("TIM2+F+G8",    "TIM2_F_G8",    8),
    ("TIM2+F+I+G4",  "TIM2_F_I_G4",  5),  # invariant + 4 gamma
    ("TIM2+F+I+G8",  "TIM2_F_I_G8",  9),  # invariant + 8 gamma
    ("TIM2+F+R2",    "TIM2_F_R2",    2),
    ("TIM2+F+R4",    "TIM2_F_R4",    4),
    ("TIM2+F+I+R2",  "TIM2_F_I_R2",  3),  # invariant + 2 free rate
    ("TIM2+F+I+R4",  "TIM2_F_I_R4",  5),  # invariant + 4 free rate
]

methods = ["scf", "scfl", "scflg", "scfl_safe"]
method_labels = {"scf": "sCF (parsimony)", "scfl": "sCFL (likelihood)", "scflg": "sCFLG (lik+gaps)", "scfl_safe": "sCFL (--safe)"}

outdir = "tests/model_tests"

# Parse all results
rows = []
for model_str, mname, ncat in models:
    for method in methods:
        stat_file = os.path.join(outdir, f"{mname}_{method}.cf.stat")
        if not os.path.exists(stat_file):
            print(f"MISSING: {stat_file}")
            continue
        with open(stat_file) as f:
            for line in f:
                if line.startswith("#") or line.startswith("ID"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 10:
                    continue
                rows.append({
                    "model": model_str,
                    "ncat": ncat,
                    "method": method,
                    "method_label": method_labels[method],
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

df = pd.DataFrame(rows)

# ============================================================
# Summary table
# ============================================================
print("=" * 120)
print("SUMMARY TABLE: Mean sCF and Mean sN across all branches")
print("=" * 120)
summary = df.groupby(["model", "ncat", "method"]).agg(
    mean_sCF=("sCF", "mean"),
    mean_sN=("sN", "mean"),
    min_sN=("sN", "min"),
    max_sN=("sN", "max"),
).reset_index()
summary = summary.sort_values(["ncat", "model", "method"])

print(f"{'Model':<20s} {'ncat':>4s} {'Method':<6s} {'mean_sCF':>10s} {'mean_sN':>10s} {'min_sN':>10s} {'max_sN':>10s}")
print("-" * 80)
for _, row in summary.iterrows():
    print(f"{row['model']:<20s} {row['ncat']:>4d} {row['method']:<6s} {row['mean_sCF']:>10.2f} {row['mean_sN']:>10.2f} {row['min_sN']:>10.2f} {row['max_sN']:>10.2f}")

# ============================================================
# Detailed per-branch table for scfl
# ============================================================
print("\n" + "=" * 120)
print("PER-BRANCH sCFL VALUES (the problematic method)")
print("=" * 120)
pivot = df[df["method"] == "scfl"].pivot_table(
    index="branch_id", columns="model", values="sCF"
)
# Sort columns by ncat
model_order = [m[0] for m in sorted(models, key=lambda x: x[2])]
pivot = pivot[[c for c in model_order if c in pivot.columns]]
print(pivot.to_string(float_format="%.1f"))

print("\n" + "=" * 120)
print("PER-BRANCH sN for sCFL (number of informative sites)")
print("=" * 120)
pivot_sn = df[df["method"] == "scfl"].pivot_table(
    index="branch_id", columns="model", values="sN"
)
pivot_sn = pivot_sn[[c for c in model_order if c in pivot_sn.columns]]
print(pivot_sn.to_string(float_format="%.1f"))

# ============================================================
# Figure 1: sCF by ncat, one line per branch, faceted by method
# ============================================================
fig, axes = plt.subplots(1, 4, figsize=(24, 7), sharey=True)
branch_ids = sorted(df["branch_id"].unique())
cmap = plt.cm.tab20(np.linspace(0, 1, len(branch_ids)))

for ax_idx, method in enumerate(methods):
    ax = axes[ax_idx]
    sub = df[df["method"] == method]

    for i, bid in enumerate(branch_ids):
        bdata = sub[sub["branch_id"] == bid].sort_values("ncat")
        ax.plot(bdata["ncat"], bdata["sCF"], marker="o", markersize=3,
                color=cmap[i], alpha=0.7, linewidth=1, label=f"Branch {bid}")

    ax.set_xlabel("Number of rate categories (ncat_mix)", fontsize=11)
    ax.set_title(method_labels[method], fontsize=13, fontweight="bold")
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 100)
    ax.grid(True, alpha=0.3)
    if ax_idx == 0:
        ax.set_ylabel("sCF (%)", fontsize=11)

axes[3].legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7, ncol=1)
fig.suptitle("Site Concordance Factor vs. Number of Rate Categories", fontsize=14, fontweight="bold")
plt.tight_layout()
plt.savefig("tests/model_tests/fig1_scf_by_ncat.png", dpi=150, bbox_inches="tight")
print(f"\nSaved: tests/model_tests/fig1_scf_by_ncat.png")

# ============================================================
# Figure 2: sN by ncat, one line per branch, faceted by method
# ============================================================
fig, axes = plt.subplots(1, 4, figsize=(24, 7), sharey=True)

for ax_idx, method in enumerate(methods):
    ax = axes[ax_idx]
    sub = df[df["method"] == method]

    for i, bid in enumerate(branch_ids):
        bdata = sub[sub["branch_id"] == bid].sort_values("ncat")
        ax.plot(bdata["ncat"], bdata["sN"], marker="o", markersize=3,
                color=cmap[i], alpha=0.7, linewidth=1, label=f"Branch {bid}")

    ax.set_xlabel("Number of rate categories (ncat_mix)", fontsize=11)
    ax.set_title(method_labels[method], fontsize=13, fontweight="bold")
    ax.set_xlim(0, 10)
    ax.grid(True, alpha=0.3)
    if ax_idx == 0:
        ax.set_ylabel("sN (informative sites)", fontsize=11)

axes[3].legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7, ncol=1)
fig.suptitle("Number of Informative Sites vs. Number of Rate Categories", fontsize=14, fontweight="bold")
plt.tight_layout()
plt.savefig("tests/model_tests/fig2_sN_by_ncat.png", dpi=150, bbox_inches="tight")
print(f"Saved: tests/model_tests/fig2_sN_by_ncat.png")

# ============================================================
# Figure 3: Direct comparison scfl vs scflg vs scf, one row per model
# ============================================================
fig, axes = plt.subplots(len(models), 3, figsize=(15, 4 * len(models)), sharey=True, sharex=True)

for row_idx, (model_str, mname, ncat) in enumerate(models):
    scf_data = df[(df["model"] == model_str) & (df["method"] == "scf")].set_index("branch_id")["sCF"]
    scfl_data = df[(df["model"] == model_str) & (df["method"] == "scfl")].set_index("branch_id")["sCF"]
    scflg_data = df[(df["model"] == model_str) & (df["method"] == "scflg")].set_index("branch_id")["sCF"]
    scfl_safe_data = df[(df["model"] == model_str) & (df["method"] == "scfl_safe")].set_index("branch_id")["sCF"]

    # Panel A: sCFL vs sCF
    ax = axes[row_idx, 0]
    common = scf_data.index.intersection(scfl_data.index)
    ax.scatter(scf_data[common], scfl_data[common], s=25, alpha=0.7, color="#F44336")
    if len(common) > 1:
        z = np.polyfit(scf_data[common], scfl_data[common], 1)
        xfit = np.linspace(0, 100, 50)
        ax.plot(xfit, np.polyval(z, xfit), alpha=0.5, linewidth=1.5, color="#F44336")
    ax.plot([0, 100], [0, 100], "k--", alpha=0.3, linewidth=1)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.grid(True, alpha=0.3)
    ax.set_ylabel(f"{model_str}\n(ncat={ncat})", fontsize=9, fontweight="bold")
    if row_idx == 0:
        ax.set_title("sCFL vs sCF", fontsize=12, fontweight="bold")
    if row_idx == len(models) - 1:
        ax.set_xlabel("sCF (parsimony) %")

    # Panel B: sCFLG vs sCF
    ax = axes[row_idx, 1]
    common = scf_data.index.intersection(scflg_data.index)
    ax.scatter(scf_data[common], scflg_data[common], s=25, alpha=0.7, color="#4CAF50")
    if len(common) > 1:
        z = np.polyfit(scf_data[common], scflg_data[common], 1)
        xfit = np.linspace(0, 100, 50)
        ax.plot(xfit, np.polyval(z, xfit), alpha=0.5, linewidth=1.5, color="#4CAF50")
    ax.plot([0, 100], [0, 100], "k--", alpha=0.3, linewidth=1)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.grid(True, alpha=0.3)
    if row_idx == 0:
        ax.set_title("sCFLG vs sCF", fontsize=12, fontweight="bold")
    if row_idx == len(models) - 1:
        ax.set_xlabel("sCF (parsimony) %")

    # Panel C: sCFL --safe vs sCF
    ax = axes[row_idx, 2]
    common = scf_data.index.intersection(scfl_safe_data.index)
    if len(common) > 0:
        ax.scatter(scf_data[common], scfl_safe_data[common], s=25, alpha=0.7, color="#FF9800")
        if len(common) > 1:
            z = np.polyfit(scf_data[common], scfl_safe_data[common], 1)
            xfit = np.linspace(0, 100, 50)
            ax.plot(xfit, np.polyval(z, xfit), alpha=0.5, linewidth=1.5, color="#FF9800")
    ax.plot([0, 100], [0, 100], "k--", alpha=0.3, linewidth=1)
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.grid(True, alpha=0.3)
    if row_idx == 0:
        ax.set_title("sCFL --safe vs sCF", fontsize=12, fontweight="bold")
    if row_idx == len(models) - 1:
        ax.set_xlabel("sCF (parsimony) %")

fig.suptitle("Concordance Factor Method Comparison\n(each row = one model)", fontsize=14, fontweight="bold")
plt.tight_layout()
plt.savefig("tests/model_tests/fig3_method_comparison.png", dpi=150, bbox_inches="tight")
print(f"Saved: tests/model_tests/fig3_method_comparison.png")

# ============================================================
# Figure 5: sCFL vs sCFL --safe (should be identical with same seed)
# ============================================================
fig, ax = plt.subplots(figsize=(8, 8))
scfl_all = df[df["method"] == "scfl"][["model", "branch_id", "sCF"]].rename(columns={"sCF": "sCFL"})
safe_all = df[df["method"] == "scfl_safe"][["model", "branch_id", "sCF"]].rename(columns={"sCF": "sCFL_safe"})
merged = scfl_all.merge(safe_all, on=["model", "branch_id"])
if len(merged) > 0:
    ax.scatter(merged["sCFL"], merged["sCFL_safe"], s=20, alpha=0.6, c="#333")
    ax.plot([0, 100], [0, 100], "r--", alpha=0.5, linewidth=1)
    max_diff = (merged["sCFL"] - merged["sCFL_safe"]).abs().max()
    ax.set_title(f"sCFL vs sCFL --safe (max |diff| = {max_diff:.2f})", fontsize=13, fontweight="bold")
else:
    ax.set_title("sCFL vs sCFL --safe (no data)", fontsize=13, fontweight="bold")
ax.set_xlabel("sCFL (%)", fontsize=11)
ax.set_ylabel("sCFL --safe (%)", fontsize=11)
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.grid(True, alpha=0.3)
ax.set_aspect("equal")
plt.tight_layout()
plt.savefig("tests/model_tests/fig5_scfl_vs_safe.png", dpi=150, bbox_inches="tight")
print(f"Saved: tests/model_tests/fig5_scfl_vs_safe.png")

# ============================================================
# Figure 4: Mean sN by model, grouped by method
# ============================================================
fig, ax = plt.subplots(figsize=(12, 6))

x = np.arange(len(models))
width = 0.2
colors = {"scf": "#2196F3", "scfl": "#F44336", "scflg": "#4CAF50", "scfl_safe": "#FF9800"}

for i, method in enumerate(methods):
    means = []
    for model_str, mname, ncat in models:
        sub = df[(df["model"] == model_str) & (df["method"] == method)]
        means.append(sub["sN"].mean())
    bars = ax.bar(x + i * width, means, width, label=method_labels[method],
                  color=colors[method], alpha=0.8)

ax.set_xlabel("Model", fontsize=11)
ax.set_ylabel("Mean sN (informative sites)", fontsize=11)
ax.set_title("Mean Informative Sites by Model and Method", fontsize=14, fontweight="bold")
ax.set_xticks(x + width)
labels = [f"{m[0]}\n(ncat={m[2]})" for m in models]
ax.set_xticklabels(labels, fontsize=7, rotation=45, ha="right")
ax.legend()
ax.grid(True, alpha=0.3, axis="y")
plt.tight_layout()
plt.savefig("tests/model_tests/fig4_mean_sN_by_model.png", dpi=150, bbox_inches="tight")
print(f"Saved: tests/model_tests/fig4_mean_sN_by_model.png")

print("\nDone!")
