#!/usr/bin/env python3
"""Inject random gaps into PHYLIP alignments at various proportions."""

import os
import random

random.seed(99)

ALNDIR = "tests/scfl_simulation/alignments"
SCF_LEVELS = [100, 90, 80, 70, 60, 50, 40, 30, 20, 10, 0]
GAP_LEVELS = [0, 5, 10, 15, 20]


def read_phylip(path):
    """Read relaxed PHYLIP alignment, return list of (name, sequence)."""
    with open(path) as f:
        lines = f.readlines()
    header = lines[0].strip().split()
    ntax, nchar = int(header[0]), int(header[1])
    seqs = []
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        parts = line.split(None, 1)
        if len(parts) == 2:
            seqs.append((parts[0], list(parts[1].replace(" ", ""))))
    assert len(seqs) == ntax, f"Expected {ntax} taxa, got {len(seqs)}"
    assert all(len(s) == nchar for _, s in seqs), "Sequence length mismatch"
    return ntax, nchar, seqs


def write_phylip(path, ntax, nchar, seqs):
    """Write relaxed PHYLIP alignment."""
    with open(path, "w") as f:
        f.write(f"{ntax} {nchar}\n")
        for name, seq in seqs:
            f.write(f"{name}  {''.join(seq)}\n")


for X in SCF_LEVELS:
    nogap_file = os.path.join(ALNDIR, f"scf{X}_nogaps.phy")
    if not os.path.exists(nogap_file):
        print(f"WARNING: Missing {nogap_file}, skipping")
        continue

    ntax, nchar, seqs = read_phylip(nogap_file)

    for G in GAP_LEVELS:
        if G == 0:
            # Just copy the no-gaps alignment
            outpath = os.path.join(ALNDIR, f"scf{X}_gap{G}.phy")
            write_phylip(outpath, ntax, nchar, seqs)
            continue

        prob = G / 100.0
        new_seqs = []
        for name, seq in seqs:
            new_seq = [("?" if random.random() < prob else c) for c in seq]
            new_seqs.append((name, new_seq))

        outpath = os.path.join(ALNDIR, f"scf{X}_gap{G}.phy")
        write_phylip(outpath, ntax, nchar, new_seqs)

    print(f"  sCF={X}%: injected gaps at {GAP_LEVELS}")

print("Gap injection complete.")
