#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 11 22:44:26 2026

@author: ghilp
"""

from pathlib import Path
import json
import shlex
import numpy as np
import pandas as pd
from collections import defaultdict

# Source files — change to desired file paths

cif_path = Path(
    "ENTER .cif PATH HERE"  # update with your path to single .cif file
)

json_path = Path(
    "ENTER .json PATH HERE"  # update with your path to single .json file
)

# User-defined protein residue ranges
protein_ranges = {
    "TIMP-1": (START, END),  # enter TIMP-1 residues
    "MMP-3cd": (START, END),  # enter MMP-3cd residues
}

# Load AF3 pLDDT values

with open(json_path) as f:
    af3 = json.load(f)

atom_plddts = np.array(af3["atom_plddts"])

print("AF3 atom pLDDTs:", len(atom_plddts))

# Read .cif file

with open(cif_path, "r") as f:
    lines = f.readlines()


# Find atom_site section
for i, line in enumerate(lines):
    if line.strip() == "_atom_site.group_PDB":
        atom_header_start = i
        break


# Get atom_site headers
headers = []
j = atom_header_start

while j < len(lines):
    line = lines[j].strip()

    if line.startswith("_atom_site."):
        headers.append(line)
        j += 1
    else:
        break


# Parse atoms
atoms = []

for line in lines[j:]:
    line = line.strip()

    if not line or line.startswith("#"):
        continue

    if line.startswith("_") or line.startswith("loop_"):
        break

    values = shlex.split(line)

    if len(values) == len(headers):
        atoms.append(dict(zip(headers, values)))


print("CIF atoms:", len(atoms))

# Assign atom-level pLDDTs to residues

residue_plddts = defaultdict(list)

for atom, plddt in zip(atoms, atom_plddts):

    chain = atom["_atom_site.label_asym_id"]
    residue = atom["_atom_site.label_seq_id"]
    resname = atom["_atom_site.label_comp_id"]

    if residue in [".", "?"]:
        continue

    key = (chain, int(residue), resname)

    residue_plddts[key].append(plddt)

    # Calculate mean pLDDT per residue

results = []

for (chain, residue, resname), values in residue_plddts.items():

    results.append({
        "chain": chain,
        "residue": residue,
        "resname": resname,
        "n_atoms": len(values),
        "pLDDT": np.mean(values)
    })


df = pd.DataFrame(results)

df = df.sort_values(
    ["chain", "residue"]
).reset_index(drop=True)

# Calculate mean pLDDT for user-defined residue ranges
for protein, (start, end) in protein_ranges.items():
    subset = df[df["residue"].between(start, end)]

    print(f"\n{protein} ({start}–{end})")
    print("Residues:", len(subset))
    print("Mean pLDDT:", subset["pLDDT"].mean())
