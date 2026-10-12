#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import numpy as np

# User-defined TIMP-1 and MMP-3 positions
TIMP1_RANGE = (START, END)  # enter TIMP-1 residues
MMP3_RANGE = (START, END)  # enter MMP-3cd residues

# Path to your AF3 confidences JSON file
JSON_FILE = "ENTER .json PATH HERE"

# Load AF3 output

with open(JSON_FILE, "r") as f:
    data = json.load(f)

pae = np.array(data["pae"])
atom_plddts = np.array(data["atom_plddts"])
token_res_ids = np.array(data["token_res_ids"])

# Convert residue ranges

t1_start, t1_end = TIMP1_RANGE
m3_start, m3_end = MMP3_RANGE

t1_idx = np.arange(t1_start - 1, t1_end)
m3_idx = np.arange(m3_start - 1, m3_end)

# TIMP-1 / MMP-3cd PAE
# Forward direction: rows = TIMP-1, columns = MMP-3cd
t1_m3_pae = pae[np.ix_(t1_idx, m3_idx)]

# Reverse direction: rows = MMP-3cd, columns = TIMP-1
m3_t1_pae = pae[np.ix_(m3_idx, t1_idx)]

# Combine both directions
inter_region_pae = np.concatenate(
    [t1_m3_pae.flatten(), m3_t1_pae.flatten()]
)

print("\n--- TIMP-1 ↔ MMP-3cd PAE ---")
print(f"Mean: {np.mean(inter_region_pae):.2f} Å")
print(f"Median: {np.median(inter_region_pae):.2f} Å")
print(f"Min: {np.min(inter_region_pae):.2f} Å")
print(f"Max: {np.max(inter_region_pae):.2f} Å")

# Check PAE matrix shapes
print("\nPAE matrix shapes:")
print(f"TIMP-1 rows × MMP-3cd columns: {t1_m3_pae.shape}")
print(f"MMP-3cd rows × TIMP-1 columns: {m3_t1_pae.shape}")

# PAE threshold summary
print("\n--- PAE Thresholds ---")
for threshold in [5, 10, 15, 20]:
    fraction = np.mean(inter_region_pae < threshold) * 100
    print(
        f"PAE < {threshold:2d} Å: "
        f"{fraction:.1f}% of TIMP-1/MMP-3cd residue pairs"
    )
