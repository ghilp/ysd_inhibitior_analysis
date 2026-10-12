#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 11 22:38:55 2026

@author: ghilp
"""
# Setup
import os
import glob
import re
import gzip
import shutil
import numpy as np
import pandas as pd


class Config():
    PATH = "./FASTQ"  # update with your path to file or single FASTQ file
    TRANSLATE_FLAG = True
    MOTIF = None  # if searching for specific motif
    # if searchign for loop, update with motifs flanking loop or interest
    FLANKED_BY = ["MOTIF_1", "MOTIF_2"]


config = Config()


for file in glob.glob(os.path.join(config.PATH, "**", "*.gz"), recursive=True):
    print(file)
    # Open and decompress the gzipped file
    with gzip.open(file, 'rb') as f_in:
        print(f"Unzipped {file}")
        with open(file[:-3], 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

# Convert FASTQ data to CSV file with all 6 frame translation
# Codoon transtaion


def translate(codon):
    codon_dict = {
        "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
        "TCA": "S", "TCG": "S", "TCC": "S", "TCT": "S",
        "TAT": "Y", "TAC": "Y", "TAA": "-", "TAG": "-",
        "TGC": "C", "TGT": "C", "TGA": "-", "TGG": "W",
        "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
        "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
        "ACA": "T", "ACT": "T", "ACC": "T", "ACG": "T",
        "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
        "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G"
    }
    return codon_dict.get(codon, "")


def complement(base):
    complement_dict = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return complement_dict.get(base, base)


def reverse_complement(sequence):
    return ''.join(complement(base) for base in sequence[::-1])


def translate_sequence(sequence, frame=0):
    return ''.join(translate(sequence[i:i+3]) for i in range(frame, len(sequence) - 2, 3))

# FASTQ reading and processing


def read_fastq(fastq_path):
    """
    Reads a FASTQ file and returns a list of (ID, Sequence, Quality).
    """
    records = []
    with open(fastq_path, 'r') as fq:
        while True:
            id_line = fq.readline().strip()
            if not id_line:
                break
            seq = fq.readline().strip()
            fq.readline()  # skip '+'
            qual = fq.readline().strip()
            seq_id = id_line[1:] if id_line.startswith('@') else id_line
            records.append((seq_id, seq, qual))
    return records


def phred_score(qual_str):
    """
    Convert ASCII-encoded Phred33 string to average quality score.
    """
    scores = [ord(c) - 33 for c in qual_str]
    return np.mean(scores) if scores else np.nan


def process_fastq(fastq_path, translate_flag=True):
    """
    Process FASTQ --> DataFrame with translation & quality metrics.
    """
    records = read_fastq(fastq_path)
    if not records:
        print(f"No reads found in {fastq_path}")
        return None

    df = pd.DataFrame(records, columns=["ID", "Sequence", "Quality"])
    df["Seq_Length"] = df["Sequence"].str.len()
    df["Avg_Quality"] = df["Quality"].apply(phred_score)

    # Remove ambiguous sequences
    df = df[~df["Sequence"].str.contains('N')]

    if translate_flag:
        print("Translating sequences (6 frames)...")
        for i in range(1, 4):
            df[f'Frame_fwd_{i}'] = df["Sequence"].apply(
                lambda s: translate_sequence(s, frame=i-1))
            df[f'Frame_rev_{i}'] = df["Sequence"].apply(
                lambda s: translate_sequence(reverse_complement(s), frame=i-1))
    else:
        print("Skipping translation (--no-translate enabled)")

    return df


def save_to_csv(df, output_path):
    """
    Save DataFrame to plain CSV.
    """
    df.to_csv(output_path + ".csv", index=False)
    print(f"Saved to {output_path}.csv")


def batch_convert_fastq(input_folder, translate_flag=True):
    """
    Convert all FASTQ files in a folder.
    """
    fastq_files = glob.glob(os.path.join(input_folder, "**", "*.fastq"), recursive=True) + \
        glob.glob(os.path.join(input_folder, "**", "*.fq"), recursive=True)

    if not fastq_files:
        print("No FASTQ files found in the folder.")
        return

    for fastq_file in fastq_files:
        print(f"Processing {fastq_file}...")
        df = process_fastq(fastq_file, translate_flag)
        if df is not None:
            output_path = os.path.splitext(fastq_file)[0]
            save_to_csv(df, output_path)


translate_flag = config.TRANSLATE_FLAG

if os.path.isdir(config.PATH):
    batch_convert_fastq(config.PATH, translate_flag)
else:
    df = process_fastq(config.PATH, translate_flag)
    if df is not None:
        output_path = os.path.splitext(config.PATH)[0]
        save_to_csv(df, output_path)

# Look for motifs or flanking regions


def find_hits(csv_path, motif=None, flanked_by=None, output_dir=None):
    """
    Find reads containing a given motif or a loop flanked by specified sequences.

    Args:
        csv_path: Input CSV (from translated FASTQ data)
        motif: Motif string to search for (optional)
        flanked_by: Tuple or list [left_flank, right_flank] (optional)
        output_dir: Optional path for output folder
    """
    base_dir = os.path.dirname(csv_path)
    base_name = os.path.splitext(os.path.basename(csv_path))[0]

    # --- Output directory setup ---
    if output_dir is None:
        output_dir = os.path.join(base_dir, "outputs")
    os.makedirs(output_dir, exist_ok=True)

    print(f"\nAnalyzing {os.path.basename(csv_path)}...")
    df = pd.read_csv(csv_path)
    frame_cols = [col for col in df.columns if col.startswith("Frame_")]
    matches = []

    # Build regex for loop search or motif search
    if flanked_by and len(flanked_by) == 2:
        left, right = flanked_by
        regex = re.compile(f"{left}([A-Z]+?){right}")
        search_mode = "flanked"
        print(f"Searching for loops flanked by '{left}' and '{right}'")
    elif motif:
        regex = re.compile(re.escape(motif))
        search_mode = "motif"
        print(f"Searching for motif '{motif}'")
    else:
        print("Must provide either --motif or --flanked_by arguments.")
        return

    # Iterate through all frames and sequences
    for _, row in df.iterrows():
        for frame_col in frame_cols:
            aa_seq = str(row[frame_col])

            # Flanked loop search
            if search_mode == "flanked":
                for match in regex.finditer(aa_seq):
                    loop_seq = match.group(1)
                    matches.append({
                        "Read_ID": row["ID"],
                        "Frame": frame_col,
                        "Sequence": row["Sequence"],
                        "Translation": aa_seq,
                        "Flank_Left": left,
                        "Loop_Seq": loop_seq,
                        "Flank_Right": right,
                        "Loop_Length": len(loop_seq),
                        "Avg_Quality": row.get("Avg_Quality", None),
                        "Seq_Length": row.get("Seq_Length", None)
                    })

            # Simple motif search
            elif search_mode == "motif" and regex.search(aa_seq):
                matches.append({
                    "Read_ID": row["ID"],
                    "Frame": frame_col,
                    "Sequence": row["Sequence"],
                    "Translation": aa_seq,
                    "Motif": motif,
                    "Avg_Quality": row.get("Avg_Quality", None),
                    "Seq_Length": row.get("Seq_Length", None)
                })
                break  # Only first matching frame is relevant

    if not matches:
        print(f"No matches found in {csv_path}")
        return None

    matches_df = pd.DataFrame(matches)

    # --- Create summaries ---
    if search_mode == "flanked":
        # Summary 1: Frame-level stats
        frame_summary = (
            matches_df.groupby("Frame")
            .agg(
                Count=("Read_ID", "count"),
                Mean_Quality=("Avg_Quality", "mean"),
                Mean_Seq_Length=("Seq_Length", "mean")
            )
            .reset_index()
        )

        # Summary 2: Unique loop sequences
        loop_summary = (
            matches_df.groupby("Loop_Seq")
            .agg(
                Count=("Read_ID", "count"),
                Mean_Quality=("Avg_Quality", "mean"),
                Mean_Seq_Length=("Seq_Length", "mean"),
                Loop_Length=("Loop_Length", "first")
            )
            .reset_index()
            .sort_values(by="Count", ascending=False)
        )

        # Save outputs inside the output folder
        matches_csv = os.path.join(
            output_dir, f"{base_name}_loops_{left}_{right}.csv")
        frame_summary_csv = os.path.join(
            output_dir, f"{base_name}_loops_{left}_{right}_frame_summary.csv")
        loop_summary_csv = os.path.join(
            output_dir, f"{base_name}_loops_{left}_{right}_loop_summary.csv")

        matches_df.to_csv(matches_csv, index=False)
        frame_summary.to_csv(frame_summary_csv, index=False)
        loop_summary.to_csv(loop_summary_csv, index=False)

        print(f"Matches saved to {matches_csv}")
        print(f"Frame summary saved to {frame_summary_csv}")
        print(f"Loop summary saved to {loop_summary_csv}")

    else:
        # Motif mode summary
        motif_summary = (
            matches_df.groupby("Frame")
            .agg(
                Count=("Read_ID", "count"),
                Mean_Quality=("Avg_Quality", "mean"),
                Mean_Seq_Length=("Seq_Length", "mean")
            )
            .reset_index()
        )

        matches_csv = os.path.join(output_dir, f"{base_name}_{motif}_hits.csv")
        summary_csv = os.path.join(
            output_dir, f"{base_name}_{motif}_summary.csv")

        matches_df.to_csv(matches_csv, index=False)
        motif_summary.to_csv(summary_csv, index=False)

        print(f"Matches saved to {matches_csv}")
        print(f"Summary saved to {summary_csv}")

    return matches_df


def batch_find_hits(input_folder, motif=None, flanked_by=None):
    """
    Apply motif or loop search to all CSVs in a folder.
    """
    csv_files = glob.glob(os.path.join(input_folder, "*.csv"))
    if not csv_files:
        print("No CSV files found in the folder.")
        return

    # Create a single output folder in batch mode
    output_dir = os.path.join(input_folder, "outputs")
    os.makedirs(output_dir, exist_ok=True)

    for csv_file in csv_files:
        find_hits(csv_file, motif=motif,
                  flanked_by=flanked_by, output_dir=output_dir)


if os.path.isdir(config.PATH):
    batch_find_hits(config.PATH, motif=config.MOTIF,
                    flanked_by=config.FLANKED_BY)
else:
    find_hits(config.PATH, motif=config.MOTIF,
              flanked_by=config.FLANKED_BY, output_dir=None)
