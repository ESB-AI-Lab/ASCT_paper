#!/usr/bin/env python3
"""
==============================================================================
 00_preprocess.py — Pairwise Alignments + Tree Parsing
==============================================================================
 Takes the extracted contigs from extract_target_contigs.py and:
   1. Parses the species tree (MrBayes Nexus or Newick) to determine
      phylogenetic ordering
   2. Runs pairwise nucleotide alignments (minimap2 or BLASTn) between
      phylogenetically adjacent species
   3. Outputs alignment files ready for the synteny plot

 No OrthoFinder, no protein extraction — the ASCT BED files are the
 only annotation needed.

 Usage:
   python 00_preprocess.py --config config.yaml
==============================================================================
"""

import argparse
import os
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import yaml


def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)


# ─────────────────────────────────────────────────────────────────────────────
# TREE PARSING
# ─────────────────────────────────────────────────────────────────────────────

def parse_tree(tree_path: str, tip_map: dict = None) -> list:
    """
    Parse a MrBayes gene tree and derive species ordering.

    This tree has GENE COPY tips (e.g. Mcal_copy1, Rphi_copy2), not species.
    We:
      1. Parse the MrBayes Nexus translate block (integer → taxon name)
      2. Resolve integers in the Newick string to taxon names
      3. Ladderize the tree
      4. Walk tips in display order, extract species code from each
         (Mcal_copy1 → Mcal), and record first-appearance order

    Returns list of unique species codes in phylogenetic display order.
    """
    from Bio import Phylo
    import io
    import re

    tree_path = str(tree_path)
    print(f"  Tree file: {tree_path}")

    # ── Parse the Nexus file manually to get translate block + tree ──
    translate_map = {}   # int_id → taxon_name
    newick_str = None
    in_trees = False
    in_translate = False

    with open(tree_path) as f:
        for line in f:
            stripped = line.strip()

            if stripped.lower().startswith("begin trees"):
                in_trees = True
                continue

            if in_trees and stripped.lower().startswith("translate"):
                in_translate = True
                continue

            if in_translate:
                # Lines like: 1       Fhep_copy1,
                # or last:    40      Ofus_copy1
                # or:         ;  (end of translate)
                if stripped == ";":
                    in_translate = False
                    continue
                # Remove trailing comma or semicolon
                clean = stripped.rstrip(",").rstrip(";").strip()
                if clean:
                    parts = clean.split()
                    if len(parts) >= 2:
                        int_id = parts[0]
                        taxon = parts[1]
                        translate_map[int_id] = taxon
                continue

            if in_trees and stripped.lower().startswith("tree "):
                # Extract the Newick string after '='
                idx = stripped.find("=")
                if idx >= 0:
                    newick_str = stripped[idx+1:].strip()
                    # Remove [&U] or [&R] prefix
                    newick_str = re.sub(r'^\s*\[&[UR]\]\s*', '', newick_str)
                    if not newick_str.endswith(";"):
                        newick_str += ";"
                break

            if stripped.lower() == "end;":
                in_trees = False

    if newick_str is None:
        print("  [ERROR] Could not find tree string in Nexus file")
        return []

    print(f"  Translate block: {len(translate_map)} taxa")

    # ── Replace integer IDs with taxon names in the Newick string ──
    if translate_map:
        # Replace integers that appear as leaf labels (not inside branch lengths)
        # Strategy: match integers that are preceded by '(' or ',' and followed
        # by '[' (annotation), ':' (branch length), ',' or ')'
        def replace_id(match):
            int_id = match.group(1)
            return match.group(0).replace(int_id, translate_map.get(int_id, int_id))

        # Replace all integer leaf references
        # Pattern: after ( or , find a bare integer before [ or : or , or )
        newick_str = re.sub(
            r'(?<=[(,])(\d+)(?=[\[:),])',
            lambda m: translate_map.get(m.group(1), m.group(1)),
            newick_str
        )

    # ── Remove MrBayes annotations [&prob=...,length_mean=...] ──
    newick_str = re.sub(r'\[&[^\]]*\]', '', newick_str)

    # ── Parse the cleaned Newick ──
    tree = Phylo.read(io.StringIO(newick_str), "newick")
    tree.ladderize()

    # ── Extract tips ──
    tips = [tip.name for tip in tree.get_terminals()]
    print(f"  Tips ({len(tips)}): {', '.join(tips[:8])}...")

    # ── Apply tip_map if provided ──
    if tip_map:
        tips = [tip_map.get(t, t) for t in tips]

    # ── Derive species tree from gene tree using LCA distances ──
    # For each species pair, compute minimum tree distance between any
    # of their copies. Then use hierarchical clustering to get a proper
    # species-level ordering.

    species_tips = defaultdict(list)
    for tip in tree.get_terminals():
        name = tip.name
        if tip_map:
            name = tip_map.get(name, name)
        sp = name.rsplit("_copy", 1)[0] if "_copy" in name else name
        species_tips[sp].append(tip)

    species_list = sorted(species_tips.keys())
    n = len(species_list)
    print(f"  {n} species, computing pairwise LCA distances...")

    # Build species distance matrix (min pairwise copy distance)
    dist = [[0.0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            min_d = min(
                tree.distance(ta, tb)
                for ta in species_tips[species_list[i]]
                for tb in species_tips[species_list[j]]
            )
            dist[i][j] = min_d
            dist[j][i] = min_d

    # Hierarchical clustering (UPGMA) to get species ordering
    # Simple implementation — no scipy needed
    species_order = _upgma_order(species_list, dist)

    print(f"  Species order ({len(species_order)}): {', '.join(species_order)}")

    return species_order


def _upgma_order(labels, dist_matrix):
    """
    UPGMA clustering returning leaf order matching the dendrogram.
    Pure Python, no scipy needed.
    """
    n = len(labels)
    if n <= 1:
        return labels[:]
    if n == 2:
        return labels[:]

    # Active cluster indices, sizes, leaf orders, and distance matrix
    active = list(range(n))
    sizes = [1] * n
    leaf_orders = [[i] for i in range(n)]
    # Full distance matrix (will add rows as we merge)
    D = [row[:] for row in dist_matrix]

    next_id = n

    while len(active) > 1:
        # Find closest pair among active clusters
        min_d = float("inf")
        mi, mj = active[0], active[1]
        for idx_i, ci in enumerate(active):
            for idx_j, cj in enumerate(active):
                if idx_i >= idx_j:
                    continue
                if D[ci][cj] < min_d:
                    min_d = D[ci][cj]
                    mi, mj = ci, cj

        # Create new merged cluster
        new_id = next_id
        next_id += 1
        new_size = sizes[mi] + sizes[mj]
        new_leaf_order = leaf_orders[mi] + leaf_orders[mj]

        # Compute distances from new cluster to all other active clusters
        new_row = [0.0] * new_id
        for ck in active:
            if ck == mi or ck == mj:
                continue
            d_new = (D[mi][ck] * sizes[mi] + D[mj][ck] * sizes[mj]) / new_size
            new_row[ck] = d_new

        # Extend D matrix to include new cluster
        # Add new_row as a new row, and extend all existing rows
        for row in D:
            row.append(0.0)
        new_row.append(0.0)  # distance to self
        D.append(new_row)
        # Set distances from existing clusters to new cluster
        for ck in active:
            if ck == mi or ck == mj:
                continue
            D[ck][new_id] = new_row[ck]

        sizes.append(new_size)
        leaf_orders.append(new_leaf_order)

        # Update active: remove mi, mj; add new_id
        active = [c for c in active if c != mi and c != mj]
        active.append(new_id)

    # Final leaf order
    final_order = leaf_orders[active[0]]
    return [labels[i] for i in final_order]


def write_phylo_order(tip_order: list, out_path: Path):
    """Write species ordering to a TSV for downstream scripts."""
    with open(out_path, "w") as f:
        f.write("phylo_order\tspecies\n")
        for i, sp in enumerate(tip_order, 1):
            f.write(f"{i}\t{sp}\n")
    print(f"  Phylo order written: {out_path}")


# ─────────────────────────────────────────────────────────────────────────────
# PAIRWISE ALIGNMENTS
# ─────────────────────────────────────────────────────────────────────────────

def find_extracted_fasta(species: str, extracted_dir: Path):
    """Find the extracted FASTA for a species."""
    for suffix in [".target_regions.fa", ".target_contigs.fa"]:
        p = extracted_dir / f"{species}{suffix}"
        if p.exists():
            return p
    return None


def get_minimap2_params(pct_divergence: int) -> list:
    """
    Map a percent divergence estimate to minimap2 parameters.

    IMPORTANT: minimap2 preset names are misleading about what divergence
    they actually handle (documented inconsistencies in issues #560, #638,
    #1188). Calibrated here to real-world effective ranges:

      asm5  actually handles ~0.1-1% divergence
      asm10 actually handles ~1-5%
      asm20 actually handles ~5-10%
      >10% requires custom parameters with smaller k, lower penalties

    Profiles below are mapped to ACTUAL divergence tolerance.
    """
    profiles = {
        1:  ["-x", "asm5"],      # <1% — asm5 is appropriate
        5:  ["-x", "asm10"],     # 1-5% — asm10 is appropriate
        10: ["-x", "asm20"],     # 5-10% — asm20's real ceiling
        15: ["-k", "17", "-w", "8",  "-A", "1", "-B", "3", "-O", "4,18",
             "-E", "2,1", "-s", "150", "-z", "200,100", "-m", "25", "-n", "2"],
        20: ["-k", "15", "-w", "7",  "-A", "1", "-B", "3", "-O", "4,18",
             "-E", "2,1", "-s", "100", "-z", "200,100", "-m", "20", "-n", "2"],
        30: ["-k", "11", "-w", "5",  "-A", "1", "-B", "2", "-O", "3,13",
             "-E", "2,1", "-s", "50",  "-z", "200,50",  "-m", "15", "-n", "2"],
        40: ["-k", "9",  "-w", "3",  "-A", "1", "-B", "2", "-O", "2,10",
             "-E", "2,1", "-s", "40",  "-z", "150,50",  "-m", "10", "-n", "2"],
        50: ["-k", "7",  "-w", "3",  "-A", "1", "-B", "2", "-O", "2,10",
             "-E", "2,1", "-s", "30",  "-z", "150,50",  "-m", "10", "-n", "2"],
    }

    # Pick the closest profile at or above the requested divergence
    thresholds = sorted(profiles.keys())
    chosen = thresholds[-1]
    for t in thresholds:
        if pct_divergence <= t:
            chosen = t
            break

    return profiles[chosen]


def run_minimap2(fa_a: Path, fa_b: Path, out_paf: Path, threads: int = 4,
                 minimap2_bin: str = "minimap2", pct_divergence: int = 20):
    """Run minimap2 with parameters tuned for the expected divergence level."""
    params = get_minimap2_params(pct_divergence)
    cmd = [
        minimap2_bin,
        *params,
        "--secondary=yes",
        "-N", "50",
        "-p", "0.3",
        "-t", str(threads),
        str(fa_a), str(fa_b)
    ]
    subprocess.run(cmd, check=True,
                   stdout=open(out_paf, "w"), stderr=subprocess.DEVNULL)


def run_blastn(fa_a: Path, fa_b: Path, out_blast: Path):
    """Run BLASTn between two FASTAs, output tabular."""
    db_prefix = str(fa_a) + ".blastdb"
    subprocess.run([
        "makeblastdb", "-in", str(fa_a),
        "-dbtype", "nucl", "-out", db_prefix
    ], check=True, capture_output=True)
    subprocess.run([
        "blastn", "-query", str(fa_b), "-db", db_prefix,
        "-outfmt", "6", "-evalue", "1e-10",
        "-max_target_seqs", "5",
        "-out", str(out_blast)
    ], check=True)
    # Cleanup db files
    for ext in [".nhr", ".nin", ".nsq", ".ndb", ".not", ".ntf", ".nto"]:
        p = Path(db_prefix + ext)
        if p.exists():
            p.unlink()


def run_pairwise_alignments(species_order: list, extracted_dir: Path,
                            align_dir: Path, cfg: dict):
    """
    Run pairwise alignments between phylogenetically adjacent species.
    """
    align_dir.mkdir(parents=True, exist_ok=True)
    aligner = cfg.get("aligner", "minimap2")
    threads = cfg.get("threads", 4)

    pairs_run = 0
    pairs_skipped = 0

    for i in range(len(species_order) - 1):
        sp_a = species_order[i]
        sp_b = species_order[i + 1]
        fa_a = find_extracted_fasta(sp_a, extracted_dir)
        fa_b = find_extracted_fasta(sp_b, extracted_dir)

        if fa_a is None or fa_b is None:
            missing = sp_a if fa_a is None else sp_b
            print(f"    [SKIP] {sp_a} vs {sp_b} — no FASTA for {missing}")
            pairs_skipped += 1
            continue

        pair_name = f"{sp_a}_vs_{sp_b}"
        print(f"    {pair_name}...", end=" ", flush=True)

        if aligner == "minimap2":
            out_file = align_dir / f"{pair_name}.paf"
            minimap2_bin = cfg.get("minimap2_path", "minimap2")
            pct_div = cfg.get("pct_divergence", 20)
            run_minimap2(fa_a, fa_b, out_file, threads, minimap2_bin, pct_div)
        else:
            out_file = align_dir / f"{pair_name}.blast6"
            run_blastn(fa_a, fa_b, out_file)

        # Count alignments
        n_alns = sum(1 for _ in open(out_file) if not _.startswith("#"))
        print(f"{n_alns} alignments")
        pairs_run += 1

    print(f"\n    {pairs_run} pairs aligned, {pairs_skipped} skipped")


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Preprocess: tree + pairwise alignments")
    parser.add_argument("--config", required=True)
    parser.add_argument("--skip-alignments", action="store_true")
    args = parser.parse_args()

    cfg = load_config(args.config)
    base_dir = Path(cfg["output_dir"])
    extracted_dir = base_dir / "00_extracted"
    align_dir = base_dir / "00_alignments"

    print("=" * 70)
    print("STEP 0: PREPROCESS — Tree + Pairwise Alignments")
    print("=" * 70)

    # ── Parse species tree ──
    tree_path = cfg.get("species_tree")
    tip_map = cfg.get("tree_tip_map")

    if tree_path and os.path.exists(tree_path):
        print(f"\n  Parsing species tree...")
        species_order = parse_tree(tree_path, tip_map)
    else:
        # Fallback: alphabetical order
        print(f"\n  [WARN] No species tree found, using alphabetical order")
        species_codes = sorted(cfg.get("species_map", {}).values())
        species_order = species_codes

    # Filter to species that have extracted FASTAs
    available = []
    for sp in species_order:
        if find_extracted_fasta(sp, extracted_dir) is not None:
            available.append(sp)
        else:
            print(f"    [INFO] {sp} has no extracted FASTA, will appear as absent")

    print(f"\n  Phylogenetic order ({len(available)} with data, "
          f"{len(species_order) - len(available)} absent):")
    for i, sp in enumerate(species_order, 1):
        marker = "✓" if sp in available else "✗"
        print(f"    {i:2d}. {sp} {marker}")

    # Save full order (including absent species — they show in the plot)
    write_phylo_order(species_order, base_dir / "phylo_order.tsv")

    # ── Pairwise alignments ──
    if not args.skip_alignments:
        print(f"\n  Running pairwise alignments ({cfg.get('aligner', 'minimap2')})...")
        run_pairwise_alignments(available, extracted_dir, align_dir, cfg)
    else:
        print("\n  Skipping alignments (--skip-alignments)")

    print(f"\n{'=' * 70}")
    print("PREPROCESS COMPLETE")
    print(f"{'=' * 70}")
    print(f"  Next: python 01_synteny_plot.py --config config.yaml")


if __name__ == "__main__":
    main()
