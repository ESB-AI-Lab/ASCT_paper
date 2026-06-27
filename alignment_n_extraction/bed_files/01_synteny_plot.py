#!/usr/bin/env python3
"""
==============================================================================
 01_synteny_plot.py — ASCT Synteny Figure (pyGenomeViz)
==============================================================================
 Produces the main manuscript figure:
   Panel A: Gene-arrow synteny diagram with alignment links, ordered by
            species phylogeny. ASCT copies shown as colored arrows; absent
            copies are visually absent. Alignment ribbons connect homologous
            regions between adjacent species.
   Panel B: Presence/absence heatmap showing copy count per species alongside
            the species tree.

 Reads directly from:
   - {species}_ASCT_mRNA.bed files (gene positions)
   - Pairwise .paf or .blast6 files (alignment links)
   - phylo_order.tsv (tree-derived species ordering)
   - config.yaml (colors, layout settings)

 Usage:
   python 01_synteny_plot.py --config config.yaml
   python 01_synteny_plot.py --config config.yaml --identity-shading
==============================================================================
"""

import argparse
import os
import sys
from collections import defaultdict
from pathlib import Path

import yaml

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
import numpy as np

try:
    from pygenomeviz import GenomeViz
except ImportError:
    print("ERROR: pip install pygenomeviz")
    sys.exit(1)


def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)


# ─────────────────────────────────────────────────────────────────────────────
# LOAD DATA
# ─────────────────────────────────────────────────────────────────────────────

def load_phylo_order(base_dir: Path) -> list:
    """Load species order from phylo_order.tsv."""
    path = base_dir / "phylo_order.tsv"
    order = []
    with open(path) as f:
        next(f)  # header
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) >= 2:
                order.append(parts[1])
    return order


def parse_asct_bed(bed_path: str) -> list:
    """Parse _ASCT_mRNA.bed → list of copy dicts."""
    copies = []
    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            copies.append({
                "scaffold": parts[0],
                "start": int(parts[1]),
                "end": int(parts[2]),
                "name": parts[3],
                "strand": parts[5] if len(parts) > 5 else "+",
            })
    return copies


def load_all_copies(bed_dir: str, species_order: list) -> dict:
    """Load ASCT copies for all species. Returns {species: [copies]}."""
    all_copies = {}
    for sp in species_order:
        bed = Path(bed_dir) / f"{sp}_ASCT_mRNA.bed"
        if bed.exists():
            all_copies[sp] = parse_asct_bed(str(bed))
        else:
            all_copies[sp] = []  # absent
    return all_copies


def parse_paf(paf_path: str, min_len: int = 500, min_ident: float = 50) -> list:
    """Parse minimap2 PAF into link dicts."""
    links = []
    if not os.path.exists(paf_path):
        return links
    with open(paf_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 12:
                continue
            aln_len = int(parts[10])
            matches = int(parts[9])
            identity = (matches / aln_len * 100) if aln_len > 0 else 0
            if aln_len < min_len or identity < min_ident:
                continue
            links.append({
                "q_name": parts[0], "q_len": int(parts[1]),
                "q_start": int(parts[2]), "q_end": int(parts[3]),
                "strand": parts[4],
                "t_name": parts[5], "t_len": int(parts[6]),
                "t_start": int(parts[7]), "t_end": int(parts[8]),
                "identity": identity, "aln_len": aln_len,
            })
    return links


def parse_blast6(blast_path: str, min_len: int = 500, min_ident: float = 50) -> list:
    """Parse BLAST6 tabular into link dicts."""
    links = []
    if not os.path.exists(blast_path):
        return links
    with open(blast_path) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 12:
                continue
            identity = float(parts[2])
            aln_len = int(parts[3])
            if aln_len < min_len or identity < min_ident:
                continue
            links.append({
                "q_name": parts[0], "t_name": parts[1],
                "identity": identity, "aln_len": aln_len,
                "q_start": int(parts[6]), "q_end": int(parts[7]),
                "t_start": int(parts[8]), "t_end": int(parts[9]),
                "strand": "+" if int(parts[8]) < int(parts[9]) else "-",
            })
    return links


def load_alignments(align_dir: Path, sp_a: str, sp_b: str, cfg: dict) -> list:
    """Load alignment links between two species (try both name orders)."""
    min_len = cfg.get("min_alignment_length", 500)
    min_ident = cfg.get("min_identity", 50)

    for name_a, name_b in [(sp_a, sp_b), (sp_b, sp_a)]:
        paf = align_dir / f"{name_a}_vs_{name_b}.paf"
        if paf.exists():
            return parse_paf(str(paf), min_len, min_ident)
        blast = align_dir / f"{name_a}_vs_{name_b}.blast6"
        if blast.exists():
            return parse_blast6(str(blast), min_len, min_ident)
    return []


def merge_links(links: list, merge_distance: int = 10_000) -> list:
    """
    Merge nearby alignment links into larger blocks.

    Two links are merged if they are:
      - Same strand
      - Within merge_distance bp of each other on BOTH target and query
      - Collinear (same relative order on target and query)

    Merged link spans the union of both, identity is length-weighted average.
    This dramatically reduces visual clutter without losing information.
    """
    if not links or merge_distance <= 0:
        return links

    # Separate by strand
    plus = [l for l in links if l.get("strand", "+") == "+"]
    minus = [l for l in links if l.get("strand", "+") == "-"]

    def merge_group(group):
        if not group:
            return []
        # Sort by target start
        group.sort(key=lambda l: l["t_start"])

        merged = []
        current = dict(group[0])  # copy first link
        current_weight = current.get("aln_len", current["t_end"] - current["t_start"])

        for link in group[1:]:
            t_gap = link["t_start"] - current["t_end"]
            q_gap = link["q_start"] - current["q_end"]

            # Check if close enough on both axes and collinear
            if (0 <= t_gap <= merge_distance and
                abs(q_gap) <= merge_distance and
                q_gap >= -merge_distance):
                # Merge: extend current to cover both
                link_weight = link.get("aln_len", link["t_end"] - link["t_start"])
                total_weight = current_weight + link_weight
                # Weighted average identity
                current["identity"] = (
                    current["identity"] * current_weight +
                    link["identity"] * link_weight
                ) / total_weight
                current["t_end"] = max(current["t_end"], link["t_end"])
                current["q_end"] = max(current["q_end"], link["q_end"])
                current["aln_len"] = current["t_end"] - current["t_start"]
                current_weight = total_weight
            else:
                merged.append(current)
                current = dict(link)
                current_weight = link.get("aln_len", link["t_end"] - link["t_start"])

        merged.append(current)
        return merged

    result = merge_group(plus) + merge_group(minus)
    return result


# ─────────────────────────────────────────────────────────────────────────────
# COMPUTE REGION EXTENTS — reads actual extracted FASTA for correct coords
# ─────────────────────────────────────────────────────────────────────────────

def read_extracted_fasta_regions(extracted_dir: Path, species: str) -> list:
    """
    Read the extracted FASTA headers to get the actual extraction boundaries.
    Returns list of (seq_name, seq_length) tuples.

    Headers from extract_target_contigs.py look like:
      >NW_026854201.1:146995-354324 Rphi_NW_026854201.1
    Or for full-contig extraction:
      >NW_026854201.1  (no range, just the seq name)
    """
    for suffix in [".target_regions.fa", ".target_contigs.fa"]:
        fa_path = extracted_dir / f"{species}{suffix}"
        if fa_path.exists():
            seqs = []
            current_name = None
            current_len = 0
            with open(fa_path) as f:
                for line in f:
                    if line.startswith(">"):
                        if current_name is not None:
                            seqs.append((current_name, current_len))
                        current_name = line[1:].strip().split()[0]
                        current_len = 0
                    else:
                        current_len += len(line.strip())
            if current_name is not None:
                seqs.append((current_name, current_len))
            return seqs
    return []


def get_species_region(copies: list, flank: int, extracted_seqs: list = None,
                       max_track_size: int = 500_000):
    """
    From a species' ASCT copies and extracted FASTA regions, build a
    virtual track that concatenates ALL extracted regions across ALL
    scaffolds. Every copy is shown regardless of which scaffold it's on.
    """
    if not copies:
        return None

    # Build segments from ALL extracted FASTA sequences (not just one scaffold)
    gap = 5000  # visual gap between non-contiguous segments
    segments = []
    track_offset = 0

    if extracted_seqs:
        for seq_name, seq_len in extracted_seqs:
            if ":" in seq_name:
                parts = seq_name.split(":")
                scaffold_name = parts[0]
                range_parts = parts[1].split("-")
                seq_start = int(range_parts[0])
                seq_end = int(range_parts[1])
            else:
                scaffold_name = seq_name
                seq_start = 0
                seq_end = seq_len

            segments.append({
                "scaffold": scaffold_name,
                "seq_name": seq_name,
                "seq_start": seq_start,
                "seq_end": seq_end,
                "seq_len": seq_len,
                "track_offset": track_offset,
            })
            track_offset += seq_len + gap

    if not segments:
        # Fallback: single segment from BED coordinates
        gene_start = min(c["start"] for c in copies)
        gene_end = max(c["end"] for c in copies)
        display_start = max(0, gene_start - flank)
        track_size = (gene_end + flank) - display_start
        segments.append({
            "scaffold": copies[0]["scaffold"],
            "seq_name": copies[0]["scaffold"],
            "seq_start": display_start,
            "seq_end": display_start + track_size,
            "seq_len": track_size,
            "track_offset": 0,
        })
        track_offset = track_size

    track_size = track_offset - gap if len(segments) > 1 else track_offset

    return {
        "track_size": track_size,
        "segments": segments,
        "copies": copies,  # ALL copies, not just one scaffold
    }


def genomic_to_track(pos: int, segments: list):
    """Convert a genomic position to track coordinate using segment mapping."""
    for seg in segments:
        if seg["seq_start"] <= pos <= seg["seq_end"]:
            return seg["track_offset"] + (pos - seg["seq_start"])
    return None


def paf_to_track(paf_pos: int, seq_name: str, segments: list):
    """Convert a PAF coordinate (0-based within extracted seq) to track coordinate."""
    # Exact name match
    for seg in segments:
        if seg["seq_name"] == seq_name:
            return seg["track_offset"] + paf_pos
    # Partial match: PAF name might be a substring or vice versa
    for seg in segments:
        if seq_name in seg["seq_name"] or seg["seq_name"] in seq_name:
            return seg["track_offset"] + paf_pos
    # Scaffold match: PAF might just have scaffold name without range
    for seg in segments:
        if seg["scaffold"] == seq_name or seq_name.startswith(seg["scaffold"]):
            if paf_pos <= seg["seq_len"]:
                return seg["track_offset"] + paf_pos
    # Single segment fallback
    if len(segments) == 1 and paf_pos <= segments[0]["seq_len"]:
        return segments[0]["track_offset"] + paf_pos
    return None


# ─────────────────────────────────────────────────────────────────────────────
# PANEL A: SYNTENY PLOT
# ─────────────────────────────────────────────────────────────────────────────

def build_synteny_plot(species_order, all_copies, align_dir, cfg, out_dir):
    """Build the pyGenomeViz synteny figure."""
    flank = cfg.get("gene_flank_bp", cfg.get("flank_bp", 50_000))
    copy_colors = cfg.get("copy_colors", ["#E63946", "#457B9D", "#2A9D8F",
                                            "#E9C46A", "#F4A261", "#6A0572"])
    min_len = cfg.get("min_alignment_length", 500)
    min_ident = cfg.get("min_identity", 50)
    identity_shading = cfg.get("identity_shading", False)
    extracted_dir = Path(cfg["output_dir"]) / "00_extracted"
    max_track = cfg.get("max_track_size", 500_000)

    # Compute regions per species, using actual extracted FASTA boundaries
    regions = {}
    for sp in species_order:
        copies = all_copies.get(sp, [])
        extracted_seqs = read_extracted_fasta_regions(extracted_dir, sp)
        region = get_species_region(copies, flank, extracted_seqs)
        regions[sp] = region
        if region:
            n_segs = len(region["segments"])
            scaffolds = set(s["scaffold"] for s in region["segments"])
            print(f"    {sp}: track_size={region['track_size']:,} "
                  f"({n_segs} segment{'s' if n_segs > 1 else ''}) "
                  f"scaffolds={','.join(scaffolds)}")

    # Species with data (for alignment links)
    present_species = [sp for sp in species_order if regions[sp] is not None]

    # Determine a default track size for absent species
    present_sizes = [r["track_size"] for r in regions.values() if r is not None]
    avg_size = int(np.mean(present_sizes)) if present_sizes else 100_000

    # ── Auto-detect track orientation ──
    # Count forward vs inverted links per pair to suggest flips
    flip_species = set(cfg.get("flip_species", []))
    if cfg.get("auto_orient", True):
        print(f"\n  Auto-detecting track orientation...")
        # For each adjacent pair, count forward vs inverted links
        vote_flip = defaultdict(int)  # species → net votes for flipping
        for i in range(len(present_species) - 1):
            sp_a = present_species[i]
            sp_b = present_species[i + 1]
            links = load_alignments(align_dir, sp_a, sp_b, cfg)
            if not links:
                continue
            fwd = sum(1 for l in links if l.get("strand", "+") == "+")
            inv = sum(1 for l in links if l.get("strand", "+") != "+")
            if inv > fwd:
                # One of them should flip — vote for sp_b (keep sp_a fixed)
                vote_flip[sp_b] += (inv - fwd)
                print(f"    {sp_a} ↔ {sp_b}: {fwd} fwd, {inv} inv → suggest flip {sp_b}")
            else:
                print(f"    {sp_a} ↔ {sp_b}: {fwd} fwd, {inv} inv → ok")
        
        # Apply auto-flips
        for sp, votes in vote_flip.items():
            if votes > 0 and sp not in flip_species:
                flip_species.add(sp)
    
    if flip_species:
        print(f"  Flipping tracks: {', '.join(sorted(flip_species))}")

    # ── Build figure ──
    n_total = len(species_order)
    gv = GenomeViz(
        fig_track_height=cfg.get("track_height", 0.6),
        feature_track_ratio=0.25,
        link_track_ratio=0.5,
        track_align_type="left",
        fig_width=cfg.get("figure_width", 22),
    )

    track_names = {}
    # Store copy positions per species for ribbon coloring
    copy_track_positions = {}  # sp → [(track_start, track_end, color_idx)]

    for sp in species_order:
        region = regions[sp]

        if region is not None:
            size = region["track_size"]
            n = len(region["copies"])
            track_label = f"{sp}  ({n} {'copy' if n == 1 else 'copies'})"
            track = gv.add_feature_track(track_label, size, labelsize=14)
            track_names[sp] = track_label

            positions = []
            is_flipped = sp in flip_species
            for idx, copy in enumerate(region["copies"]):
                color = copy_colors[idx % len(copy_colors)]
                strand = 1 if copy["strand"] == "+" else -1
                rel_start = genomic_to_track(copy["start"], region["segments"])
                rel_end = genomic_to_track(copy["end"], region["segments"])
                if rel_start is None or rel_end is None:
                    continue
                # Flip coordinates if this species needs inverting
                if is_flipped:
                    rel_start, rel_end = size - rel_end, size - rel_start
                    strand = -strand
                positions.append((rel_start, rel_end, idx))
                track.add_feature(
                    start=rel_start, end=rel_end, strand=strand,
                    plotstyle="bigarrow",
                    fc=color, ec="black", lw=0.6,
                    label=copy["name"],
                    text_kws=dict(size=8, rotation=30),
                )
            copy_track_positions[sp] = positions
        else:
            # Absent species — empty track
            track_label = f"{sp}  (absent)"
            track = gv.add_feature_track(track_label, avg_size, labelsize=14)
            track_names[sp] = track_label
            track.add_feature(
                start=avg_size // 4, end=avg_size * 3 // 4, strand=0,
                plotstyle="box", fc="#F5F5F5", ec="#CCCCCC", lw=0.5,
                label="no ASCT detected",
                text_kws=dict(size=6, color="#999999"),
            )

    # ── Add alignment links between adjacent present species ──
    print(f"\n  Adding alignment links...")
    for i in range(len(present_species) - 1):
        sp_a = present_species[i]
        sp_b = present_species[i + 1]

        links = load_alignments(align_dir, sp_a, sp_b, cfg)
        if not links:
            print(f"    {sp_a} ↔ {sp_b}: no links")
            continue

        # Merge nearby fragments into larger blocks
        merge_dist = cfg.get("merge_distance", 10_000)
        raw_count = len(links)
        links = merge_links(links, merge_dist)

        name_a = track_names[sp_a]
        name_b = track_names[sp_b]
        segs_a = regions[sp_a]["segments"]
        segs_b = regions[sp_b]["segments"]
        size_a = regions[sp_a]["track_size"]
        size_b = regions[sp_b]["track_size"]

        added = 0
        skipped_range = 0
        first_error = None
        for link in links:
            # Map PAF coordinates to track coordinates via segments
            # PAF coords are 0-based within each extracted sequence
            t_start = paf_to_track(link["t_start"], link.get("t_name", ""), segs_a)
            t_end = paf_to_track(link["t_end"], link.get("t_name", ""), segs_a)
            q_start = paf_to_track(link["q_start"], link.get("q_name", ""), segs_b)
            q_end = paf_to_track(link["q_end"], link.get("q_name", ""), segs_b)

            if any(v is None for v in [t_start, t_end, q_start, q_end]):
                skipped_range += 1
                continue
            if t_end < 0 or t_start > size_a or q_end < 0 or q_start > size_b:
                skipped_range += 1
                continue
            t_start = max(0, min(t_start, size_a))
            t_end = max(0, min(t_end, size_a))
            q_start = max(0, min(q_start, size_b))
            q_end = max(0, min(q_end, size_b))

            # Apply flip for inverted species
            if sp_a in flip_species:
                t_start, t_end = size_a - t_end, size_a - t_start
            if sp_b in flip_species:
                q_start, q_end = size_b - q_end, size_b - q_start

            # Color ribbon by nearest gene copy
            ident = link["identity"]
            t_mid = (t_start + t_end) / 2
            
            # Find which copy the target midpoint is closest to
            color = "#AAAAAA"  # default grey
            best_dist = float("inf")
            positions_a = copy_track_positions.get(sp_a, [])
            for cp_start, cp_end, cp_idx in positions_a:
                cp_mid = (cp_start + cp_end) / 2
                d = abs(t_mid - cp_mid)
                if d < best_dist:
                    best_dist = d
                    color = copy_colors[cp_idx % len(copy_colors)]

            try:
                gv.add_link(
                    (name_a, t_start, t_end),
                    (name_b, q_start, q_end),
                    color=color, inverted_color="#D55E00",
                    curve=True, v=ident,
                    vmin=cfg.get("min_identity", 30),
                    vmax=100,
                    size=0.8,
                )
                added += 1
            except Exception as e:
                if first_error is None:
                    first_error = (e, link)

        if first_error and added == 0:
            err, lnk = first_error
            print(f"    {sp_a} ↔ {sp_b}: 0 links drawn — ERROR: {err}")
            print(f"      sample link: t={lnk['t_start']}-{lnk['t_end']} "
                  f"q={lnk['q_start']}-{lnk['q_end']} "
                  f"track_a_size={size_a} track_b_size={size_b}")
        elif added == 0 and skipped_range > 0:
            # All links skipped — show debug info
            sample = links[0]
            seg_names_a = [s["seq_name"][:40] for s in segs_a]
            seg_names_b = [s["seq_name"][:40] for s in segs_b]
            print(f"    {sp_a} ↔ {sp_b}: 0 links ({skipped_range} coord mismatch)")
            print(f"      PAF t_name: {sample.get('t_name', '?')[:50]}")
            print(f"      seg_a names: {seg_names_a}")
            print(f"      PAF q_name: {sample.get('q_name', '?')[:50]}")
            print(f"      seg_b names: {seg_names_b}")
        else:
            extra = f" ({skipped_range} out of range)" if skipped_range else ""
            merge_info = f" (merged from {raw_count})" if raw_count != len(links) else ""
            print(f"    {sp_a} ↔ {sp_b}: {added} links drawn{merge_info}{extra}")

    # ── Render ──
    fig = gv.plotfig()

    # Add legend for copy colors
    max_copies = max((len(all_copies.get(sp, [])) for sp in species_order), default=0)
    legend_patches = []
    for i in range(max_copies):
        color = copy_colors[i % len(copy_colors)]
        legend_patches.append(mpatches.Patch(
            facecolor=color, edgecolor="black", linewidth=0.5,
            label=f"Copy {i+1}"
        ))
    legend_patches.append(mpatches.Patch(
        facecolor="#F5F5F5", edgecolor="#CCCCCC", linewidth=0.5, label="Absent"
    ))

    if identity_shading:
        ident_colors = cfg.get("identity_colormap", {})
        legend_patches.extend([
            mpatches.Patch(facecolor=ident_colors.get("high", "#2A9D8F"), label="≥90% identity"),
            mpatches.Patch(facecolor=ident_colors.get("mid", "#E9C46A"), label="70-90%"),
            mpatches.Patch(facecolor=ident_colors.get("low", "#E63946"), label="<70%"),
        ])

    fig.legend(handles=legend_patches, loc="lower right", fontsize=6,
               title="ASCT Copies", title_fontsize=7, frameon=True,
               bbox_to_anchor=(0.98, 0.02))

    # Save
    pdf_path = out_dir / "panel_A_synteny.pdf"
    svg_path = out_dir / "panel_A_synteny.svg"
    fig.savefig(str(pdf_path), dpi=300, bbox_inches="tight")
    fig.savefig(str(svg_path), bbox_inches="tight")
    plt.close(fig)
    print(f"\n  Panel A saved: {pdf_path}")
    print(f"                 {svg_path}")

    return pdf_path


# ─────────────────────────────────────────────────────────────────────────────
# PANEL B: PRESENCE/ABSENCE HEATMAP WITH TREE
# ─────────────────────────────────────────────────────────────────────────────

def build_pa_heatmap(species_order, all_copies, cfg, out_dir):
    """
    Build a presence/absence + copy count heatmap alongside the species tree.
    Uses matplotlib directly (no R dependencies).
    """
    from Bio import Phylo
    import io, re

    tree_path = cfg.get("species_tree")
    copy_colors = cfg.get("copy_colors", ["#E63946", "#457B9D", "#2A9D8F",
                                            "#E9C46A", "#F4A261", "#6A0572"])

    n_species = len(species_order)
    copy_counts = [len(all_copies.get(sp, [])) for sp in species_order]
    max_copies = max(copy_counts) if copy_counts else 1

    # ── Parse tree manually (MrBayes Nexus with translate block) ──
    tree = None
    if tree_path and os.path.exists(tree_path):
        try:
            # Use the same manual parser that works in 00_preprocess.py
            translate_map = {}
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
                        if stripped == ";":
                            in_translate = False
                            continue
                        clean = stripped.rstrip(",").rstrip(";").strip()
                        if clean:
                            parts = clean.split()
                            if len(parts) >= 2:
                                translate_map[parts[0]] = parts[1]
                        continue
                    if in_trees and stripped.lower().startswith("tree "):
                        idx = stripped.find("=")
                        if idx >= 0:
                            newick_str = stripped[idx+1:].strip()
                            newick_str = re.sub(r'^\s*\[&[UR]\]\s*', '', newick_str)
                            if not newick_str.endswith(";"):
                                newick_str += ";"
                        break

            if newick_str and translate_map:
                # Replace integer IDs with taxon names
                newick_str = re.sub(
                    r'(?<=[(,])(\d+)(?=[\[:),])',
                    lambda m: translate_map.get(m.group(1), m.group(1)),
                    newick_str
                )
                # Strip MrBayes annotations
                newick_str = re.sub(r'\[&[^\]]*\]', '', newick_str)
                tree = Phylo.read(io.StringIO(newick_str), "newick")
                tree.ladderize()
        except Exception as e:
            print(f"    [WARN] Tree parse failed: {e}")
            tree = None

    # ── Figure layout ──
    has_tree = tree is not None
    n_cols = 2 if has_tree else 1
    width_ratios = [1.2, 2] if has_tree else [1]

    fig, axes = plt.subplots(
        1, n_cols,
        figsize=(max(8, max_copies * 1.5 + (5 if has_tree else 0)),
                 max(5, n_species * 0.4)),
        gridspec_kw={"width_ratios": width_ratios, "wspace": 0.1},
        squeeze=False,
    )
    axes = axes[0]
    ax_idx = 0

    # ── Tree panel ──
    if has_tree:
        ax_tree = axes[ax_idx]
        ax_idx += 1
        try:
            Phylo.draw(
                tree, axes=ax_tree, do_show=False,
                label_func=lambda c: c.name if c.name else "",
                label_colors=lambda name: "black",
            )
            ax_tree.set_ylabel("")
            ax_tree.set_xlabel("")
            ax_tree.set_title("Gene Tree", fontsize=9, fontweight="bold")
            ax_tree.spines["top"].set_visible(False)
            ax_tree.spines["right"].set_visible(False)
        except Exception as e:
            ax_tree.text(0.5, 0.5, f"Tree render failed:\n{e}",
                        ha="center", va="center", fontsize=7,
                        transform=ax_tree.transAxes)
            ax_tree.axis("off")

    # ── Heatmap panel ──
    ax_heat = axes[ax_idx]

    matrix = np.zeros((n_species, max_copies))
    for i, sp in enumerate(species_order):
        for j in range(copy_counts[i]):
            matrix[i, j] = 1

    for i in range(n_species):
        for j in range(max_copies):
            if matrix[i, j] == 1:
                color = copy_colors[j % len(copy_colors)]
                ax_heat.add_patch(plt.Rectangle(
                    (j, n_species - 1 - i), 1, 1,
                    facecolor=color, edgecolor="white", linewidth=1.5
                ))
            else:
                ax_heat.add_patch(plt.Rectangle(
                    (j, n_species - 1 - i), 1, 1,
                    facecolor="#F0F0F0", edgecolor="white", linewidth=1.5
                ))

    ax_heat.set_xlim(0, max_copies)
    ax_heat.set_ylim(0, n_species)
    ax_heat.set_xticks([x + 0.5 for x in range(max_copies)])
    ax_heat.set_xticklabels([f"Copy {i+1}" for i in range(max_copies)],
                             fontsize=8, rotation=45, ha="right")
    ax_heat.set_yticks([n_species - 1 - i + 0.5 for i in range(n_species)])
    ax_heat.set_yticklabels(species_order, fontsize=8)
    ax_heat.set_title("ASCT Copy Presence / Absence", fontsize=10, fontweight="bold")
    ax_heat.set_aspect("equal")

    for i, sp in enumerate(species_order):
        y = n_species - 1 - i + 0.5
        ax_heat.text(max_copies + 0.2, y, str(copy_counts[i]),
                    va="center", ha="left", fontsize=8, fontweight="bold")

    fig.subplots_adjust(bottom=0.12)
    pdf_path = out_dir / "panel_B_presence_absence.pdf"
    svg_path = out_dir / "panel_B_presence_absence.svg"
    fig.savefig(str(pdf_path), dpi=300, bbox_inches="tight")
    fig.savefig(str(svg_path), bbox_inches="tight")
    plt.close(fig)
    print(f"  Panel B saved: {pdf_path}")

    return pdf_path


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="ASCT synteny plot")
    parser.add_argument("--config", required=True)
    parser.add_argument("--identity-shading", action="store_true",
                        help="Color links by percent identity")
    args = parser.parse_args()

    cfg = load_config(args.config)
    if args.identity_shading:
        cfg["identity_shading"] = True

    base_dir = Path(cfg["output_dir"])
    bed_dir = cfg["bed_dir"]
    align_dir = base_dir / "00_alignments"
    out_dir = base_dir / "01_figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("STEP 1: ASCT SYNTENY FIGURE")
    print("=" * 70)

    # Load phylo order
    species_order = load_phylo_order(base_dir)
    print(f"\n  {len(species_order)} species in phylogenetic order")

    # Load all ASCT copies
    all_copies = load_all_copies(bed_dir, species_order)
    present = sum(1 for v in all_copies.values() if v)
    absent = sum(1 for v in all_copies.values() if not v)
    total_copies = sum(len(v) for v in all_copies.values())
    print(f"  {present} species with ASCT, {absent} without")
    print(f"  {total_copies} total copies")

    # Panel A: synteny
    print(f"\n  Building Panel A (synteny diagram)...")
    build_synteny_plot(species_order, all_copies, align_dir, cfg, out_dir)

    # Panel B: presence/absence heatmap
    print(f"\n  Building Panel B (presence/absence)...")
    build_pa_heatmap(species_order, all_copies, cfg, out_dir)

    print(f"\n{'=' * 70}")
    print("STEP 1 COMPLETE")
    print(f"{'=' * 70}")
    print(f"  Outputs: {out_dir}")
    print(f"  Next: python 02_composite_figure.py --config config.yaml")


if __name__ == "__main__":
    main()
