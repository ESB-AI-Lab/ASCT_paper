#!/usr/bin/env python3
"""
All-vs-all pairwise comparison of ASCT gene copies across species.
Produces:
  1. Clustered heatmap of percent identity
  2. Ribbon plot showing relationships between gene copies
"""

import subprocess
import tempfile
import os
import sys
from pathlib import Path
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

WORKDIR = Path(__file__).parent

SPECIES_COLORS = {
    "Asen": "#e6194b",
    "Mcal": "#3cb44b",
    "Mcor": "#4363d8",
    "Medi": "#f58231",
    "Mgal": "#911eb4",
    "Pvir": "#42d4f4",
}


def read_fasta(path):
    seqs = {}
    name = None
    lines = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(lines)
                name = line[1:].strip().split()[0]
                lines = []
            else:
                lines.append(line.strip())
    if name:
        seqs[name] = "".join(lines)
    return seqs


def pairwise_identity_mafft(seq1, seq2, name1, name2):
    """Align two sequences with MAFFT and return percent identity."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as tmp:
        tmp.write(f">{name1}\n{seq1}\n>{name2}\n{seq2}\n")
        tmp_path = tmp.name

    try:
        result = subprocess.run(
            ["mafft", "--quiet", "--auto", tmp_path],
            capture_output=True, text=True, timeout=60
        )
        aln_seqs = {}
        current = None
        lines = []
        for line in result.stdout.strip().split("\n"):
            if line.startswith(">"):
                if current:
                    aln_seqs[current] = "".join(lines)
                current = line[1:].strip().split()[0]
                lines = []
            else:
                lines.append(line.strip())
        if current:
            aln_seqs[current] = "".join(lines)

        s1 = aln_seqs[name1].upper()
        s2 = aln_seqs[name2].upper()
        matches = sum(1 for a, b in zip(s1, s2) if a == b and a != "-")
        aln_len = sum(1 for a, b in zip(s1, s2) if not (a == "-" and b == "-"))
        return (matches / aln_len * 100) if aln_len > 0 else 0.0
    finally:
        os.unlink(tmp_path)


def build_identity_matrix(genes):
    """Compute all-vs-all percent identity matrix."""
    names = list(genes.keys())
    n = len(names)
    matrix = np.zeros((n, n))
    np.fill_diagonal(matrix, 100.0)

    pairs = list(combinations(range(n), 2))
    total = len(pairs)

    for idx, (i, j) in enumerate(pairs):
        pct = pairwise_identity_mafft(
            genes[names[i]], genes[names[j]], names[i], names[j]
        )
        matrix[i, j] = pct
        matrix[j, i] = pct
        if (idx + 1) % 20 == 0 or idx + 1 == total:
            print(f"  {idx + 1}/{total} pairs done")

    return pd.DataFrame(matrix, index=names, columns=names)


def plot_heatmap(df, outpath):
    """Clustered heatmap of pairwise percent identity."""
    species_labels = [n.split("_ASCT")[0] for n in df.index]
    row_colors = [SPECIES_COLORS.get(s, "#999999") for s in species_labels]

    dist = 100.0 - df.values
    np.fill_diagonal(dist, 0)
    condensed = squareform(dist, checks=False)
    linkage_mat = linkage(condensed, method="average")

    g = sns.clustermap(
        df,
        row_linkage=linkage_mat,
        col_linkage=linkage_mat,
        cmap="RdYlGn",
        vmin=30,
        vmax=100,
        figsize=(14, 12),
        row_colors=row_colors,
        col_colors=row_colors,
        annot=True,
        fmt=".1f",
        annot_kws={"size": 7},
        linewidths=0.5,
        cbar_kws={"label": "Percent Identity (%)"},
    )
    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("")

    legend_patches = [mpatches.Patch(color=c, label=s) for s, c in SPECIES_COLORS.items()]
    g.ax_heatmap.legend(
        handles=legend_patches, loc="upper left",
        bbox_to_anchor=(1.06, 1.22), frameon=True, title="Species"
    )

    g.fig.suptitle("ASCT Gene Copy Pairwise Identity", y=1.02, fontsize=16, fontweight="bold")
    g.savefig(outpath, dpi=200, bbox_inches="tight")
    pdf_path = outpath.with_suffix(".pdf")
    g.savefig(pdf_path, bbox_inches="tight")
    plt.close()
    print(f"Heatmap saved: {outpath} + {pdf_path.name}")


def parse_bed_for_species(species):
    """Parse BED file and return entries for a species."""
    bed_path = WORKDIR / f"{species}_ASCT_mRNA.bed"
    if not bed_path.exists():
        return []
    entries = []
    with open(bed_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            entries.append({
                "scaffold": cols[0],
                "start": int(cols[1]),
                "name": cols[3],
                "strand": cols[5],
            })
    return entries


def get_display_order(species):
    """
    Determine display order from BED data:
    - Group by scaffold, sort within scaffold by start position
    - If all genes on a scaffold are minus-strand, reverse their order
    - Arrange scaffolds with ASCT1.x (clade A) scaffold first, then others
    Returns: list of (gene_name, break_before) tuples
    """
    entries = parse_bed_for_species(species)
    if not entries:
        return []

    scaffolds = {}
    for e in entries:
        scaffolds.setdefault(e["scaffold"], []).append(e)

    for scaf in scaffolds:
        scaffolds[scaf].sort(key=lambda e: e["start"])
        names = [e["name"] for e in scaffolds[scaf]]
        strands = [e["strand"] for e in scaffolds[scaf]]
        has_clade_a = any("_ASCT1." in n for n in names)
        if has_clade_a and "_ASCT1." not in names[0]:
            scaffolds[scaf].reverse()
        elif not has_clade_a and all(s == "-" for s in strands):
            scaffolds[scaf].reverse()

    def _gene_clade_priority(name):
        """ASCT1.x (clade A) = 0, ASCT2a (clade B) = 1, ASCT2b (clade C) = 2."""
        part = name.split("_ASCT")[-1] if "_ASCT" in name else ""
        if part.startswith("1."):
            return 0
        if part.startswith("2a"):
            return 1
        if part.startswith("2b"):
            return 2
        return 9

    def scaffold_priority(scaf):
        best = min(_gene_clade_priority(e["name"]) for e in scaffolds[scaf])
        return (best, scaf)

    ordered_scaffolds = sorted(scaffolds.keys(), key=scaffold_priority)

    result = []
    for si, scaf in enumerate(ordered_scaffolds):
        for ei, entry in enumerate(scaffolds[scaf]):
            break_before = (si > 0 and ei == 0)
            result.append((entry["name"], break_before))
    return result


def draw_break_symbol(ax, x, y, h=0.4):
    """Draw a publication-style zigzag break symbol."""
    zw = 0.12
    zh = h / 2
    n_teeth = 3
    tooth_h = zh / n_teeth
    xs = []
    ys = []
    yy = y + zh
    xs.append(x)
    ys.append(yy)
    for i in range(n_teeth):
        xs.append(x + zw)
        ys.append(yy - tooth_h * 0.5)
        yy -= tooth_h
        xs.append(x - zw)
        ys.append(yy + tooth_h * 0.5)
    xs.append(x)
    ys.append(y - zh)
    ax.plot(xs, ys, color="black", linewidth=1.5, solid_capstyle="round", zorder=5)


def plot_ribbon(df, outpath):
    """
    Ribbon plot: species as rows, gene copies as boxes ordered by
    genomic position from BED files. Contig breaks shown with zigzag.
    If all genes on a contig are minus-strand, display order is reversed.
    ASCT3's contig placed leftmost.
    """
    species_order = ["Asen", "Mcal", "Mcor", "Medi", "Mgal", "Pvir"]

    species_display = {}
    for sp in species_order:
        species_display[sp] = get_display_order(sp)

    species_genes = {sp: [g for g, _ in items] for sp, items in species_display.items() if items}

    fig, ax = plt.subplots(figsize=(18, 10))

    y_positions = {}
    box_positions = {}
    row_height = 1.0
    box_h = 0.5
    max_slots = max(
        sum(1 for _ in items) + sum(1 for _, brk in items if brk) * 0.4
        for items in species_display.values() if items
    )
    box_w = 2.2
    gap = 0.3
    break_gap = 0.5

    for row_idx, sp in enumerate(species_order):
        if sp not in species_genes:
            continue
        y = (len(species_order) - 1 - row_idx) * (row_height + 0.8)
        y_positions[sp] = y
        display = species_display[sp]

        total_width = 0
        for gi, (gname, brk) in enumerate(display):
            if gi > 0:
                total_width += break_gap if brk else gap
            total_width += box_w

        max_total = max(
            sum(box_w + (break_gap if brk and gi > 0 else gap if gi > 0 else 0)
                for gi, (_, brk) in enumerate(d))
            for d in species_display.values() if d
        )
        start_x = (max_total - total_width) / 2

        cur_x = start_x
        for gi, (gname, brk) in enumerate(display):
            if gi > 0:
                if brk:
                    draw_break_symbol(ax, cur_x + break_gap / 2 - gap / 2, y, box_h)
                    cur_x += break_gap
                else:
                    cur_x += gap

            box_positions[gname] = (cur_x, y)

            rect = FancyBboxPatch(
                (cur_x, y - box_h / 2), box_w, box_h,
                boxstyle="round,pad=0.05",
                facecolor=SPECIES_COLORS.get(sp, "#cccccc"),
                edgecolor="black", linewidth=1.2, alpha=0.9
            )
            ax.add_patch(rect)
            ax.text(cur_x + box_w / 2, y, gname,
                    ha="center", va="center", fontsize=7, fontweight="bold", color="white")
            cur_x += box_w

        ax.text(-1.5, y, sp, ha="right", va="center", fontsize=12,
                fontweight="bold", color=SPECIES_COLORS.get(sp, "#333333"))

    cmap = plt.cm.RdYlGn
    norm = plt.Normalize(vmin=30, vmax=100)

    ribbons = []
    for i, name_i in enumerate(df.index):
        sp_i = name_i.split("_ASCT")[0]
        for j, name_j in enumerate(df.columns):
            sp_j = name_j.split("_ASCT")[0]
            if sp_i >= sp_j:
                continue
            if name_i not in box_positions or name_j not in box_positions:
                continue
            pct = df.iloc[i, j]
            ribbons.append((name_i, name_j, pct))

    ribbons.sort(key=lambda r: r[2])

    top_n = 60
    if len(ribbons) > top_n:
        threshold = sorted([r[2] for r in ribbons], reverse=True)[top_n - 1]
        ribbons = [r for r in ribbons if r[2] >= threshold]

    from matplotlib.path import Path as MPath
    from collections import defaultdict

    box_ribbons = defaultdict(list)
    for idx, (name_i, name_j, pct) in enumerate(ribbons):
        box_ribbons[name_i].append((idx, name_j))
        box_ribbons[name_j].append((idx, name_i))

    margin = 0.15
    usable = box_w - 2 * margin
    ribbon_x = {}
    for gene, connections in box_ribbons.items():
        bx, by = box_positions[gene]
        connections.sort(key=lambda c: box_positions[c[1]][0])
        n = len(connections)
        for slot, (ribbon_idx, _) in enumerate(connections):
            if n == 1:
                x_pos = bx + box_w / 2
            else:
                x_pos = bx + margin + slot * usable / (n - 1)
            ribbon_x[(ribbon_idx, gene)] = x_pos

    for idx, (name_i, name_j, pct) in enumerate(ribbons):
        x1 = ribbon_x[(idx, name_i)]
        x2 = ribbon_x[(idx, name_j)]
        _, y1 = box_positions[name_i]
        _, y2 = box_positions[name_j]

        color = cmap(norm(pct))
        alpha = 0.15 + 0.6 * ((pct - 30) / 70)
        lw = 1 + 3 * ((pct - 30) / 70)

        mid_y = (y1 - box_h / 2 + y2 + box_h / 2) / 2
        verts = [
            (x1, y1 - box_h / 2),
            (x1, mid_y),
            (x2, mid_y),
            (x2, y2 + box_h / 2),
        ]
        codes = [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4]
        path = MPath(verts, codes)
        patch = mpatches.PathPatch(
            path, facecolor="none", edgecolor=color,
            linewidth=lw, alpha=alpha
        )
        ax.add_patch(patch)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.5, pad=0.02)
    cbar.set_label("Percent Identity (%)", fontsize=11)

    all_box_xs = [bx for bx, _ in box_positions.values()]
    ax.set_xlim(min(all_box_xs) - 3, max(all_box_xs) + box_w + 1.5)
    all_ys = list(y_positions.values())
    ax.set_ylim(min(all_ys) - 1.5, max(all_ys) + 1.5)
    ax.set_aspect("auto")
    ax.axis("off")
    ax.set_title("ASCT Gene Copy Relationships Across Species",
                 fontsize=16, fontweight="bold", pad=20)

    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    pdf_path = outpath.with_suffix(".pdf")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close()
    print(f"Ribbon plot saved: {outpath} + {pdf_path.name}")


def report_extremes(df):
    """Print most and least similar pairs."""
    names = list(df.index)
    pairs = []
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            pairs.append((names[i], names[j], df.iloc[i, j]))

    pairs.sort(key=lambda x: x[2], reverse=True)

    print("\n=== TOP 10 MOST SIMILAR PAIRS ===")
    for a, b, pct in pairs[:10]:
        print(f"  {a} <-> {b}: {pct:.1f}%")

    print("\n=== TOP 10 MOST DISSIMILAR PAIRS ===")
    for a, b, pct in pairs[-10:]:
        print(f"  {a} <-> {b}: {pct:.1f}%")

    print("\n=== PER-SPECIES SUMMARY ===")
    for sp in sorted(SPECIES_COLORS.keys()):
        sp_genes = [n for n in names if n.startswith(sp + "_")]
        if len(sp_genes) < 2:
            continue
        within = []
        for i, a in enumerate(sp_genes):
            for b in sp_genes[i+1:]:
                within.append((a, b, df.loc[a, b]))
        within.sort(key=lambda x: x[2], reverse=True)
        print(f"\n  {sp} (within-species):")
        for a, b, pct in within:
            print(f"    {a} <-> {b}: {pct:.1f}%")


def main():
    genes = {}
    for fasta in sorted(WORKDIR.glob("*_ASCT*_gene.fasta")):
        seqs = read_fasta(fasta)
        genes.update(seqs)

    print(f"Loaded {len(genes)} gene copies")
    print("Computing pairwise alignments with MAFFT...")

    df = build_identity_matrix(genes)
    df.to_csv(WORKDIR / "pairwise_identity.csv")
    print(f"Identity matrix saved: pairwise_identity.csv")

    report_extremes(df)

    print("\nGenerating heatmap...")
    plot_heatmap(df, WORKDIR / "ASCT_identity_heatmap.png")

    print("Generating ribbon plot...")
    plot_ribbon(df, WORKDIR / "ASCT_ribbon_plot.png")

    print("\nDone!")


if __name__ == "__main__":
    main()
