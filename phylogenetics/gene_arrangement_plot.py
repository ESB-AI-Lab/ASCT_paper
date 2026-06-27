#!/usr/bin/env python3
"""
Gene arrangement plot: colored blocks for each ASCT gene copy,
with intergenic distances shown as log-scaled connecting line lengths.
"""

import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
from pathlib import Path

WORKDIR = Path(__file__).parent

SPECIES_ORDER = ["Pvir", "Mcor", "Medi", "Mgal", "Mcal", "Asen"]

SPECIES_FULL = {
    "Pvir": "Perna viridis",
    "Mcor": "Mytilus coruscus",
    "Medi": "Mytilus edulis",
    "Mgal": "Mytilus galloprovincialis",
    "Mcal": "Mytilus californianus",
    "Asen": "Arcuatula senhousia",
}

GENE_COLORS = {
    "ASCT1": "#4363d8",
    "ASCT2": "#42d4f4",
    "ASCT3": "#3cb44b",
}


def get_gene_color(name):
    for key in ("ASCT3", "ASCT2", "ASCT1"):
        if key in name:
            return GENE_COLORS[key]
    return "#808080"


def parse_bed(path):
    entries = []
    with open(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            entries.append({
                "scaffold": cols[0],
                "start": int(cols[1]),
                "end": int(cols[2]),
                "name": cols[3],
                "strand": cols[5],
            })
    return entries


def format_distance(bp):
    if bp >= 1_000_000:
        return f"{bp / 1_000_000:.1f} Mb"
    if bp >= 1_000:
        return f"{bp / 1_000:.1f} kb"
    return f"{bp} bp"


def gap_to_line_length(distance_bp):
    if distance_bp <= 0:
        return 0.4
    return max(0.4, min(2.5, 0.35 * (math.log10(distance_bp) - 3)))


def draw_break_symbol(ax, x, y, h=0.3):
    zw, n_teeth = 0.08, 3
    zh = h / 2
    tooth_h = zh / n_teeth
    xs, ys = [x], [y + zh]
    yy = y + zh
    for _ in range(n_teeth):
        xs.append(x + zw)
        ys.append(yy - tooth_h * 0.5)
        yy -= tooth_h
        xs.append(x - zw)
        ys.append(yy + tooth_h * 0.5)
    xs.append(x)
    ys.append(y - zh)
    ax.plot(xs, ys, color="black", lw=1.5, solid_capstyle="round", zorder=5)


def main():
    n_species = len(SPECIES_ORDER)
    fig, ax = plt.subplots(figsize=(18, n_species * 1.2 + 1.5))

    row_height = 1.4
    box_h = 0.5
    box_w = 1.6
    break_gap = 0.7
    max_x = 0

    for row_idx, sp in enumerate(SPECIES_ORDER):
        bed_path = WORKDIR / f"{sp}_ASCT_mRNA.bed"
        if not bed_path.exists():
            continue
        entries = parse_bed(bed_path)
        if not entries:
            continue

        y = (n_species - 1 - row_idx) * row_height

        # Group by scaffold, sort within scaffold by start
        scaffolds = {}
        for e in entries:
            scaffolds.setdefault(e["scaffold"], []).append(e)
        for scaf in scaffolds:
            scaffolds[scaf].sort(key=lambda e: e["start"])

        # Scaffold with ASCT3 first, then others by smallest start
        def scaffold_priority(scaf):
            has_asct3 = any("ASCT3" in e["name"] for e in scaffolds[scaf])
            return (0 if has_asct3 else 1, min(e["start"] for e in scaffolds[scaf]))

        scaffold_order = sorted(scaffolds.keys(), key=scaffold_priority)

        cur_x = 4.5

        for si, scaf in enumerate(scaffold_order):
            genes = scaffolds[scaf]

            if si > 0:
                cur_x += 0.25
                draw_break_symbol(ax, cur_x + break_gap / 2, y, box_h)
                cur_x += break_gap + 0.25

            for gi, gene in enumerate(genes):
                if gi > 0:
                    prev = genes[gi - 1]
                    distance = gene["start"] - prev["end"]
                    if distance < 0:
                        distance = 0
                    line_len = gap_to_line_length(distance)
                    line_start = cur_x
                    cur_x += line_len

                    ax.plot([line_start, cur_x], [y, y],
                            color='black', lw=1, zorder=1)
                    ax.text((line_start + cur_x) / 2, y + box_h / 2 + 0.08,
                            format_distance(distance),
                            fontsize=6, ha='center', va='bottom',
                            color='#555555', style='italic')

                color = get_gene_color(gene["name"])
                rect = FancyBboxPatch(
                    (cur_x, y - box_h / 2), box_w, box_h,
                    boxstyle="round,pad=0.03",
                    facecolor=color, edgecolor='black',
                    lw=1.0, alpha=0.9, zorder=2,
                )
                ax.add_patch(rect)

                display = gene["name"].replace(f"{sp}_", "")
                ax.text(cur_x + box_w / 2, y, display,
                        fontsize=7.5, ha='center', va='center',
                        color='white', fontweight='bold', zorder=3)

                arrow_y = y - box_h / 2 - 0.18
                if gene["strand"] == "+":
                    ax.annotate("",
                                xy=(cur_x + box_w * 0.8, arrow_y),
                                xytext=(cur_x + box_w * 0.2, arrow_y),
                                arrowprops=dict(arrowstyle="-|>", color='black', lw=1.5))
                else:
                    ax.annotate("",
                                xy=(cur_x + box_w * 0.2, arrow_y),
                                xytext=(cur_x + box_w * 0.8, arrow_y),
                                arrowprops=dict(arrowstyle="-|>", color='black', lw=1.5))

                cur_x += box_w

        max_x = max(max_x, cur_x)

        full_name = SPECIES_FULL.get(sp, sp)
        ax.text(0.3, y, full_name, ha='left', va='center',
                fontsize=11, fontweight='bold', style='italic')

    legend_patches = [
        mpatches.Patch(color=GENE_COLORS["ASCT3"], label="ASCT3"),
        mpatches.Patch(color=GENE_COLORS["ASCT1"], label="ASCT1"),
        mpatches.Patch(color=GENE_COLORS["ASCT2"], label="ASCT2"),
    ]
    ax.legend(handles=legend_patches, loc='lower right',
              frameon=True, fontsize=9, title="Gene type")

    ax.text(0.98, 0.01,
            'Intergenic distances shown as log-scaled line lengths\nZigzag = different scaffold',
            transform=ax.transAxes, fontsize=8, style='italic',
            ha='right', va='bottom', color='#555555')

    ax.set_xlim(-0.5, max_x + 1.5)
    ax.set_ylim(-1.5, n_species * row_height + 0.5)
    ax.axis('off')
    ax.set_title("ASCT Gene Arrangement Across Species",
                 fontsize=16, fontweight='bold', pad=20)

    plt.tight_layout()
    fig.savefig(WORKDIR / "ASCT_gene_arrangement.png", dpi=200, bbox_inches='tight')
    fig.savefig(WORKDIR / "ASCT_gene_arrangement.pdf", bbox_inches='tight')
    plt.show()
    print("Saved: ASCT_gene_arrangement.png + .pdf")


if __name__ == "__main__":
    main()
