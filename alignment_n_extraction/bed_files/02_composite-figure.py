#!/usr/bin/env python3
"""
==============================================================================
 02_composite_figure.py — Assemble Manuscript Figure
==============================================================================
 Combines Panel A (synteny diagram) and Panel B (PA heatmap) into a
 single multi-panel manuscript figure. Also generates a LaTeX template.

 Usage:
   python 02_composite_figure.py --config config.yaml
==============================================================================
"""

import argparse
import os
import sys
from pathlib import Path

import yaml
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.image as mpimg


def load_config(path):
    with open(path) as f:
        return yaml.safe_load(f)


def pdf_to_image(pdf_path, dpi=200):
    """Convert PDF first page to numpy image array."""
    # Try pdf2image (poppler)
    try:
        from pdf2image import convert_from_path
        import numpy as np
        images = convert_from_path(str(pdf_path), dpi=dpi, first_page=1, last_page=1)
        return np.array(images[0])
    except ImportError:
        pass
    # Try PyMuPDF
    try:
        import fitz
        import numpy as np
        doc = fitz.open(str(pdf_path))
        page = doc[0]
        pix = page.get_pixmap(matrix=fitz.Matrix(dpi/72, dpi/72))
        img = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.h, pix.w, pix.n)
        doc.close()
        return img
    except ImportError:
        pass
    # Try PNG fallback
    png = str(pdf_path).replace(".pdf", ".png")
    if os.path.exists(png):
        return mpimg.imread(png)
    return None


def main():
    parser = argparse.ArgumentParser(description="Assemble composite manuscript figure")
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    cfg = load_config(args.config)
    base_dir = Path(cfg["output_dir"])
    fig_dir = base_dir / "01_figures"
    out_dir = base_dir / "02_composite"
    out_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 70)
    print("STEP 2: COMPOSITE FIGURE ASSEMBLY")
    print("=" * 70)

    panel_a = fig_dir / "panel_A_synteny.pdf"
    panel_b = fig_dir / "panel_B_presence_absence.pdf"

    # ── Try rendering composite ──
    panels = {}
    for label, path in [("A", panel_a), ("B", panel_b)]:
        if path.exists():
            img = pdf_to_image(path)
            if img is not None:
                panels[label] = img
                print(f"  Panel {label}: loaded from {path.name}")
            else:
                print(f"  Panel {label}: exists but no renderer (install pdf2image or PyMuPDF)")
        else:
            print(f"  Panel {label}: NOT FOUND at {path}")

    if panels:
        fig = plt.figure(figsize=(16, 14))
        gs = GridSpec(2, 1, figure=fig, height_ratios=[2, 1], hspace=0.1)

        for i, (label, gs_slot) in enumerate([("A", gs[0]), ("B", gs[1])]):
            if label in panels:
                ax = fig.add_subplot(gs_slot)
                ax.imshow(panels[label], aspect="auto")
                ax.axis("off")
                ax.text(-0.02, 1.02, label, transform=ax.transAxes,
                       fontsize=18, fontweight="bold", va="bottom", ha="right")

        composite_path = out_dir / "Figure_ASCT_synteny.pdf"
        fig.savefig(str(composite_path), dpi=300, bbox_inches="tight",
                    facecolor="white")
        plt.close(fig)
        print(f"\n  Composite figure: {composite_path}")
    else:
        print("\n  Could not render composite (no PDF renderer installed)")
        print("  The individual panels are still available:")
        print(f"    {panel_a}")
        print(f"    {panel_b}")

    # ── LaTeX template ──
    latex = r"""\begin{figure}[htbp]
\centering
\begin{subfigure}[b]{\textwidth}
    \centering
    \includegraphics[width=\textwidth]{results/01_figures/panel_A_synteny.pdf}
    \caption{Synteny of the ASCT locus across species. Gene arrows represent
    individual ASCT copies colored by copy number. Gray ribbons indicate
    nucleotide-level sequence homology between phylogenetically adjacent species.
    Species lacking detectable ASCT copies are marked as absent.}
    \label{fig:synteny-a}
\end{subfigure}
\vspace{0.5cm}
\begin{subfigure}[b]{0.75\textwidth}
    \centering
    \includegraphics[width=\textwidth]{results/01_figures/panel_B_presence_absence.pdf}
    \caption{ASCT gene copy number variation mapped onto the species phylogeny.
    Colored cells indicate presence; gray cells indicate absence.}
    \label{fig:synteny-b}
\end{subfigure}
\caption{\textbf{Comparative synteny analysis reveals ASCT copy number variation
across molluscan and annelid genomes.}
(A)~Locus-level synteny with gene arrows and alignment links.
(B)~Phylogenetic distribution of ASCT gene copies.}
\label{fig:synteny}
\end{figure}
"""
    latex_path = out_dir / "figure_synteny.tex"
    with open(latex_path, "w") as f:
        f.write(latex)
    print(f"  LaTeX template: {latex_path}")

    print(f"\n{'=' * 70}")
    print("DONE")
    print(f"{'=' * 70}")
    print(f"  All outputs: {out_dir}")
    print(f"  SVGs for Inkscape/Illustrator editing: {fig_dir}")


if __name__ == "__main__":
    main()
