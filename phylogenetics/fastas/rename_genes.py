#!/usr/bin/env python3
"""
Rename ASCT gene copies: copy names → Mgal-anchored → paper names.
Already applied — kept as documentation of the full rename chain.
"""

import shutil
from pathlib import Path

WORKDIR = Path(__file__).parent

# Phase 1: copy names → Mgal-anchored orthology (already applied)
COPY_TO_ORTHOLOGY = {
    "Asen_copy1": "Asen_ASCT2",  "Asen_copy2": "Asen_ASCT1",  "Asen_copy3": "Asen_ASCT3",
    "Mcal_copy1": "Mcal_ASCT2",  "Mcal_copy2": "Mcal_ASCT1",  "Mcal_copy3": "Mcal_ASCT3",
    "Mcor_copy1": "Mcor_ASCT1",  "Mcor_copy2": "Mcor_ASCT3a", "Mcor_copy3": "Mcor_ASCT3b",
    "Medi_copy1": "Medi_ASCT2a", "Medi_copy2": "Medi_ASCT1b", "Medi_copy3": "Medi_ASCT1a",
    "Medi_copy4": "Medi_ASCT3",  "Medi_copy5": "Medi_ASCT2b",
    "Mgal_copy1": "Mgal_ASCT1",  "Mgal_copy2": "Mgal_ASCT2",  "Mgal_copy3": "Mgal_ASCT3",
    "Pvir_copy1": "Pvir_ASCT1",  "Pvir_copy2": "Pvir_ASCT2",  "Pvir_copy3": "Pvir_ASCT3",
}

# Phase 2: Mgal-anchored → paper names (already applied)
ORTHOLOGY_TO_PAPER = {
    "Asen_ASCT1":  "Asen_ASCT2b2",  "Asen_ASCT2":  "Asen_ASCT2a1",  "Asen_ASCT3":  "Asen_ASCT1.3",
    "Mcal_ASCT1":  "Mcal_ASCT2b1",  "Mcal_ASCT2":  "Mcal_ASCT2a2",  "Mcal_ASCT3":  "Mcal_ASCT1.3",
    "Mcor_ASCT1":  "Mcor_ASCT2b1",  "Mcor_ASCT3a": "Mcor_ASCT1.3a", "Mcor_ASCT3b": "Mcor_ASCT1.3b",
    "Medi_ASCT1a": "Medi_ASCT2b3",  "Medi_ASCT1b": "Medi_ASCT2b2",  "Medi_ASCT2a": "Medi_ASCT2a1",
    "Medi_ASCT2b": "Medi_ASCT2a5",  "Medi_ASCT3":  "Medi_ASCT1.4",
    "Mgal_ASCT1":  "Mgal_ASCT2b1",  "Mgal_ASCT2":  "Mgal_ASCT2a2",  "Mgal_ASCT3":  "Mgal_ASCT1.3",
    "Pvir_ASCT1":  "Pvir_ASCT2b1",  "Pvir_ASCT2":  "Pvir_ASCT2a2",  "Pvir_ASCT3":  "Pvir_ASCT1.3",
}

# Current active map: copy names → paper names (composite)
RENAME_MAP = {k: ORTHOLOGY_TO_PAPER[v] for k, v in COPY_TO_ORTHOLOGY.items()}


def rename_fasta_files():
    """Rename FASTA files and update headers."""
    for old_name, new_name in RENAME_MAP.items():
        old_path = WORKDIR / f"{old_name}_gene.fasta"
        new_path = WORKDIR / f"{new_name}_gene.fasta"

        if not old_path.exists():
            print(f"  SKIP: {old_path.name} not found")
            continue

        with open(old_path) as f:
            content = f.read()

        content = content.replace(f">{old_name}", f">{new_name}", 1)

        with open(new_path, "w") as f:
            f.write(content)

        old_path.unlink()
        print(f"  {old_name}_gene.fasta -> {new_name}_gene.fasta")


def update_bed_files():
    """Update the name column in BED files."""
    for bed_path in sorted(WORKDIR.glob("*_ASCT_mRNA.bed")):
        lines = bed_path.read_text().splitlines()
        updated = []
        changed = False
        for line in lines:
            if line.startswith("#"):
                updated.append(line)
                continue
            cols = line.split("\t")
            old_name = cols[3]
            if old_name in RENAME_MAP:
                cols[3] = RENAME_MAP[old_name]
                changed = True
            updated.append("\t".join(cols))
        if changed:
            bed_path.write_text("\n".join(updated) + "\n")
            print(f"  Updated {bed_path.name}")


def main():
    print("=== Renaming FASTA files ===")
    rename_fasta_files()

    print("\n=== Updating BED files ===")
    update_bed_files()

    print("\n=== Rename complete ===")
    print("\nNew files:")
    for f in sorted(WORKDIR.glob("*_ASCT*_gene.fasta")):
        print(f"  {f.name}")


if __name__ == "__main__":
    main()
