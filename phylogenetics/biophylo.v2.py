import matplotlib.pyplot as plt
from Bio import Phylo

# Create publication-quality tree
tree = Phylo.read("dendropy_clean.newick", "newick")


# Store confidence values then clear them so Bio.Phylo doesn't draw them separately
saved_conf = {}
for clade in tree.find_clades():
    if not clade.is_terminal() and clade.confidence is not None:
        saved_conf[id(clade)] = clade.confidence
        clade.confidence = None


def support_label(clade):
    if clade.is_terminal():
        return clade.name
    conf = saved_conf.get(id(clade))
    if conf is not None:
        pct = conf * 100 if conf <= 1.0 else conf
        if pct >= 95.0:
            return "*"
        if pct < 50.0:
            return "-"
        return f"{pct:.0f}"
    return ""


fig, ax = plt.subplots(1, 1, figsize=(12, 8))
Phylo.draw(tree, axes=ax, do_show=False, label_func=support_label)

# Professional styling
ax.set_title('Phylogenetic Analysis of ASCT Gene Copies\nacross Mytilus Species Complex',
            fontsize=16, fontweight='bold', pad=20)

for spine in ax.spines.values():
    spine.set_visible(False)
ax.set_xticks([])
ax.set_yticks([])

ax.text(0.98, 0.02, 'Branch lengths = substitutions per site\nMrBayes consensus tree (50% majority rule)\n* = posterior probability ≥ 0.95',
        transform=ax.transAxes, fontsize=9, style='italic',
        ha='right', va='bottom')

plt.tight_layout()
plt.savefig("ASCT_phylogeny.png", dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig("ASCT_phylogeny.pdf", bbox_inches='tight')
plt.show()
