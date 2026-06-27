import matplotlib.pyplot as plt
from Bio import Phylo

# Use the dendropy-cleaned file
tree = Phylo.read("dendropy_clean.newick", "newick")

# Check names first
print("\nTerminal names in Bio.Phylo:")
for clade in tree.get_terminals():
    print(f"  {clade.name}")


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


# Plot
fig, ax = plt.subplots(1, 1, figsize=(12, 8))
Phylo.draw(tree, axes=ax, do_show=False, label_func=support_label)
ax.set_title("MrBayes Consensus Tree - ASCT Gene Copies", fontsize=14)
plt.tight_layout()
plt.savefig("tree_biophylo.png", dpi=300, bbox_inches='tight')
#plt.savefig("tree_biophylo.pdf", bbox_inches='tight')
plt.show()
