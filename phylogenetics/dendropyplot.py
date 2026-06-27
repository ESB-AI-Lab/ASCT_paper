import dendropy
import re

# Read the original nexus file (dendropy handles it properly)
tns = dendropy.TaxonNamespace()
tree = dendropy.Tree.get(
    path="output.best.nex.con.tre",
    schema="nexus",
    taxon_namespace=tns,
    preserve_underscores=True
)

# Check the names are correct
print("Taxon names:")
for taxon in tree.taxon_namespace:
    print(f"  {taxon.label}")

# Extract posterior probabilities and set as node labels
for node in tree.preorder_node_iter():
    if not node.is_leaf():
        try:
            prob_val = float(node.annotations.get_value('prob'))
            node.label = f"{prob_val:.2f}"
        except (TypeError, ValueError):
            pass

# Save as clean newick with support values preserved as node labels
tree.write_to_path("dendropy_clean.newick", schema="newick",
                   suppress_annotations=True)

# Strip any remaining bracket annotations (e.g., [&U])
with open("dendropy_clean.newick", "r") as f:
    nwk = f.read()
nwk_clean = re.sub(r'\[&[^\]]*\]', '', nwk)
with open("dendropy_clean.newick", "w") as f:
    f.write(nwk_clean)
