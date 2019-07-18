from fragmenter import fragment, chemi
import json

"""
This is an example script to generate fragments from a molecule
Note: The current fragmentation scheme is not optimal yet so it can be slow
for larger molecules.
"""

# Create an oemol from a SMILES
oemol = chemi.smiles_to_oemol('OC1(CN(C1)C(=O)C1=C(NC2=C(F)C=C(I)C=C2)C(F)=C(F)C=C1)[C@@H]1CCCCN1',
                              name='Cobimetinib')

# Instantiate a fragmenter engine
frag_engine = fragment.WBOFragmenter(oemol)
# Use default options to fragment molecule.
frag_engine.fragment()
# Generate PDF with fragments
frag_engine.depict_fragments(fname='example_fragments.pdf')

# Generate input for torsiondrive jobs. These jobs will drive the central rotatable bond in the fragment
td_inputs = frag_engine.to_torsiondrive_json(max_confs=10)
with open('example_td_inputs.json', 'w') as f:
    json.dump(td_inputs, f, indent=2, sort_keys=True)

# If you want to only generate starting conformations for other calculations, use `to_qcschema_mols`
qcschema_mols = frag_engine.to_qcschema_mols(max_confs=10)
with open('example_qcschema_mols.json', 'w') as f:
    json.dump(qcschema_mols, f, indent=2, sort_keys=True)