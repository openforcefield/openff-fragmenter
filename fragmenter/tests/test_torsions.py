""" Tests generating files for qm torsion scan """

import fragmenter.torsions as torsions
from fragmenter.chemi import smiles_to_molecule


def test_find_torsion_around_bond():

    mol = smiles_to_molecule("CCCC", add_atom_map=True)
    dihedral = torsions.find_torsion_around_bond(mol, (3, 4))
    assert dihedral == (0, 2, 3, 1)
