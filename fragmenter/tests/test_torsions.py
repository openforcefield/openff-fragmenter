""" Tests generating files for qm torsion scan """

import fragmenter.torsions as torsions
from fragmenter import chemi


def test_find_torsion_around_bond():
    mol = chemi.smiles_to_oemol('CCCC', add_atom_map=True)
    dihedral = torsions.find_torsion_around_bond(mol, (3, 4))
    assert dihedral == [0, 2, 3, 1]
