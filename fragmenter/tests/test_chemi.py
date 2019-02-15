"""Test chemi module"""

from fragmenter import chemi
from fragmenter.tests.utils import has_openeye
import pytest
import numpy as np
from .utils import using_openeye
import cmiles

@pytest.fixture
def mapped_molecule():
    return chemi.smiles_to_oemol('[H:5][C:1]([H:6])([H:7])[C:3]([H:11])([H:12])[C:4]([H:13])([H:14])[C:2]([H:8])([H:9])[H:10]')


def test_connectivity_table(mapped_molecule):
    """Test generate connectivity table"""
    expected_table = np.array([[0, 2, 1],
                        [1, 3, 1],
                        [2, 3, 1],
                        [0, 4, 1],
                        [0, 5, 1],
                        [0, 6, 1],
                        [1, 7, 1],
                        [1, 8, 1],
                        [1, 9, 1],
                        [2, 10, 1],
                        [2, 11, 1],
                        [3, 12, 1],
                        [3, 13, 1]])
    connectivity_table = chemi.get_mapped_connectivity_table(mapped_molecule)

    for bond in connectivity_table:
        xi = np.isin(expected_table, bond[:2])
        match = np.where(np.array([i[:2].sum() for i in xi]) == 2)[0]
        # assert that a match was found and only one was found
        assert len(match) == 1
        # assert that bond order is the same
        assert expected_table[match][0][-1] == bond[-1]

@using_openeye
def test_to_mapped_xyz():
    from openeye import oechem
    smiles = 'HC(H)(C(H)(H)OH)OH'
    mapped_smiles = '[H:5][C:1]([H:6])([C:2]([H:7])([H:8])[O:4][H:10])[O:3][H:9]'
    mol = cmiles.utils.load_molecule(smiles)
    mapped_mol = cmiles.utils.load_molecule(mapped_smiles)

    with pytest.raises(ValueError):
        chemi.to_mapped_xyz(mapped_mol)
    # generate conformer
    mol = chemi.generate_conformers(mol, max_confs=1)
    mapped_mol = chemi.generate_conformers(mapped_mol, max_confs=1)
    atom_map = cmiles.utils.get_atom_map(mol, mapped_smiles)

    with pytest.raises(ValueError):
        chemi.to_mapped_xyz(mol)

    xyz_1 = chemi.to_mapped_xyz(mol, atom_map)
    xyz_2 = chemi.to_mapped_xyz(mapped_mol)
    assert xyz_1 == xyz_2

    xyz_1 = sorted(xyz_1.split('\n')[2:-1])
    xyz_2 = sorted(xyz_2.split('\n')[2:-1])
    assert xyz_1 == xyz_2

@using_openeye
def teat_qcschema_to_xyz():
    smiles = 'HC(H)(C(H)(H)OH)OH'
    mapped_smiles = '[H:5][C:1]([H:6])([C:2]([H:7])([H:8])[O:4][H:10])[O:3][H:9]'
    mol = cmiles.utils.load_molecule(smiles)
    mapped_mol = cmiles.utils.load_molecule(mapped_smiles)

    dihedrals = [(2, 0, 1, 3), (0, 1, 3, 9), (1, 0, 2, 8)]
    intervals = [90, 90, 90]
    mult_conf = chemi.generate_grid_conformers(mapped_mol, dihedrals, intervals)

    # generate list of qcschema molecules
    qcschema_molecules = [cmiles.utils.mol_to_map_ordered_qcschema(conf, mol_id) for conf in multi_conf.GetConfs()]


@using_openeye
def test_grid_multi_conformers():
    "Test generating grid multiconformer"
    smiles = 'HC(H)(C(H)(H)OH)OH'
    mapped_smiles = '[H:5][C:1]([H:6])([C:2]([H:7])([H:8])[O:4][H:10])[O:3][H:9]'
    mol = cmiles.utils.load_molecule(smiles)
    mapped_mol = cmiles.utils.load_molecule(mapped_smiles)

    dihedrals = [(2, 0, 1, 3), (0, 1, 3, 9), (1, 0, 2, 8)]
    intervals = [60, 60, 60]
    with pytest.raises(ValueError):
        chemi.generate_grid_conformers(mol, dihedrals, intervals)

    mult_conf = chemi.generate_grid_conformers(mapped_mol, dihedrals, intervals)
    assert mult_conf.GetMaxConfIdx() == 216

    intervals = [90, 90, 90]
    mult_conf = chemi.generate_grid_conformers(mapped_mol, dihedrals, intervals)
    assert mult_conf.GetMaxConfIdx() == 64

@using_openeye
def test_remove_atom_map():
    from openeye import oechem
    mapped_smiles = '[H:5][C:1]([H:6])([C:2]([H:7])([H:8])[O:4][H:10])[O:3][H:9]'
    mapped_mol = oechem.OEMol()
    oechem.OESmilesToMol(mapped_mol, mapped_smiles)

    chemi.remove_map(mapped_mol)
    assert oechem.OEMolToSmiles(mapped_mol) == 'C(CO)O'

    chemi.restore_map(mapped_mol)
    assert oechem.OEMolToSmiles(mapped_mol) == mapped_smiles

@using_openeye
def test_remove_clashes():
    pass

@using_openeye
def test_resolve_clashes():
    pass