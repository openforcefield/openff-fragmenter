"""Test chemi module"""

from fragmenter import chemi, torsions
import pytest
import numpy as np
from .utils import using_openeye
import cmiles

@pytest.fixture
def mapped_molecule():
    return chemi.smiles_to_oemol('[H:5][C:1]([H:6])([H:7])[C:3]([H:11])([H:12])[C:4]([H:13])([H:14])[C:2]([H:8])([H:9])[H:10]')

@pytest.mark.parametrize('keep_confs, output', [(None, 1), (1, 1),  (2, 2),  (-1, 3),  (-2, 'error')])
def test_get_charges(keep_confs, output):
    mol = chemi.smiles_to_oemol('CCCCCCC')
    if output == 'error':
        with pytest.raises(ValueError):
            chemi.get_charges(mol, keep_confs=keep_confs)
    else:
        charged = chemi.get_charges(mol, keep_confs=keep_confs)
        for i, c in enumerate(charged.GetConfs()):
            i +=1
        assert i == output

def test_generate_conformers():
    mol = chemi.smiles_to_oemol('CCCCCCC')
    confs = chemi.generate_conformers(mol, max_confs=1)
    assert confs.GetMaxConfIdx() == 1

    confs = chemi.generate_conformers(mol)
    assert confs.GetMaxConfIdx() == 3

def generate_grid_conformers():

    mol = chemi.smiles_to_oemol('CCCC')
    dihedrals = [(0, 2, 3, 1), (3, 2, 0, 4)]
    intervals = [30, 90]
    grid = chemi.generate_grid_conformers(mol, dihedrals=dihedrals, intervals=intervals)
    assert grid.GetMaxConfIdx() == 48

@pytest.mark.parametrize('smiles, charge', [('CCCC', 0),
                                             ('CC(=O)([O-])', -1),
                                             ('C[N+](C)(C)[H]', +1)])
def test_get_charge(smiles, charge):
    mol = chemi.smiles_to_oemol(smiles)
    assert chemi.get_charge(mol) == charge

def test_smiles_to_oemol():
    from openeye import oechem
    mol = chemi.smiles_to_oemol('CCCC')
    assert isinstance(mol, oechem.OEMol)
    assert oechem.OEMolToSmiles(mol) == 'CCCC'
    assert mol.GetTitle() == 'butane'

    mol = chemi.smiles_to_oemol('CCCC', normalize=False)
    assert mol.GetTitle() == ''

    mol = chemi.smiles_to_oemol('CCCC', add_atom_map=True)
    assert oechem.OEMolToSmiles(mol) == '[H:5][C:1]([H:6])([H:7])[C:3]([H:11])([H:12])[C:4]([H:13])([H:14])[C:2]([H:8])([H:9])[H:10]'

def test_normalize_molecule():
    from openeye import oechem
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, 'CCCC')
    assert mol.GetTitle() == ''

    normalized_mol = chemi.normalize_molecule(mol)
    assert normalized_mol.GetTitle() == 'butane'

def test_smiles_to_smi():
    """Test writing out list of SMILES to smi file"""
    pass

def test_file_to_oemol():
    """Test read file to oemol list"""
    pass

def test_oemols_to_smiles():
    """Test write oemols to list of SMILES"""
    pass

def test_file_to_smiles_list():
    """Test write out file to list SMILES"""
    pass



#@using_openeye
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

#@using_openeye
# def test_qcschema_to_xyz():
#     smiles = 'HC(H)(C(H)(H)OH)OH'
#     mapped_smiles = '[H:5][C:1]([H:6])([C:2]([H:7])([H:8])[O:4][H:10])[O:3][H:9]'
#     mol = cmiles.utils.load_molecule(smiles)
#     mapped_mol = cmiles.utils.load_molecule(mapped_smiles)
#
#     dihedrals = [(2, 0, 1, 3), (0, 1, 3, 9), (1, 0, 2, 8)]
#     intervals = [90, 90, 90]
#     mult_conf = chemi.generate_grid_conformers(mapped_mol, dihedrals, intervals)
#
#     # generate list of qcschema molecules
#     qcschema_molecules = [cmiles.utils.mol_to_map_ordered_qcschema(conf, mol_id) for conf in multi_conf.GetConfs()]


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
    # assert mult_conf.GetMaxConfIdx() == 216

    # intervals = [90, 90, 90]
    # mult_conf = chemi.generate_grid_conformers(mapped_mol, dihedrals, intervals)
    # assert mult_conf.GetMaxConfIdx() == 64

@using_openeye
def test_remove_clashes():
    pass

@using_openeye
def test_resolve_clashes():
    pass

@using_openeye
@pytest.mark.parametrize('smiles', ['OCCO', 'C(CO)O', '[H]C([H])(C([H])([H])O[H])O[H]',
                                    '[H:5][C:1]([H:6])([C:2]([H:7])([H:8])[O:4][H:10])[O:3][H:9]',
                                   '[H][O][C]([H])([H])[C]([H])([H])[O][H]',
                                    '[O:1]([C:3]([C:4]([O:2][H:6])([H:9])[H:10])([H:7])[H:8])[H:5]'])
def canonical_order_conformer(smiles):
    """Test that geometry is ordered the same way every time no matter the SMILES used to create the molecule"""
    import cmiles
    mapped_smiles = '[H:5][C:1]([H:6])([C:2]([H:7])([H:8])[O:4][H:10])[O:3][H:9]'
    mol_id_oe = cmiles.to_molecule_id(mapped_smiles, canonicalization='openeye')
    oemol = cmiles.utils.load_molecule(mapped_smiles, toolkit='openeye')
    # Generate canonical geometry
    conf = chemi.generate_conformers(oemol, can_order=True, max_confs=1)
    mapped_symbols, mapped_geometry = cmiles._cmiles_oe.get_map_ordered_geometry(conf, mapped_smiles)
    # #mapped_symbols = ['C', 'C', 'O', 'O', 'H', 'H', 'H', 'H', 'H', 'H']
    # mapped_geometry = [-1.6887193912042044, 0.8515190939276903, 0.8344587822904272, -4.05544806361675, -0.3658269566455062,
    #                    -0.22848169646448416, -1.6111611950422127, 0.4463128276938808, 3.490617694146934, -3.97756355964586,
    #                    -3.0080934853087373, 0.25948499322223956, -1.6821252026076652, 2.891135395246369, 0.4936556190978574,
    #                    0.0, 0.0, 0.0, -4.180315034973438, -0.09210893239246959, -2.2748227320305525, -5.740516456782416,
    #                    0.4115539217904015, 0.6823267491485907, -0.07872657410528058, 1.2476492272884379, 4.101615944163073,
    #                    -5.514569080545831, -3.7195945404657222, -0.4441653010509862]

    mol = cmiles.utils.load_molecule(smiles, toolkit='openeye')
    # if not cmiles.utils.has_explicit_hydrogen(mol):
    #     mol = utils.add_explicit_hydrogen(mol)
    atom_map = cmiles.utils.get_atom_map(mol, mapped_smiles=mapped_smiles)
    # use the atom map to add coordinates to molecule. First reorder mapped geometry to order in molecule
    mapped_coords = np.array(mapped_geometry, dtype=float).reshape(int(len(mapped_geometry)/3), 3)
    coords = np.zeros((mapped_coords.shape))
    for m in atom_map:
        coords[atom_map[m]] = mapped_coords[m-1]
    # flatten
    coords = coords.flatten()
    # convert to Angstroms
    coords = coords*cmiles.utils.BOHR_2_ANGSTROM
    # set coordinates in oemol
    mol.SetCoords(coords)
    mol.SetDimension(3)

    # Get new atom map
    atom_map = cmiles.utils.get_atom_map(mol, mapped_smiles)
    symbols, geometry = cmiles._cmiles_oe.get_map_ordered_geometry(mol, mapped_smiles)
    assert geometry == mapped_geometry
    assert symbols == mapped_symbols
