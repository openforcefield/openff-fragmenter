"""
Unit and regression test for the fragmenter package.
"""

# Import package, test suite, and other packages as needed
import fragmenter
from fragmenter import chemi
from cmiles.utils import mol_to_smiles, remove_atom_map
import sys
import unittest
from fragmenter.tests.utils import get_fn, has_openeye, using_openeye

import pytest


def test_fragmenter_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "fragmenter" in sys.modules


class TestFragment(unittest.TestCase):

    @unittest.skipUnless(has_openeye, 'Cannot test without OpenEye')
    def test_stereo_parent(self):
        pass

    @unittest.skipUnless(has_openeye, 'Cannot test without OpenEye')
    def test_expand_protonation_states(self):
        """Test expand protonation states"""
        smiles = 'C5=C(C1=CN=CC=C1)N=C(NC2=C(C=CC(=C2)NC(C3=CC=C(C=C3)CN4CCN(CC4)C)=O)C)N=C5'
        molecule = chemi.smiles_to_oemol(smiles)
        protonation = fragmenter.fragment._expand_states(molecule)
        protonation_1 = {'Cc1ccc(cc1Nc2nccc(n2)c3ccc[nH+]c3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C',
                         'Cc1ccc(cc1Nc2nccc(n2)c3ccc[nH+]c3)NC(=O)c4ccc(cc4)CN5CC[NH+](CC5)C',
                         'Cc1ccc(cc1Nc2nccc(n2)c3ccc[nH+]c3)NC(=O)c4ccc(cc4)C[NH+]5CCN(CC5)C',
                         'Cc1ccc(cc1Nc2nccc(n2)c3ccc[nH+]c3)NC(=O)c4ccc(cc4)C[NH+]5CC[NH+](CC5)C',
                         'Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C',
                         'Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CC[NH+](CC5)C',
                         'Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)C[NH+]5CCN(CC5)C',
                         'Cc1ccc(cc1Nc2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)C[NH+]5CC[NH+](CC5)C',
                         'Cc1ccc(cc1[N-]c2nccc(n2)c3ccc[nH+]c3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C',
                         'Cc1ccc(cc1[N-]c2nccc(n2)c3ccc[nH+]c3)NC(=O)c4ccc(cc4)CN5CC[NH+](CC5)C',
                         'Cc1ccc(cc1[N-]c2nccc(n2)c3ccc[nH+]c3)NC(=O)c4ccc(cc4)C[NH+]5CCN(CC5)C',
                         'Cc1ccc(cc1[N-]c2nccc(n2)c3ccc[nH+]c3)NC(=O)c4ccc(cc4)C[NH+]5CC[NH+](CC5)C',
                         'Cc1ccc(cc1[N-]c2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CCN(CC5)C',
                         'Cc1ccc(cc1[N-]c2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)CN5CC[NH+](CC5)C',
                         'Cc1ccc(cc1[N-]c2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)C[NH+]5CCN(CC5)C',
                         'Cc1ccc(cc1[N-]c2nccc(n2)c3cccnc3)NC(=O)c4ccc(cc4)C[NH+]5CC[NH+](CC5)C'}
        protonation_2 = set()
        for mol in protonation:
            protonation_2.add(mol_to_smiles(mol, mapped=False, explicit_hydrogen=False, isomeric=True))

        intersection = protonation_1.intersection(protonation_2)
        self.assertEqual(len(intersection), len(protonation_1))
        self.assertEqual(len(intersection), len(protonation_2))


    @unittest.skipUnless(has_openeye, 'Cannot test without OpenEye')
    def test_expand_tautomers(self):
        """Test expand tautomer"""
        smiles_1 ='c1ccc2c(c1)C=CCC2=O'
        smiles_2 = 'c1ccc2c(c1)cccc2O'
        molecule_1 = chemi.smiles_to_oemol(smiles_1)
        molecule_2 = chemi.smiles_to_oemol(smiles_2)
        tautomers_1 = fragmenter.fragment.expand_states(molecule_1, protonation=False, tautomers=True, stereoisomers=False)
        tautomers_2 = fragmenter.fragment.expand_states(molecule_2, protonation=False, tautomers=True, stereoisomers=False)

        self.assertEqual(tautomers_1, tautomers_2)


    @unittest.skipUnless(has_openeye, 'Cannot test without OpenEye')
    def test_expand_enantiomers(self):
        smiles = 'CN(C)C/C=C/C(=O)NC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC(=C(C=C3)F)Cl)O[C@H]4CCOC4'
        molecule = chemi.smiles_to_oemol(smiles)
        stereoisomers = fragmenter.fragment._expand_states(molecule, enumerate='stereoisomers')

        stereoisomers_1 = {'CN(C)C/C=C/C(=O)Nc1cc2c(cc1O[C@@H]3CCOC3)ncnc2Nc4ccc(c(c4)Cl)F',
                        'CN(C)C/C=C/C(=O)Nc1cc2c(cc1O[C@H]3CCOC3)ncnc2Nc4ccc(c(c4)Cl)F',
                        'CN(C)C/C=C\\C(=O)Nc1cc2c(cc1O[C@@H]3CCOC3)ncnc2Nc4ccc(c(c4)Cl)F',
                        'CN(C)C/C=C\\C(=O)Nc1cc2c(cc1O[C@H]3CCOC3)ncnc2Nc4ccc(c(c4)Cl)F'}

        stereoisomers_2 = set()
        for mol in stereoisomers:
            stereoisomers_2.add(mol_to_smiles(mol, mapped=False, explicit_hydrogen=False, isomeric=True))
        intersection = stereoisomers_1.intersection(stereoisomers_2)
        self.assertEqual(len(intersection), len(stereoisomers_1))
        self.assertEqual(len(intersection), len(stereoisomers_2))
        self.assertEqual(len(stereoisomers_1), len(stereoisomers_2))

@using_openeye
def test_keep_track_of_map():
    from openeye import oechem
    mapped_smiles = '[H:45][c:8]1[cH:18][c:10]([c:19]([cH:17][c:7]1[H:44])[N:35]([H:67])[c:21]2[n:32][cH:20][c:9]([c:12]([n:31]2)[H:49])[H:46])[H:47]'
    mapped_mol = oechem.OEMol()
    oechem.OESmilesToMol(mapped_mol, mapped_smiles)

    frags = fragmenter.fragment.CombinatorialFragmenter(mapped_mol)
    frags.fragment()
    #frags.fragment_all_bonds_not_in_ring_systems()
    #frags.combine_fragments(min_rotors=1, max_rotors=frags.n_rotors+1, restore_maps=True)

    keys = list(frags.fragments.keys())
    assert oechem.OEMolToSmiles(frags.fragments[keys[0]][0]) == '[H:45][c:8]1[cH:18][c:10]([c:19]([cH:17][c:7]1[H:44])[NH:35][H:67])[H:47]'
    assert oechem.OEMolToSmiles(frags.fragments[keys[1]][0]) == '[H:46][c:9]1[cH:20][n:32][c:21]([n:31][c:12]1[H:49])[NH:35][H:67]'


@using_openeye
def test_tag_fgroups():
    from openeye import oechem
    import itertools
    smiles = '[H:40][c:3]1[c:8]([c:20]2[n:30][c:12]([c:14]([n:32]2[n:31][c:11]1[H:48])[C:2]#[C:1][c:13]3[c:9]([c:15]([c:4]([c:5]([c:16]3[C:26]' \
             '([H:58])([H:59])[H:60])[H:42])[H:41])[C:21](=[O:36])[N:35]([H:66])[c:19]4[c:7]([c:6]([c:17]([c:18]([c:10]4[H:47])[C:29]([F:37])([F:38])' \
             '[F:39])[C:28]([H:64])([H:65])[N:34]5[C:24]([C:22]([N:33]([C:23]([C:25]5([H:56])[H:57])([H:52])[H:53])[C:27]([H:61])([H:62])[H:63])([H:50])' \
             '[H:51])([H:54])[H:55])[H:43])[H:44])[H:46])[H:49])[H:45]'
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    frags = fragmenter.fragment.CombinatorialFragmenter(mol)
    fgroups = {}
    fgroups['alkyne_0'] = [1, 2]
    fgroups['carbonyl_0'] = [21, 36]
    fgroups['amide_0'] = [35]
    fgroups['tri_halide_0'] = [29, 37, 38, 39]
    for group in fgroups:
        for i in fgroups[group]:
            a = frags.molecule.GetAtom(oechem.OEHasMapIdx(i))
            assert a.GetData('fgroup') == group
    for group in fgroups:
        atoms = [frags.molecule.GetAtom(oechem.OEHasMapIdx(i)) for i in fgroups[group]]
        for atom in itertools.combinations(atoms, 2):
            # Check for bond
            b = frags.molecule.GetBond(atom[0], atom[1])
            if b:
                assert b.GetData('fgroup') == group

@using_openeye
def test_rotor_wbo():
    from openeye import oechem
    smiles ='[H:5][C:1]([H:6])([H:7])[C:3]([H:11])([H:12])[C:4]([H:13])([H:14])[C:2]([H:8])([H:9])[H:10]'
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    f = fragmenter.fragment.WBOFragmenter(mol)
    assert f.rotors_wbo == {}
    f._get_rotor_wbo()
    assert list(f.rotors_wbo.keys()) == [(3, 4)]
    assert round(f.rotors_wbo[(3, 4)], ndigits=3) ==  0.986


@using_openeye
def test_get_bond():
    from openeye import oechem
    smiles ='[H:5][C:1]([H:6])([H:7])[C:3]([H:11])([H:12])[C:4]([H:13])([H:14])[C:2]([H:8])([H:9])[H:10]'
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    f = fragmenter.fragment.WBOFragmenter(mol)
    bond = f.get_bond(bond_tuple=(3, 4))
    assert bond.IsRotor()

def test_build_fragment():
    from openeye import oechem
    smiles = 'CCCCCC'
    mol = chemi.smiles_to_oemol(smiles)
    f = fragmenter.fragment.WBOFragmenter(mol)
    f.calculate_wbo()
    f._get_rotor_wbo()
    for bond in f.rotors_wbo:
        f.build_fragment(bond)
    assert len(f.fragments) == 3
    for bond in f.fragments:
        remove_atom_map(f.fragments[bond])
        assert oechem.OEMolToSmiles(f.fragments[bond]) == 'CCCC'


@using_openeye
def test_build_WBOfragment():
    """ Test build fragment"""
    from openeye import oechem
    smiles = 'CCCCC'
    mol = chemi.smiles_to_oemol(smiles)
    oechem.OESmilesToMol(mol, smiles)
    f = fragmenter.fragment.WBOFragmenter(mol)
    f.fragment()
    assert len(f.fragments) == len(f.rotors_wbo)
    assert f.fragments.keys() == f.rotors_wbo.keys()


def test_to_atom_bond_set():
    from openeye import oechem
    smiles = '[H:38][c:1]1[c:2]([c:14]([n:28][c:5]([c:8]1[C:25]([H:64])([H:65])[N:33]2[C:17]([C:19]([N:34]([C:20]([C:18]2([H:46])[H:47])([H:50])[H:51])[C:26]([H:66])([H:67])[C:22]([H:55])([H:56])[H:57])([H:48])[H:49])([H:44])[H:45])[H:42])[N:35]([H:69])[c:16]3[n:29][c:6]([c:12]([c:13]([n:31]3)[c:7]4[c:3]([c:10]5[c:9]([c:11]([c:4]4[H:41])[F:36])[n:30][c:15]([n:32]5[C:27]([H:68])([C:23]([H:58])([H:59])[H:60])[C:24]([H:61])([H:62])[H:63])[C:21]([H:52])([H:53])[H:54])[H:40])[F:37])[H:43])[H:39]'
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    f = fragmenter.fragment.Fragmenter(mol)
    atoms = {17, 18, 19, 20, 22, 26, 33, 34, 66, 67}
    bonds = {(17, 19), (17, 33), (18, 20), (18, 33), (19, 34), (20, 34), (22, 26), (26, 34),  (26, 66), (26, 67)}
    atom_bond_set = f._to_atom_bond_set(atoms=atoms, bonds=bonds)
    atoms_2 = set([a.GetMapIdx() for a in atom_bond_set.GetAtoms()])
    assert atoms == atoms_2
    bonds_2 = set()
    for b in atom_bond_set.GetBonds():
        a1 = b.GetBgn().GetMapIdx()
        a2 = b.GetEnd().GetMapIdx()
        bonds_2.add(tuple(sorted((a1, a2))))
    assert bonds == bonds_2

def test_atom_bond_set_to_mol():
    from openeye import oechem
    smiles = '[H:38][c:1]1[c:2]([c:14]([n:28][c:5]([c:8]1[C:25]([H:64])([H:65])[N:33]2[C:17]([C:19]([N:34]([C:20]([C:18]2([H:46])[H:47])([H:50])[H:51])[C:26]([H:66])([H:67])[C:22]([H:55])([H:56])[H:57])([H:48])[H:49])([H:44])[H:45])[H:42])[N:35]([H:69])[c:16]3[n:29][c:6]([c:12]([c:13]([n:31]3)[c:7]4[c:3]([c:10]5[c:9]([c:11]([c:4]4[H:41])[F:36])[n:30][c:15]([n:32]5[C:27]([H:68])([C:23]([H:58])([H:59])[H:60])[C:24]([H:61])([H:62])[H:63])[C:21]([H:52])([H:53])[H:54])[H:40])[F:37])[H:43])[H:39]'
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    f = fragmenter.fragment.Fragmenter(mol)
    atoms = {17, 18, 19, 20, 22, 26, 33, 34, 66, 67}
    bonds = {(17, 19), (17, 33), (18, 20), (18, 33), (19, 34), (20, 34), (22, 26), (26, 34),  (26, 66), (26, 67)}
    atom_bond_set = f._to_atom_bond_set(atoms=atoms, bonds=bonds)
    mol = f.atom_bond_set_to_mol(atom_bond_set)
    for b in mol.GetBonds():
        a1 = b.GetBgn()
        a2 = b.GetEnd()
        if not a1.IsHydrogen() and not a2.IsHydrogen():
            assert tuple(sorted((a1.GetMapIdx(), a2.GetMapIdx()))) in bonds

def test_calculate_wbo():
    smiles = 'CCCC'
    oemol = chemi.smiles_to_oemol(smiles, name='butane')
    f = fragmenter.fragment.WBOFragmenter(oemol)
    mol = f.calculate_wbo()
    assert not mol
    for bond in f.molecule.GetBonds():
        assert 'WibergBondOrder' in bond.GetData()

    mol = f.calculate_wbo(f.molecule)
    assert mol
    for bond in mol.GetBonds():
        assert 'WibergBondOrder' in bond.GetData()

def test_compare_wbo():
    from openeye import oechem
    smiles ='[H:5][C:1]([H:6])([H:7])[C:3]([H:11])([H:12])[C:4]([H:13])([H:14])[C:2]([H:8])([H:9])[H:10]'
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    f = fragmenter.fragment.WBOFragmenter(mol)
    f.calculate_wbo()
    f._get_rotor_wbo()
    assert f.compare_wbo(fragment=mol, bond_tuple=(3, 4)) == 0.0

@pytest.mark.parametrize('input, output', [('CCCC', 0),
                                           ('c1ccccc1', 1),
                                           ('c1ccccc1C', 1)])

def test_find_ring_systems(input, output):
    mol = chemi.smiles_to_oemol(input)
    f = fragmenter.fragment.WBOFragmenter(mol)
    f._find_ring_systems()
    assert len(f.ring_systems) == output

@pytest.mark.parametrize('input, output', [(True, 7),
                                           (False, 6)])
def test_keep_non_rotor(input, output):
    mol = chemi.smiles_to_oemol('c1ccccc1C')
    f = fragmenter.fragment.WBOFragmenter(mol)
    f._find_ring_systems(keep_non_rotor_ring_substituents=input)
    assert len(f.ring_systems[1][0]) == output

def test_find_ortho_substituent():
    from openeye import oechem
    smiles ="[H:34][c:1]1[c:2]([c:6]([c:7]([c:8]([c:3]1[H:36])[Cl:33])[N:28]([H:57])[C:14](=[O:30])[c:9]2[c:5]([n:23][c:13]([s:32]2)[N:29]([H:58])[c:11]3[c:4]([c:10]([n:24][c:12]([n:25]3)[C:20]([H:50])([H:51])[H:52])[N:26]4[C:15]([C:17]([N:27]([C:18]([C:16]4([H:41])[H:42])([H:45])[H:46])[C:21]([H:53])([H:54])[C:22]([H:55])([H:56])[O:31][H:59])([H:43])[H:44])([H:39])[H:40])[H:37])[H:38])[C:19]([H:47])([H:48])[H:49])[H:35]"
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    f = fragmenter.fragment.WBOFragmenter(mol)
    f._find_ring_systems(keep_non_rotor_ring_substituents=False)
    ortho = f._find_ortho_substituent(ring_idx=1, rot_bond=(7, 28))
    assert len(ortho[0]) == 2
    assert len(ortho[1]) == 2
    assert ortho[0] == set((19, 33))
    assert ortho[1] == set(((19, 6), (33, 8)))

def test_find_rotatable_bonds():
    from openeye import oechem
    smiles = '[H:38][c:1]1[c:2]([c:14]([n:28][c:5]([c:8]1[C:25]([H:64])([H:65])[N:33]2[C:17]([C:19]([N:34]([C:20]([C:18]2([H:46])[H:47])([H:50])[H:51])[C:26]([H:66])([H:67])[C:22]([H:55])([H:56])[H:57])([H:48])[H:49])([H:44])[H:45])[H:42])[N:35]([H:69])[c:16]3[n:29][c:6]([c:12]([c:13]([n:31]3)[c:7]4[c:3]([c:10]5[c:9]([c:11]([c:4]4[H:41])[F:36])[n:30][c:15]([n:32]5[C:27]([H:68])([C:23]([H:58])([H:59])[H:60])[C:24]([H:61])([H:62])[H:63])[C:21]([H:52])([H:53])[H:54])[H:40])[F:37])[H:43])[H:39]'
    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, smiles)
    # Add map
    f = fragmenter.fragment.WBOFragmenter(mol)
    rot_bonds = f._find_rotatable_bonds()
    assert len(rot_bonds) == 7
    expected_rot_bonds = [(14, 35), (8, 25), (25, 33), (34, 26), (35, 16), (13, 7), (32, 27)]
    for bond in rot_bonds:
        assert bond in expected_rot_bonds

def test_add_substituent():
    smiles = 'CCCCCC'
    mol = chemi.smiles_to_oemol(smiles)
    f = fragmenter.fragment.WBOFragmenter(mol)
    f.fragment()
    for bond in f.fragments:
        assert mol_to_smiles(f.fragments[bond], mapped=False, explicit_hydrogen=False) == 'CCCC'

    mol = f.fragments[(3, 5)]
    atoms = set()
    bonds = set()
    for a in mol.GetAtoms():
        if a.IsHydrogen():
            continue
        atoms.add(a.GetMapIdx())
    for b in mol.GetBonds():
        a1 = b.GetBgn()
        a2 = b.GetEnd()
        if a1.IsHydrogen() or a2.IsHydrogen():
            continue
        bonds.add((a1.GetMapIdx(), a2.GetMapIdx()))

    mol = f._add_next_substituent(atoms, bonds, target_bond=(3, 5))

    assert mol_to_smiles(mol, mapped=False, explicit_hydrogen=False) == 'CCCCC'

def test_to_json():
    smiles = 'CCCCCC'
    mol = chemi.smiles_to_oemol(smiles)
    f = fragmenter.fragment.WBOFragmenter(mol)
    f.fragment()
    json_dict = f.to_json()
    assert len(json_dict) == 1
    assert list(json_dict.keys())[0] == 'CCCC'
    assert 'provenance' in json_dict['CCCC']
    assert 'cmiles_identifiers' in json_dict['CCCC']

def test_to_qcscheme_mol():
    smiles = 'CCCCCC'
    mol = chemi.smiles_to_oemol(smiles)
    f = fragmenter.fragment.WBOFragmenter(mol)
    f.fragment()

    qcschema_mol = f._to_qcschema_mol(f.fragments[(3, 5)])
    assert 'initial_molecule' in qcschema_mol
    assert 'geometry' in qcschema_mol['initial_molecule'][0]
    assert 'symbols' in qcschema_mol['initial_molecule'][0]
    assert 'connectivity' in qcschema_mol['initial_molecule'][0]
    assert 'identifiers' in qcschema_mol
    assert 'provenance' in qcschema_mol

def test_to_qcschema_mols():
    smiles = 'CCCCCC'
    mol = chemi.smiles_to_oemol(smiles)
    f = fragmenter.fragment.WBOFragmenter(mol)
    f.fragment()

    qcschema_mol = f.to_qcschema_mols()
    assert len(qcschema_mol) == 1

def test_td_inputs():
    smiles = 'CCCCCC'
    mol = chemi.smiles_to_oemol(smiles)
    f = fragmenter.fragment.WBOFragmenter(mol)
    f.fragment()

    td_inputs = f.to_torsiondrive_json()
    assert len(td_inputs) == 1
