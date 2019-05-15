"""
Unit and regression test for the fragmenter package.
"""

# Import package, test suite, and other packages as needed
import fragmenter
from fragmenter import chemi
from cmiles.utils import mol_to_smiles, restore_atom_map
import sys
import unittest
from fragmenter.tests.utils import get_fn, has_openeye, using_openeye


def test_fragmenter_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "fragmenter" in sys.modules


class TestFragment(unittest.TestCase):

    @unittest.skipUnless(has_openeye, 'Cannot test without OpenEye')
    def test_stereo_parent(self):
        """Test non isomeric and isomeric parent molecule SMILES"""
        smiles = 'NC(C)(F)C(=O)O'
        isomeric_smiles_r = 'N[C@](C)(F)C(=O)O'
        isomeric_smiles_s = 'N[C@@](C)(F)C(=O)O'
        mol_1 = chemi.smiles_to_oemol(smiles)
        mol_2 = chemi.smiles_to_oemol(isomeric_smiles_r)
        mol_3 = chemi.smiles_to_oemol(isomeric_smiles_s)


        self.assertFalse(fragmenter.fragment._generate_fragments(mol_1))
        fragmenter.fragment._generate_fragments(mol_1, strict_stereo=False)
        fragmenter.fragment._generate_fragments(mol_2)
        fragmenter.fragment._generate_fragments(mol_3)


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

    frags = fragmenter.fragment.Fragmenter(mapped_mol)
    frags.fragment_all_bonds_not_in_ring_systems()
    frags.combine_fragments(min_rotors=1, max_rotors=frags.n_rotors+1, restore_maps=True)

    keys = list(frags.fragment_combinations.keys())
    assert oechem.OEMolToSmiles(frags.fragment_combinations[keys[0]][0]) == '[H:45][c:8]1[cH:18][c:10]([c:19]([cH:17][c:7]1[H:44])[NH:35][H:67])[H:47]'
    assert oechem.OEMolToSmiles(frags.fragment_combinations[keys[1]][0]) == '[H:46][c:9]1[cH:20][n:32][c:21]([n:31][c:12]1[H:49])[NH:35][H:67]'


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
    frags = fragmenter.fragment.Fragmenter(mol)
    restore_atom_map(frags.molecule)
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
    f = fragmenter.fragment.Fragmenter(mol)
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
    f = fragmenter.fragment.Fragmenter(mol)
    bond = f.get_bond(bond_tuple=(3, 4))
    assert bond.IsRotor()

@using_openeye
def test_build_fragment():
    pass

def test_to_atom_bond_set():
    pass

def test_frag_to_mol():
    pass

def test_compare_wbo():
    pass

def test_find_ortho_substituent():
    pass
