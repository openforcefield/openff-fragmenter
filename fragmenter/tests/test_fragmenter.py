"""
Unit and regression test for the fragmenter package.
"""

# Import package, test suite, and other packages as needed
import fragmenter
import sys
import unittest
from openmoltools import openeye
from fragmenter.tests.utils import get_fn, has_openeye


def test_fragmenter_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "fragmenter" in sys.modules


class TestFragment(unittest.TestCase):

    def test_stereo_parent(self):
        """Test non isomeric and isomeric parent molecule SMILES"""
        smiles = 'NC(C)(F)C(=O)O'
        isomeric_smiles_r = 'N[C@](C)(F)C(=O)O'
        isomeric_smiles_s = 'N[C@@](C)(F)C(=O)O'
        mol_1 = openeye.smiles_to_oemol(smiles)
        mol_2 = openeye.smiles_to_oemol(isomeric_smiles_r)
        mol_3 = openeye.smiles_to_oemol(isomeric_smiles_s)


        self.assertFalse(fragmenter.fragment._generate_fragments(mol_1))
        fragmenter.fragment._generate_fragments(mol_1, strict_stereo=False)
        fragmenter.fragment._generate_fragments(mol_2)
        fragmenter.fragment._generate_fragments(mol_3)

    def test_expand_protonation_states(self):
        """Test expand protonation states"""
        smiles = 'C5=C(C1=CN=CC=C1)N=C(NC2=C(C=CC(=C2)NC(C3=CC=C(C=C3)CN4CCN(CC4)C)=O)C)N=C5'
        molecule = openeye.smiles_to_oemol(smiles)
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
            protonation_2.add(fragmenter.utils.create_mapped_smiles(mol, tagged=False, explicit_hydrogen=False))

        intersection = protonation_1.intersection(protonation_2)
        self.assertEqual(len(intersection), len(protonation_1))
        self.assertEqual(len(intersection), len(protonation_2))

    def test_expand_tautomers(self):
        """Test expand tautomer"""
        smiles_1 ='c1ccc2c(c1)C=CCC2=O'
        smiles_2 = 'c1ccc2c(c1)cccc2O'
        molecule_1 = fragmenter.utils.smiles_to_oemol(smiles_1)
        molecule_2 = fragmenter.utils.smiles_to_oemol(smiles_2)
        tautomers_1 = fragmenter.fragment.expand_states(molecule_1, protonation=False, tautomers=True, stereoisomers=False)
        tautomers_2 = fragmenter.fragment.expand_states(molecule_2, protonation=False, tautomers=True, stereoisomers=False)

        self.assertEqual(tautomers_1, tautomers_2)

    def test_expand_enantiomers(self):
        smiles = 'CN(C)C/C=C/C(=O)NC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC(=C(C=C3)F)Cl)O[C@H]4CCOC4'
        molecule = openeye.smiles_to_oemol(smiles)
        stereoisomers = fragmenter.fragment._expand_states(molecule, enumerate='stereoisomers')

        stereoisomers_1 = {'CN(C)C/C=C/C(=O)Nc1cc2c(cc1O[C@@H]3CCOC3)ncnc2Nc4ccc(c(c4)Cl)F',
                        'CN(C)C/C=C/C(=O)Nc1cc2c(cc1O[C@H]3CCOC3)ncnc2Nc4ccc(c(c4)Cl)F',
                        'CN(C)C/C=C\\C(=O)Nc1cc2c(cc1O[C@@H]3CCOC3)ncnc2Nc4ccc(c(c4)Cl)F',
                        'CN(C)C/C=C\\C(=O)Nc1cc2c(cc1O[C@H]3CCOC3)ncnc2Nc4ccc(c(c4)Cl)F'}

        stereoisomers_2 = set()
        for mol in stereoisomers:
            stereoisomers_2.add(fragmenter.utils.create_mapped_smiles(mol, tagged=False, explicit_hydrogen=False))
        intersection = stereoisomers_1.intersection(stereoisomers_2)
        self.assertEqual(len(intersection), len(stereoisomers_1))
        self.assertEqual(len(intersection), len(stereoisomers_2))
        self.assertEqual(len(stereoisomers_1), len(stereoisomers_2))



