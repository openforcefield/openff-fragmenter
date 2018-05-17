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

    @unittest.skipUnless(has_openeye, "Cannot test without openeye")
    def test_canoicalization_details(self):
        """ test canonicalization detail"""
        input_smi = get_fn('butane.smi')
        fragments = fragmenter.generate_fragments(inputf=input_smi)

        provenance = fragments['provenance']
        canon_detail = provenance['canonicalization_details']
        isomeric_smiles = canon_detail['canonical_isomeric_SMILES']
        for flag in isomeric_smiles['Flags']:
            self.assertTrue(flag in ['ISOMERIC', 'Isotopes', 'AtomStereo', 'BondStereo', 'Canonical', 'AtomMaps', 'RGroups'])
        canonical_smiles = canon_detail['canonical_SMILES']
        for flag in canonical_smiles['Flags']:
            print(flag)
            self.assertTrue(flag in ['DEFAULT', 'AtomMaps', 'Canonical', 'RGroups'])

    def test_stereo_parent(self):
        """Test non isomeric and isomeric parent molecule SMILES"""
        smiles = 'NC(C)(F)C(=O)O'
        isomeric_smiles_r = 'N[C@](C)(F)C(=O)O'
        isomeric_smiles_s = 'N[C@@](C)(F)C(=O)O'
        mol_1 = openeye.smiles_to_oemol(smiles)
        mol_2 = openeye.smiles_to_oemol(isomeric_smiles_r)
        mol_3 = openeye.smiles_to_oemol(isomeric_smiles_s)

        with self.assertRaises(RuntimeError):
            fragmenter.fragment._generate_fragments(mol_1)
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
        smiles ='C1=C(SC(=N1)NC2=CC(=NC(=N2)C)N3CCN(CC3)CCO)C(=O)NC4=C(C=CC=C4Cl)C'
        molecule = openeye.smiles_to_oemol(smiles)
        tautomers = fragmenter.fragment._expand_states(molecule, enumerate='tautomers')
        tautomers_1 = {'Cc1cccc(c1NC(=C2C=NC(=[NH+]c3cc(nc(n3)C)N4CCN(CC4)CCO)S2)[O-])Cl',
                        'Cc1cccc(c1NC(=O)c2cnc(s2)Nc3cc(nc(n3)C)N4CCN(CC4)CCO)Cl',
                        'Cc1cccc(c1[NH+]=C(c2cnc(s2)Nc3cc(nc(n3)C)N4CCN(CC4)CCO)[O-])Cl'}
        tautomers_2 = set()
        for mol in tautomers:
            tautomers_2.add(fragmenter.utils.create_mapped_smiles(mol, tagged=False, explicit_hydrogen=False))

        intersection = tautomers_1.intersection(tautomers_2)
        self.assertEqual(len(tautomers_1), len(intersection))
        self.assertEqual(len(tautomers_2), len(intersection))
        self.assertEqual(len(tautomers_1), len(tautomers_2))

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



