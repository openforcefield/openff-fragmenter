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
        for flag in isomeric_smiles['Flags_set_to_True']:
            self.assertTrue(flag in ['ISOMERIC', 'Isotopes', 'AtomStereo', 'BondStereo', 'Canonical', 'AtomMaps', 'RGroups'])
        canonical_smiles = canon_detail['canonical_SMILES']
        for flag in canonical_smiles['Flags_set_to_True']:
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



