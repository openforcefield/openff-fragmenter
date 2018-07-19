"""Test util functions"""

import unittest
import json
from fragmenter.tests.utils import get_fn, has_openeye
from fragmenter import utils
from openmoltools import openeye

class TesTorsions(unittest.TestCase):

    @unittest.skipUnless(has_openeye, 'Cannot test without openeye')
    def test_is_mapped(self):
        """Test checking if atom map exists"""
        smiles = 'CCCC'
        tagged_smiles = '[H:5][C:1]([H:6])([H:7])[C:2]([H:8])([H:9])[C:3]([H:10])([H:11])[C:4]([H:12])([H:13])[H:14]'
        molecule = openeye.smiles_to_oemol(smiles)
        tagged_molecule = openeye.smiles_to_oemol(tagged_smiles)

        self.assertTrue(utils.is_mapped(tagged_molecule))
        self.assertFalse(utils.is_mapped(molecule))

        # Add tags
        tagged_smiles = utils.create_mapped_smiles(molecule)
        self.assertTrue(utils.is_mapped(molecule))

    def test_check_molecule(self):
        """Test check moelcule"""
        pass

    def test_formal_charge(self):
        """Test formal charge"""

        mol_1 = utils.smiles_to_oemol('c1cc(c[nH+]c1)c2ccncn2')
        charge = utils.get_charge(mol_1)
        self.assertEqual(charge, 1)

        mol_2 = utils.smiles_to_oemol('C[NH+]1CC[NH+](CC1)Cc2ccccc2')
        charge = utils.get_charge(mol_2)
        self.assertEqual(charge, 2)

        mol_3 = utils.smiles_to_oemol('CCC(C)(C)C(=O)[O-]')
        charge = utils.get_charge(mol_3)
        self.assertEqual(charge, -1)



