"""Test util functions"""

import unittest
import json
from fragmenter.tests.utils import get_fn, has_openeye
from fragmenter import utils
from openmoltools import openeye
from openeye import oechem


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
        tagged_mol = utils.smiles_to_oemol(tagged_smiles)
        self.assertTrue(utils.is_mapped(tagged_mol))

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

    def test_create_mapped_smiles(self):
        """Test mapped SMILES"""
        smiles = 'CC(c1c(ccc(c1Cl)F)Cl)OC'
        mol = oechem.OEMol()
        oechem.OESmilesToMol(mol, smiles)

        canonical_smiles = utils.create_mapped_smiles(mol, tagged=False, explicit_hydrogen=False, isomeric=False)
        self.assertEqual(canonical_smiles, 'CC(c1c(ccc(c1Cl)F)Cl)OC')
        canonical_isomeric_smiles = utils.create_mapped_smiles(mol, tagged=False, explicit_hydrogen=False, isomeric=True)
        self.assertEqual(canonical_isomeric_smiles, 'C[C@@H](c1c(ccc(c1Cl)F)Cl)OC')
        canonical_explicit_h_smiles = utils.create_mapped_smiles(mol, tagged=False, explicit_hydrogen=True, isomeric=False)
        self.assertEqual(canonical_explicit_h_smiles, '[H]c1c(c(c(c(c1F)Cl)C([H])(C([H])([H])[H])OC([H])([H])[H])Cl)[H]')
        canonical_isomeric_explicit_h_smiles = utils.create_mapped_smiles(mol, tagged=False, explicit_hydrogen=True, isomeric=True)
        self.assertEqual(canonical_isomeric_explicit_h_smiles,
                         '[H]c1c(c(c(c(c1F)Cl)[C@]([H])(C([H])([H])[H])OC([H])([H])[H])Cl)[H]')

    def test_has_conformer(self):
        """Test has conformer"""
        infile = get_fn('butane.pdb')
        ifs = oechem.oemolistream(infile)
        molecule_with_conf = oechem.OEMol()
        oechem.OEReadMolecule(ifs, molecule_with_conf)

        self.assertTrue(utils.has_conformer(molecule_with_conf))

        molecule_without_conf = utils.smiles_to_oemol('CCCC')
        self.assertFalse(utils.has_conformer(molecule_without_conf))

    def test_2D_conformation(self):
        """Test checking for 2D conformation"""
        from fragmenter import fragment, chemi
        mol = utils.smiles_to_oemol('CCCC')
        states = fragment.expand_states(mol, return_molecules=True)
        for state in states:
            self.assertFalse(utils.has_conformer(state, check_two_dimension=True))

        conf = chemi.generate_conformers(mol, max_confs=1)
        self.assertTrue(utils.has_conformer(conf, check_two_dimension=True))





