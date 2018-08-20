""" Tests generating files for qm torsion scan """

import unittest
import json
from fragmenter.tests.utils import get_fn, has_openeye
import fragmenter.torsions as torsions
from fragmenter import utils, chemi
from cmiles import to_canonical_smiles_oe
import warnings


class TestTorsions(unittest.TestCase):

    @unittest.skipUnless(has_openeye, 'Cannot test without openeye')
    def test_generate_torsions(self):
        """ Tests finding torsion to drive """
        from openeye import oechem
        infile = get_fn('butane.pdb')
        ifs = oechem.oemolistream(infile)
        inp_mol = oechem.OEMol()
        oechem.OEReadMolecule(ifs, inp_mol)
        needed_torsion_scans = torsions.find_torsions(molecule=inp_mol)
        self.assertEqual(len(needed_torsion_scans['internal']), 1)
        self.assertEqual(len(needed_torsion_scans['terminal']), 2)
        self.assertEqual(needed_torsion_scans['internal']['torsion_0'], (14, 10, 7, 4))
        self.assertEqual(needed_torsion_scans['terminal']['torsion_0'], (10, 7, 4, 3))
        self.assertEqual(needed_torsion_scans['terminal']['torsion_1'], (7, 10, 14, 13))

    @unittest.skipUnless(has_openeye, 'Cannot test without OpenEye')
    def test_tagged_smiles(self):
        """Test index-tagges smiles"""
        from openeye import oechem
        inf = get_fn('ethylmethylidyneamonium.mol2')
        ifs = oechem.oemolistream(inf)
        inp_mol = oechem.OEMol()
        oechem.OEReadMolecule(ifs, inp_mol)

        tagged_smiles = to_canonical_smiles_oe(inp_mol, isomeric=True, mapped=True, explicit_hydrogen=True)

        # Tags should always be the same as mol2 molecule ordering
        self.assertEqual(tagged_smiles, '[H:5][C:1]#[N+:4][C:3]([H:9])([H:10])[C:2]([H:6])([H:7])[H:8]')

    @unittest.skipUnless(has_openeye, "Cannot test without OpneEye")
    def test_atom_map(self):
        """Test get atom map"""
        from openeye import oechem
        tagged_smiles = '[H:5][C:1]#[N+:4][C:3]([H:9])([H:10])[C:2]([H:6])([H:7])[H:8]'
        mol_1 = chemi.smiles_to_oemol('CC[N+]#C')
        inf = get_fn('ethylmethylidyneamonium.mol2')
        ifs = oechem.oemolistream(inf)
        mol_2 = oechem.OEMol()
        oechem.OEReadMolecule(ifs, mol_2)

        mol_1, atom_map = chemi.get_atom_map(tagged_smiles, mol_1)

        for i, mapping in enumerate(atom_map):
            atom_1 = mol_1.GetAtom(oechem.OEHasAtomIdx(atom_map[mapping]))
            atom_1.SetAtomicNum(i+1)
            atom_2 = mol_2.GetAtom(oechem.OEHasAtomIdx(mapping-1))
            atom_2.SetAtomicNum(i+1)
            self.assertEqual(oechem.OECreateCanSmiString(mol_1), oechem.OECreateCanSmiString(mol_2))

        # Test aromatic molecule
        tagged_smiles = '[H:10][c:4]1[c:3]([c:2]([c:1]([c:6]([c:5]1[H:11])[H:12])[C:7]([H:13])([H:14])[H:15])[H:8])[H:9]'
        mol_1 = chemi.smiles_to_oemol('Cc1ccccc1')
        inf = get_fn('toluene.mol2')
        ifs = oechem.oemolistream(inf)
        mol_2 = oechem.OEMol()
        oechem.OEReadMolecule(ifs, mol_2)

        mol_1, atom_map = chemi.get_atom_map(tagged_smiles, mol_1)
        for i, mapping in enumerate(atom_map):
            atom_1 = mol_1.GetAtom(oechem.OEHasAtomIdx(atom_map[mapping]))
            atom_1.SetAtomicNum(i+1)
            atom_2 = mol_2.GetAtom(oechem.OEHasAtomIdx(mapping-1))
            atom_2.SetAtomicNum(i+1)
            self.assertEqual(oechem.OECreateCanSmiString(mol_1), oechem.OECreateCanSmiString(mol_2))

    @unittest.skipUnless(has_openeye, "Cannot test without OpenEye")
    def test_atom_map_order(self):
        """Test atom map"""
        from openeye import oechem
        tagged_smiles = '[H:5][C:1]#[N+:4][C:3]([H:9])([H:10])[C:2]([H:6])([H:7])[H:8]'
        mol_from_tagged_smiles = chemi.smiles_to_oemol(tagged_smiles)
        mol_1, atom_map = chemi.get_atom_map(tagged_smiles, mol_from_tagged_smiles)

        # Compare atom map to tag
        for i in range(1, len(atom_map) +1):
            atom_1 = mol_from_tagged_smiles.GetAtom(oechem.OEHasAtomIdx(atom_map[i]))
            self.assertEqual(i, atom_1.GetMapIdx())

    @unittest.skipUnless(has_openeye, "Cannot test without OpneEye")
    def test_mapped_xyz(self):
        """Test writing out mapped xyz"""
        from openeye import oechem, oeomega
        tagged_smiles = '[H:10][c:4]1[c:3]([c:2]([c:1]([c:6]([c:5]1[H:11])[H:12])[C:7]([H:13])([H:14])[H:15])[H:8])[H:9]'
        mol_1 = chemi.smiles_to_oemol('Cc1ccccc1')
        inf = get_fn('toluene.mol2')
        ifs = oechem.oemolistream(inf)
        mol_2 = oechem.OEMol()
        oechem.OEReadMolecule(ifs, mol_2)

        mol_1, atom_map = chemi.get_atom_map(tagged_smiles, mol_1)
        for i, mapping in enumerate(atom_map):
            atom_1 = mol_1.GetAtom(oechem.OEHasAtomIdx(atom_map[mapping]))
            atom_1.SetAtomicNum(i+1)
            atom_2 = mol_2.GetAtom(oechem.OEHasAtomIdx(mapping-1))
            atom_2.SetAtomicNum(i+1)

        xyz_1 = chemi.to_mapped_xyz(mol_1, atom_map)
        # molecule generated from mol2 should be in the right order.
        atom_map_mol2 = {1:0, 2:1, 3:2, 4:3, 5:4, 6:5, 7:6, 8:7, 9:8, 10:9, 11:10, 12:11, 13:12, 14:13, 15:14}
        xyz_2 = chemi.to_mapped_xyz(mol_2, atom_map_mol2)

        for ele1, ele2 in zip(xyz_1.split('\n')[:-1], xyz_2.split('\n')[:-1]):
            self.assertEqual(ele1.split(' ')[2], ele2.split(' ')[2])

    @unittest.skipUnless(has_openeye, "Cannot test without OpenEye")
    def test_to_mapped_geometry(self):
        """Test mapped geometry"""
        from openeye import oechem

        infile = get_fn('butane.pdb')
        ifs = oechem.oemolistream(infile)
        molecule = oechem.OEMol()
        oechem.OEReadMolecule(ifs, molecule)
        tagged_smiles = to_canonical_smiles_oe(molecule, isomeric=True, explicit_hydrogen=True, mapped=True)
        molecule, atom_map = chemi.get_atom_map(tagged_smiles, molecule)
        mapped_geometry = chemi.to_mapped_QC_JSON_geometry(molecule, atom_map)

        f = open(infile)
        line = f.readline()
        symbols = []
        geometry = []
        while line.strip():
            if line.startswith('ATOM'):
                line = line.split()
                symbols.append(line[2][0])
                geometry.append(line[5:8])
            line = f.readline()
        f.close()
        geometry = sum(geometry, [])
        self.assertEqual(symbols, mapped_geometry['symbols'])
        for x, y in zip(geometry, mapped_geometry['geometry']):
            self.assertAlmostEqual(float(x), y, 3)

    @unittest.skipUnless(has_openeye, "Cannot test without OpenEye")
    def test_crank_specs(self):
        """Test crank job JSON"""

        test_crank = {'canonical_isomeric_SMILES': 'CCCC',
                      'needed_torsion_drives': {'internal': {
                                                        'torsion_0': (1, 2, 3, 4)},
                                                'terminal': {
                                                             'torsion_0': (3, 2, 1, 5),
                                                             'torsion_1': (2, 3, 4, 12)}},
                      'provenance': {'canonicalization': 'openeye v2017.Oct.1',
                                     'package': 'fragmenter',
                                     'parent_molecule': 'CCCC',
                                     'routine': 'fragmenter.fragment.generate_fragments',
                                     'routine_options': {'MAX_ROTORS': 2,
                                                         'combinatorial': True,
                                                         'generate_visualization': False,
                                                         'json_filename': None,
                                                         'remove_map': True,
                                                         'strict_stereo': True},
                                     'user': 'chayastern',
                                     'version': '0.0.0+29.g7d02c4d.dirty'},
                      'tagged_SMARTS': '[H:5][C:1]([H:6])([H:7])[C:2]([H:8])([H:9])[C:3]([H:10])([H:11])[C:4]([H:12])([H:13])[H:14]'}

        crank_job = torsions.define_crank_job(test_crank, terminal_torsion_resolution=30)
        self.assertEqual(crank_job['crank_torsion_drives']['crank_job_0']['internal_torsions']['torsion_0'], 30)
        self.assertEqual(crank_job['crank_torsion_drives']['crank_job_1']['terminal_torsions']['torsion_0'], 30)
        self.assertEqual(crank_job['crank_torsion_drives']['crank_job_1']['terminal_torsions']['torsion_1'], 30)
