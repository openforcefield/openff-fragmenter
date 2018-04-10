""" Tests generating files for qm torsion scan """

import unittest
from fragmenter.tests.utils import get_fn, has_openeye
import fragmenter.torsion_scan as qmscan
from fragmenter import utils
from openmoltools import openeye
import tempfile
import os
from fnmatch import fnmatch
import shutil


class TestQmscan(unittest.TestCase):

    @unittest.skipUnless(has_openeye, 'Cannot test without openeye')
    def test_generat_torsions(self):
        """ Tests finding torsion to drive """
        from openeye import oechem
        infile = get_fn('butane.pdb')
        ifs = oechem.oemolistream(infile)
        inp_mol = oechem.OEMol()
        oechem.OEReadMolecule(ifs, inp_mol)
        outfile_path = tempfile.mkdtemp()[1]
        qmscan.generate_torsions(inp_mol=inp_mol, output_path=outfile_path, interval=30, tar=False)
        input_files = []
        pattern = '*.pdb'
        for path, subdir, files in os.walk(outfile_path):
            for name in files:
                if fnmatch(name, pattern):
                    input_files.append(os.path.join(path, name))

        contents = open(input_files[0]).read()
        pdb = get_fn('butane_10_7_4_3_0.pdb')
        compare_contents = open(pdb).read()
        self.assertEqual(contents, compare_contents )

        shutil.rmtree(outfile_path)

    def test_generate_input(self):
        """Test generate psi4 input files"""
        root = get_fn('torsion_scan/10_7_4_3')
        qmscan.generate_scan_input(root, 'pdb', 'butane', ['MP2'], ['aug-cc-pvtz'], symmetry='C1')

        contents = open(get_fn('torsion_scan/10_7_4_3/0/butane_10_7_4_3_0.dat')).read()
        compare_content = open(get_fn('butane_10_7_4_3_0.dat')).read()
        self.assertEqual(contents, compare_content)

    @unittest.skipUnless(has_openeye, 'Cannot test without OpenEye')
    def test_tagged_smiles(self):
        """Test index-tagges smiles"""
        from openeye import oechem
        inf = get_fn('ethylmethylidyneamonium.mol2')
        ifs = oechem.oemolistream(inf)
        inp_mol = oechem.OEMol()
        oechem.OEReadMolecule(ifs, inp_mol)

        tagged_smiles = utils.create_mapped_smiles(inp_mol)

        # Tags should always be the same as mol2 molecule ordering
        self.assertEqual(tagged_smiles, '[H:5][C:1]#[N+:4][C:3]([H:9])([H:10])[C:2]([H:6])([H:7])[H:8]')

    @unittest.skipUnless(has_openeye, "Cannot test without OpneEye")
    def test_atom_map(self):
        """Test get atom map"""
        from openeye import oechem
        tagged_smiles = '[H:5][C:1]#[N+:4][C:3]([H:9])([H:10])[C:2]([H:6])([H:7])[H:8]'
        mol_1 = openeye.smiles_to_oemol('CC[N+]#C')
        inf = get_fn('ethylmethylidyneamonium.mol2')
        ifs = oechem.oemolistream(inf)
        mol_2 = oechem.OEMol()
        oechem.OEReadMolecule(ifs, mol_2)

        mol_1, atom_map = utils.get_atom_map(tagged_smiles, mol_1)

        for i, mapping in enumerate(atom_map):
            atom_1 = mol_1.GetAtom(oechem.OEHasAtomIdx(atom_map[mapping]))
            atom_1.SetAtomicNum(i+1)
            atom_2 = mol_2.GetAtom(oechem.OEHasAtomIdx(mapping-1))
            atom_2.SetAtomicNum(i+1)
            self.assertEqual(oechem.OECreateCanSmiString(mol_1), oechem.OECreateCanSmiString(mol_2))

        # Test aromatic molecule
        tagged_smiles = '[H:10][c:4]1[c:3]([c:2]([c:1]([c:6]([c:5]1[H:11])[H:12])[C:7]([H:13])([H:14])[H:15])[H:8])[H:9]'
        mol_1 = openeye.smiles_to_oemol('Cc1ccccc1')
        inf = get_fn('toluene.mol2')
        ifs = oechem.oemolistream(inf)
        mol_2 = oechem.OEMol()
        oechem.OEReadMolecule(ifs, mol_2)

        mol_1, atom_map = utils.get_atom_map(tagged_smiles, mol_1)
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
        mol_from_tagged_smiles = openeye.smiles_to_oemol(tagged_smiles)
        mol_1, atom_map = utils.get_atom_map(tagged_smiles, mol_from_tagged_smiles)

        # Compare atom map to tag
        for i in range(1, len(atom_map) +1):
            atom_1 = mol_from_tagged_smiles.GetAtom(oechem.OEHasAtomIdx(atom_map[i]))
            self.assertEqual(i, atom_1.GetMapIdx())


    @unittest.skipUnless(has_openeye, "Cannot test without OpneEye")
    def test_mapped_xyz(self):
        """Test writing out mapped xyz"""
        from openeye import oechem, oeomega
        tagged_smiles = '[H:10][c:4]1[c:3]([c:2]([c:1]([c:6]([c:5]1[H:11])[H:12])[C:7]([H:13])([H:14])[H:15])[H:8])[H:9]'
        mol_1 = openeye.smiles_to_oemol('Cc1ccccc1')
        inf = get_fn('toluene.mol2')
        ifs = oechem.oemolistream(inf)
        mol_2 = oechem.OEMol()
        oechem.OEReadMolecule(ifs, mol_2)

        mol_1, atom_map = utils.get_atom_map(tagged_smiles, mol_1)
        for i, mapping in enumerate(atom_map):
            atom_1 = mol_1.GetAtom(oechem.OEHasAtomIdx(atom_map[mapping]))
            atom_1.SetAtomicNum(i+1)
            atom_2 = mol_2.GetAtom(oechem.OEHasAtomIdx(mapping-1))
            atom_2.SetAtomicNum(i+1)

        xyz_1 = utils.to_mapped_xyz(mol_1, atom_map)
        # molecule generated from mol2 should be in the right order.
        atom_map_mol2 = {1:0, 2:1, 3:2, 4:3, 5:4, 6:5, 7:6, 8:7, 9:8, 10:9, 11:10, 12:11, 13:12, 14:13, 15:14}
        xyz_2 = utils.to_mapped_xyz(mol_2, atom_map_mol2)

        for ele1, ele2 in zip(xyz_1.split('\n')[:-1], xyz_2.split('\n')[:-1]):
            self.assertEqual(ele1.split(' ')[2], ele2.split(' ')[2])

