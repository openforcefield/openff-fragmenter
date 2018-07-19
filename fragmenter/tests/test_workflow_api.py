""" Test launch fragmenter """

import unittest
import fragmenter
from fragmenter import workflow_api
from fragmenter.tests.utils import get_fn
import json


class TestWorkflow(unittest.TestCase):

    def test_get_provenance(self):
        """Test get provenance"""
        provenance = workflow_api._get_provenance(routine='enumerate_states')
        canonicalization_details = {'canonical_isomeric_SMILES': {'Flags': ['ISOMERIC',
                                                                            'Isotopes',
                                                                            'AtomStereo',
                                                                            'BondStereo',
                                                                            'Canonical',
                                                                            'AtomMaps',
                                                                            'RGroups'],
                                    'oe_function': 'openeye.oechem.OEMolToSmiles(molecule)'},
                                    'package': 'openeye-v2017.Oct.1'}
        self.assertEqual(provenance['canonicalization_details'], canonicalization_details)
        default_kewyords = {'carbon_hybridization': True,
                            'level': 0,
                            'max_states': 200,
                            'protonation': True,
                            'reasonable': True,
                            'stereoisomers': True,
                            'suppress_hydrogen': True,
                            'tautomers': False}
        self.assertEqual(provenance['routine']['enumerate_states']['keywords'], default_kewyords)
        self.assertEqual(provenance['routine']['enumerate_states']['version'], fragmenter.__version__)

        options = get_fn('options.yaml')
        provenance = workflow_api._get_provenance(routine='enumerate_states', options=options)
        self.assertEqual(provenance['canonicalization_details'], canonicalization_details)
        self.assertTrue(provenance['routine']['enumerate_states']['keywords']['tautomers'])
        self.assertFalse(provenance['routine']['enumerate_states']['keywords']['stereoisomers'])

        provenance = workflow_api._get_provenance(routine='enumerate_fragments', options=options)
        self.assertEqual(provenance['canonicalization_details'], workflow_api._canonicalization_details)
        self.assertTrue(provenance['routine']['enumerate_fragments']['keywords']['generate_visualization'])
        self.assertFalse(provenance['routine']['enumerate_fragments']['keywords']['strict_stereo'])

    def test_load_options(self):
        """Test load options"""
        options = workflow_api._load_options('enumerate_states')
        default_options_wf = workflow_api._default_options['enumerate_states']
        default_options = {'carbon_hybridization': True,
                            'level': 0,
                            'max_states': 200,
                            'protonation': True,
                            'reasonable': True,
                            'stereoisomers': True,
                            'suppress_hydrogen': True,
                            'tautomers': False}
        self.assertEqual(options, default_options_wf)
        self.assertEqual(default_options_wf, default_options)
        self.assertEqual(options, default_options)

        user_options = get_fn('options.yaml')
        options = workflow_api._load_options(routine='enumerate_states', load_path=user_options)
        self.assertTrue(options['tautomers'])
        self.assertFalse(options['stereoisomers'])

        options = workflow_api._load_options('enumerate_fragments')
        default_options_wf = workflow_api._default_options['enumerate_fragments']
        default_options = {'strict_stereo': True,
                            'combinatorial': True,
                            'MAX_ROTORS': 2,
                            'remove_map': True}

        self.assertEqual(options, default_options_wf)
        self.assertEqual(default_options_wf, default_options)
        self.assertEqual(options, default_options)

        options = workflow_api._load_options(routine='enumerate_fragments', load_path=user_options)
        self.assertFalse(options['strict_stereo'])
        self.assertTrue(options['generate_visualization'])

        options = workflow_api._load_options('generate_crank_jobs')
        default_options_wf = workflow_api._default_options['generate_crank_jobs']
        default_options = {'max_conf': 1,
                            'terminal_torsion_resolution': 30,
                            'internal_torsion_resolution': 30,
                            'scan_internal_external_combination': 0,
                            'scan_dimension': 2,
                            'options':{
                            'qc_program': 'psi4',
                            'method': 'B3LYP',
                            'basis': 'aug-cc-pVDZ'}}

        self.assertEqual(options, default_options_wf)
        self.assertEqual(default_options_wf, default_options)
        self.assertEqual(options, default_options)

        options = workflow_api._load_options(routine='generate_crank_jobs', load_path=user_options)
        self.assertEqual(options['terminal_torsion_resolution'], 0)
        self.assertEqual(options['internal_torsion_resolution'], 15)

    def test_remove_extraneous_options(self):
        """Test remove extraneous options"""

        default_options = workflow_api._default_options
        options = get_fn('options.yaml')
        user_options = workflow_api._load_options(routine='enumerate_states', load_path=options)
        needed_options = workflow_api._remove_extraneous_options(user_options, 'enumerate_states')

        with self.assertRaises(KeyError):
            self.assertTrue(needed_options['verbose'])
            self.assertTrue(default_options['verbose'])

        self.assertFalse(user_options['verbose'])

        user_options = workflow_api._load_options(routine='enumerate_fragments', load_path=options)
        needed_options = workflow_api._remove_extraneous_options(user_options, 'enumerate_fragments')
        with self.assertRaises(KeyError):
            self.assertTrue(needed_options['generate_visualization'])

    def test_enumerate_states(self):
        """Test enumerate states"""

        states = workflow_api.enumerate_states('CCC(C)(C)C(=O)O')
        smiles = {'states': {'CCC(C)(C)C(=O)O', 'CCC(C)(C)C(=O)[O-]'}}
        self.assertEqual(len(states['states']), 2)
        self.assertEqual(states['provenance']['routine']['enumerate_states']['keywords'], workflow_api._default_options['enumerate_states'])
        self.assertEqual(states['provenance']['routine']['enumerate_states']['version'], fragmenter.__version__)

        states.pop('provenance')
        self.assertEqual(states, smiles)

    def test_enumerate_fragents(self):
        """Test enumerate fragments"""

        mol_smiles = 'CCCCC'
        fragments = workflow_api.enumerate_fragments(mol_smiles)
        self.assertEqual(len(fragments), 1)

        molecule = {'geometry': [0.31914281845092773,
                                  -1.093637466430664,
                                  -1.5644147396087646,
                                  0.09283685684204102,
                                  -0.7512494325637817,
                                  -0.10052239894866943,
                                  -0.09279406070709229,
                                  0.7513599395751953,
                                  0.1004934310913086,
                                  -0.3191012144088745,
                                  1.0937272310256958,
                                  1.564411997795105,
                                  0.4511583745479584,
                                  -2.172018527984619,
                                  -1.699379324913025,
                                  -0.5405434370040894,
                                  -0.7790045738220215,
                                  -2.1644840240478516,
                                  1.208341360092163,
                                  -0.5867938995361328,
                                  -1.9529553651809692,
                                  0.9486593008041382,
                                  -1.1027858257293701,
                                  0.48745453357696533,
                                  -0.7917245626449585,
                                  -1.287994146347046,
                                  0.26173627376556396,
                                  -0.9482856392860413,
                                  1.1030991077423096,
                                  -0.48761388659477234,
                                  0.7921697497367859,
                                  1.287682056427002,
                                  -0.26162096858024597,
                                  -0.4494825005531311,
                                  2.173631429672241,
                                  1.6855350732803345,
                                  0.5340865850448608,
                                  0.7838987112045288,
                                  2.176231622695923,
                                  -1.2162054777145386,
                                  0.5980985760688782,
                                  1.9490993022918701],
                                 'molecular_charge': 0,
                                 'molecular_multiplicity': 1,
                                 'symbols': ['C',
                                  'C',
                                  'C',
                                  'C',
                                  'H',
                                  'H',
                                  'H',
                                  'H',
                                  'H',
                                  'H',
                                  'H',
                                  'H',
                                  'H',
                                  'H']}
        self.assertEqual(fragments['CCCC']['molecule'],
                         molecule)

        mol_smiles_iso = 'N[C@H](C)CCF'
        frags_iso = fragmenter.workflow_api.enumerate_fragments(mol_smiles_iso)

        self.assertEqual(len(frags_iso.keys()), 2)
        iso_frag = frags_iso['CC[C@@H](C)N']
        self.assertEqual(iso_frag['SMILES']['canonical_SMILES'], 'CCC(C)N')
        self.assertEqual(iso_frag['SMILES']['canonical_isomeric_explicit_hydrogen_SMILES'],
                         '[H][C@@](C([H])([H])[H])(C([H])([H])C([H])([H])[H])N([H])[H]')
        self.assertEqual(iso_frag['SMILES']['canonical_explicit_hydrogen_SMILES'],
                         '[H]C([H])([H])C([H])([H])C([H])(C([H])([H])[H])N([H])[H]')

    def test_generate_crank_jobs(self):
        """Test generate crank jobs"""

        fragment = json.load(open(get_fn('CCCC.json'), 'r'))

        crank_jobs = workflow_api.generate_crank_jobs(fragment['CCCC'])

        self.assertEqual(len(crank_jobs), 2)
        self.assertEqual(crank_jobs['crank_job_0']['dihedrals'][0], [0, 1, 2, 3])
        self.assertEqual(crank_jobs['crank_job_0']['grid_spacing'], [30])
        self.assertEqual(crank_jobs['crank_job_1']['dihedrals'], [[2, 1, 0, 4], [1, 2, 3, 11]])
        self.assertEqual(crank_jobs['crank_job_1']['grid_spacing'], [30, 30])

        with self.assertRaises(Warning):
            workflow_api.generate_crank_jobs(fragment['CCCC'], options=get_fn('options.yaml'))

    def test_workflow(self):
        """Test workflow"""

        smiles_list = ['CCCC', 'CCCCCC']
        crank_jobs = workflow_api.workflow(smiles_list, write_json_crank_job=False)

        self.assertEqual(len(crank_jobs.keys()), 1)
        self.assertEqual(list(crank_jobs.keys())[0], 'CCCC')
        self.assertEqual(len(crank_jobs['CCCC'].keys()), 2)
        self.assertEqual(len(crank_jobs['CCCC']['crank_job_0']['dihedrals']), 1)
        self.assertEqual(len(crank_jobs['CCCC']['crank_job_1']['dihedrals']), 2)
        self.assertEqual(crank_jobs['CCCC']['crank_job_0']['provenance'], crank_jobs['CCCC']['crank_job_1']['provenance'])
        self.assertEqual(len(crank_jobs['CCCC']['crank_job_0']['provenance']['SMILES']), 5)
        self.assertEqual(crank_jobs['CCCC']['crank_job_0']['provenance']['SMILES']['canonical_isomeric_explicit_hydrogen_SMILES'],
                         crank_jobs['CCCC']['crank_job_0']['provenance']['SMILES']['canonical_explicit_hydrogen_SMILES'])







