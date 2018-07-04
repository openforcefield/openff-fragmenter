""" Test launch fragmenter """

import unittest
import fragmenter
from fragmenter import workflow_api
from fragmenter.tests.utils import get_fn


class TestWorkflow(unittest.TestCase):

    def test_get_provenance(self):
        """Test get provenance"""
        provenance = workflow_api.get_provenance(routine='enumerate_states')
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
        provenance = workflow_api.get_provenance(routine='enumerate_states', options=options)
        self.assertEqual(provenance['canonicalization_details'], canonicalization_details)
        self.assertTrue(provenance['routine']['enumerate_states']['keywords']['tautomers'])
        self.assertFalse(provenance['routine']['enumerate_states']['keywords']['stereoisomers'])

    def test_load_options(self):
        """Test load options"""
        options = workflow_api.load_options('enumerate_states')
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
        options = workflow_api.load_options(routine='enumerate_states', load_path=user_options)
        self.assertTrue(options['tautomers'])
        self.assertFalse(options['stereoisomers'])

    def test_remove_extraneous_options(self):
        """Test remove extraneous options"""

        default_options = workflow_api._default_options
        options = get_fn('options.yaml')
        user_options = workflow_api.load_options(routine='enumerate_states', load_path=options)
        needed_options = workflow_api.remove_extraneous_options(user_options, 'enumerate_states')

        with self.assertRaises(KeyError):
            self.assertTrue(needed_options['verbose'])
            self.assertTrue(default_options['verbose'])

        self.assertFalse(user_options['verbose'])

    def test_enumerate_states(self):
        """Test enumerate states"""

        states = workflow_api.enumerate_states('CCC(C)(C)C(=O)O')
        smiles = {'CCC(C)(C)C(=O)O': {'CCC(C)(C)C(=O)O', 'CCC(C)(C)C(=O)[O-]'}}
        self.assertEqual(len(states), 2)
        self.assertEqual(states['provenance']['routine']['enumerate_states']['keywords'], workflow_api._default_options['enumerate_states'])
        self.assertEqual(states['provenance']['routine']['enumerate_states']['version'], fragmenter.__version__)

        states.pop('provenance')
        self.assertEqual(states, smiles)

        pass


    # def test_launch_default(self):
    #     """Test default launch fragmenter"""
    #     fragments = launch.launch_fragmenter('CCC(C)(C)C(=O)O')
    #     self.assertEqual(len(fragments), 4)
    #     for frag in fragments:
    #         self.assertEqual(fragments[frag]['provenance']['routine_options'], launch._default_options)
    #         dimension = 0
    #         crank_dimension = 0
    #         torsion_type = ['mid', 'terminal']
    #         for torsion in torsion_type:
    #             dimension += len(fragments[frag]['needed_torsion_drives'][torsion]) - 1
    #             crank_dimension += len(fragments[frag]['crank_torsion_drives']['crank_job_0']['{}_torsions'.format(torsion)])
    #         self.assertEqual(dimension, crank_dimension)
    #
    # def test_omit_terminal_torsions(self):
    #     """Test skip terminal torions"""
    #     options = get_fn('options.yaml')
    #     fragments = launch.launch_fragmenter(molecule='CCC(C)(C)C(=O)O', options=options)
    #     self.assertEqual(len(fragments), 2)
    #     for frag in fragments:
    #         dimension = 0
    #         crank_dimension = 0
    #         torsion_type = ['mid', 'terminal']
    #         mid_d = len(fragments[frag]['needed_torsion_drives']['mid']) - 1
    #         term_d = len(fragments[frag]['needed_torsion_drives']['terminal']) - 1
    #         for torsion in torsion_type:
    #             dimension += len(fragments[frag]['needed_torsion_drives'][torsion]) - 1
    #             crank_dimension += len(fragments[frag]['crank_torsion_drives']['crank_job_0']['{}_torsions'.format(torsion)])
    #         self.assertEqual(dimension-term_d, crank_dimension)
    #         self.assertEqual(crank_dimension, mid_d)
    #         self.assertEqual(fragments[frag]['needed_torsion_drives']['mid']['grid_spacing'], 15)
    #         self.assertIsNone(fragments[frag]['needed_torsion_drives']['terminal']['grid_spacing'])
    #
    # def test_1D_scans(self):
    #     """Test 1D torsion scan grid"""
    #     options = get_fn('options2.yaml')
    #     fragments = launch.launch_fragmenter(molecule='CCC(C)(C)C(=O)O', options=options)
    #     self.assertEqual(len(fragments), 2)
    #     for frag in fragments:
    #         needed_t = fragments[frag]['needed_torsion_drives']
    #         jobs = len(needed_t['mid']['grid_spacing'])
    #         self.assertIsNone(needed_t['terminal']['grid_spacing'])
    #         crank_jobs = fragments[frag]['crank_torsion_drives']
    #         self.assertEqual(jobs, len(crank_jobs))
    #         for i, job in enumerate(crank_jobs):
    #             self.assertEqual(crank_jobs[job]['terminal_torsions'], {})
    #             print(job)
    #             print(crank_jobs[job]['mid_torsions']['torsion_{}'.format(str(i))])
    #             self.assertEqual(crank_jobs[job]['mid_torsions']['torsion_{}'.format(str(i))], 15)
    #
    # def test_crank_initial_state(self):
    #     """ Test generate crank initial state"""
    #     jsonfile = open(get_fn('butane_crankjob.json'), 'r')
    #     test_crank_job = json.load(jsonfile)
    #     jsonfile.close()
    #
    #     crank_initial_state = launch.get_initial_crank_state(test_crank_job['CCCC'])
    #     for dihedral in crank_initial_state['crank_job_0']['dihedrals']:
    #         self.assertTrue(dihedral in [[2, 1, 0, 4], [0, 1, 2, 3], [1, 2, 3, 11]])
    #     self.assertEqual(crank_initial_state['crank_job_0']['grid_spacing'], [30, 30, 30])
    #     self.assertFalse(crank_initial_state['crank_job_0']['grid_status'])





