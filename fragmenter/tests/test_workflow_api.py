""" Test launch fragmenter """

import pytest
import qcfractal.interface as portal
import fragmenter
from fragmenter import workflow_api, chemi
from fragmenter.tests.utils import get_fn, has_crank, has_openeye
import json
import copy

# Import QCFractal and bounce if not available
qcfractal = pytest.importorskip("qcfractal")
from qcfractal import testing
from qcfractal.testing import fractal_compute_server

@testing.using_rdkit
@testing.using_geometric
@testing.using_torsiondrive
def test_workflow(fractal_compute_server):
    """Fragmenter regression test"""

    pytest.importorskip("torsiondrive")
    pytest.importorskip("geometric")
    pytest.importorskip("rdkit")

    client = portal.FractalClient(fractal_compute_server)
    smiles = 'CCCCC'
    workflow_id = 'example'
    workflow_json = get_fn('workflows.json')

    workflow = workflow_api.WorkFlow(client=client, workflow_json=workflow_json, workflow_id=workflow_id)
    workflow.workflow(molecules_smiles=smiles, molecule_titles=['pentane'])

    assert len(workflow.qcfractal_jobs) == 1
    key = list(workflow.qcfractal_jobs.keys())[0]
    assert len(workflow.qcfractal_jobs[key]) == 3
    assert len(workflow.qcfractal_jobs[key]['optimization_input']) == 0
    assert len(workflow.qcfractal_jobs[key]['torsiondrive_input']) == 3
    assert len(workflow.qcfractal_jobs[key]['torsiondrive_input']['(0, 2, 3, 1)']['initial_molecule']) == 8
# class TestWorkflow(unittest.TestCase):
#
#     def test_get_provenance(self):
#         """Test get provenance"""
#         provenance = workflow_api._get_provenance(workflow_id='workflow_1', routine='enumerate_states')
#
#         self.assertIn('workflow_id', provenance)
#         self.assertIn('enumerate_states', provenance['routine'])
#
#         provenance = workflow_api._get_provenance(workflow_id='workflow_1', routine='enumerate_fragments')
#         self.assertIn('workflow_id', provenance)
#         self.assertIn('enumerate_fragments', provenance['routine'])
#
#
#         # default_kewyords = {'carbon_hybridization': True,
#         #                     'level': 0,
#         #                     'max_states': 200,
#         #                     'protonation': True,
#         #                     'reasonable': True,
#         #                     'stereoisomers': True,
#         #                     'suppress_hydrogen': True,
#         #                     'tautomers': False}
#         # self.assertEqual(provenance['routine']['enumerate_states']['keywords'], default_kewyords)
#         # self.assertEqual(provenance['routine']['enumerate_states']['version'], fragmenter.__version__)
#         #
#         # #options = get_fn('options.yaml')
#         # provenance = workflow_api._get_provenance(workflow_id='workflow_1', routine='enumerate_states')
#         # self.assertTrue(provenance['routine']['enumerate_states']['keywords']['tautomers'])
#         # self.assertFalse(provenance['routine']['enumerate_states']['keywords']['stereoisomers'])
#         #
#         # provenance = workflow_api._get_provenance(workflow_id='workflow_1', routine='enumerate_fragments')
#         # self.assertTrue(provenance['routine']['enumerate_fragments']['keywords']['generate_visualization'])
#         # self.assertFalse(provenance['routine']['enumerate_fragments']['keywords']['strict_stereo'])
#
#     def test_load_options(self):
#         """Test load options"""
#         options_1 = workflow_api._get_options(workflow_id='workflow_1', all=True)
#
#         fn = get_fn('workflows.json')
#         f = open(fn)
#         options_2 = json.load(f)['workflow_1']['fragmenter']
#
#         self.assertEqual(options_1, options_2)
#         # default_options_wf = workflow_api._default_options['enumerate_states']
#         # default_options = {'carbon_hybridization': True,
#         #                     'level': 0,
#         #                     'max_states': 200,
#         #                     'protonation': True,
#         #                     'reasonable': True,
#         #                     'stereoisomers': True,
#         #                     'suppress_hydrogen': True,
#         #                     'tautomers': False}
#         # self.assertEqual(options, default_options_wf)
#         # self.assertEqual(default_options_wf, default_options)
#         # self.assertEqual(options, default_options)
#         #
#         # user_options = get_fn('options.yaml')
#         # options = workflow_api._load_options(routine='enumerate_states', load_path=user_options)
#         # self.assertTrue(options['tautomers'])
#         # self.assertFalse(options['stereoisomers'])
#         #
#         # options = workflow_api._load_options('enumerate_fragments')
#         # default_options_wf = workflow_api._default_options['enumerate_fragments']
#         # default_options = {'strict_stereo': True,
#         #                     'combinatorial': True,
#         #                     'MAX_ROTORS': 2,
#         #                     'remove_map': True}
#         #
#         # self.assertEqual(options, default_options_wf)
#         # self.assertEqual(default_options_wf, default_options)
#         # self.assertEqual(options, default_options)
#         #
#         # options = workflow_api._load_options(routine='enumerate_fragments', load_path=user_options)
#         # self.assertFalse(options['strict_stereo'])
#         # self.assertTrue(options['generate_visualization'])
#         #
#         # options = workflow_api._load_options('generate_crank_jobs')
#         # default_options_wf = workflow_api._default_options['generate_crank_jobs']
#         # default_options = {'max_conf': 1,
#         #                     'terminal_torsion_resolution': 30,
#         #                     'internal_torsion_resolution': 30,
#         #                     'scan_internal_external_combination': 0,
#         #                     'scan_dimension': 2,
#         #                     'options':{
#         #                     'qc_program': 'psi4',
#         #                     'method': 'B3LYP',
#         #                     'basis': 'aug-cc-pVDZ'}}
#         #
#         # self.assertEqual(options, default_options_wf)
#         # self.assertEqual(default_options_wf, default_options)
#         # self.assertEqual(options, default_options)
#         #
#         # options = workflow_api._load_options(routine='generate_crank_jobs', load_path=user_options)
#         # self.assertEqual(options['terminal_torsion_resolution'], 0)
#         # self.assertEqual(options['internal_torsion_resolution'], 15)
#
#     # def test_remove_extraneous_options(self):
#     #     """Test remove extraneous options"""
#     #
#     #     default_options = workflow_api._default_options
#     #     options = get_fn('options.yaml')
#     #     user_options = workflow_api._load_options(routine='enumerate_states', load_path=options)
#     #     needed_options = workflow_api._remove_extraneous_options(user_options, 'enumerate_states')
#     #
#     #     with self.assertRaises(KeyError):
#     #         self.assertTrue(needed_options['verbose'])
#     #         self.assertTrue(default_options['verbose'])
#     #
#     #     self.assertFalse(user_options['verbose'])
#     #
#     #     user_options = workflow_api._load_options(routine='enumerate_fragments', load_path=options)
#     #     needed_options = workflow_api._remove_extraneous_options(user_options, 'enumerate_fragments')
#     #     with self.assertRaises(KeyError):
#     #         self.assertTrue(needed_options['generate_visualization'])
#
#     def test_enumerate_states(self):
#         """Test enumerate states"""
#
#         states = workflow_api.enumerate_states(molecule='CCC(C)(C)C(=O)O', workflow_id='workflow_1')
#         smiles = {'states': {'CCC(C)(C)C(=O)O', 'CCC(C)(C)C(=O)[O-]'}}
#         self.assertEqual(len(states['states']), 2)
#         self.assertEqual(states['provenance']['routine']['enumerate_states']['parent_molecule'], 'CCC(C)(C)C(=O)O')
#         self.assertEqual(states['provenance']['routine']['enumerate_states']['version'], fragmenter.__version__)
#
#         states.pop('provenance')
#         self.assertEqual(states, smiles)
#
#     def test_enumerate_fragments(self):
#         """Test enumerate fragments"""
#
#         mol_smiles = 'CCCCC'
#         fragments = workflow_api.enumerate_fragments(molecule=mol_smiles, workflow_id='workflow_1')
#         self.assertEqual(len(fragments), 1)
#
#         molecule = {'geometry': [0.31914281845092773,
#                                   -1.093637466430664,
#                                   -1.5644147396087646,
#                                   0.09283685684204102,
#                                   -0.7512494325637817,
#                                   -0.10052239894866943,
#                                   -0.09279406070709229,
#                                   0.7513599395751953,
#                                   0.1004934310913086,
#                                   -0.3191012144088745,
#                                   1.0937272310256958,
#                                   1.564411997795105,
#                                   0.4511583745479584,
#                                   -2.172018527984619,
#                                   -1.699379324913025,
#                                   -0.5405434370040894,
#                                   -0.7790045738220215,
#                                   -2.1644840240478516,
#                                   1.208341360092163,
#                                   -0.5867938995361328,
#                                   -1.9529553651809692,
#                                   0.9486593008041382,
#                                   -1.1027858257293701,
#                                   0.48745453357696533,
#                                   -0.7917245626449585,
#                                   -1.287994146347046,
#                                   0.26173627376556396,
#                                   -0.9482856392860413,
#                                   1.1030991077423096,
#                                   -0.48761388659477234,
#                                   0.7921697497367859,
#                                   1.287682056427002,
#                                   -0.26162096858024597,
#                                   -0.4494825005531311,
#                                   2.173631429672241,
#                                   1.6855350732803345,
#                                   0.5340865850448608,
#                                   0.7838987112045288,
#                                   2.176231622695923,
#                                   -1.2162054777145386,
#                                   0.5980985760688782,
#                                   1.9490993022918701],
#                                  'molecular_charge': 0,
#                                  'molecular_multiplicity': 1,
#                                  'symbols': ['C',
#                                   'C',
#                                   'C',
#                                   'C',
#                                   'H',
#                                   'H',
#                                   'H',
#                                   'H',
#                                   'H',
#                                   'H',
#                                   'H',
#                                   'H',
#                                   'H',
#                                   'H']}
#         #self.assertEqual(fragments['CCCC']['molecule']['geometry'],
#         #                 molecule['geometry'])
#
#         mol_smiles_iso = 'N[C@H](C)CCF'
#         frags_iso = fragmenter.workflow_api.enumerate_fragments(molecule=mol_smiles_iso, workflow_id='workflow_1')
#
#         self.assertEqual(len(frags_iso.keys()), 2)
#         iso_frag = frags_iso['CC[C@@H](C)N']
#         self.assertEqual(iso_frag['identifiers']['canonical_smiles'], 'CCC(C)N')
#         self.assertEqual(iso_frag['identifiers']['canonical_isomeric_explicit_hydrogen_smiles'],
#                          '[H][C@@](C([H])([H])[H])(C([H])([H])C([H])([H])[H])N([H])[H]')
#         self.assertEqual(iso_frag['identifiers']['canonical_explicit_hydrogen_smiles'],
#                          '[H]C([H])([H])C([H])([H])C([H])(C([H])([H])[H])N([H])[H]')
#
#     def test_generate_crank_jobs(self):
#         """Test generate crank jobs"""
#
#         fragment = json.load(open(get_fn('CCCC.json'), 'r'))
#
#         crank_jobs = workflow_api.generate_torsiondrive_input(fragment['CCCC'], workflow_id='workflow_1')
#
#         key = '[H:5][C:1]([H:6])([H:7])[C:3]([H:11])([H:12])[C:4]([H:13])([H:14])[C:2]([H:8])([H:9])[H:10]'
#         expected_dih = [(0, 2, 3, 1), (3, 2, 0, 4), (2, 3, 1, 7)]
#         all_dih = crank_jobs[key]['torsiondrive_input']['job_0']['dihedrals'] + crank_jobs[key]['torsiondrive_input']['job_1']['dihedrals']
#         for dih in all_dih:
#             self.assertIn(dih, expected_dih)
#
#         self.assertEqual(len(crank_jobs), 1)
#         self.assertEqual(len(crank_jobs[key]['torsiondrive_input']), 2)
#         self.assertEqual(crank_jobs[key]['torsiondrive_input']['job_0']['grid_spacing'][0], 30)
#         self.assertEqual(crank_jobs[key]['torsiondrive_input']['job_1']['grid_spacing'][0], 30)
#
#     def test_workflow(self):
#         """Test workflow"""
#
#         smiles_list = ['CCCC', 'CCCCCC']
#         crank_jobs = workflow_api.workflow(smiles_list, write_json_intermediate=False, workflow_id='workflow_1')
#
#         key = '[H:5][C:1]([H:6])([H:7])[C:3]([H:11])([H:12])[C:4]([H:13])([H:14])[C:2]([H:8])([H:9])[H:10]'
#
#         self.assertEqual(len(crank_jobs.keys()), 1)
#         self.assertEqual(list(crank_jobs.keys())[0], key)
#         self.assertEqual(len(crank_jobs[key].keys()), 2)
#
#         #self.assertEqual(crank_jobs[key]['torsiondrive_input']['job_0']['provenance'], crank_jobs[key]['torsiondrive_input']['job_1']['provenance'])
#         self.assertEqual(len(crank_jobs[key]['torsiondrive_input']['job_0']['initial_molecule']['identifiers']), 7)
#         self.assertEqual(crank_jobs[key]['torsiondrive_input']['job_0']['initial_molecule']['identifiers']['canonical_isomeric_explicit_hydrogen_smiles'],
#                          crank_jobs[key]['torsiondrive_input']['job_0']['initial_molecule']['identifiers']['canonical_explicit_hydrogen_smiles'])
#
#     # @unittest.skipUnless(has_crank, 'Cannot test without crank')
#     # def test_crank(self):
#     #     """Test fragmenter interfacing with crank"""
#     #
#     #     from crank import crankAPI
#     #     crank_jobs = workflow_api.workflow(['CCCC'], write_json_crank_job=False)
#     #
#     #     for mol in crank_jobs:
#     #         for job in crank_jobs[mol]:
#     #             state = copy.deepcopy(crank_jobs[mol][job])
#     #             next_job = crankAPI.next_jobs_from_state(state)
#     #             crankAPI.update_state(state, next_job)
#     #
#     #             self.assertEqual(crank_jobs[mol][job]['dihedrals'], state['dihedrals'])
#     #             self.assertEqual(state['grid_status'], next_job)
#
#     @unittest.skipUnless(has_openeye, 'Cannot test without OpenEye')
#     def test_dihedral_numbering(self):
#         """Test dihedral indices correspond to attached atoms"""
#
#         from openeye import oechem
#         crank_jobs = workflow_api.workflow(['CCCC'], workflow_id='workflow_1')
#         key = '[H:5][C:1]([H:6])([H:7])[C:3]([H:11])([H:12])[C:4]([H:13])([H:14])[C:2]([H:8])([H:9])[H:10]'
#
#         mol_with_map = chemi.smiles_to_oemol(key)
#
#         for job in crank_jobs[key]['torsiondrive_input']:
#             dihedrals = crank_jobs[key]['torsiondrive_input'][job]['dihedrals']
#             for dihedral in dihedrals:
#                 prev_atom = mol_with_map.GetAtom(oechem.OEHasMapIdx(dihedral[0]+1))
#                 for d in dihedral[1:]:
#                     atom = mol_with_map.GetAtom(oechem.OEHasMapIdx(d+1))
#                     bond = mol_with_map.GetBond(prev_atom, atom)
#                     self.assertIsNotNone(bond)
#                     prev_atom = atom
#








