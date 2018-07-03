""" Test launch fragmenter """

import unittest
from fragmenter import launch, utils
from fragmenter.tests.utils import get_fn


class TestLaunch(unittest.TestCase):

    def test_launch_default(self):
        """Test default launch fragmenter"""
        fragments = launch.launch_fragmenter('CCC(C)(C)C(=O)O')
        self.assertEqual(len(fragments), 4)
        for frag in fragments:
            self.assertEqual(fragments[frag]['provenance']['routine_options'], launch._default_options)
            dimension = 0
            crank_dimension = 0
            torsion_type = ['mid', 'terminal']
            for torsion in torsion_type:
                dimension += len(fragments[frag]['needed_torsion_drives'][torsion]) - 1
                crank_dimension += len(fragments[frag]['crank_torsion_drives']['crank_job_0']['{}_torsions'.format(torsion)])
            self.assertEqual(dimension, crank_dimension)

    def test_omit_terminal_torsions(self):
        """Test skip terminal torions"""
        options = get_fn('options.yaml')
        fragments = launch.launch_fragmenter(molecule='CCC(C)(C)C(=O)O', options=options)
        self.assertEqual(len(fragments), 2)
        for frag in fragments:
            dimension = 0
            crank_dimension = 0
            torsion_type = ['mid', 'terminal']
            mid_d = len(fragments[frag]['needed_torsion_drives']['mid']) - 1
            term_d = len(fragments[frag]['needed_torsion_drives']['terminal']) - 1
            for torsion in torsion_type:
                dimension += len(fragments[frag]['needed_torsion_drives'][torsion]) - 1
                crank_dimension += len(fragments[frag]['crank_torsion_drives']['crank_job_0']['{}_torsions'.format(torsion)])
            self.assertEqual(dimension-term_d, crank_dimension)
            self.assertEqual(crank_dimension, mid_d)
            self.assertEqual(fragments[frag]['needed_torsion_drives']['mid']['grid_spacing'], 15)
            self.assertIsNone(fragments[frag]['needed_torsion_drives']['terminal']['grid_spacing'])

    def test_1D_scans(self):
        """Test 1D torsion scan grid"""
        options = get_fn('options2.yaml')
        fragments = launch.launch_fragmenter(molecule='CCC(C)(C)C(=O)O', options=options)
        self.assertEqual(len(fragments), 2)
        for frag in fragments:
            needed_t = fragments[frag]['needed_torsion_drives']
            jobs = len(needed_t['mid']['grid_spacing'])
            self.assertIsNone(needed_t['terminal']['grid_spacing'])
            crank_jobs = fragments[frag]['crank_torsion_drives']
            self.assertEqual(jobs, len(crank_jobs))
            for i, job in enumerate(crank_jobs):
                self.assertEqual(crank_jobs[job]['terminal_torsions'], {})
                print(job)
                print(crank_jobs[job]['mid_torsions']['torsion_{}'.format(str(i))])
                self.assertEqual(crank_jobs[job]['mid_torsions']['torsion_{}'.format(str(i))], 15)

    def test_initial_crank_job(self):
        """Test initial crank job"""






