"""
Unit and regression test for the fragmenter package.
"""

# Import package, test suite, and other packages as needed
import fragmenter
import pytest
import sys
import unittest
from fragmenter.tests.utils import get_fn, has_openeye


def test_fragmenter_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "fragmenter" in sys.modules


class TestFragment(unittest.TestCase):

    @unittest.skipUnless(has_openeye, "Cannot test without openeye")
    def test_generate_fragments(self):
        """ generate fragments regression test"""
        input_smi = get_fn('butane.smi')
        fragments = fragmenter.generate_fragments(inputf=input_smi)

        provenance = fragments['provenance']
        canon_detail = provenance['canonicalization_details']
        self.assertTrue(canon_detail['AtomMaps'])
        self.assertTrue(canon_detail['AtomStereo'])
        self.assertTrue(canon_detail['BondStereo'])
        self.assertTrue(canon_detail['Canonical'])
        self.assertTrue(canon_detail['Isotopes'])
        self.assertTrue(canon_detail['RGroups'])
        self.assertTrue(canon_detail['ISOMERIC'])
        self.assertFalse(canon_detail['DEFAULT'])

        with self.assertRaises(Warning):
            fragmenter.generate_fragments(input_smi, OESMILESFlag='DEFAULT')

        fragments = fragmenter.generate_fragments(input_smi, strict_stereo=False, OESMILESFlag='DEFAULT')
        provenance = fragments['provenance']
        canon_detail = provenance['canonicalization_details']
        self.assertTrue(canon_detail['AtomMaps'])
        self.assertFalse(canon_detail['AtomStereo'])
        self.assertFalse(canon_detail['BondStereo'])
        self.assertTrue(canon_detail['Canonical'])
        self.assertFalse(canon_detail['Isotopes'])
        self.assertTrue(canon_detail['RGroups'])
        self.assertFalse(canon_detail['ISOMERIC'])
        self.assertTrue(canon_detail['DEFAULT'])
