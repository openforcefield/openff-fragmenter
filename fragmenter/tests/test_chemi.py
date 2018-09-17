"""Test chemi module"""

from fragmenter import chemi
from fragmenter.tests.utils import has_openeye
import pytest
import numpy as np

@pytest.fixture
def mapped_molecule():
    return chemi.smiles_to_oemol('[H:5][C:1]([H:6])([H:7])[C:3]([H:11])([H:12])[C:4]([H:13])([H:14])[C:2]([H:8])([H:9])[H:10]')


def test_connectivity_table(mapped_molecule):
    """Test generate connectivity table"""
    expected_table = np.array([[0, 2, 1],
                        [1, 3, 1],
                        [2, 3, 1],
                        [0, 4, 1],
                        [0, 5, 1],
                        [0, 6, 1],
                        [1, 7, 1],
                        [1, 8, 1],
                        [1, 9, 1],
                        [2, 10, 1],
                        [2, 11, 1],
                        [3, 12, 1],
                        [3, 13, 1]])
    connectivity_table = chemi.get_mapped_connectivity_table(mapped_molecule)

    for bond in connectivity_table:
        xi = np.isin(expected_table, bond[:2])
        match = np.where(np.array([i[:2].sum() for i in xi]) == 2)[0]
        # assert that a match was found and only one was found
        assert len(match) == 1
        # assert that bond order is the same
        assert expected_table[match][0][-1] == bond[-1]




