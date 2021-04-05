from contextlib import nullcontext

import pytest
from openff.toolkit.topology import Molecule

from fragmenter.utils import get_fgroup_smarts, get_fgroup_smarts_comb, get_map_index


def test_get_fgroup_smarts():
    """Tests that the `get_fgroup_smarts` utility returns correctly."""

    smarts = get_fgroup_smarts()

    assert "hydrazine" in smarts
    assert smarts["hydrazine"] == "[NX3][NX3]"

    assert "phosphon" not in smarts

    assert len(smarts) == 21


def test_get_fgroup_smarts_comb():
    """Tests that the `get_fgroup_smarts_comb` utility returns correctly."""

    smarts = get_fgroup_smarts_comb()

    assert "phosphon" in smarts
    assert smarts["phosphon"] == "[PX4](=[OX1])([OX2])([OX2])"

    assert "amide_2" not in smarts

    assert len(smarts) == 12


def test_get_map_index():

    molecule = Molecule.from_smiles("[C:5]([H:1])([H:2])([H:3])([H:4])")
    assert get_map_index(molecule, 0) == 5


@pytest.mark.parametrize(
    "raise_error, expected_raises",
    [
        (False, nullcontext()),
        (True, pytest.raises(KeyError, match="is not in the atom map")),
    ],
)
def test_get_map_index_error(raise_error, expected_raises):

    molecule = Molecule.from_smiles("C")

    with expected_raises:
        map_index = get_map_index(molecule, 0, error_on_missing=raise_error) == 5
        assert map_index == 0
