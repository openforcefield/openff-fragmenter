import importlib
from collections import defaultdict
from typing import TypeVar

import pytest
from openff.fragmenter.fragment import BondTuple
from openff.fragmenter.utils import get_map_index
from openff.toolkit.topology import Molecule

T = TypeVar("T")

try:
    importlib.import_module("openeye.oechem")

    has_openeye = True
except ImportError:
    has_openeye = False

using_openeye = pytest.mark.skipif(not has_openeye, reason="Cannot run without OpenEye")


def smarts_set_to_map_indices(
    input_set: set[str], molecule: Molecule
) -> set[int | BondTuple]:
    """A utility function which maps a set of the form ``{"SMARTS_1", "SMARTS_2", ...}``
    to one of the form ``{map_index_1, map_index_2, ...}`` when the SMARTS defines a
    single atom or ``{(map_index_1_a, map_index_1_b), ...}`` when the SMARTS defines a
    bond.

    The is useful for tests which need to define expected values for atoms / bonds
    with specific map indices in a way that is agnostic of how the molecule used in
    the test is ordered (i.e. whether it was canonically ordered with RDKit or OE.
    """

    return_value = set()

    for smarts in input_set:
        matches = list(
            {
                tuple(sorted(match))
                for match in molecule.chemical_environment_matches(smarts)
            }
        )
        assert len(matches) == 1

        match = matches[0]
        assert 0 < len(match) <= 2

        if len(match) == 1:
            return_value.add(get_map_index(molecule, match[0]))
        else:
            return_value.add(tuple(get_map_index(molecule, i) for i in match))

    return return_value


def key_smarts_to_map_indices(
    input_dictionary: dict[str, T], molecule: Molecule
) -> dict[int | BondTuple, T]:
    """A utility function which maps a dictionary of the form ``value["SMARTS"] = x``
    to one of the form ``value[map_index] = x`` when the SMARTS defines a single atom
    or ``value[(map_index_1, map_index_2)] = x`` when the SMARTS defines a bond.

    The is useful for tests which need to define expected values for atoms / bonds
    with specific map indices in a way that is agnostic of how the molecule used in
    the test is ordered (i.e. whether it was canonically ordered with RDKit or OE.
    """

    return_value = {}

    for smarts, expected in input_dictionary.items():
        matches = list(
            {
                tuple(sorted(match))
                for match in molecule.chemical_environment_matches(smarts)
            }
        )
        assert len(matches) == 1

        match = matches[0]
        assert 0 < len(match) <= 2

        if len(match) == 1:
            return_value[get_map_index(molecule, match[0])] = expected
        else:
            return_value[tuple(get_map_index(molecule, i) for i in match)] = expected

    return return_value


def value_smarts_to_map_indices(
    input_dictionary: dict[str, set[str]], molecule: Molecule
) -> dict[str, set[int | BondTuple]]:
    """A utility function which maps a dictionary of the form
    ``value[key] = {"SMARTS_1", "SMARTS_2", ...}`` to one of the form
    ``value[key] = {map_index_1, map_index_2, ...}`` when the SMARTS defines a single
    atom or ``value[key] = {(map_index_1_a, map_index_1_b), ...}`` when the SMARTS
    defines a bond.

    The is useful for tests which need to define expected values for atoms / bonds
    with specific map indices in a way that is agnostic of how the molecule used in
    the test is ordered (i.e. whether it was canonically ordered with RDKit or OE.
    """

    return_value = defaultdict(set)

    for key, values in input_dictionary.items():
        for smarts in values:
            matches = list(
                {
                    tuple(sorted(match))
                    for match in molecule.chemical_environment_matches(smarts)
                }
            )

            for match in matches:
                assert 0 < len(match) <= 2

                if len(match) == 1:
                    return_value[key].add(get_map_index(molecule, match[0]))
                else:
                    return_value[key].add(
                        tuple(get_map_index(molecule, i) for i in match)
                    )

    return {**return_value}
