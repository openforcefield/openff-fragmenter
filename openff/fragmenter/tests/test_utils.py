from contextlib import nullcontext

import pytest
from openff.fragmenter.utils import (
    default_functional_groups,
    get_atom_index,
    get_map_index,
    global_toolkit_registry,
)
from openff.toolkit.topology import Molecule
from openff.toolkit.utils import (
    GLOBAL_TOOLKIT_REGISTRY,
    ToolkitRegistry,
    ToolkitWrapper,
)


def test_default_functional_groups():
    """Tests that the `default_functional_groups` utility returns correctly."""

    smarts = default_functional_groups()

    assert "hydrazine" in smarts
    assert smarts["hydrazine"] == "[NX3:1][NX3:2]"

    assert "phosphon" not in smarts

    assert len(smarts) == 21


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


def test_get_atom_index():
    molecule = Molecule.from_smiles("[C:5]([H:1])([H:2])([H:3])([H:4])")
    assert get_atom_index(molecule, 5) == 0


def test_global_toolkit_registry():
    class DummyToolkitWrapper(ToolkitWrapper):
        def from_smiles(self, *args, **kwargs):
            return type(self)

    original_toolkits = GLOBAL_TOOLKIT_REGISTRY.registered_toolkits

    with global_toolkit_registry(DummyToolkitWrapper()):
        return_value = Molecule.from_smiles("C")

    assert return_value == DummyToolkitWrapper

    # Make sure the registry is returned to it's previous state.
    assert isinstance(Molecule.from_smiles("C"), Molecule)

    assert all(
        type(original) is type(current)
        for original, current in zip(
            original_toolkits, GLOBAL_TOOLKIT_REGISTRY.registered_toolkits
        )
    )

    with global_toolkit_registry(ToolkitRegistry([DummyToolkitWrapper])):
        assert Molecule.from_smiles("C") == DummyToolkitWrapper
