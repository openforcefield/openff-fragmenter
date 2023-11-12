"""Test chemi module"""
import numpy
import pytest
from openff.fragmenter import chemi
from openff.fragmenter._tests.utils import using_openeye
from openff.fragmenter.chemi import (
    _extract_oe_fragment,
    _extract_rd_fragment,
    _find_oe_stereocenters,
    _find_rd_stereocenters,
    assign_elf10_am1_bond_orders,
    extract_fragment,
    find_ring_systems,
    find_stereocenters,
    smiles_to_molecule,
)
from openff.fragmenter.utils import global_toolkit_registry
from openff.toolkit.topology import Molecule
from openff.toolkit.utils import (
    AmberToolsToolkitWrapper,
    OpenEyeToolkitWrapper,
    ToolkitRegistry,
)
from openff.units import unit
from openff.utilities import MissingOptionalDependencyError


def test_assign_elf10_am1_bond_orders():
    molecule = assign_elf10_am1_bond_orders(Molecule.from_smiles("CCCC"))

    for bond in molecule.bonds:
        assert bond.fractional_bond_order is not None


@using_openeye
def test_assign_elf10_am1_bond_orders_simple_parity():
    with global_toolkit_registry(OpenEyeToolkitWrapper()):
        molecule = assign_elf10_am1_bond_orders(Molecule.from_smiles("C"))
        oe_bond_orders = [bond.fractional_bond_order for bond in molecule.bonds]

    with global_toolkit_registry(
        ToolkitRegistry([AmberToolsToolkitWrapper(), OpenEyeToolkitWrapper()])
    ):
        molecule = assign_elf10_am1_bond_orders(Molecule.from_smiles("C"))
        at_bond_orders = [bond.fractional_bond_order for bond in molecule.bonds]

    assert numpy.allclose(oe_bond_orders, at_bond_orders, atol=1.0e-2)


@pytest.mark.parametrize("smiles, max_confs", [("CCCCCCC", 1), ("CCCCCCC", 3)])
def test_generate_conformers(smiles, max_confs):
    returned_molecule = chemi._generate_conformers(
        Molecule.from_smiles(smiles), max_confs=max_confs
    )

    assert 0 < returned_molecule.n_conformers <= max_confs


def test_generate_conformers_ordering():
    original_molecule = Molecule.from_smiles("CCCC")

    returned_molecule = chemi._generate_conformers(original_molecule, max_confs=1)
    assert returned_molecule.n_conformers == 1

    # Make sure the atom ordering did not change.
    _, atom_map = Molecule.are_isomorphic(
        original_molecule, returned_molecule, return_atom_map=True
    )

    assert all(i == j for i, j in atom_map.items())


def test_generate_conformers_canonical_check():
    original_molecule = Molecule.from_smiles("CCCCl").canonical_order_atoms()
    original_molecule = chemi._generate_conformers(original_molecule, max_confs=1)

    # Generate a conformer using a molecule with permuted atom orderings.
    atom_map = numpy.arange(original_molecule.n_atoms)
    numpy.random.shuffle(atom_map)

    remapped_molecule = original_molecule.remap(
        {int(i): int(j) for i, j in enumerate(atom_map)}
    )
    remapped_molecule = chemi._generate_conformers(remapped_molecule, max_confs=1)

    original_conformer = original_molecule.conformers[0].m_as(unit.angstrom)
    remapped_conformer = remapped_molecule.conformers[0].m_as(unit.angstrom)

    assert numpy.allclose(original_conformer[:4], remapped_conformer[atom_map][:4])


@using_openeye
@pytest.mark.parametrize(
    "smiles",
    [
        "c1ccccc1",
        "C1CC2CCC1C2",
        "C12C3C4C1C5C2C3C45",
        "c1ccc2ccccc2c1",
        "c1ccc(cc1)C2CCCC2",
        "c1ccc(cc1)C2CCC3C2CC3",
        "c1cc2ccc3cccc4c3c2c(c1)cc4",
        "C1CCC2(CC1)CCCCC2",
        "c1ccc(cc1)Cc2ccccc2",
    ],
)
def test_find_ring_systems(smiles):
    from openeye import oechem

    molecule = Molecule.from_smiles(smiles, allow_undefined_stereo=True)

    oe_molecule = molecule.to_openeye()
    _, expected = oechem.OEDetermineRingSystems(oe_molecule)

    ring_systems = find_ring_systems(molecule)

    for i, ring_system in enumerate(expected):
        if ring_system > 0:
            assert ring_systems[i] == ring_system

        else:
            assert i not in ring_systems


@pytest.mark.parametrize("add_atom_map", [True, False])
def test_smiles_to_molecule(add_atom_map):
    molecule = smiles_to_molecule("CCCC", add_atom_map=add_atom_map)

    assert isinstance(molecule, Molecule)
    assert ("atom_map" in molecule.properties) == add_atom_map


@pytest.mark.parametrize(
    "smiles, expected_atoms, expected_bonds",
    [
        ("[C:1]([H:6])([H:7])([H:8])[C:5]([Br:3])([Cl:4])([H:2])", [4], []),
        ("[C:3]([Cl:1])([H:2])=[C:6]([H:4])([Cl:5])", [], [(2, 5)]),
        (
            "[C:3]([Cl:1])([H:2])=[C:6]([H:4])[C:5]([Br:7])([Cl:8])([H:9])",
            [4],
            [(2, 5)],
        ),
    ],
)
@pytest.mark.parametrize(
    "find_method", [_find_oe_stereocenters, _find_rd_stereocenters, find_stereocenters]
)
def test_find_stereocenters(smiles, expected_atoms, expected_bonds, find_method):
    molecule = Molecule.from_mapped_smiles(smiles, allow_undefined_stereo=True)

    stereogenic_atoms, stereogenic_bonds = None, None

    try:
        stereogenic_atoms, stereogenic_bonds = find_method(molecule)
    except (ModuleNotFoundError, MissingOptionalDependencyError) as e:
        pytest.skip(str(e))

    assert stereogenic_atoms == expected_atoms
    assert stereogenic_bonds == expected_bonds


@pytest.mark.parametrize(
    "smiles, atoms, bonds, expected",
    [
        (
            "[C:1]([H:4])([H:5])([H:6])"
            "[C:2]([H:7])([H:8])"
            "[C:3]([H:9])([H:10])([H:11])",
            {1, 2},
            {(1, 2)},
            "CC",
        ),
        (
            "[C:1]([H:4])([H:5])([H:6])"
            "[C:2]([H:7])([H:8])"
            "[C:3]([H:9])([H:10])([H:11])",
            {1, 2, 4},
            {(1, 2), (1, 4)},
            "CC",
        ),
        (
            r"[H:6]/[C:1](=[C:2](\[C:3]([H:7])([H:8])[H:9])/[Cl:5])/[Cl:4]",
            {1, 2, 4, 5},
            {(1, 2), (1, 4), (2, 5)},
            r"Cl\C=C/Cl",
        ),
        (
            "[H:7][C:1]([H:8])([C@:2]([F:3])([Cl:5])[Br:6])[Cl:4]",
            {1, 2, 3, 5, 6},
            {(1, 2), (2, 3), (2, 5), (2, 6)},
            r"C[C@](F)(Cl)Br",
        ),
    ],
)
@pytest.mark.parametrize(
    "extract_method", [_extract_rd_fragment, _extract_oe_fragment, extract_fragment]
)
def test_extract_fragment(smiles, atoms, bonds, expected, extract_method):
    molecule = Molecule.from_mapped_smiles(smiles)
    molecule.properties["atom_map"] = {i: i + 1 for i in range(molecule.n_atoms)}

    fragment = None

    try:
        fragment = extract_method(molecule, atoms, bonds)
    except (ModuleNotFoundError, MissingOptionalDependencyError) as e:
        pytest.skip(str(e))

    expected_fragment = Molecule.from_smiles(expected)

    assert Molecule.are_isomorphic(
        fragment, expected_fragment, bond_stereochemistry_matching=False
    )[0]


def test_extract_fragment_bonds_in_atoms():
    """Tests that an exception is raised when the bonds set contains atoms not in the
    atoms set."""

    molecule = Molecule.from_smiles("[H:1][C:2]#[C:3][H:4]")

    with pytest.raises(ValueError, match="set includes atoms not in the"):
        extract_fragment(molecule, {2}, {(2, 3)})


@using_openeye
def test_extract_fragment_disconnected_fragment_warning():
    molecule = Molecule.from_smiles(
        "[C:1]([H:3])([H:4])([H:5])[C:2]([H:6])([H:7])([H:8])"
    )

    with pytest.raises(AssertionError, match="An atom that is not bonded"):
        extract_fragment(molecule, {1, 2}, set())
