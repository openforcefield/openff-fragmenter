"""Test chemi module"""
import numpy
import pytest
from openff.toolkit.topology import Molecule
from openff.toolkit.utils import (
    AmberToolsToolkitWrapper,
    OpenEyeToolkitWrapper,
    RDKitToolkitWrapper,
    ToolkitRegistry,
)
from simtk import unit

from fragmenter import chemi
from fragmenter.chemi import assign_elf10_am1_bond_orders, find_ring_systems
from fragmenter.tests.utils import global_toolkit_wrapper


@pytest.mark.parametrize(
    "toolkit_wrapper",
    [
        OpenEyeToolkitWrapper(),
        ToolkitRegistry([AmberToolsToolkitWrapper(), RDKitToolkitWrapper()]),
    ],
)
def test_assign_elf10_am1_bond_orders(toolkit_wrapper):

    with global_toolkit_wrapper(toolkit_wrapper):

        oe_molecule = Molecule.from_smiles("CCCC").to_openeye()
        oe_molecule = assign_elf10_am1_bond_orders(oe_molecule)

    for bond in oe_molecule.GetBonds():
        assert "WibergBondOrder" in bond.GetData()


def test_assign_elf10_am1_bond_orders_simple_parity():

    with global_toolkit_wrapper(OpenEyeToolkitWrapper()):

        oe_molecule = Molecule.from_smiles("C").to_openeye()
        oe_molecule = assign_elf10_am1_bond_orders(oe_molecule)

        oe_orders = [bond.GetData("WibergBondOrder") for bond in oe_molecule.GetBonds()]

    with global_toolkit_wrapper(
        ToolkitRegistry([AmberToolsToolkitWrapper(), OpenEyeToolkitWrapper()])
    ):

        oe_molecule = Molecule.from_smiles("C").to_openeye()
        oe_molecule = assign_elf10_am1_bond_orders(oe_molecule)

        at_orders = [bond.GetData("WibergBondOrder") for bond in oe_molecule.GetBonds()]

    assert numpy.allclose(oe_orders, at_orders, atol=1.0e-2)


@pytest.mark.parametrize("smiles, max_confs", [("CCCCCCC", 1), ("CCCCCCC", 3)])
@pytest.mark.parametrize(
    "toolkit_wrapper", [OpenEyeToolkitWrapper(), RDKitToolkitWrapper()]
)
def test_generate_conformers(smiles, max_confs, toolkit_wrapper):

    with global_toolkit_wrapper(toolkit_wrapper):

        returned_molecule = chemi.generate_conformers(
            Molecule.from_smiles(smiles), max_confs=max_confs
        )

    assert 0 < returned_molecule.n_conformers <= max_confs


@pytest.mark.parametrize(
    "toolkit_wrapper", [OpenEyeToolkitWrapper(), RDKitToolkitWrapper()]
)
def test_generate_conformers_ordering(toolkit_wrapper):

    with global_toolkit_wrapper(toolkit_wrapper):

        original_molecule = Molecule.from_smiles("CCCC")

        returned_molecule = chemi.generate_conformers(original_molecule, max_confs=1)
        assert returned_molecule.n_conformers == 1

        # Make sure the atom ordering did not change.
        _, atom_map = Molecule.are_isomorphic(
            original_molecule, returned_molecule, return_atom_map=True
        )

        assert all(i == j for i, j in atom_map.items())


@pytest.mark.parametrize(
    "toolkit_wrapper", [OpenEyeToolkitWrapper(), RDKitToolkitWrapper()]
)
def test_generate_conformers_canonical_check(toolkit_wrapper):

    with global_toolkit_wrapper(toolkit_wrapper):

        original_molecule = Molecule.from_smiles("CCCCl").canonical_order_atoms()
        original_molecule = chemi.generate_conformers(original_molecule, max_confs=1)

        # Generate a conformer using a molecule with permuted atom orderings.
        atom_map = numpy.arange(original_molecule.n_atoms)
        numpy.random.shuffle(atom_map)

        remapped_molecule = original_molecule.remap(
            {int(i): int(j) for i, j in enumerate(atom_map)}
        )
        remapped_molecule = chemi.generate_conformers(remapped_molecule, max_confs=1)

        original_conformer = original_molecule.conformers[0].value_in_unit(
            unit.angstrom
        )
        remapped_conformer = remapped_molecule.conformers[0].value_in_unit(
            unit.angstrom
        )

        assert numpy.allclose(original_conformer[:4], remapped_conformer[atom_map][:4])


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
@pytest.mark.parametrize(
    "toolkit_wrapper", [OpenEyeToolkitWrapper(), RDKitToolkitWrapper()]
)
def test_find_ring_systems(smiles, toolkit_wrapper):

    from openeye import oechem

    molecule = Molecule.from_smiles(
        smiles, allow_undefined_stereo=True, toolkit_registry=toolkit_wrapper
    )

    oe_molecule = molecule.to_openeye()
    _, expected = oechem.OEDetermineRingSystems(oe_molecule)

    with global_toolkit_wrapper(toolkit_wrapper):
        ring_systems = find_ring_systems(molecule)

    for i, ring_system in enumerate(expected):

        if ring_system > 0:
            assert ring_systems[i] == ring_system

        else:
            assert i not in ring_systems


@pytest.mark.parametrize(
    "add_atom_map, expected",
    [
        (False, "CCCC"),
        (
            True,
            "[H:5][C:1]([H:6])([H:7])[C:3]([H:11])([H:12])"
            "[C:4]([H:13])([H:14])[C:2]([H:8])([H:9])[H:10]",
        ),
    ],
)
@pytest.mark.parametrize(
    "toolkit_wrapper", [OpenEyeToolkitWrapper(), RDKitToolkitWrapper()]
)
def test_smiles_to_oemol(add_atom_map, expected, toolkit_wrapper):

    from openeye import oechem

    with global_toolkit_wrapper(toolkit_wrapper):
        mol = chemi.smiles_to_oemol("CCCC", add_atom_map=add_atom_map)

    assert isinstance(mol, oechem.OEMol)
    assert oechem.OEMolToSmiles(mol) == expected
