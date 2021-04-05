"""Test chemi module"""

import pytest
from openff.toolkit.topology import Molecule
from openff.toolkit.utils import OpenEyeToolkitWrapper, RDKitToolkitWrapper

from fragmenter import chemi
from fragmenter.chemi import find_ring_systems
from fragmenter.tests.utils import global_toolkit_wrapper


@pytest.mark.parametrize(
    "keep_confs, output", [(None, 1), (1, 1), (2, 2), (-1, 3), (-2, "error")]
)
def test_get_charges(keep_confs, output):
    mol = chemi.smiles_to_oemol("CCCCCCC")
    if output == "error":
        with pytest.raises(ValueError):
            chemi.get_charges(mol, keep_confs=keep_confs)
    else:
        charged = chemi.get_charges(mol, keep_confs=keep_confs)
        for i, c in enumerate(charged.GetConfs()):
            i += 1
        assert i == output


def test_generate_conformers():
    mol = chemi.smiles_to_oemol("CCCCCCC")
    confs = chemi.generate_conformers(mol, max_confs=1)
    assert confs.GetMaxConfIdx() == 1

    confs = chemi.generate_conformers(mol)
    assert confs.GetMaxConfIdx() == 3


def test_smiles_to_oemol():
    from openeye import oechem

    mol = chemi.smiles_to_oemol("CCCC")
    assert isinstance(mol, oechem.OEMol)
    assert oechem.OEMolToSmiles(mol) == "CCCC"
    assert mol.GetTitle() == "butane"

    mol = chemi.smiles_to_oemol("CCCC", normalize=False)
    assert mol.GetTitle() == ""

    mol = chemi.smiles_to_oemol("CCCC", add_atom_map=True)
    assert (
        oechem.OEMolToSmiles(mol)
        == "[H:5][C:1]([H:6])([H:7])[C:3]([H:11])([H:12])[C:4]([H:13])([H:14])[C:2]([H:8])([H:9])[H:10]"
    )


def test_normalize_molecule():
    from openeye import oechem

    mol = oechem.OEMol()
    oechem.OESmilesToMol(mol, "CCCC")
    assert mol.GetTitle() == ""

    normalized_mol = chemi.normalize_molecule(mol)
    assert normalized_mol.GetTitle() == "butane"


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
