"""
Unit and regression test for the fragmenter package.
"""

# Import package, test suite, and other packages as needed
import logging

import numpy
import pytest
from openff.toolkit.topology import Molecule

from fragmenter.chemi import smiles_to_molecule
from fragmenter.fragment import PfizerFragmenter, WBOFragmenter
from fragmenter.tests.utils import (
    key_smarts_to_map_indices,
    smarts_set_to_map_indices,
    value_smarts_to_map_indices,
)
from fragmenter.utils import get_atom_index, get_map_index


@pytest.mark.parametrize(
    "smiles, output",
    [
        (
            "OC1(CN(C1)C(=O)C1=C(NC2=C(F)C=C(I)C=C2)C(F)=C(F)C=C1)[C@@H]1CCCCN1",
            {"[#7r6]-[#6:1]-[#6!r6]": "S"},
        ),
        (
            "OC1(CN(C1)C(=O)C1=C(NC2=C(F)C=C(I)C=C2)C(F)=C(F)C=C1)[C@H]1CCCCN1",
            {"[#7r6]-[#6:1]-[#6!r6]": "R"},
        ),
        (r"C\C=C\C", {"[#6:1]=[#6:2]": "E"}),
        (r"C/C=C\C", {"[#6:1]=[#6:2]": "Z"}),
    ],
)
def test_find_stereo(smiles, output):

    fragmenter = WBOFragmenter(Molecule.from_smiles(smiles))
    fragmenter._find_stereo()

    expected_stereo = key_smarts_to_map_indices(output, fragmenter.molecule)
    assert fragmenter.stereo == expected_stereo


@pytest.mark.parametrize(
    "smiles, fragment_smiles, output, warning",
    [
        ("C[C@@](F)(Cl)I", "C[C@@](F)(Cl)I", True, None),
        ("C[C@@](F)(Cl)I", "C[C@](F)(Cl)I", False, "Stereochemistry for atom "),
        ("C[C@@](F)(Cl)F", "C[C@@](F)(Cl)I", False, "A new stereocenter formed at"),
        (r"C\C=C\C", r"C\C=C\C", True, None),
        (r"C/C=C\C", r"C\C=C\C", False, "Stereochemistry for bond"),
        ("CC=C(C)C", r"C\C=C(C)\CC", False, "A new chiral bond formed at"),
    ],
)
def test_check_stereo(smiles, fragment_smiles, output, warning, caplog):

    fragmenter = WBOFragmenter(Molecule.from_smiles(smiles))
    fragmenter._find_stereo()

    fragment = smiles_to_molecule(fragment_smiles, add_atom_map=True)

    with caplog.at_level(logging.WARNING):
        assert fragmenter._check_stereo(fragment) == output

    if warning is None:
        assert len(caplog.records) == 0
    else:
        assert len(caplog.records) == 1
        assert caplog.records[0].message.startswith(warning)


@pytest.mark.parametrize(
    "smiles, fragment_smiles",
    [
        ("C[C@@](F)(Cl)I", "C[C@@](F)(Cl)I"),
        ("C[C@@](F)(Cl)I", "C[C@](F)(Cl)I"),
        (r"C\C=C\C", r"C\C=C\C"),
        (r"C/C=C\C", r"C\C=C\C"),
    ],
)
def test_fix_stereo(smiles, fragment_smiles):

    fragmenter = WBOFragmenter(Molecule.from_smiles(smiles))
    fragmenter._find_stereo()

    fragment = smiles_to_molecule(fragment_smiles, add_atom_map=True)
    fixed = fragmenter._fix_stereo(fragment)

    assert fragmenter._check_stereo(fixed) is True


def test_tag_functional_groups(potanib):

    fragmenter = WBOFragmenter(
        potanib,
        functional_groups={
            # Taken from the old `fgroup_smarts_comb.yml` file.
            "amide": "[NX3R0:1][CX3:2](=[OX1:3])",
            "tri_halide": "[#6:1]([#9:2])([#9:3])([#9:4])",
            "carbonyl": "[CX3R0:1]=[OX1:2]",
            "alkyne": "[CX2:1]#[CX2:2]",
        },
    )

    expected_functional_groups = value_smarts_to_map_indices(
        {
            "alkyne_0": {"[#6r6]-[#6:1]#[#6]", "[#6r5]-[#6:1]#[#6]"},
            "carbonyl_0": {"[#6:1]=[#8]", "[#6]=[#8:1]"},
            "amide_0": {"[#6:1]=[#8]", "[#6]=[#8:1]", "[#7:1]-[#6]=[#8]"},
            "tri_halide_0": {"[#6:1]-[#9]", "[#6]-[#9:1]"},
        },
        fragmenter.molecule,
    )
    actual_functional_groups = {
        name: groups[0] for name, groups in fragmenter.functional_groups.items()
    }

    assert actual_functional_groups == expected_functional_groups

    for group in fragmenter.functional_groups:

        assert all(
            index in expected_functional_groups[group]
            for bond in fragmenter.functional_groups[group][1]
            for index in bond
        )


def test_find_rotatable_bonds(abemaciclib):

    fragmenter = WBOFragmenter(abemaciclib)
    rotatable_bonds = fragmenter._find_rotatable_bonds()

    assert len(rotatable_bonds) == 7

    expected_rotatable_bonds = smarts_set_to_map_indices(
        {
            "[#6ar6:1]-[#6ar6:2]",
            "[#6ar6:1]-[#6H2:2]",
            "[#7r6:1]-[#6H2:2]-[#6a]",
            "[#7r6:1]-[#6H2:2]-[#6H3]",
            "[#6ar6:1]-[#7H1:2]-[#6a](:[#7a])(:[#7a])",
            "[#7r5:1]-[#6:2](-[#6H3])(-[#6H3])",
            "[#6ar6:1]-[#7H1:2]-[#6a]:[#6a]",
        },
        fragmenter.molecule,
    )

    assert {*rotatable_bonds} == expected_rotatable_bonds


def test_atom_bond_set_to_mol(abemaciclib):

    fragmenter = WBOFragmenter(abemaciclib)

    atoms = {
        get_map_index(fragmenter.molecule, atom_index)
        for match in fragmenter.molecule.chemical_environment_matches(
            "[C:1][C:2][N:3]1[C:4][C:5][N:6][C:7][C:8]1"
        )
        for atom_index in match
    }

    bonds = {
        (
            get_map_index(fragmenter.molecule, bond.atom1_index),
            get_map_index(fragmenter.molecule, bond.atom2_index),
        )
        for bond in fragmenter.molecule.bonds
        if get_map_index(fragmenter.molecule, bond.atom1_index) in atoms
        and get_map_index(fragmenter.molecule, bond.atom2_index) in atoms
    }

    fragment = fragmenter._atom_bond_set_to_mol(atoms=atoms, bonds=bonds)

    for bond in fragment.bonds:

        if bond.atom1.atomic_number == 1 or bond.atom2.atomic_number == 1:
            continue

        map_index_1 = get_map_index(fragment, bond.atom1_index)
        map_index_2 = get_map_index(fragment, bond.atom2_index)

        assert tuple(sorted((map_index_1, map_index_2))) in bonds


@pytest.mark.parametrize(
    "input_smiles, expected_n_atoms, expected_n_bonds", [("CCCCCCC", 14, 20)]
)
def test_get_torsion_quartet(input_smiles, expected_n_atoms, expected_n_bonds):

    fragmenter = PfizerFragmenter(Molecule.from_smiles(input_smiles))

    bond_match = fragmenter.molecule.chemical_environment_matches("C[C:1][C:2]CCCC")[0]

    atoms, bonds = fragmenter._get_torsion_quartet(
        (
            get_map_index(fragmenter.molecule, bond_match[0]),
            get_map_index(fragmenter.molecule, bond_match[1]),
        )
    )

    # This also includes explicit hydrogen
    assert len(atoms) == expected_n_atoms
    assert len(bonds) == expected_n_bonds


@pytest.mark.parametrize(
    "input_smiles, n_ring_systems", [("CCCC", 0), ("c1ccccc1", 1), ("c1ccccc1C", 1)]
)
def test_find_ring_systems(input_smiles, n_ring_systems):

    fragmenter = WBOFragmenter(Molecule.from_smiles(input_smiles))
    fragmenter._find_ring_systems()

    assert len(fragmenter.ring_systems) == n_ring_systems


@pytest.mark.parametrize(
    "keep_non_rotor_ring_substituents, n_output", [(True, 7), (False, 6)]
)
def test_keep_non_rotor(keep_non_rotor_ring_substituents, n_output):

    fragmenter = WBOFragmenter(Molecule.from_smiles("c1ccccc1C"))
    fragmenter._find_ring_systems(
        keep_non_rotor_ring_substituents=keep_non_rotor_ring_substituents
    )

    assert len(fragmenter.ring_systems[1][0]) == n_output


@pytest.mark.parametrize(
    "input_smiles, bond_smarts, expected",
    [
        ("CCCCCCC", "C[C:1][C:2]CCCC", True),
        ("CCCCc1ccccc1", "[#6H3]-[#6H2:1]-[#6H2:2]", False),
    ],
)
def test_get_ring_and_fgroup(input_smiles, bond_smarts, expected):

    fragmenter = PfizerFragmenter(Molecule.from_smiles(input_smiles))
    # noinspection PyTypeChecker
    atoms, bonds = fragmenter._get_torsion_quartet(
        tuple(
            get_map_index(fragmenter.molecule, i)
            for i in fragmenter.molecule.chemical_environment_matches(bond_smarts)[0]
        )
    )

    bonds = {tuple(sorted(bond)) for bond in bonds}

    l_atoms = len(atoms)
    l_bonds = len(bonds)

    atoms_2, bonds_2 = fragmenter._get_ring_and_fgroups(atoms, bonds)

    assert (l_atoms == len(atoms_2)) == expected
    assert (l_bonds == len(bonds_2)) == expected


@pytest.mark.parametrize(
    "input_smiles, bond_smarts, expected_pattern",
    [
        (
            "CCc1ccccc1C=O",
            "C[CH2:2][c:1]1ccccc1C=O",
            "[CH3:8][CH2:9][c:6]1[cH:4][cH:2][cH:1][cH:3][c:5]1[CH:7]=[O:10]",
        ),
        (
            "c1ccc(cc1)c2ccccc2c3ccccc3",
            "c1cc[c:1](cc1)[c:2]2ccccc2c3ccccc3",
            "[cH:1]1[cH:3][cH:9][c:15]([cH:10][cH:4]1)[c:17]2[cH:13][cH:7][cH:8]"
            "[cH:14][c:18]2[c:16]3[cH:11][cH:5][cH:2][cH:6][cH:12]3",
        ),
    ],
)
def test_get_ring_and_fgroup_ortho(input_smiles, bond_smarts, expected_pattern):
    """Ensure that FGs and rings attached to ortho groups are correctly
    detected.

    The expected values were generated using fragmenter=0.0.7
    """

    fragmenter = PfizerFragmenter(Molecule.from_smiles(input_smiles))

    bond = tuple(
        get_map_index(fragmenter.molecule, i)
        for i in fragmenter.molecule.chemical_environment_matches(bond_smarts)[0]
    )

    atoms, bonds = fragmenter._get_torsion_quartet(bond)
    atoms, bonds = fragmenter._get_ring_and_fgroups(atoms, bonds)

    actual_atoms = {
        map_index
        for map_index in atoms
        if fragmenter.molecule.atoms[
            get_atom_index(fragmenter.molecule, map_index)
        ].atomic_number
        != 1
    }
    expected_atoms = {
        get_map_index(fragmenter.molecule, atom_index)
        for match in fragmenter.molecule.chemical_environment_matches(expected_pattern)
        for atom_index in match
    }

    assert actual_atoms == expected_atoms


def test_find_ortho_substituents(dasatanib):

    # TODO: add test to test adding substituent directly bonded to rotatable bond
    fragmenter = WBOFragmenter(dasatanib)
    fragmenter._find_ring_systems(keep_non_rotor_ring_substituents=False)

    ortho_atoms, ortho_bonds = fragmenter._find_ortho_substituents(
        bonds=smarts_set_to_map_indices(
            {"[#6ar6:1]-[#7:2]-[#6]=[#8]"}, fragmenter.molecule
        )
    )

    assert ortho_atoms == smarts_set_to_map_indices(
        {"[#6ar6]-[#7:1]-[#6]=[#8]", "[#6ar6]-[#17:1]", "[#6H3:1]-[#6r6]:[#6r6]"},
        fragmenter.molecule,
    )
    assert ortho_bonds == smarts_set_to_map_indices(
        {"[#6ar6:1]-[#7:2]-[#6]=[#8]", "[#6ar6:1]-[#17:2]", "[#6H3:2]-[#6r6:1]:[#6r6]"},
        fragmenter.molecule,
    )


def test_cap_open_valance():

    fragmenter = WBOFragmenter(Molecule.from_smiles("CNCCc1ccccc1"))

    expected_atom = get_map_index(
        fragmenter.molecule,
        fragmenter.molecule.chemical_environment_matches("[#7]-[#6H3:1]")[0][0],
    )

    # noinspection PyTypeChecker
    atoms, bonds = fragmenter._get_torsion_quartet(
        tuple(
            get_map_index(fragmenter.molecule, i)
            for i in fragmenter.molecule.chemical_environment_matches(
                "[#6a:1]-[#6H2:2]"
            )[0]
        )
    )
    atoms, bonds = fragmenter._get_ring_and_fgroups(atoms, bonds)

    # Remove the cap atom from the current list to make sure it gets included during
    # capping.
    atoms -= {expected_atom}

    atoms, _ = fragmenter._cap_open_valence(atoms, bonds)

    # Check that carbon bonded to N was added
    assert expected_atom in atoms


# def test_to_qcscheme_mol():
#
#     fragmenter = WBOFragmenter(Molecule.from_smiles("CCCCCC"))
#     fragmenter.fragment()
#
#     qcschema_mol = fragmenter._to_qcschema_mol(fragmenter.fragments[(3, 5)])
#
#     assert "initial_molecule" in qcschema_mol
#     assert "geometry" in qcschema_mol["initial_molecule"][0]
#     assert "symbols" in qcschema_mol["initial_molecule"][0]
#     assert "connectivity" in qcschema_mol["initial_molecule"][0]
#     assert "identifiers" in qcschema_mol
#     assert "provenance" in qcschema_mol
#
#
# def test_td_inputs():
#
#     fragmenter = WBOFragmenter(Molecule.from_smiles("CCCCCC"))
#     fragmenter.fragment()
#
#     td_inputs = fragmenter.to_torsiondrive_json()
#     assert len(td_inputs) == 2


def test_wbo_fragment():
    """ Test build fragment"""

    fragmenter = WBOFragmenter(Molecule.from_smiles("CCCCC"))
    fragmenter.fragment()

    assert len(fragmenter.fragments) == len(fragmenter.rotors_wbo)
    assert fragmenter.fragments.keys() == fragmenter.rotors_wbo.keys()


def test_keep_track_of_map():

    fragments = WBOFragmenter(Molecule.from_smiles("CCCCC"))
    fragments.fragment()

    assert all(
        "atom_map" in fragment.properties for fragment in fragments.fragments.values()
    )


def test_calculate_wbo():

    fragmenter = WBOFragmenter(Molecule.from_smiles("CCCC"))

    molecule = fragmenter.calculate_wbo()
    assert not molecule

    for bond in fragmenter.molecule.bonds:
        assert bond.fractional_bond_order is not None

    molecule = fragmenter.calculate_wbo(fragmenter.molecule)
    assert molecule

    for bond in molecule.bonds:
        assert bond.fractional_bond_order is not None


def test_get_rotor_wbo():

    fragmenter = WBOFragmenter(Molecule.from_smiles("CCCC"))
    assert fragmenter.rotors_wbo == {}

    expected_bonds = {
        (
            get_map_index(fragmenter.molecule, match[0]),
            get_map_index(fragmenter.molecule, match[1]),
        )
        for match in fragmenter.molecule.chemical_environment_matches("[#6:1]-[#6:2]")
    }

    fragmenter._get_rotor_wbo()

    assert len(fragmenter.rotors_wbo) == 1

    rotor_index = next(iter(fragmenter.rotors_wbo))

    assert rotor_index in expected_bonds
    assert numpy.isclose(fragmenter.rotors_wbo[rotor_index], 0.986, atol=0.001)


def test_compare_wbo():

    fragmenter = WBOFragmenter(Molecule.from_smiles("CCCC"))
    fragmenter.calculate_wbo()
    fragmenter._get_rotor_wbo()

    fragment = smiles_to_molecule("CCCC", add_atom_map=True)

    for bond_tuple in fragmenter.rotors_wbo:

        assert numpy.isclose(
            fragmenter._compare_wbo(fragment=fragment, bond_tuple=bond_tuple),
            0.0,
            atol=1.0e-6,
        )


def test_build_fragment():

    fragmenter = WBOFragmenter(Molecule.from_smiles("CCCCCC"))
    fragmenter.calculate_wbo()
    fragmenter._get_rotor_wbo()

    fragmenter.threshold = 0.05

    for bond in fragmenter.rotors_wbo:
        fragmenter._build_fragment(bond)

    assert len(fragmenter.fragments) == 3

    assert (
        fragmenter.fragments[(3, 5)].to_smiles(explicit_hydrogens=False, mapped=False)
        == "CCCCC"
    )
    assert (
        fragmenter.fragments[(4, 6)].to_smiles(explicit_hydrogens=False, mapped=False)
        == "CCCCC"
    )
    assert (
        fragmenter.fragments[(5, 6)].to_smiles(explicit_hydrogens=False, mapped=False)
        == "CCCCCC"
    )


@pytest.mark.parametrize(
    "input_smiles, n_output",
    [
        (
            "[H:42][N:19]([C@H:11]1[C@@H:13]2[N:18]([C:7]1=[O:22])[CH2:12][C:14]"
            "([S:26]2)([CH3:15])[CH3:16])[C:9](=[O:24])[CH3:17]",
            8,
        )
    ],
)
def test_ring_fgroups(input_smiles, n_output):

    fragmenter = WBOFragmenter(Molecule.from_smiles(input_smiles))
    fragmenter._find_ring_systems()

    assert len(fragmenter.ring_systems[1][0]) == n_output


def test_add_substituent():

    fragmenter = WBOFragmenter(Molecule.from_smiles("CCCCCC"))
    fragmenter.fragment()

    assert (
        fragmenter.fragments[(3, 5)].to_smiles(mapped=False, explicit_hydrogens=False)
        == "CCCCC"
    )

    fragment = fragmenter.fragments[(3, 5)]

    atoms = set(
        get_map_index(fragment, i)
        for i in range(fragment.n_atoms)
        if fragment.atoms[i].atomic_number != 1
    )

    bonds = set(
        (
            get_map_index(fragment, bond.atom1_index),
            get_map_index(fragment, bond.atom2_index),
        )
        for bond in fragment.bonds
        if bond.atom1.atomic_number != 1 and bond.atom2.atomic_number != 1
    )

    fragment = fragmenter._add_next_substituent(
        fragmenter.molecule, atoms, bonds, target_bond=(3, 5)
    )

    assert fragment.to_smiles(mapped=False, explicit_hydrogens=False) == "CCCCCC"


@pytest.mark.parametrize("input_smiles, n_output", [("CCCCCCC", 4)])
def test_pfizer_fragmenter(input_smiles, n_output):

    fragmenter = PfizerFragmenter(Molecule.from_smiles(input_smiles))
    fragmenter.fragment()

    assert len(fragmenter.fragments) == n_output
