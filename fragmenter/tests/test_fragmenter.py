"""
Unit and regression test for the fragmenter package.
"""

# Import package, test suite, and other packages as needed
import logging

import numpy
import pytest
from openff.toolkit.topology import Molecule

from fragmenter.chemi import smiles_to_molecule
from fragmenter.fragment import Fragmenter, PfizerFragmenter, WBOFragmenter
from fragmenter.tests.utils import using_openeye
from fragmenter.utils import get_fgroup_smarts, get_fgroup_smarts_comb, get_map_index


class DummyFragmenter(Fragmenter):
    """A mock fragmenter class used for testing"""

    def fragment(self):
        return


def test_keep_track_of_map():

    fragments = WBOFragmenter(Molecule.from_smiles("c1ccc(cc1)Nc2ncccn2"))
    fragments.fragment()

    assert all(
        "atom_map" in fragment.properties for fragment in fragments.fragments.values()
    )


def test_tag_fgroups(potanib):

    functional_groups = get_fgroup_smarts_comb()

    frags = DummyFragmenter(potanib, functional_groups=functional_groups)

    functional_groups = {
        "alkyne_0": {1, 2},
        "carbonyl_0": {21, 36},
        "amide_0": {35, 21, 36},
        "tri_halide_0": {29, 37, 38, 39},
    }

    for group, atoms in functional_groups.items():
        assert frags.functional_groups[group][0] == atoms

        assert len(frags.functional_groups[group][1]) > 0

        assert all(
            index in atoms
            for bond in frags.functional_groups[group][1]
            for index in bond
        )


def test_rotor_wbo(butane):

    fragmenter = WBOFragmenter(butane)
    assert fragmenter.rotors_wbo == {}

    fragmenter._get_rotor_wbo()

    assert list(fragmenter.rotors_wbo.keys()) == [(3, 4)]
    assert round(fragmenter.rotors_wbo[(3, 4)], ndigits=3) == 0.986


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


def test_wbo_fragment():
    """ Test build fragment"""

    fragmenter = WBOFragmenter(Molecule.from_smiles("CCCCC"))
    fragmenter.fragment()

    assert len(fragmenter.fragments) == len(fragmenter.rotors_wbo)
    assert fragmenter.fragments.keys() == fragmenter.rotors_wbo.keys()


def test_atom_bond_set_to_mol(abemaciclib):

    fragmenter = DummyFragmenter(abemaciclib, get_fgroup_smarts())

    atoms = {17, 18, 19, 20, 22, 26, 33, 34, 66, 67}
    bonds = {
        (17, 19),
        (17, 33),
        (18, 20),
        (18, 33),
        (19, 34),
        (20, 34),
        (22, 26),
        (26, 34),
        (26, 66),
        (26, 67),
    }

    fragment = fragmenter._atom_bond_set_to_mol(atoms=atoms, bonds=bonds)

    for bond in fragment.bonds:

        if bond.atom1.atomic_number == 1 or bond.atom2.atomic_number == 1:
            continue

        map_index_1 = get_map_index(fragment, bond.atom1_index)
        map_index_2 = get_map_index(fragment, bond.atom2_index)

        assert tuple(sorted((map_index_1, map_index_2))) in bonds


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


def test_compare_wbo(butane):

    fragmenter = WBOFragmenter(butane)
    fragmenter.calculate_wbo()
    fragmenter._get_rotor_wbo()

    assert numpy.isclose(
        fragmenter._compare_wbo(fragment=butane, bond_tuple=(3, 4)), 0.0
    )


@pytest.mark.parametrize(
    "input_smiles, n_output", [("CCCC", 0), ("c1ccccc1", 1), ("c1ccccc1C", 1)]
)
def test_find_ring_systems(input_smiles, n_output):

    fragmenter = WBOFragmenter(Molecule.from_smiles(input_smiles))
    fragmenter._find_ring_systems()

    assert len(fragmenter.ring_systems) == n_output


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


def test_find_ortho_substituents(dasatanib):

    # TODO: add test to test adding substituent directly bonded to rotatable bond
    fragmenter = WBOFragmenter(dasatanib)
    fragmenter._find_ring_systems(keep_non_rotor_ring_substituents=False)

    ortho_atoms, ortho_bonds = fragmenter._find_ortho_substituents(bonds={(6, 28)})

    assert ortho_atoms == {19, 28, 33}
    assert ortho_bonds == {(5, 19), (6, 28), (7, 33)}


def test_find_rotatable_bonds(abemaciclib):

    fragmenter = WBOFragmenter(abemaciclib)
    rotatable_bonds = fragmenter._find_rotatable_bonds()

    assert len(rotatable_bonds) == 7

    expected_rotatable_bonds = {
        (7, 13),
        (8, 25),
        (25, 33),
        (26, 34),
        (14, 35),
        (27, 32),
        (15, 35),
    }

    assert {*rotatable_bonds} == expected_rotatable_bonds


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


@pytest.mark.parametrize(
    "smiles, output",
    [
        (
            "OC1(CN(C1)C(=O)C1=C(NC2=C(F)C=C(I)C=C2)C(F)=C(F)C=C1)[C@@H]1CCCCN1",
            {20: "S"},
        ),
        (
            "OC1(CN(C1)C(=O)C1=C(NC2=C(F)C=C(I)C=C2)C(F)=C(F)C=C1)[C@H]1CCCCN1",
            {20: "R"},
        ),
        (r"C\C=C\C", {(1, 2): "E"}),
        (r"C/C=C\C", {(1, 2): "Z"}),
    ],
)
def test_find_stereo(smiles, output):

    fragmenter = WBOFragmenter(Molecule.from_smiles(smiles))
    fragmenter._find_stereo()

    assert fragmenter.stereo == output


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


@using_openeye
@pytest.mark.parametrize("input_smiles, n_output", [("CCCCCCC", 4)])
def test_pfizer_fragmenter(input_smiles, n_output):

    fragmenter = PfizerFragmenter(Molecule.from_smiles(input_smiles))
    fragmenter.fragment()

    assert len(fragmenter.fragments) == n_output


@using_openeye
@pytest.mark.parametrize("input_smiles, output", [("CCCCCCC", (14, 20))])
def test_get_torsion_quartet(input_smiles, output):

    fragmenter = PfizerFragmenter(Molecule.from_smiles(input_smiles))
    atoms, bonds = fragmenter._get_torsion_quartet((3, 5))

    # This also includes explicit hydrogen
    assert len(atoms) == output[0]
    assert len(bonds) == output[1]


@using_openeye
@pytest.mark.parametrize(
    "input_smiles, bond, output",
    [("CCCCCCC", (3, 5), True), ("CCCCc1ccccc1", (9, 10), False)],
)
def test_get_ring_and_fgroup(input_smiles, bond, output):

    fragmenter = PfizerFragmenter(Molecule.from_smiles(input_smiles))
    atoms, bonds = fragmenter._get_torsion_quartet(bond)

    bonds = {tuple(sorted(bond)) for bond in bonds}

    l_atoms = len(atoms)
    l_bonds = len(bonds)

    atoms_2, bonds_2 = fragmenter._get_ring_and_fgroups(atoms, bonds)

    assert (l_atoms == len(atoms_2)) == output
    assert (l_bonds == len(bonds_2)) == output


@using_openeye
@pytest.mark.parametrize(
    "input_smiles, bond, expected_atoms, expected_bonds",
    [
        (
            "[H:11][c:1]1[c:2]([c:4]([c:6]([c:5]([c:3]1[H:13])[C:7](=[O:10])[H:15])"
            "[C:9]([H:19])([H:20])[C:8]([H:16])([H:17])[H:18])[H:14])[H:12]",
            (6, 9),
            {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 14, 16, 17, 18, 19, 20},
            # fmt: off
            {
                (1, 2), (6, 9), (1, 3), (4, 6), (9, 20), (5, 6), (9, 19), (5, 7),
                (4, 14), (8, 9), (8, 18), (7, 10), (8, 17), (8, 16), (2, 4), (3, 5)
            }
        ),
        (
            "[H:11][c:1]1[c:2]([c:4]([c:6]([c:5]([c:3]1[H:13])[C:7](=[O:10])[H:15])"
            "[C:9]([H:19])([H:20])[C:8]([H:16])([H:17])[H:18])[H:14])[H:12]",
            (5, 7),
            {1, 2, 3, 4, 5, 6, 7, 9, 10, 13, 15},
            # fmt: off
            {
                (1, 2), (6, 9), (1, 3), (4, 6), (7, 15), (3, 13), (5, 6), (5, 7),
                (7, 10), (2, 4), (3, 5)
            }
        ),
        (
            "[H:19][c:1]1[c:3]([c:9]([c:15]([c:10]([c:4]1[H:22])[H:28])[c:17]2[c:13]"
            "([c:7]([c:8]([c:14]([c:18]2[c:16]3[c:11]([c:5]([c:2]([c:6]([c:12]3[H:30])"
            "[H:24])[H:20])[H:23])[H:29])[H:32])[H:26])[H:25])[H:31])[H:27])[H:21]",
            (15, 17),
            {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 27, 28, 31},
            # fmt: off
            {
                (11, 16), (1, 3), (10, 15), (13, 17), (12, 16), (2, 5), (8, 14),
                (9, 15), (14, 18), (16, 18), (4, 10), (2, 6), (5, 11), (1, 4), (3, 9),
                (6, 12), (10, 28), (15, 17), (17, 18), (9, 27), (7, 13), (13, 31),
                (7, 8)
            }
        ),
    ],
)
def test_get_ring_and_fgroup_ortho(input_smiles, bond, expected_atoms, expected_bonds):
    """Ensure that FGs and rings attached to ortho groups are correctly
    detected.

    The expected values were generated using fragmenter=0.0.7
    """

    fragmenter = PfizerFragmenter(Molecule.from_smiles(input_smiles))

    atoms, bonds = fragmenter._get_torsion_quartet(bond)
    atoms, bonds = fragmenter._get_ring_and_fgroups(atoms, bonds)

    assert atoms == expected_atoms
    assert bonds == expected_bonds


@pytest.mark.parametrize("input_smiles, bond, output", [("CNCCc1ccccc1", (6, 8), 7)])
def test_cap_open_valance(input_smiles, bond, output):

    # TODO: This test doesn't test anything?

    fragmenter = PfizerFragmenter(Molecule.from_smiles(input_smiles))

    atoms, bonds = fragmenter._get_torsion_quartet(bond)
    atoms, bonds = fragmenter._get_ring_and_fgroups(atoms, bonds)

    fragmenter._cap_open_valence(atoms, bonds)

    # Check that carbon bonded to N was added
    assert get_map_index(fragmenter.molecule, output)
