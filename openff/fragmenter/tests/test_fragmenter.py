"""
Unit and regression test for the fragmenter.fragmenter module.
"""

import logging

import numpy
import pytest
from openff.fragmenter.chemi import assign_elf10_am1_bond_orders, smiles_to_molecule
from openff.fragmenter.fragment import (
    FragmentationResult,
    Fragmenter,
    PfizerFragmenter,
    WBOFragmenter,
)
from openff.fragmenter.tests import does_not_raise
from openff.fragmenter.tests.utils import (
    key_smarts_to_map_indices,
    smarts_set_to_map_indices,
    value_smarts_to_map_indices,
)
from openff.fragmenter.utils import (
    default_functional_groups,
    get_atom_index,
    get_map_index,
)
from openff.toolkit.topology import Molecule
from openff.toolkit.utils import (
    GLOBAL_TOOLKIT_REGISTRY,
    RDKitToolkitWrapper,
    ToolkitRegistry,
)


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

    molecule = smiles_to_molecule(smiles, add_atom_map=True)

    actual_stereo = Fragmenter._find_stereo(molecule)
    expected_stereo = key_smarts_to_map_indices(output, molecule)

    assert actual_stereo == expected_stereo


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

    parent = smiles_to_molecule(smiles, add_atom_map=True)
    fragment = smiles_to_molecule(fragment_smiles, add_atom_map=True)

    parent_stereo = Fragmenter._find_stereo(parent)

    with caplog.at_level(logging.WARNING):
        assert Fragmenter._check_stereo(fragment, parent_stereo) == output

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

    parent_stereo = Fragmenter._find_stereo(smiles_to_molecule(smiles, True))

    fragment = smiles_to_molecule(fragment_smiles, add_atom_map=True)
    fixed = Fragmenter._fix_stereo(fragment, parent_stereo)

    assert Fragmenter._check_stereo(fixed, parent_stereo) is True


def test_find_functional_groups(potanib):

    found_groups = Fragmenter._find_functional_groups(
        molecule=potanib,
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
        potanib,
    )
    actual_functional_groups = {
        name: groups[0] for name, groups in found_groups.items()
    }

    assert actual_functional_groups == expected_functional_groups

    for group in found_groups:

        assert all(
            index in expected_functional_groups[group]
            for bond in found_groups[group][1]
            for index in bond
        )


def test_find_rotatable_bonds_default(abemaciclib):

    rotatable_bonds = Fragmenter.find_rotatable_bonds(abemaciclib, None)
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
        abemaciclib,
    )

    assert {*rotatable_bonds} == expected_rotatable_bonds


@pytest.mark.parametrize(
    "smarts, expected_raises",
    [
        (["[#6:1]-[#6:2]"], does_not_raise()),
        (
            ["[#6]-[#6]"],
            pytest.raises(ValueError, match="The `target_bond_smarts` pattern "),
        ),
    ],
)
def test_find_rotatable_bonds_custom(smarts, expected_raises):

    ethane = Molecule.from_smiles("CC")
    ethane.properties["atom_map"] = {i: i + 1 for i in range(ethane.n_atoms)}

    with expected_raises:

        rotatable_bonds = Fragmenter.find_rotatable_bonds(ethane, smarts)
        assert len(rotatable_bonds) == 1

        assert rotatable_bonds == [(1, 2)]


def test_atom_bond_set_to_mol(abemaciclib):

    molecule = smiles_to_molecule(abemaciclib.to_smiles(mapped=False), True)

    atoms = {
        get_map_index(molecule, atom_index)
        for match in molecule.chemical_environment_matches(
            "[C:1][C:2][N:3]1[C:4][C:5][N:6][C:7][C:8]1"
        )
        for atom_index in match
    }

    bonds = {
        (
            get_map_index(molecule, bond.atom1_index),
            get_map_index(molecule, bond.atom2_index),
        )
        for bond in molecule.bonds
        if get_map_index(molecule, bond.atom1_index) in atoms
        and get_map_index(molecule, bond.atom2_index) in atoms
    }

    fragment = Fragmenter._atom_bond_set_to_mol(molecule, {}, atoms=atoms, bonds=bonds)

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

    molecule = smiles_to_molecule(input_smiles, True)

    bond_match = molecule.chemical_environment_matches("C[C:1][C:2]CCCC")[0]

    atoms, bonds = Fragmenter._get_torsion_quartet(
        molecule,
        (
            get_map_index(molecule, bond_match[0]),
            get_map_index(molecule, bond_match[1]),
        ),
    )

    # This also includes explicit hydrogen
    assert len(atoms) == expected_n_atoms
    assert len(bonds) == expected_n_bonds


@pytest.mark.parametrize(
    "input_smiles, n_ring_systems", [("CCCC", 0), ("c1ccccc1", 1), ("c1ccccc1C", 1)]
)
def test_find_ring_systems(input_smiles, n_ring_systems):

    ring_systems = Fragmenter._find_ring_systems(
        smiles_to_molecule(input_smiles, True), {}
    )

    assert len(ring_systems) == n_ring_systems


@pytest.mark.parametrize(
    "keep_non_rotor_ring_substituents, n_output", [(True, 7), (False, 6)]
)
def test_keep_non_rotor(keep_non_rotor_ring_substituents, n_output):

    ring_systems = Fragmenter._find_ring_systems(
        smiles_to_molecule("c1ccccc1C", True),
        {},
        keep_non_rotor_ring_substituents=keep_non_rotor_ring_substituents,
    )

    assert len(ring_systems[1][0]) == n_output


@pytest.mark.parametrize(
    "input_smiles, bond_smarts, expected",
    [
        ("CCCCCCC", "C[C:1][C:2]CCCC", True),
        ("CCCCc1ccccc1", "[#6H3]-[#6H2:1]-[#6H2:2]", False),
    ],
)
def test_get_ring_and_fgroup(input_smiles, bond_smarts, expected):

    molecule, _, functional_groups, ring_systems = Fragmenter._prepare_molecule(
        smiles_to_molecule(input_smiles, True), default_functional_groups(), False
    )

    # noinspection PyTypeChecker
    atoms, bonds = Fragmenter._get_torsion_quartet(
        molecule,
        tuple(
            get_map_index(molecule, i)
            for i in molecule.chemical_environment_matches(bond_smarts)[0]
        ),
    )

    bonds = {tuple(sorted(bond)) for bond in bonds}

    l_atoms = len(atoms)
    l_bonds = len(bonds)

    atoms_2, bonds_2 = Fragmenter._get_ring_and_fgroups(
        molecule, functional_groups, ring_systems, atoms, bonds
    )

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

    molecule, _, functional_groups, ring_systems = Fragmenter._prepare_molecule(
        smiles_to_molecule(input_smiles, True), default_functional_groups(), False
    )

    bond = tuple(
        get_map_index(molecule, i)
        for i in molecule.chemical_environment_matches(bond_smarts)[0]
    )

    # noinspection PyTypeChecker
    atoms, bonds = Fragmenter._get_torsion_quartet(molecule, bond)
    atoms, bonds = Fragmenter._get_ring_and_fgroups(
        molecule, functional_groups, ring_systems, atoms, bonds
    )

    actual_atoms = {
        map_index
        for map_index in atoms
        if molecule.atoms[get_atom_index(molecule, map_index)].atomic_number != 1
    }
    expected_atoms = {
        get_map_index(molecule, atom_index)
        for match in molecule.chemical_environment_matches(expected_pattern)
        for atom_index in match
    }

    assert actual_atoms == expected_atoms


def test_find_ortho_substituents(dasatanib):

    # TODO: add test to test adding substituent directly bonded to rotatable bond
    ortho_atoms, ortho_bonds = Fragmenter._find_ortho_substituents(
        dasatanib,
        bonds=smarts_set_to_map_indices({"[#6ar6:1]-[#7:2]-[#6]=[#8]"}, dasatanib),
    )

    assert ortho_atoms == smarts_set_to_map_indices(
        {"[#6ar6]-[#7:1]-[#6]=[#8]", "[#6ar6]-[#17:1]", "[#6H3:1]-[#6r6]:[#6r6]"},
        dasatanib,
    )
    assert ortho_bonds == smarts_set_to_map_indices(
        {"[#6ar6:1]-[#7:2]-[#6]=[#8]", "[#6ar6:1]-[#17:2]", "[#6H3:2]-[#6r6:1]:[#6r6]"},
        dasatanib,
    )


def test_cap_open_valance():

    molecule, _, functional_groups, ring_systems = Fragmenter._prepare_molecule(
        smiles_to_molecule("CNCCc1ccccc1", True), default_functional_groups(), False
    )

    expected_atom = get_map_index(
        molecule,
        molecule.chemical_environment_matches("[#7]-[#6H3:1]")[0][0],
    )

    # noinspection PyTypeChecker
    atoms, bonds = Fragmenter._get_torsion_quartet(
        molecule,
        tuple(
            get_map_index(molecule, i)
            for i in molecule.chemical_environment_matches("[#6a:1]-[#6H2:2]")[0]
        ),
    )
    atoms, bonds = Fragmenter._get_ring_and_fgroups(
        molecule, functional_groups, ring_systems, atoms, bonds
    )

    # Remove the cap atom from the current list to make sure it gets included during
    # capping.
    atoms -= {expected_atom}

    atoms, _ = Fragmenter._cap_open_valence(molecule, functional_groups, atoms, bonds)

    # Check that carbon bonded to N was added
    assert expected_atom in atoms


def test_prepare_molecule():

    molecule, stereo, functional_groups, ring_systems = Fragmenter._prepare_molecule(
        Molecule.from_smiles("C[C@H](O)C1CCCC1"), {"alcohol": "[#6:1]-[#8X2H1:2]"}, True
    )

    assert "atom_map" in molecule.properties

    assert {*stereo.values()} == {"S"}
    assert len(stereo) == 1

    assert "alcohol_0" in functional_groups
    assert len(functional_groups["alcohol_0"][0]) == 2
    assert len(functional_groups["alcohol_0"][1]) == 1

    assert len(ring_systems) == 1
    assert len(ring_systems[1][0]) == 5


@pytest.mark.parametrize(
    "toolkit_registry, expected_provenance",
    [
        (
            None,
            [
                toolkit.__class__.__name__
                for toolkit in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits
            ],
        ),
        (RDKitToolkitWrapper(), ["RDKitToolkitWrapper"]),
        (ToolkitRegistry([RDKitToolkitWrapper]), ["RDKitToolkitWrapper"]),
    ],
)
def test_fragmenter_provenance(toolkit_registry, expected_provenance):
    class DummyFragmenter(Fragmenter):
        def _fragment(
            self, molecule: Molecule, target_bond_smarts: str
        ) -> FragmentationResult:

            return FragmentationResult(
                parent_smiles="[He:1]", fragments=[], provenance={}
            )

    result = DummyFragmenter().fragment(
        Molecule.from_smiles("[He]"), ["[*:1]~[*:2]"], toolkit_registry
    )

    assert "toolkits" in result.provenance
    assert [name for name, _ in result.provenance["toolkits"]] == expected_provenance

    assert "options" in result.provenance
    assert result.provenance["options"]["target_bond_smarts"] == ["[*:1]~[*:2]"]


def test_wbo_fragment():
    """Test build fragment"""

    result = WBOFragmenter().fragment(Molecule.from_smiles("CCCCC"))

    assert len(result.fragments) == 2
    assert {fragment.bond_indices for fragment in result.fragments} == {(3, 5), (4, 5)}


def test_keep_track_of_map():

    result = WBOFragmenter().fragment(Molecule.from_smiles("CCCCC"))

    assert all(
        "atom_map" in fragment.molecule.properties for fragment in result.fragments
    )


def test_get_rotor_wbo():

    molecule = smiles_to_molecule("CCCC", True)

    for bond in molecule.bonds:
        bond.fractional_bond_order = 0.986

    expected_bonds = {
        (
            get_map_index(molecule, match[0]),
            get_map_index(molecule, match[1]),
        )
        for match in molecule.chemical_environment_matches("[#6:1]-[#6:2]")
    }

    rotors_wbo = WBOFragmenter._get_rotor_wbo(
        molecule, WBOFragmenter.find_rotatable_bonds(molecule, None)
    )

    assert len(rotors_wbo) == 1

    rotor_index = next(iter(rotors_wbo))

    assert rotor_index in expected_bonds
    assert numpy.isclose(rotors_wbo[rotor_index], 0.986, atol=0.001)


def test_compare_wbo():

    parent = assign_elf10_am1_bond_orders(smiles_to_molecule("CCCC", add_atom_map=True))
    rotors_wbo = WBOFragmenter._get_rotor_wbo(
        parent, WBOFragmenter.find_rotatable_bonds(parent, None)
    )

    fragment = smiles_to_molecule("CCCC", add_atom_map=True)

    for bond_tuple, value in rotors_wbo.items():

        assert numpy.isclose(
            WBOFragmenter._compare_wbo(fragment, bond_tuple, value),
            0.0,
            atol=1.0e-6,
        )
        assert numpy.isclose(
            WBOFragmenter._compare_wbo(fragment, bond_tuple, value + 1.0),
            1.0,
            atol=1.0e-6,
        )


def test_build_fragment():

    parent = assign_elf10_am1_bond_orders(smiles_to_molecule("CCCCCC", True))
    rotors_wbo = WBOFragmenter._get_rotor_wbo(
        parent, WBOFragmenter.find_rotatable_bonds(parent, None)
    )

    fragments = {
        bond: WBOFragmenter._build_fragment(
            parent, {}, {}, {}, bond, parent_wbo, threshold=0.05
        )
        for bond, parent_wbo in rotors_wbo.items()
    }

    assert len(fragments) == 3

    assert (
        fragments[(3, 5)].to_smiles(explicit_hydrogens=False, mapped=False) == "CCCCC"
    )
    assert (
        fragments[(4, 6)].to_smiles(explicit_hydrogens=False, mapped=False) == "CCCCC"
    )
    assert (
        fragments[(5, 6)].to_smiles(explicit_hydrogens=False, mapped=False) == "CCCCCC"
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

    parent = smiles_to_molecule(input_smiles, True)
    parent_groups = Fragmenter._find_functional_groups(
        parent, default_functional_groups()
    )

    parent_rings = Fragmenter._find_ring_systems(parent, parent_groups)

    assert len(parent_rings[1][0]) == n_output


def test_add_substituent():

    result = WBOFragmenter().fragment(Molecule.from_smiles("CCCCCC"))

    fragment = result.fragments_by_bond[(3, 5)].molecule

    assert fragment.to_smiles(mapped=False, explicit_hydrogens=False) == "CCCCC"

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

    fragment = WBOFragmenter._add_next_substituent(
        result.parent_molecule, {}, {}, {}, atoms, bonds, target_bond=(3, 5)
    )

    assert fragment.to_smiles(mapped=False, explicit_hydrogens=False) == "CCCCCC"


@pytest.mark.parametrize("input_smiles, n_output", [("CCCCCCC", 4)])
def test_pfizer_fragmenter(input_smiles, n_output):

    result = PfizerFragmenter().fragment(Molecule.from_smiles(input_smiles))
    assert len(result.fragments) == n_output
