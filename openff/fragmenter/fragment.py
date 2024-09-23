import abc
import logging
import warnings
from collections import defaultdict
from typing import Any, Literal, Union

import networkx
import openff.fragmenter
from openff.fragmenter.chemi import (
    assign_elf10_am1_bond_orders,
    extract_fragment,
    find_ring_systems,
    find_stereocenters,
)
from openff.fragmenter.states import _enumerate_stereoisomers
from openff.fragmenter.utils import (
    default_functional_groups,
    get_atom_index,
    get_map_index,
    global_toolkit_registry,
)
from openff.toolkit.topology import Atom, Molecule
from openff.toolkit.utils import (
    GLOBAL_TOOLKIT_REGISTRY,
    ToolkitRegistry,
    ToolkitWrapper,
)
from openff.toolkit.utils.exceptions import AtomMappingWarning

from pydantic.v1 import BaseModel, Field

logger = logging.getLogger(__name__)

BondTuple = tuple[int, int]
AtomAndBondSet = tuple[set[int], set[BondTuple]]

Stereochemistries = dict[Union[int, BondTuple], str]
RingSystems = dict[int, AtomAndBondSet]

FunctionalGroups = dict[str, AtomAndBondSet]

Heuristic = Literal["path_length", "wbo"]


class Fragment(BaseModel):
    """An object which stores minimal information about a molecules fragment."""

    smiles: str = Field(
        ...,
        description="A mapped SMILES pattern describing the fragment. The map indices "
        "assigned to each atom in the pattern will correspond to the map index of the "
        "corresponding parent atom. If an atom does not have a map index it is likely "
        "that the atom was added (either H, or C) to ensure every atom in the fragment "
        "has a sensible valence.",
    )

    bond_indices: tuple[int, int] = Field(
        ...,
        description="The map indices of the atoms involved in the bond that the "
        "fragment was built around.",
    )

    @property
    def molecule(self) -> Molecule:
        """The fragment represented as an OpenFF molecule object."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", AtomMappingWarning)

            return Molecule.from_smiles(self.smiles, allow_undefined_stereo=True)


class FragmentationResult(BaseModel):
    """An object which stores the results of fragmenting a molecule."""

    parent_smiles: str = Field(
        ...,
        description="A mapped SMILES pattern describing the parent molecule that was "
        "fragmented.",
    )

    fragments: list[Fragment] = Field(..., description="The generated fragments.")

    provenance: dict[str, Any] = Field(
        ...,
        description="A dictionary storing provenance information about how the "
        "fragments were generated.",
    )

    @property
    def parent_molecule(self) -> Molecule:
        """The parent molecule represented as an OpenFF molecule object."""
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", AtomMappingWarning)

            return Molecule.from_smiles(self.parent_smiles)

    @property
    def fragment_molecules(self) -> dict[BondTuple, Molecule]:
        """A dictionary of the fragment molecules represented as OpenFF molecule
        objects, indexed by the map indices of the bond that each fragment was built
        around."""
        return {fragment.bond_indices: fragment.molecule for fragment in self.fragments}

    @property
    def fragments_by_bond(self) -> dict[BondTuple, Fragment]:
        """Returns a dictionary of fragments indexed by the bond (defined in terms of
        the map indices of the atoms that form it) that the fragment was built around.
        """
        return {fragment.bond_indices: fragment for fragment in self.fragments}


class Fragmenter(BaseModel, abc.ABC):
    """The base class that all fragmentation engines should inherit from."""

    functional_groups: dict[str, str] = Field(
        default_factory=default_functional_groups,
        description="A dictionary of SMARTS of functional groups that should not be "
        "fragmented, indexed by an informative name, e.g. 'alcohol': '[#6]-[#8X2H1]'.",
    )

    @classmethod
    def _find_stereo(cls, molecule: Molecule) -> Stereochemistries:
        """Find chiral atoms and bonds, store the chirality.

        Notes
        -----
            * This is needed to check if a fragment has flipped chirality. Currently
              this can happen and it is a bug.

        Parameters
        ----------
        molecule
            The molecule to search for stereochemistry.

        Returns
        -------
            The stereochemistry associated with atom and bond stereocenters
        """

        atom_stereo = {
            get_map_index(molecule, atom.molecule_atom_index): atom.stereochemistry
            for atom in molecule.atoms
            if atom.stereochemistry is not None
        }

        bond_stereo = {
            (
                get_map_index(molecule, bond.atom1_index),
                get_map_index(molecule, bond.atom2_index),
            ): bond.stereochemistry
            for bond in molecule.bonds
            if bond.stereochemistry is not None
        }

        return {**atom_stereo, **bond_stereo}

    @classmethod
    def _check_stereo(
        cls, fragment: Molecule, parent_stereo: Stereochemistries
    ) -> bool:
        """Checks if the stereochemistry of a fragment is different to the
        stereochemistry of the parent.

        Parameters
        ----------
        fragment
            The fragment whose stereo should be compared to the parent.
        parent_stereo
            The stereochemistry of the parent molecule.

        Returns
        -------
            Whether the fragment has the same stereochemistry as the parent.
        """

        atom_stereocenters, bond_stereocenters = find_stereocenters(fragment)

        # Check for new / flipped chiral centers.
        for atom_index in atom_stereocenters:
            map_index = get_map_index(fragment, atom_index)

            if map_index not in parent_stereo:
                logger.warning(f"A new stereocenter formed at atom {map_index}")
                return False

            fragment_stereo = fragment.atoms[atom_index].stereochemistry

            if fragment_stereo != parent_stereo[map_index]:
                logger.warning(
                    f"Stereochemistry for atom {map_index} flipped from "
                    f"{parent_stereo[map_index]} to {fragment_stereo}"
                )

                return False

        for index_tuple in bond_stereocenters:
            map_tuple = tuple(get_map_index(fragment, i) for i in index_tuple)

            map_tuple = (
                map_tuple if map_tuple in parent_stereo else tuple(reversed(map_tuple))
            )

            if map_tuple not in parent_stereo:
                logger.warning(f"A new chiral bond formed at bond {map_tuple}")
                return False

            fragment_stereo = fragment.get_bond_between(*index_tuple).stereochemistry

            if fragment_stereo != parent_stereo[map_tuple]:
                logger.warning(
                    f"Stereochemistry for bond {map_tuple} flipped from "
                    f"{parent_stereo[map_tuple]} to {fragment_stereo}"
                )

                return False

        return True

    @classmethod
    def _fix_stereo(
        cls, fragment: Molecule, parent_stereo: Stereochemistries
    ) -> Molecule | None:
        """Flip all stereocenters and find the stereoisomer that matches the parent

        Parameters
        ----------
        fragment
            The fragment whose stereochemistry should be flipped to match the parent.
        parent_stereo
            The stereochemistry of the parent molecule.

        Returns
        -------
            A fragment with the same stereochemistry as parent molecule if possible else ``None``.
        """

        for stereoisomer in _enumerate_stereoisomers(fragment, 200, True):
            if not cls._check_stereo(stereoisomer, parent_stereo):
                continue

            return stereoisomer

        return None

    @classmethod
    def _find_functional_groups(
        cls, molecule: Molecule, functional_groups: dict[str, str]
    ) -> FunctionalGroups:
        """Find the atoms and bonds involved in the functional groups specified by
        ``functional_groups``.

        Parameters
        ----------
        molecule
            The molecule to search for function groups.
        functional_groups
            A dictionary of SMARTS of functional groups that should not be fragmented
            indexed by a friendly string label, e.g. 'alcohol: [#6:1]-[#8H1X2:2]'

        Returns
        -------
            The atoms and bonds in the found function groups stored in a dictionary
            indexed by a unique key associated with each functional group.
        """

        found_groups = {}

        for functional_group, smarts in functional_groups.items():
            unique_matches = {
                tuple(sorted(match))
                for match in molecule.chemical_environment_matches(smarts)
            }

            for i, match in enumerate(unique_matches):
                atoms = {get_map_index(molecule, index) for index in match}
                bonds = {
                    (
                        get_map_index(molecule, bond.atom1_index),
                        get_map_index(molecule, bond.atom2_index),
                    )
                    for bond in molecule.bonds
                    if bond.atom1_index in match and bond.atom2_index in match
                }

                found_groups[f"{functional_group}_{i}"] = (atoms, bonds)

        return found_groups

    @classmethod
    def find_rotatable_bonds(
        cls, molecule: Molecule, target_bond_smarts: list[str] | None
    ) -> list[BondTuple]:
        """Finds the rotatable bonds in a molecule *including* rotatable double
        bonds.

        Parameters
        ----------
        molecule
            The molecule to search for rotatable bonds.
        target_bond_smarts
            An optional list of SMARTS patterns that should be used to identify the bonds
            within the parent molecule to grow fragments around. Each SMARTS pattern
            should include **two** indexed atoms that correspond to the two atoms
            involved in the central bond.

            If no pattern is provided fragments will be constructed around all 'rotatable
            bonds'. A 'rotatable bond' here means any bond matched by a
            `[!$(*#*)&!D1:1]-,=;!@[!$(*#*)&!D1:2]` SMARTS pattern with the added
            constraint that the **heavy** degree (i.e. the degree excluding hydrogen) of
            both atoms in the bond must be >= 2.

        Returns
        -------
            A list of the **map** indices of the atoms that form the rotatable
            bonds, ``[(m1, m2),...]``.
        """

        if target_bond_smarts is None:
            matches = molecule.chemical_environment_matches(
                "[!$(*#*)&!D1:1]-,=;!@[!$(*#*)&!D1:2]"
            )

        else:
            matches = [
                match
                for smarts in target_bond_smarts
                for match in molecule.chemical_environment_matches(smarts)
            ]

            if not all(len(match) == 2 for match in matches):
                raise ValueError(
                    f"The `target_bond_smarts` pattern ({target_bond_smarts}) "
                    f"must define exactly two indexed atoms to match."
                )

        unique_matches = {tuple(sorted(match)) for match in matches}

        if target_bond_smarts is None:
            # Drop bonds without a heavy degree of at least 2 on each end to avoid
            # finding terminal bonds
            def heavy_degree(atom_index: int) -> int:
                atom = molecule.atoms[atom_index]
                return sum(1 for atom in atom.bonded_atoms if atom.atomic_number != 1)

            unique_matches = {
                match
                for match in unique_matches
                if all(heavy_degree(i) > 1 for i in match)
            }

        return [
            (
                get_map_index(molecule, match[0]),
                get_map_index(molecule, match[1]),
            )
            for match in unique_matches
        ]

    @classmethod
    def _atom_bond_set_to_mol(
        cls,
        parent: Molecule,
        parent_stereo: Stereochemistries,
        atoms: set[int],
        bonds: set[BondTuple],
    ) -> tuple[Molecule, bool]:
        """Extracts a subset of a molecule based on a set of atom and bond indices.

        Parameters
        ----------
        parent
            The parent molecule to slice.
        parent_stereo
            The stereochemistry of the parent. This will be used to ensure the returned
            subset molecule retains the correct stereochemistry.
        atoms
            set of map indices
        bonds
            set of bond tuples (m1, m2)

        Returns
        -------
            The subset molecule and a flag to indicate whether any new stereocenters are present in the fragment.
        """

        fragment = extract_fragment(parent, atoms, bonds)

        if not cls._check_stereo(fragment, parent_stereo):
            fixed_fragment = cls._fix_stereo(fragment, parent_stereo)
            if fixed_fragment is None:
                return fragment, True
            else:
                return fixed_fragment, False

        return fragment, False

    @classmethod
    def _get_torsion_quartet(
        cls, molecule: Molecule, bond: BondTuple
    ) -> AtomAndBondSet:
        """Get all atoms bonded to the torsion quartet around rotatable bond

        Parameters
        ----------
        molecule
            The molecule containing the rotatable bond.
        bond
            map indices of atoms in bond

        Returns
        -------
            The map indices of atoms in quartet and the bonds in quartet.
        """

        atom_map_indices = {*bond}
        bond_map_indices = {bond}

        atoms = [
            molecule.atoms[i]
            for i, j in molecule.properties["atom_map"].items()
            if j in bond
        ]

        for atom in atoms:
            map_index = get_map_index(molecule, atom.molecule_atom_index)

            for neighbor in atom.bonded_atoms:
                neighbour_map_index = get_map_index(
                    molecule, neighbor.molecule_atom_index
                )

                atom_map_indices.add(neighbour_map_index)
                bond_map_indices.add((map_index, neighbour_map_index))

                for next_neighbour in neighbor.bonded_atoms:
                    next_neighbour_map_index = get_map_index(
                        molecule, next_neighbour.molecule_atom_index
                    )

                    atom_map_indices.add(next_neighbour_map_index)
                    bond_map_indices.add(
                        (neighbour_map_index, next_neighbour_map_index)
                    )

        return atom_map_indices, bond_map_indices

    @classmethod
    def _find_ring_systems(
        cls,
        molecule: Molecule,
        functional_groups: FunctionalGroups,
        keep_non_rotor_ring_substituents: bool = False,
    ) -> RingSystems:
        """This function finds all ring systems in a molecule.

        Parameters
        ----------
        molecule
            The molecule to search for ring systems.
        functional_groups
            A dictionary of the functional groups on the molecule which should not
            be fragmented.
        keep_non_rotor_ring_substituents
            If True, keep all non rotatable ring substituents. According to the
            benchmark, it is not necessary.

        Returns
        -------
            Any found ring systems.
        """

        atom_to_ring_indices = find_ring_systems(molecule)

        # Find the map indices of the atoms involved in each ring system.
        ring_system_atoms = {
            ring_index: {
                get_map_index(molecule, i)
                for i in atom_to_ring_indices
                if atom_to_ring_indices[i] == ring_index
            }
            for ring_index in {*atom_to_ring_indices.values()}
        }

        # Find the map indices of the bonds involved in each ring system.
        ring_system_bonds = defaultdict(set)

        for bond in molecule.bonds:
            ring_index_1 = atom_to_ring_indices.get(bond.atom1_index, -1)
            ring_index_2 = atom_to_ring_indices.get(bond.atom2_index, -2)

            if ring_index_1 != ring_index_2:
                continue

            ring_system_bonds[ring_index_1].add(
                (
                    get_map_index(molecule, bond.atom1_index),
                    get_map_index(molecule, bond.atom2_index),
                )
            )

        # Scan the neighbours of the ring system atoms for any functional groups
        # / non-rotor substituents which should be included in the ring systems.
        for ring_index in ring_system_atoms:
            # If any atoms are part of a functional group, include the other atoms in the
            # group in the ring system lists
            ring_functional_groups = {
                functional_group
                for map_index in ring_system_atoms[ring_index]
                for functional_group in functional_groups
                if map_index in functional_groups[functional_group][0]
            }

            ring_system_atoms[ring_index].update(
                map_index
                for functional_group in ring_functional_groups
                for map_index in functional_groups[functional_group][0]
            )
            ring_system_bonds[ring_index].update(
                map_tuple
                for functional_group in ring_functional_groups
                for map_tuple in functional_groups[functional_group][1]
            )

            if not keep_non_rotor_ring_substituents:
                continue

            non_rotor_atoms, non_rotor_bonds = cls._find_non_rotor_ring_substituents(
                molecule, ring_system_atoms[ring_index]
            )

            ring_system_atoms[ring_index].update(non_rotor_atoms)
            ring_system_bonds[ring_index].update(non_rotor_bonds)

        ring_systems = {
            ring_index: (
                ring_system_atoms[ring_index],
                ring_system_bonds[ring_index],
            )
            for ring_index in ring_system_atoms
        }

        return ring_systems

    @classmethod
    def _find_non_rotor_ring_substituents(
        cls, molecule: Molecule, ring_system_atoms: set[int]
    ) -> AtomAndBondSet:
        """Find the non-rotor substituents attached to a particular ring system.

        Parameters
        ----------
        molecule
            The molecule to search for non-rotor ring substituents.
        ring_system_atoms
            The map indices of the atoms in the ring system of interest.

        Returns
        -------
            The map indices of the atoms and bonds involved in any found
            functional groups.
        """

        rotatable_bonds = molecule.find_rotatable_bonds()

        def heavy_degree(atom: Atom) -> int:
            return sum(1 for atom in atom.bonded_atoms if atom.atomic_number != 1)

        rotor_bonds = [
            bond
            for bond in rotatable_bonds
            if heavy_degree(bond.atom1) >= 2 and heavy_degree(bond.atom2) >= 2
        ]

        non_rotor_atoms = set()
        non_rotor_bonds = set()

        for bond in molecule.bonds:
            # Check if the bond is a rotor.
            if bond in rotor_bonds:
                continue

            if bond.atom1.atomic_number == 1 or bond.atom2.atomic_number == 1:
                continue

            map_index_1 = get_map_index(molecule, bond.atom1_index)
            map_index_2 = get_map_index(molecule, bond.atom2_index)

            in_system_1 = map_index_1 in ring_system_atoms
            in_system_2 = map_index_2 in ring_system_atoms

            if (in_system_1 and in_system_2) or (not in_system_1 and not in_system_2):
                continue

            non_rotor_atoms.update((map_index_1, map_index_2))
            non_rotor_bonds.add((map_index_1, map_index_2))

        return non_rotor_atoms, non_rotor_bonds

    @classmethod
    def _get_ring_and_fgroups(
        cls,
        parent: Molecule,
        parent_groups: FunctionalGroups,
        parent_rings: RingSystems,
        atoms: set[int],
        bonds: set[BondTuple],
    ) -> AtomAndBondSet:
        """Adds the atom and bond indices of groups ortho to those already in
        a fragment, such that they are retained during fragmentation.

        Parameters
        ----------
        parent
            The molecule being fragmented.
        parent_groups
            A dictionary of the functional groups on the molecule which should not
            be fragmented.
        parent_rings
            A dictionary of the ring systems in the molecule which should not
            be fragmented.
        atoms
            The map indices of the atoms in the fragment.
        bonds
            The map indices of the bonds in the fragment.

        Returns
        -------
            The updated set of atom and bond map indices to retain.
        """

        # Find the sets of atoms which are located ortho to one of the bonds being
        # fragmented.
        ortho_atoms, ortho_bonds = cls._find_ortho_substituents(parent, bonds)

        atoms.update(ortho_atoms)
        bonds.update(ortho_bonds)

        # Include the rings systems and functional groups connected to the current
        # atom sets.
        new_atoms = set()
        new_bonds = set()

        fragment_groups = {
            group
            for group in parent_groups
            if any(atom in parent_groups[group][0] for atom in atoms)
        }

        for functional_group in fragment_groups:
            new_atoms.update(parent_groups[functional_group][0])
            new_bonds.update(parent_groups[functional_group][1])

        fragment_rings = {
            ring_index
            for ring_index in parent_rings
            if any(atom in parent_rings[ring_index][0] for atom in atoms)
        }

        for ring_system in fragment_rings:
            new_atoms.update(parent_rings[ring_system][0])
            new_bonds.update(parent_rings[ring_system][1])

        atoms.update(new_atoms)
        bonds.update(new_bonds)

        # Ensure the matched bonds doesn't include duplicates.
        bonds = {tuple(sorted(bond)) for bond in bonds}

        return atoms, bonds

    @classmethod
    def _find_ortho_substituents(
        cls, parent: Molecule, bonds: set[BondTuple]
    ) -> AtomAndBondSet:
        """Find ring substituents that are ortho to one of the rotatable bonds specified
        in a list of bonds.

        Parameters
        ----------
        parent
            The parent molecule being fragmented.
        bonds
            The map indices of the rotatable bonds.

        Returns
        -------
            The set of map indices of atoms in ortho group and of bond tuples in ortho
            group.
        """

        matched_atoms = set()
        matched_bonds = set()

        for match in parent.chemical_environment_matches(
            "[!#1:1]~&!@[*:2]@[*:3]~&!@[!#1*:4]"
        ):
            map_tuple = tuple(get_map_index(parent, i) for i in match)

            if map_tuple[:2] not in bonds and map_tuple[:2][::-1] not in bonds:
                continue

            matched_atoms.update(map_tuple[::3])
            matched_bonds.update((map_tuple[i], map_tuple[i + 1]) for i in [0, 2])

        # Ensure the matched bonds doesn't include duplicates.
        matched_bonds = {tuple(sorted(bond)) for bond in matched_bonds}

        return matched_atoms, matched_bonds

    @classmethod
    def _cap_open_valence(
        cls,
        parent: Molecule,
        parent_groups: FunctionalGroups,
        atoms: set[int],
        bonds: set[BondTuple],
    ) -> AtomAndBondSet:
        """Cap with methyl for fragments that ends with N, O or S. Otherwise cap with H

        Parameters
        ----------
        parent
            The molecule being fragmented.
        parent_groups
            A dictionary of the functional groups on the molecule which should not
            be fragmented.
        atoms
            The map indices of the atoms in the fragment being constructed.
        bonds
            The map indices of the bonds in the fragment being constructed.
        """

        map_index_to_functional_group = {
            map_index: functional_group
            for functional_group in parent_groups
            for map_index in parent_groups[functional_group][0]
        }

        atoms_to_add = set()
        bonds_to_add = set()

        for map_index in atoms:
            atom_index = get_atom_index(parent, map_index)
            atom = parent.atoms[atom_index]

            if (
                atom.atomic_number not in (7, 8, 16)
                and map_index not in map_index_to_functional_group
            ):
                continue

            # If atom is N, O or S, it needs to be capped
            should_cap = False

            for neighbour in atom.bonded_atoms:
                neighbour_map_index = get_map_index(
                    parent, neighbour.molecule_atom_index
                )

                if neighbour.atomic_number == 1 or neighbour_map_index in atoms:
                    continue

                should_cap = True
                break

            if not should_cap:
                continue

            for neighbour in atom.bonded_atoms:
                if neighbour.atomic_number != 6:
                    continue

                neighbour_map_index = get_map_index(
                    parent, neighbour.molecule_atom_index
                )

                atoms_to_add.add(neighbour_map_index)
                bonds_to_add.add((map_index, neighbour_map_index))

        atoms.update(atoms_to_add)
        bonds.update(bonds_to_add)

        return atoms, bonds

    @classmethod
    def _prepare_molecule(
        cls,
        molecule: Molecule,
        functional_groups: dict[str, str],
        keep_non_rotor_ring_substituents: bool,
    ) -> tuple[Molecule, Stereochemistries, FunctionalGroups, RingSystems]:
        """Prepare a molecule for fragmentation.

        This involves canonically ordering the molecule, determining the stereochemistry
        of any stereocenters, detecting any functional groups which should be preserved,
        and finding any ring systems which should be preserved.

        Parameters
        ----------
        molecule
            The parent molecule that should be fragmented.
        functional_groups:
            A dictionary of SMARTS of functional groups that should not be fragmented.
        keep_non_rotor_ring_substituents:
            If True, will always keep all non rotor substituents on ring.

        Returns
        -------
            The prepared molecule to fragment, its stereochemistry, and the functional
            groups and ring systems it contains that should not be fragmented.
        """

        # Canonically order the molecule to try and make the fragmentation more
        # deterministic.
        molecule: Molecule = molecule.canonical_order_atoms()
        molecule.properties["atom_map"] = {i: i + 1 for i in range(molecule.n_atoms)}

        # Keep track of stereo to make sure it does not flip
        stereo = cls._find_stereo(molecule)

        # Find the functional groups and ring systems which should not be fragmented.
        found_functional_groups = cls._find_functional_groups(
            molecule, functional_groups
        )
        found_ring_systems = cls._find_ring_systems(
            molecule, found_functional_groups, keep_non_rotor_ring_substituents
        )

        return molecule, stereo, found_functional_groups, found_ring_systems

    @abc.abstractmethod
    def _fragment(
        self,
        molecule: Molecule,
        target_bond_smarts: list[str] | None,
    ) -> FragmentationResult:
        """The internal implementation of ``fragment``.

        Parameters
        ----------
        molecule
            The molecule to fragment.
        target_bond_smarts
            An optional SMARTS pattern that should be used to identify the bonds within
            the parent molecule to grow fragments around.

        Returns
        -------
            The results of the fragmentation including the fragments and provenance
            about the fragmentation.
        """

        raise NotImplementedError()

    def fragment(
        self,
        molecule: Molecule,
        target_bond_smarts: list[str] | None = None,
        toolkit_registry: ToolkitRegistry | ToolkitWrapper | None = None,
    ) -> FragmentationResult:
        """Fragments a molecule according to this class' settings.

        Notes
        -----
        * This method is currently *not* guaranteed to be thread safe as it uses and
          modifies the OpenFF toolkits' ``GLOBAL_TOOLKIT_REGISTRY``.

        Parameters
        ----------
        molecule
            The molecule to fragment.
        target_bond_smarts
            An optional list of SMARTS patterns that should be used to identify the bonds
            within the parent molecule to grow fragments around. Each SMARTS pattern
            should include **two** indexed atoms that correspond to the two atoms
            involved in the central bond.

            If no pattern is provided fragments will be constructed around all 'rotatable
            bonds'. A 'rotatable bond' here means any bond matched by a
            `[!$(*#*)&!D1:1]-,=;!@[!$(*#*)&!D1:2]` SMARTS pattern with the added
            constraint that the **heavy** degree (i.e. the degree excluding hydrogen) of
            both atoms in the bond must be >= 2. Note this will not find terminal
            rotatable bonds such as -OH, -NH2 -CH3.
        toolkit_registry
            The underlying cheminformatics toolkits to use for things like conformer
            generation, WBO computation etc. If no value is provided, the current
            ``GLOBAL_TOOLKIT_REGISTRY`` will be used. See the OpenFF toolkit
            documentation for more information.

        Returns
        -------
            The results of the fragmentation including the fragments and provenance
            about the fragmentation.
        """

        if toolkit_registry is None:
            toolkit_registry = GLOBAL_TOOLKIT_REGISTRY

        with global_toolkit_registry(toolkit_registry):
            result = self._fragment(molecule, target_bond_smarts)

            result.provenance["toolkits"] = [
                (toolkit.__class__.__name__, toolkit.toolkit_version)
                for toolkit in GLOBAL_TOOLKIT_REGISTRY.registered_toolkits
            ]

        if "options" not in result.provenance:
            result.provenance["options"] = {}

        if target_bond_smarts is not None:
            result.provenance["options"]["target_bond_smarts"] = target_bond_smarts

        return result

    def _default_provenance(self) -> dict[str, Any]:
        """Returns a dictionary containing default provenance information."""

        provenance = {
            "creator": openff.fragmenter.__package__,
            "version": openff.fragmenter.__version__,
            "options": self.dict(),
        }

        return provenance


class WBOOptions(BaseModel):
    """A set of options for controlling how Wiberg Bond Orders are computed."""

    method: Literal["am1-wiberg-elf10"] = Field(
        "am1-wiberg-elf10", description="The method to use when computing the WBOs."
    )

    max_conformers: int = Field(
        800, description="The maximum number of conformers to average the WBOs over."
    )
    rms_threshold: float = Field(
        1.0,
        description="The minimum RMS value [Angstrom] at which two conformers are "
        "considered redundant and one is deleted.",
    )


class WBOFragmenter(Fragmenter):
    """Fragment engine for fragmenting molecules using Wiberg Bond Order."""

    scheme: Literal["WBO"] = "WBO"

    wbo_options: WBOOptions = Field(
        WBOOptions(), description="The options to use when computing the WBOs."
    )

    threshold: float = Field(
        0.03,
        description="The threshold for the central bond WBO. If the fragment WBO is "
        "below this threshold, fragmenter will grow out the fragment one bond at a "
        "time via the path specified by the heuristic option",
    )
    heuristic: Heuristic = Field(
        "path_length",
        description="The path fragmenter should take when fragment needs to be grown "
        "out. The options are ``['wbo', 'path_length']``.",
    )

    keep_non_rotor_ring_substituents: bool = Field(
        False,
        description="Whether to always keep all non rotor substituents on rings. If "
        "``False``, rotor substituents on rings will only be retained if they are "
        "ortho to the central bond or if it's needed for WBO to be within the "
        "threshold.",
    )

    def _fragment(
        self, molecule: Molecule, target_bond_smarts: list[str] | None
    ) -> FragmentationResult:
        """Fragments a molecule in such a way that the WBO of the bond that a fragment
        is being built around does not change beyond the specified threshold.
        """

        (
            molecule,
            stereochemistry,
            functional_groups,
            ring_systems,
        ) = self._prepare_molecule(
            molecule, self.functional_groups, self.keep_non_rotor_ring_substituents
        )

        # Calculate WBO for molecule
        if self.wbo_options.method != "am1-wiberg-elf10":
            raise NotImplementedError(
                "WBOs can currently only be computed using 'am1-wiberg-elf10'."
            )

        molecule = assign_elf10_am1_bond_orders(
            molecule, self.wbo_options.max_conformers, self.wbo_options.rms_threshold
        )

        rotatable_bonds = self.find_rotatable_bonds(molecule, target_bond_smarts)
        wbo_rotor_bonds = self._get_rotor_wbo(molecule, rotatable_bonds)

        fragments = {
            bond: self._build_fragment(
                molecule,
                stereochemistry,
                functional_groups,
                ring_systems,
                bond,
                wbo_rotor_bonds[bond],
                threshold=self.threshold,
                heuristic=self.heuristic,
            )
            for bond in wbo_rotor_bonds
        }

        return FragmentationResult(
            parent_smiles=molecule.to_smiles(mapped=True),
            fragments=[
                Fragment(smiles=fragment.to_smiles(mapped=True), bond_indices=bond)
                for bond, fragment in fragments.items()
            ],
            provenance=self._default_provenance(),
        )

    @classmethod
    def _get_rotor_wbo(
        cls, molecule: Molecule, rotor_bonds: list[BondTuple]
    ) -> dict[BondTuple, float]:
        """Cache the WBO of each bond in a specific set of rotor bonds..

        Parameters
        ----------
            molecule
                The molecule containing the rotors.
            rotor_bonds
                The map indices of the rotor bonds to return the WBOs of.

        Returns
        -------
            The WBO of each rotor bond.
        """

        if any(bond.fractional_bond_order is None for bond in molecule.bonds):
            raise RuntimeError(
                "WBO was not calculated for this molecule. Calculating WBO..."
            )

        rotors_wbo = {}

        for bond_indices in rotor_bonds:
            bond = molecule.get_bond_between(
                get_atom_index(molecule, bond_indices[0]),
                get_atom_index(molecule, bond_indices[1]),
            )

            rotors_wbo[bond_indices] = bond.fractional_bond_order

        return rotors_wbo

    @classmethod
    def _compare_wbo(
        cls, fragment: Molecule, bond_tuple: BondTuple, parent_wbo: float, **kwargs
    ) -> float:
        """Compare Wiberg Bond order of rotatable bond in a fragment to the parent.

        Parameters
        ----------
        fragment
            The fragment containing the rotatable bond.
        bond_tuple
            The map indices of the rotatable bond.
        parent_wbo
            The WBO of the parent bond with map indices matching ``bond_tuple``.

        Returns
        -------
            The absolute difference between the fragment and parent WBOs.
        """

        # Create new fragment object because sometimes the molecule created from atom
        # bond set is wonky and then the WBOs are not reproducible

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", AtomMappingWarning)

            fragment = Molecule.from_smiles(
                fragment.to_smiles(mapped=True), allow_undefined_stereo=True
            )

        fragment_map = fragment.properties.pop("atom_map", None)

        try:
            fragment = assign_elf10_am1_bond_orders(fragment, **kwargs)

        except RuntimeError:
            # Most of the time it fails because it is either missing parameters or a
            # functional group that should not be fragmented was fragmented
            logger.warning(
                f"Cannot calculate WBO for fragment {fragment.to_smiles()}. Continue "
                f"growing fragment"
            )

            # TODO: handle different kinds of failures instead of just continuing to
            #      grow until the failure goes away. Some fail because there are
            #      functional groups that should not be fragmented.
            return 1.0

        if fragment_map is not None:
            fragment.properties["atom_map"] = fragment_map

        bond = fragment.get_bond_between(
            get_atom_index(fragment, bond_tuple[0]),
            get_atom_index(fragment, bond_tuple[1]),
        )

        fragment_wbo = bond.fractional_bond_order

        return abs(parent_wbo - fragment_wbo)

    @classmethod
    def _build_fragment(
        cls,
        parent: Molecule,
        parent_stereo: Stereochemistries,
        parent_groups: FunctionalGroups,
        parent_rings: RingSystems,
        bond_tuple: BondTuple,
        parent_wbo: float,
        threshold: float,
        heuristic: Heuristic = "path_length",
        cap: bool = True,
        **kwargs,
    ) -> Molecule:
        """Build a fragment around a specified bond.

        Parameters
        ----------
        parent
            The original molecule being fragmented.
        parent_stereo
            The stereochemistry of the parent molecule.
        parent_groups
            A dictionary of the functional groups on the molecule which should not
            be fragmented.
        parent_rings
            A dictionary of the ring systems in the molecule which should not
            be fragmented.
        bond_tuple
            The map indices specifying which bond to build the fragment around.
        parent_wbo
            The WBO of the parent bond with map indices matching ``bond_tuple``.
        threshold
            The threshold for the central bond WBO. If the fragment WBO is below this
            threshold, fragmenter will grow out the fragment one bond at a time via the
            path specified by the heuristic option
        heuristic
            The heuristic to use when building the fragment.
        cap
            Whether to cap open valences.

        Returns
        -------
            The built fragment molecule.
        """

        atoms, bonds = cls._get_torsion_quartet(parent, bond_tuple)
        atoms, bonds = cls._get_ring_and_fgroups(
            parent, parent_groups, parent_rings, atoms, bonds
        )

        # Cap open valence
        if cap:
            atoms, bonds = cls._cap_open_valence(parent, parent_groups, atoms, bonds)

        fragment, has_new_stereocenter = cls._atom_bond_set_to_mol(
            parent, parent_stereo, atoms, bonds
        )

        if has_new_stereocenter:
            wbo_difference = threshold + 1.0
        else:
            wbo_difference = cls._compare_wbo(
                fragment, bond_tuple, parent_wbo, **kwargs
            )

        while fragment is not None and wbo_difference > threshold:
            fragment, has_new_stereocenter = cls._add_next_substituent(
                parent,
                parent_stereo,
                parent_groups,
                parent_rings,
                atoms,
                bonds,
                target_bond=bond_tuple,
                heuristic=heuristic,
            )

            if fragment is None:
                break

            if has_new_stereocenter:
                # For now keep growing the fragment as it is not yet clear how to handle such cases robustly.
                wbo_difference = threshold + 1.0

            else:
                wbo_difference = cls._compare_wbo(
                    fragment, bond_tuple, parent_wbo, **kwargs
                )

        # A work around for a known bug where if stereochemistry changes or gets removed,
        # the WBOs can change more than the threshold (this will sometimes happen if a
        # very small threshold is chosen) and even the parent will have a WBO difference
        # greater than the threshold. In this case, return the molecule
        if fragment is None:
            fragment, _ = cls._atom_bond_set_to_mol(parent, parent_stereo, atoms, bonds)

        return fragment

    @classmethod
    def _select_neighbour_by_path_length(
        cls, molecule: Molecule, atoms: set[int], target_bond: BondTuple
    ) -> tuple[int, BondTuple] | None:
        atom_indices = {get_atom_index(molecule, atom) for atom in atoms}

        atoms_to_add = [
            (atom_index, neighbour.molecule_atom_index)
            for atom_index in atom_indices
            for neighbour in molecule.atoms[atom_index].bonded_atoms
            if neighbour.atomic_number != 1
            and neighbour.molecule_atom_index not in atom_indices
        ]
        map_atoms_to_add = [
            (
                get_map_index(molecule, j),
                (get_map_index(molecule, i), get_map_index(molecule, j)),
            )
            for i, j in atoms_to_add
        ]

        # Compute the distance from each neighbouring atom to each of the atoms in the
        # target bond.
        nx_molecule = molecule.to_networkx()

        target_indices = [get_atom_index(molecule, atom) for atom in target_bond]

        path_lengths = list(
            zip(
                *(
                    (
                        networkx.shortest_path_length(
                            nx_molecule, target_index, neighbour_index
                        )
                        for target_index in target_indices
                    )
                    for atom_index, neighbour_index in atoms_to_add
                )
            )
        )

        if len(path_lengths) == 0:
            return None
        path_lengths_1, path_lengths_2 = path_lengths

        reverse = False

        min_path_length_1 = min(path_lengths_1)
        min_path_length_2 = min(path_lengths_2)

        if min_path_length_1 < min_path_length_2:
            sort_by = path_lengths_1
        elif min_path_length_2 < min_path_length_1:
            sort_by = path_lengths_2

        else:
            # If there are multiple neighbouring atoms the same path length away
            # from the target bond fall back to sorting by the WBO.
            map_atoms_to_add = [
                map_tuple
                for map_tuple, *path_length_tuple in zip(
                    map_atoms_to_add, path_lengths_1, path_lengths_2
                )
                if min_path_length_1 in path_length_tuple
            ]

            sort_by = [
                molecule.get_bond_between(
                    get_atom_index(molecule, neighbour_bond[0]),
                    get_atom_index(molecule, neighbour_bond[1]),
                ).fractional_bond_order
                for _, neighbour_bond in map_atoms_to_add
            ]

            reverse = True

        sorted_atoms = [
            a for _, a in sorted(zip(sort_by, map_atoms_to_add), reverse=reverse)
        ]

        return None if len(sorted_atoms) == 0 else sorted_atoms[0]

    @classmethod
    def _select_neighbour_by_wbo(
        cls, molecule: Molecule, atoms: set[int]
    ) -> tuple[int, BondTuple] | None:
        """A function which return those atoms which neighbour those in the ``atoms``
        list sorted by the WBO of the bond between the input atom and the neighbouring
        atom from largest to smallest.

        Parameters
        ----------
        molecule
            The original molecule being fragmented.
        atoms
            The map indices of the atoms currently in the fragment.

        Returns
        -------
            The indices of the atoms to be added to the fragment sorted into the
            order that they should be added in.
        """

        map_bond_orders = {
            (
                get_map_index(molecule, bond.atom1_index),
                get_map_index(molecule, bond.atom2_index),
            ): bond.fractional_bond_order
            for bond in molecule.bonds
            if bond.atom1.atomic_number != 1 and bond.atom2.atomic_number != 1
        }

        neighbour_bond_orders = {
            (bond_order, (map_tuple[1 - i], map_tuple))
            for i in range(2)
            for map_tuple, bond_order in map_bond_orders.items()
            if map_tuple[i] in atoms and map_tuple[1 - i] not in atoms
        }

        sorted_atoms = [
            atom_to_add
            for _, atom_to_add in sorted(
                neighbour_bond_orders, key=lambda x: x[0], reverse=True
            )
        ]

        return None if len(sorted_atoms) == 0 else sorted_atoms[0]

    @classmethod
    def _add_next_substituent(
        cls,
        parent: Molecule,
        parent_stereo: Stereochemistries,
        parent_groups: FunctionalGroups,
        parent_rings: RingSystems,
        atoms: set[int],
        bonds: set[BondTuple],
        target_bond: BondTuple,
        heuristic: Heuristic = "path_length",
    ) -> tuple[Molecule | None, bool]:
        """Expand the fragment to include the next set of substituents / ring systems.

        Parameters
        ----------
        parent
            The original molecule being fragmented.
        parent_stereo
            The stereochemistry of the parent molecule.
        parent_groups
            A dictionary of the functional groups on the molecule which should not
            be fragmented.
        parent_rings
            A dictionary of the ring systems in the molecule which should not
            be fragmented.
        atoms
            The map indices of the atoms currently in the fragment.
        bonds
            The map indices of the bonds currently in the fragment.
        target_bond
            The bond the fragment is being built around.
        heuristic
            How to add substituents. The choices are `'path_length'` or `'wbo'`

        Returns
        -------
            The expanded fragment, or None if the fragment already includes the
            entire parent.
        """

        # Select the next atom neighbour (and the groups / rings that it is part of)
        # that should be added to the fragment.
        if heuristic == "wbo":
            neighbour_atom_and_bond = cls._select_neighbour_by_wbo(parent, atoms)
        elif heuristic == "path_length":
            neighbour_atom_and_bond = cls._select_neighbour_by_path_length(
                parent, atoms, target_bond
            )
        else:
            raise NotImplementedError(
                "Only `'wbo'` and `'path_length'` are supported heuristics."
            )

        if neighbour_atom_and_bond is None:
            return None, False

        neighbour_atom, neighbour_bond = neighbour_atom_and_bond

        # If the neighbour to include is part of a functional group / ring system
        # the entire group should be included in the fragment.
        for group, group_atoms in parent_groups.items():
            if neighbour_atom not in group_atoms[0]:
                continue

            atoms.update(parent_groups[group][0])
            bonds.update(parent_groups[group][1])

        for ring_index, ring_atoms in parent_rings.items():
            if neighbour_atom not in ring_atoms[0]:
                continue

            atoms.update(parent_rings[ring_index][0])
            bonds.update(parent_rings[ring_index][1])

        atoms.add(neighbour_atom)
        bonds.add(neighbour_bond)

        # Check new WBO
        return cls._atom_bond_set_to_mol(parent, parent_stereo, atoms, bonds)


class PfizerFragmenter(Fragmenter):
    """Fragment engine for fragmenting molecules using Pfizer's protocol
    (doi: 10.1021/acs.jcim.9b00373)
    """

    scheme: Literal["Pfizer"] = "Pfizer"

    def _fragment(
        self, molecule: Molecule, target_bond_smarts: list[str] | None
    ) -> FragmentationResult:
        """Fragments a molecule according to Pfizer protocol."""

        (
            parent,
            parent_stereo,
            parent_groups,
            parent_rings,
        ) = self._prepare_molecule(molecule, self.functional_groups, False)

        rotatable_bonds = self.find_rotatable_bonds(parent, target_bond_smarts)

        fragments = {}

        for bond in rotatable_bonds:
            atoms, bonds = self._get_torsion_quartet(parent, bond)

            atoms, bonds = self._get_ring_and_fgroups(
                parent, parent_groups, parent_rings, atoms, bonds
            )

            atoms, bonds = self._cap_open_valence(parent, parent_groups, atoms, bonds)

            fragments[bond], _ = self._atom_bond_set_to_mol(
                parent, parent_stereo, atoms, bonds
            )
        return FragmentationResult(
            parent_smiles=parent.to_smiles(mapped=True),
            fragments=[
                Fragment(smiles=fragment.to_smiles(mapped=True), bond_indices=bond)
                for bond, fragment in fragments.items()
            ],
            provenance=self._default_provenance(),
        )
